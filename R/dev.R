library(gulf.data)
library(gulf.graphics)
library(TMB)

year <- 2018
sex <- 1

# Compile TMB program:
clc(); compile("dev.cpp")
dyn.load(dynlib("dev"))

# Load biological data:
b <- read.scsbio(year, survey = "regular")
b <- b[which((b$sex == sex) & (b$carapace.width > 0) & (b$carapace.width < 150) & !is.na(b$carapace.width)), ]
b$missing.legs[b$missing.legs == "*********"] <- "**********"
b$year <- year(b)
b$missing.legs <- gsub("[2-7]", "*", b$missing.legs)
b$maturity <- as.numeric(is.mature(b))
b <- b[!is.na(b$maturity) & !is.na(b$carapace.width), ]
b <- b[which(b$carapace.width >= 40), ]

# Missing leg matrix:
m <- matrix(NA, nrow = nrow(b), ncol = 10)
for (i in 1:10) m[,i] <- as.numeric(substr(b$missing.legs,i,i) == "1")

# Define data input:
data <- list(z = m,
             maturity = b$maturity,
             size = round(b$carapace.width))

# Define initial parameter values:
parameters <- list(alpha = -3.7,
                   leg_effect = c(-0.39, 0.38, 0.169, 0.082, -0.047, -0.496,  0.42, 0.20, 0.089, 0.083),
                   log_sigma_leg = -1.52,
                   beta_maturity = 1.09,
                   size_effect = c(0.31,0.31,0.32,0.32,0.33,0.33,0.34,0.34,0.34,0.35,0.35,0.36,0.36,0.37,
                                   0.37,0.38,0.38,0.39,0.39,0.4,0.4,0.41,0.41,0.42,0.43,0.43,0.44,0.44,
                                   0.45,0.46,0.46,0.47,0.47,0.48,0.49,0.49,0.5,0.51,0.51,0.52,0.53,0.59,
                                   0.54,0.53,0.55,0.62,0.66,0.64,0.49,0.38,0.34,0.33,0.39,0.45,0.51,0.58,
                                   0.7,0.77,0.88,0.96,0.97,1.02,1.11,1.13,1.16,1.19,1.27,1.35,1.35,1.39,
                                   1.47,1.5,1.47,1.46,1.54,1.6,1.6,1.61,1.59,1.55,1.52,1.52,1.53,1.41,1.34,
                                   1.32,1.22,1.23,1.19,1.18,1.21,1.22,1.14,1.04,0.94,0.83,0.78,0.7,0.58,0.5,
                                   0.45,0.41,0.32,0.21,0.19,0.16,0.14,0.12,0.12,0.14,0.16,0.21,0.23,0.25,0.29,
                                   0.31,0.34,0.39,0.44,0.51,0.61,0.67,0.69,0.72,0.75,0.78,0.82,0.86,0.88,0.88,
                                   0.87,0.87,0.86,0.86,0.85,0.84,0.84,0.83,0.82,0.8,0.79,0.78,0.77,0.76,0.75,
                                   0.74,0.73,0.72,0.71,0.7),
                   beta_size = c(-0.00031, 0.0893),
                   xp_size = 70.609,
                   logit_rho_size = 2.01,
                   crab_effect = rep(0, nrow(m)),
                   log_sigma_crab = -0.141)


# Create TMB object:
obj <- MakeADFun(data = data, parameters = parameters, DLL = "dev", random = c("leg_effect", "size_effect", "crab_effect"))

# Fit model:
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 1500, trace = 3))
obj$par <- theta$par

p <- obj$report()$p
rep  <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")


ix <- b$carapace.width > 60
plot(random[grep("crab", rownames(random)),1][ix], col = rainbow(5)[b$shell.condition[ix]],
     cex = 0.3 * sqrt(exp(random[grep("crab", rownames(random)),1][ix])))
legend("topleft", legend = 1:5, pch = 21, pt.bg = rainbow(5))

# Model predictions on a grid:
grid <- expand.grid(cw = 1:150, maturity = 0:1, leg = 1:10)
leg_effect <- as.numeric(random[grep("leg_effect", rownames(random)), 1])
size_effect <- as.numeric(random[grep("size_effect", rownames(random)), 1])
crab_effect <- as.numeric(random[grep("crab_effect", rownames(random)), 1])
mu_grid <- fixed["alpha",1] +
      fixed["beta_maturity",1] * grid$maturity +
      leg_effect[grid$leg] + size_effect[grid$cw]
p_grid <- 1 / (1 + exp(-mu_grid))

# Model predictions of obsservations:
mu_obs <- NULL
for (i in 1:10){
   mu_obs <- cbind(mu_obs, fixed["alpha",1] +
                           fixed["beta_maturity",1] * data$maturity +
                           leg_effect[i] + fixed["beta_size",1] * size_effect[data$size])# +
                           #crab_effect)
}
p_obs <- 1 / (1 + exp(-mu_obs))

# ========================================= Diagnostic plots ========================================

# Missing legs by size by maturity:
ix <- data$maturity == 0
tmp <- aggregate(apply(m[ix,], 1, sum), by = list(size = data$size[ix]), mean)
gbarplot(tmp$x, tmp$size)
tmp <- aggregate(apply(p[ix,], 1, sum, na.rm = TRUE), by = list(size = data$size[ix]), mean)
lines(tmp$size, tmp$x, lwd = 2, col = "blue")

# Missing legs by shell condition:
ix <- data$maturity == 1
tmp <- aggregate(apply(m[ix,], 1, sum), by = list(size = b$shell.condition[ix]), mean)
gbarplot(tmp$x, tmp$size)
tmp <- aggregate(apply(p[ix,], 1, sum, na.rm = TRUE), by = list(size = b$shell.condition[ix]), mean)
lines(tmp$size, tmp$x, lwd = 2, col = "blue")

# Missing legs by leg position and maturity:
ix <- data$maturity == 1
tmp <- apply(m[ix,], 2, mean)
gbarplot(tmp, 1:10)
tmp <- apply(p[ix,], 2, mean, na.rm = TRUE)
lines(1:10, tmp, lwd = 2, col = "blue")

# Missing legs by size for each leg position:
ix <- data$maturity == 1
ll <- kronecker(matrix(1:5), matrix(1, nrow = 6, ncol = 6))
ll <- rbind(0, cbind(0,ll,0), 0, 0, 0)
layout(ll)
par(mar = c(0,0,0,0))
for (i in 1:5){
   rm <- aggregate(m[ix, c(0,5) + i], by = list(size = data$size[ix]), mean)
   rp <- aggregate(p[ix, c(0,5) + i], by = list(size = data$size[ix]), mean)
   #robs <- aggregate(p_obs[ix, c(0,5) + i], by = list(size = data$size[ix]), mean)

   plot(c(40, 120), c(0, 0.20), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
   grid()
   lines(rm$size, rm[,2], lwd = 2, col = "blue")
   lines(rp$size, rp[,2], lwd = 2, col = "blue", lty = "dashed")
   lines(rm$size, rm[,3], lwd = 2, col = "green")
   lines(rp$size, rp[,3], lwd = 2, col = "green", lty = "dashed")

   #lines(robs$size, rp[,2], lwd = 2, col = "red", lty = "solid")
   #lines(robs$size, rp[,3], lwd = 2, col = "red", lty = "dashed")

   if (i == 3) mtext("Missing legs (#/crab)", 2,  2.75, cex = 1.25)
   mtext(i, 4, 1.5, cex = 1.25, las = 1)

   axis(2, at = seq(0, 0.15, by = 0.05), cex.axis = 1.25)
   axis(2, at = seq(0.05, 0.15, by = 0.05), cex.axis = 1.25)
   box()
}
axis(1, cex.axis = 1.25)
mtext("Carapace width (mm)", 1, 2.5, cex = 1.25)

#


