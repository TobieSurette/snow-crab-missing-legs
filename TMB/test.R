library(gulf.data)
library(gulf.graphics)
library(TMB)

clm()
clg()

# Read snow crab survey biological data:
b <- read.scsbio(1988:2021)
b$year <- year(b)
b$maturity <- is.mature(b)
b <- b[which((b$carapace.width > 0) & (b$carapace.width <= 150)), ]
b <- b[which(b$carapace.width >= 50 & b$sex == 1), ]

# Format durometer readings:
b$durometer <- gsub("[*]", "", b$durometer)
b$durometer <- as.numeric(b$durometer)
b$durometer[which((b$durometer <= 0) & (b$durometer > 100))] <- NA

# Fix missing legs:
b$missing.legs[b$missing.legs == "*********"] <- "**********"

# Define variables:
years <- sort(unique(b$year))

# Define missing leg patterns:
m <- matrix(0, ncol = 10, nrow = nrow(b))
for (i in 1:10) m[, i] <- as.numeric(substr(b$missing.legs, i, i))
m[is.na(m)] <- 0

# Compile and load TMB program:
clc(); compile("test.cpp")
dyn.load(dynlib("test"))

# Define data:
ix <- sort(sample(1:nrow(m), 20000))
data = list(z = m)

# Define initial parameters:
parameters <- list(leg_effect = rep(0, ncol(m)), # Leg effect parameters.
                   log_sigma_leg = -1,    # Leg effect error parameter.
                   L_eps = rep(0, 0.5 * ncol(m) * (ncol(m)-1)),
                   log_scale_eps = rep(0, ncol(m)),
                   log_sigma_L_eps = -1,
                   log_sigma_scale_eps = -1)

# Initialize TMB object:
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 random = c("leg_effect", "L_eps", "log_scale_eps"),
                 DLL = "test")

# Fit model:
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 500, trace = 3))
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")


# Compile rates of missing legs by size and year:
res <- matrix(NA, nrow = 150, ncol = length(years))
rownames(res) <- 1:150
colnames(res) <- years
for (i in 1:length(years)){
   ix <- which(b$maturity & (b$year == years[i]))
   r <- (m[ix,j] == 1) | (m[ix,j+5] == 1)
   r <- apply(m[ix, ] == 1, 1, function(x) sum(x) / 10)
   tmp <- aggregate(r, by = list(cw = round(b$carapace.width[ix])), mean)
   res[as.character(tmp$cw), as.character(years[i])] <- tmp$x
}
res[, "1996"] <- NA

# Plot rates:
image(1:150, years, res, xlab = "", ylab = "", col = colorRampPalette(c("white", "black"))(10),
      zlim = c(0, 0.15), xlim = c(60, 120))
mtext("Carapace width (mm)", 1, 2.5, cex = 1.5)
mtext("Year", 2, 2.5, cex = 1.5)
box()


logit <- function(x) log(x/(1-x))

# Plot rate anomalies:
cols <- c(rev(colorRampPalette(c("white", "red"))(10)[-1]), "white", colorRampPalette(c("white", "black"))(10)[-1])
image(1:150, years, logit(res) + 2.6, xlab = "", ylab = "",
      col = cols,
      breaks = seq(-1, 1, len = length(cols)+1),
      xlim = c(60, 110))
mtext("Carapace width (mm)", 1, 2.5, cex = 1.5)
mtext("Year", 2, 2.5, cex = 1.5)
box()


