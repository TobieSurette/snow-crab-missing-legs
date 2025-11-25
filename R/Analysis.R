library(gulf)

rm(list = ls())
graphics.off()
setwd("U:/Snow Crab/Missing Legs")
load("Survey Missing Leg Data.Rdata")

# Define Gaussian likelihood function:
log.likelihood <- function(theta, x, fun, index = NULL, theta.fixed = NULL){
   if (!is.null(index) & !is.null(theta.fixed)){
      temp <- index * NA
      temp[index] <- theta
      temp[!index] <- theta.fixed
      theta <- temp
   }

   sigma <- theta[length(theta)]
   n <- dim(x)[1]
   res <- (log(x$weight) - fun(theta, x))
   ll <- - 0.5 * log(2*pi) - log(sigma) - ((res)^2 / (2*sigma^2))
   ll <- sum(ll)

   return(ll)
}

# Allometric model with missing leg predictors only:
fun <- function(theta, x){
   beta <- theta[1]
   alpha <- theta[2:7]
   
   mu <- log(alpha[1] + (alpha[2] * x$m1) + (alpha[3] * x$m2) + (alpha[4] * x$m3) + (alpha[5] * x$m4) + (alpha[6] * x$m5))
   mu <- beta * log(x$carapace.width) + mu

   return(mu)
}

# Allometric model with missing leg and regenerated legs predictors:
fun2 <- function(theta, x){
   beta <- theta[1]
   alpha <- theta[c(2, 3:7, 8:12)]

   mu <- alpha[1] + (alpha[2] * x$m1) + (alpha[3] * x$m2) + (alpha[4] * x$m3) +
                    (alpha[5] * x$m4) + (alpha[6] * x$m5)
   mu <- mu + (alpha[7] * x$r1) + (alpha[8] * x$r2) + (alpha[9] * x$r3) +
                    (alpha[10] * x$r4) + (alpha[11] * x$r5)
   mu <- beta * log(x$carapace.width) + log(mu)

   return(mu)
}

# Allometric model with missing leg , regenerated leg and region predictors:
fun3 <- function(theta, x){
   beta <- theta[1]
   alpha <- theta[c(2, 3:7, 8:12, 13)]

   mu <- alpha[1] + (alpha[2] * x$m1) + (alpha[3] * x$m2) + (alpha[4] * x$m3) +
                    (alpha[5] * x$m4) + (alpha[6] * x$m5)
   mu <- mu + (alpha[7] * x$r1) + (alpha[8] * x$r2) + (alpha[9] * x$r3) +
                    (alpha[10] * x$r4) + (alpha[11] * x$r5)
   mu <- mu + alpha[12]*(x$area == "CHETICAMP NS        ")
   mu <- beta * log(x$carapace.width) + log(mu)

   return(mu)
}

# Fit missing leg model to data:
theta <- c(3, 0.0003870332, 0, 0, 0, 0, 0, 0.2)
fun(theta, x)
log.likelihood(theta, x, fun)
l1 <- optim(theta, log.likelihood, control = list(maxit = 5000, trace = 1, fnscale = -1),
            x = x, fun = fun)
l1 <- optim(l1$par, log.likelihood, control = list(maxit = 5000, trace = 1, fnscale = -1),
            x = x, fun = fun)

# Fit missing leg and regenerated leg model to data:
theta <- c(l1$par[1], l1$par[2], l1$par[3:7], 0, 0, 0, 0, 0, l1$par[length(l1$par)])
fun2(theta, x)
log.likelihood(theta, x, fun2)
l2 <- optim(theta, log.likelihood, control = list(maxit = 5000, trace = 1, fnscale = -1),
            x = x, fun = fun2)
l2 <- optim(l2$par, log.likelihood, control = list(maxit = 5000, trace = 1, fnscale = -1),
            x = x, fun = fun2)
            



theta <- c(l2$par[1], l2$par[2], l2$par[3:12], 0, l2$par[length(l2$par)])
fun3(theta, x)
log.likelihood(theta, x, fun3)
l3 <- optim(theta, log.likelihood, control = list(maxit = 5000, trace = 1, fnscale = -1),
            x = x, fun = fun3)
l3 <- optim(l3$par, log.likelihood, control = list(maxit = 5000, trace = 1, fnscale = -1),
            x = x, fun = fun3)
            
# Plot coeffcient values:
windows()
dbarplot(abs(l3$par[3:13]))

windows()
res <- (log(x$weight) - fun3(l3$par, x)) / l3$par[length(l3$par)]
index <- abs(res) < 4
hist(res[abs(res) < 4], n = 100)
x[rownames(x)[abs(res) >= 4], ]

windows()
layout(matrix(1:5, ncol =1))
for (i in 1:4){
   # A numerical vector of the form c(bottom, left, top, right) c(5, 4, 4, 2) + 0.1.
   par(mar = (c(2, 4, 2, 2) + 0.1))
   temp <- index & (x$total == (i - 1))
   hist(res[temp], n = 100, xlab = NA, main = NA, xlim = c(-4, 4))
}
temp <- apply(m == 2, 1, sum) > 0
hist(res[temp], n = 100, xlab = NA, main = NA, xlim = c(-4, 4))

# Residuals trends for #missing legs
windows()
layout(matrix(1:5, ncol =1))
for (i in 1:5){
   # A numerical vector of the form c(bottom, left, top, right) c(5, 4, 4, 2) + 0.1.
   par(mar = (c(2, 4, 2, 2) + 0.1))
   temp <- index & (x$total == (i - 1))
   plot(x$carapace.width[temp], res[temp], xlab = NA, main = NA,
        xlim = c(60, 150), ylim = c(-4, 4), pch = 21, bg = "blue", cex = 0.5)
   lines(par("usr")[c(1, 2)], c(0, 0), col = "red", lwd = 2)
   text(65, 2, paste("Missing legs =", i-1))
}

# Residuals trends for regenerated legs:
windows()
temp <- index & (apply(m==2, 1, sum) == 1)
plot(x$carapace.width[temp], res[temp],
     xlab = "Carapace width(mm)", ylab = "Standardized residuals",
     main = "Regenerated leg residual plot",
        xlim = c(60, 150), ylim = c(-4, 4), pch = 21, bg = "blue", cex = 1)
grid()
temp <- index & (apply(m==2, 1, sum) == 2)
points(x$carapace.width[temp], res[temp], pch = 21, bg = "green", cex = 1)
lines(par("usr")[c(1, 2)], c(0, 0), col = "red", lwd = 2)
legend("topleft", pch = 21, pt.bg = c("blue", "green"),
       legend = c("1 regenerated leg", "2 regenerated legs"))


windows()
index <- (x$total == 0)
plot(log(x$carapace.width[index]), log(x$weight[index]),
     pch = 21, bg = "red", cex = 0.3, xlim = c(4, 5),
     xlab = "ln(carapace width)", ylab = "ln(weight)")
grid()

log.cw <- seq(4, 5, len = 100)
log.w <- l3$par[1]*log.cw + log(l3$par[2])
lines(log.cw, log.w, lwd = 2, col = "red")

index <- (x$total == 1)
points(log(x$carapace.width[index]), log(x$weight[index]), pch = 21, bg = "green", cex = 0.5)
for (i in 3:7){
   log.w <- l3$par[1]*log.cw + log(l3$par[2] + l3$par[i])
   lines(log.cw, log.w, lwd = 1, col = "green")
}

index <- (x$total == 2)
points(log(x$carapace.width[index]), log(x$weight[index]), pch = 21, bg = "blue", cex = 0.5)
index <- (x$total == 3)
points(log(x$carapace.width[index]), log(x$weight[index]), pch = 21, bg = "magenta", cex = 0.7)

legend("topleft", pt.bg = c("red", "green", "blue", "magenta"),
       pch = 21, legend = paste(0:3, "missing legs"))
       



windows()
layout(matrix(1:3, ncol = 1))
leg <- 1
for (i in 3:5){
   temp <- aggregate(data.frame(n = apply((m == 1)[x$shell.condition == i,c(leg, leg+5)], 1, mean)),
                     by = x[x$shell.condition == i, "len", drop = FALSE], mean)
   rownames(temp) <- temp[,1]
   temp <- temp[,2, drop = FALSE]
   dbarplot(temp, xlim = c(80, 145), ylim = c(0, 0.12))
   box()
}

# Plot carapace width versus pereopod weight:
windows()
beta <- l3$par[1]
alpha <- l3$par[3:7]
w <- seq(0, 150, len = 500)
plot(c(0, 150), c(0, 150), type = "n", xlab = "Leg weight(g)", ylab = "Carapace width(mm)")
abline(h = seq(0, 150, by = 10), v = seq(0, 150, by = 10), col = "lightgray", lty = "dotted",)
colours <- c("blue", "green", "purple", "red", "magenta")
for (i in 1:5){
   temp <- exp((log(w) - log(-alpha[i])) / beta)
   lines(w, temp, col = colours[i], lwd = 2)
}
legend("bottomright", legend = paste("Leg", 1:5), lwd = 2, col = colours)
box()

library(mgcv)
model <- gam(m1 ~ s(log.length), data = x, family = binomial, weights = rep(2, dim(x)[1]))

