rm(list = ls())
graphics.off()

# setwd("C:/Documents and Settings/SuretteTJ.ENT/Desktop/Snow Crab Missing Legs")
load("Missing Leg Data.Rdata")

# Define vector of colours to used:
cols <- c("blue", "green", "orange", "red", "purple")

# Carapace width distribution:
windows()
hist(x$carapace.width, n = 100,
     main = "", xlab = "Carapace width(mm)", ylab = "Frequency")

# Durometer distribution:
windows()
hist(x$durometer[x$durometer <= 100], n = 100, main = "", xlab = "Durometer", ylab = "Frequency")

# Plot morphometric measurements:
windows()
plot(c(60, 150), c(10, 45), type = "n", xlab = "Carapace width(mm)", ylab = "Chela height(mm)",)
grid()
points(x$carapace.width, x$chela.height, pch = 21, bg = cols[1], cex = 0.5)

# Frequency table of crab having specified number of missing legs:
windows()
dbarplot(table(x$total), xlab = "Number of missing legs", ylab = "Number of crab observed")

# Missing leg code frequency table:
temp <- table(unlist(strsplit(substr(data, 166, 184), " ")))
windows()
dbarplot(temp[setdiff(names(temp), "*")],
         main = "Number of missing legs by code",
         xlab = "Missing leg code", ylab = "Observed number")

# Length-weight relationship on regular scale:
windows()
plot(c(60, 150), c(0, 1300), type = "n", xlab = "Carapace width(mm)", ylab = "Weight(g)")
grid()
for (i in 0:3){
   index <- (x$total == i)
   points(x$carapace.width[index], x$weight[index], pch = 21, bg = cols[i+1], cex = 0.3+i/10)
}
legend("topleft", pt.bg = cols,
       pch = 21, legend = paste(0:3, "missing legs"))

# Length-weight relationship on log scale:
windows()
plot(c(4.2, 5), c(4.8, 7.3), type = "n", xlab = "ln(carapace width)", ylab = "ln(weight)")
grid()
for (i in 0:3){
   index <- (x$total == i)
   points(log(x$carapace.width[index]), log(x$weight[index]), pch = 21, bg = cols[i+1], cex = 0.3+i/10)
}
legend("topleft", pt.bg = cols,
       pch = 21, legend = paste(0:3, "missing legs"))

# Plot de-trended data, assuming a beta of 3:
windows()
f <- function(x, y) return(log(y) - 3 * log(x))
plot(c(60, 150), c(-8.1, -7.55), type = "n", xlab = "Carapace width(mm)", ylab = "Alpha")
grid()
for (i in 0:3){
   index <- x$total == i
   xx <- x$carapace.width[index]
   yy <- f(x$carapace.width[index], x$weight[index])
   points(xx, yy, pch = 21, bg = cols[i+1], cex = 0.5+i/10)
   lines(par("usr")[1:2], c(mean(yy), mean(yy)), col = cols[i+1], lwd = 2)
   # Regression line:
   # beta <- coef(lm(yy ~ xx))
   # abline(beta[1], beta[2], , col = cols[i+1], lwd = 2)
}
legend("topleft", pt.bg = cols, lwd = 2, col = cols,
       pch = 21, legend = paste(0:3, "missing legs"))

# Mean number of missing legs versus cw:
windows()
plot(c(80, 145), c(0, 0.16),
     ylab = "Mean number of missing legs",
     xlab = "Carapace width(mm)", type = "n")
grid()
for (i in 1:5){
   tmp <- aggregate(data.frame(n = apply((m == 1)[,c(i, i+5)], 1, mean)),
                     by = x["len"], mean)
   lines(temp[, 1], temp[, 2], col = cols[i], lwd = 2)
}
legend("topright", legend = paste("Leg", 1:5), lwd = 2, col = cols)

# Mean number of regenerated legs versus cw:
windows()
plot(c(80, 145), c(0, 0.05),
     ylab = "Mean number of regenerated legs",
     xlab = "Carapace width(mm)", type = "n")
grid()
for (i in 1:5){
   temp <- aggregate(data.frame(n = apply((m == 2)[,c(i, i+5)], 1, mean)),
                     by = x["len"], mean)
   lines(temp[, 1], temp[, 2], col = cols[i], lwd = 2)
}
legend("topright", legend = paste("Leg", 1:5), lwd = 2, col = cols)

