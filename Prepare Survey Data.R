
x <- read.gulf(year = 1987:2012, survey = "sc", card = "bio")
temp <- x
xx <- x

# Define vector of colours to used:
cols <- c("blue", "green", "orange", "red", "purple")

rm(data, year, files, i, index, m)

# save(list = ls(), file = "Survey Missing Leg Data.Rdata")
temp <- rbind(x, data)

m <- matrix(0, ncol = 10, nrow = dim(x)[1])
for (i in 1:10){
   m[, i] <- as.numeric(substr(x$missing.leg, i, i))
}
m[is.na(m)] <- 0

aggregate(m[,1]+m[,6], by = x[c("year", "shell.condition")], mean)

# Check for parity:
temp <- aggregate(m==1, by = x[c("year")], mean)

windows()
plot(c(2000, 2012), c(0.01, 0.065),
     xlab = "Year", ylab = "Mean number of missing legs", type = "n")
grid()
for(i in 1:5){
   lines(temp[, 1], temp[, i+1], lwd = 2, col = cols[i], lty = "solid")
   lines(temp[, 1], temp[, i+6], lwd = 2, col = cols[i], lty = "dashed")
}
legend("topright", lwd = 2, col = rep(cols, rep(2, 5)), lty = c("solid", "dashed"), legend = paste(paste("Leg", rep(1:5, rep(2, 5))), c("left", "right")))


index <- (x$month %in% 7:9) & (x$sex %in% 1) & (x$carapace.width < 80)
index <- (x$month %in% 7:9) & (x$sex %in% 2)
temp <- aggregate(m[index, ]==1, by = x[index, c("year"), drop = FALSE], mean)
windows()
plot(c(2000, 2012), c(min(temp[, 2:11]), max(temp[, 2:11])),
     xlab = "Year", ylab = "Mean number of missing legs", type = "n")
grid()
for(i in 1:5){
   lines(temp[, 1], temp[, i+1], lwd = 2, col = cols[i], lty = "solid")
   lines(temp[, 1], temp[, i+6], lwd = 2, col = cols[i], lty = "dashed")
}
legend("topright", lwd = 2, col = rep(cols, rep(2, 5)), lty = c("solid", "dashed"), legend = paste(paste("Leg", rep(1:5, rep(2, 5))), c("left", "right")))


