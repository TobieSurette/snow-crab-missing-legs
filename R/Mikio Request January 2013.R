library(gulf)

rm(list = ls())
graphics.off()
setwd("U:/Snow Crab/Missing Legs")
load("Survey Missing Leg Data.Rdata")

data$maturity <- is.mature(data)

# Prepare matrix of code values:
m <- matrix(0, ncol = 10, nrow = dim(data)[1])
for (i in 1:10){
   m[, i] <- as.numeric(substr(data$missing.leg, i, i))
}
m[is.na(m)] <- 0
data$total <- apply(m == 1, 1, sum, na.rm = TRUE)


year <- sort(unique(data$year))
p <- rep(NA, length(year))
for (i in 1:length(year)){
   index <- (data$year == year[i]) & (data$carapace.width >= 95) & (data$carapace.width < 100) & !is.na(data$carapace.width)
   a <- sum(data$total[index & (data$sex == 1) & (data$maturity) & !is.na(data$maturity)] > 0)
   b <- sum(data$total[index ] > 0)
   p[i] <- a / b
}


