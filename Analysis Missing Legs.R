library(gulf)

# Define data paths:
files <- "W:/Crab/Offshore Crab Common/Fishing Year 2010/Miscellaneous Research Data/Missing Legs Project/2010 Missing Legs.txt"
files <- c(files, "W:/Crab/Offshore Crab Common/Fishing Year 2011/Miscellaneous Research Data/Missing Legs Project/2011 Missing Legs.txt")

# Load lines as text:
data <- NULL
for (i in 1:2) data <- c(data, readLines(files[i]))

# Parse data:
x <- data.frame(day = as.numeric(substr(data, 2, 3)))
x$month <- as.numeric(substr(data, 4, 5))
x$year <- as.numeric(substr(data, 6, 9))
x$zone <- substr(data, 11, 12)
x$area <- substr(data, 14, 33)
x$vessel <- substr(data, 35, 54)
x$observer <- substr(data, 56, 95)
x$crab.number <- as.numeric(substr(data, 97, 100))
x$sex <- as.numeric(substr(data, 102, 102))
x$carapace.width <- as.numeric(substr(data, 104, 110))
x$chela.height <- as.numeric(substr(data, 112, 117))
x$shell.condition <- as.numeric(substr(data, 119, 119))
x$shell.condition.mossy <- substr(data, 120, 120)
x$durometer <- as.numeric(substr(data, 122, 124))
x$weight <- as.numeric(substr(data, 126, 133))
x$comment <- substr(data, 135, 164)
x$missing.legs <- substr(data, 166, 184)
x$missing.legs <- gsub(" ", "", x$missing.legs)
x$log.length <- log(x$carapace.width)
x$log.weight <- log(x$weight)

# Statistics on missing legs codes:
table(unlist(strsplit(substr(data, 166, 184), " ")))

m <- matrix(0, ncol = 10, nrow = dim(x)[1])
for (i in 1:10){
   m[, i] <- as.numeric(substr(x$missing.leg, i, i))
}
m[is.na(m)] <- 0
index <- (m == 3) | (m == 4) | (m == 5) | (m == 6) | (m == 8)
index <- apply(index, 1, sum) == 0
x <- x[index, ]
m <- m[index, ]

#=========================================================================================================
x$total <- apply(m, 1, sum)
x$len <- round(x$carapace.width / 5) * 5
x$m1 <- apply(m[,c(1,6)] == 1, 1, sum)
x$m2 <- apply(m[,c(2,7)] == 1, 1, sum)
x$m3 <- apply(m[,c(3,8)] == 1, 1, sum)
x$m4 <- apply(m[,c(4,9)] == 1, 1, sum)
x$m5 <- apply(m[,c(5,10)] == 1, 1, sum)
x$r1 <- apply(m[,c(1,6)] == 2, 1, sum)
x$r2 <- apply(m[,c(2,7)] == 2, 1, sum)
x$r3 <- apply(m[,c(3,8)] == 2, 1, sum)
x$r4 <- apply(m[,c(4,9)] == 2, 1, sum)
x$r5 <- apply(m[,c(5,10)] == 2, 1, sum)
x$d1 <- apply(m[,c(1,6)] == 7, 1, sum)
x$d2 <- apply(m[,c(2,7)] == 7, 1, sum)
x$d3 <- apply(m[,c(3,8)] == 7, 1, sum)
x$d4 <- apply(m[,c(4,9)] == 7, 1, sum)
x$d5 <- apply(m[,c(5,10)] == 7, 1, sum)

#=========================================================================================================
m0   <- lm(log.weight ~ log.length, data = x)
m1.1 <- lm(log.weight ~ log.length + m1, data = x)
m1.2 <- lm(log.weight ~ log.length + m2, data = x)
m1.3 <- lm(log.weight ~ log.length + m3, data = x)
m1.4 <- lm(log.weight ~ log.length + m4, data = x)
m1.5 <- lm(log.weight ~ log.length + m5, data = x)
m2   <- lm(log.weight ~ log.length + m1 + m2 + m3 + m4 + m5, data = x)
m2.1 <- lm(log.weight ~ log.length + m1 + m2 + I(m3 + m4) + m5, data = x)
m2.2 <- lm(log.weight ~ log.length + m1 + I(m2 + m3 + m4) + m5, data = x)
m2.3 <- lm(log.weight ~ log.length + m1 + I(m2 + m3 + m4), data = x)
m3.2 <- lm(log.weight ~ log.length + m1 + I(m2 + m3 + m4) + m5 + r1 + I(r2 + r3 + r4) + r5, data = x)
m4.2 <- lm(log.weight ~ log.length + m1 + I(m2 + m3 + m4) + m5 + r1 + I(r2 + r3 + r4) + r5 + d1 + I(d2 + d3 + d4) + d5 , data = x)

m3   <- lm(log.weight ~ log.length + total, data = x)

#=========================================================================================================
index <- (x$total == 0)
plot(log(x$carapace.width[index]), log(x$weight[index]), pch = 21, bg = "red", cex = 0.3, xlim = c(4, 5))
index <- (x$total == 1)
points(log(x$carapace.width[index]), log(x$weight[index]), pch = 21, bg = "green", cex = 0.3)
index <- (x$total == 2)
points(log(x$carapace.width[index]), log(x$weight[index]), pch = 21, bg = "blue", cex = 0.5)
index <- (x$total == 3)
points(log(x$carapace.width[index]), log(x$weight[index]), pch = 21, bg = "magenta", cex = 0.7)

