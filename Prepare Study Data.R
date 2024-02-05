library(gulf)

rm(list = ls())

# Define data paths:
files <- "W:/Crab/Offshore Crab Common/Fishing Year 2010/Miscellaneous Research Data/Missing Legs Project/2010 Missing Legs.txt"
files <- c(files, "W:/Crab/Offshore Crab Common/Fishing Year 2011/Miscellaneous Research Data/Missing Legs Project/2011 Missing Legs.txt")

# Load lines as text:
data <- NULL
for (i in 1:2){
   data <- c(data, readLines(files[i]))
}

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
x$durometer[x$durometer <= 0] <- NA
x$weight <- as.numeric(substr(data, 126, 133))
x$comment <- substr(data, 135, 164)
x$missing.legs <- substr(data, 166, 184)
x$missing.legs <- gsub(" ", "", x$missing.legs)
x$log.length <- log(x$carapace.width)
x$log.cw <- log(x$carapace.width)
x$log.weight <- log(x$weight)
x$abdomen.width <- NA
x$maturity <- is.mature.scbio(x)

# Remove NA weight data:
x <- x[!is.na(x$weight), ]

m <- matrix(0, ncol = 10, nrow = dim(x)[1])
for (i in 1:10){
   m[, i] <- as.numeric(substr(x$missing.leg, i, i))
}
m[is.na(m)] <- 0
index <- (m == 3) | (m == 4) | (m == 5) | (m == 6) | (m == 8)
#index <- (apply(index, 1, sum) == 0 ) & (!is.na(x$weight)) & !is.na(x$chela.height) &
#          !(rownames(x) %in% c("128","237","275","290","348","450","1444","1870","2367","2423","2705","2870","3078","3166","3240","3245"))
model <- lm(log.weight ~ log.cw, data = x)
index <- (apply(index, 1, sum) == 0 ) & (!is.na(x$weight)) & !is.na(x$chela.height) &
          (abs(residuals(model)) < 0.3)
x <- x[index, ]
m <- m[index, ]

model <- lm(log.weight ~ log.cw, data = x)

x$total <- apply(m == 1, 1, sum)
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

rm(files, i, index, model)

save(list = ls(), file = "Missing Leg Data.Rdata")
