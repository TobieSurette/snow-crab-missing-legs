library(gulf.data)
library(gulf.graphics)
clm()
clg()

# Read snow crab survey biological data:
b <- read.scsbio(1988:2025, survey = "regular")
b$year <- year(b)
b <- b[which((b$carapace.width > 0) & (b$carapace.width <= 150)), ]
b$tow.id <- tow.id(b)

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

# Load tow data:
x <- read.scsset(years, survey = "regular")
x$year <- year(x)
import(x) <- aggregate(list(n = apply(m == 1, 1, sum, na.rm = TRUE)), by = b[key(x)], mean)

library(gulf.spatial)

map.new()
points(longitude(x), latitude(x), cex = 1.5*x$n)

# Female time series of relative missing leg frequencies:

windows()
temp <- (m[, 1:5] == 1) + (m[, 6:10] == 1)
sex = 1
p <- matrix(NA, nrow = length(years), ncol = 5)
for (i in 1:length(years)){
   ix <- (b$year == years[i]) & !is.mature(b) & (b$sex == sex) & (b$carapace.width <= 80)  & (!is.na(b$carapace.width))
   p[i, ] <- apply(temp[ix, ], 2, sum, na.rm = TRUE)
}
p <- p / repvec(apply(p, 1, sum), ncol = 5)
dimnames(p) <- list(years, paste("Leg", 1:5) )
plot(c(1987, 2012), c(0.05, 0.35), type = "n",
     xlab = "Year", ylab = "Relative Missing Leg Frequency")
title(main = ifelse(sex == 1, "Males", "Females"))
grid()
cols <- rainbow(5)
for (i in 1:5){
   print(i)
   lines(years, p[, i], col = cols[i], lwd = 2)
}
legend("topright", lwd = 2, col = cols, legend = paste("Leg", 1:5), bg = "white")

# Relative frequencies versus carapace width:
sex <- 2
year <- 1990:2012
maturity <- TRUE
b$bincw <- round(b$carapace.width)
if (maturity){
   temp <- b[(b$sex %in% sex) & (b$year %in% year) & is.mature(b), ]
}else{
   temp <- b[(b$sex %in% sex) & (b$year %in% years) & !is.mature(b), ]
}
# Prepare matrix of code values:
mm <- matrix(0, ncol = 10, nrow = dim(temp)[1])
for (i in 1:10){
   mm[, i] <- as.numeric(substr(temp$missing.leg, i, i))
}
mm[is.na(mm)] <- 0
t <- table(temp$bincw)
t <- as.numeric(names(t[t > 100]))
bins <- seq(min(t), max(t), by = unique(diff(t))[1])
windows()
f <- (mm[, 1:5] == 1) + (mm[, 6:10] == 1)
p <- matrix(NA, nrow = length(bins), ncol = 5)
for (i in 1:length(bins)){
   p[i, ] <- apply(f[(temp$bincw == bins[i]), ], 2, sum, na.rm = TRUE)
}
p <- p / repvec(apply(p, 1, sum), ncol = 5)
dimnames(p) <- list(bins, paste("Leg", 1:5) )
plot(c(min(bins), max(bins)), c(0.00, 0.35), type = "n",
     xlab = "Carapace width(mm)", ylab = "Relative Missing Leg Frequency")
str <- ifelse(sex == 1, "Males", "Females")
str <- paste(str, ifelse(maturity, "Mature", "Immature"))
title(main = str)
grid()
for (i in 1:5){
   print(i)
   lines(bins, p[, i], col = cols[i], lwd = 2)
}
legend("topleft", lwd = 2, col = cols, legend = paste("Leg", 1:5), bg = "white")

# Time series of relative missing leg frequencies:
temp <- (m[, 1:5] == 1) + (m[, 6:10] == 1)
sex = 1
p <- matrix(NA, nrow = length(year), ncol = 5)
for (i in 1:length(year)){
   index <- (data$year == year[i]) & (data$sex == sex) & (data$carapace.width <= 80)  & (!is.na(data$carapace.width))
   p[i, ] <- apply(temp[index, ], 2, sum)
}
p <- p / repvec(apply(p, 1, sum), ncol = 5)
dimnames(p) <- list(year, paste("Leg", 1:5) )
dbarplot(p, xlab = "Year")

windows()
plot(c(1987, 2012), c(0.05, 0.35), type = "n",
     xlab = "Year", ylab = "Relative Missing Leg Frequency")
title(main = ifelse(sex == 1, "Males", "Females"))
grid()
for (i in 1:5){
   print(i)
   lines(year, p[, i], col = cols[i], lwd = 2)
}
legend("topright", lwd = 2, col = cols, legend = paste("Leg", 1:5))


data$cwbin <- round(data$carapace.width/5)*5

aggregate(m[,1]+m[,6], by = data[c("year", "shell.condition")], mean)

# Check for parity:
#index <- (data$month %in% 7:9) & (data$sex %in% 1) & (data$carapace.width < 80)
index <- (data$month %in% 7:9) & (data$sex %in% 2)
                                                          temp <- aggregate(m[index, ] ==1, by = data[index, c("year"), drop = FALSE], mean)

windows()
plot(c(1987, 2012), c(0.01, 0.065),
     xlab = "Year", ylab = "Mean number of missing legs", type = "n")
grid()
for(i in 1:5){
   lines(temp[, 1], temp[, i+1], lwd = 2, col = cols[i], lty = "solid")
   lines(temp[, 1], temp[, i+6], lwd = 2, col = cols[i], lty = "dashed")
}
legend("topright", lwd = 2, col = rep(cols, rep(2, 5)), lty = c("solid", "dashed"), legend = paste(paste("Leg", rep(1:5, rep(2, 5))), c("left", "right")))

# Mean number of missing legs versus cw:
windows()
index <- (data$month %in% 7:9) & (data$sex %in% 1)
plot(c(40, 140), c(0, 0.075),
     ylab = "Mean number of missing legs",
     xlab = "Carapace width(mm)", type = "n")
grid()
for (i in 1:5){
   temp <- aggregate(data.frame(n = apply((m[index, ] == 1)[,c(i, i+5)], 1, mean)),
                     by = data[index, "cwbin", drop = FALSE], mean)
   lines(temp[, 1], temp[, 2], col = cols[i], lwd = 2)
}
legend("topright", legend = paste("Leg", 1:5), lwd = 2, col = cols)

# Mean number of missing legs versus cw:
windows()
index <- (data$month %in% 7:9) & (data$sex %in% 1)
plot(c(40, 140), c(0, 0.075),
     ylab = "Mean number of missing legs",
     xlab = "Carapace width(mm)", type = "n")
grid()
leg <- 1:5
for (j in 1:length(year)){
   index2 <- index & (data$year == year[j])
   temp <- aggregate(data.frame(n = apply((m[index2, ] == 1)[,c(leg, leg+5)], 1, mean)),
                     by = data[index2, "cwbin", drop = FALSE], mean)
   lines(temp[, 1], temp[, 2], col = blues9[(j %% 9)+1], lwd = (j/10)+1)
}
#legend("topright", legend = paste("Leg", 1:5), lwd = 2, col = cols)

# Mean number of missing legs versus cw:
windows()
index <- (data$month %in% 7:9) & (data$sex %in% 2)
plot(c(40, 80), c(0, 0.075),
     ylab = "Mean number of missing legs",
     xlab = "Carapace width(mm)", type = "n")
grid()
leg <- 1:5
for (j in 1:length(year)){
   index2 <- index & (data$year == year[j])
   temp <- aggregate(data.frame(n = apply((m[index2, ] == 1)[,c(leg, leg+5)], 1, mean)),
                     by = data[index2, "cwbin", drop = FALSE], mean)
   lines(temp[, 1], temp[, 2], col = blues9[(j %% 9)+1], lwd = (j/10)+1)
}
#legend("topright", legend = paste("Leg", 1:5), lwd = 2, col = cols)

# Bubble plot:
windows()
index <- (data$month %in% 7:9) & (data$sex %in% 1) & (data$carapace.width < 130) & (data$shell.condition == 1)
plot(c(1987, 2012), c(40, 130),
     xlab = "Year",
     ylab = "Carapace width(mm)",
     main = "Mean number of missing legs in males",
     type = "n")
grid()
leg <- 1:5
for (j in 1:length(year)){
   index2 <- index & (data$year == year[j])
   if (year[j] != 1996){
      temp <- aggregate(data.frame(n = apply((m[index2, ] == 1)[,c(leg, leg+5)], 1, mean)),
                        by = data[index2, "cwbin", drop = FALSE], mean)
      points(rep(year[j], dim(temp)[1]), temp[, 1], pch = 21, cex = 40*temp[, 2], bg = "red")
   }
}

# Bubble plot:
windows()
index <- (data$month %in% 7:9) & (data$sex %in% 2) & (data$carapace.width < 80)
plot(c(1987, 2012), c(40, 80),
     xlab = "Year",
     ylab = "Carapace width(mm)",
     main = "Mean number of missing legs in females",
     type = "n")
grid()
leg <- 1:5
for (j in 1:length(year)){
   index2 <- index & (data$year == year[j])
   if (year[j] != 1996){
      temp <- aggregate(data.frame(n = apply((m[index2, ] == 1)[,c(leg, leg+5)], 1, mean)),
                        by = data[index2, "cwbin", drop = FALSE], mean)
      points(rep(year[j], dim(temp)[1]), temp[, 1], pch = 21, cex = 40*temp[, 2], bg = "red")
   }
}
