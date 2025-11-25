library(gulf.data)
library(gulf.graphics)

b <- read.scsbio(1988:2025, survey = "regular")

b <- b[which((b$sex == 1) & (b$carapace.width > 0) & (b$carapace.width < 150) & !is.na(b$carapace.width)), ]
b$missing.legs[b$missing.legs == "*********"] <- "**********"
b$year <- year(b)
b$missing.legs <- gsub("[2-7]", "*", b$missing.legs)
b$maturity <- is.mature(b)

# Missing leg matrix:
m <- matrix(NA, nrow = nrow(b), ncol = 10)
for (i in 1:10) m[,i] <- as.numeric(substr(b$missing.legs,i,i) == "1")

ix <- which(m[,1] == 1 & m[,6] == 1 & !is.na(b$chela.height))
table(b$year[ix])
ix <- which(m[,2] == 1 & m[,7] == 1)
table(b$year[ix])
ix <- which(m[,3] == 1 & m[,8] == 1)
table(b$year[ix])
ix <- which(m[,4] == 1 & m[,9] == 1)
table(b$year[ix])
ix <- which(m[,5] == 1 & m[,10] == 1)
table(b$year[ix])

ix <- which(apply(m[,1:5] ==  m[,6:10], 1, all) & b$missing.legs != "**********")
gbarplot(table(b$year[ix]))

# Define missing leg patterns:
m <- matrix(0, ncol = 10, nrow = nrow(b))
for (i in 1:10) m[, i] <- as.numeric(substr(b$missing.legs, i, i))
m[is.na(m)] <- 0

# Proportions of crab having no missing legs by size and year:
t <- table(b$year, round(b$carapace.width), apply(m == 1, 1, sum) == 0)
t <- t[,,"TRUE"] / (t[,,"FALSE"] + t[,,"TRUE"])
logit <- function(x) return(log(x / (1-x)))
image(as.numeric(dimnames(t)[[1]]), as.numeric(dimnames(t)[[2]]),
      logit(t), xlab = "", ylab = "", zlim = c(-3, 7), ylim = c(30, 120), col = terrain.colors(100))

# Average number of missing legs by specific group:
ix <- which(b$maturity & b$carapace.width > 60 & b$carapace.width < 70)
r <- aggregate(apply(m[ix, ], 1, sum), by = b[ix, "year", drop = FALSE], function(x) sum(x)/length(x))
gbarplot(r[,2], r[,1], ylim = c(0, 1.25))

r <- matrix(NA, nrow = 5, ncol = 5)
for (i in 1:5){
   for (j in 1:5){
      r[i,j] <- sum((m[,i] == 1) & (m[,j+5] == 1))
   }
}

p <- apply(m, 2, sum) / nrow(m)
nrow(m) * t(t(p[1:5])) %*% t(p[6:10])

# Generate all missing leg patterns:
a <- substr(as.character(intToBits(0:31)), 2,2)
dim(a) <- c(32, 32)
a <- t(a)[, 1:5]
a <- apply(a, 1, function(x) paste(x, collapse = ""))
a <- gsub("0", "*", a)


ll <- t(kronecker(matrix(1:2), matrix(1, ncol = 5, nrow = 5)))
ll <- rbind(0, cbind(0, 0, ll, 0), 0)
layout(ll)

bb <- b[which(b$maturity & b$carapace.width >= 60 & b$carapace.width < 95), ]

# Tabulate frequency of each missing leg pattern by year:
years <- sort(unique(bb$year))
r <- matrix(NA, nrow = length(a), ncol = length(years))
dimnames(r) <- list(pattern = a, year = years)
l <- r
for (i in 1:length(years)){
   ix <- which(bb$year == years[i])
   t <- substr(bb$missing.legs[ix], 1, 5)
   t <- table(t)
   l[names(t),i] <- as.numeric(t)
   t <- substr(bb$missing.legs[ix], 6, 10)
   t <- table(t)
   r[names(t),i] <- as.numeric(t)
}
l <- l / repvec(apply(l, 2, sum, na.rm = TRUE), nrow = nrow(r))
r <- r / repvec(apply(r, 2, sum, na.rm = TRUE), nrow = nrow(r))

layout(ll)
par(mar = c(0,0,0,0))

l <- l[-1, ]
ix <- order(nchar(gsub("[*]", "", a)))
ix <- rev(order(apply(l, 1, mean, na.rm = TRUE)))

image(years,  1:32, t(log(l[ix,])), ylab = "", yaxt = "n", xlab = "", zlim = c(-8, 0), col = terrain.colors(100))
axis(2, at = 0.5 + (1:31), labels = rownames(l)[ix], las = 2)
mtext("Pattern", 2, 4.0, cex = 1.5)
mtext("Year", 1, 3.0, cex = 1.5)

r <- r[-1, ]
ix <- order(nchar(gsub("[*]", "", a)))
ix <- rev(order(apply(r, 1, mean, na.rm = TRUE)))

image(years,  1:32, t(log(r[ix,])), ylab = "", yaxt = "n", xlab = "", zlim = c(-8, 0), col = terrain.colors(100))
mtext("Year", 1, 3.0, cex = 1.5)

# Number of missing legs by shell condition:
boxplot(apply(m, 1, mean) ~ b$shell.condition + b$year)
r <- aggregate(apply(m, 1, sum), by = b[c("year", "shell.condition")], mean)
plot(range(years), c(0, 1.8), type = "n", yaxs = "i")
grid()
for (i in 1:5){
   ix <- r$shell.condition == i
   lines(r$year[ix], r$x[ix], col = rainbow(5)[i], lwd = 2)
}


b$cw <- round(b$carapace.width)
ix <- which(b$maturity)
r <- aggregate(m[ix,], by = b[ix, c("cw"), drop = FALSE], mean)
plot(c(40, 110), c(0, 0.1), type = "n", yaxs = "i")
grid()
for (i in 1:5){
   lines(r$cw, r[,i+1], col = rainbow(5)[i], lwd = 2)
   lines(r$cw, r[,i+6], col = rainbow(5)[i], lwd = 2, lty = "dashed")
}


b$cw <- round(b$carapace.width)
ix <- which(!b$maturity)
r <- aggregate(m[ix,], by = b[ix, c("year", "cw"), drop = FALSE], mean)
plot(c(40, 110), c(0, 0.4), type = "n", yaxs = "i")
grid()
res <- matrix(NA, nrow = 150, ncol = length(years))
dimnames(res) <- list(cw = 1:150, year = years)
for (i in 1:length(years)){
   ix <- r$year == years[i]
   xx <- r$cw[ix]
   yy <- r[ix,3] + r[ix,3+5]
   tmp <- aggregate(yy, list(5*round(xx/5)), mean)

   lines(tmp[,1], tmp[,2], col = rainbow(length(years))[i], lwd = 2)

   res[as.character(r$cw[ix]), i] <- apply(r[ix,3:12], 1, sum)
}
image(years, as.numeric(rownames(res)),t(log(res)))

