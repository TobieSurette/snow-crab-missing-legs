library(gulf)

setwd("C:/Users/SuretteTJ/Desktop")

# Load missing legs data set:
b <- read.scbio(year = 1989:2017)
b <- b[which(b$sex == 2), ]
b$maturity <- is.mature(b, probability = TRUE)
b$cw <- round(b$carapace.width)
b <- b[!is.na(b$carapace.width), ]

# Fix weird missing cheliped errors:
index <- which((substr(b$missing.legs,1,1) %in% c("1", "2")) & (substr(b$missing.legs,6,6) %in% c("1", "2")) & !is.na(b$chela.height))
ii <- c(1, 6)[(runif(length(index)) > 0.5) + 1]
temp <- b$missing.legs[index]
for (i in 1:length(ii)) substr(temp, ii[i], ii[i]) <- "*"
b$missing.legs[index] <- temp

# Missing leg pattern inmatrix form:
M <- NULL
for (i in 1:10) M <- cbind(M, substr(b$missing.legs, i, i) == "1")
M <- M + 1 - 1

# Missing leg pattern inmatrix form:
R <- NULL
for (i in 1:10) R <- cbind(R, substr(b$missing.legs, i, i) == "2")
R <- R + 1 - 1

years <- sort(unique(b$year))


# ========== Variation along carapace widths of missing chelipeds ==============
k <- array(0, dim = c(length(years), 140, 5))
dimnames(k) <- list(year = years, width = 1:140, missing = 1:5)
n <- k
for (i in 1:5){
   for (j in 1:length(years)){
      index <- (b$year == years[j]) & (b$maturity > 0.5) & is.new(b)
      res <- aggregate(list(k = (M[index, i] == 1) | (M[index, i+5] == 1)), by = b[index, "cw", drop = FALSE], sum)
      res$n <- aggregate(list(n = (M[index, i] == 1) | (M[index, i+5] == 1)), by = b[index, "cw", drop = FALSE], length)$n
      res <- res[res$cw <= 140, ]
      tmp <- colorRamp(c("black", cols[i]))(0.75)
      tmp <- rgb(tmp[1,1]/255, tmp[1,2]/255, tmp[1,3]/255)
      k[j, res$cw, i] <- res$k
      n[j, res$cw, i] <- res$n        
   }
}

#windows(width = 11, height = 8.5)
jpeg(filename = "Missing Legs - Leg versus width - females.jpeg",
     width = 12*480, height = 11*480, quality = 75, res = 600) 
plot(c(40, 80), c(0, 0.25), type = "n", xlab = "Carapace width (mm)", ylab = "# Missing / crab", xaxs = "i", yaxs = "i", cex.lab = 1.5, xaxt = "n", cex.axis = 1.25)
axis(1, at = seq(40, 140, by = 5), cex.axis = 1.25) 
grid()
cols <- rainbow(5)
for (i in 1:5){
   ratio <- apply(k[as.character(2008:2017), , i] / n[as.character(2008:2017), , i], 2, mean, na.rm = TRUE)
   r <- NULL
   for (j in 2:(length(ratio)-1)) r[j-1] <- mean(ratio[(j-1):(j+1)], na.rm = TRUE)
   names(r) <- names(ratio)[2:(length(ratio)-1)]
   lines(as.numeric(names(r)), r, lwd = 3, col = cols[i])
}
box()
legend("topleft", legend = c("Cheliped", paste("Walking Leg", 1:4)), lwd = 3, col = cols, bg = "white", cex = 1.25)

# Immature
k <- array(0, dim = c(length(years), 140, 5))
dimnames(k) <- list(year = years, width = 1:140, missing = 1:5)
n <- k
for (i in 1:5){
   for (j in 1:length(years)){
      index <- (b$year == years[j]) & (b$maturity < 0.5) & is.new(b)
      res <- aggregate(list(k = (M[index, i] == 1) | (M[index, i+5] == 1)), by = b[index, "cw", drop = FALSE], sum)
      res$n <- aggregate(list(n = (M[index, i] == 1) | (M[index, i+5] == 1)), by = b[index, "cw", drop = FALSE], length)$n
      res <- res[res$cw <= 140, ]
      tmp <- colorRamp(c("black", cols[i]))(0.75)
      tmp <- rgb(tmp[1,1]/255, tmp[1,2]/255, tmp[1,3]/255)
      k[j, res$cw, i] <- res$k
      n[j, res$cw, i] <- res$n        
   }
}

for (i in 1:5){
   ratio <- apply(k[as.character(2008:2017), , i] / n[as.character(2008:2017), , i], 2, mean, na.rm = TRUE)
   r <- NULL
   for (j in 2:(length(ratio)-1)) r[j-1] <- mean(ratio[(j-1):(j+1)], na.rm = TRUE)
   names(r) <- names(ratio)[2:(length(ratio)-1)]
   r <- r[as.character(40:65)]
   lines(as.numeric(names(r)), r, lwd = 3, col = cols[i], lty = "dashed")
}
box()
dev.off()


# Annual rates of missing chelipeds
windows(width = 11, height = 8.5)
plot(range(years), c(0.00, 0.21), type = "n", xlab = "Year", ylab = "# Missing / crab", yaxs = "i", cex.lab = 1.5, xaxt = "n", cex.axis = 1.25)
cols <- rainbow(5)  
grid()

for (i in 1:5){
   index <- (b$maturity >= 0.5) #& (b$carapace.width >= 95) 
   res <- aggregate(list(k = (M[index, i] == 1) | (M[index, i+5] == 1)), by = b[index, "year", drop = FALSE], sum)
   res$n <- aggregate(list(n = (M[index, i] == 1) | (M[index, i+5] == 1)), by = b[index, "year", drop = FALSE], length)$n
   tmp <- colorRamp(c("black", cols[i]))(0.75)
   tmp <- rgb(tmp[1,1]/255, tmp[1,2]/255, tmp[1,3]/255)

   # Binomial confidence intervals:
   p <- res$k / res$n
   delta <- qnorm(0.975) * sqrt((1/1000)*p*(1-p))
   ci <- cbind(p - delta, p + delta)
   for (j in 1:nrow(res)){
      #lines(c(res$year[j], res$year[j]), ci[j,], lwd = 2, col = cols[i])
   }
   
   lines(res$year, res$k / res$n, col = tmp, lwd = 4)
   lines(res$year, res$k / res$n, col = cols[i], lwd = 2)     
}
axis(1, at = seq(1990, 2017, by = 2), cex.axis = 1.25)
legend("topright", legend = c("Cheliped", paste("Leg", 2:5)), lwd = 3, col = cols, bg = "white", cex = 1.25)
box()

# Variation along carapace widths of missing chelipeds
k <- array(0, dim = c(length(years), 140, 5))
dimnames(k) <- list(year = years, width = 1:140, missing = 1:5)
n = k
cols <- rainbow(5)  
for (i in 1:5){
   for (j in 1:length(years)){
      index <- (b$year == years[j]) & (b$maturity > 0.5) & is.new(b)
      res <- aggregate(list(k = (M[index, i] == 1) | (M[index, i+5] == 1)), by = b[index, "cw", drop = FALSE], sum)
      res$n <- aggregate(list(n = (M[index, i] == 1) | (M[index, i+5] == 1)), by = b[index, "cw", drop = FALSE], length)$n
      res <- res[res$cw <= 140, ]
      tmp <- colorRamp(c("black", cols[i]))(0.75)
      tmp <- rgb(tmp[1,1]/255, tmp[1,2]/255, tmp[1,3]/255)
      k[j, res$cw, i] <- res$k
      n[j, res$cw, i] <- res$n        
   }
}

windows(width = 11, height = 8.5)
plot(c(40, 130), c(0, 0.30), type = "n", xlab = "Carapace width (mm)", ylab = "# Missing / crab", xaxs = "i", yaxs = "i", cex.lab = 1.5, xaxt = "n", cex.axis = 1.25)
cols <- rainbow(5)  
grid()

#I <- apply(k / n, c(1, 2), sum, na.rm = TRUE)
#I <- k[,,2] / n[,,2]
#I <- I[, 60:140]

I <- (k[,,1] / n[,,1])[, 40:90]

# Image version:
windows(width = 11, height = 8.5)
image(40:90, years, t(I), col = colorRampPalette(c("white", "black"))(100), 
      breaks = c(seq(0, 0.50, len = 100), 1.1),
      xlab = "Carapace width (mm)", ylab = "Years", cex.lab = 1.25,
      xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
      xlim = c(39.5, 90.5), ylim = c(1989.5, 2017.5))
axis(2, at = seq(1990, 2018, by = 2))
axis(1, at = seq(40, 90, by = 5))
box()
lens <- 40:131
for (j in 1:length(years)) lines(par("usr")[1:2], rep(years[j],2)+0.5, lwd = 0.5, col = "grey60")
for (i in 1:length(lens)) lines(rep(lens[i],2)+0.5, par("usr")[3:4], lwd = 0.5, col = "grey60")
box()


# Bubble plot version:
plot(c(59.5, 130.5), c(1989.5, 2017.5), type = "n", xlab = "Carapace width (mm)", ylab = "Years", cex.lab = 1.25, xaxs = "i", yaxs = "i")
for (i in 1:nrow(I)){
   for (j in 1:ncol(I)){
      points(as.numeric(colnames(I)[j]), as.numeric(rownames(I)[i]), pch = 21, bg = "grey", cex = 2.2 * sqrt(I[i,j]))
   }
}

# Single leg correlation:
years <- sort(unique(b$year))
res <- NULL
ratio <- matrix(NA, nrow = length(years), ncol = 5)
for (i in 1:length(years)){
   index <- b$year == years[i]
   for (j in 1:5){
      p <- c(sum((M[index,j] == 0) & (M[index,j+5] == 0)), sum(xor(M[index,j] == 1, M[index,j+5] == 1)), sum((M[index,j] == 1) & (M[index,j+5] == 1)))
      p <- c(years[i], p, round(((p[2] / sum(p)) ^ 2) * sum(p), 1))
      ratio[i,j] <- round(p[length(p)-1] / p[length(p)], 2)
   } 
}
rownames(ratio) <- years
colnames(ratio) <- 1:5
ratio <- as.data.frame(ratio)

windows()
plot(range(years), c(0, 5), type = "n", xlab = "Year", ylab = "Ratio", cex.lab = 1.25, yaxs = "i")
dbarplot(ratio[, 1], years, add = TRUE)
abline(1, 0, lwd = 2, col = "red")
cols <- rainbow(4)
for (i in 2:5){
   lines(years, ratio[, i], lwd = 2, col = cols[i-1])
}
legend("topright", legend = 2:5, col = cols, lwd = 2)

# Spatial plots:
mif <- read.mif("U:/Snow Crab/Trawl Survey/2017/Generate Survey Stations 2017/grids2017.mif", mid = FALSE)
for (i in 1:length(mif)){
   tmp <- km2deg(mif[[i]]$x, mif[[i]]$y)
   mif[[i]]$longitude <- tmp$longitude
   mif[[i]]$latitude <- tmp$latitude
}
n <- length(mif) # Number of stations, which defines the sampling grid scheme.]

s <- read.scset(year = 1989:2017)
s <- s[s$season == "fall", ]
s <- s[s$valid == 1, ]
index <- match(b[c("year", "tow.id")], s[c("year", "tow.id")])
b <- b[-which(is.na(index)), ]

# Missing leg pattern inmatrix form:
M <- NULL
for (i in 1:10) M <- cbind(M, substr(b$missing.legs, i, i) == "1")
M <- M + 1 - 1

s$grid <- NA
for (i in 1:length(mif)){
   p <- as.polygon(mif[[i]]$longitude, mif[[i]]$latitude)
   s$grid[in.polygon(p, longitude(s), latitude(s))] <- i
}
index <- match(b[c("year", "tow.id")], s[c("year", "tow.id")])
b$grid <- s$grid[index]


index <- (b$year > 2007) & (b$year <= 2017) & (b$maturity > 0.5) & (b$carapace.width >= 60)
res <- aggregate(list(k = apply(M[index, ], 1, sum)), by = b[index, "grid", drop = FALSE], sum)
res$n <- aggregate(list(n = M[index, 1]), by = b[index, "grid", drop = FALSE], function(x) return(length(x)))$n
res$rate <- res$k / res$n

upper <- 1.5
windows(width = 11, height = 10)
#jpeg(filename = "Missing Legs - Mature Males 2013-2017.jpeg",
#     width = 12*480, height = 11*480, quality = 75, res = 600) 
     
gulf.map(land = FALSE, sea = TRUE)
for (i in 1:length(mif)){
   index <- which(res$grid == i)
   if (length(index) == 1 ){
      if (res$n[index] > 10){
         col <- colorRamp(c("white", "black"))(min(c(res$rate[index] / upper, 1)))
         col <- rgb(col[, 1] / 255, col[, 2] / 255 , col[, 3] / 255)
         polygon(mif[[i]]$longitude, mif[[i]]$latitude, col = col)
      }
   }
}
coastline(col = "grey80", border = "lightsteelblue4", lwd = 1.5)
box()
t <- seq(0, upper, by = 0.3)
cols <- (colorRamp(c("white", "black"))(t/upper)) / 255
cols <- rgb(cols[, 1], cols[, 2], cols[, 3])  
legend("bottomleft", legend = c(t[1:(length(t)-1)], paste0(t[length(t)], "+")), pch = 22, pt.bg = cols, pt.cex = 4, cex = 1.5, bg = "grey95")
#bathymetry(dem = FALSE)
wind.rose()
#dev.off()

# Cross leg correlation:
years <- sort(unique(b$year))
res <- NULL
ratio <- matrix(NA, nrow = length(years), ncol = 10)
y <- c(1, 2)
for (i in 1:length(years)){
   index <- b$year == years[i]
   for (j in 2:10){
      p <- c(sum((M[index,1] == 0) & (M[index,j] == 0)), sum((M[index,1] == 1) & (M[index,j] == 0)), sum((M[index,1] == 0) & (M[index,j] == 1)), sum((M[index,1] == 1) & (M[index,j] == 1)))
      ratio[i,j] <- round(p[4] / (sum(p) * (p[2] / sum(p)) * (p[3] / sum(p))), 2)
   } 
}
rownames(ratio) <- years
colnames(ratio) <- c(paste0("L", 1:5), paste0("R", 1:5))
ratio <- as.data.frame(ratio)

windows()
plot(range(years), c(0, 8), type = "n", xlab = "Year", ylab = "Ratio", cex.lab = 1.25, yaxs = "i")
dbarplot(ratio[, 1], years, add = TRUE)
abline(1, 0, lwd = 2, col = "red")
cols <- rainbow(ncol(ratio)-1)
for (i in 2:10) lines(years, ratio[, i], lwd = 2, col = cols[i-1])
legend("topright", legend = 2:5, col = cols, lwd = 2)


# Contigency table analysis:
p <- c(sum((M[,2] == 0) & (M[,7] == 0)), sum(xor(M[,2] == 1, M[,7] == 1)), sum((M[,2] == 1) & (M[,7] == 1)))
((p[2] / sum(p)) ^ 2) * sum(p)


# Left side:
t <- t(M[,1:5]) %*% M[,1:5]
p <- apply(M[,1:5],2, sum) / nrow(M)
P <- t(t(p )) %*% t(p)
diag(P) <- p
round(t / round(P * nrow(M)), 1)

# Right side:
t <- t(M[,6:10]) %*% M[,6:10]
p <- apply(M[,6:10],2, sum) / nrow(M)
P <- t(t(p )) %*% t(p)
diag(P) <- p
round(t / round(P * nrow(M)), 1)

# Cross-side:
t <- t(M[,1:5]) %*% M[,6:10]
p <- apply(M,2,sum) / nrow(M)
P <- repvec(p[1:5], ncol = 5) * repvec(p[6:10], nrow = 5)
round(t / round(P * nrow(M)), 1)

# Cross-side:
t <- t(M[,6:10]) %*% M[,1:5]
p <- apply(M,2,sum) / nrow(M)
P <- repvec(p[6:10], ncol = 5) * repvec(p[1:5], nrow = 5)
round(t / round(P * nrow(M)), 1)




loglike <- function(theta, M){
   # Base probabilities:
   beta <- matrix(0, nrow = 10, ncol = 10)
  # diag(beta) <- c(theta[1:5], theta[1:5])
   mu <- repvec(c(theta[1:5], theta[1:5]), nrow = nrow(M))
   
   p <- 1 / (1 + exp(-mu))
   
   ll <- -sum(M * log(p) + (1-M) * log(1-p))
   
   return(ll) 
}

t(combn(10, 2))

loglike(1:5, M)
p <- optim(runif(5), loglike, M = M, control = list(trace = 3, maxit = 5000))


index <-  t(combn(10, 2))
T <- NULL
for (i in 1:10){
   for (j in 1:i){
      print(c(i,j))
      T <- cbind(T, M[,i] * M[,j])
   }
}



for (i in 1:10){
   for (j in 1:i){
      mu[i,] <- mu[,i] + beta[i,j] * M[,i] * M[,j]
   }
}
 
loglike <- function(theta, M){
   # Base probabilities:
   beta <- matrix(NA, nrow = 10, ncol = 10)

   # Diagonal elements:
   k <- 1
   for (i in 1:5){
      for (j in 1:i){
         beta[i,j] <- theta[k]     # Left side.
         beta[i+5,j+5] <- theta[k] # Right side.
         beta[i+5,j] <- theta[k] + theta[k+15]
         k <- k + 1
      }
   }
   
   # Symmetry mapping lower quadrant: 
   for (i in 2:5){
      for (j in 1:(i-1)){
         beta[j+5,i] <- beta[i+5,j] 
      }
   }

   # Symmetry mapping:
   for (i in 2:10){ 
      for (j in 1:(i-1)){ 
         beta[j,i] <- beta[i,j]
      }
   }
   
   # Calculate all two-way interaction contributions:
   mu <- matrix(0, nrow = nrow(M), ncol = 10)
   for (i in 1:10){
      for (j in 1:10){
         mu[,i] <- mu[,i] + beta[i,j] * M[,i] * M[,j]
      }
   }
   
   # Logistic transform:
   p <- 1 / (1 + exp(-mu))
   
   # Bernouilli log-likelihood:
   ll <- -sum(M * log(p) + (1-M) * log(1-p))
   
   return(ll)
}

loglike(1:5, M)
p <- optim(runif(30), loglike, M = M, control = list(trace = 3, maxit = 5000))


mapping <- rbind(c(1, 6, 7, 9, 12, 16, 21, 22, 24, 27)

loglike <- function(theta, M){
   # Base probabilities:
   p <- exp(theta[1:4]) / (1 + sum(exp(theta[1:4])))
   p <- c(1-sum(p), p)
   
   
   theta[1] <- 0
   
   alpha <- c(theta[c(1, 3, 6, 10, 15)], theta[c(1, 3, 6, 10, 15)])
   alpha <- repvec(alpha, nrow = nrow(M) )
   beta <- matrix(NA, nrow = 10, ncol = 10)

   # Diagonal elements:
   k <- 1
   for (i in 1:5){
      for (j in 1:i){
         beta[i,j] <- theta[k]
         beta[i+5,j] <- beta[i,j]
         beta[i+5,j+5] <- beta[i,j]
         beta[i+5,j] <- theta[k] + theta[k+15]
         k <- k + 1
      }
   }
   diag(beta) <- 0
   
   # Symmetry mapping low quadrant: 
   for (i in 2:5){
      for (j in 1:(i-1)){
         beta[j+5,i] <- beta[i+5,j] 
      }
   }

   # Symmetry mapping:
   for (i in 2:10){ 
      for (j in 1:(i-1)){ 
         beta[j,i] <- beta[i,j]
      }
   }
   
   #print(beta)
   mu <- (M %*% beta) + alpha
   mu <- alpha
   p <- 1 / (1 + exp(-mu))
   
   return(-sum(M * log(p))) 
}

# Marginal probabilities:
logit <- function(x) return(log(x) - log(1-x))
p <- 1 / (1 + exp(-theta[1:5])
p <- c(p, p)
R <- kronecker(t(p), t(t(p)))

t(M) %*% M

> t(M) %*% M

      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
 [1,]  500  107   75   62   48   96   73   51   45    45
 [2,]  107  809  167  120   90   86  118   75   72    67
 [3,]   75  167  727  141   83   68  103   70   63    53
 [4,]   62  120  141  641  117   46   74   68   82    52
 [5,]   48   90   83  117  558   47   72   58   67    71
 [6,]   96   86   68   46   47  609  130  101   77    64
 [7,]   73  118  103   74   72  130  882  202  123   114
 [8,]   51   75   70   68   58  101  202  725  123    85
 [9,]   45   72   63   82   67   77  123  123  686   135
[10,]   45   67   53   52   71   64  114   85  135   567

M[,1] * t(M[,2])


loglike <- function(theta, M){
   # Base probabilities:
   alpha <- theta[1:5]
   alpha <- c(alpha, alpha)
   mu <- NA * M
   
   beta <- matrix()
   
   k <- 6
   for (i in 1:10){
      mu[, i] <- alpha[i] + M[, setdiff(1:10, i)] %*% t(t(theta[setdiff(1:10, i)]))
   } 
   
   
   theta[1] <- 0
   
   alpha <- c(theta[c(1, 3, 6, 10, 15)], theta[c(1, 3, 6, 10, 15)])
   alpha <- repvec(alpha, nrow = nrow(M) )
   beta <- matrix(NA, nrow = 10, ncol = 10)

   # Diagonal elements:
   k <- 1
   for (i in 1:5){
      for (j in 1:i){
         beta[i,j] <- theta[k]
         beta[i+5,j] <- beta[i,j]
         beta[i+5,j+5] <- beta[i,j]
         beta[i+5,j] <- theta[k] + theta[k+15]
         k <- k + 1
      }
   }
   diag(beta) <- 0
   
   # Symmetry mapping low quadrant: 
   for (i in 2:5){
      for (j in 1:(i-1)){
         beta[j+5,i] <- beta[i+5,j] 
      }
   }

   # Symmetry mapping:
   for (i in 2:10){ 
      for (j in 1:(i-1)){ 
         beta[j,i] <- beta[i,j]
      }
   }
   
   #print(beta)
   mu <- (M %*% beta) + alpha
   mu <- alpha
   p <- 1 / (1 + exp(-mu))
   
   return(-sum(M * log(p))) 
}


