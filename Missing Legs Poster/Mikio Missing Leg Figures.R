library(gulf.data)
library(gulf.spatial)

# Load missing legs data set:
b <- read.scsbio(year = 1989:2022, survey = "regular")
b <- b[which(b$sex == 1), ]
b$maturity <- is.mature(b, probability = TRUE)
b$cw <- round(b$carapace.width)
b <- b[!is.na(b$carapace.width), ]
b$year <- year(b)

# Fix weird missing cheliped errors:
index <- which((substr(b$missing.legs,1,1) %in% c("1", "2")) & (substr(b$missing.legs,6,6) %in% c("1", "2")) & !is.na(b$chela.height))
ii <- c(1, 6)[(runif(length(index)) > 0.5) + 1]
temp <- b$missing.legs[index]
for (i in 1:length(ii)) substr(temp, ii[i], ii[i]) <- "*"
b$missing.legs[index] <- temp

# ============================= Spatial plots =================================
years <- sort(unique(b$year))
s <- read.scsset(year = years, survey = "regular", valid = 1)
s$year <- year(s)
index <- match(b[c("year", "tow.id")], s[c("year", "tow.id")])
b <- b[-which(is.na(index)), ]

# Missing leg pattern in matrix form:
M <- NULL
for (i in 1:10) M <- cbind(M, substr(b$missing.legs, i, i) == "1")
M <- M + 1 - 1

s$grid <- grid.scs(lon(s), lat(s))
index <- match(b[c("year", "tow.id")], s[c("year", "tow.id")])
b$grid <- s$grid[index]


index <- (b$year > 2012) & (b$year <= 2022) & (b$maturity > 0.5) & (b$carapace.width >= 60)
res <- aggregate(list(k = apply(M[index, ], 1, sum)), by = b[index, "grid", drop = FALSE], sum)
res$n <- aggregate(list(n = M[index, 1]), by = b[index, "grid", drop = FALSE], function(x) return(length(x)))$n
res$rate <- res$k / res$n

grids <- read.gulf.spatial("scs grids")
for (i in 1:length(grids)){
   tmp <- km2deg(grids[[i]]$x, grids[[i]]$y)
   grids[[i]]$longitude <- tmp[,1]
   grids[[i]]$latitude <- tmp[,2]
}

upper <- 1.5
#windows(width = 11, height = 10)
jpeg(filename = "Missing Legs - Mature Males 2013-2017.jpeg",
     width = 12*480, height = 11*480, quality = 75, res = 600)

map.new()
map("bathymetry")
for (i in 1:length(grids)){
   ix <- which(res$grid == i)
   if (length(ix) == 1 ){
      if (res$n[ix] > 10){
         col <- colorRamp(c("white", "black"))(min(c(res$rate[ix] / upper, 1)))
         col <- rgb(col[, 1] / 255, col[, 2] / 255 , col[, 3] / 255)
         polygon(grids[[i]]$longitude, grids[[i]]$latitude, col = col)
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
dev.off()

years <- sort(unique(b$year))

# Commercial remaining correlation with missing leg rate of intermediate-sized new matures:
index <- (b$maturity > 0.5) & !is.new.shell(b) & (b$carapace.width >= 95)
res <- aggregate(list(rate = apply(M[index,], 1, sum)), by = b[index, c("year", "tow.id")], mean)
res$n <- aggregate(list(n = apply(M[index,], 1, sum)), by = b[index, c("year", "tow.id")], length)$n
index <- (b$maturity > 0.5) & is.new.shell(b) & (b$carapace.width >= 60) & (b$carapace.width <= 95)
new <- aggregate(list(rate = apply(M[index,], 1, sum)), by = b[index, c("year", "tow.id")], mean)
new$n <- aggregate(list(n = apply(M[index,], 1, sum)), by = b[index, c("year", "tow.id")], length)$n
index <- match(new[c("year", "tow.id")], res[c("year", "tow.id")])
new$res <- res$n[index]
new$res[is.na(new$res)] <- 0
boxplot(log(new$rate) ~ new$res)

index <- ((b$maturity > 0.5) & (b$carapace.width >= 95)) + 1 - 1#& !is.new(b)
res <- aggregate(list(n = index), by = b[c("year"), drop = FALSE], mean)
index <- (b$maturity > 0.5) & is.new.shell(b) & (b$carapace.width >= 60) & (b$carapace.width <= 95)
new <- aggregate(list(rate = apply(M[index, c(1,6)], 1, sum)), by = b[index, c("year"), drop = FALSE], mean)
new$n <- aggregate(list(n = apply(M[index, c(1,6)], 1, sum)), by = b[index, c("year"), drop = FALSE], length)$n
index <- match(new[c("year")], res[c("year")])
new$res <- res$n[index]
new$res[is.na(new$res)] <- 0

jpeg(filename = "Missing Legs - Commercial versus Intermediate Rate.jpeg",
     width = 12*480, height = 12*480, quality = 75, res = 600)
windows()
plot(new$res, new$rate, pch = 21, cex = 1.5, xlim = c(0, 0.35), xaxs = "i", ylim = c(0, 0.22), yaxs = "i", bg = "grey", xlab = "Commercial density", ylab = "# Missing / Crab", cex.lab = 1.5)
text(new$res, new$rate, new$year, pos = 2)
m <- lm(new$rate ~ new$res)
abline(m, col = "red", lwd = 2)
dev.off()

lines(new$res, new$rate)


plot(new$res[1:(length(new$res)-1)], new$rate[2:length(new$res)], pch = 21, cex = 1.5, bg = "grey")
text(new$res, new$rate, new$year, pos = 2)
lines(new$res, new$rate)


# ========== Variation along carapace widths of missing chelipeds ==============
k <- array(0, dim = c(length(years), 140, 5))
dimnames(k) <- list(year = years, width = 1:140, missing = 1:5)
n <- k
for (i in 1:5){
   for (j in 1:length(years)){
      index <- (b$year == years[j]) & (b$maturity > 0.5) & is.new.shell(b)
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
#jpeg(filename = "Missing Legs - Leg versus width.jpeg", width = 24*480, height = 11*480, quality = 75, res = 600)

tiff(file = paste0("Missing Legs Poster/Missing legs by cw.tiff"), compression = "lzw", units = "in", res = 300, height = 7.5, width = 7.5)
plot(c(40, 130), c(0, 0.25), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", cex.lab = 1.0, xaxt = "n", cex.axis = 1.0)
axis(1, at = seq(40, 140, by = 10), cex.axis = 1.0)
axis(1, at = seq(45, 140, by = 10), cex.axis = 1.0)
grid()
cols <- rainbow(5)
for (i in 1:5){
   ratio <- apply(k[as.character(2008:2017), , i] / n[as.character(2008:2017), , i], 2, mean, na.rm = TRUE)
   r <- NULL
   for (j in 2:(length(ratio)-1)) r[j-1] <- mean(ratio[(j-1):(j+1)], na.rm = TRUE)
   names(r) <- names(ratio)[2:(length(ratio)-1)]
   lines(as.numeric(names(r)), r, lwd = 2, col = cols[i])
}
mtext("Carapace width (mm)", 1, 2.5, cex = 1.25)
mtext("Missing leg rate (#/crab)", 2, 2.5, cex = 1.25)
box()
legend("topright", legend = c("Cheliped", paste("Walking Leg", 1:4)), lwd = 2, col = cols, bg = "white", cex = 1.25)

# Immature
k <- array(0, dim = c(length(years), 140, 5))
dimnames(k) <- list(year = years, width = 1:140, missing = 1:5)
n <- k
for (i in 1:5){
   for (j in 1:length(years)){
      index <- (b$year == years[j]) & (b$maturity < 0.5) & is.new.shell(b)
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
   r <- r[as.character(40:110)]
   lines(as.numeric(names(r)), r, lwd = 2, col = cols[i], lty = "dashed")
}
box()
dev.off()

#===============================================================================
k <- array(0, dim = c(length(years), 140, 5))
dimnames(k) <- list(year = years, width = 1:140, missing = 1:5)
n = k
for (i in 1:5){
   for (j in 1:length(years)){
      index <- (b$year == years[j]) & (b$maturity > 0.5) & is.new.shell(b)

      res <- aggregate(list(k = (M[index, i] == 1) | (M[index, i+5] == 1)), by = b[index, "cw", drop = FALSE], sum)
      res$n <- aggregate(list(n = (M[index, i] == 1) | (M[index, i+5] == 1)), by = b[index, "cw", drop = FALSE], length)$n
      res <- res[res$cw <= 140, ]
      tmp <- colorRamp(c("black", cols[i]))(0.75)
      tmp <- rgb(tmp[1,1]/255, tmp[1,2]/255, tmp[1,3]/255)
      k[j, res$cw, i] <- res$k
      n[j, res$cw, i] <- res$n
   }
}


I <- (k[,,1] / n[,,1])[, 39:139]
II <- NULL
for (i in 2:(ncol(I)-1)) II <- cbind(II, apply(I[, (i-1):(i+1)], 1, mean, na.rm = TRUE))
colnames(II) <- 40:138
II["1996", ] <- NA

# Image version:
tiff(file = paste0("Missing Legs Poster/Missing legs by cw x year.tiff"), compression = "lzw", units = "in", res = 300, height = 7.5, width = 7.5)
image(40:138, years, t(II), col = colorRampPalette(c("white", "black"))(100),
      breaks = c(seq(0, 0.25, len = 100), 1.1),
      xlab = "", ylab = "", cex.lab = 1.25,
      xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
      xlim = c(39.5, 130.5), ylim = c(1989.5, 2022.5))
axis(2, at = seq(1990, 2022, by = 4), cex.axis = 1, las = 2)
axis(2, at = seq(1992, 2022, by = 4), cex.axis = 1, las = 2)
axis(1, at = seq(40, 140, by = 10))
axis(1, at = seq(45, 140, by = 10))
mtext("Carapace width (mm)", 1, 2.25, cex = 1.25)
mtext("Years", 2, 3.25, cex = 1.25)
box()
lens <- 60:131
for (j in 1:length(years)) lines(par("usr")[1:2], rep(years[j],2)+0.5, lwd = 0.5, col = "grey60")
for (i in 1:length(lens)) lines(rep(lens[i],2)+0.5, par("usr")[3:4], lwd = 0.5, col = "grey60")
box()
dev.off()

# Image version with biomass time series:
windows(width = 11, height = 7)
jpeg(filename = "Missing Legs - Year x width immature.jpeg",
     width = 12*480, height = 6.25*480, quality = 75, res = 600)
image(40:138, years, t(II), col = colorRampPalette(c("white", "black"))(100),
      breaks = c(seq(0, 0.25, len = 100), 1.1),
      xlab = "", ylab = "", cex.lab = 1.3, cex.axis = 1.1,
      xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
      xlim = c(49.5, 110.5), ylim = c(1989.5, 2017.5))
axis(2, at = seq(1990, 2018, by = 2))
axis(1, at = seq(60, 140, by = 5))
box()
mtext("Carapace width (mm)", 1, 2.5, cex = 1.35)
mtext("Years", 2, 2.5, cex = 1.35)
lens <- 60:131
for (j in 1:length(years)) lines(c(59.5, par("usr")[2]), rep(years[j],2)+0.5, lwd = 0.5, col = "grey60")
for (i in 1:length(lens)) lines(rep(lens[i],2)-0.5, par("usr")[3:4], lwd = 0.5, col = "grey60")
box()

# Calculate biomass:
res <- aggregate(list(n = (b$cw >= 95) & (b$maturity >= 0.5)), b[c("year", "tow.id")], function(x) return(sum(x)))
res <- sort(res, by = c("year", "tow.id"))
index <- match(res[c("year", "tow.id")], s[c("year", "tow.id")])
s$n <- 0
s$n[index] <- res$n
res <- aggregate(list(density = 1000000 * s$n / s$swept.area), by = s[c("year")], mean)
res$density[res$year == 1996] <- NA
lines(50 + 10 * (1 - (res$density / 6000)), res$year, lwd = 2)
lines(c(60, 60)-0.5, par("usr")[3:4])
mtext("Density", side = 1, at = 55, cex = 1.25, padj = 0.5)
dev.off()



