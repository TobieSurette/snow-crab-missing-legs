library(gulf.data)
library(gulf.graphics)
clm()
clg()

# Define variables:
years <- 1988:2025

# Read snow crab survey biological data:
b <- read.scsbio(years, survey = "regular")
b$year <- year(b)
b <- b[which((b$carapace.width > 0) & (b$carapace.width <= 150)), ]
b$tow.id <- tow.id(b)
b$cw <- round(b$carapace.width)
b$cw2 <- round(b$carapace.width/2)*2

# Format durometer readings:
b$durometer <- gsub("[*]", "", b$durometer)
b$durometer <- as.numeric(b$durometer)
b$durometer[which((b$durometer <= 0) & (b$durometer > 100))] <- NA

# Fix missing legs:
b$missing.legs[b$missing.legs == "*********"] <- "**********"

# Define missing leg patterns:
m <- matrix(0, ncol = 10, nrow = nrow(b))
for (i in 1:10) m[, i] <- as.numeric(substr(b$missing.legs, i, i))
m[is.na(m)] <- 0

# Time series by maturity and new/old-shelled:
png(file = "results/figures/old-shelled immature rates.png", res = 500, width = 7, height = 9, units = "in")
l <- rbind(0, cbind(0, kronecker(1:2, matrix(1, ncol = 5, nrow = 5)), 0), 0)
layout(l)
cols <- fade(c("red", "orange", "gold", "green3", "blue"), 0.5)
for (j in 1:2){
   ix <- which((b$sex == j) & !is.new.shell(b) & !is.mature(b))

   r <- aggregate((m == 1)[ix, ], by = b[ix, "year", drop = FALSE], mean)

   par(mar = c(0,0,0,0))
   plot(range(years), c(0, 12), type = "n", yaxs = "i", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
   grid()
   for (i in 1:ncol(m)){
      if (i >= 6) lty <- "dashed" else lty <- "solid"
      lines(r$year, 100 * r[,i+1], col = cols[((i-1) %% 5) + 1], lwd = 2, lty = lty)
   }
   if (j == 1){
      legend("topleft",
             legend = rev(c("Cheliped", paste0("Leg ", 2:5))),
             lwd = 2, cex = 1.25,
             col = rev(cols),
             bg = NA,
             box.col = "grey50", box.lwd = 0.5)

      mtext("Rate (%)", 2, 2.5, at = 0, font = 2, cex = 1.25)
      axis(2)
   }
   mtext(sex(j), 4, 1.5, font = 2, cex = 1.25)

   if (j == 2) axis(2, at = seq(0, 10, by = 2))

   box(col = "grey50")
}
mtext("Year", 1, 2.75, font = 2, cex = 1.25)
axis(1, at = seq(1990, 2025, by = 5))
dev.off()

# Time series by maturity and shell condition:
leg <- 1
png(file = paste0("results/figures/Immature leg ", leg, " loss rates.png"), res = 500, width = 7, height = 9, units = "in")
l <- rbind(0, cbind(0, kronecker(1:2, matrix(1, ncol = 5, nrow = 5)), 0), 0)
layout(l)
ylim = c(0, 24)
cols <- fade(c("green3", "gold", "orange", "red"), 0.5)
if (leg == 1) title.str <- "Cheliped loss" else title.str <- paste0("Walking leg ", leg, " loss")
for (j in 1:2){
   # Compile stats:
   legs <- c(leg, leg+5)
   r <- data.frame(year = years, SC12 = NA, SC3 = NA, SC4 = NA, SC5 = NA)
   ix <- (b$sex == j) & !is.mature(b)
   tmp <- aggregate((m == 1)[ix & (b$shell.condition %in% 1:2), legs], by = b[ix & (b$shell.condition %in% 1:2), "year", drop = FALSE], mean)
   r$SC12[match(tmp$year, r$year)] <- apply(tmp[, 2:3], 1, mean)
   tmp <- aggregate((m == 1)[ix & (b$shell.condition %in% 3), legs], by = b[ix & (b$shell.condition %in% 3), "year", drop = FALSE], mean)
   r$SC3[match(tmp$year, r$year)] <- apply(tmp[, 2:3], 1, mean)
   tmp <- aggregate((m == 1)[ix & (b$shell.condition %in% 4), legs], by = b[ix & (b$shell.condition %in% 4), "year", drop = FALSE], mean)
   r$SC4[match(tmp$year, r$year)] <- apply(tmp[, 2:3], 1, mean)
   tmp <- aggregate((m == 1)[ix & (b$shell.condition %in% 5), legs], by = b[ix & (b$shell.condition %in% 5), "year", drop = FALSE], mean)
   r$SC5[match(tmp$year, r$year)] <- apply(tmp[, 2:3], 1, mean)

   par(mar = c(0,0,0,0))
   plot(range(years), ylim, type = "n", yaxs = "i", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
   grid()
   for (i in 1:4){
      lines(r$year, 100 * r[,i+1], col = cols[((i-1) %% 5) + 1], lwd = 2, lty = "solid")
   }
   if (j == 1){
      legend("topleft",
             legend = rev(c("SC1&2", "SC3", "SC4", "SC5")),
             lwd = 2, cex = 1.25,
             col = rev(cols),
             bg = NA,
             box.col = "grey50", box.lwd = 0.5)

      mtext("Rate (%)", 2, 2.5, at = 0, font = 2, cex = 1.25)
      axis(2, at = seq(ylim[1], ylim[2], by = 2))

      mtext(title.str, 3, 0.5, font = 2, cex = 1.25)
   }
   mtext(sex(j), 4, 1.0, font = 2, cex = 1.25)

   if (j == 2) axis(2, at = seq(0, ylim[2]-2, by = 2))

   box(col = "grey50")
}
mtext("Year", 1, 2.75, font = 2, cex = 1.25)
axis(1, at = seq(1990, 2025, by = 5))
dev.off()

# Time series by maturity and carapace width:
leg <- 1
png(file = paste0("results/figures/Immature leg ", leg, " vs cw loss rates.png"), res = 500, width = 7, height = 9, units = "in")
l <- rbind(0, cbind(0, kronecker(c(1, 1, 2), matrix(1, ncol = 5, nrow = 5)), 0), 0, 0)
layout(l)

cols <- fade(c("green3", "gold", "orange", "red"), 0.5)
if (leg == 1) title.str <- "Cheliped loss" else title.str <- paste0("Walking leg ", leg, " loss")
logit <- function(x) return(log(x/(1-x)))
for (j in 1:2){
   # Compile stats:
   legs <- c(leg, leg+5)
   ix <- (b$sex == j) & !is.mature(b)
   r <- matrix(NA, nrow = length(years), ncol = 76)
   dimnames(r) <- list(year = years, cw = seq(0, 150, by = 2))

   tmp <- aggregate(m[ix, legs], by = b[ix, c("year", "cw2")], mean)
   tmp$value <- (tmp$V1 + tmp$V2) / 2
   for (i in 1:length(years)){
      r[as.character(years[i]), as.character(tmp$cw2[tmp$year == years[i]])] <- tmp$value[tmp$year == years[i]]
   }
   r[r == 0 | r == 1] <- NA
   r <- (logit(r) - mean(logit(r), na.rm = TRUE)) / sd(logit(r), na.rm = TRUE)
   breaks <- seq(-3, 3, by = 0.1)
   cols <- fade(colorRampPalette(c("red", "orange", "yellow", "green2", "blue"))(length(breaks)-1), 0.8)
   cols <- fade(colorRampPalette(c("red", "white", "blue"))(length(breaks)-1), 0.8)

   par(mar = c(0,0,0,0))
   #breaks = c(-7, -6 , -5, -4.5, -4, -3.75, -3.5, -3.25, -3, -2.5, -2, -1, 0)
   if (j == 1) ylim <- c(40, 120) else ylim <- c(40, 80)
   image(years, as.numeric(colnames(r)), r, ylim = ylim, xaxt = "n", yaxt = "n",
         breaks = breaks, col = cols)

   if (j == 1) mtext("Carapace width (mm)", 2, 2.5, at = ylim[1] +20, font = 2, cex = 1.25)
   if (j == 1) axis(2, at = seq(ylim[1], ylim[2], by = 10))
   if (j == 2) axis(2, at = seq(ylim[1], ylim[2]-10, by = 10))
   if (j == 1) mtext(title.str, 3, 0.5, font = 2, cex = 1.25)
   mtext(sex(j), 4, 1.0, font = 2, cex = 1.25)

   box(col = "grey50")
}
mtext("Year", 1, 2.75, font = 2, cex = 1.25)
axis(1, at = seq(1990, 2025, by = 5))
dev.off()
