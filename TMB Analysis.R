library(gulf)
library(TMB)

setwd("U:/Snow Crab/Missing Legs")
source("U:/TMB/TMB utilities.R")

years <- 1989:2019

# Read survey data:
s <- read.scset(year = years, valid = 1)
s$station <- as.numeric(substr(s$tow.id, 3, 5))

# Read biological observations:
b <- read.scbio(year = years)
b <- b[b$sex == 1, ]
b$carapace.width <- round(b$carapace.width)
b$maturity <- (is.mature(b, prob = TRUE) >= 0.5) + 1 - 1
b <- b[!is.na(b$carapace.width) & !is.na(b$shell.condition) & !is.na(b$maturity), ]
b$station <- as.numeric(substr(b$tow.id, 3, 5))
 
# Load survey grids:
mif <- read.mif("grids2017.mif", mid = FALSE)
for (i in 1:length(mif)){
   tmp <- km2deg(mif[[i]]$x, mif[[i]]$y)
   mif[[i]]$x <- tmp$longitude
   mif[[i]]$y <- tmp$latitude
}
mif <- as.polygon(mif)

# Assign grids to each survey station:
s$grid <- NA
for (i in 1:nrow(s)){
   if ((i %% 25) == 0) print(i)
   flag <- FALSE
   j <- 1
   while (!flag & (j <= length(mif))){
      if (in.polygon(mif[[j]], longitude(s[i,]), latitude(s[i,]))){
         s$grid[i] <- j
         flag <- TRUE
      }
      j <- j + 1
   }
}

index <- match(b[c("year", "tow.id")], s[c("year", "tow.id")])
b$grid <- s$grid[index]
b <- b[!is.na(b$grid), ]

# Missing leg matrix:
M <- matrix(NA, ncol = 10, nrow = nrow(b))
parity <- M
position <- M
for (i in 1:10){
   M[,i] <- (substr(b$missing.legs, i, i) == "1") + 1 - 1
   parity[,i] <- (i > 5) + 1 - 1
   position[,i] <- ((i-1) %% 5) + 1
}

res <- aggregate(apply(M, 1, sum), by = b[c("year", "grid")], mean)

clg()
years <- sort(unique(res$year))
for (j in 26:31){ #length(years)){
   index <- which(res$year == years[j])
   windows()
   gulf.map()
   v <- res[index,3] / 1.25
   v[v >= 1] <- 1
   for (i in 1:length(index)){
      plot(as.polygon(mif[res[index[i],2]]), col = grey(1 - v[i]))
   }
   mtext(years[j], 3, 2.0, cex = 1.5)
}

clg()
years <- sort(unique(res$year))
for (j in 11:21){ #length(years)){
   index <- which(res$year == years[j])
   windows()
   gulf.map()
   for (i in 1:length(index)){
      plot(as.polygon(mif[res[index[i],2]]), col = grey(1 - (res[index[i],3]/max(res[index,3]))))
   }
}

clg()
res <- aggregate(list(commercial = (b$carapace.width >= 95) & (b$maturity == 1)), by = b[c("year", "grid")], sum)
res$mature <- aggregate(list(mature = b$maturity), by = b[c("year", "grid")], sum)$mature
res$p <- res$commercial / res$mature
years <- sort(unique(res$year))
for (j in 26:31){ #length(years)){
   index <- which(res$year == years[j] & !is.na(res$p))
   windows()
   gulf.map()
   for (i in 1:length(index)){
      plot(as.polygon(mif[res[index[i],2]]), col = grey(1 - res$p[index[i]]))
   }
   mtext(years[j], 3, 2.0, cex = 1.5)
}


res <- aggregate(list(missing = apply(M, 1, sum)), by = b["year"], sum)
res$regenerated <- aggregate(list(missing = apply(M, 1, sum)), by = b["year"], sum)[, 2]
res$human <- aggregate(list(missing = apply(M, 1, sum)), by = b["year"], sum)[, 2]
res$total <- aggregate(list(b$year), by = b["year"], length)[, 2]


# Expand variables:
maturity <- repvec(b$maturity, ncol = 10)
sc <- repvec(b$shell.condition, ncol = 10)
cw <- repvec(b$carapace.width, ncol = 10)
grid <- repvec(b$grid, ncol = 10)
station <- repvec(b$station, ncol = 10)

# Define data:
data <- data.frame(z = as.numeric(M),                # Indicator vector of missing pereopods. 
                   width = as.numeric(cw),           # Vector of observed carapace widths IDs.
                   parity = as.numeric(parity),      # Vector indicating the side of the crab. 
                   position = as.numeric(position),  # Vector indicating the position of the pereopod.
                   maturity = as.numeric(maturity),  # Vector of crab maturity IDs.
                   carapace = as.numeric(sc),        # Vector of carapace condition categories.
                   grid = as.numeric(grid))          # Vector of sampling grid IDs.        

# Sub-sample:
#index <- sort(sample(1:length(data$z), 25000))

#data <- lapply(data, function(x) x[index])
   // Define leg-position cross-correlation parameter matrix:
   matrix<Type> alpha_position(n_position, n_position);
   for (int i = 0; i < n_position; i++){
      // Define diagonal elements (5 parameters):
      alpha_position(i,i) = position_effect[i];
      alpha_position(i + n_position, i + n_position) = position_effect[i];

      // Read cross-probabilities left side (10 parameters):
      k = 5;
      for (int j = 0; j < i; j++){
         alpha_position(i,j) = position_effect[k];
         alpha_position(j,i) = alpha_position(i,j);         
         k = k + 1;
      }
      // Read cross-probabilities right side (10 parameters):
      for (int j = 0; j < i; j++){ 
         alpha_position(i+n_position, j+n_position) = position_effect[k];
         alpha_position(j+n_position, i+n_position) = alpha_position(i+n_position, j+n_position);
         k = k + 1;
      }
   }
   
# Translate into factors:
#data$width <- match(data$width, sort(unique(data$width)))
#data$grid <- match(data$grid, sort(unique(data$grid)))

# Define initial parameter values:
parameters = list(log_alpha = 0,                                         # Log-scale average missing pereopod rate.
                  width_effect = rep(0, max(data$width)),                # Crab size effect effect.
                  position_effect = rep(0, 5),                           # Pereopod position effect.
                  parity_effect = 0,                                     # Crab side effect.
                  maturity_effect = 0,                                   # Maturity effect.
                  carapace_effect = rep(0,5),                            # Carapace condition effect.
                  grid_effect = rep(0, max(data$grid)),                  # Sampling station effect.
                  log_sigma_width = 0,                                   # Log-scale error parameter for length effects.
                  log_sigma_position = 0,                                # Log-scale error parameter for position effects.
                  log_sigma_carapace = 0,                                # Log-scale error parameter for carapace condition effects.
                  log_sigma_grid = 0,                                    # Log-scale error parameter for sampling grid effects.
                  width_position_effect = rep(0, max(data$width) * 5),   # Width x position interaction term.
                  log_sigma_width_position = 0,                          # Log-scale error for width x position interaction term.
                  width_maturity_effect = rep(0, max(data$width)),       # Width x position interaction term.
                  log_sigma_width_maturity = 0,                          # Log-scale error for width x position interaction term.
                  grid_maturity_effect = rep(0, max(data$grid)),         # Grid x position interaction term.
                  log_sigma_grid_maturity = 0,                           # Log-scale error for grid x position interaction term.                  
                  position_maturity_effect = rep(0, 5),                  # Position x maturity interaction term.
                  log_sigma_position_maturity = 0,                       # Log-scale error for position x maturity interaction term.    
                  grid_position_effect = rep(0, max(data$grid) * 5),     # Grid x position interaction term.
                  log_sigma_grid_position = 0)                           # Log-scale error for grid x position interaction term.
                                  

# Function to update parameters:
update.parameters <- function(x, fixed, random){
   if (!missing(fixed)){
      for (i in 1:length(x)){
         if (names(x[i]) %in% rownames(fixed)){
            x[[i]] <- as.numeric(fixed[rownames(fixed) == names(x[i]), 1])
         }
      }
   }
    
   if (!missing(random)){
      for (i in 1:length(x)){
         if (names(x[i]) %in% rownames(random)){
            x[[i]] <- as.numeric(random[names(x[i]) == rownames(random), 1])
         }
      }
   }
   
   return(x)   
}

# m <- gam(z ~  as.factor(grid), data = data, family = binomial)

# Compile and load 
cpp.file <- "Missing_Pereopod.cpp" 
compile(cpp.file)
dll.file <- gsub(".cpp", "", cpp.file)  
dyn.load(dynlib(dll.file))

# Define TMB object:
active.parameters <- c("log_alpha", "width_effect", "log_sigma_width")
obj <- MakeADFun(data = data[data.cpp(cpp.file)], 
                 parameters = parameters[parameters.cpp(cpp.file)], 
                 map = control.parameters(parameters, active.parameters),
                 random = "width_effect",
                 DLL = dll.file)

theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))

rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")

parameters <- update.parameters(parameters, fixed, random)

active.parameters <- unique(c(active.parameters, "position_effect", "log_sigma_position", "carapace_effect", "log_sigma_carapace"))
obj <- MakeADFun(data = data[data.cpp(cpp.file)], 
                 parameters = parameters[parameters.cpp(cpp.file)],
                 map = control.parameters(parameters, active.parameters),
                 random = c("width_effect", "position_effect", "carapace_effect"),
                 DLL = dll.file)
               
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
#theta <- optim(theta$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))

rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")

parameters <- update.parameters(parameters, fixed, random)


active.parameters <- c(active.parameters, "maturity_effect", "parity_effect")
obj <- MakeADFun(data = data[data.cpp(cpp.file)], 
                 parameters = parameters[parameters.cpp(cpp.file)],
                 map = control.parameters(parameters, active.parameters),
                 random = c("width_effect", "position_effect", "carapace_effect"),
                 DLL = dll.file)
               
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 500, trace = 3))
#theta <- optim(theta$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))

rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")

parameters <- update.parameters(parameters, fixed, random)

active.parameters <- unique(c(active.parameters, "grid_effect", "log_sigma_grid"))
obj <- MakeADFun(data = data[data.cpp(cpp.file)], 
                 parameters = parameters[parameters.cpp(cpp.file)],
                 map = control.parameters(parameters, active.parameters),
                 random = c("width_effect", "position_effect", "carapace_effect", "grid_effect"),
                 DLL = dll.file)
               
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 1000, trace = 3))
#theta <- optim(theta$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))

rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")

parameters <- update.parameters(parameters, fixed, random)

#===================================================================================================
active.parameters <- unique(c(active.parameters, "width_position_effect", "log_sigma_width_position", "width_maturity_effect", "log_sigma_width_maturity"))
obj <- MakeADFun(data = data[data.cpp(cpp.file)], 
                 parameters = parameters[parameters.cpp(cpp.file)],
                 map = control.parameters(parameters, active.parameters),
                 random = c("width_effect", "position_effect", "carapace_effect", "grid_effect", "width_position_effect", "width_maturity_effect"),
                 DLL = dll.file)
               
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 1000, trace = 3))
#theta <- optim(theta$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))

rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")

parameters <- update.parameters(parameters, fixed, random)

#===================================================================================================
active.parameters <- unique(c(active.parameters, "grid_maturity_effect", "log_sigma_grid_maturity"))
obj <- MakeADFun(data = data[data.cpp(cpp.file)], 
                 parameters = parameters[parameters.cpp(cpp.file)],
                 random = c("width_effect", "position_effect", "carapace_effect", "grid_effect", "width_position_effect", "width_maturity_effect", "grid_maturity_effect"),
                 DLL = dll.file)
               
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
#theta <- optim(theta$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))

rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")

parameters <- update.parameters(parameters, fixed, random)

grid_effects <- parameters$grid_effects
dbarplot(grid_effect)



dbarplot(random[rownames(random) == "width_effect", 1], width = 1)
dbarplot(random[rownames(random) == "width_maturity_effect", 1], width = 1)
m <- parameters$width_position_effect
dim(m) <- c(5, length(parameters$width_effect))
plot(c(0, 140), c(0.5, 5.5), type = "n", xlab = "Carapace width (mm)", ylab = "Position")
for (i in 1:5){
   mm <- m[i,]
   index <- mm >= 0
   points((1:ncol(m))[index], rep(i,sum(index)), cex = 30 * sqrt(abs(mm[index])), pch = 21, bg = "white")
   points((1:ncol(m))[!index], rep(i,sum(!index)), cex = 30 * sqrt(abs(mm[!index])), pch = 21, bg = "black")
}

dbarplot(random[rownames(random) == "position_effect", 1])
dbarplot(random[rownames(random) == "grid_effect", 1])
dbarplot(random[rownames(random) == "grid_maturity_effect", 1])

dbarplot(random[grep("width_effect", rownames(random)) , 1])
dbarplot(random[grep("width_maturity_effect", rownames(random)) , 1])
dbarplot(random[grep("grid_effect", rownames(random)) , 1])

s$effect <- NA
s$effect <- parameters$station_effect[s$station]

# Grid_effect:
gulf.map(sea = FALSE)
g <- parameters$grid_effect
for (i in 1:length(g)){
   if (g[i] >= 0) col = colorRamp(c("white", "black"))(g[i]/max(abs(g)))
   if (g[i] < 0) col = colorRamp(c("white", "red"))(-g[i]/max(abs(g)))
   col <- col / 255
   col <- rgb(col[1,1], col[1,2], col[1,3])
   polygon(mif[[i]]$x, mif[[i]]$y, col = col)
}
bathymetry(dem = FALSE)

# Grid x maturity effect:
gulf.map(sea = FALSE)
g <- parameters$grid_maturity_effect
for (i in 1:length(g)){
   if (g[i] >= 0) col = colorRamp(c("white", "black"))(g[i]/max(abs(g)))
   if (g[i] < 0) col = colorRamp(c("white", "red"))(-g[i]/max(abs(g)))
   col <- col / 255
   col <- rgb(col[1,1], col[1,2], col[1,3])
   polygon(mif[[i]]$x, mif[[i]]$y, col = col)
}
bathymetry(dem = FALSE)
points(longitude(s), latitude(s), pch = 21, bg = "grey")



index <- parameters$station_effect < 0
points(longitude(s)[index], latitude(s)[index], pch = 21, bg = "grey", cex = 3*sqrt(abs(s$effect[index])))
index <- parameters$station_effect >= 0
points(longitude(s)[index], latitude(s)[index], pch = 21, bg = "black", cex = 3*sqrt(abs(s$effect[index])))

