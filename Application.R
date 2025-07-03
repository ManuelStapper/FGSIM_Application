library("surveillance")
library("sf")
library("scoringutils")
library("spdep")
library("ggplot2")
library("dplyr")
library("geosphere")

# All of the code below needs to be reviewed after changing the sampling

# Setting the working directory here

source("RawWeights.R")
source("DerivativeHelper.R")
source("Powerlaw.R")
source("W_fgsim.R")
source("W_centroid.R")
source("RawWeightsMTL.R")
source("WeightsFromFit.R")


### Reading in some data

# Population data
popMat = as.matrix(read.csv("Data/population.csv", header = T, sep = ";")[, -1])
# Dates observed
dates = as.matrix(read.csv("Data/dates.csv", header = T))

# Transform population data
# Idea: Assume constant total population but split across districts
# Also: Interpolate to obtain weekly data
# Reference data: End of year
# Do rescaling in place
for(i in 1:20){
  popMat[, i] = popMat[, i] / sum(popMat[, i]) * sum(popMat[, 12])
}
POP = matrix(NA,nrow = 1036, ncol = 401)
POP[1, ] = popMat[, 1]

for(i in 2:20){
  r = min(which(dates[, 1] == 2000 + i))
  POP[r, ] = popMat[, i]
}
ind = which(!is.na(POP[, 1]))
for(r in 993:1036){
  POP[r, ] = POP[992, ]
}

for(i in 1:19){
  for(j in 1:401){
    i0 = ind[i]
    i1 = ind[i+1]
    POP[i0:i1, j] = seq(POP[i0, j], POP[i1, j], length = i1 - i0 + 1)
  }
}


# Infection country
cases = t(as.matrix(read.csv("Data/cases.csv", header = T)))
# Shapefiles for administrative districts
# Could be replaced by the data used in getLNpars.R?
Ger = read_sf(dsn = "Data/Shapefile", layer = "LKGermany")

ID2Geom = read.csv("Data/IDlookup.csv", header = T)
IDs = sort(unique(ID2Geom[, 1]))
dateName = paste(dates[, 1], dates[, 2], sep = "-")

shp = Ger$geometry[ID2Geom$GeomID1]
shp[IDs == 3159] = st_union(Ger$geometry[c(306, 372)])

A = nb2mat(poly2nb(shp, snap = 1)) > 0
O = nbOrder(A)

colnames(A) = IDs
rownames(A) = IDs
colnames(O) = IDs
rownames(O) = IDs
rownames(popMat) = IDs
colnames(POP) = IDs
rownames(POP) = dateName
colnames(cases) = IDs
rownames(cases) = dateName

########################
### Data preparation ###
########################

cols = c(rgb(253, 185, 19, maxColorValue = 255),
         rgb(0, 160, 226, maxColorValue = 255),
         rgb(89, 18, 68, maxColorValue = 255),
         rgb(0, 191, 111, maxColorValue = 255),
         rgb(32, 154, 141, maxColorValue = 255),
         rgb(0, 49, 81, maxColorValue = 255))

# Restrict to southern germany states
indS = which(floor(IDs/1000) %in% c(8, 9))
indN = which(floor(IDs/1000) %in% c(1, 2, 3, 4))
indE = which(floor(IDs/1000) %in% c(11, 12, 13, 14, 15, 16))
indW = which(floor(IDs/1000) %in% c(5, 6, 7, 10))
ind = list(indN, indE, indS, indW)

### Caution: Neighbourhood order depends on the observation region!
### Overwrite the ngb order matrix!

O[indN, indN] = nbOrder(nb2mat(poly2nb(shp[indN], snap = 1)) > 0)
O[indE, indE] = nbOrder(nb2mat(poly2nb(shp[indE], snap = 1)) > 0)
O[indS, indS] = nbOrder(nb2mat(poly2nb(shp[indS], snap = 1)) > 0)
O[indW, indW] = nbOrder(nb2mat(poly2nb(shp[indW], snap = 1)) > 0)

#############################################################
### Interlude: Examples why NGB order is not always great ###
#############################################################

cols = c("#FDB913", "#00A0E2", "#591244", "#00BF6F", "#209A8D", "#003151")

# 1)
# District of Göttingen = Göttingen + Osterode am Harz

# png("Plots/Example1.png", res = 350)
# plot(shp[O[26, ] == 1], col = c(cols[1], rep(cols[2], 3), cols[1], rep(cols[2], 2)), border = "white", lwd = 2)
# plot(Ger$geometry[c(306, 372)], border = "red", add = T, col = cols[2], lwd = 2)
# plot(shp[IDs == 3159], add = T, lwd = 2)
# dev.off()

png("Plots/Example1.png", 
    width = 6*1.1, height = 6, units = "in", res = 350)
# Minimal margins for maps
par(mar = c(0, 0, 0, 0))
plot(shp[O[26, ] == 1], col = c(cols[1], rep(cols[5], 3), cols[1], rep(cols[5], 2)), 
     border = "white", lwd = 2, axes = FALSE)  # Remove axes for cleaner look
plot(Ger$geometry[c(306, 372)], border = "red", add = T, col = cols[5], lwd = 2)
plot(shp[IDs == 3159], add = T, lwd = 2)
dev.off()

# 2)

# Find good example for connecting district
OW = nbOrder(nb2mat(poly2nb(shp[indW], snap = 1)) > 0)
OT = nbOrder(nb2mat(poly2nb(shp, snap = 1)) > 0)[indW, indW]
indTemp = which((OW - OT) == max(OW - OT), arr.ind = T)
indTemp = sort(unique(c(which(OW[33, ] <= 1), which(OW[40, ] <= 2))))

# png("Plots/Example2.png")
# plot(shp[indW[indTemp]], col = cols[4], border = "white", lwd = 2)
# plot(shp[c(48, 58)], add = T, col = cols[2], border = "white", lwd = 2)
# plot(shp[indW[33]], col = cols[4], border = "black", lwd = 2, add = T)
# plot(shp[indW[40]], col = cols[4], border = "black", lwd = 2, add = T)
# dev.off()

png("Plots/Example2.png",
    width = 8.96, height = 6, units = "in", res = 350)
par(mar = c(0, 0, 0, 0))
plot(shp[c(indW[indTemp], 48, 58)], border = "white")
plot(shp[indW[indTemp]], col = cols[5], border = "white", lwd = 2, add = T)
plot(shp[c(48, 58)], add = T, col = cols[1], border = "white", lwd = 2)
plot(shp[indW[33]], col = cols[5], border = "black", lwd = 2, add = T)
plot(shp[indW[40]], col = cols[5], border = "black", lwd = 2, add = T)
dev.off()

# 3)

load("Parameters.RData")

D = parsBs_E$mu
exp(parsB_E$u[22, 27])
exp(parsB_E$u[59, 71])
D[22, 27]
D[59, 71]

png("Plots/Example3.png",
    width = 15.52, height = 6, units = "in", res = 350)
par(mar = c(0, 0, 0, 0))
plot(shp[indE[c(22, 27)]], col = cols[5], border = "white", lwd = 2, xlim = c(9.9, 13.8))
shpTemp <- st_geometry(shp[indE[c(59, 71)]]) + c(-1.2, 2.5)
plot(shpTemp, add = T, col = cols[5], border = "white", lwd = 2)
rect(xleft = 9.9, ybottom = 53.29, xright = 10.5, ytop = 53.62)
dev.off()

# 4) 

pairMin = c(0, 0)
currMin = Inf
pairMax = c(0, 0)
currMax = 0
for(i in 1:401){
  indTemp = which(O[i, ] == 1)
  for(j in indTemp){
    temp = as.numeric(st_length(st_intersection(st_boundary(shp[i]), st_boundary(shp[j]))))
    if(min(temp[temp > 1000]) < currMin){
      currMin = min(temp[temp > 0])
      pairMin = c(i, j)
    }
    if(max(temp) > currMax){
      currMax = max(temp)
      pairMax = c(i, j)
    }
  }
}

png("Plots/Example4.png",
    width = 8.69, height = 6, units = "in", res = 350)
par(mar = c(0.5, 0.5, 1, 0.5))
plot(shp[c(351, 338, 346, 350)], col = cols[5], border = "white", lwd = 2)
dev.off()

# 5)

which(ID2Geom$Name == "Mainz, kfS")
which(ID2Geom$Name == "Wiesbaden, kfS")

plot(shp[c(120, 163)])

which(ID2Geom$Name == "Uckermark")
which(ID2Geom$Name == "Vorpommern-Greifswald")

plot(shp[c(343, 350)], col = cols[1], border = "white", lwd = 2)

library("raster")
popDens = raster("PopDistrict/DEU2011.tif")
vals = values(popDens)
vals[!is.na(vals)] = log(vals[!is.na(vals)] + 1)
values(popDens) = vals


png("Plots/Example5a.png",
    width = 6, height = 6, units = "in", res = 350)
par(bty = "n", mar = c(0.5, 0.5, 1, 0.5))
plot(mask(crop(popDens, as_Spatial(shp[c(120, 163)])), as_Spatial(shp[c(120, 163)])), legend = F, axes = F)
plot(shp[c(120, 163)], add = T, border = "black", lwd = 2)
dev.off()

png("Plots/Example5b.png",
    width = 6, height = 6, units = "in", res = 350)
par(bty = "n", mar = c(0.5, 0.5, 1, 0.5))
plot(mask(crop(popDens, as_Spatial(shp[c(343, 350)])), as_Spatial(shp[c(343, 350)])), legend = F, axes = F)
plot(shp[c(343, 350)], add = T, border = "black", lwd = 2)
dev.off()


png("Plots/Example5.png",
    width = 8, height = 6, units = "in", res = 350)
layout(matrix(c(1, 2), nrow = 1, ncol = 2))
par(bty = "n", mar = c(0, 0, 0, 0))

plot(mask(crop(popDens, as_Spatial(shp[c(120, 163)])), as_Spatial(shp[c(120, 163)])), 
     legend = F, axes = F)
plot(shp[c(120, 163)], add = T, border = "black", lwd = 2)

# Right panel
plot(mask(crop(popDens, as_Spatial(shp[c(343, 350)])), as_Spatial(shp[c(343, 350)])), 
     legend = F, axes = F)
plot(shp[c(343, 350)], add = T, border = "black", lwd = 2)
dev.off()

#####################
### Interlude end ###
#####################

# Read in parameter data here
load("Parameters.RData")

centroids = st_centroid(shp)
mus = list()
popRatio = list()

for(r in 1:4){
  n = length(ind[[r]])
  popVec = colMeans(POP[, ind[[r]]])
  
  mu = matrix(0, nrow = n, ncol = n)
  pr = matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in i:n){
      mu[i, j] = distGeo(centroids[i][[1]][1:2], centroids[j][[1]][1:2])/1000
      mu[j, i] = mu[i, j]
      pr[i, j] = popVec[j]/popVec[i]
      pr[j, i] = 1/pr[i, j]
    }
  }
  mus[[r]] = mu
  popRatio[[r]] = pr
}

parsBc_N = list(mu = mus[[1]])
parsBc_E = list(mu = mus[[2]])
parsBc_S = list(mu = mus[[3]])
parsBc_W = list(mu = mus[[4]])

parsGc_N = list(mu1 = mus[[1]], mu2 = popRatio[[1]])
parsGc_E = list(mu1 = mus[[2]], mu2 = popRatio[[2]])
parsGc_S = list(mu1 = mus[[3]], mu2 = popRatio[[3]])
parsGc_W = list(mu1 = mus[[4]], mu2 = popRatio[[4]])


temp = colMeans(POP[, indN])
parsB_N$pop = parsT_N$pop = parsG_N$pop = parsC_N$pop = parsR_N$pop = diag(temp/exp(mean(log(temp))))
parsBs_N$pop = parsTs_N$pop = parsGs_N$pop = parsCs_N$pop = parsRs_N$pop = diag(temp/exp(mean(log(temp))))
parsBc_N$pop = parsGc_N$pop = diag(temp/exp(mean(log(temp))))

temp = colMeans(POP[, indE])
parsB_E$pop = parsT_E$pop = parsG_E$pop = parsC_E$pop = parsR_E$pop = diag(temp/exp(mean(log(temp))))
parsBs_E$pop = parsTs_E$pop = parsGs_E$pop = parsCs_E$pop = parsRs_E$pop = diag(temp/exp(mean(log(temp))))
parsBc_E$pop = parsGc_E$pop = diag(temp/exp(mean(log(temp))))

temp = colMeans(POP[, indS])
parsB_S$pop = parsT_S$pop = parsG_S$pop = parsC_S$pop = parsR_S$pop = diag(temp/exp(mean(log(temp))))
parsBs_S$pop = parsTs_S$pop = parsGs_S$pop = parsCs_S$pop = parsRs_S$pop = diag(temp/exp(mean(log(temp))))
parsBc_S$pop = parsGc_S$pop = diag(temp/exp(mean(log(temp))))

temp = colMeans(POP[, indW])
parsB_W$pop = parsT_W$pop = parsG_W$pop = parsC_W$pop = parsR_W$pop = diag(temp/exp(mean(log(temp))))
parsBs_W$pop = parsTs_W$pop = parsGs_W$pop = parsCs_W$pop = parsRs_W$pop = diag(temp/exp(mean(log(temp))))
parsBc_W$pop = parsGc_W$pop = diag(temp/exp(mean(log(temp))))

# Quick plot
colsTemp = sapply(1:401, function(i) ifelse(i %in% indN, 1, ifelse(i %in% indE, 2, ifelse(i %in% indS, 3, 4))))

shp_sf = st_as_sf(shp)
shp_sf$colsTemp = as.factor(c("North", "East", "South", "West")[colsTemp])
fillAlpha = 1

p = ggplot(data = shp_sf) +
  geom_sf(aes(fill = colsTemp), lwd = 0.2, show.legend = F, color = "white") +
  scale_fill_manual(values = alpha(cols, fillAlpha)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
p = p + annotate("text", x = 6.5, y = 54.50, label = "North", color = alpha(cols[2], fillAlpha), size = 8)
p = p + annotate("text", x = 14, y = 55.0, label = "East", color = alpha(cols[1], fillAlpha), size = 8)
p = p + annotate("text", x = 14.00, y = 50.00, label = "South", color = alpha(cols[3], fillAlpha), size = 8)
p = p + annotate("text", x = 6.5, y = 48.50, label = "West", color = alpha(cols[4], fillAlpha), size = 8)
p
ggsave("Plots/MapRegions.png", plot = p, width = 6, height = 8, dpi = 300)


# Define sts object for hhh4
influenza = lapply(ind, function(xx){
  sts(observed = cases[, xx],
      start = c(2001, 1),
      frequency = 52,
      population = POP[, xx],
      neighbourhood = O[xx, xx])
})

names(influenza) = c("N", "E", "S", "W")

# tInd = setdiff(2:1036, 434:485) # If run on complete data
tInd = setdiff(2:928, 434:485) # If run on only training data

model = list(
  ar = list(NULL),
  ne = list(f = addSeason2formula(~-1 + fe(1, unitSpecific = TRUE), period = 52, S = 1),
            weights = W_powerlaw(maxlag = 7, from0 = T),
            scale = NULL, normalize = T),
  end = list(f = addSeason2formula(~-1 + fe(1, unitSpecific = TRUE), period = 52, S = 1)),
  family = "NegBin1", subset = tInd,
  optimizer = list(stop = list(tol=1e-5, niter=100000),
                   regression = list(method="nlminb", iter.max = 1000),
                   variance = list(method="nlminb", iter.max = 1000)))

models = list()
fits = list()

for(i in 1:4){
  models[[i]] = list()
  fits[[i]] = list()
}
names(models) = names(fits) = c("N", "E", "S", "W")

for(i in 1:4){
  for(j in 1:11){
    models[[i]][[j]] = model
  }
}

# Not pretty but does the job. Problem with just running a loop:
# parameters are forwarded to weight function by importing the argument in "pars"
# Thus, weight function would import all parameters first whenever being evaluated

models$N[[1]]$ne$weights = diag(rep(1, length(indN)))
models$N[[2]]$ne$weights = W_powerlaw2(pars = parsB_N, maxlag = 50, from0 = T, initial = 2.5, popScale = T)
models$N[[3]]$ne$weights = W_fgsim(pars = parsB_N, truncPL = F, popScale = T, initial = 1.5)
models$N[[4]]$ne$weights = W_fgsim(pars = parsBs_N, truncPL = F, popScale = T, initial = 1.5, simplify = T)
models$N[[5]]$ne$weights = W_centroid(pars = parsBc_N, popScale = T, initial = c(1, 0))
models$N[[6]]$ne$weights = W_fgsim(pars = parsT_N, truncPL = F, popScale = T, initial = 2)
models$N[[7]]$ne$weights = W_fgsim(pars = parsG_N, truncPL = F, popScale = T, initial = c(1.5, 0), gravity = T)
models$N[[8]]$ne$weights = W_fgsim(pars = parsGs_N, truncPL = F, popScale = T, initial = c(1.5, 0), gravity = T, simplify = T)
models$N[[9]]$ne$weights = W_centroid(pars = parsGc_N, popScale = T, gravity = T, initial = c(1, 0, 0))
models$N[[10]]$ne$weights = W_fgsim(pars = parsC_N, truncPL = F, popScale = T, initial = 1)
models$N[[11]]$ne$weights = W_fgsim(pars = parsR_N, truncPL = F, popScale = T, initial = -1)

models$E[[1]]$ne$weights = diag(rep(1, length(indE)))
models$E[[2]]$ne$weights = W_powerlaw2(pars = parsB_E, maxlag = 50, from0 = T, initial = 2.5, popScale = T)
models$E[[3]]$ne$weights = W_fgsim(pars = parsB_E, truncPL = F, popScale = T, initial = 1.5)
models$E[[4]]$ne$weights = W_fgsim(pars = parsBs_E, truncPL = F, popScale = T, initial = 1.5, simplify = T)
models$E[[5]]$ne$weights = W_centroid(pars = parsBc_E, popScale = T, initial = c(1, 0))
models$E[[6]]$ne$weights = W_fgsim(pars = parsT_E, truncPL = F, popScale = T, initial = 2)
models$E[[7]]$ne$weights = W_fgsim(pars = parsG_E, truncPL = F, popScale = T, initial = c(1.5, 0), gravity = T)
models$E[[8]]$ne$weights = W_fgsim(pars = parsGs_E, truncPL = F, popScale = T, initial = c(1.5, 0), gravity = T, simplify = T)
models$E[[9]]$ne$weights = W_centroid(pars = parsGc_E, popScale = T, gravity = T, initial = c(1, 0, 0))
models$E[[10]]$ne$weights = W_fgsim(pars = parsC_E, truncPL = F, popScale = T, initial = 1)
models$E[[11]]$ne$weights = W_fgsim(pars = parsR_E, truncPL = F, popScale = T, initial = -1)

models$S[[1]]$ne$weights = diag(rep(1, length(indS)))
models$S[[2]]$ne$weights = W_powerlaw2(pars = parsB_S, maxlag = 50, from0 = T, initial = 2.5, popScale = T)
models$S[[3]]$ne$weights = W_fgsim(pars = parsB_S, truncPL = F, popScale = T, initial = 1.5)
models$S[[4]]$ne$weights = W_fgsim(pars = parsBs_S, truncPL = F, popScale = T, initial = 1.5, simplify = T)
models$S[[5]]$ne$weights = W_centroid(pars = parsBc_S, popScale = T, initial = c(1, 0))
models$S[[6]]$ne$weights = W_fgsim(pars = parsT_S, truncPL = F, popScale = T, initial = 2)
models$S[[7]]$ne$weights = W_fgsim(pars = parsG_S, truncPL = F, popScale = T, initial = c(1.5, 0), gravity = T)
models$S[[8]]$ne$weights = W_fgsim(pars = parsGs_S, truncPL = F, popScale = T, initial = c(1.5, 0), gravity = T, simplify = T)
models$S[[9]]$ne$weights = W_centroid(pars = parsGc_S, popScale = T, gravity = T, initial = c(1, 0, 0))
models$S[[10]]$ne$weights = W_fgsim(pars = parsC_S, truncPL = F, popScale = T, initial = 1)
models$S[[11]]$ne$weights = W_fgsim(pars = parsR_S, truncPL = F, popScale = T, initial = -1)

models$W[[1]]$ne$weights = diag(rep(1, length(indW)))
models$W[[2]]$ne$weights = W_powerlaw2(pars = parsB_W, maxlag = 50, from0 = T, initial = 2.5, popScale = T)
models$W[[3]]$ne$weights = W_fgsim(pars = parsB_W, truncPL = F, popScale = T, initial = 1.5)
models$W[[4]]$ne$weights = W_fgsim(pars = parsBs_W, truncPL = F, popScale = T, initial = 1.5, simplify = T)
models$W[[5]]$ne$weights = W_centroid(pars = parsBc_W, popScale = T, initial = c(1, 0))
models$W[[6]]$ne$weights = W_fgsim(pars = parsT_W, truncPL = F, popScale = T, initial = 2)
models$W[[7]]$ne$weights = W_fgsim(pars = parsG_W, truncPL = F, popScale = T, initial = c(1.5, 0), gravity = T)
models$W[[8]]$ne$weights = W_fgsim(pars = parsGs_W, truncPL = F, popScale = T, initial = c(1.5, 0), gravity = T, simplify = T)
models$W[[9]]$ne$weights = W_centroid(pars = parsGc_W, popScale = T, gravity = T, initial = c(1, 1, 2))
models$W[[10]]$ne$weights = W_fgsim(pars = parsC_W, truncPL = F, popScale = T, initial = 1)
models$W[[11]]$ne$weights = W_fgsim(pars = parsR_W, truncPL = F, popScale = T, initial = -1)


#############################################################
##### Stage 1: Simple fit, complete period, no forecast #####
#############################################################

# Sloppy way to handle sensitivity to initial values

fitFun = function(infl, model){
  if(length(model$ne$weights) != 4){
    return(hhh4(infl, model))
  }
  
  dInit <- model$ne$weights$initial[1]
  
  for(i in 0:30){
    model$ne$weights$initial[1] = dInit + i * 0.1
    
    # Try to fit the model
    res = tryCatch({
      return(hhh4(infl, model))
    }, error = function(e){
      print(paste("Bumping initial value to", dInit + (i + 1) * 0.1))
      return(NULL)
    })
    
    if(!is.null(res)){
      return(res)
    }
  }
  
  return(0)
}

for(r in 1:4){
  for(i in 1:11){
    fits[[r]][[i]] = fitFun(influenza[[r]], models[[r]][[i]])
    print(c(r, i))
    print(BIC(fits[[r]][[i]]))
    if(i > 1){
      temp = coefficients(fits[[r]][[i]])
      print(temp[startsWith(names(temp), "neweights")])
    }
  }
}

models$W[[3]]$ne$weights$initial = 0.5
fits$W[[3]] = fitFun(influenza$W, models$W[[3]])

unlist(lapply(fits$N, BIC))
unlist(lapply(fits$E, BIC))
unlist(lapply(fits$S, BIC))
unlist(lapply(fits$W, BIC))


# Use overall fit as starting values for the estimation below

for(r in 1:4){
  for(i in 1:11){
    temp = fits[[r]][[i]]$coefficients
    indWgt = startsWith(names(temp), "neweights")
    indFE = startsWith(names(temp), "ne.1.") | startsWith(names(temp), "end.1.")
    models[[r]][[i]]$start = list(temp[indFE])
  }
  if(i > 1){
    models[[r]][[i]]$ne$weights$initial = temp[indWgt]
  }
}

# Select seasonality order

Smat = array(0, dim = c(4, 11, 2))

for(r in 1:4){
  for(m in 1:11){
    S = c(1, 1)
    models[[r]][[m]]$ne$f = addSeason2formula(~-1 + fe(1, unitSpecific = TRUE), period = 52.1775, S = S[1])
    models[[r]][[m]]$end$f = addSeason2formula(~-1 + fe(1, unitSpecific = TRUE), period = 52.1775, S = S[2])
    tempFit = fitFun(influenza[[r]], models[[r]][[m]])
    if(!is.list(tempFit)){
      models[[r]][[m]]$start = NULL
      tempFit = fitFun(influenza[[r]], models[[r]][[m]])
    }
    BICbest = BIC(tempFit)
    if(is.na(BICbest)) BICbest = Inf
    foundBest = F
    
    while(!foundBest){
      # Check increase of Epidemic component
      if(S[1] + 1 <= 16){
        models[[r]][[m]]$ne$f = addSeason2formula(~-1 + fe(1, unitSpecific = TRUE), period = 52.1775, S = S[1] + 1)
        models[[r]][[m]]$end$f = addSeason2formula(~-1 + fe(1, unitSpecific = TRUE), period = 52.1775, S = S[2])
        tempFit = fitFun(influenza[[r]], models[[r]][[m]])
        if(!is.list(tempFit)){
          models[[r]][[m]]$start = NULL
          tempFit = fitFun(influenza[[r]], models[[r]][[m]])
        }
        bicEpi = BIC(tempFit)
        if(is.na(bicEpi)){
          bicEpi = Inf
          print(paste0("No convergence at ", S[1] + 1, " and ", S[2]))
        }
      }else{
        bicEpi = Inf
      }
      
      if(S[2] + 1 <= 16){
        models[[r]][[m]]$ne$f = addSeason2formula(~-1 + fe(1, unitSpecific = TRUE), period = 52.1775, S = S[1])
        models[[r]][[m]]$end$f = addSeason2formula(~-1 + fe(1, unitSpecific = TRUE), period = 52.1775, S = S[2] + 1)
        tempFit = fitFun(influenza[[r]], models[[r]][[m]])
        if(!is.list(tempFit)){
          models[[r]][[m]]$start = NULL
          tempFit = fitFun(influenza[[r]], models[[r]][[m]])
        }
        bicEnd = BIC(tempFit)
        if(is.na(bicEnd)){
          bicEnd = Inf
          print(paste0("No convergence at ", S[1], " and ", S[2] + 1))
        }
      }else{
        bicEnd = Inf
      }
      
      if((bicEnd < BICbest) | (bicEpi < BICbest)){
        if(bicEnd < bicEpi){
          print(paste0("Region: ", r, " Model: ", m, " Increase Endemic to ", S[2] + 1))
          print(bicEnd)
          S[2] = S[2] + 1
          BICbest = bicEnd
        }else{
          print(paste0("Region: ", r, " Model: ", m, " Increase Epidemic to ", S[1] + 1))
          print(bicEpi)
          S[1] = S[1] + 1
          BICbest = bicEpi
        }
      }else{
        print(paste0("Region: ", r, " Model: ", m, " Best Model: ", S[1], " ", S[2]))
        Smat[r, m, ] = S
        foundBest = T
      }
    }
  }
}

Smat[, , 1]
Smat[, , 2]

# Using seasonality order and using full data set
for(r in 1:4){
  for(m in 1:11){
    models[[r]][[m]]$ne$f = addSeason2formula(~-1 + fe(1, unitSpecific = TRUE), period = 52.1775, S = Smat[r, m, 1])
    models[[r]][[m]]$end$f = addSeason2formula(~-1 + fe(1, unitSpecific = TRUE), period = 52.1775, S = Smat[r, m, 2])
    models[[r]][[m]]$subset = setdiff(2:1036, 434:485)
    # models[[r]][[m]]$subset = setdiff(2:928, 434:485)
  }
}

for(r in 1:4){
  for(i in 1:11){
    fits[[r]][[i]] = fitFun(influenza[[r]], models[[r]][[i]])
    print(c(i, r))
  }
}

unlist(lapply(fits$N, BIC))
unlist(lapply(fits$E, BIC))
unlist(lapply(fits$S, BIC))
unlist(lapply(fits$N, BIC))


round(sapply(1:4, function(r) unlist(lapply(fits[[r]], BIC))))

modelNames = c("None", "NGB", "FG-B", "FGs-B", "Centroid-B", "FG-T", "FG-G", "FGs-G", "Centroid-G", "FG-C", "FG-R")

tab = paste0(modelNames, " & ", paste0("(", Smat[1, , 2], ",", Smat[1, , 1], ")"))
tab = paste0(tab, " & ", round(unlist(lapply(fits[[1]], BIC))))
tab = paste0(tab, " & ", paste0("(", Smat[2, , 2], ",", Smat[2, , 1], ")"))
tab = paste0(tab, " & ", round(unlist(lapply(fits[[2]], BIC))))
tab = paste0(tab, " & ", paste0("(", Smat[3, , 2], ",", Smat[3, , 1], ")"))
tab = paste0(tab, " & ", round(unlist(lapply(fits[[3]], BIC))))
tab = paste0(tab, " & ", paste0("(", Smat[4, , 2], ",", Smat[4, , 1], ")"))
tab = paste0(tab, " & ", round(unlist(lapply(fits[[4]], BIC))))

paste0(tab, collapse = "\\")

plotBIC = function(r, legendPos = "bottomright"){
  bics = unlist(lapply(fits[[r]], BIC))
  decrease = bics[1] - bics[2:11]
  
  plot(decrease[c(1, 2, 5, 6, 9, 10)], xaxt = "n", ylim = c(min(decrease), max(decrease)),
       pch = c(17, rep(16, 5)), cex = 3, xlab = "Model", ylab = "BIC Decrease - Relative to None",
       main = c("North", "East", "South", "West")[r], col = alpha(cols[c(2, 1, 3, 4)[r]], 0.5),
       xlim = c(0.5, 6.5))
  points(decrease[c(1, 2, 5, 6, 9, 10)], pch = c(2, rep(1, 5)), cex = 3, col = cols[c(2, 1, 3, 4)[r]], lwd = 2)
  axis(1, 1:6, labels = c("NGB", "Beeline", "Time", "Gravity", "Circle", "Radiation"))
  points(c(2, 2), decrease[c(3, 4)], pch = c(17, 15), cex = 3, col = alpha(cols[c(2, 1, 3, 4)[r]], 0.5))
  points(c(4, 4), decrease[c(7, 8)], pch = c(17, 15), cex = 3, col = alpha(cols[c(2, 1, 3, 4)[r]], 0.5))
  points(c(2, 2), decrease[c(3, 4)], pch = c(2, 0), cex = 3, col = cols[c(2, 1, 3, 4)[r]], lwd = 2)
  points(c(4, 4), decrease[c(7, 8)], pch = c(2, 0), cex = 3, col = cols[c(2, 1, 3, 4)[r]], lwd = 2)
  orders = paste0("(", Smat[1, 2:11, 2], ",", Smat[1, 2:11, 1], ")")
  text(c(1, 2, 2, 2, 3, 4, 4, 4, 5, 6) + 0.4, decrease, label = orders)
  legend(legendPos, pch = c(16, 17, 15), legend = c("FG", "FGs", "Centroid"), bty = "n",
         pt.cex = 2, y.intersp = 0.3)
}

par(mfrow = c(2, 2))
plotBIC(1, "topleft")
plotBIC(2)
plotBIC(3, "bottomleft")
plotBIC(4)
par(mfrow = c(1, 1))

getWpars = function(fitt){
  temp = fitt$coefficients
  temp[startsWith(names(temp), "neweights")]
}

getSE = function(fitt){
  temp = fitt$se
  temp[startsWith(names(temp), "neweights")]
}

decay = matrix(0, ncol = 8, nrow = 13)
u = qnorm(0.975)

decay[, 1] = unlist(lapply(fits$N, getWpars))[c(1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 14, 15)]
decay[, 2] = u*unlist(lapply(fits$N, getSE))[c(1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 14, 15)]

decay[, 3] = unlist(lapply(fits$E, getWpars))[c(1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 14, 15)]
decay[, 4] = u*unlist(lapply(fits$E, getSE))[c(1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 14, 15)]

decay[, 5] = unlist(lapply(fits$S, getWpars))[c(1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 14, 15)]
decay[, 6] = u*unlist(lapply(fits$S, getSE))[c(1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 14, 15)]

decay[, 7] = unlist(lapply(fits$W, getWpars))[c(1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 14, 15)]
decay[, 8] = u*unlist(lapply(fits$W, getSE))[c(1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 14, 15)]

decay = paste(round(decay, digits = 3))
decay = matrix(decay, nrow = 13)
decay[, 2] = paste0("$\\pm$ ", decay[, 2])
decay[, 4] = paste0("$\\pm$ ", decay[, 4])
decay[, 6] = paste0("$\\pm$ ", decay[, 6])
decay[, 8] = paste0("$\\pm$ ", decay[, 8])

paste0(apply(decay, 1, function(x) paste0(x, collapse = " & ")), collapse = "\\")

########################################
##### Stage 2: One season forecast #####
########################################

predicthhh4 = function(fit, h, nChains, p, true, addCase = c(0, 1)){
  temp = fit$control$ne$f
  temp2 = strsplit(as.character(temp[2]), "cos")[[1]]
  temp3 = strsplit(temp2[length(temp2)], " * ")[[1]][1]
  Sepi = as.numeric(paste0(strsplit(temp3, "")[[1]][-1], collapse = ""))/2
  temp = fit$control$end$f
  temp2 = strsplit(as.character(temp[2]), "cos")[[1]]
  temp3 = strsplit(temp2[length(temp2)], " * ")[[1]][1]
  Send = as.numeric(paste0(strsplit(temp3, "")[[1]][-1], collapse = ""))/2
  
  n = ncol(fit$fitted.values)
  beta = coefficients(fit)
  # Structure parameters
  betaLambda = beta[1:(Sepi*2)]
  alphaLambda = beta[(Sepi*2 + 1):(Sepi*2 + n)]
  betaNu = beta[(Sepi*2 + n + 1):(Sepi*2 + n + Send*2)]
  alphaNu = beta[(Sepi*2 + n + Send*2 + 1):(Sepi*2 + 2*n + Send*2)]
  # psi = 1/exp(-tail(beta, ifelse(fit$control$family == "NegBin1", 1, n)))
  psi = 1/tail(beta, ifelse(fit$control$family == "NegBin1", 1, n))
  
  # Get weight matrix, transposing to make computations below easier!
  W = t(getWeights(fit))
  # Find last observations
  iLast = tail(fit$control$subset, 1)
  Ylast = observed(fit$stsObj)[iLast, ]
  Ylast[addCase[2]] = Ylast[addCase[2]] + addCase[1]
  Y = array(0, dim = c(n, h+1, nChains))
  for(i in 1:nChains) Y[, 1, i] = Ylast
  # Prepare log-linear predictors
  tSeq = iLast:(iLast + h)
  zLambda = matrix(0, nrow = h+1, ncol = Sepi*2)
  zNu = matrix(0, nrow = h+1, ncol = Send*2)
  
  for(i in 1:Sepi){
    zLambda[, 2*i - 1] = sin(2*i*pi*tSeq/52)
    zLambda[, 2*i] = cos(2*i*pi*tSeq/52)
  }
  for(i in 1:Send){
    zNu[, 2*i - 1] = sin(2*i*pi*tSeq/52)
    zNu[, 2*i] = cos(2*i*pi*tSeq/52)
  }

  seasLambda = (zLambda%*%betaLambda)[, 1]
  seasNu = (zNu%*%betaNu)[, 1]
  
  lambda = exp(outer(alphaLambda, seasLambda, "+"))
  nu = exp(outer(alphaNu, seasNu, "+"))
  
  for(i in 1:nChains){
    for(hh in 1:h){
      Wtemp = (W%*%Y[, hh, i])[, 1]
      mu = nu[, hh+1] + lambda[, hh+1] * Wtemp
      Y[, hh+1, i] = rnbinom(n, size = psi, mu = mu)
    }
  }
  
  # Prepare output
  alphaSeq = 2*p[p < 0.5]
  
  m = matrix(0, nrow = n + 1, ncol = h)
  Q = array(0, dim = c(n + 1, h, length(p)))
  WIS = matrix(0, nrow = n + 1, ncol = h)
  WIS_D = matrix(0, nrow = n + 1, ncol = h)
  WIS_U = matrix(0, nrow = n + 1, ncol = h)
  WIS_O = matrix(0, nrow = n + 1, ncol = h)
  LIS = matrix(0, nrow = n + 1, ncol = h)
  LIS_D = matrix(0, nrow = n + 1, ncol = h)
  LIS_U = matrix(0, nrow = n + 1, ncol = h)
  LIS_O = matrix(0, nrow = n + 1, ncol = h)
  
  RPS = matrix(0, nrow = n + 1, ncol = h)
  
  M = max(true)
  
  for(hh in 1:h){
    for(i in 1:n){
      y = Y[i, hh+1, ]
      m[i, hh] = mean(y)
      Q[i, hh, ] = sapply(p, function(pp) quantile(y, pp))
      WIS[i, hh] = abs(quantile(y, 0.5) - true[hh, i])/2
      for(j in 1:length(alphaSeq)){
        a = alphaSeq[j]
        l = quantile(y, a/2)
        u = quantile(y, 1 - a/2)
        tt = true[hh, i]
        temp = scoringutils:::interval_score(tt, l, u, (1 - a)*100, T, T)
        WIS[i, hh] = WIS[i, hh] + temp$interval_score
        WIS_O[i, hh] = WIS_O[i, hh] + temp$overprediction
        WIS_U[i, hh] = WIS_U[i, hh] + temp$underprediction
        WIS_D[i, hh] = WIS_D[i, hh] + temp$dispersion
        # WIS[i, hh] = WIS[i, hh] + ((u - l) + 2/a*(l - tt)*(tt < l) + 2/a*(tt - u)*(tt > u))*wSeq[j]
        
        l = log(l + 1)
        u = log(u + 1)
        tt = log(tt + 1)
        temp = scoringutils:::interval_score(tt, l, u, (1 - a)*100, T, T)
        LIS[i, hh] = LIS[i, hh] + temp$interval_score
        LIS_O[i, hh] = LIS_O[i, hh] + temp$overprediction
        LIS_U[i, hh] = LIS_U[i, hh] + temp$underprediction
        LIS_D[i, hh] = LIS_D[i, hh] + temp$dispersion
        # LIS[i, hh] = LIS[i, hh] + (u - l) + 2/a*(l - tt)*(tt < l) + 2/a*(tt - u)*(tt > u)
      }

      RPS[i, hh] = scoringutils::crps_sample(observed = c(true[hh, i]), predicted = matrix(y, nrow = 1))
    }
  }
  
  return(list(mean = m, quantiles = Q, WIS = WIS, WIS_D = WIS_D, WIS_U = WIS_U, WIS_O = WIS_O,
              LIS = LIS, LIS_D = LIS_D, LIS_U = LIS_U, LIS_O = LIS_O, RPS = RPS))
}


# Let's look at the number of components here
bar_df = data.frame(region = character(100),
                    measure = character(100),
                    nComponent = numeric(100),
                    value = numeric(100))
bar_df[, 1:3] = expand.grid(c("North", "East", "South", "West"), c("B", "T", "G", "C", "R"), 1:5)

addVal = 1/(length(indN)^2/2 + length(indN)/2)
for(i in 1:length(indN)){
  for(j in i:length(indN)){
    nc = sum(parsB_N$w[i, j, ] != 0)
    indTemp = (bar_df$region == "North") & (bar_df$measure == "B") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
    
    nc = sum(parsT_N$w[i, j, ] != 0)
    indTemp = (bar_df$region == "North") & (bar_df$measure == "T") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
    
    nc = sum(parsG_N$w[i, j, ] != 0)
    indTemp = (bar_df$region == "North") & (bar_df$measure == "G") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
    
    nc = sum(parsC_N$w[i, j, ] != 0)
    indTemp = (bar_df$region == "North") & (bar_df$measure == "C") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
    
    nc = sum(parsR_N$w[i, j, ] != 0)
    indTemp = (bar_df$region == "North") & (bar_df$measure == "R") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
  }
}

addVal = 1/(length(indE)^2/2 + length(indE)/2)
for(i in 1:length(indE)){
  for(j in i:length(indE)){
    nc = sum(parsB_E$w[i, j, ] != 0)
    indTemp = (bar_df$region == "East") & (bar_df$measure == "B") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
    
    nc = sum(parsT_E$w[i, j, ] != 0)
    indTemp = (bar_df$region == "East") & (bar_df$measure == "T") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
    
    nc = sum(parsG_E$w[i, j, ] != 0)
    indTemp = (bar_df$region == "East") & (bar_df$measure == "G") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
    
    nc = sum(parsC_E$w[i, j, ] != 0)
    indTemp = (bar_df$region == "East") & (bar_df$measure == "C") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
    
    nc = sum(parsR_E$w[i, j, ] != 0)
    indTemp = (bar_df$region == "East") & (bar_df$measure == "R") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
  }
}

addVal = 1/(length(indS)^2/2 + length(indS)/2)
for(i in 1:length(indS)){
  for(j in i:length(indS)){
    nc = sum(parsB_S$w[i, j, ] != 0)
    indTemp = (bar_df$region == "South") & (bar_df$measure == "B") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
    
    nc = sum(parsT_S$w[i, j, ] != 0)
    indTemp = (bar_df$region == "South") & (bar_df$measure == "T") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
    
    nc = sum(parsG_S$w[i, j, ] != 0)
    indTemp = (bar_df$region == "South") & (bar_df$measure == "G") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
    
    nc = sum(parsC_S$w[i, j, ] != 0)
    indTemp = (bar_df$region == "South") & (bar_df$measure == "C") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
    
    nc = sum(parsR_S$w[i, j, ] != 0)
    indTemp = (bar_df$region == "South") & (bar_df$measure == "R") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
  }
}

addVal = 1/(length(indW)^2/2 + length(indW)/2)
for(i in 1:length(indW)){
  for(j in i:length(indW)){
    nc = sum(parsB_W$w[i, j, ] != 0)
    indTemp = (bar_df$region == "West") & (bar_df$measure == "B") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
    
    nc = sum(parsT_W$w[i, j, ] != 0)
    indTemp = (bar_df$region == "West") & (bar_df$measure == "T") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
    
    nc = sum(parsG_W$w[i, j, ] != 0)
    indTemp = (bar_df$region == "West") & (bar_df$measure == "G") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
    
    nc = sum(parsC_W$w[i, j, ] != 0)
    indTemp = (bar_df$region == "West") & (bar_df$measure == "C") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
    
    nc = sum(parsR_W$w[i, j, ] != 0)
    indTemp = (bar_df$region == "West") & (bar_df$measure == "R") & (bar_df$nComponent == nc)
    bar_df$value[indTemp] = bar_df$value[indTemp] + addVal
  }
}

levels(bar_df$measure) = c("Beeline", "Time", "Gravity", "Circle", "Radiation")

colsTemp = cols[c(2, 1, 3, 4)]

p = ggplot(bar_df, aes(x = nComponent, y = value, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +  
  facet_wrap(~ measure, scales = "fixed") +                              
  scale_fill_manual(values = colsTemp) +                                 
  labs(x = "Number of components", y = "Value", fill = "Region") +
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.15),                                     
    legend.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank(),                               
    strip.text = element_text(size = 12),                                
    axis.text.x = element_text(angle = 45, hjust = 1)                    
  )
p
ggsave("Plots/Components.png", plot = p, width = 8, height = 6, dpi = 300)


sapply(unique(bar_df$measure), function(measure){
  sum(bar_df[bar_df$measure == measure, ]$nComponent * bar_df[bar_df$measure == measure, ]$value)/4
})




###########################################
##### Stage 3: 8-weeks-ahead forecast #####
###########################################

# Start with 2:928 and predict 8 weeks
# Then 3:929, ...

h = 8
tShiftMax = 100
p = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)
n = length(shp)

# That should run over night?
forecast2 = matrix(0, nrow = n*(tShiftMax + 1)*(length(p) + 10)*h*length(models$N), ncol = 7)

counter = 1
t1 = Sys.time()
for(r in 1:4){
  n = length(ind[[r]])
  for(tShift in 0:tShiftMax){
    # Check if the block has not been computed yet
    checkInd = counter + length(models$N)*n*h*(length(p) + 10) - 1
    if(forecast2[checkInd, 1] == 0){
      for(mm in 1:length(models$N)){
        m = mm
        modelTemp = models[[r]][[m]]
        modelTemp$subset = setdiff((2:928) + tShift,  434:485)
        modelTemp$start = as.list(coefficients(fits[[r]][[m]]))
        if(mm > 1){
          indTemp = which(startsWith(names(fits[[r]][[m]]$coefficients), "neweights"))
          modelTemp$ne$weights$initial = fits[[r]][[m]]$coefficients[indTemp]
        }
        
        
        fitTemp = tryCatch({
          fitFun(influenza[[r]], modelTemp)
        }, error = function(e){
          modelTemp$start = NULL
          if(mm > 1){
            modelTemp$ne$weights$initial = modelTemp$ne$weights$initial / 2
          }
          fitFun(influenza[[r]], modelTemp)
        })
        
        fitTemp$control$subset = modelTemp$subset
        
        set.seed(06042024 + m + tShift * 15 + r)
        temp = predicthhh4(fitTemp, h, nChains = 1000, p = p, cases[(929 + tShift):(936 + tShift), ind[[r]]])
        for(i in 1:n){
          for(hh in 1:h){
            forecast2[counter, ] = c(i, 928 + tShift, hh, 0, mm, temp$mean[i, hh], r)
            counter = counter + 1
            for(pp in 1:length(p)){
              forecast2[counter, ] = c(i, 928 + tShift, hh, p[pp], mm, temp$quantiles[i, hh, pp], r)
              counter = counter + 1
            }
            forecast2[counter, ] = c(i, 928 + tShift, hh, 2, mm, temp$WIS[i, hh], r)
            counter = counter + 1
            forecast2[counter, ] = c(i, 928 + tShift, hh, 3, mm, temp$WIS_D[i, hh], r)
            counter = counter + 1
            forecast2[counter, ] = c(i, 928 + tShift, hh, 4, mm, temp$WIS_U[i, hh], r)
            counter = counter + 1
            forecast2[counter, ] = c(i, 928 + tShift, hh, 5, mm, temp$WIS_O[i, hh], r)
            counter = counter + 1
            
            forecast2[counter, ] = c(i, 928 + tShift, hh, 6, mm, temp$LIS[i, hh], r)
            counter = counter + 1
            forecast2[counter, ] = c(i, 928 + tShift, hh, 7, mm, temp$LIS_D[i, hh], r)
            counter = counter + 1
            forecast2[counter, ] = c(i, 928 + tShift, hh, 8, mm, temp$LIS_U[i, hh], r)
            counter = counter + 1
            forecast2[counter, ] = c(i, 928 + tShift, hh, 9, mm, temp$LIS_O[i, hh], r)
            counter = counter + 1
            
            forecast2[counter, ] = c(i, 928 + tShift, hh, 10, mm, temp$RPS[i, hh], r)
            counter = counter + 1
          }
        }
        print(c(m, tShift))
      }
      print(paste("Estimated finish:", t1 + (Sys.time() - t1)/counter*nrow(forecast2)))
    }else{
      # If it has already been computed, increase counter
      counter = checkInd + 1
    }
  }
}


# Fix weird scoringutils update
forecast2 = as.data.frame(forecast2)
names(forecast2) = c("district", "tEnd", "h", "type", "model", "value", "region")

# write.csv(x = forecast2, file = "forecast2.csv")
# forecast2 = read.csv("forecast2.csv", header = T)


#####################
### Check results ###
#####################

# Find the WIS, lWIS and RPS for 1 week ahead forecasts
# Split by region and model
forecast1wk = forecast2[forecast2$h == 1, ]
forecast1wk = forecast1wk[forecast1wk$type %in% c(2, 6, 10), ]

tempMat = matrix(nrow = 12, ncol = 11)
counter = 1
for(r in 1:4){
  for(tt in c(2, 6, 10)){
    tempMat[counter, ] = sapply(1:11, function(m){
      temp = (forecast1wk$region == r) & (forecast1wk$type == tt) & (forecast1wk$model == m)
      sum(forecast1wk$value[temp])/length(ind[[r]])
    })
    counter = counter + 1
  }
}

# Fix weird weighting
tempMat[c(1, 2, 4, 5, 7, 8, 10, 11), ] = tempMat[c(1, 2, 4, 5, 7, 8, 10, 11), ]/2
tempMat = t(tempMat)

relWIS1 = matrix(0, nrow = 10, ncol = 12)
for(i in 2:11){
  relWIS1[i-1, ] = tempMat[i, ]/tempMat[1, ]
}

paste0(apply(round(relWIS1, digits = 3), 1, function(x){
  paste0(x, collapse = " & ")
}), collapse = "\\")

tempMat2 = matrix(nrow = 10, ncol = 12)

for(cc in 1:12){
  theta = matrix(0, nrow = 11, ncol = 11)
  for(i in 1:11) for(j in 1:11) theta[i, j] = tempMat[i, cc] / tempMat[j, cc]
  theta = apply(theta, 1, function(x) prod(x)^(1/10))
  tempMat2[, cc] = theta[2:11]/theta[1]
}

round(tempMat2*100, digits = 1)
paste0(apply(round(tempMat2*100, digits = 1), 1, function(x){
  paste0(x, collapse = " & ")
}), collapse = "\\")


# For regions 1-4, horizons c(1, 4, 8)
forecast148 = forecast2[forecast2$h %in% c(1, 4, 8) & (forecast2$type == 6), ]

tempMat2 = matrix("", nrow = 33, ncol = 8)
for(h in 1:3){
  for(r in 1:4){
    temp = sapply(1:11, function(m){
      temp = (forecast148$region == r) & (forecast148$h == c(1, 4, 8)[h]) & (forecast148$model == m)
      sum(forecast148$value[temp])/length(ind[[r]])/2
    })
    
    nms = c("None", "NGB", "FG-B", "FGs-B", "Cent-B", "FG-T", "FG-G", "FGs-G", "Cent-G", "FG-C", "FG-R")
    tempS = sort(temp)
    tempMat2[(h*11 - 10):(h*11), (r*2 - 1):(r*2)] = cbind(nms[order(temp)], round(tempS, digits = 2))
  }
}

paste0(apply(tempMat2, 1, function(x){
  paste0(x, collapse = " & ")
}), collapse = "\\")


# Risk mapping on a fine grid

# Select the week of 2019 with the most cases: Index 947
current = cases[947, indS]
modelTemp = models$S[[3]]
modelTemp$subset = modelTemp$subset[modelTemp$subset <= 947]
fittedTemp = fitFun(infl = influenza$S, model = modelTemp)
dEst = tail(fittedTemp$coefficients, 2)[1]

# Set wd to population data folder 
# Download population data
library(raster)
library(geosphere)
pop = getPOP("DEU", 2011, 1000) # Needs PopBoundaryAPI.R from other GitHub repo
# set wd back to project folder

popS = mask(crop(pop, as_Spatial(shp[indS])), as_Spatial(shp[indS]))
popS10 = aggregate(x = popS, fact = 10, FUN = sum)
riskRaster = popS10
popSmat = as.matrix(popS10)
# dim(popSmat)
# plot(popS10)

risks = matrix(0, nrow = prod(dim(popSmat)), ncol = length(current))
risksW = matrix(0, nrow = prod(dim(popSmat)), ncol = length(current))
coord = coordinates(popS10)

# Find the risk map of one (!) infected person in a certain district.
# From there, aggregate the risks according to the observed counts!

for(i in 1:length(current)){
  popTemp = mask(popS10, as_Spatial(shp[indS[i]]))
  indDistrict = !is.na(values(popTemp)) & (values(popTemp) > 0)
  coordOrigin = coord[indDistrict, ]
  popOrigin = values(popS10)[indDistrict]
  
  for(j in 1:nrow(coordOrigin)){
    distTemp = (distGeo(coordOrigin[j, 1:2], coord) + 500)/1000
    risk = distTemp^(-dEst)
    risks[, i] = risks[, i] + risk
    risksW[, i] = risksW[, i] + risk * popOrigin[j]
  }
  print(i)
}


plot(rowSums(cases[940:950, indS]))
current = cases[947, indS]

# t = 8 (week of the year)
beta = fits$S[[4]]$coefficients
beta1 = beta[startsWith(names(beta), "end.sin")]
beta2 = beta[startsWith(names(beta), "end.cos")]
beta3 = beta[startsWith(names(beta), "ne.sin")]
beta4 = beta[startsWith(names(beta), "ne.cos")]

endFE = beta[startsWith(names(beta), "end.1")]
neFE = beta[startsWith(names(beta), "ne.1")]

# Endemic component at district level:
nu = exp(sum(beta1*sin((1:length(beta1))*2*pi*9/52.1775)) + sum(beta2*cos((1:length(beta2))*2*pi*9/52.1775)) + endFE)

# Log-linear predictor for epidemic compoenent
lambda = exp(sum(beta3*sin((1:length(beta3))*2*pi*9/52.1775)) + sum(beta4*cos((1:length(beta4))*2*pi*9/52.1775)) + neFE)

# Split endemic component to grid cell level
# weighted by population density

riskEndemic = rep(0, length(popS10))
for(i in 1:length(current)){
  temp = values(mask(popS10, as_Spatial(shp[indS[i]])))
  temp[is.na(temp)] = 0
  temp = temp/sum(temp)
  riskEndemic = riskEndemic + nu[i] * temp
}
riskEndemic = riskEndemic * (values(popS10) > -1)

riskEndemicRaster = popS10
values(riskEndemicRaster) = riskEndemic
plot(riskEndemicRaster)

# Next epidemic component


# From risksW split to grid cells
# risksW: For each grid cell and each district of origin
#         gives the contact intensity of a randomly selected person
#         from origin to grid cell
#         Normalised instead of row normalisation
risksWnormalised = risksW
for(i in 1:140){
  risksWnormalised[, i] = risksW[, i]/sum(risksW[, 1])
}
# risksWnormalised gives contact probability of a randomly selected infected
# with someone from grid cell

riskEpidemic = rep(0, length(values(popS10)))

# Collect all transmission risks from each district to each grid cell
# For the grid cells, important what the district effect is!
# Either iterate through all districts, or create a vector of effects?

neFElong = rep(0, length(values(popS10)))
for(i in 1:140){
  indDistrict = which(!is.na(values(mask(popS10, as_Spatial(shp[indS[i]])))))
  neFElong[indDistrict] = lambda[i]
}

for(i in 1:140){
  riskEpidemic = riskEpidemic + risksWnormalised[, i] * current[i]
}


riskEpidemic[is.na(values(popS10))] = NA
riskEpidemicRaster = popS10
values(riskEpidemicRaster) = log(riskEpidemic * neFElong + riskEndemic)

riskEpidemicRaster1 = riskEpidemicRaster
values(riskEpidemicRaster1) = log(exp(values(riskEpidemicRaster)) / values(popS10))

temp = values(riskEpidemicRaster1)
indTemp = which(!is.finite(temp) & (!is.na(temp)))
indTemp2 = which(is.finite(temp))
temp[indTemp] = min(temp[indTemp2])

values(riskEpidemicRaster1) = temp

par(mfrow = c(1, 2))
plot(riskEpidemicRaster)
plot(riskEpidemicRaster1)

values(riskEpidemicRaster1)
raster1Mat = as.matrix(riskEpidemicRaster1)
raster1Mat[is.na(raster1Mat)] = 0
raster1MatSmooth = raster1Mat
for(i in 2:(nrow(raster1MatSmooth) - 1)){
  for(j in 2:(ncol(raster1MatSmooth) - 1)){
    raster1MatSmooth[i, j] = mean(raster1Mat[(i-1):(i+1), (j-1):(j+1)])
  }
}
riskEpidemicRaster1Smooth = riskEpidemicRaster1
values(riskEpidemicRaster1Smooth) = c(t(raster1MatSmooth))*(values(popS10) > -100)
plot(riskEpidemicRaster1)
plot(riskEpidemicRaster1Smooth)


raster_df <- as.data.frame(rasterToPoints(riskEpidemicRaster), stringsAsFactors = FALSE)
colnames(raster_df) <- c("x", "y", "value")
raster_df$value = exp(raster_df$value)
raster1_df <- as.data.frame(rasterToPoints(riskEpidemicRaster1Smooth), stringsAsFactors = FALSE)
colnames(raster1_df) <- c("x", "y", "value") 
raster1_df$value = exp(raster1_df$value)
raster1_df$value = raster1_df$value/100


library(latex2exp)

ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = value)) +
  geom_sf(data = shp[indS, ], fill = NA, color = "black", size = 0.5) +
  scale_fill_gradientn(
    name = TeX("Expected\nCases per $km^2$"),                # Legend title
    trans = "log",                            # Log scale for colors
    colors = c("lightgreen", "darkgreen", "orange", "darkred"),  # Define your gradient
    breaks = c(0.0001, 0.001, 0.01, 0.1),     # Example original values
    labels = c(0.0001, 0.001, 0.01, 0.1)    # Corresponding labels
  ) +
  theme_minimal() +
  theme(
    legend.background = element_rect(fill = NA, color = NA), # White background for legend
    panel.grid = element_blank(),                                # Remove all grid lines
    strip.text = element_text(size = 12),                        # Customize facet label size
    axis.text = element_blank(),            # Rotate x-axis labels for readability
    axis.title = element_blank(),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  )
ggsave("Plots/RiskMap.png", width = 8, height = 6, dpi = 350)

ggplot() +
  geom_raster(data = raster1_df, aes(x = x, y = y, fill = value)) +
  geom_sf(data = shp[indS, ], fill = NA, color = "black", size = 0.5) +
  scale_fill_gradientn(
    name = "Individual\nRisk",                # Legend title
    trans = "log",                            # Log scale for colors
    colors = c("lightgreen", "darkgreen", "orange", "darkred"),  # Define your gradient
    breaks = c(0.0001, 0.001, 0.01),     # Example original values
    labels = c("1/10000", "1/1000", "1/100")    # Corresponding labels
  ) +
  theme_minimal() +
  theme(
    legend.background = element_rect(fill = NA, color = NA), # White background for legend
    panel.grid = element_blank(),                                # Remove all grid lines
    strip.text = element_text(size = 12),                        # Customize facet label size
    axis.text = element_blank(),                                 # Remove axis text
    axis.title = element_blank(),                                 # Remove axis titles
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  )
ggsave("Plots/RiskMapIndiv.png", width = 8, height = 6, dpi = 350)
