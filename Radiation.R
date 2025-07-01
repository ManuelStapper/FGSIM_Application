rm(list = ls())

library("sf")
library("sp")
library("geosphere")
library("Rcpp")
library("rnaturalearth")
library("units")

# Setting working directory here

source("wspsample.R")
source("PopBoundaryAPI.R")
source("ParameterPrep/Fitting.R")
# sourceCpp("ParameterPrep/Radiation.cpp")
sourceCpp("ParameterPrep/Radiation2.cpp")

# Needs a table IDlookup.csv mapping shapefile and IDs
ID2Geom = read.csv("Data/IDlookup.csv", header = T, fileEncoding = "UTF-8")
IDs = ID2Geom$ID

countryCode = "DEU"
level = 3
year = 2011

# Download data
# Setting working directory for data here

Ger = getJSON(countryCode, level, timeout = 1000)
GerS = as_Spatial(Ger$geometry[ID2Geom$JsonNr])


# Download population data
pop = getPOP(countryCode, year, 1000)
# Setting working directory back to project folder

samO = list()
samD = list()

for(i in 1:401){
  set.seed(i)
  samO[[i]] = wspsample(GerS[i], 2000, pop, weighted = T)
  set.seed(401 + i)
  samD[[i]] = wspsample(GerS[i], 2000, pop, weighted = T)
}


indN = which(floor(IDs/1000) %in% c(1, 2, 3, 4))
indE = which(floor(IDs/1000) %in% c(11, 12, 13, 14, 15, 16))
indS = which(floor(IDs/1000) %in% c(8, 9))
indW = which(floor(IDs/1000) %in% c(5, 6, 7, 10))
ind = list(indN, indE, indS, indW)


# We need four population data sets for the four regions
# Ideally, they are of suitable size (covering the observation region, but not too large)


# Maximum centroid distances:
# (NESW): 360, 443, 443, 367

# Find min and max coordinates and check with countries overlap
# N: (0.5, 17.8) & (47.5, 58.8)
# E: (2.6, 22.4) & (45.6, 59.3)
# S: (0.7, 20.9) & (42.6, 55.2)
# W: (0.0, 16.2) & (45.2, 56.3)

regions = sapply(ind, function(xx) st_union(Ger$geometry[ID2Geom$JsonNr][xx]))
bufferRegion = lapply(ind, function(xx){
  region = st_union(Ger$geometry[ID2Geom$JsonNr][xx])
  st_buffer(region, dist = set_units(500, km))
})

world = ne_countries(scale = "medium", returnclass = "sf")
europe = world[world$continent == "Europe", ]

# plot(europe$geometry, xlim = c(0, 22.4), ylim = c(42.6, 59.3))
# plot(bufferRegion[[1]], col = rgb(1, 0, 0, alpha = 0.2), add = T)
# plot(bufferRegion[[2]], col = rgb(0, 1, 0, alpha = 0.2), add = T)
# plot(bufferRegion[[3]], col = rgb(0, 0, 1, alpha = 0.2), add = T)
# plot(bufferRegion[[4]], col = rgb(1, 0, 1, alpha = 0.2), add = T)

GerAll = st_union(Ger$geometry[unlist(ind)])
GerBuffer = st_buffer(GerAll, dist = set_units(500, km))

europe$name

# plot(europe$geometry, xlim = c(0, 22.4), ylim = c(42.6, 59.3))
# plot(GerBuffer, add = T, col = rgb(1, 0, 0, alpha = 0.2))

# Load worldpop data and save everything that is in the MBR around the buffers
ngbCodes = c("GBR", "CHE", "SWE", "SVK", "SVN", "POL", "NOR", "NLD", "LUX",
             "LIE", "ITA", "HUN", "FRA", "DNK", "CZE", "BIH", "BEL", "AUT", "HRV")
# Setting working directory to population data folder
popNGB = list()
for(i in 1:length(ngbCodes)){
  temp = getPOP(ngbCodes[i], year, 1000)
  temp = crop(temp, as(GerBuffer, "Spatial"))
  popNGB[[i]] = temp
  print(i)
}
# Setting working directory back to project folder

popNGB[[20]] = pop
popAll = do.call(merge, popNGB)
popAll = mask(popAll, as(GerBuffer, "Spatial"))

popAllLog = aggregate(popAll, fact = 10, sum)
values(popAllLog) = log(values(popAllLog))

# library(ggspatial)
# popAllLog_df <- as.data.frame(rasterToPoints(popAllLog), xy = TRUE)
# p = ggplot() +
#   geom_raster(data = popAllLog_df, aes(x = x, y = y, fill = layer)) +
#   geom_sf(data = europe$geometry, fill = "transparent", color = "black", size = 0.5) +
#   geom_sf(data = GerBuffer, fill = "transparent", color = "black") +
#   
#   theme_minimal() +
#   theme(
#     panel.grid.major = element_line(color = "grey90"),   # Light grid lines
#     axis.title = element_blank(),                        # Remove axis titles
#     axis.text = element_blank(),                         # Remove axis text
#     axis.ticks = element_blank(),                        # Remove axis ticks
#     plot.background = element_rect(fill = "white", color = NA), # Clean white background
#     legend.position = "right",                           # Adjust legend placement
#     legend.title = element_text(size = 10),
#     legend.text = element_text(size = 8)
#   ) +
#   
#   scale_fill_viridis_c(name = "Log(Population)") +
#   coord_sf(xlim = c(-5, 27.41), ylim = c(41, 61), expand = FALSE)
# ggsave("Plots/PopDensity.png", plot = p, width = 6, height = 8, dpi = 300)

popMat = as.matrix(popAll)
pop10 = aggregate(popAll, fact = 10, fun = sum)
popMat10 = as.matrix(pop10)

for(i in 1:ncol(popMat)){
  popMat[is.na(popMat[, i]), i] = 0
}
popMat10[is.na(popMat10)] = 0

replaceIndices = function(coords, rangeLat, rangeLong, dims){
  i = dims[1] - floor((coords[, 2] - rangeLong[1])/(rangeLong[2] - rangeLong[1])*dims[1])
  j = ceiling((coords[, 1] - rangeLat[1])/(rangeLat[2] - rangeLat[1])*dims[2])
  return(cbind(i, j))
}

ext1 = extent(popAll)
ext10 = extent(pop10)

for(i in 1:401){
  I1 = replaceIndices(samO[[i]][, 1:2],
                      rangeLong = c(ymin(ext1), ymax(ext1)),
                      rangeLat = c(xmin(ext1), xmax(ext1)),
                      dims = dim(popMat))
  I10 = replaceIndices(samO[[i]][, 1:2],
                       rangeLong = c(ymin(ext10), ymax(ext10)),
                       rangeLat = c(xmin(ext10), xmax(ext10)),
                       dims = dim(popMat10))
  samO[[i]] = cbind(samO[[i]], I1, I10)
  
  I1 = replaceIndices(samD[[i]][, 1:2],
                      rangeLong = c(ymin(ext1), ymax(ext1)),
                      rangeLat = c(xmin(ext1), xmax(ext1)),
                      dims = dim(popMat))
  I10 = replaceIndices(samD[[i]][, 1:2],
                       rangeLong = c(ymin(ext10), ymax(ext10)),
                       rangeLat = c(xmin(ext10), xmax(ext10)),
                       dims = dim(popMat10))
  samD[[i]] = cbind(samD[[i]], I1, I10)
}


getPop = function(origin, destination, popMat){
  n = nrow(origin)
  out = rep(0, n)
  dims = dim(popMat)
  
  for(i in 1:n){
    ind1 = origin[i, 1:2]
    ind2 = destination[i, 1:2]
    r = ceiling(sqrt(sum((ind1 - ind2)^2)))
    out[i] = computePop(ind1[1], ind1[2], r, popMat)
  }
  out
}

D = array(NA, dim = c(401, 401, 2000))

n = unlist(lapply(ind, length))
N = sum(n^2)

counter = 0
t1 = Sys.time()
for(i in 1:401){
  iRegion = which(unlist(lapply(ind, function(xx) i %in% xx)))
  for(j in ind[[iRegion]]){
    if(i == j){
      D[i, j, ] = getPop(samO[[i]][, 4:5], samD[[j]][, 4:5], popMat)
    }else{
      D[i, j, ] = getPop(samO[[i]][, 6:7], samD[[j]][, 6:7], popMat10)
    }
    counter = counter + 1
    if(counter/100 == floor(counter/100)){
      print(paste("Counter:", counter))
      estTime = t1 + (Sys.time() - t1)/counter*N
      print(paste("Estimated time:", estTime))
    }
  }
}

n = unlist(lapply(ind, length))
N = sum(n^2)
counterOuter = 1
t1 = Sys.time()
for(r in 1:4){
  indTemp = ind[[r]]
  n = length(indTemp)
  RES = matrix(nrow = n^2, ncol = 2 + 3*5 + 3)
  counter = 1
  for(i in 1:n){
    for(j in 1:n){
      x = D[indTemp[i], indTemp[j], ]/1000
      x = log(x[!is.na(x)])
      bicMin = Inf
      kMin = 0
      for(k in 1:5){
        est = estMTN(x = x, k = k, l = min(x), u = max(x), nIter = 100)
        if(!is.na(est$bic) & is.finite(est$bic)){
          if(est$bic < bicMin){
            bicMin = est$bic
            kMin = k
          }
        } else{
          break
        }
      }
      est = estMTN(x = x, k = kMin, l = min(x), u = max(x), nIter = 100)
      
      buffer = rep(0, 5 - kMin)
      RES[counter, ] = c(indTemp[i], indTemp[j], est$m, buffer, est$s, buffer + 1, est$w, buffer, est$bic, min(x), max(x))
      counter = counter + 1
      counterOuter = counterOuter + 1
      if(counterOuter/100 == floor(counterOuter/100)){
        print(paste("Counter:", counterOuter))
        estTime = t1 + (Sys.time() - t1)/counterOuter*N
        print(paste("Estimated time:", estTime))
      }
    }
  }
  RES = data.frame(RES)
  colnames(RES) = c("Origin", "Destination", paste0("mu", 1:5), paste0("sigma", 1:5), paste0("w", 1:5), "bic", "l", "u")
  
  write.csv(x = RES, file = paste("ParameterPrep/Results/Circle", c("N", "E", "S", "W")[r], ".csv", sep = ""))
}

counterOuter = 1
t1 = Sys.time()
for(r in 1:4){
  indTemp = ind[[r]]
  n = length(indTemp)
  RES = matrix(nrow = n^2, ncol = 2 + 3*5 + 3)
  counter = 1
  for(i in 1:n){
    for(j in 1:n){
      Sij = D[indTemp[i], indTemp[j], ]
      indKeep = which(!is.na(Sij))
      
      Sij = Sij[indKeep]
      Ni = samO[[indTemp[i]]][indKeep, 3]
      Nj = samD[[indTemp[j]]][indKeep, 3]
      x = log((Ni*Nj)/((Ni + Sij)*(Ni + Nj + Sij))*1000000)
      
      bicMin = Inf
      kMin = 0
      for(k in 1:5){
        est = estMTN(x = x, k = k, l = min(x), u = max(x), nIter = 100)
        if(!is.na(est$bic) & is.finite(est$bic)){
          if(est$bic < bicMin){
            bicMin = est$bic
            kMin = k
          }
        } else{
          break
        }
      }
      est = estMTN(x = x, k = kMin, l = min(x), u = max(x), nIter = 100)
      
      buffer = rep(0, 5 - kMin)
      RES[counter, ] = c(indTemp[i], indTemp[j], est$m, buffer, est$s, buffer + 1, est$w, buffer, est$bic, min(x), max(x))
      counter = counter + 1
      counterOuter = counterOuter + 1
      if(counterOuter/100 == floor(counterOuter/100)){
        print(paste("Counter:", counterOuter))
        estTime = t1 + (Sys.time() - t1)/counterOuter*N
        print(paste("Estimated time:", estTime))
      }
    }
  }
  RES = data.frame(RES)
  colnames(RES) = c("Origin", "Destination", paste0("mu", 1:5), paste0("sigma", 1:5), paste0("w", 1:5), "bic", "l", "u")
  
  write.csv(x = RES, file = paste("ParameterPrep/Results/Radiation", c("N", "E", "S", "W")[r], ".csv", sep = ""))
}
