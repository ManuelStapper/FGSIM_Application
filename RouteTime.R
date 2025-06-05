rm(list = ls())

library("sf")
library("geosphere")
library("osrm")

setwd("/Users/lshms101/Desktop/Projects/FGSIM")

source("wspsample.R")
source("PopBoundaryAPI.R")
source("ParameterPrep/Fitting.R")


ID2Geom = read.csv("Data/IDlookup.csv", header = T, fileEncoding = "UTF-8")
IDs = ID2Geom$ID

countryCode = "DEU"
level = 3
year = 2011

# Download data
setwd("/Users/lshms101/Desktop/Projects/FGSIM/ParameterPrep/PopDistrict")

Ger = getJSON(countryCode, level, timeout = 1000)
GerS = as_Spatial(Ger$geometry[ID2Geom$JsonNr])


# Download population data
pop = getPOP(countryCode, year, 1000)
setwd("/Users/lshms101/Desktop/Projects/FGSIM")

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

getTime = function(origin, destination, batchSize = 100){
  n = nrow(origin)
  batch = floor(n/batchSize)
  nDist = batchSize^2/2 - batchSize/2
  out = rep(0, batch*nDist)
  for(b in 1:batch){
    ind = ((b-1)*batchSize + 1):(b*batchSize)
    temp = osrmTable(origin[ind, ], destination[ind, ], osrm.profile = "car", measure = "duration")$durations
    out[((b-1)*nDist + 1):(b*nDist)] = temp[upper.tri(temp)]
  }
  out
}

t1 = Sys.time()
for(r in 1:4){
  indTemp = ind[[r]]
  n = length(indTemp)
  RES = matrix(nrow = n^2/2 + n/2, ncol = 2 + 3*5 + 3)
  counter = 1
  for(i in 1:n){
    for(j in i:n){
      x = getTime(samO[[i]][1:100, 1:2], samD[[j]][1:100, 1:2], 100)
      x = log(x + runif(length(x))/10)
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
    }
  }
  RES = data.frame(RES)
  colnames(RES) = c("Origin", "Destination", paste0("mu", 1:5), paste0("sigma", 1:5), paste0("w", 1:5), "bic", "l", "u")
  
  write.csv(x = RES, file = paste("ParameterPrep/Results/RouteTime", c("N", "E", "S", "W")[r], ".csv", sep = ""))
}


