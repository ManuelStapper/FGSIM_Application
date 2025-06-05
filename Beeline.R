library("sf")
library("geosphere")

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

D = array(NA, dim = c(401, 401, 2000))
for(i in 1:401){
  iRegion = which(unlist(lapply(ind, function(xx) i %in% xx)))
  for(j in ind[[iRegion]]){
    if(i <= j){
      D[i, j, ] = D[j, i, ] = distGeo(samO[[i]][, 1:2], samD[[j]][, 1:2])
    }
  }
}

n = unlist(lapply(ind, length))
N = sum(n^2/2 + n/2)
counterOuter = 1
t1 = Sys.time()
for(r in 1:4){
  indTemp = ind[[r]]
  n = length(indTemp)
  RES = matrix(nrow = n^2/2 + n/2, ncol = 2 + 3*5 + 3)
  counter = 1
  for(i in 1:n){
    for(j in i:n){
      x = log(D[indTemp[i], indTemp[j], ]/1000)
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
        print(paste("Estimated Time:", estTime))
      }
    }
  }
  RES = data.frame(RES)
  colnames(RES) = c("Origin", "Destination", paste0("mu", 1:5), paste0("sigma", 1:5), paste0("w", 1:5), "bic", "l", "u")
  
  write.csv(x = RES, file = paste("ParameterPrep/Results/Beeline", c("N", "E", "S", "W")[r], ".csv", sep = ""))
}
