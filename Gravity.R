rm(list = ls())

library("sf")
library("geosphere")

# Set wd to project folder here

source("wspsample.R")
source("PopBoundaryAPI.R")
source("ParameterPrep/Fitting.R")


ID2Geom = read.csv("Data/IDlookup.csv", header = T, fileEncoding = "UTF-8")
IDs = ID2Geom$ID

countryCode = "DEU"
level = 3
year = 2011

# Download data
# Setting wd to population data folder

Ger = getJSON(countryCode, level, timeout = 1000)
GerS = as_Spatial(Ger$geometry[ID2Geom$JsonNr])


# Download population data
pop = getPOP(countryCode, year, 1000)
# Setting wd back to project folder

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

t1 = Sys.time()
for(r in 1:4){
  indTemp = ind[[r]]
  n = length(indTemp)
  RES = matrix(nrow = n^2/2 + n/2, ncol = 6*5 + 3)
  counter = 1
  for(i in 1:n){
    for(j in i:n){
      D = distGeo(samO[[indTemp[i]]][, 1:2], samD[[indTemp[j]]][, 1:2])
      x = cbind(log(D/1000), log(samO[[indTemp[i]]][, 3]/samD[[indTemp[i]]][, 3]))
      bicMin = Inf
      kMin = 0
      for(k in 1:5){
        est = estMN2(x = x, k = k, nIter = 100)
        if(!is.na(est$bic) & is.finite(est$bic)){
          if(est$bic < bicMin){
            bicMin = est$bic
            kMin = k
          }
        } else{
          break
        }
      }
      est = estMN2(x = x, k = kMin, nIter = 100)
      
      buffer = rep(0, 5 - kMin)
      RES[counter, ] = c(i, j, est$m[, 1], buffer, est$m[, 2], buffer, est$S[, 1],
                         buffer + 1, est$S[, 2], buffer, est$S[, 4], buffer + 1,
                         est$w, buffer, est$bic)
      counter = counter + 1
    }
  }
  RES = data.frame(RES)
  colnames(RES) = c("Origin", "Destination", paste0("mu1", 1:5),
                    paste0("mu2", 1:5), paste0("sigma1", 1:5), paste0("sigma12", 1:5),
                    paste0("sigma2", 1:5), paste0("w", 1:5), "bic")
  
  write.csv(x = RES, file = paste("ParameterPrep/Results/Gravity", c("N", "E", "S", "W")[r], ".csv", sep = ""))
}

