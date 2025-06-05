# Include weights wrapper here
# Double checking the transformations

setwd("/Users/lshms101/Desktop/Projects/FGSIM/ParameterPrep/Results")

classes = c("Beeline", "RouteTime", "Gravity", "Circle", "Radiation")
classesAbbr = c("B", "T", "G", "C", "R")


# For Beeline, RouteTime, Circle and Radiation:
# list of three arrays (mu, sigma, weights)
# One matrix for lower ends, one matrix for upper

# For Gravity:
# list of six arrays (mu1, mu2, sigma1, sigma2, sigma12, weights)

### Caution: Only for gravity model, sigmas give variances, not stds!

for(cl in c(1, 2, 4, 5)){
  for(reg in 1:4){
    nmWrite = paste("pars", classesAbbr[cl], "_", c("N", "E", "S", "W")[reg], sep = "")
    nmRead = paste(classes[cl], c("N", "E", "S", "W")[reg], ".csv", sep = "")
    RES = read.csv(nmRead, header = T)[, -1]
    ind = sort(unique(RES[, 1]))
    n = length(ind)
    pars = list(mu = array(0, dim = c(n, n, 5)),
                sigma = array(0, dim = c(n, n, 5)),
                w = array(0, dim = c(n, n, 5)),
                l = matrix(0, nrow = n, ncol = n),
                u = matrix(0, nrow = n, ncol = n))
    
    for(i in 1:nrow(RES)){
      o = which(ind == RES[i, 1])
      d = which(ind == RES[i, 2])
      
      for(k in 1:5){
        pars$mu[d, o, k] = RES[i, 2 + k]
        pars$sigma[d, o, k] = RES[i, 7 + k]
        pars$w[d, o, k] = RES[i, 12 + k]
        
        if(cl %in% c(1, 2)){
          pars$mu[o, d, k] = pars$mu[d, o, k]
          pars$sigma[o, d, k] = pars$sigma[d, o, k]
          pars$w[o, d, k] = pars$w[d, o, k]
        }
      }
      pars$l[d, o] = RES[i, 19]
      pars$u[d, o] = RES[i, 20]
      if(cl %in% c(1, 2)){
        pars$l[o, d] = pars$l[d, o]
        pars$u[o, d] = pars$u[d, o]
      }
    }
    assign(nmWrite, pars)
  }
}

for(reg in 1:4){
  nmWrite = paste("parsG", "_", c("N", "E", "S", "W")[reg], sep = "")
  nmRead = paste("Gravity", c("N", "E", "S", "W")[reg], ".csv", sep = "")
  RES = read.csv(nmRead, header = T)[, -1]
  ind = sort(unique(RES[, 1]))
  n = length(ind)
  pars = list(mu1 = array(0, dim = c(n, n, 5)),
              mu2 = array(0, dim = c(n, n, 5)),
              sigma1 = array(0, dim = c(n, n, 5)),
              sigma12 = array(0, dim = c(n, n, 5)),
              sigma2 = array(0, dim = c(n, n, 5)),
              w = array(0, dim = c(n, n, 5)))
  
  for(i in 1:nrow(RES)){
    o = which(ind == RES[i, 1])
    d = which(ind == RES[i, 2])
    
    for(k in 1:5){
      pars$mu1[d, o, k] = RES[i, 2 + k]
      pars$mu2[d, o, k] = RES[i, 7 + k]
      pars$sigma1[d, o, k] = RES[i, 12 + k]
      pars$sigma12[d, o, k] = RES[i, 17 + k]
      pars$sigma2[d, o, k] = RES[i, 22 + k]
      pars$w[d, o, k] = RES[i, 27 + k]
      
      pars$mu1[o, d, k] = pars$mu1[d, o, k]
      pars$mu2[o, d, k] = -pars$mu2[d, o, k] # Sign inverted!
      pars$sigma1[o, d, k] = pars$sigma1[d, o, k]
      pars$sigma12[o, d, k] = pars$sigma12[d, o, k]
      pars$sigma2[o, d, k] = pars$sigma2[d, o, k]
      pars$w[o, d, k] = pars$w[d, o, k]
    }
    
  }
  assign(nmWrite, pars)
}

# Translate expectations for simplified estimation:
# E(D)^{-d} instead of E(D^{-d})

# For univariate measures
translatePars1 = function(x){
  nn = dim(x$mu)[1]
  out = matrix(0, nrow = nn, ncol = nn)
  for(i in 1:nn){
    for(j in 1:nn){
      m = x$mu[i, j, ]
      s = x$sigma[i, j, ]
      l = x$l[i, j]
      u = x$u[i, j]
      w = x$w[i, j, ]
      scl = (pnorm((u - m)/s - s) - pnorm((l - m)/s - s)) / (pnorm((u - m)/s) - pnorm((l - m)/s))
      scl[is.na(scl)] = 1.0
      out[i, j] = sum(w * exp(m + s^2/2) * scl)
    }
  }
  list(mu = out)
}

# For gravity model
translatePars2 = function(x){
  nn = dim(x$mu1)[1]
  out1 = matrix(0, nrow = nn, ncol = nn) # For distances
  out2 = matrix(0, nrow = nn, ncol = nn) # For population ratio
  for(i in 1:nn){
    for(j in 1:nn){
      m1 = x$mu1[i, j, ]
      m2 = x$mu2[i, j, ]
      s1 = x$sigma1[i, j, ]
      s2 = x$sigma2[i, j, ]
      w = x$w[i, j, ]
      out1[i, j] = sum(w * exp(m1 + s1/2))
      out2[i, j] = sum(w * exp(m2 + s2/2))
    }
  }
  list(mu1 = out1, mu2 = out2)
}

for(r in 1:4){
  for(cl in c(1, 2, 4, 5)){
    nmIn = paste0("pars", classesAbbr[cl], "_", c("N", "E", "S", "W")[r])
    nmOut = paste0("pars", classesAbbr[cl], "s_", c("N", "E", "S", "W")[r])
    parsTemp = get(nmIn)
    assign(nmOut, translatePars1(parsTemp))
  }
  
  nmIn = paste0("pars", "G_", c("N", "E", "S", "W")[r])
  nmOut = paste0("pars", "Gs_", c("N", "E", "S", "W")[r])
  parsTemp = get(nmIn)
  assign(nmOut, translatePars2(parsTemp))
}




keep = c(outer(paste("pars", classesAbbr, "_", sep = ""), c("N", "E", "S", "W"), FUN = paste0))
keep = c(keep, c(outer(paste("pars", classesAbbr, "s_", sep = ""), c("N", "E", "S", "W"), FUN = paste0)))

rm(list = setdiff(ls(), keep))

setwd("/Users/lshms101/Desktop/Projects/FGSIM")
save.image("Parameters.RData")







