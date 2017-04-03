library(plyr)
library(scales)
library(flux)
library(infotheo)
library(evd)
library(tseries)
library(forecast)
library(fitdistrplus)
library(stats)
library(TSdist)
library(quantmod)



decomposition <- function(matrix){
  
  result <- matrix(0, nrow = 2*nrow(matrix), ncol = ncol(matrix))
  j <- 1
  
  
  for(i in 1:nrow(matrix)){
    
    temp <- decompose(ts(matrix[i,], frequency = 24, start = c(1,1)))
    
    noise <- temp$random
    seasonal <- temp$seasonal
    
    result[j, ] <- noise
    result[j+1, ] <- seasonal
    
    j <- j+2
    
  }
  
  return(result)
}

getSeasonality <- function(matrix){
  
  ris <- matrix(0, nrow = nrow(matrix), ncol = 96)
  
  for(i in 1:12){
    
    var <- matrix[i,]
    
    temp <- fft(var[1:96])
    
    ris[i, ] <- temp
    
  }
  
  return(ris)
  
}

probCalculator <- function(array, breaks){
  
  z <- length(breaks)-1
  l <- length(array)
  p <- matrix(0, ncol = z, nrow = z)
  
  bMatrix <- breakIdentifier(array, breaks)
  
  for( b in 1:z){
    
    extrA <- breaks[b]
    extrB <- breaks[b+1]
    x <- 1
    
    for (i in 1:l){
      
      x <- i+1
      
      if(extrA <= array[i] && array[i] <= extrB && x<=l){
        
        index <- bMatrix[x, 2]
        
        p[b, index] <- p[b, index] + 1
        
      }
      
    }
    
  }
  
  for(s in 1:z){
    
    div <- sum(p[s,])
    
    if(div == 0 ){
      
      p[s, ] <- rep(1/length(p[s,]), length(p[s,]))
      
    } 
    
    else {
      
      for(v in 1:z){
        
        p[s,v] <- p[s,v]/div
        
      }
      
    }
  }
  
  return(p)
  
}

tsTrendExtraction <- function(m){
  
  r <- nrow(m)
  c <- ncol(m)
  j <- 1
  i <- 2
  
  tsm <- matrix(0, nrow = (c/2), ncol = r)
  
  while(i <= c){
    
    var <- decompose(ts(m[, i], frequency = 24, start = c(12,2016)))
    
    tsm[j, ] <- var$trend
    
    j <- j+1
    i <- i+2
  }
  
  return(tsm)
  
}

breakIdentifier <- function(array, breaks){
  
  l <- length(array)
  z <- length(breaks)-1
  p <- matrix(0, ncol = 2, nrow = l)
  
  for(i in 1:l){
    
    k <- 1
    dd <- FALSE
    
    while(k<z && dd == FALSE){
      
      a <- breaks[k]
      b <- breaks[k+1]
      
      if(array[i] >= a && array[i] <= b){
        
        p[i,] <- c(array[i], k)
        dd <- TRUE
        
      }
      
      k <- k+1
      
    }
    
  }
  
  return(p)
  
}

amplitudeCalculator <- function(trendTs, debug){
  
  results <- 0
  
  real <- trendTs
  real <- real[!is.na(real)]
  
  maxAmpl <- 0
  ind <- 0
  
  a <- 1
  b <- length(real)
  
  for(i in a:b){
    
    j <- seq(i+1, i+100, 1)
    
    max <- 0
    k <- 1
    flag <- FALSE
    mode <- 0
    
    while(k <= (length(j)-1) && !flag && i+k <= b){
      
      c <- abs(real[j[k]] - real[i])
      
      flag1 <- c >= mode
      
      if(c > max){
        
        max <- c
        ind <- seq(from = i, to = j[k], by = 1)
        
      }
      
      k <- k+1
      mode <- c
    }
    
    if(max > maxAmpl){
      
      maxAmpl <- max
      mInd <- ind
      
    }
    
  }
  results <- maxAmpl
  
  if(debug == 1){
    
    plot(real[1:length(real)], type = "l")
    segments(x0 = mInd[1], y0 = real[mInd[1]], x1 = mInd[length(mInd)], y1 = real[mInd[length(mInd)]], col = "red")
    print(c(p, results[p]))
    
  }
  
  return(results)
}

#----------------------------#

debug <- 0

datapath1 <- "Dataset 28-02/AVT_"
topDir <- "v5_ModelBuilder2017-03-23"
sensors <- c(12, 32, 32, 12, 32, 12, 32, 32, 32, 32)

for(pt in 1:10){
  
  if(pt != 10){
  datapath2 <- paste(datapath1, "0", pt, ".csv", sep = "")
  } else {
    datapath2 <- paste(datapath1, pt, ".csv", sep = "")
  }
  
  print(paste("STATISTICHE MACCHINA AVT_", pt, sep = ""))
  
  topDir2 <- paste(topDir, "/Sim", pt, sep = "")
  
  data <- data.matrix(read.csv(datapath2))    # Historical data
  data <- data[-1,-1]
  
  time <- ncol(data)
  
  # TREND MODELING
  
  trendTs <- tsTrendExtraction(data)
  dimC <- length(trendTs[1,][!is.na(trendTs[1,])])
  dimR <- nrow(trendTs)
  temp <- matrix(0, nrow = dimR, ncol = dimC)
  for(i in 1:dimR){
    
    temp[i,] <- trendTs[i,][!is.na(trendTs[i,])]
    
  }
  
  trendTs <- matrix(0, nrow = dimR, ncol = dimC)
  trendTs <- temp
  
  cc <- cor(t(trendTs))
  
  write.csv(cc, paste(topDir2,"/correlation.csv", sep = ""))

  for(i in 1: sensors[pt]){
    
    print(paste("SENSORE ", i, sep = ""))
    
    path <- paste(topDir2, "/",i, "/", sep = "")
    
    minimum <- min(trendTs[i,], na.rm = TRUE)
    maximum <- max(trendTs[i,], na.rm = TRUE)
    mean <- mean(trendTs[i,], na.rm = TRUE)
    sd <- sd(trendTs[i,], na.rm = TRUE)
    amplitude <- amplitudeCalculator(trendTs[i,][!is.na(trendTs[i,])], debug)
    
    toPrint <- data.frame(minimum, maximum, mean, sd, amplitude)
    write.csv(toPrint, paste(path, "/stats.csv", sep = ""))
    
    peaks <- findPeaks(trendTs[i,],thresh = mean(trendTs[i,], na.rm = TRUE))
    write.csv(peaks, paste(path, "/peaks.csv", sep = ""))
    
    if(length(unique(trendTs[i,][!is.na(trendTs[i,])])) != 1){
      acf <- acf(trendTs[i,][!is.na(trendTs[i,])])
      toPrint <- data.frame(acf$acf, acf$lag)
      write.csv(toPrint, paste(path, "/acf.csv", sep = ""))}
    
    CCorDist <- 0
    
    for(j in 1:nrow(trendTs)){
      
      score <- CCorDistance(trendTs[i,], trendTs[j,])
      CCorDist[j] <- score
      
    }
    
    write.csv(CCorDist, paste(path, "/ccordist.csv", sep = ""))
    
  }
  
  
}

