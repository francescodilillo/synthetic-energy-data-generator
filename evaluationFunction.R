#       EVALUATION FUNCTION           #
# ----------------------------------- #

# LIB
# ----------------------------------- #
library(plyr)
library(scales)
library(flux)
library(infotheo)
library(evd)
library(tseries)
library(forecast)
library(fitdistrplus)
library(stats)
library(pracma)
library(TSdist)


# FUNC
# ----------------------------------- #

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

evaluateStep1 <- function(synth, index, pathTR){
  
  path1 <- paste(pathTR, "/",index,"/stats.csv", sep = "")
  path3 <- paste(pathTR, "/",index,"/acf.csv", sep = "")
  
  stats <- read.csv(path1)
  stats <- stats[2:length(stats)]
  
  ACFscore <- 0
  
  if(length(unique(synth)) > 1){
    rACF <- read.csv(path3)
    rACF <- data.matrix(rACF[, c(2:length(rACF))])
    mLag <- rACF[nrow(rACF),2]
    synthACF <- acf(synth, lag.max = mLag, plot = FALSE)
    for(i in 1:mLag){
      
      ACFscore <- ACFscore + (as.numeric(rACF[i,1]) - as.numeric(synthACF$acf[i]))^2
      
    }
    }
  
  minScore <- stats$minimum - min(synth, na.rm = TRUE)
  maxScore <- stats$maximum - max(synth, na.rm = TRUE)
  meanScore <- stats$mean - mean(synth, na.rm = TRUE)
  sdScore <- stats$sd - sd(synth, na.rm = TRUE)
  ampScore <- stats$amplitude - amplitudeCalculator(synth, 0)

  result <- data.frame(minScore, maxScore, meanScore, sdScore, ampScore, ACFscore)
  result <- data.matrix(result)
  
  return(result)
}

evaluateStep2 <- function(syntheticDataset, pathTR){
  
  path <- paste(pathTR, "/correlation.csv", sep = "")
  
  correlationMatrix <- read.csv(path)
  correlationMatrix <- data.matrix(correlationMatrix[,c(2:ncol(correlationMatrix))])

  sCorr <- cor(t(syntheticDataset))

  result <- correlationMatrix - sCorr
  
  return(result)
 
}

evaluateStep3 <- function(synth, dsynth, dim, ind, pathTR){
  
  path <- paste(pathTR, "/",ind,"/ccordist.csv", sep = "")
  
  ccordist <- read.csv(path)
  ccordist <- ccordist[2:length(ccordist)]
  difference <- 0
  result <- data.frame(ccordist)
  names(result) <- "CCF"
  
  for(i in 1:dim){
    
    dist1 <- CCorDistance(synth, dsynth[i,])
    difference[i] <- abs(ccordist$x[i]-dist1)
    
  }
  
  difference <- data.frame(difference)
  names(difference) <- "DIFF_CCF"
  
  result <- cbind(result$CCF[1:dim], difference$DIFF_CCF[1:dim])
  
  return(result)
  
}

