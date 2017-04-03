#_________________________________________________________________________________________________#
#                                                                                                 #
#                                     MACRO INTERVAL TREND                                        #
#                                       MODEL GENERATOR                                           #
#                                                                                                 #
#                                     Francesco Di Lillo                                          #
#                             Karlsruher Institut für Technologie                                 #
#                       Scuola Politecnica e delle Scienze di Base Federico II                    #
#_________________________________________________________________________________________________#


# Libraries
#----------------------------#
library(plyr)
library(scales)
library(flux)
library(infotheo)
library(evd)
library(tseries)
library(forecast)
library(fitdistrplus)
library(stats)

#----------------------------#

# Utility
#----------------------------#

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
      
      if(extrA <= array[i] & array[i] <= extrB & x<=l){
        
        index <- bMatrix[x, 2]
        
        p[b, index] <- p[b, index] + 1
        
      }
      
    }
    
  }
  
  for(s in 1:z){
    
    div <- sum(p[s,])
    
    for(v in 1:z){
      
      p[s,v] <- p[s,v]/div
      
    }
    
  }
  
  return(p)
  
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
      
      if(array[i] >= a & array[i] <= b){
        
        p[i,] <- c(array[i], k)
        dd <- TRUE
        
      }
      
      k <- k+1
      
    }
    
  }
  
  return(p)
  
}

probCalib <- function(matrixX){
  
  dim <- nrow(matrixX)
  dim2 <- ncol(matrixX)
  
  for(i in 1:dim){
    
    rrow <- matrixX[i,]
    
    if(all(is.na(rrow))){
      
      rrow <- rep((1/dim2), dim2)
      
    }
    
    matrixX[i,] <- rrow
    
  }
  
  return(matrixX)
  
}

# Function
#----------------------------#

macroIntervalGenerator <- function(path1, trainingData){
  
  path2 <- paste(path1,"/macroInterval", sep = "")
  dir.create(path2, showWarnings = FALSE)
  
  histograms <- 0
  nval <- 0
  
  # --> HISTOGRAM GENERATOR
    
    temp <- hist(trainingData, breaks = 51, plot = FALSE)
    temp2 <- data.frame(temp$breaks, c(temp$density,0))
    write.csv(temp2, paste(path2,"/histogram.csv", sep =""))
    histograms<- temp
    
    nval <- length(histograms$breaks)-1
    state <- seq(from = 1, to = nval, by = 1)
    l <- length(state)
    
  # --> STATE LIST GENERATOR      
    maxlength <- l
    stateList <- 0
    stateList <- state
  
    write.csv(stateList, paste(path2,"/stateList.csv", sep = ""))
    
  # --> TRANSITION MATRIX
    
    transitions <- 0
      
    temp2 <- probCalculator(trainingData, histograms$breaks)
    transitions <- temp2
    if(any(is.na(temp2))>0){
  
      temp2 <- probCalib(temp2)
      
    }
    write.csv(transitions, paste(path2,"/transitions.csv", sep = ""))
  
}