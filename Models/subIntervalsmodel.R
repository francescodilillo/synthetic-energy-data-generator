#_________________________________________________________________________________________________#
#                                                                                                 #
#                                    SUB INTERVALS MODEL                                          #
#                                      MODEL GENERATOR                                            #
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
      
      if(!is.na(array[i])){
        
        if(extrA <= array[i] & array[i] <= extrB & x<=l){
          
          index <- bMatrix[x, 2]
          
          p[b, index] <- p[b, index] + 1
          
        }}
      
    }
    
  }
  
  # if(debugG == 1){
  #   
  #   print(p)
  #   write.csv(p, paste("probCalculatorError",runif(1),".csv"))
  #   
  # }
  
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
      if(!is.na(array[i])){
        if(array[i] >= a & array[i] <= b){
          
          p[i,] <- c(array[i], k)
          dd <- TRUE
          
        }}
      
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

divisors <- function(x){
  #  Vector of numberes to test against
  y <- seq_len(x)
  #  Modulo division. If remainder is 0 that number is a divisor of x so return it
  y[ x%%y == 0 ]
}

intervalIdentifier <- function(trainingData, limitV){
  
  l <- 1
  intervals <- 0
  nchange <- 1
  intervals[nchange] <- 1
  cond <- FALSE
  
  while(l < length(trainingData)){
    
    if(trainingData[l] <= limitV && cond == FALSE){
      
      cond <- TRUE
      nchange <- nchange + 1
      intervals[nchange] <- l
      l <- l+1
      
    } else if(trainingData[l] >= limitV && cond == TRUE){
      
      cond <- FALSE
      nchange <- nchange + 1
      intervals[nchange] <- l
      l <- l+1
      
    } else{
      
      l <- l+1
      
    }
  }
  
  return(intervals)
  
}

intervalSummarizer <- function(intervals, window){
  
  temp <- intervals
  k <- 1
  
  for(i in 2:length(intervals)){
    
    past <- temp[k]
    
    if(intervals[i]-past > window){
      
      k <- k+1
      temp[k] <- intervals[i]
      
    }
    
  }
  
  intervals <- temp[1:k]
  intervals[length(intervals)+1] <- length(subject)
  
  return(intervals)
}

summaryGenerator <- function(subject, intervals, limitV){
  
  if(debugG == 1){
    
    print(limitV)
    print(intervals)
    
  }
  
  synthesis <- matrix(0, nrow = 3, ncol = length(intervals))
  
  for(i in 1 : length(intervals)){
    
    synthesis[1,i] <- intervals[i]
    
    if(i <= (length(intervals)-1)){
      
      index <- (intervals[i+1] - intervals[i])%/%2
      
      if(subject[(intervals[i]+index)] >= limitV){
        
        synthesis[2,i] <- FALSE
        
      } else {
        
        synthesis[2,i] <- TRUE
        
      }
      
    } else {
      
      synthesis[2,i] <- FALSE
      
    }
    
  }
  
  for(i in 1:(length(synthesis[1,]) -1)){
    
    x <- sort(table(subject[intervals[i]:intervals[i+1]]),decreasing=TRUE)[1:3]
    x <- data.frame(x)
    
    if(length(x$Var1)==0){
      
      xx <- as.numeric(x$x[1])
      
    } else {
      
      xx <- as.numeric(as.matrix(x$Var1)[1])
      
    }
    
    synthesis[3,i] <- xx
    
  }
  
  return(synthesis)
  
}

subIntervalIdentifier <- function(s2, index, m2, cond3){
  
  subIntervals <- 0
  subIntervals[1] <- 1
  nchange <- 2
  counter <- index[1]

  stopC <- FALSE
  
  if(!is.na(counter)){
  
    subIntervals[nchange] <- counter
    nchange <- 3
      
    while(i < length(s2) && !stopC){
      
      if(cond3 == TRUE){
        
        index2 <- which(s2[counter:length(s2)] < m2)
        if(length(index2) == 0){
          
          stopC <- TRUE
          
        } else {
          
          subIntervals[nchange] <- index2[1] + counter
          counter <- counter+index2[1]
          nchange <- nchange+1
          cond3 <- FALSE
          
        }
      }else {
        #print(c(counter, length(s2), m2))
        index2 <- which(s2[counter:length(s2)] > m2)
        if(length(index2) == 0){
          
          stopC <- TRUE
          
        } else {
          
          subIntervals[nchange] <- index2[1] + counter
          counter <- counter+index2[1]
          nchange <- nchange+1
          cond3 <- TRUE
          
        }
      }
      
    }
    
    subIntervals <- unique(subIntervals)
  
  } else {
    
    subIntervals[nchange] <- length(s2)
    
  }
  
  #plot(s2, type = "l")
  #abline(v = subIntervals, col = "blue")
  
  print(subIntervals)
  
  return(subIntervals)
  
}

subIntervalSummarizer <- function(subIntervals){
  
  save <- rep(FALSE, length(subIntervals))
  save[1] <- save[length(subIntervals)] <- TRUE
  
  for(j in 1:(length(subIntervals)-1)){
    
    dist2 <- subIntervals[j+1] - subIntervals[j]
    
    if(dist2 > 20){
      
      save[j+1] <- TRUE
      
    }
    
  }
  
  return(save)
  
}

reshape <- function(subIntervalsMatrix){
  
  rIndex <- 0
  cIndex <- 0
  k2k <- 1
  p2p <- 1
  
  for(i in 1:nrow(subIntervalsMatrix)){
    
    if(sum(subIntervalsMatrix[i,])==0){
      
      rIndex[k2k] <- i
      k2k <- k2k + 1
      
    }
    
  }
  
  for(j in 1:ncol(subIntervalsMatrix)){
    
    if(sum(subIntervalsMatrix[,j])==0){
      
      cIndex[p2p] <- j
      p2p <- p2p + 1
      
    }
    
  }
  
  matrixxx <- subIntervalsMatrix[-rIndex, -cIndex]
  
  return(matrixxx)
  
}

stateListGenerator <- function(intervals, subject, indexes, path, rows){
  
  nval <- 0
  stateMatrix <- list()
  histograms <- list()
  
  for(i in 1:rows){
    
    if(rows != 1){
      
      
      oSubject <- intervals[i,][!is.na(intervals[i,])]
      oSubject <- oSubject[oSubject>0]
      oSubject <- oSubject[!is.na(oSubject)]
      cond <- TRUE
      
    } else {
      
      
      oSubject <- intervals[intervals > 0]
      oSubject <- oSubject[!is.na(oSubject)]
      cond <- FALSE
      
    }
    
    dimO <- length(oSubject) - 1
    
    for(j in 1:dimO){
      if(i > 1){
        trendUS <- subject[(oSubject[j]+indexes[i]):(oSubject[j+1]+indexes[i])]
      }
      else{
        #print(c(oSubject[j],oSubject[j+1]))
        trendUS <- subject[oSubject[j]:oSubject[j+1]]
      }
      #plot(trendUS, type = "l")
      temp <- hist(trendUS, breaks = 51, plot = FALSE)
      temp2 <- data.frame(temp$breaks, c(temp$density,0))
      
      if(!cond){
        
        write.csv(temp2, paste(path, "/histogram",j,".csv", sep =""))
        
      } else {
        
        write.csv(temp2, paste(path, "/State",i,"histogram",j,".csv", sep =""))
        
      }
      
      histograms[[j]]<- temp
      
      nval[j] <- length(histograms[[j]]$breaks)-1
      state <- seq(from = 1, to = nval[j], by = 1)
      l <- length(state)
      
      if(j==1){
        
        maxlength <- l
        stateList <- matrix(0, nrow = dimO, ncol = maxlength)
        stateList[1,] <- state
        
      } else{
        
        if(l<maxlength){
          
          diff <- maxlength-l
          z <- rep(0,diff)
          state <- c(state,z)
          stateList[j,] <- state
          
        } else{
          
          maxlength <- l
          tempMatrix <- matrix(0, nrow = dimO, ncol = maxlength)
          index <- j-1
          
          for(jj in 1:index){
            
            pState <- stateList[jj,]
            pL <- length(pState)
            diff <- maxlength - pL
            z <- rep(0, diff)
            pState <- c(pState, z)
            tempMatrix[jj,] <- pState
            
          }
          
          tempMatrix[j,] <- state
          stateList <- tempMatrix
          
        }
        
      }
      
    }
    
    if(rows>1){
      
      stateMatrix[[i]] <- stateList
      stateListSaver(path, 2, stateMatrix)
      
    } else {
      
      stateListSaver(path, 1 , stateList)
      
    }
    
  }
  
  if(rows > 1){
    
    return(stateMatrix)
    
  } else {
    
    return(stateList)
    
  }
  
}

stateListSaver <- function(path, op, stateList){
  
  if(op == 1){
    
    write.csv(stateList, paste(path,"/stateList.csv", sep = ""))
    
  }
  
  else{
    
    iter <- length(stateList)
    
    for(i in 1:iter){
      
      write.csv(stateList[[i]], paste(path,"/subStateList",i,".csv", sep = ""))
      
    }
    
  }
  
  
}

transitionMatrixBuilder <- function(intervals, subject, path, rows){
  
  r <- rows
  transition <- list()
  cond <- TRUE
  
  for(i in 1:r){
    
    if(r > 1){
      
      oSubject <- intervals[i,][!is.na(intervals[i,])]
      oSubject <- oSubject[oSubject>0]
      oSubject <- oSubject[!is.na(oSubject)]
      cond <- TRUE
      
    } else {
      
      oSubject <- intervals[intervals > 0]
      oSubject <- oSubject[!is.na(oSubject)]
      cond <- FALSE
    }
    
    
    dimO <- length(oSubject)-1
    
    for(j in 1:dimO){
      
      
      if(cond == FALSE){
        
        histograms <- read.csv(paste(path,"/histogram",j,".csv", sep = ""))
        
      } else{
        
        histograms <- read.csv(paste(path,"/State",i,"histogram",j,".csv", sep = ""))
        
      }
      histograms <- histograms[, c(2,3)]
      names(histograms) <- c("breaks", "density")
      
      if(length(histograms$breaks)!= 2){
        
        temp2 <- probCalculator(subject[oSubject[j]:oSubject[j+1]], histograms$breaks)
        if(any(is.na(temp2))>0){
          
          # if(debugG == 1){
          # 
          #   print(c("NULL in probCalculator, iter",j))
          #   print(c(oSubject[j], oSubject[j+1]))
          #   write.csv(temp2, paste(path, "/errormatrix.csv", sep =""))
          #   
          # }
          
          temp2 <- probCalib(temp2)
          
        }
      } else { 
        
        temp2 <- matrix(1, nrow = 1, ncol = 1) 
        
      }
      
      transition[[j]] <- temp2
      
      if(cond == TRUE){
        
        write.csv(transition[[j]], paste(path, "/State",i,"subState",j,"transitions.csv", sep = "")) 
        
      }
      else {
        
        write.csv(transition[[j]], paste(path, "/transitions",j,".csv", sep = "")) 
        
      }
    }
    
  }
  
  return(transition)
  
}

# Function
#----------------------------#

subIntervalsGenerator <- function(path3, trainingData, intervals){

  lint <- length(intervals)

    
  for(il in 1:2){
  
  # SUB-INTERVALS IDENTIFIER
  path1 <- paste(path3, "/", il, sep = "")
  dir.create(path1, showWarnings = FALSE)
  subIntervalsMatrix <- matrix(0, nrow = lint, ncol = 50)
  
  for(jj in 1:(lint-1)){
    
    a <- intervals[jj]
    b <- intervals[jj+1]
    cond3 <- FALSE
    s2 <- trainingData[a:b]
    incr <- 0
    fIndex <- 0
    rIndex <- 0
    
    limitA <- 0
    
    limitA[1] <- mean(s2, na.rm= TRUE)
    limitA[2] <- limitA[1] - abs(min((s2 - mean(s2, na.rm= TRUE))))/2
    
    if(s2[1] < limitA[il]){
      
      cond3 <- TRUE
      index <- which(s2 > limitA[il])
      
    } else {
      
      cond3 <- FALSE
      index <- which(s2 < limitA[il])
      
    }
    
    v <- sort(table(s2[index]),decreasing=TRUE)[1]
    v <- as.numeric(names(v))
    # v <- data.frame(v)
    # 
    # if(class(v$Var1) != "NULL"){
    #   
    #   vv <- as.numeric(as.matrix(v$Var1)[1])
    #   
    # } else {
    #   
    #   vv <- v$v[1]
    #   
    # }
    # 
    #kk <- 1
    limitV <- v #vv
    
    if(debugG != 1){
      
      plot(s2, type = "l")
      abline(h = limitA[il], col = "blue")
      
    }
    
    subIntervals <- subIntervalIdentifier(s2, index, limitA[il], cond3)
    
    if(subIntervals[length(subIntervals)] != length(s2)){
    
    subIntervals[length(subIntervals)+1] <- length(s2)}
    
    # SUB-INTERVALS SUMMARIZER
    
    saveIndex <- subIntervalSummarizer(subIntervals)
    
    diffIndex <- ncol(subIntervalsMatrix) - length(subIntervals[saveIndex])
    
    if(diffIndex >0){
    
    subIntervalsMatrix[jj,] <- c(subIntervals[saveIndex],rep(0,diffIndex))
    
    } else {
      
      if(length(subIntervals[saveIndex])==ncol(subIntervalsMatrix)){
      subIntervalsMatrix[jj,] <- subIntervals[saveIndex]
      } else {
        
        subIntervalsMatrix[jj,] <- subIntervals[saveIndex][1:ncol(subIntervalsMatrix)]
        
      }
      
    }
  }
  
  subIntervalsMatrix <- reshape(subIntervalsMatrix)
  
  # modSubIntervalsMatrix <- subIntervalsMatrix
  # 
  # for(i in 1:nrow(modSubIntervalsMatrix)){ 
  #   
  #   l <- length(modSubIntervalsMatrix[i,][modSubIntervalsMatrix[i,] > 0])
  #   
  #   modSubIntervalsMatrix[i,1] <- modSubIntervalsMatrix[i,1] + 20
  #   
  #   modSubIntervalsMatrix[i,l] <- modSubIntervalsMatrix[i,l] - 20
  #   
  #   
  # }
  
  np1 <- paste(path1, "/subIntervalsMatrix.csv", sep = "")
  write.csv(subIntervalsMatrix, np1)
  
  # np2 <- paste(path1, "/modSubIntervalsMatrix.csv", sep = "")
  # write.csv(modSubIntervalsMatrix, np2)
  
  rows <- nrow(subIntervalsMatrix)
  
  stateList <- stateListGenerator(subIntervalsMatrix, trainingData, intervals, path1, rows)
  
  transition <- transitionMatrixBuilder(subIntervalsMatrix, trainingData, path1, rows)
  
  }

}