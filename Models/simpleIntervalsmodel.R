#_________________________________________________________________________________________________#
#                                                                                                 #
#                                   SIMPLE INTERVALS TREND                                        #
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

# External Sources
#----------------------------#
source("subIntervalsmodel.R")

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

intervalSummarizer <- function(intervals, window, upperlimit){
  
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
  intervals[length(intervals)+1] <- upperlimit
  
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
    
    x <- sort(table(subject[intervals[i]:intervals[i+1]]),decreasing=TRUE)[1]
    #x <- data.frame(x)
    x <- as.numeric(names(x))
    
    # if(length(x$Var1)==0){
    #   
    #   xx <- as.numeric(x$x[1])
    #   
    # } else {
    #   
    # xx <- as.numeric(as.matrix(x$Var1)[1])
    # 
    # }
    # 
    # synthesis[3,i] <- xx
    
    synthesis[3,i] <- x
    
  }
  
  return(synthesis)
  
}

model1 <- function(trainingData, intervals){
  
  dim <- length(intervals)
  
  xcoord <- 0
  
  for(i in 1:(dim-1)){
    
    s <- trainingData[intervals[i]:intervals[i+1]]
      
    lim <- min(subject[intervals[i]:intervals[i+1]])
    
    tempp <- which(s == lim)  
    xcoord[i] <- tempp[1]
    xcoord[i] <- xcoord[i]+intervals[i]
    
  }
  
  newInt <- 0
  
  for(i in 1:length(xcoord)){
    
    dim2 <- length(intervals)
    
    index <- which(intervals < xcoord[i])
    index <- index[length(index)]
    
    temp1 <- intervals[1:index]
    temp2 <- intervals[index+1: dim2]
    newInt <- c(temp1, xcoord[i], temp2)
    
  }
  
  return(newInt)
  
}

model2 <- function(trainingData, intervals){
  
  dim <- length(intervals)
  iD <- 0
  hh <- 0
  
  for(i in 1:(dim-1)){
    
    hh[i] <- round(((intervals[i+1]-intervals[i])*5)/100)
    
    if( i == 1 ){
      
      a <- intervals[i]+hh
      
    } else {
      
      a <- intervals[i]
      
    }
    
    b <- intervals[i+1]
    
    s <- trainingData[a:b]
    
    temp <- which(s == max(s))
    
    iD[i] <- round(temp[1]) + a
    
  }
  
  newInt <- 0
  newInt[1] <- 1
  
  j <- 2
  
  for(i in 1:length(iD)){
    
    
    newInt[j] <- iD[i]-hh[i]
    newInt[j+1] <- iD[i]+hh[i]
    j <- j+2
    
  }
  
  newInt[j] <- length(trainingData)
  
  return(newInt)
  
}

model3 <- function(trainingData, intervals){
  
  limitSD <- sd(trainingData) - (sd(trainingData)*5)/100
  dim <- length(intervals)
  vals <- 0
  k <- 1
  
  for(i in 1:(dim-1)){
    
    localsd <- sd(trainingData[intervals[i]:intervals[i+1]])
    
    if(localsd > limitSD){
      
      min1 <- min(trainingData[intervals[i]:intervals[i+1]], na.rm = TRUE)
      ex1 <- which(trainingData[intervals[i]:intervals[i+1]] == min1)
      max1 <- max(trainingData[intervals[i]:intervals[i+1]], na.rm = TRUE)
      ex2 <- which(trainingData[intervals[i]:intervals[i+1]] == max1)
      
      if(length(ex1) > 1){
        
        ex1 <- ex1[length(ex1)]
        
      }
      
      if(length(ex2) > 1){
        
        ex2 <- ex2[length(ex2)]
        
      }
      
      t1 <- sd(subject[intervals[i]: (ex1+intervals[i])], na.rm = TRUE)
      t2 <- sd(subject[(ex1+intervals[i]): intervals[i+1]], na.rm = TRUE)
      if(!is.na(t1) && !is.na(t2)){
      s1 <- t1 - t2} else {s1 <- 100}
      
      t3 <- sd(subject[intervals[i]: (ex2+intervals[i])], na.rm = TRUE)
      t4 <- sd(subject[(ex2+intervals[i]): intervals[i+1]], na.rm = TRUE)
      if(!is.na(t3) && !is.na(t4)){
      s2 <- t3 - t4}else {s2 <- 100}
      
      
      score <- min(abs(s1), abs(s2))
      
      if(score == s1){
        
        vals[k] <- ex1+intervals[i]
        
      } else {

        vals[k] <- ex2+intervals[i]
        
      }
      
    }
  }
  
  newInt <- 0
  
  if(vals[1]!=0){
  
    for(i in 1:length(vals)){
      
      jj <- which(intervals <= vals[i])
      jj <- jj[length(jj)]
      temp <- intervals[1: jj]
      temp2 <- intervals[(jj+1):length(intervals)]
      
      newInt <- c(temp, vals[i], temp2)
      
    }
  } else {
      
    newInt <- intervals
    
    }
  
  return(newInt)
  
}

model4 <- function(trainingData, intervals, path1){
  
  m <- mean(trainingData, na.rm= TRUE)
  min1 <- min(trainingData, na.rm = TRUE)
  index <- which(trainingData > mean(trainingData))
  
  limit1 <- m
    
  limit2 <- (max(trainingData,na.rm = TRUE)-min(trainingData, na.rm=TRUE))/2 + min(trainingData, na.rm=TRUE)

  limits <- 0
  limits[1] <- limit1
  limits[2] <- limit2
  
  #print(limits)
  
  intt <- list()
  
  for(i in 1:2){
  
    temp <- 0
    tpath <- paste(path1,"/", i, sep = "")
    dir.create(tpath, showWarnings = FALSE)
    
    # INTERVALS IDENTIFIER
    temp <- intervalIdentifier(trainingData, limits[i])
    
    # INTERVALS SUMMARIZER
    temp <- intervalSummarizer(temp, 30, length(trainingData))
    write.csv(temp, paste(tpath, "/intervals.csv", sep = ""))
    
    # SUMMARY MATRIX
    synthesis <- summaryGenerator(trainingData, temp, limits[i])
    write.csv(synthesis, paste(tpath, "/synthesis.csv", sep = ""))
  
    intt[[i]] <- temp
    #print(intt[[i]])
    
  }
  
  return(intt)
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

simpleIntervalGenerator <- function(path1, trainingData){
  
  path2 <- paste(path1,"/simpleIntervals", sep = "")
  dir.create(path2, showWarnings = FALSE)
  
  ls <- length(trainingData)
  div <- divisors(ls)
  i <- 1
  temp <- div[i]
  qt <- ls
  
  while(qt > 10){
    
    temp <- div[i]
    qt <- ls/temp
    i <- i+1
    
  }
  
  temp <- div[i-1]
  div <- temp
  qt <- ls/div
  
  intervals <- 0
  
  if(length(unique(trainingData))>1){
      
      for( i in 1:qt){
        
        intervals[i] <- (i-1)*div + 1
        
      }
      
      intervals <- round(intervals)
      intervals[i+1] <- ls
      
      if(debugG == 1){ print("Model1")}
      
      path21 <- paste(path2, "/model1", sep = "")
      dir.create(path21, showWarnings = FALSE)
      mod1int <- model1(trainingData, intervals)
      write.csv(mod1int, paste(path21,"/intervals.csv", sep = ""))
      stateList <- stateListGenerator(mod1int, trainingData, mod1int, path21, 1)
      transition <- transitionMatrixBuilder(mod1int, trainingData, path21, 1)
      
      if(debugG == 1){ print("Model2")}
      
      path22 <- paste(path2, "/model2", sep = "")
      dir.create(path22, showWarnings = FALSE)
      mod2int <- model2(trainingData, intervals)
      write.csv(mod2int, paste(path22,"/intervals.csv", sep = ""))
      stateList <- stateListGenerator(mod2int, trainingData, mod2int, path22, 1)
      transition <- transitionMatrixBuilder(mod2int, trainingData, path22, 1)
      
      if(debugG == 1){ print("Model3")}
      
      path23 <- paste(path2, "/model3", sep = "")
      dir.create(path23, showWarnings = FALSE)
      mod3int <- model3(trainingData, intervals)
      write.csv(mod3int, paste(path23,"/intervals.csv", sep = ""))
      stateList <- stateListGenerator(mod3int, trainingData, mod3int, path23, 1)
      transition <- transitionMatrixBuilder(mod3int, trainingData, path23, 1)
      
      if(debugG == 1){ print("Model4")}

      mod4int <- list()
      path24 <- paste(path2, "/model4", sep = "")
      dir.create(path24, showWarnings = FALSE)
      mod4int <- model4(trainingData, intervals, path24)
      for(i in 1:2){

        tpath <- paste(path24,"/", i, sep = "")

        stateList <- stateListGenerator(mod4int[[i]], trainingData, mod4int[[i]], tpath, 1)
        transition <- transitionMatrixBuilder(mod4int[[i]], trainingData, tpath, 1)

      }
      
      if(debugG == 1){ print("SubIntervals")}
      
      path25 <- paste(path2, "/subIntervals", sep = "")
      dir.create(path25, showWarnings = FALSE)
      
      v <- sort(table(trainingData),decreasing=TRUE)[1]
      v <- as.numeric(names(v))
      # v <- data.frame(v)
      # if(ncol(v) == 2){
      # names(v) <- c("Var1", "Freq")
      # } else {
      #   names(v) <- "Var1"
      # }
      # vv <- 1
      # 
      # if(class(v$Var1) != "NULL"){
      #   
      #   vv <- as.numeric(as.matrix(v$Var1)[1])
      #   
      # }
      
      #if(length(trainingData[trainingData>vv]) > 250){
      if(length(trainingData[trainingData>v]) > 250){
        subIntervalsGenerator(path25, trainingData, mod4int[[1]])
        
      } else {
        
        write.csv(0, paste(path25, "/nl.csv", sep = ""))
        
      }
     
  } else {
    
    intervals <- c(1, length(trainingData))
    path20 <- paste(path2, "/model0", sep = "")
    dir.create(path20, showWarnings = FALSE)
    write.csv(intervals, paste(path20, "/intervals.csv", sep = ""))
    
    synthesis <- matrix(0, nrow = 3, ncol = 2)
    synthesis[1,] <- intervals
    synthesis[2,] <- c(0,0)
    synthesis[3,] <- trainingData[1]
    
    write.csv(synthesis, paste(path20, "/synthesis.csv", sep = ""))
    
    
    
  }
  
}