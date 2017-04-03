#_________________________________________________________________________________________________#
#                                                                                                 #
#                                 SYNTHETIC TREND GENERATOR                                       #
#                                                                                                 #
#                                     Francesco Di Lillo                                          #
#                             Karlsruher Institut für Technologie                                 #
#                       Scuola Politecnica e delle Scienze di Base Federico II                    #
#_________________________________________________________________________________________________#


# Utility
#----------------------------#

macroGenerator <- function(mStateList, mHistogram, mTransitions){
  
  densities <- mHistogram$density
  densities <- densities[1:(length(densities)-1)]
  syntheticDataset <- 0
  
  s <- mStateList[which(mStateList>0)]
  prob <- densities/sum(densities)
  
  if(length(s) > 1){
  
    currentState <- sample(s, 1, prob = prob, replace = TRUE)
    
    syntheticDataset[1] <- mHistogram$breaks[currentState + 1]
    
    
    for(j in 2:time){
      
      probNextState <- mTransitions[currentState,]
      
      if(any(is.na(probNextState))){
        npezz <- length(probNextState)
        probNextState <- rep(1/npezz, npezz)
        
      }
      
      currentState <- sample(s, 1, prob=probNextState, replace = TRUE)
      
      if(length(s) > 1){
        
        extrA <- mHistogram$breaks[currentState]
        extrB <- mHistogram$breaks[currentState + 1]
        nextSeq <- seq(from = extrA, to = extrB, length.out = 4)
        nextValue <- sample(nextSeq, 1, replace = FALSE)
        
      }
      
      else {
        
        nextValue <- mHistogram$breaks[currentState + 1]
        
      }
      
      syntheticDataset[j] <- nextValue
      
    }
  
  } else {
    
    for(i in 1:time){
      
      syntheticDataset[i] <- mHistogram$breaks[2]
      
    }
    
  }
  
  return(syntheticDataset)
  
}
hLoader <- function(splitIntervalFpath, li, op, macroState){
  
  histograms <- list()
  dim <- li-1
  
  if(op == 1){
    for(i in 1:dim){
      
      histograms[[i]] <- read.csv(paste(splitIntervalFpath,"/histogram",i,".csv", sep = ""))
      histograms[[i]] <- histograms[[i]][, c(2,3)]
      names(histograms[[i]]) <- c("breaks", "density")
      
      
    }
  } else {
    
    for(j in 1:dim){
      
      hPath <- paste(splitIntervalFpath, "/State",macroState,"histogram",j,".csv", sep = "")
      
      histograms[[j]] <- read.csv(hPath)
      histograms[[j]] <- histograms[[j]][, c(2,3)]
      names(histograms[[j]]) <- c("breaks", "density")
      
    }
    
  }
  
  return(histograms)
  
}
tLoader <- function(splitIntervalFpath, li, op, macroState){
  
  transitions <- list()
  dim <- li-1
  
  if(op == 1){
  
    for(i in 1:dim){
      
      m <- data.matrix(read.csv(paste(splitIntervalFpath,"/transitions",i,".csv", sep = "")))
      m <- m[, c(2: ncol(m))]
      #m <- probModifier(m)
      
      transitions[[i]] <- m
      
    }
    
  } else {
    
    for(j in 1:dim){
      
      tPath <- paste(splitIntervalFpath, "/State",macroState,"subState",j,"transitions.csv", sep = "")
      
      m <- data.matrix(read.csv(tPath))
      m <- m[, c(2: ncol(m))]
      ddeb <- m
      
      transitions[[j]] <- m
      
    }
    
  }
  
  return(transitions)
  
}
splitterGenerator <- function(sStateList, sIntervals, splitIntervalFpath){
  
  li <- length(sIntervals)
  syntheticDataset <- 0
  histograms <- list()  
  transitionMatrix <- list()
  
  histograms <- hLoader(splitIntervalFpath, li, 1, 1)
  transitionMatrix <- tLoader(splitIntervalFpath, li, 1, 1)
  
  currentState <- 0
  switch <- FALSE
  macroState <- 1
  model <- 1
    
    for(i in 1:time){
      
      cond <- i%%time
      
      if(cond == 0) {
        cond <- time
      }
      
      tModel <- which(sIntervals <= cond)
      tModel <- tModel[length(tModel)]
      
      if(tModel != model){
        
        model <- tModel
        switch <- TRUE
        
        if(model == length(sIntervals)){
          model <- length(sIntervals) - 1
        }
        
        #if(debugG == 1) {print(model)}
        
      }
      
      densities <- histograms[[model]]$density
      densities <- densities[1:(length(densities)-1)]
      transitions <- transitionMatrix[model]
      transitions <- as.matrix(transitions[[1]])
      
      s <- sStateList[model,][which(sStateList[model,]>0)]
      prob <- densities/sum(densities)
      
      if(i == 1){
        
        currentState <- as.numeric(sample(s, 1, prob = prob, replace = TRUE))
        syntheticDataset[i] <- histograms[[model]]$breaks[currentState + 1]
        
      } else {
        
        if(switch == TRUE){
          
          if(currentState > ncol(transitions)) {
            
            currentState <- as.numeric(sample(s, 1, prob = prob, replace = TRUE))  
            
          }
          
          switch <- FALSE
          
        }
        
        probNextState <- transitions[currentState,]
        currentState <- as.numeric(sample(s, 1, prob=probNextState, replace = TRUE))
        
        extrA <- histograms[[model]]$breaks[currentState]
        extrB <- histograms[[model]]$breaks[currentState + 1]
        nextSeq <- seq(from = extrA, to = extrB, length.out = 4)
        nextValue <- sample(nextSeq, 1, replace = FALSE)
        
        syntheticDataset[i] <- nextValue
        
      }
      
    }
  
  return(syntheticDataset)
  
}
basicGenerator <- function(synthesis){
  
  intervals <- as.numeric(synthesis[1,])
  syntheticDataset <- 0
  macroState <- 1
  cond2 <- FALSE
  
  for(i in 1:time){
    
    cond <- i%%time
    if(cond == 0) {
      cond <- time
    }
    
    newMState <- which(intervals <= cond)
    newMState <- as.numeric(newMState[length(newMState)])
    
    if(newMState != macroState){
      
      macroState <- newMState
      
      if(macroState == length(intervals)){
        
        macroState <- length(synthesis[1,]) - 1
        
      }
      
    }
    
    syntheticDataset[i] <- synthesis[3,macroState]
    
    
  }
  
  return(syntheticDataset)
  
}
subStateGenerator <- function(synthesisIntervalFpath, synthesisSubIntervalpath){
  
  macroState <- 1
  subState <- 1
  start <- TRUE
  randomPick <- 0
  syntheticDataset <- 0
  
  sPath <- paste(synthesisIntervalFpath,"/synthesis.csv", sep = "")
  
  sLpath <- paste(synthesisSubIntervalpath, "/subStateList", macroState,".csv", sep = "")
  sSiPath <- paste(synthesisSubIntervalpath, "/subIntervalsMatrix.csv", sep = "")
  
  synthesis <- read.csv(sPath)
  synthesis <- data.matrix(synthesis[,c(2:ncol(synthesis))])
  
  subIntervalsMatrix <- read.csv(sSiPath)
  subIntervalsMatrix <- data.matrix(subIntervalsMatrix[,c(2:ncol(subIntervalsMatrix))])
  
  subStateList <- read.csv(sLpath)
  subStateList <- data.matrix(subStateList[,c(2:ncol(subStateList))])
  
  intervals <- as.numeric(synthesis[1,])
  subIntervals <- as.numeric(subIntervalsMatrix[macroState,])
  subIntervals <- subIntervals[subIntervals > 0]
  ls <- length(subIntervals) 
  
  stateList <- read.csv(sLpath)
  stateList <- data.matrix(stateList[,c(2:ncol(stateList))])
  
  histograms <- list()
  transitionMatrix <- list()
  histograms <- hLoader(synthesisSubIntervalpath, ls, 2, macroState)
  transitionMatrix <- tLoader(synthesisSubIntervalpath, ls, 2, macroState)
  print(c(length(histograms), length(transitionMatrix)))
  
  counter <- 0
  
  for(i in 1:time){
    
    cond <- i%%intervals[length(intervals)]
    if(cond == 0) {
      cond <- intervals[length(intervals)]
    }
    
    newMState <- which(intervals <= cond)
    newMState <- as.numeric(newMState[length(newMState)])
    
    if(newMState != macroState){
      
      switch <- TRUE
      start <- TRUE
      macroState <- newMState
      subState <- 1
      
      if(macroState == length(intervals)){
        
        macroState <- length(intervals) - 1
        
      }
      
      sLpath <- paste(synthesisSubIntervalpath, "/subStateList", macroState,".csv", sep = "")
      
      subStateList <- read.csv(sLpath)
      subStateList <- data.matrix(subStateList[,c(2:ncol(subStateList))])
      subIntervals <- as.numeric(subIntervalsMatrix[macroState,])
      subIntervals <- subIntervals[subIntervals > 0]
      ls <- length(subIntervals)
      
      histograms <- list()
      transitionMatrix <- list()
      histograms <- hLoader(synthesisSubIntervalpath, ls, 2, macroState)
      transitionMatrix <- tLoader(synthesisSubIntervalpath, ls, 2, macroState)
      print(c(length(histograms), length(transitionMatrix)))
    }
    
    else {
     
      newSState <- which((subIntervals+intervals[macroState] - 1) <= cond)
      
      newSState <- as.numeric(newSState[length(newSState)])
      
      if(newSState != subState){
        
        if(randomPick!=1){
          subState <- newSState
        } else {
          
          subState <- sample(1:length(subIntervals), 1, replace=TRUE)
          
        }
        
        switch <- TRUE
        
        if(subState == length(subIntervals)){
          
          subState <- length(subIntervals) - 1
          
        }
        
        
        
      }
      
    }
    
    densities <- histograms[[subState]]$density
    densities <- densities[1:(length(densities)-1)]
    transitions <- transitionMatrix[[subState]]
    transitions <- as.matrix(transitions)
    
    s <- as.numeric(subStateList[subState,][which(subStateList[subState,]>0)])
    prob <- densities/sum(densities)
    
    if(start == TRUE){
      
      start <- FALSE
      switch <- FALSE
      maxP <- max(densities)
      currentState <- which(histograms[[subState]]$density == maxP)[1]
      syntheticDataset[i] <- histograms[[subState]]$breaks[currentState[1] + 1]
      
      #print(c("START",i))
      #print(c(macroState, subState))
      
    } else {
      
      if(switch == TRUE){
        
        counter <- counter+1
        
        if(currentState > nrow(transitions)) {
          
          maxP <- max(densities)
          currentState <- which(histograms[[subState]]$density == maxP)
          
        }
        
        
        #print(c("SWITCH",i))
        switch <- FALSE
        
      }
      
      
      probNextState <- transitions[currentState,]
      if(class(probNextState) == "matrix"){ if(nrow(probNextState)>1){probNextState <- probNextState[1,]}}
      currentState <- as.numeric(sample(s, 1, prob=probNextState, replace = TRUE))
      
      extrA <- histograms[[subState]]$breaks[currentState[1]]
      extrB <- histograms[[subState]]$breaks[currentState[1] + 1]
      
      if(length(histograms[[subState]]$breaks) == 2){
        nextSeq <- extrA
      }
      else {
        nextSeq <- seq(from = extrA, to = extrB, length.out = 4)
      }
      
      
      nextValue <- sample(nextSeq, 1, replace = FALSE)
      
      syntheticDataset[i] <- nextValue
      
    }
    
  }
  
  return(syntheticDataset)
  
}


# External Sources
#----------------------------#

source("evaluationFunction.R")

# Fixed Variables
#----------------------------#

time <- 9755
debugG <- 1

# Script
#----------------------------#

finalDataset <- list()
evaluationSet <- list()
topFolder <- paste("v5_ModelBuilder", Sys.Date(), sep = "")
machinenumber <- length(list.dirs(path = topFolder, full.names = TRUE, recursive = FALSE))

for(machine in 1:machinenumber){
  
  if(debugG == 1){print(paste("Machine ", machine, sep = ""))}
  
  simFolder <- paste(topFolder,"/Sim", machine, sep = "")
  subFolderlength <- length(list.dirs(path = simFolder, full.names = TRUE, recursive = FALSE))
  bDS <- matrix(0, nrow = subFolderlength, ncol = time)
  
  for(sensors in 1:subFolderlength){
    
    if(debugG == 1){print(paste("Sensor ", sensors, sep = ""))}
    
    sensorFolder <- paste(simFolder, "/", sensors, sep = "")
    splitIntervalpath <- paste(sensorFolder,"/simpleIntervals", sep = "")
    model0 <- paste(splitIntervalpath,"/model0", sep = "")
    
    if(file.exists(model0)){
      
      if(debugG == 1){print("MODEL 0")}
      
      synthesisIntervalpath <- paste(splitIntervalpath,"/model0/synthesis.csv", sep = "")
      synthesis <- read.csv(synthesisIntervalpath)
      synthesis <- data.matrix(synthesis[,c(2:ncol(synthesis))])
      
      bDS[sensors,] <- basicGenerator(synthesis)
      if(debugG == 1){plot(bDS[sensors,], type = "l", main = "Model 0", ylab = paste("Best Dataset[",sensors,"]", sep = ""))}
      
    } else{
    
    syntheticDataset <- matrix(0, nrow = 10, ncol = time)
    k <- 1
    
    # MACRO INTERVAL GENERATOR
    if(debugG == 1){print(c("MACRO INTERVAL", k))}
    macroIntervalpath <- paste(sensorFolder,"/macroInterval", sep = "")
    mTransitions <- read.csv(paste(macroIntervalpath, "/transitions.csv", sep = ""))
    mTransitions <- data.matrix(mTransitions[,c(2:ncol(mTransitions))])
    mStateList <- read.csv(paste(macroIntervalpath, "/stateList.csv", sep = ""))
    mStateList <- mStateList$x
    mHistogram <- read.csv(paste(macroIntervalpath, "/histogram.csv", sep = ""))
    mHistogram <- mHistogram[, c(2,3)]
    names(mHistogram) <- c("breaks", "density")
    
    syntheticDataset[k,] <- macroGenerator(mStateList, mHistogram, mTransitions)
    
    if(debugG == 1){plot(syntheticDataset[1,], type = "l", main = "Macro Interval", ylab = "syntheticDataset")}
    k <- k+1
  
    # SPLITTER CRITERIA INTERVAL GENERATOR
    for(i in 1:3){
      
      model <- paste("/model", i, sep = "")
      if(debugG == 1){print(c(paste("MODEL ", i)), k)}
      splitIntervalFpath <- paste(splitIntervalpath, model, sep = "")
      sStateList <- read.csv(paste(splitIntervalFpath, "/stateList.csv", sep = ""))
      sStateList <- data.matrix(sStateList[,c(2:ncol(sStateList))])
      sIntervals <- read.csv(paste(splitIntervalFpath, "/intervals.csv", sep = ""))
      sIntervals <- sIntervals$x[!is.na(sIntervals$x)]
      
      syntheticDataset[k,] <- splitterGenerator(sStateList, sIntervals, splitIntervalFpath)
      
      if(debugG == 1){plot(syntheticDataset[k,], type = "l", main = paste("Model ",i,sep =""), ylab = paste("syntheticDataset[",k,"]", sep = ""))}
      
      k <- k+1
    }
    
    # SYNTHESIS APPROACHES

    for(i in 1:2){

      if(debugG == 1){print(c(paste("MODEL 4 -", i), k))}
      
      synthesisIntervalpath <- paste(splitIntervalpath,"/model4", sep = "")
      criteria <- i
      synthesisIntervalFpath <- paste(synthesisIntervalpath, "/", criteria, sep = "")
      synthesis <- read.csv(paste(synthesisIntervalFpath, "/synthesis.csv", sep = ""))
      synthesis <- data.matrix(synthesis[,c(2:ncol(synthesis))])

      syntheticDataset[k,] <- basicGenerator(synthesis)
      if(debugG == 1){plot(syntheticDataset[k,], type = "l", main = paste("Model ", k,sep = ""), ylab = paste("syntheticDataset[",k,"]", sep = ""))}

      synStateList <- read.csv(paste(synthesisIntervalFpath, "/stateList.csv", sep = ""))
      synStateList <- data.matrix(synStateList[,c(2:ncol(synStateList))])
      synIntervals <- read.csv(paste(synthesisIntervalFpath, "/intervals.csv", sep = ""))
      synIntervals <- synIntervals$x[!is.na(synIntervals$x)]

      syntheticDataset[k+1,] <- splitterGenerator(synStateList, synIntervals, synthesisIntervalFpath)
      if(debugG == 1){plot(syntheticDataset[k+1,], type = "l", main = paste("Model ", k+1 ,sep = ""), ylab = paste("syntheticDataset[",k+1,"]", sep = ""))}

      k <- k+2


    }
    
    nlpath <- paste(splitIntervalpath,"/subIntervals/nl.csv", sep = "")
      
    if(!file.exists(nlpath)){
      
      for(i in 1:2){
      
        if(debugG == 1){print(paste("SUBINTERVALS ", i))}
          
        synthesisSubIntervalpath <- paste(splitIntervalpath,"/subIntervals/", i, sep = "") 
        syntheticDataset[k,] <- subStateGenerator(paste(synthesisIntervalpath, "/1", sep = ""),synthesisSubIntervalpath)
        if(debugG == 1){plot(syntheticDataset[k,], type = "l", main = "SubIntervals", ylab = paste("syntheticDataset[",k,"]", sep = ""))}
        k <- k+1
      
      }
      
    } else {
      
      temp <- syntheticDataset[c(1:k),]
      syntheticDataset <- matrix(0, nrow = k, ncol = time)
      syntheticDataset <- temp
      
    }
    
    print("EVALUATION")
    
    evalFramework <- matrix(0, nrow = nrow(syntheticDataset), ncol = 6)
    bestACF <- 0
    bestAmp <- 0
    
    bInd <- -1
    for(j in 1:nrow(syntheticDataset)){
      
      evalFramework[j,] <- evaluateStep1(syntheticDataset[j,], sensors, simFolder)
      
    }
    
    evaluationSet[[sensors]] <- evalFramework
    
    oRthreshold <- sort(evalFramework[,6])
    if(length(oRthreshold) >= 7){
    threshold <- oRthreshold[7]}
    else{
      threshold <- oRthreshold[length(oRthreshold)]
    }
    startInd <- which(evalFramework[,6]==threshold)
    startInd <- startInd[1]
    bestStats <- rep(99, 6)
    
    for(j in 1:nrow(evalFramework)){
      
      print(j)
      
        if(j == 1){
          
          bestStats <- evalFramework[startInd,]
          bInd <- startInd
          bDS[sensors, ] <- syntheticDataset[startInd, ]
          
        } else {
          
          if(j != startInd && j != bInd){
         
            if(evalFramework[j,6] <= threshold){
            
            cond <- FALSE
            k <- 5
            
            while(cond == FALSE && k > 0){
              
              print(c(k, abs(evalFramework[j, k]), abs(bestStats[k])))
              
              if(abs(evalFramework[j, k]) < abs(bestStats[k])){
                
                bestStats <- evalFramework[j, ]
                bInd <- j
                bDS[sensors, ] <- syntheticDataset[j, ]
                cond <- TRUE
                
              } else if( abs(evalFramework[j, k]) == abs(bestStats[k]) ){
                
                k <- k-1
                
              } else {
                
                cond <- TRUE
                
              }
              
            } 
            
            
            }
            
          
          }
          
        
      }
      
      
    }
     
    if(debugG == 1){
      
      print(evalFramework[bInd,])
      plot(bDS[sensors,], type = "l")
      
    }
     
    }
    
  }
     
    finalDataset[[machine]] <- bDS
  
}