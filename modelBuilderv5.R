#_________________________________________________________________________________________________#
#                                                                                                 #
#                                 MODEL GENERATOR FRAMEWORK                                       #
#                                                                                                 #
#                                     Francesco Di Lillo                                          #
#                             Karlsruher Institut für Technologie                                 #
#                       Scuola Politecnica e delle Scienze di Base Federico II                    #
#_________________________________________________________________________________________________#


# Utility
#----------------------------#

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

# External Sources
#----------------------------#
source("custompath/macroIntervalmodel.R")
source("custompath/simpleIntervalsmodel.R")

# Fixed Variables
#----------------------------#

debugG <- 1

# Script
#----------------------------#

sensors <- c(12, 32, 32, 12, 32, 12, 32, 32, 32, 32)
topFolder <- paste("v5_ModelBuilder", Sys.Date(), sep = "")
dir.create(topFolder, showWarnings = FALSE)
startName <- "path"

for(mac in 1:10){
  
  machine <- mac
  nzero = "0"
  if(mac == 10){nzero = ""}
  machineName <- paste("AVT_",nzero,machine, sep = "")
  fileName <- paste(startName,"/",machineName,".csv", sep = "")
  nome <- paste(topFolder,"/Sim", machine, sep = "")
  dir.create(nome, showWarnings = FALSE)
  
  if(debugG == 1){print(paste("Machine ", machine, sep = ""))}
  
  data <- data.matrix(read.csv(fileName))                          # Historical data
  data <- data[-1,-1]

  rows <- ncol(data)
  cols <- nrow(data)
  time <- cols
  
  # TREND MODELING
  
  trendTs <- tsTrendExtraction(data)
  debug <- 1
  dim <- nrow(trendTs)
  
  for(i in 1:dim){
    
    if(debugG == 1) {print(paste("Sensor ", i, sep = ""))}
    
    ind <- i
    
    path1 <- paste(nome,"/",ind, sep = "")
    dir.create(path1, showWarnings = FALSE)
    
    subject <- trendTs[ind,][!is.na(trendTs[ind,])]
    
    macroIntervalGenerator(path1, subject)
    simpleIntervalGenerator(path1, subject)
    
  }
  
}