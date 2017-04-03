#_________________________________________________________________________________________________#
#                                                                                                 #
#                                       FILE FILTER v3                                            #
#                                                                                                 #
#                                     Francesco Di Lillo                                          #
#                             Karlsruher Institut für Technologie                                 #
#                       Scuola Politecnica e delle Scienze di Base Federico II                    #
#_________________________________________________________________________________________________#

# Libraries
#----------------------------#
library(data.table)

#----------------------------#

# Utility
#----------------------------#

# ///

# Script
#----------------------------#


path1 <- "C:/Users/dilillo/Desktop/R/AVT_export/export-all-2017-02-28.csv"
limit <- c(12, 32, 32, 12, 32, 12, 32, 32, 32, 32)

path <- path1
AVT_export <- read.csv(path, sep=",", header = FALSE)
table1 <- as.data.table(AVT_export)

for(m in 1: 10){
  
  machine <- paste("AVT_0",m,sep = "")
  
  if(m == 10){
    
    machine <- "AVT_10"
    
  }
  
  print(machine)
  
  preFiltered <- table1[table1$V1 == machine]
  
  result <- preFiltered
  
  pathZeros <- "00"
  finalDF <- data.frame(0, 0)
  
  for(j in 1:limit[m]){
    
    tempStore <- data.frame(0, 0, 0, 0)
    names(tempStore) <- c("V1", "V2", "V3", "V4")
    
    if(j == 10){
      
      pathZeros <- "0"
      
    }
    
    sensor <- paste(pathZeros,j, sep ="")    
    
    print(sensor)
    
    cond <- FALSE
    for(o in 1:nrow(result)){
      
      temp <- toString(result$V2[o])
      temp <- substring(temp, 1,3)
      
      if(temp == sensor){
        
        
        if(ncol(tempStore) == 0){
          
          tempStore <- data.frame(toString(result$V1[o]), toString(result$V2[o]), toString(result$V3[o]),
                                  toString(result$V4[o]))
          names(tempStore) <- c("V1", "V2", "V3", "V4")
          
        }else{
          
          tST <- data.frame(toString(result$V1[o]), toString(result$V2[o]), toString(result$V3[o]),
                            toString(result$V4[o]))
          names(tST) <- c("V1", "V2", "V3", "V4")
          tempStore <- rbind(tempStore,tST)
          
        }
      }
      
      
    }
    
    
    order.tempStore <- tempStore[order(tempStore$V3),]
    
    
    if(j == 1){
      
      finalDF <- data.frame(order.tempStore$V3, order.tempStore$V4)
      names(finalDF) <- c(paste("timestampSensor ",j, sep = ""), paste("Sensor ",j, sep = ""))
      
    } else{
      
      tempDF <- data.frame(order.tempStore$V3, order.tempStore$V4)
      names(tempDF) <- c(paste("timestampSensor ",j, sep = ""), paste("Sensor ",j, sep = ""))
      finalDF <- data.frame(finalDF, tempDF)
      
    }
    
    print(paste("Fine elab ", sensor, sep = ""))
    
  }
  
  write.csv(finalDF, paste("C:/Users/dilillo/Desktop/R/Dataset 28-02/",machine,".csv", sep = ""))
  
}