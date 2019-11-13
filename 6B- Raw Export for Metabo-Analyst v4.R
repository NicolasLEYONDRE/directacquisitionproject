## This script is intended to extract the adduct-corrected peak lists from the preVKs files in 
## order to send them to the MetaboAnalyst platform
## https://www.metaboanalyst.ca/faces/ModuleView.xhtml

library(readr)

# Set Directory
input_dir <- 'C:/Users/Simon Ollivier/Desktop/working directory/1- deisopeaklist'
setwd(input_dir)
input_files <- list.files(input_dir, pattern = "deisopeaklist[.]CSV")

# Select files
for(i in 1:length(input_files)){
  csv_active<-read_csv(input_files[[i]])
  
  basenames <- c("mz", "into")
  Xtract <- matrix(ncol=length(basenames), nrow = nrow(csv_active))
  colnames(Xtract) <- c(basenames)
  
  Xtract[,"mz"]=cbind(csv_active$"mz")
  Xtract[,"into"]=cbind(csv_active$"into")
  
  Xtractfile <- paste(substr(input_files[[i]],1,nchar(input_files[[i]])-17), "_rawforMA.csv",sep ="")
  write.csv(Xtract, file = Xtractfile, row.names = FALSE)
  }