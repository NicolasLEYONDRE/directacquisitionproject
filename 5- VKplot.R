                      #######################################################################
                      #   THE FIRST PART OF THIS SCRIPT IS SIMILAR TO CHEMQUANT,            #
                      #   IF YOU HAVE ALREADY GENERATED THE DATASET W/ CHEMQUANT            #
                      #   PLEASE LOAD THE .RDATA FILE AND GO DIRECTLY TO THE PLOTTING PART  #
                      #######################################################################

library(readr)
library(dplyr)
library(tidyr)
library(plot3D)
library(plot3Drgl)
library(ggplot2)
library(plotly)
                      
# Set Directory
input_dir <- 'C:/Users/Simon Ollivier/Desktop/working directory/3- preVKs'
setwd(input_dir)
temp <- list.files(input_dir, pattern = "preVKs[.]csv")

# Group files

list2env(lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), read.csv), envir = .GlobalEnv)
dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))
data <- bind_rows(dfs)
data$chem <- NA

# Compute chemical family according to Van Krevelen coordinates for all data
for(i in (1:nrow(data))){
  HC = data[i,"H.C"]
  OC = data[i,"O.C"]
  
  ## Correlates index w/ chemical family (rectangles approx.)
  
  ### Unsaturated hydrocarbons
  if((between(HC,0.85,1.5)==TRUE)&(between(OC,0,0.1)==TRUE)){
    data$chem[i] = "Unsaturated hydrocarbons"
  }
  
  ### Condensed aromatic compounds
  if((between(HC,0,0.65)==TRUE)&(between(OC,0,0.65)==TRUE)){
    data$chem[i] = "Condensed aromatic compounds"
  }
  if((between(HC,0,0.5)==TRUE)&(between(OC,0,0.7)==TRUE)){
    data$chem[i] = "Condensed aromatic compounds"
  }
  if((between(HC,0,0.45)==TRUE)&(between(OC,0,0.75)==TRUE)){
    data$chem[i] = "Condensed aromatic compounds"
  }
  if((between(HC,0,0.4)==TRUE)&(between(OC,0,0.85)==TRUE)){
    data$chem[i] = "Condensed aromatic compounds"
  }
  
  ### (Poly)phenolic compounds
  if((between(HC,0.65,0.75)==TRUE)&(between(OC,0.45,0.8)==TRUE)){
    data$chem[i] = "Phenolic compounds"
  }
  if((between(HC,0.75,0.85)==TRUE)&(between(OC,0.4,0.85)==TRUE)){
    data$chem[i] = "Phenolic compounds"
  }
  if((between(HC,0.85,1.0)==TRUE)&(between(OC,0.35,0.9)==TRUE)){
    data$chem[i] = "Phenolic compounds"
  }
  if((between(HC,1.0,1.15)==TRUE)&(between(OC,0.385,0.9)==TRUE)){
    data$chem[i] = "Phenolic compounds"
  }
  if((between(HC,1.15,1.30)==TRUE)&(between(OC,0.385,0.85)==TRUE)){
    data$chem[i] = "Phenolic compounds"
  }
  
  ### Terpenoids
  if((between(HC,0.85,1.0)==TRUE)&(between(OC,0.225,0.35)==TRUE)){
    data$chem[i] = "Terpenoids"
  }
  if((between(HC,1.0,1.15)==TRUE)&(between(OC,0.20,0.385)==TRUE)){
    data$chem[i] = "Terpenoids"
  }
  if((between(HC,1.15,1.5)==TRUE)&(between(OC,0.20,0.4)==TRUE)){
    data$chem[i] = "Terpenoids"
  }
  if((between(HC,1.5,1.6)==TRUE)&(between(OC,0.225,0.395)==TRUE)){
    data$chem[i] = "Terpenoids"
  }
  if((between(HC,1.6,1.75)==TRUE)&(between(OC,0.25,0.385)==TRUE)){
    data$chem[i] = "Terpenoids"
  }
  
  ### Nucleic acids
  if((between(HC,1.30,1.50)==TRUE)&(between(OC,0.40,1.05)==TRUE)){
    data$chem[i] = "Nucleic acids"
  }
  if((between(HC,1.50,1.60)==TRUE)&(between(OC,0.60,0.90)==TRUE)){
    data$chem[i] = "Nucleic acids"
  }
  
  ### Steroids
  if((between(HC,1.50,2.00)==TRUE)&(between(OC,0,0.10)==TRUE)){
    data$chem[i] = "Sterol-likes"
  }
  
  ## Fatty acids
  if((between(HC,1.85,2.00)==TRUE)&(between(OC,0.1,0.35)==TRUE)){
    data$chem[i] = "Fatty acids"
  }
  if((between(HC,2.00,2.30)==TRUE)&(between(OC,0,0.35)==TRUE)){
    data$chem[i] = "Fatty acids"
  }
  
  ## Small acids
  if((between(HC,1.5,1.6)==TRUE)&(between(OC,0.1,0.225)==TRUE)){
    data$chem[i] = "Small acids"
  }
  if((between(HC,1.6,1.85)==TRUE)&(between(OC,0.1,0.25)==TRUE)){
    data$chem[i] = "Small acids"
  }
  
  ## Amino acids
  if((between(HC,1.5,1.6)==TRUE)&(between(OC,0.395,0.60)==TRUE)){
    data$chem[i] = "Amino acids"
  }
  if((between(HC,1.6,1.75)==TRUE)&(between(OC,0.385,0.70)==TRUE)){
    data$chem[i] = "Amino acids"
  }
  if((between(HC,1.75,1.9)==TRUE)&(between(OC,0.32,0.79)==TRUE)){
    data$chem[i] = "Amino acids"
  }
  if((between(HC,1.9,2.6)==TRUE)&(between(OC,0.35,0.79)==TRUE)){
    data$chem[i] = "Amino acids"
  }
  
  ## Carbohydrates
  if((between(HC,1.50,1.60)==TRUE)&(between(OC,0.90,1.20)==TRUE)){
    data$chem[i] = "Carbohydrates"
  }
  if((between(HC,1.60,2.50)==TRUE)&(between(OC,0.79,1.20)==TRUE)){
    data$chem[i] = "Carbohydrates"
  }  
}

                  ###########################################################################
                  #  FROM THIS POINT ONWARDS THE CODE HAS NOT BEEN OPTIMISED AND IS ONLY    #
                  #  INTENDED FOR OUR EXPERIMENTAL CONDITIONS, PLEASE MAKE THE APPROPRIATE  #
                  #  CHANGES IF NECESSARY (e.g. "Acetone" condition)                        #
                  ###########################################################################

# Sort all the data into variables & coordinates for plotting

data$x <- NA
data$y <- NA

for (j in 1:nrow(data)){
  
  ## y axis
  if(grepl("Condensed aromatic compounds", data$chem[j])==TRUE){
    data$y[j] = 1
  }
  if(grepl("Phenolic compounds", data$chem[j])==TRUE){
    data$y[j] = 2
  }
  if(grepl("Terpenoids", data$chem[j])==TRUE){
    data$y[j] = 3
  } 
  if(grepl("Sterol-likes", data$chem[j])==TRUE){
    data$y[j] = 4
  }
  if(grepl("Unsaturated hydrocarbons", data$chem[j])==TRUE){
    data$y[j] = 5
  }
  if(grepl("Fatty acids", data$chem[j])==TRUE){
    data$y[j] = 6
  }
  if(grepl("Small acids", data$chem[j])==TRUE){
    data$y[j] = 7
  }
  if(grepl("Nucleic acids", data$chem[j])==TRUE){
    data$y[j] = 8
  }
  if(grepl("Amino acids", data$chem[j])==TRUE){
    data$y[j] = 9
  }
  if(grepl("Carbohydrates", data$chem[j])==TRUE){
    data$y[j] = 10
  }
  
  ## x axis
  if(grepl("DART", data$Ionisation[j])==TRUE){
    data$x[j] = 1
  }
  if(grepl("ASAP", data$Ionisation[j])==TRUE){
    data$x[j] = 2
  }
  if((grepl("ESI",data$Ionisation[j])==TRUE)&(grepl("Acétone", data$Condition[j])==TRUE)){
    data$x[j] = 3
  }
  if((grepl("ESI",data$Ionisation[j])==TRUE)&(grepl("Méthanol", data$Condition[j])==TRUE)){
    data$x[j] = 4
  }
}

######################################################################################################################
# Creating the plotting function

plotVanKrevelen <- function(data2){
  
  # Correlate cex to frequency of Van Krevelen coordinates
  width <- table(data2$O.C,data2$H.C)

scatter2D (data2$O.C, data2$H.C, colvar = data2$x,
           col = NULL, NAcol = "white", breaks = NULL,
           colkey = NULL, clim = NULL, clab = NULL, 
           alpha=0.25,
           pch = 20, cex = width+0.5,
           xlab="O/C ratio", ylab="H/C ratio",
           CI = NULL, add = FALSE, plot = TRUE)
}
######################################################################################################################

# Apply desired filters to your data before plotting, if no 'chem' filtering is desired, please skip the 'computing'
# part of the script to reduce calculation time

#test <- filter(data, Lichen=="Rf" & Polarité=="POS" & Nature=="Multi 1" & Ionisation=="ESI" & Condition=="Acétone")
test <- filter(data, Lichen=="Rf" & Polarité=="POS" & Nature=="Broyat" & Ionisation=="ASAP")
plotVanKrevelen(test)
