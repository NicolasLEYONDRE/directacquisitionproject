library(readr)
library(dplyr)
library(tidyr)
library(plot3D)
library(plot3Drgl)
library(ggplot2)
library(plotly)

# Set Directory
input_dir <- 'C:/Users/Simon Ollivier/Desktop/working directory'
setwd(input_dir)
temp <- list.files(input_dir, pattern = "preVKs[.]csv")

# Find the sample codes
#check<-unique(na.omit(as.numeric(unlist(strsplit(unlist(substr(temp,1,10)), "[^0-9]+")))))

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

# Remove rows without a chem tag
dataset1 <- subset(data, !is.na(data$chem)==TRUE)

##########################################################################################################################
# adapted from http://www.sthda.com/english/wiki/impressive-package-for-3d-and-4d-graph-r-software-and-data-visualization
hist3D_fancy<- function(x, y, dataset, break.func = c("Sturges", "scott", "FD"), breaks = NULL,
                        colvar = NULL, col="white", clab=NULL, phi = 25, theta = 45, ...){
  
  # Compute the number of classes for a histogram
  break.func <- break.func [1]
  if(is.null(breaks)){
    x.breaks <- switch(break.func,
                       Sturges = nclass.Sturges(x),
                       scott = nclass.scott(x),
                       FD = nclass.FD(x))
    y.breaks <- switch(break.func,
                       Sturges = nclass.Sturges(y),
                       scott = nclass.scott(y),
                       FD = nclass.FD(y))
  } else x.breaks <- y.breaks <- breaks
  
  # Cut x and y variables in bins for counting
  x.bin <- seq(min(x)-1, max(x)+0.5, length.out = x.breaks)
  y.bin <- seq(0, max(y), length.out = y.breaks)
  xy <- table(cut(x, x.bin), cut(y, y.bin))
  z <- prop.table(xy,1)
  
  xmid <- 0.5*(x.bin[-1] + x.bin[-length(x.bin)])
  ymid <- 0.5*(y.bin[-1] + y.bin[-length(y.bin)])
  oldmar <- par("mar")
  par (mar = par("mar") + c(0, 0, 0, 0.5))
  
  H3D<- hist3D(x = xmid, y = ymid, z = z, ...,
         zlim = c(0, 1), zlab = "frequency",
         xlab="", ylab="", bty= "g", 
         phi = phi, theta = theta,
         col=ramp.col(c("magenta1", "purple", "midnightblue")), colkey=FALSE,
         shade = 0.2, border = "black",
         d = 1, ticktype = "detailed", alpha = 1,
         axes=FALSE) 
  
  #Set limit for the axes
  max.x <- max(x)+0.5
  x.axis <- min(x):max.x
  min.x <- min(x)-1
  y.axis <- seq(0.5,9.5,by=1)
  min.y <- 0
  max.y <- 10
  z.axis <- seq(0, 1, by=0.2)
  min.z <- 0
  max.z <- 100

  # Draw the ticks on the axes
  ##X
  tick.start <- trans3d(x.axis-0.125, min.y, min.z, H3D)
  tick.end <- trans3d(x.axis-0.125, (min.y - 0.50), min.z, H3D)
  segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)
  
  ##Y
  tick.start <- trans3d(max.x, y.axis, min.z, H3D)
  tick.end <- trans3d(max.x + 0.05*max(x), y.axis, min.z, H3D)
  segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)
  
  ##Z
  tick.start <- trans3d(min.x, min.y, z.axis, H3D)
  tick.end <- trans3d(min.x, (min.y - 0.50), z.axis, H3D)
  segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)

  # Label the axes
  ##X
  labels <- c('DART-MS', 'ASAP-MS', 'ESI-MS \n in acetone', 'ESI-MS \n in methanol')
  label.pos <- trans3d(x.axis-0.5, (min.y - 1.85), min.z, H3D)
  text(label.pos$x, label.pos$y, labels=labels[min(x):max(x)], adj=c(0, NA), cex=0.75, las=1)
  
  ##Y
  labels <- c('Condensed aromatic compounds', 'Phenolic compounds', 'Terpenoids', 'Sterol-like compounds',
              'Unsaturated hydrocarbons', 'Fatty acids', 'Small acids', 'Nucleic acids', 'Amino-acids', 'Carbohydrates')
  label.pos <- trans3d((max.x + 0.05*max(x) + 0.05), y.axis, min.z, H3D)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), cex=0.7)
  
  ##Z
  labels <- c('0','20','40','60','80','100')
  label.pos <- trans3d(min.x, (min.y - 0.5), z.axis, H3D)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(1, NA), cex=0.75)
  
  ##Z axis: title
  z.title <- c('Frequency (%)')
  z.title.pos <- trans3d(min.x, (min.y-1.7), 0.75, H3D)
  text(z.title.pos$x, z.title.pos$y, labels=z.title, adj=c(1, NA), cex=1, srt=102)
  
  
  # Heatmap
  if(is.null(breaks)){
    y.breaks2 <- switch(break.func,
                       Sturges = nclass.Sturges(dataset$y),
                       scott = nclass.scott(dataset$y),
                       FD = nclass.FD(dataset$y))
  } else y.breaks2 <- breaks+1
  y.bin2 <- seq(0, max(dataset$y), length.out = 1.75*y.breaks2)
  Dx <- table(cut(dataset$x,x.bin),cut(dataset$y,y.bin2))
  ymid2 <- 0.5*(y.bin2[-1] + y.bin2[-length(y.bin2)])

  par(mar=c(4,4,4,4))
image3D(z = 1, x=xmid, y=ymid2, colvar = Dx, add = TRUE,
         colkey = list(length = 0.5, width = 0.4, shift = 0,
         cex.axis = 0.8, cex.clab = 0.8), 
         alpha =0.4,
         clab = c("","Number of ionized molecules","(Absolute value)"), 
         plot = TRUE)
  
  par(mar = oldmar)
}
##########################################################################################################################
#ep_pos_t <- filter(dataset1, Polarité=="POS" & (Nature=="Thalle"|Nature=="Broyat"))
ep_pos_t <- filter(dataset1, Polarité=="NEG" & Nature=="Multi 1")
#histo<- 
  hist3D_fancy(ep_pos_t$x, ep_pos_t$y, ep_pos_t, breaks = 11)
#plotrgl(histo)
