
library(xcms)
library(msdata)
library(MSnbase)
library(CAMERA)
library(readr)
library(dplyr)
library(tidyr)
library(plot3D)
library(plot3Drgl)
library(ggplot2)
library(plotly)

                      


# Set Directory

#input_dir <- # your path beware with / 
	input_dir <- 'C:/Users/CRMPO/Documents/2022/Simon.rapport.DATA/VK/TSV'
	setwd(input_dir)
	temp <- list.files(input_dir, pattern = "dataset_solvents_comp_publi1.tsv")
	temp2 <- list.files(input_dir, pattern = "DOJ.DB.MOD.simp.csv")
#reading function for tsv, csv
#take care of special characters like accent (we are French)
	data <- read.csv2(file = temp, sep = '\t', header = TRUE)
	DB.DOJ <- read.csv2(file = temp2, sep = ';', header = TRUE)
	data <- (rbind(data, DB.DOJ)	
#depending on file charateristics
#data <- read.csv2(file = temp, sep = ';', header = TRUE)
#data <- read.csv2(file = temp, sep = ',', header = TRUE)
#if the file is from OSF, it should be processed and saved with the right (chem, x, y) => goto "sort"

#Intitialisation of chem, y and x
	data$chem <- NA
	data$y <- NA
	data$x <- NA
# Compute chemical family according to Van Krevelen coordinates for all data
# one run for all 3
	for(i in (1:nrow(data)))
	{
	  HC = data[i,"H.C"]
	  OC = data[i,"O.C"]  
	  ## Correlates index w/ chemical family (rectangles approx.)
		### Condensed aromatic compounds y <- 1
	  if((between(HC,0.2,0.65)==TRUE)&(between(OC,0,0.65)==TRUE)){
		data$chem[i] = "Condensed aromatic compounds"
		data$y[i] = 1
	  } 
	   if((between(HC,0.65,0.85)==TRUE)&(between(OC,0,0.40)==TRUE)){
		data$chem[i] = "Condensed aromatic compounds"
		data$y[i] = 1
	  }  
	  if((between(HC,0.65,0.75)==TRUE)&(between(OC,0.40,0.45)==TRUE)){
		data$chem[i] = "Condensed aromatic compounds"
		data$y[i] = 1
	  }   

	 ### Polyphenols and derivatives y <- 2
	   if((between(HC,0.75,0.85)==TRUE)&(between(OC,0.40,0.45)==TRUE)){
		data$chem[i] = "Polyphenols and derivatives"
		data$y[i] = 2
	  }   
	   if((between(HC,0.65,0.75)==TRUE)&(between(OC,0.45,0.8)==TRUE)){
		data$chem[i] = "Polyphenols and derivatives"
		data$y[i] = 2
	  }
	  if((between(HC,0.75,0.85)==TRUE)&(between(OC,0.42,0.85)==TRUE)){
		 data$chem[i] = "Polyphenols and derivatives"
		 data$y[i] = 2
	  }
	  if((between(HC,0.85,1.0)==TRUE)&(between(OC,0.35,1)==TRUE)){
		data$chem[i] = "Polyphenols and derivatives"
		data$y[i] = 2
	  }
	  if((between(HC,1.0,1.15)==TRUE)&(between(OC,0.385,0.9)==TRUE)){
		data$chem[i] = "Polyphenols and derivatives"
		data$y[i] = 2
	  }
	  if((between(HC,1.15,1.30)==TRUE)&(between(OC,0.4,0.85)==TRUE)){
		data$chem[i] = "Polyphenols and derivatives"
		data$y[i] = 2
	  }
	  
	  ### Benzenoids or atypic terpenoids y <- 3
	  if((between(HC,0.85,1.0)==TRUE)&(between(OC,0.225,0.35)==TRUE)){
		data$chem[i] = "Benzenoids or atypic terpenoids"
		data$y[i] = 3
	  }
	  if((between(HC,1.0,1.15)==TRUE)&(between(OC,0.20,0.385)==TRUE)){
		data$chem[i] = "Benzenoids or atypic terpenoids"
		data$y[i] = 3
	  }
	  if((between(HC,1.15,1.5)==TRUE)&(between(OC,0.20,0.4)==TRUE)){
		data$chem[i] = "Benzenoids or atypic terpenoids"
		data$y[i] = 3
	  }
	  if((between(HC,1.5,1.6)==TRUE)&(between(OC,0.225,0.395)==TRUE)){
		data$chem[i] = "Benzenoids or atypic terpenoids"
		data$y[i] = 3
	  }
	  if((between(HC,1.6,1.75)==TRUE)&(between(OC,0.25,0.385)==TRUE)){
		data$chem[i] = "Benzenoids or atypic terpenoids"
		data$y[i] = 3
	  }
		 
	  ### Unsaturated hydrocarbons y <- 5
	  if((between(HC,0.85,1.5)==TRUE)&(between(OC,0,0.1)==TRUE)){
		data$chem[i] = "Unsaturated hydrocarbons"
		data$y[i] = 5
	  }
		
	 ### Terpenes y <- 4
	  if((between(HC,1.50,2.00)==TRUE)&(between(OC,0,0.10)==TRUE)){
		data$chem[i] = "Terpenes"
		data$y[i] = 4
	  }  
	  
	  ## Fatty acyls y <- 6
	  if((between(HC,1.75,1.8)==TRUE)&(between(OC,0.255,0.32)==TRUE)){
		data$chem[i] = "Fatty acyls"
		data$y[i] = 6
	  }
	  if((between(HC,1.8,2.00)==TRUE)&(between(OC,0.1,0.35)==TRUE)){
		data$chem[i] = "Fatty acyls"
		data$y[i] = 6
	  }
	  if((between(HC,2.00,2.30)==TRUE)&(between(OC,0,0.35)==TRUE)){
		data$chem[i] = "Fatty acyls"
		data$y[i] = 6
	  }
	  
	  ## Prenol derivatives y <- 7
	  if((between(HC,1.5,1.6)==TRUE)&(between(OC,0.1,0.225)==TRUE)){
		data$chem[i] = "Prenol derivatives"
		data$y[i] = 7
	  }
	  if((between(HC,1.6,1.85)==TRUE)&(between(OC,0.1,0.25)==TRUE)){
		data$chem[i] = "Prenol derivatives"
		data$y[i] = 7
	  }  
	  
	  ### Nucleic acids y<- 8
	  if((between(HC,1.3,1.50)==TRUE)&(between(OC,0.40,1.05)==TRUE)){
		data$chem[i] = "Nucleic acids"
		data$y[i] = 8
	  }
	  if((between(HC,1.50,1.60)==TRUE)&(between(OC,0.60,0.90)==TRUE)){
		data$chem[i] = "Nucleic acids"
		data$y[i] = 8
	  } 
	 
	  ## Amino acids y <- 9
	  if((between(HC,1.5,1.6)==TRUE)&(between(OC,0.395,0.60)==TRUE)){
		data$chem[i] = "Amino acids"
		data$y[i] = 9
	  }
	  if((between(HC,1.6,1.75)==TRUE)&(between(OC,0.385,0.70)==TRUE)){
		data$chem[i] = "Amino acids"
		data$y[i] = 9
	  }
	  if((between(HC,1.75,1.9)==TRUE)&(between(OC,0.32,0.82)==TRUE)){
		data$chem[i] = "Amino acids"
		data$y[i] = 9
	  }
	  if((between(HC,1.9,2.6)==TRUE)&(between(OC,0.35,0.82)==TRUE)){
		data$chem[i] = "Amino acids"
		data$y[i] = 9
	  }
	  
	  ## Carbohydrates y <- 10
	  if((between(HC,1.50,1.60)==TRUE)&(between(OC,0.90,1.20)==TRUE)){
		data$chem[i] = "Carbohydrates"
		data$y[i] = 10
	  }
	  if((between(HC,1.60,2.50)==TRUE)&(between(OC,0.82,1.20)==TRUE)){
		data$chem[i] = "Carbohydrates"
		data$y[i] = 10
	  } 
	  
	  #with all data
	  ## type DART, ASAP, ESI Acétone, ESI Méthanol 
	  if(grepl("ASAP", data$Ionisation[i])==TRUE){
		data$x[i] = 1
	  }
	  if((grepl("ESI",data$Ionisation[i])==TRUE)&(grepl("Acétone", data$Condition[i])==TRUE)){
		data$x[i] = 2
	  }
	  if((grepl("ESI",data$Ionisation[i])==TRUE)&(grepl("Méthanol", data$Condition[i])==TRUE)){
		data$x[i] = 2
	  }
	  ## type DataBase 
	  if(grepl("DataBase", data$Nature[i])==TRUE){
		data$x[i] = 3
	  }
	}

# sort

######################################################################################################################
# the data.frame contains everything
# Apply desired filters to your data before plotting, if no 'chem' filtering is desired, please skip the 'computing'
#
#########################################################################################################
#tricks
# use str(data) for checking the whole content
	str(data)
# use table(data$x) for checking if the right ionizations are available
	table(data$x)
# Remove rows without a chem tag
	data2 <- subset(data, y!="NA")
# choose the right filters

# in this work the need of ChemQuant and its plotting

#  publication 1 figure 4 a,b,c & d, SI figure S2, S5 a,b,c & d
#  publication 2 figure 3 a, figure 3 b, SI figure S1 a,b,c & d

#  array of focalisation

# with plot for Publication 1

#  publication 1 / [publication 2 supplementary info]
#	4a / S1a : Nature=="Multi 1" & Polarité =="POS"
	data3 <- subset(data2, ((Polarité=="POS")&(Nature=="Multi 1"))==TRUE|((x==4)==TRUE))   
#	4b / S1b : Nature=="Multi 1" & Polarité =="NEG"
	data3 <- subset(data2, ((Polarité=="NEG")&(Nature=="Multi 1"))==TRUE|((x==4)==TRUE))
#	4c / S1c : Nature=="Multi 2" & Polarité =="POS"
	data3 <- subset(data2, ((Polarité=="POS")&(Nature=="Multi 2"))==TRUE|((x==4)==TRUE))
#	4d / S1d : Nature=="Multi 2" & Polarité =="NEG"
	data3 <- subset(data2, ((Polarité=="NEG")&(Nature=="Multi 2"))==TRUE|((x==4)==TRUE))
#  publication 1 supplementary info
#	S2  : Lichen=="Rf" & Nature=="Multi 1" & Polarité =="POS"
	data3 <- subset(data2, (((Polarité=="POS")&(Nature=="Multi 1"))&(Lichen=="Rf"))==TRUE|((x==4)==TRUE))
#	S5a : Nature=="Acétone" & Polarité =="POS"	
	data3 <- subset(data2, ((Polarité=="POS")&(Nature=="Acétone"))==TRUE|((x==4)==TRUE))
#	S5b : Nature=="Acétone" & Polarité =="NEG" 
	data3 <- subset(data2, ((Polarité=="NEG")&(Nature=="Acétone"))==TRUE|((x==4)==TRUE))
#	S5c : Nature=="Méthanol 1" & Polarité =="POS"	
	data3 <- subset(data2, ((Polarité=="POS")&(Nature=="Méthanol 1"))==TRUE|((x==4)==TRUE))
#	S5d : Nature=="Méthanol 1" & Polarité =="NEG"
	data3 <- subset(data2, ((Polarité=="NEG")&(Nature=="Méthanol 1"))==TRUE|((x==4)==TRUE))
	
# 27/04/22 simplification esi methanol + acetone = ESI

#  publication 1 / [publication 2 supplementary info]
#	4a / S1a : Nature=="Multi 1" & Polarité =="POS"
	data3 <- subset(data2, ((Polarité=="POS")&(Nature=="Multi 1"))==TRUE|((x==3)==TRUE))   
#	4b / S1b : Nature=="Multi 1" & Polarité =="NEG"
	data3 <- subset(data2, ((Polarité=="NEG")&(Nature=="Multi 1"))==TRUE|((x==3)==TRUE))
#	4c / S1c : Nature=="Multi 2" & Polarité =="POS"
	data3 <- subset(data2, ((Polarité=="POS")&(Nature=="Multi 2"))==TRUE|((x==3)==TRUE))
#	4d / S1d : Nature=="Multi 2" & Polarité =="NEG"
	data3 <- subset(data2, ((Polarité=="NEG")&(Nature=="Multi 2"))==TRUE|((x==3)==TRUE))
#  publication 1 supplementary info
#	S2  : Lichen=="Rf" & Nature=="Multi 1" & Polarité =="POS"
	data3 <- subset(data2, (((Polarité=="POS")&(Nature=="Multi 1"))&(Lichen=="Rf"))==TRUE|((x==3)==TRUE))
#	S5a : Nature=="Acétone" & Polarité =="POS"	
	data3 <- subset(data2, ((Polarité=="POS")&(Nature=="Acétone"))==TRUE|((x==3)==TRUE))
#	S5b : Nature=="Acétone" & Polarité =="NEG" 
	data3 <- subset(data2, ((Polarité=="NEG")&(Nature=="Acétone"))==TRUE|((x==3)==TRUE))
#	S5c : Nature=="Méthanol 1" & Polarité =="POS"	
	data3 <- subset(data2, ((Polarité=="POS")&(Nature=="Méthanol 1"))==TRUE|((x==3)==TRUE))
#	S5d : Nature=="Méthanol 1" & Polarité =="NEG"
	data3 <- subset(data2, ((Polarité=="NEG")&(Nature=="Méthanol 1"))==TRUE|((x==3)==TRUE))




######################################################################################################################
#
# 										Creating the plot of histograms
#
#									   			Publication 1
#
###############################################################################################
#
# adapted from http://www.sthda.com/english/wiki/impressive-package-for-3d-and-4d-graph-r-software-and-data-visualization
#
######################################################################################################################

	hist3D_fancy<- function(x, y, dataset, break.func = c("Sturges", "scott", "FD"), breaks = NULL,
							colvar = NULL, col="white", clab=NULL, phi = 30, theta = 45, ...){
	  
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
	  #x.bin <- seq(min(x)-1, max(x)+0.5, length.out = 14) # modification NLY 12.04.2022 match with x => 1,2, & 4
	  x.bin <- c(0, 0.35, 0.7, 1, 1.35, 1.7, 2, 2.35, 2.7, 3, 3.35, 3.7, 4, 4.25)
	  y.bin <- seq(0, max(y), length.out = y.breaks)
	  xy <- table(cut(x, x.bin), cut(y, y.bin))
	  z <- prop.table(xy,1)
	  xmid <- c(0.2, 0.5, 1, 1.2, 1.5, 2, 2.2, 2.5, 3, 3.2, 3.5, 4, 4.1)
	  #xmid <- 0.5*(x.bin[-1] + x.bin[-length(x.bin)])
	  
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
	  max.x <- max(x)+0.2
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
	  
	  labels <- c( 'ASAP-MS', 'ESI-MS \n in acetone', 'ESI-MS \n in methanol', 'Database')
	  label.pos <- trans3d(x.axis-0.5, (min.y - 1.85), min.z, H3D)
	  text(label.pos$x, label.pos$y, labels=labels[min(x):max(x)], adj=c(0, NA), cex=0.75, las=1)
	  
	  ##Y
	  labels <- c('Condensed aromatic compounds', 'Polyphenols and derivatives', 'Benzenoids', 'Terpenes',
				  'Unsaturated hydrocarbons', 'Fatty acyls', 'Prenol derivatives', 'Nucleic acids', 'Amino-acids', 'Carbohydrates')
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
	  #x.bin <- c(0.75,1,1.75,2,2.75,3,3.75,4,4.4,4.5)	  
	   x.bin <- c(0, 0.35, 0.7, 1, 1.35, 1.7, 2, 2.35, 2.7, 3, 3.35, 3.7, 4, 4.5)
	  x.bin2 <- x.bin
	  y.bin2 <- seq(0, 10.5, by=0.5)# modification NLY 12.04.2022
	  #y.bin2 <- seq(0, max(dataset$y), length.out = 1.75*y.breaks2)
	  Dx <- table(cut(dataset$x,x.bin2),cut(dataset$y,y.bin2))
	  Dx <- 4*Dx # modification NLY 12.04.2022 (trouble of intensity in colkey)
	  xmid2 <- c(0.2, 0.5, 1, 1.2, 1.5, 2, 2.2, 2.5, 3, 3.2, 3.5, 4, 4.1)
	  #xmid2 <- 0.5*(x.bin2[-1] + x.bin2[-length(x.bin2)])
	  ymid2 <- 0.5*(y.bin2[-1] + y.bin2[-length(y.bin2)])
	  xmid2[1]<-0
	  par(mar=c(4,4,4,4))
	image3D(z = 1, x=xmid2, y=ymid2, colvar = Dx, add = TRUE,
			 colkey = list(length = 0.5, width = 0.4, shift = 0,
			 cex.axis = 0.8, cex.clab = 0.8), 
			 alpha =1,
			 clab = c("","Number of ionized molecules","(Absolute value)"), 
			 plot = TRUE)  
	  par(mar = oldmar)
	}



######################################################################################################################
#
# 										Creating the plot of histograms
#
#									   			Publication 1 // 3 Items (27/04/22)
#
###############################################################################################
#
# adapted from http://www.sthda.com/english/wiki/impressive-package-for-3d-and-4d-graph-r-software-and-data-visualization
#
######################################################################################################################

	hist3D_fancy<- function(x, y, dataset, break.func = c("Sturges", "scott", "FD"), breaks = NULL,
							colvar = NULL, col="white", clab=NULL, phi = 30, theta = 45, ...){
	  
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
	   if(max(x)==3){
	  x.bin <- c(0,0.33,0.75,1,1.33,1.75,2,2.33,2.75,3,3.5) # modification NLY 12.04.2022 match with x => 1,2,3
	  }
	  else 	  x.bin <- seq(min(x)-1, max(x)+0.5, length.out = x.breaks)
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
	  
	  labels <- c( 'ASAP-MS','DI-ESI-MS', 'Database')
	  label.pos <- trans3d(x.axis-0.5, (min.y - 1.85), min.z, H3D)
	  text(label.pos$x, label.pos$y, labels=labels[min(x):max(x)], adj=c(0, NA), cex=0.75, las=1)
	  
	  ##Y
	  labels <- c('Condensed aromatic compounds', 'Polyphenols and derivatives', 'Benzenoids', 'Terpenes',
				  'Unsaturated hydrocarbons', 'Fatty acyls', 'Prenol derivatives', 'Nucleic acids', 'Amino-acids', 'Carbohydrates')
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
	  x.bin <- seq(min(x)-1, max(x)+0.5, length.out = x.breaks)	
	  #x.bin <- seq(-0.5, 3.5, by=0.5)
	  x.bin2 <- x.bin
	  y.bin2 <- seq(0, 10.5, by=0.5)# modification NLY 12.04.2022
	  #y.bin2 <- seq(0, max(dataset$y), length.out = 1.75*y.breaks2)
	  Dx <- table(cut(dataset$x,x.bin2),cut(dataset$y,y.bin2))
	  Dx <- 4*Dx # modification NLY 12.04.2022 (trouble of intensity in colkey)
	  	  
	  xmid2 <- 0.5*(x.bin2[-1] + x.bin2[-length(x.bin2)])
	  ymid2 <- 0.5*(y.bin2[-1] + y.bin2[-length(y.bin2)])
	   xmid2[1]<-0
	  par(mar=c(4,4,4,4))
	image3D(z = 1, x=xmid2, y=ymid2, colvar = Dx, add = TRUE,
			 colkey = list(length = 0.5, width = 0.4, shift = 0,
			 cex.axis = 0.8, cex.clab = 0.8), 
			 alpha =1,
			 clab = c("","Number of ionized molecules","(Absolute value)"), 
			 plot = TRUE)  
	  par(mar = oldmar)
	}





######################################################################################################################
# 										applying to data3
######################################################################################################################

	hist3D_fancy(data3$x, data3$y, data3, breaks = 11)

###########################################
# for the stats easy to use
	table(data3$x,data3$y)

######################################################################################################################
# 										debuging lines
######################################################################################################################

x <- data3$x
y <- data3$y
breaks <- 11
dataset <- data3
 x.breaks <- y.breaks <- breaks
	  
	  # Cut x and y variables in bins for counting
	  x.bin <- seq(min(x)-1, max(x)+0.5, length.out = x.breaks) # modification NLY 12.04.2022 match with x => 1,2, & 4
	  y.bin <- seq(0, max(y), length.out = y.breaks)
	  xy <- table(cut(x, x.bin), cut(y, y.bin))
	  z <- prop.table(xy,1)
	  #xmid <- c(0.25,1,1.75,2,2.75,3,3.75,4,4.4)
	  xmid <- 0.5*(x.bin[-1] + x.bin[-length(x.bin)])	  
	  ymid <- 0.5*(y.bin[-1] + y.bin[-length(y.bin)])
	  
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
y.breaks2 <- breaks+1
	  #x.bin <- c(0.75,1,1.75,2,2.75,3,3.75,4,4.4,4.5)	  
	  x.bin <- seq(min(x)-1, max(x)+0.5, length.out = x.breaks)	
	  x.bin2 <- x.bin
	  y.bin2 <- seq(0, 10.5, by=0.5)# modification NLY 12.04.2022
	  #y.bin2 <- seq(0, max(dataset$y), length.out = 1.75*y.breaks2)
	  Dx <- table(cut(dataset$x,x.bin2),cut(dataset$y,y.bin2))
	  Dx <- 4*Dx # modification NLY 12.04.2022 (trouble of intensity in colkey)
	  
	  xmid2 <- 0.5*(x.bin2[-1] + x.bin2[-length(x.bin2)])
	  ymid2 <- 0.5*(y.bin2[-1] + y.bin2[-length(y.bin2)])
	  xmid2[1]<-0
