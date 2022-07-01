                      #######################################################################
                      #   THE FIRST PART OF THIS SCRIPT IS SIMILAR TO CHEMQUANT,            #
                      #   IF YOU HAVE ALREADY GENERATED THE DATASET W/ CHEMQUANT            #
                      #   PLEASE LOAD THE .RDATA FILE AND GO DIRECTLY TO THE PLOTTING PART  #
                      #######################################################################

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
library(tcltk)
                      
# Set Directory

#input_dir <- # your path beware with / 
	input_dir <- 'C:/Users/CRMPO/Documents/2022/Simon.rapport.DATA/VK/TSV'
	setwd(input_dir)
	temp <- list.files(input_dir, pattern = "dataset_solvents_comp_publi1.tsv")
#temp <- list.files(input_dir, pattern = "dataset_pre-extract_comp_publi2.tsv")
#reading function for tsv, csv
#take care of special characters like accent (we are French)
	data  <- read.csv2(file = temp, header = TRUE, sep ="\t")
#depending on file charateristics
#data <- read.csv2(file = temp, sep = '\t', header = TRUE)
#data <- read.csv2(file = temp, sep = ',', header = TRUE)
#if the file is from OSF, it should be processed and saved with the right (chem, x, y) => goto "sort"

#Intitialisation of chem, y and x
	data$chem <- NA
	data$y <- NA
	data$x <- NA
# Compute chemical family according to Van Krevelen coordinates for all data
# one run for all 3
	for(i in (1:nrow(data))){
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
	  # publication 1 VK diagramm only with ASAP & ESI
	  ## type DART, ASAP, ESI Acétone, ESI Méthanol 
	  #"if(grepl("DART", data$Ionisation[i])==TRUE){
	   #	data$x[i] = 1
	  #}"
	  if(grepl("ASAP", data$Ionisation[i])==TRUE){
		data$x[i] = 2
	  }
	  if((grepl("ESI",data$Ionisation[i])==TRUE)&(grepl("Acétone", data$Condition[i])==TRUE)){
		data$x[i] = 1
	  }
	  #27/04/22 - Simplification ESI is Acétone and Méthanol
	  if((grepl("ESI",data$Ionisation[i])==TRUE)&(grepl("Méthanol", data$Condition[i])==TRUE)){
		data$x[i] = 3
	  }
	}

# sort

######################################################################################################################
# the data.frame contains everything
# Apply desired filters to your data before plotting, if no 'chem' filtering is desired, please skip the 'computing'
# part of the script to reduce calculation time beware for plotting frequency data are really big some 100 Mo up to 1 or 2 Go
###########################################################################################################
# tricks
###########################################################################################################
# use str(data) for checking the whole content
#	str(data)
# use table(data$x) for checking if the right ionizations are available
#	table(data$x)
###########################################################################################################

# data of interest have y!="NA"
# Remove rows without a chem tag
	data <- subset(data, y!="NA")
# one should first sort data
	data2 <- data[order(data[,"x"],data[,"O.C"],data[,"H.C"]), ]
# choose the right filters

# in this work the need of VKplot and its plotting

#  publication 1 figure 5 a & b SI figure S4 a & b
#  publication 2 figure 3 c, SI figure S2 => 5b- VKplot.TSV.publi2.R

#  array of focalisation

#  publication 1
#	5a	: Nature=="Multi 1" & y==3
	data3 <- subset(data2, ((y==3)&(Nature=="Multi 1"))==TRUE)
#	5b  : Nature=="Acétone" & y==3
	data3 <- subset(data2, ((y==3)&(Nature=="Acétone"))==TRUE)

#  publication 1 supplementary info with y==NA please
#	S4a : Nature=="Multi 1" & Polarité =="POS"	
	data3 <- subset(data2, ((Polarité=="POS")&(Nature=="Multi 1"))==TRUE)
#	S4b : Nature=="Multi 1" & Polarité =="NEG" 
	data3 <- subset(data2, ((Polarité=="NEG")&(Nature=="Multi 1"))==TRUE)

######################################################################################################################
#
# Creating the plot
###############################################################################################
# the frequence of apparition for a given couple O.C and H.C is needed
# this will correspond to the size of the colored dot or circle
# with the correlation with x for color layer.
# frequence of x, O.C, H.C
	frequence <- as.data.frame(table(data3$x, data3$O.C, data3$H.C))
# frequence is a mix of factor and classical variable
# it's not compatible with following process (scatter2d)
# small trick write and read
	write.csv2(frequence, file ="myfilefrequence.csv")
	data4 <- read.csv2(file ="myfilefrequence.csv", header=TRUE)
# always stay tidy
	unlink("myfilefrequence.csv")
# only data with existing O.C H.C are wished
	data4 <- subset(data4, Freq!=0)
	colnames(data4) <- c('X','x','O.C','H.C','Freq')
######################################################################################################################
# 
# 					addition of a dimension zcolor to qualify any possible color layer.
#					original 34 layers : 7 colors
#					A, B, C, A=B=C, A=B, B=C, A=C.
#
###############################################################################################
	data4$zcolor = 0
	i<-1
	while(i<(nrow(data4)))
	{
		if(((data4$O.C[i]==data4$O.C[i+1])&&(data4$H.C[i]==data4$H.C[i+1]))&&((data4$O.C[i]==data4$O.C[i+2])&&(data4$H.C[i]==data4$H.C[i+2])))
		{
		#   triple points found
		#  	line x   O.C  H.C    Freq
		#	i    1	 0.25 0.85   freq1
		#	i+1  2   0.25 0.85   freq2
		#	i+2  3   0.25 0.85   freq3
			if((data4$Freq[i]==data4$Freq[i+1]))
			{
			# freq1 (x=1) == freq2 (x=2) #
			# what about freq3 (x=3)
				if((data4$Freq[i+1]==data4$Freq[i+2]))
				{
				# freq1 == freq2 == fre3# 
				# same size = other color (not x=1 or x=2 or x=3)
				# the color is given to the first
				data4$zcolor[i] = 13
				i<-i+3
				}
				else
				{
				if((data4$Freq[i+1]>data4$Freq[i+2]))
				{
				# (freq1 == freq2) > fre3# 
				# same size = other color (not x=1 or x=2 or x=3) for 1 and 2
				# the color is given to the first
				# 3 stays above with its color

				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				myx1 = data4$x[i]
				myx2 = data4$x[i+1]
				myx3 = data4$x[i+2]
				# 3 1 2

				data4$Freq[i] 	= myFreq1 
				data4$Freq[i+1] = myFreq2  
				data4$Freq[i+2] = myFreq3
				data4$x[i]   = myx1
				data4$x[i+1] = myx2
				data4$x[i+2] = myx3				
				
				data4$zcolor[i] = 14
				data4$zcolor[i+2] = 17
				i<-i+3
				}
				else
				if((data4$Freq[i+1]<data4$Freq[i+2]))
				{
				#  fre3 > (freq1 == freq2) # 
				# same size = other color (not x=1 or x=2 or x=3) for 1 and 2
				# the color is given to the first
				# 3 goes under with its color (exchange with 1)

				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				myx1 = data4$x[i]
				myx2 = data4$x[i+1]
				myx3 = data4$x[i+2]
				# 3 1 2

				data4$Freq[i] 	= myFreq3 
				data4$Freq[i+1] = myFreq1  
				data4$Freq[i+2] = myFreq2
				data4$x[i]   = myx3
				data4$x[i+1] = myx1
				data4$x[i+2] = myx2

				data4$zcolor[i] = 20
				data4$zcolor[i+1] = 23
				i<-i+3
				}
				}
			}
			else
			{
			if((data4$Freq[i]==data4$Freq[i+2]))
			{
			# freq1 (x=1) == freq3 (x=3) #
			# what about freq2 (x=2)
				if((data4$Freq[i]>data4$Freq[i+1]))
				{
				# (freq1 == freq3) > fre2# 
				# same size = other color (not x=1 or x=2 or x=3) for 1 and 3
				# the color is given to the first
				# 2 stays above with its color

				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				myx1 = data4$x[i]
				myx2 = data4$x[i+1]
				myx3 = data4$x[i+2]
				# (2 3) 1

				data4$Freq[i] 	= myFreq1 
				data4$Freq[i+1] = myFreq3  
				data4$Freq[i+2] = myFreq2
				data4$x[i]   = myx1
				data4$x[i+1] = myx3
				data4$x[i+2] = myx2
				
				data4$zcolor[i] = 15	
				data4$zcolor[i+2] = 18		
				i<-i+3
				}
				else
				if((data4$Freq[i]<data4$Freq[i+1]))
				{
				#  fre2 > (freq1 == freq3) # 
				# same size = other color (not x=1 or x=2 or x=3) for 1 and 3
				# the color is given to the first
				# 2 goes under with its color (exchange with 1)

				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				myx1 = data4$x[i]
				myx2 = data4$x[i+1]
				myx3 = data4$x[i+2]
				# (2 3) 1

				data4$Freq[i] 	= myFreq2 
				data4$Freq[i+1] = myFreq1  
				data4$Freq[i+2] = myFreq3
				data4$x[i]   = myx2
				data4$x[i+1] = myx1
				data4$x[i+2] = myx3
				
				data4$zcolor[i] = 21				
				data4$zcolor[i+1] = 24
						
				i<-i+3
				}
			}
			else 
			if((data4$Freq[i+1]==data4$Freq[i+2]))
			{
			# freq2 (x=2) == freq3 (x=3) #
			# what about freq1 (x=1)
				if((data4$Freq[i]>data4$Freq[i+1]))
				{
				# freq1 > (freq3= fre2)# 
				# same size = other color (not x=1 or x=2 or x=3) for 2 and 3
				# the color is given to the first
				# 1 stays underneath with its color

				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				myx1 = data4$x[i]
				myx2 = data4$x[i+1]
				myx3 = data4$x[i+2]
				# (2 3) 1

				data4$Freq[i] 	= myFreq1 
				data4$Freq[i+1] = myFreq3  
				data4$Freq[i+2] = myFreq2
				data4$x[i]   = myx1
				data4$x[i+1] = myx3
				data4$x[i+2] = myx2
				data4$zcolor[i] = 22
				data4$zcolor[i+1] = 25				
				i<-i+3
				}
				else
				if((data4$Freq[i]<data4$Freq[i+1]))
				{
				#  (freq2 == freq3) > freq1 # 
				# same size = other color (not x=1 or x=2 or x=3) for 2 and 3
				# the color is given to the first
				# 1 goes 3, 2 goes 1 3 goes 2

				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				myx1 = data4$x[i]
				myx2 = data4$x[i+1]
				myx3 = data4$x[i+2]
				# (2 3) 1

				data4$Freq[i] 	= myFreq2 
				data4$Freq[i+1] = myFreq3  
				data4$Freq[i+2] = myFreq1
				data4$x[i]   = myx2
				data4$x[i+1] = myx3
				data4$x[i+2] = myx1				
				data4$zcolor[i] = 16
				data4$zcolor[i+2] = 19
				i<-i+3
				}
			}
			else 
			if((data4$Freq[i]>data4$Freq[i+1]))
			{
			# freq1 (x=1) > freq2 (x=2) #
			# what about freq3 (x=3)	
				if((data4$Freq[i+1]>data4$Freq[i+2]))
				{
				# freq1 > freq2 > freq3# 
				# 3 stay in position each with its color
				data4$zcolor[i]   = 26
				data4$zcolor[i+1] = 30
				data4$zcolor[i+2] = 34
				i<-i+3
				}
				else
				if((data4$Freq[i]>data4$Freq[i+2]))
					{
					# freq1 > freq3 > fre2# 
					# 1 stays 
					# 3 goes under 2 with its color (swap 2 3)

					myFreq = data4$Freq[i+1]
					data4$Freq[i+1] = data4$Freq[i+2]
					data4$Freq[i+2] = myFreq
					myx = data4$x[i+1]
					data4$x[i+1] = data4$x[i+2]
					data4$x[i+2] = myx					
					
					data4$zcolor[i]   = 26
					data4$zcolor[i+1] = 31
					data4$zcolor[i+2] = 33
					i<-i+3
					}
					else
					{
					# freq3 > freq1 > fre2# 
					# 1 stays 
					# 3 goes under 2 with its color

					myFreq1 = data4$Freq[i]
					myFreq2 = data4$Freq[i+1]
					myFreq3 = data4$Freq[i+2]
					myx1 = data4$x[i]
					myx2 = data4$x[i+1]
					myx3 = data4$x[i+2]
					# 3 1 2

					data4$Freq[i] 	= myFreq3 
					data4$Freq[i+1] = myFreq1  
					data4$Freq[i+2] = myFreq2
					data4$x[i]   = myx3
					data4$x[i+1] = myx1
					data4$x[i+2] = myx2
					data4$zcolor[i]   = 28
					data4$zcolor[i+1] = 29
					data4$zcolor[i+2] = 33	



					i<-i+3					
					}			
			}
			else
			if((data4$Freq[i+1]>data4$Freq[i]))
			{
			# freq2 (x=2) > freq1 (x=1) #
			# what about freq3 (x=3)	
				if((data4$Freq[i]>data4$Freq[i+2]))
				{
				# freq2 > freq1 > fre3# 
				# 3 stays in position
				# 2 goes under 1

					myFreq1 = data4$Freq[i]
					myFreq2 = data4$Freq[i+1]
					myFreq3 = data4$Freq[i+2]
					myx1 = data4$x[i]
					myx2 = data4$x[i+1]
					myx3 = data4$x[i+2]
					# 2 1 3

					data4$Freq[i] 	= myFreq2 
					data4$Freq[i+1] = myFreq1  
					data4$Freq[i+2] = myFreq3
					data4$x[i]   = myx2
					data4$x[i+1] = myx1
					data4$x[i+2] = myx3
					data4$zcolor[i]   = 27
					data4$zcolor[i+1] = 29
					data4$zcolor[i+2] = 34	

				i<-i+3
				}
				else{
					if((data4$Freq[i+1]>data4$Freq[i+2]))
					{
					# freq2 > freq3 > fre1# 
					# 1 goes above 
					# 2 and 3 under 1

					myFreq1 = data4$Freq[i]
					myFreq2 = data4$Freq[i+1]
					myFreq3 = data4$Freq[i+2]
					myx1 = data4$x[i]
					myx2 = data4$x[i+1]
					myx3 = data4$x[i+2]
					# 2 3 1 

					data4$Freq[i] 	= myFreq2 
					data4$Freq[i+1] = myFreq3  
					data4$Freq[i+2] = myFreq1
					data4$x[i]   = myx2
					data4$x[i+1] = myx3
					data4$x[i+2] = myx1
					data4$zcolor[i]   = 27
					data4$zcolor[i+1] = 31
					data4$zcolor[i+2] = 32	

					i<-i+3
					}	
					else
					{
					# freq3 > freq2 > fre1# 
					# 1 goes above
					# 3 goes under 2 with its color

					myFreq1 = data4$Freq[i]
					myFreq2 = data4$Freq[i+1]
					myFreq3 = data4$Freq[i+2]
					myx1 = data4$x[i]
					myx2 = data4$x[i+1]
					myx3 = data4$x[i+2]
					# 3 1 2

					data4$Freq[i] 	= myFreq3 
					data4$Freq[i+1] = myFreq2  
					data4$Freq[i+2] = myFreq1
					data4$x[i]   = myx3
					data4$x[i+1] = myx2
					data4$x[i+2] = myx1
					data4$zcolor[i]   = 28
					data4$zcolor[i+1] = 30
					data4$zcolor[i+2] = 32	

					i<-i+3					
					}
					}			
			}
			}			
		}
		else
		{
			if((data4$O.C[i]==data4$O.C[i+1])&&(data4$H.C[i]==data4$H.C[i+1]))
			{
			#   double point found
			#  	line x   O.C  H.C    Freq
			#	i    1	 0.25 0.85   freq1
			#	i+1  2   0.25 0.85   freq2
			#	i+2  1   0.26 0.85   freq3
			#	at the end the cursor is driven to i+2	
				if((data4$Freq[i]==data4$Freq[i+1]))
				{
				# freq1 == freq2 # 
				# same size = other color (not x=1 or x=2 or x=3)
				# the color is given to the first
				# 3 cases possible
					if((data4$x[i]==1)&&(data4$x[i+1]==2))
					{
					# x=1, x=2	
						data4$zcolor[i] = 4
						i<-i+2
					}
					else
					{
					if((data4$x[i]==1)&&(data4$x[i+1]==3))
					{
					# x=1, x=3	
						data4$zcolor[i] = 5
						i<-i+2
					}
					else
					if((data4$x[i]==2)&&(data4$x[i+1]==3))
					{
					# x=2, x=3
						data4$zcolor[i] = 6
						i<-i+2
					}
					}

				
				}
				else
				{
					if((data4$Freq[i]>data4$Freq[i+1]))
					{
						if((data4$x[i]==1)&&(data4$x[i+1]==2))
						{
						# freq1 (x=1) > freq2 (x=2) #
						data4$zcolor[i] = 7  # will be the color of x=1 (big and underneath)
						#data4$zcolorcode[i]= "#0000FF"
						data4$zcolor[i+1] = 8 # will be the color of x=2 (small and above)
						#data4$zcolorcode[i+1]= "#FF00FF"
						i<-i+2
						}
						else
						{
							if((data4$x[i]==1)&&(data4$x[i+1]==3))
							{
							# freq1 (x=1) > freq3 (x=3) #
							data4$zcolor[i] = 7  # will be the color of x=1 (big and underneath)
							#data4$zcolorcode[i]= "#0000FF"
							data4$zcolor[i+1] = 9 # will be the color of x=3 (small and above)
							#data4$zcolorcode[i+1]= "#FF00FF"
							i<-i+2
							}
							else
							if((data4$x[i]==2)&&(data4$x[i+1]==3))
							{
							# freq2 (x=2) > freq3 (x=3) #
							data4$zcolor[i] = 8  # will be the color of x=2 (big and underneath)
							#data4$zcolorcode[i]= "#0000FF"
							data4$zcolor[i+1] = 9 # will be the color of x=3 (small and above)
							#data4$zcolorcode[i+1]= "#FF00FF"
							i<-i+2
						}
						}

					}
					else
					{
						if((data4$x[i]==1)&&(data4$x[i+1]==2))
						{
						#myX = data4$X[i+1]
						#data4$X[i+1] = data4$X[i]
						#data4$X[i] = myX
						myFreq = data4$Freq[i+1]
						data4$Freq[i+1] = data4$Freq[i]
						data4$Freq[i] = myFreq
						myx = data4$x[i+1]
						data4$x[i+1] = data4$x[i]
						data4$x[i] = myx
						# freq1 (x=1) < freq2 (x=2) # 2 1
						data4$zcolor[i] = 11  # will be the color of x=2 (big and underneath)
						data4$zcolor[i+1] = 12 # will be the color of x=1 (small and above)

						i<-i+2
						}
						else
						{
						if((data4$x[i]==1)&&(data4$x[i+1]==3))
						{
						#myX = data4$X[i+1]
						#data4$X[i+1] = data4$X[i]
						#data4$X[i] = myX
						myFreq = data4$Freq[i+1]
						data4$Freq[i+1] = data4$Freq[i]
						data4$Freq[i] = myFreq
						myx = data4$x[i+1]
						data4$x[i+1] = data4$x[i]
						data4$x[i] = myx
						# freq1 (x=1) < freq3 (x=3) # 3 1
						data4$zcolor[i] = 10  # will be the color of x=3 (big and underneath)
						#data4$zcolorcode[i]= "#0000FF"
						data4$zcolor[i+1] = 12 # will be the color of x=1 (small and above)
						#data4$zcolorcode[i+1]= "#FF00FF"
						i<-i+2
						}
						else
						if((data4$x[i]==2)&&(data4$x[i+1]==3))
						{
						#myX = data4$X[i+1]
						#data4$X[i+1] = data4$X[i]
						#data4$X[i] = myX
						myFreq = data4$Freq[i+1]
						data4$Freq[i+1] = data4$Freq[i]
						data4$Freq[i] = myFreq
						myx = data4$x[i+1]
						data4$x[i+1] = data4$x[i]
						data4$x[i] = myx
						# freq2 (x=2) < freq3 (x=3) # 3 2
						data4$zcolor[i] = 10  # will be the color of x=3 (big and underneath)
						#data4$zcolorcode[i]= "#0000FF"
						data4$zcolor[i+1] = 11 # will be the color of x=2 (small and above)
						#data4$zcolorcode[i+1]= "#FF00FF"
						i<-i+2
						}
						}						
				
					}			
				}
			}	
			else
			{			
		#   single point found
		#  	line x      O.C H.C    Freq
		#	i    xi	   0.25 0.85   freq1
		#	i+1  xi+1  0.29 0.85   freq2
		#	i+2  xi+2  0.29 0.85   freq3
		#	at the end the cursor is driven to i+1
				if((data4$x[i]==1))
				{
				# x=1 #
				data4$zcolor[i] = 1  
				}
				if((data4$x[i]==2))
				{
				# x=2 #
				data4$zcolor[i] = 2  
							}
				if((data4$x[i]==3))
				{
				# x=3 #
				data4$zcolor[i] = 3  
				}
				i<-i+1
			}
		}
		if (i==nrow(data4))
		{
		# last line before 
			if((data4$x[i]==1))
			{
			# x=1 #
			data4$zcolor[i] = 1  
			}
			if((data4$x[i]==2))
			{
			# x=2 #
			data4$zcolor[i] = 2  
						}
			if((data4$x[i]==3))
			{
			# x=3 #
			data4$zcolor[i] = 3  
			}
			i<-i+1
		}
	}

######################################################################################################################
# 
# 				addition of a dimension zcolor to qualify any possible color layer.
#				3 levels of x, 34 levels but zcolor only 4 colors
#				A only, B only, C only  D= the rest 
#				the size of the first will be the sum of the rest
#
###############################################################################################
	
	data4$zcolor = 0
	i<-1
	while(i<(nrow(data4)))
	{
		if(((data4$O.C[i]==data4$O.C[i+1])&&(data4$H.C[i]==data4$H.C[i+1]))&&((data4$O.C[i]==data4$O.C[i+2])&&(data4$H.C[i]==data4$H.C[i+2])))
		{
		#   triple points found
		#  	line x   O.C  H.C    Freq
		#	i    1	 0.25 0.85   freq1
		#	i+1  2   0.25 0.85   freq2
		#	i+2  3   0.25 0.85   freq3
			if((data4$Freq[i]==data4$Freq[i+1]))
			{
			# freq1 (x=1) == freq2 (x=2) #
			# what about freq3 (x=3)
				if((data4$Freq[i+1]==data4$Freq[i+2]))
				{
				# freq1 == freq2 == fre3# 
				# same size = other color (not x=1 or x=2 or x=3)
				# the color is given to the first
				data4$zcolor[i] = 4
				data4$Freq[i] = 3*data4$Freq[i]
				i<-i+3
				}
				else
				{
				if((data4$Freq[i+1]>data4$Freq[i+2]))
				{
				# (freq1 == freq2) > fre3# 
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				data4$Freq[i] 	= myFreq1+ myFreq2  +myFreq3					
				data4$zcolor[i] = 4

				i<-i+3
				}
				else
				if((data4$Freq[i+1]<data4$Freq[i+2]))
				{
				#  fre3 > (freq1 == freq2) # 
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				data4$Freq[i] 	= myFreq1+ myFreq2  +myFreq3					
				data4$zcolor[i] = 4
				i<-i+3
				}
				}
			}
			else
			{
			if((data4$Freq[i]==data4$Freq[i+2]))
			{
			# freq1 (x=1) == freq3 (x=3) #
			# what about freq2 (x=2)
				if((data4$Freq[i]>data4$Freq[i+1]))
				{
				# (freq1 == freq3) > fre2# 
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				data4$Freq[i] 	= myFreq1+ myFreq2  +myFreq3					
				data4$zcolor[i] = 4	
				i<-i+3
				}
				else
				if((data4$Freq[i]<data4$Freq[i+1]))
				{
				#  fre2 > (freq1 == freq3) # 
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				data4$Freq[i] 	= myFreq1+ myFreq2  +myFreq3					
				data4$zcolor[i] = 4
						
				i<-i+3
				}
			}
			else 
			if((data4$Freq[i+1]==data4$Freq[i+2]))
			{
			# freq2 (x=2) == freq3 (x=3) #
			# what about freq1 (x=1)
				if((data4$Freq[i]>data4$Freq[i+1]))
				{
				# freq1 > (freq3= fre2)# 
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				data4$Freq[i] 	= myFreq1+ myFreq2  +myFreq3					
				data4$zcolor[i] = 4			
				i<-i+3
				}
				else
				if((data4$Freq[i]<data4$Freq[i+1]))
				{
				#  (freq2 == freq3) > freq1 # 
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				data4$Freq[i] 	= myFreq1+ myFreq2  +myFreq3					
				data4$zcolor[i] = 4
				i<-i+3
				}
			}
			else 
			if((data4$Freq[i]>data4$Freq[i+1]))
			{
			# freq1 (x=1) > freq2 (x=2) #
			# what about freq3 (x=3)	
				if((data4$Freq[i+1]>data4$Freq[i+2]))
				{
				# freq1 > freq2 > freq3# 
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				data4$Freq[i] 	= myFreq1+ myFreq2  +myFreq3					
				data4$zcolor[i] = 4
				i<-i+3
				}
				else
				if((data4$Freq[i]>data4$Freq[i+2]))
					{
					# freq1 > freq3 > fre2# 
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				data4$Freq[i] 	= myFreq1+ myFreq2  +myFreq3					
				data4$zcolor[i] = 4
					i<-i+3
					}
					else
					{
					# freq3 > freq1 > fre2# 
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				data4$Freq[i] 	= myFreq1+ myFreq2  +myFreq3					
				data4$zcolor[i] = 4



					i<-i+3					
					}			
			}
			else
			if((data4$Freq[i+1]>data4$Freq[i]))
			{
			# freq2 (x=2) > freq1 (x=1) #
			# what about freq3 (x=3)	
				if((data4$Freq[i]>data4$Freq[i+2]))
				{
				# freq2 > freq1 > fre3# 
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				data4$Freq[i] 	= myFreq1+ myFreq2  +myFreq3					
				data4$zcolor[i] = 4

				i<-i+3
				}
				else{
					if((data4$Freq[i+1]>data4$Freq[i+2]))
					{
					# freq2 > freq3 > fre1# 
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				data4$Freq[i] 	= myFreq1+ myFreq2  +myFreq3					
				data4$zcolor[i] = 4	

					i<-i+3
					}	
					else
					{
					# freq3 > freq2 > fre1# 
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				myFreq3 = data4$Freq[i+2]
				data4$Freq[i] 	= myFreq1+ myFreq2  +myFreq3					
				data4$zcolor[i] = 4

					i<-i+3					
					}
					}			
			}
			}			
		}
		else
		{
			if((data4$O.C[i]==data4$O.C[i+1])&&(data4$H.C[i]==data4$H.C[i+1]))
			{
			#   double point found
			#  	line x   O.C  H.C    Freq
			#	i    1	 0.25 0.85   freq1
			#	i+1  2   0.25 0.85   freq2
			#	i+2  1   0.26 0.85   freq3
			#	at the end the cursor is driven to i+2	
				if((data4$Freq[i]==data4$Freq[i+1]))
				{
				# freq1 == freq2 # 
				# same size = other color (not x=1 or x=2 or x=3)
				# the color is given to the first
				# 3 cases possible
					if((data4$x[i]==1)&&(data4$x[i+1]==2))
					{
					# x=1, x=2	
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				data4$Freq[i] 	= myFreq1+ myFreq2					
				data4$zcolor[i] = 4
						i<-i+2
					}
					else
					{
					if((data4$x[i]==1)&&(data4$x[i+1]==3))
					{
					# x=1, x=3	
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				data4$Freq[i] 	= myFreq1+ myFreq2					
				data4$zcolor[i] = 4
						i<-i+2
					}
					else
					if((data4$x[i]==2)&&(data4$x[i+1]==3))
					{
					# x=2, x=3
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				data4$Freq[i] 	= myFreq1+ myFreq2					
				data4$zcolor[i] = 4
						i<-i+2
					}
					}

				
				}
				else
				{
					if((data4$Freq[i]>data4$Freq[i+1]))
					{
						if((data4$x[i]==1)&&(data4$x[i+1]==2))
						{
						# freq1 (x=1) > freq2 (x=2) #
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				data4$Freq[i] 	= myFreq1+ myFreq2					
				data4$zcolor[i] = 4
						i<-i+2
						}
						else
						{
							if((data4$x[i]==1)&&(data4$x[i+1]==3))
							{
							# freq1 (x=1) > freq3 (x=3) #
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				data4$Freq[i] 	= myFreq1+ myFreq2					
				data4$zcolor[i] = 4
							i<-i+2
							}
							else
							if((data4$x[i]==2)&&(data4$x[i+1]==3))
							{
							# freq2 (x=2) > freq3 (x=3) #
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				data4$Freq[i] 	= myFreq1+ myFreq2					
				data4$zcolor[i] = 4
							i<-i+2
						}
						}

					}
					else
					{
						if((data4$x[i]==1)&&(data4$x[i+1]==2))
						{
						# freq1 (x=1) < freq2 (x=2) # 2 1
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				data4$Freq[i] 	= myFreq1+ myFreq2					
				data4$zcolor[i] = 4

						i<-i+2
						}
						else
						{
						if((data4$x[i]==1)&&(data4$x[i+1]==3))
						{
						# freq1 (x=1) < freq3 (x=3) # 3 1
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				data4$Freq[i] 	= myFreq1+ myFreq2					
				data4$zcolor[i] = 4
						i<-i+2
						}
						else
						if((data4$x[i]==2)&&(data4$x[i+1]==3))
						{
						# freq2 (x=2) < freq3 (x=3) # 3 2
				myFreq1 = data4$Freq[i]
				myFreq2 = data4$Freq[i+1]
				data4$Freq[i] 	= myFreq1+ myFreq2					
				data4$zcolor[i] = 4
						i<-i+2
						}
						}						
				
					}			
				}
			}	
			else
			{			
		#   single point found
		#  	line x      O.C H.C    Freq
		#	i    xi	   0.25 0.85   freq1
		#	i+1  xi+1  0.29 0.85   freq2
		#	i+2  xi+2  0.29 0.85   freq3
		#	at the end the cursor is driven to i+1
				if((data4$x[i]==1))
				{
				# freq1 (x=1) #
				data4$zcolor[i] = 1  
				}
				if((data4$x[i]==2))
				{
				# fre2 x=2 #
				data4$zcolor[i] = 2  
							}
				if((data4$x[i]==3))
				{
				#  freq3 (x=3) #
				data4$zcolor[i] = 3  
				}
				i<-i+1
			}
		}
		if (i==nrow(data4))
		{
		# last line before 
			if((data4$x[i]==1))
			{
			# x=1 #
			data4$zcolor[i] = 1  
			}
			if((data4$x[i]==2))
			{
			# x=2 #
			data4$zcolor[i] = 2  
						}
			if((data4$x[i]==3))
			{
			# x=3 #
			data4$zcolor[i] = 3  
			}
			i<-i+1
		}
	}




##############################################
#
#           plotting function
#              from plot3D
#				expe data color=f(zcolor)
#		
# https://rdrr.io/cran/plot3D/man/scatter.html
#############################################

	plotVanKrevelen.expe <- function(myfrequence, alpha , radius, mini, nbcol){
	c2 <- "#0000FF"
	c3 <- "#ed2b2b"
	#c2 <- "#00FFFF"
	c1 <- "#339900"	
	c4 <- "#9900CC"
	c5 <- "#00FFFF"
	c6 <- "#00FFFF"
	c7 <- "#0000FF"
	if (nbcol=="34")
	{
	mycol = c(c1,c2,c3,c4,c5,c6,c1,c2,c3,c3,c2,c1,c7,c4,c5,c6,c3,c2,c1,c3,c2,c1,c4,c5,c6,c1,c2,c3,c1,c2,c3,c1,c2,c3)
	}
	if (nbcol=="3")
	{
	mycol = c(c1,c2,c3)
	}
	if (nbcol=="4")
	{
	mycol = c(c1,c2,c3,c4)
	}
	if (nbcol=="1")
	{
	mycol = c(c4)
	}
	if (nbcol=="6")
	{
	mycol = c(c1,c2,c3,c4,c5,c6)
	}
	if (nbcol=="12")
	{
	mycol = c(c4,c5,c6,c3,c2,c1,c3,c2,c1,c4,c5,c6)
	}
	  # Correlate cex to frequency of Van Krevelen coordinates
	  # ici mafrequence est un tableau global avec O.C, H.C, x et la fréquence en variable
	  # ainsi x va fournir la couleur, Freq le cex nécessaire
	  # le cex est 50 % de la frequence : c'est améliorable 
	  # il serait intéressant de définir une granulométrie en fonction du % de représentatitivité

	scatter2D (myfrequence$O.C, myfrequence$H.C, colvar = myfrequence$zcolor,
col = ramp.col (col = mycol, n = 100, alpha = 1),#level existing complete overview
#				col = ramp.col (col = c("#0000FF","#FF00FF","#339900","#9900CC","#9900CC","#00FFFF","#CCCCFF","#CC9966","#000099","#CC33CC","#99FF33","#00FFFF","#CCCCFF","#CC9966","#000099","#CC33CC","#99FF33"), n = 100, alpha = 1),#level existing complete overview
#				col = ramp.col (col = c("#0000FF","#FF00FF","#FF00FF","#0000FF"), n = 100, alpha = 1),#level existing in both with <> frequence
#				col = ramp.col (col = c("#0000FF","#FF00FF"), n = 100, alpha = 1),#level existing in both with <> frequence
#				col = ramp.col (col = c("#0000FF","#FF00FF","#339900"), n = 100, alpha = 1),#level existing in both with <> frequence
#col = ramp.col (col = c("#0000FF","#FF00FF","#339900","#9900CC","#00FFFF","#CCCCFF","#CC9966","#000099","#CC33CC","#99FF33"), n = 100, alpha = 1),#level existing complete overview
			   NAcol = "white", breaks = NULL,
			   colkey = NULL, clim = NULL, clab = NULL, 
			   alpha=alpha,
			   pch = 20, cex = radius*(myfrequence$Freq)+mini,
			   xlab="O/C ratio", ylab="H/C ratio",
			   CI = NULL, add = FALSE, plot = TRUE)
	}



data5 <- subset(data4, zcolor>=1&zcolor<=3)
data5 <- subset(data4, zcolor>=1&zcolor<=4)
plotVanKrevelen.expe (data5, 0.4,1,0,1)


data5 <- subset(data4, zcolor>=1&zcolor<=3)
plotVanKrevelen.expe (data5, 0.4,0.5,0,3)






##############################################
#
#           plotting function
#              from plot3D
#				expe data color=f(zcolor)
#		
# https://rdrr.io/cran/plot3D/man/scatter.html
#############################################

	plotVanKrevelen.expe.onecolor <- function(myfrequence, alpha , radius, mini){
	c2 <- "#0000FF"
	c3 <- "#ed2b2b"
	#c2 <- "#00FFFF"
	c1 <- "#339900"	
	c4 <- "#9900CC"
	c5 <- "#00FFFF"
	c6 <- "#00FFFF"
	c7 <- "#0000FF"
	
	  # Correlate cex to frequency of Van Krevelen coordinates
	  # ici mafrequence est un tableau global avec O.C, H.C, x et la fréquence en variable
	  # ainsi x va fournir la couleur, Freq le cex nécessaire
	  # le cex est 50 % de la frequence : c'est améliorable 
	  # il serait intéressant de définir une granulométrie en fonction du % de représentatitivité

	scatter2D (myfrequence$O.C, myfrequence$H.C, colvar = myfrequence$zcolor,
col = "#9900CC",#level existing complete overview
#				col = ramp.col (col = c("#0000FF","#FF00FF","#339900","#9900CC","#9900CC","#00FFFF","#CCCCFF","#CC9966","#000099","#CC33CC","#99FF33","#00FFFF","#CCCCFF","#CC9966","#000099","#CC33CC","#99FF33"), n = 100, alpha = 1),#level existing complete overview
#				col = ramp.col (col = c("#0000FF","#FF00FF","#FF00FF","#0000FF"), n = 100, alpha = 1),#level existing in both with <> frequence
#				col = ramp.col (col = c("#0000FF","#FF00FF"), n = 100, alpha = 1),#level existing in both with <> frequence
#				col = ramp.col (col = c("#0000FF","#FF00FF","#339900"), n = 100, alpha = 1),#level existing in both with <> frequence
#col = ramp.col (col = c("#0000FF","#FF00FF","#339900","#9900CC","#00FFFF","#CCCCFF","#CC9966","#000099","#CC33CC","#99FF33"), n = 100, alpha = 1),#level existing complete overview
			   NAcol = "white", breaks = NULL,
			   colkey = NULL, clim = NULL, clab = NULL, 
			   alpha=alpha,
			   pch = 20, cex = radius*(myfrequence$Freq)+mini,
			   xlab="O/C ratio", ylab="H/C ratio",
			   CI = NULL, add = FALSE, plot = TRUE)
	}
data5 <- subset(data4, zcolor==4)
plotVanKrevelen.expe.onecolor (data5, 0.4,0.1,0)




#
##############################################
#
#           plotting function
#              from plot3D
# https://rdrr.io/cran/plot3D/man/scatter.html
#############################################

	plotVanKrevelen.f <- function(myfrequence){
	  
	  # Correlate cex to frequency of Van Krevelen coordinates
	  # ici mafrequence est un tableau global avec O.C, H.C, x et la fréquence en variable
	  # ainsi x va fournir la couleur, Freq le cex nécessaire
	  # le cex est 50 % de la frequence : c'est améliorable 
	  # il serait intéressant de définir une granulométrie en fonction du % de représentatitivité

	scatter2D (myfrequence$O.C, myfrequence$H.C, colvar = myfrequence$x,
			   col = NULL, NAcol = "white", breaks = NULL,
			   colkey = NULL, clim = NULL, clab = NULL, 
			   alpha=0.4,
			   pch = 20, cex = 0.075*(myfrequence$Freq),
			   xlab="O/C ratio", ylab="H/C ratio",
			   CI = NULL, add = FALSE, plot = TRUE)
	}
	plotVanKrevelen.f(data4)

#########################################################
# trash
#input_dir <- 'C:/Users/CRMPO/Documents/2022/Simon.rapport.DATA/VK/TSV'
#input_dir <- 'C:/Users/Simon Ollivier/Desktop/working directory/3- preVKs'
#input_dir <- 'C:/Users/CRMPO/Documents/2022/Simon.rapport.DATA/VK/3-preVKs'
#input_dir <- 'C:/Users/CRMPO/Documents/2022/Simon.rapport.DATA/VK/3-preVKs-2'
#input_dir <- 'C:/Users/CRMPO/Documents/2022/Simon.rapport.DATA/VK/preVKs-3'
#temp <- list.files(input_dir, pattern ="data.3-preVKs-2.sup2.3sources.csv")
#temp <- list.files(input_dir, pattern = "preVKs[.]csv")
#temp <- #load your file or files
# Find the sample codes
#check<-unique(na.omit(as.numeric(unlist(strsplit(unlist(substr(temp,1,10)), "[^0-9]+")))))

# Group files (not used with TSV)

#list2env(lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), read.csv), envir = .GlobalEnv)
#dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))
#data <- bind_rows(dfs)
#TSV's been saved as csv and will used as it is


#plotVanKrevelen(test)
#plotVanKrevelen(test)
#test <-  filter(data, Nature=="Multi 1" & y=="3")



#data2 <- data2[order(data2[,"x"],data2[,"O.C"],data2[,"H.C"]), ]
##data  <- read.csv2(file = temp, header = TRUE, sep =";")
#data  <- read.csv2(file = temp, header = TRUE, sep ="\t")

#data <- read.table(file =temp, header = true, sep =";")
#dataset1 <- subset(data, (y=="3")==TRUE)
#test <- filter(data, (Nature=="Multi 1")==TRUE & (y=="3")==TRUE)
 


#dataset1 <- subset(data, !is.na(data$chem)==TRUE)
# Remove rows without a chem tag
#dataset1 <- subset(data, (y!='NA')==TRUE)
#dataset1 <- subset(data, !is.na(data$chem)==TRUE)
#dataset1 <- subset(data, (y=="3")==TRUE)
#test <- subset(dataset1, (Nature=="Multi 1")==TRUE)
#plotVanKrevelen <- function(data2, digit){
 # plotVanKrevelen <- function(data2){
  # Correlate cex to frequency of Van Krevelen coordinates
 #width <- table(data2$O.C,data2$H.C) 
#scatter2D (data2$O.C, data2$H.C, colvar = data2$x,
#           col = NULL, NAcol = "white", breaks = NULL,
#           colkey = NULL, clim = NULL, clab = NULL, 
#           alpha=0.25,
#           pch = 20, cex = width+0.5,
#           xlab="O/C ratio", ylab="H/C ratio",
#           CI = NULL, add = FALSE, plot = TRUE)		   	   
#}


