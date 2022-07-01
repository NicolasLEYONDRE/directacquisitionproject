                      #######################################################################
                      #   THE FIRST PART OF THIS SCRIPT IS SIMILAR TO CHEMQUANT,            #
                      #   IF YOU HAVE ALREADY GENERATED THE DATASET W/ CHEMQUANT            #
                      #   PLEASE LOAD THE .RDATA FILE AND GO DIRECTLY TO THE PLOTTING PART  #
                      #######################################################################
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
	temp <- list.files(input_dir, pattern = "dataset_pre-extract_comp_publi2.tsv")
	temp2 <- list.files(input_dir, pattern = "DOJ.DB.MOD.simp.csv")

#for VKplot : 2 protocols in 1 :  experimental data (temp)  //the database (temp2)


################################################################################
#
#     generation of data for the experimental part
#
################################################################################
#reading function for tsv, csv
#take care of special characters like accent (we are French :-) )
	data <- read.csv2(file = temp, sep = '\t', header = TRUE)
#get rid of ESI data.
	data <- subset(data, x!=4)
#depending on file charateristics
#data <- read.csv2(file = temp, sep = ';', header = TRUE)
#data <- read.csv2(file = temp, sep = ',', header = TRUE)
#if the file is from OSF, it should be processed and saved with the right (chem, x, y) => goto "sort"
#Intitialisation of chem, y and x
	data$chem <- NA
	data$y <- NA
	data$x <- NA
# Compute chemical family according to Van Krevelen coordinates for all data
# one run for all 3 (chem, y and x)
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
	  if(grepl("DART", data$Ionisation[i])==TRUE){
		data$x[i] = 1
	  }
	  if(grepl("ASAP", data$Ionisation[i])==TRUE){
		data$x[i] = 2
	  }
	 # ## type DataBase 
	 # if(grepl("DataBase", data$Nature[i])==TRUE){
	 #	data$x[i] = 3
	 # }
	}


######################################################################################################################
#
#   Creating the plot
#   expe - 3 colors - > (1) alone - (2) alone - (1) and (2) 1 color.
#   Focus on y = 3

#
###############################################################################################

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

# one should first sort data
	data2 <- data[order(data[,"x"],data[,"O.C"],data[,"H.C"]), ]
	
# choose the right filters

# in this work the need of VKplot and its plotting

#  array of focalisation

#  publication 2
#	3c  : (((Nature=="Thalle")|(Nature=="Broyat"))&(y==3)) 
	data3 <- subset(data2, ((((Nature=="Thalle")|(Nature=="Broyat"))==TRUE)&(y==3)==TRUE))

#  publication 2 supplementary info
#	S2 : ((Nature=="Thalle")|(Nature=="Broyat"))
	data3 <- subset(data2, (((Nature=="Thalle")|(Nature=="Broyat"))==TRUE))	
	
# the frequence of apparition for a given couple O.C and H.C is needed
# this will correspond to the size of the colored dot or circle
# with the correlation with x for color layer.
# frequence of x, O.C, H.C
# pub 1 5 a/5 b /  pub 2 3 c
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
# addition of a dimension (present only in x==1 (red), present only in x==2 (blue), in both but 1 (3)in green and 2(4) in yellow)
	data4$zcolor = 0
	i<-1
	while(i<(nrow(data4)))
	{
		if((data4$O.C[i]==data4$O.C[i+1])&&(data4$H.C[i]==data4$H.C[i+1]))
		{
		#   double points found
		#  	line x   O.C H.C  Freq
		#	i    1	0.25 0.85   freq1
		#	i+1  2   0.25 0.85   freq2
		#	at the end the cursor is driven to i+2	
			if((data4$Freq[i]==data4$Freq[i+1]))
			{
			# freq1 == freq2 # 
			# same size = other color
			data4$zcolor[i] = 3
			data4$Freq[i] = 2*data4$Freq[i]
			i<-i+2
			}
			else
			{
				if((data4$Freq[i]>data4$Freq[i+1]))
				{
				# freq1 (x=1) > freq2 (x=2) #
				#data4$zcolor[i] = 4  # will be the color of x=1 (big and underneath)
				#data4$zcolor[i+1] = 5 # will be the color of x=2 (small and above)
				data4$zcolor[i] = 3
				data4$Freq[i] = data4$Freq[i]+data4$Freq[i+1]
				i<-i+2
				}
				else
				{
				# freq1 (x=1) < freq2 (x=2) #				
				#data4$zcolor[i] = 6 # will be the color of x=2 after swap (big and underneath)
				#data4$zcolor[i+1] = 7 #will be the color of x=1 after swap (small and above)
				data4$zcolor[i] = 3
				data4$Freq[i] = data4$Freq[i]+data4$Freq[i+1]
				# line swapping to put x=2 under x=1
				# the first ordering is still active
				# the same o/c and h/c
				# exchange X Freq and x
				#myFreq = data4$Freq[i+1]
				#data4$Freq[i+1] = data4$Freq[i]
				#data4$Freq[i] = myFreq
				#myx = data4$x[i+1]
				#data4$x[i+1] = data4$x[i]
				#data4$x[i] = myx				
				i<-i+2				
				}			
			}
		}
		else
		{
		#   single point found
		#  	line x      O.C H.C    Freq
		#	i    xi	   0.25 0.85   freq1
		#	i+1  xi+1  0.29 0.85   freq2
		#
		#	at the end the cursor is driven to i+1
			if((data4$x[i]==1))
			{
				# (x=1) #
				data4$zcolor[i] = 1
				i<-i+1
			}
			else
			{
				# (x=2) #
				data4$zcolor[i] = 2
				i<-i+1				
			}	
		}
		if (i==nrow(data4))
		{
		# last line before 
			if((data4$x[i]==1))
			{
				# (x=1) #
				data4$zcolor[i] = 1
				i<-i+1
			}
			else
			{
				# (x=2) #
				data4$zcolor[i] = 2
				i<-i+1				
			}	
		}
	}

#################################################
#
#  complete story
#
#################################################

data5 <- subset(data4, zcolor>=1&zcolor<=3)
plotVanKrevelen.expe(data5, 0.5,0.05,1, 3)

#################################################
#
#  focus on specific
#
#################################################

data5 <- subset(data4, zcolor>=1&zcolor<=2)
plotVanKrevelen.expe(data5,0.5, 0.1,0, 2)



##############################################
#
#           plotting function
#              from plot3D
#				expe data color=f(zcolor)
#             7 layers
# https://rdrr.io/cran/plot3D/man/scatter.html
#############################################

	plotVanKrevelen.expe <- function(myfrequence, alpha, radius, mini, nbcol){
	c2 <- "#0000FF"
	c3 <- "#336600"
	#c2 <- "#00FFFF"
	#c1 <- "#339900"	
	c1 <- "#FF0000"
	c4 <- "#00FFFF"
	c5 <- "#00FFFF"
	c6 <- "#00FFFF"
	c7 <- "#0000FF"
	if (nbcol=="2")
	{
	mycol = c(c1,c2)
	}
	if (nbcol=="3")
	{
	mycol = c(c1,c2,c3)
	}
	if (nbcol=="5")
	{
	mycol = c(c3,c4,c5,c6,c7)
	}
	
	if (nbcol=="7")
	{
	mycol = c(c1,c2,c3,c4,c5,c6,c7)
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





##############################################
#
#           plotting function
#              from plot3D
#				expe data color=f(zcolor)
# https://rdrr.io/cran/plot3D/man/scatter.html
#############################################

	plotVanKrevelen.expe <- function(myfrequence, radius){
	  
	  # Correlate cex to frequency of Van Krevelen coordinates
	  # ici mafrequence est un tableau global avec O.C, H.C, x et la fréquence en variable
	  # ainsi x va fournir la couleur, Freq le cex nécessaire
	  # le cex est 50 % de la frequence : c'est améliorable 
	  # il serait intéressant de définir une granulométrie en fonction du % de représentatitivité

	scatter2D (myfrequence$O.C, myfrequence$H.C, colvar = myfrequence$x,
#				col = ramp.col (col = c("#0000FF","#FF00FF","#339900","#0000FF","#FF00FF","#FF00FF","#0000FF"), n = 100, alpha = 1),#level existing complete overview
				col = ramp.col (col = c("#0000FF","#FF00FF"), n = 100, alpha = 1),#level existing in both with <> frequence
			   NAcol = "white", breaks = NULL,
			   colkey = NULL, clim = NULL, clab = NULL, 
			   alpha=0.5,
			   pch = 20, cex = radius*(myfrequence$Freq)+1,
			   xlab="O/C ratio", ylab="H/C ratio",
			   CI = NULL, add = FALSE, plot = TRUE)
	}
	plotVanKrevelen.expe(data4,0.1)




################################################################################
#
#     generation of data for the database
#
################################################################################


#input_dir <- # your path beware with / 
	input_dir <- 'C:/Users/CRMPO/Documents/2022/Simon.rapport.DATA/VK/TSV'
	setwd(input_dir)
	temp2 <- list.files(input_dir, pattern = "DOJ.DB.MOD.simp.csv")
#reading function for tsv, csv
#take care of special characters like accent (we are French :-) )
	data <- read.csv2(file = temp2, sep = ';', header = TRUE)
#depending on file charateristics
#data <- read.csv2(file = temp, sep = ';', header = TRUE)
#data <- read.csv2(file = temp, sep = ',', header = TRUE)
#if the file is from OSF, it should be processed and saved with the right (chem, x, y) => goto "sort"
#Intitialisation of chem, y and x
	data$chem <- NA
	data$y <- NA
	data$x <- NA
# Compute chemical family according to Van Krevelen coordinates for all data
# one run for all 3 (chem, y and x)
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
	 # ## type DART, ASAP, ESI Acétone, ESI Méthanol 
	 # if(grepl("DART", data$Ionisation[i])==TRUE){
	 #	data$x[i] = 1
	 # }
	 # if(grepl("ASAP", data$Ionisation[i])==TRUE){
	 #	data$x[i] = 2
	 # }
	 # type DataBase 
	  if(grepl("DataBase", data$Nature[i])==TRUE){
	 	data$x[i] = 1
	  }
	}

######################################################################################################################
#
#   Creating the plot  
#   Focus on all y
# 	db alone - 11 colors  - part 2
#
# frequence array with (x, O.C and H.C)
#
# each x-level will get an existence and a real color.
#
###############################################################################################

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

# one should first sort data
	data2 <- data[order(data[,"y"],data[,"O.C"],data[,"H.C"]), ]
	
# choose the right filters
	data2<- subset(data, y!="NA") 
	
# the frequence of apparition for a given couple O.C and H.C is needed
# this will correspond to the size of the colored dot or circle
# with the correlation with y for color layer.
# frequence of y, O.C, H.C
# pub 1 5 a/5 b /  pub 2 3 c
	frequence <- as.data.frame(table(data2$y, data2$O.C, data2$H.C))
# frequence is a mix of factor and classical variable
# it's not compatible with following process (scatter2d)
# small trick write and read
	write.csv2(frequence, file ="myfilefrequence.csv")
	data4 <- read.csv2(file ="myfilefrequence.csv", header=TRUE)
# always stay tidy
	unlink("myfilefrequence.csv")
# only data with existing O.C H.C are wished
	data4 <- subset(data4, Freq!=0)
	colnames(data4) <- c('X','y','O.C','H.C','Freq')
# logically the data is ready for plotting
	


##############################################
#
#           plotting function
#              from plot3D
#				database
# https://rdrr.io/cran/plot3D/man/scatter.html
#############################################

	plotVanKrevelen.DB <- function(myfrequence, radius){
	  
	  # Correlate cex to frequency of Van Krevelen coordinates
	  # ici mafrequence est un tableau global avec O.C, H.C, x et la fréquence en variable
	  # ainsi x va fournir la couleur, Freq le cex nécessaire
	  # le cex est 50 % de la frequence : c'est améliorable 
	  # il serait intéressant de définir une granulométrie en fonction du % de représentatitivité
	c1 <- "#8dcd48"
	c2 <- "#00ae4e"
	c3 <- "#7f7f7f"
	c4 <- "#556bef"
	c5 <- "#212121"
	c6 <- "#002ef0"
	c7 <- "#00aaf0"
	c8 <- "#ffbb00"
	c9 <- "#ff312a"
	c10 <- "#cc12ff"
	mycol = c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)

	
	scatter2D (myfrequence$O.C, myfrequence$H.C, colvar = myfrequence$y,
			   col = ramp.col (col = mycol, n = 100, alpha = 1),#level existing complete overview 
			   NAcol = "white", breaks = NULL,
			   colkey = NULL, clim = NULL, clab = NULL, 
			   alpha=1,
			   pch = 20, cex = radius*(myfrequence$Freq)+1,
			   xlab="O/C ratio", ylab="H/C ratio",
			   CI = NULL, add = FALSE, plot = TRUE)
	}
	plotVanKrevelen.DB(data4,0.02)








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


