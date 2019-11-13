library(rjson)
library(readr)

# Set Directory
input_dir <- 'C:/Users/Simon Ollivier/Desktop/working directory'
setwd(input_dir)
input_files <- list.files(input_dir, pattern = "[.]CSV")

# Select files
for(i in 1:length(input_files)){

		#Ensuring that we select the deisotoped peak list
		if (grepl('_deisopeaklist',input_files[[i]])==TRUE){
		
		#Setting parameters	depending on instrument
		  ## Mass range & max number of results
		  if (grepl('ASAP_',input_files[[i]])==TRUE){
		    massRange="0.0017"
		    max.nresults=10
		  }
		  if (grepl('DART_',input_files[[i]])==TRUE){
		    massRange="0.0010"
		    max.nresults=10
		  }
		  if (grepl('ESI_',input_files[[i]])==TRUE){
		    massRange="0.0010"
		    max.nresults=10
		  }
		  if (grepl('MALDI_',input_files[[i]])==TRUE){
		    massRange="0.0050"
		    max.nresults=20
		  }
		
		## Remaining parameters for the query (e.g. unsaturation)	
					minUnsaturation="-1"
					integerUnsaturation=0
					numberOfResultsOnly=FALSE
					typedResult=FALSE		
		
		#Applying electron mass correction according to polarity 
			##Positive ionisation (MS_01 et MS_03) 
			if ((grepl('_MS_01',input_files[[i]])==TRUE)|(grepl('_MS_03',input_files[[i]])==TRUE)){	
			correction <- +0.000549
			}
			##Negative ionisation (MS_02 et MS_04)
			if ((grepl('_MS_02',input_files[[i]])==TRUE)|(grepl('_MS_04',input_files[[i]])==TRUE)){	
			correction <- -0.000549
			}
				
		#Reading m/z information from the peak list		
				csv_active<-read_csv(input_files[[i]])
				report<-NULL
				for (j in 1:length(csv_active$mz)){
					# Correction of electron mass
					target<- (csv_active$mz[j]+correction)
					
					# Setting the elemental ranges
					## Carbon
					  ### If isotopic contributions are not detected or the number of C is too low (i.e. < 10, S/N reduction), we work w/ DNP cut-offs
					  ### Also if the intensity was too high and had to be replaced by a high value during detection
				if((is.na(csv_active$estim_nC[j])==TRUE)|((csv_active$estim_nC[j]<10)==TRUE)|((csv_active$into[j]==1e10)==TRUE)){	
					if ((target<= 200)==TRUE)	{
					  Cmin=0
					  Cmax=15
					}				
					if (((target> 200)==TRUE)&((target<= 400)==TRUE))	{
					  Cmin=0
					  Cmax=30
					}
					if (((target> 400)==TRUE)&((target<= 600)==TRUE))	{
					  Cmin=0
					  Cmax=42
					}
					if (((target> 600)==TRUE)&((target<= 800)==TRUE))	{
					  Cmin=0
					  Cmax=56
					}
					if (((target> 800)==TRUE)&((target<= 1000)==TRUE))	{
					  Cmin=0
					  Cmax=66
					}
					if (((target> 1000)==TRUE)&((target<= 1500)==TRUE))	{
					  Cmin=0
					  Cmax=100
					}
				}
				    ### Otherwise we consider the calculated contribution +- 4
				if((is.na(csv_active$estim_nC[j])==FALSE)&((csv_active$estim_nC[j] >= 10)==TRUE)&((csv_active$into[j]==1e10)==FALSE)){
				    Cmin=csv_active$estim_nC[j]-4
				    Cmax=csv_active$estim_nC[j]+4
				}	
					
					## Hydrogen (low isotopic contribution of 2D, works with DNP cut-offs)
					if (((target<= 200)==TRUE))	{
					  Hmin=0
					  Hmax=31
					}				
					if (((target> 200)==TRUE)&((target<= 400)==TRUE))	{
					  Hmin=0
					  Hmax=59
					}
					if (((target> 400)==TRUE)&((target<= 600)==TRUE))	{
					  Hmin=0
					  Hmax=87
					}
					if (((target> 600)==TRUE)&((target<= 800)==TRUE))	{
					  Hmin=0
					  Hmax=109
					}
					if (((target> 800)==TRUE)&((target<= 1000)==TRUE))	{
					  Hmin=0
					  Hmax=127
					}
					if (((target> 1000)==TRUE)&((target<= 1500)==TRUE))	{
					  Hmin=0
					  Hmax=183
					}		
					
					## Nitrogen (low isotopic contribution of 15N, we consider that not detected != not present and allow 2 N + 1 for DART-MS+ cf [M+NH4]+ adduct)
					if(is.na(csv_active$estim_nN[j])==TRUE){
					  Nmin=0
					  if (grepl('DART_',input_files[[i]])==TRUE){
					  			if ((grepl('_MS_01',input_files[[i]])==TRUE)|(grepl('_MS_03',input_files[[i]])==TRUE)){	
								Nmax=3
								}
								if ((grepl('_MS_02',input_files[[i]])==TRUE)|(grepl('_MS_04',input_files[[i]])==TRUE)){
								Nmax=2
								}
					  }
					  if (grepl('ASAP_',input_files[[i]])==TRUE){
					  Nmax=2
					  }
					  if (grepl('ESI_',input_files[[i]])==TRUE){
					  Nmax=2
					  }
					  if (grepl('ESI_',input_files[[i]])==TRUE){
					  Nmax=2
					  }
					} 
					if(is.na(csv_active$estim_nN[j])==FALSE) {
					  Nmin=0
					  Nmax=csv_active$estim_nN[j]+1
					}
					
					## Oxygen (low isotopic contribution of 18O, works with DNP cut-offs)
					if (((target<= 200)==TRUE))	{
					  Omin=0
					  Omax=7
					}				
					if (((target> 200)==TRUE)&((target<= 400)==TRUE))	{
					  Omin=0
					  Omax=14
					}
					if (((target> 400)==TRUE)&((target<= 600)==TRUE))	{
					  Omin=0
					  Omax=21
					}
					if (((target> 600)==TRUE)&((target<= 800)==TRUE))	{
					  Omin=0
					  Omax=25
					}
					if (((target> 800)==TRUE)&((target<= 1000)==TRUE))	{
					  Omin=0
					  Omax=37
					}
					if (((target> 1000)==TRUE)&((target<= 1500)==TRUE))	{
					  Omin=0
					  Omax=44
					}
					
					## Sulfur (relatively lower isotopic contribution compared to halogens, do not force Smin)
					if(is.na(csv_active$estim_nS[j])==TRUE){
					  Smin=0
					  Smax=0
					} 
					if(is.na(csv_active$estim_nS[j])==FALSE) {
					  Smin=0
					  Smax=csv_active$estim_nS[j]
					}
					
					## For halogens (high M+2 contribution, we force estim_nX)
					## Chlorine
					if((csv_active$into[j]!=1e10)==TRUE){
					if(is.na(csv_active$estim_nCl[j])==TRUE){
					  Clmin=0
					  Clmax=0
					} 
					if(is.na(csv_active$estim_nCl[j])==FALSE) {
					  Clmin=csv_active$estim_nCl[j]
					  Clmax=csv_active$estim_nCl[j]
					}}
					if((csv_active$into[j]==1e10)==TRUE){
					  Clmin=0
					  Clmax=1
					}
					
					## Bromine
					if(is.na(csv_active$estim_nBr[j])==TRUE){
					  Brmin=0
					  Brmax=0
					} 
					if(is.na(csv_active$estim_nBr[j])==FALSE) {
					  Brmin=csv_active$estim_nBr[j]
					  Brmax=csv_active$estim_nBr[j]
					}
					
					## Adduct
					## Sodium (depending on ionisation source)
					if (grepl('ESI_',input_files[[i]])==TRUE|grepl('MALDI_',input_files[[i]])==TRUE){
					  # Positive ionisation
					  if ((grepl('_MS_01',input_files[[i]])==TRUE)|(grepl('_MS_03',input_files[[i]])==TRUE)){	
					    if (((target<= 1000)==TRUE))	{
					      Namin=0
					      Namax=1
					    }
					    ## Forcing 1 Na over m/z 1000
					    if (((target> 1000)==TRUE)&((target<= 1500)==TRUE))	{
					      Namin=1
					      Namax=1
					    }
					  }
					  # Negative ionisation
					  if ((grepl('_MS_02',input_files[[i]])==TRUE)|(grepl('_MS_04',input_files[[i]])==TRUE)){	
					    Namin=0
					    Namax=0
					  }
					}
					if (grepl('ASAP_',input_files[[i]])==TRUE|grepl('DART_',input_files[[i]])==TRUE){
					  Namin=0
					  Namax=0
					}
					
					## Potassium (ESI+)
					if (grepl('ESI_',input_files[[i]])==TRUE|grepl('MALDI_',input_files[[i]])==TRUE){
					  # Positive ionisation
					  if ((grepl('_MS_01',input_files[[i]])==TRUE)|(grepl('_MS_03',input_files[[i]])==TRUE)){	
					    Kmin=0
						#Kmax=1 ## Turned off because the script returned too many incoherent formula i.e. Na+K. Please choose one of the two according to your data.
						Kmax=0
					  }
					  # Negative ionisation
					  if ((grepl('_MS_02',input_files[[i]])==TRUE)|(grepl('_MS_04',input_files[[i]])==TRUE)){	
					    Kmin=0
					    Kmax=0
					  }
					}
					if (grepl('ASAP_',input_files[[i]])==TRUE|grepl('DART_',input_files[[i]])==TRUE){
					  Kmin=0
					  Kmax=0
					}
					
					# Sending query to chemcalc.org
					mfRange=paste(paste("C",paste(Cmin,Cmax,sep="-"),sep=""),paste("H",paste(Hmin,Hmax,sep="-"),sep=""),
					              paste("N",paste(Nmin,Nmax,sep="-"),sep=""),paste("O",paste(Omin,Omax,sep="-"),sep=""),
					              paste("S",paste(Smin,Smax,sep="-"),sep=""),paste("Cl",paste(Clmin,Clmax,sep="-"),sep=""),
					              paste("Br",paste(Brmin,Brmax,sep="-"),sep=""),paste("Na",paste(Namin,Namax,sep="-"),sep=""),
								  paste("K",paste(Kmin,Kmax,sep="-"),sep=""),sep=" ")
					query <- paste("http://www.chemcalc.org/chemcalc/em?monoisotopicMass=",target,"&mfRange=", mfRange, "&massRange=", massRange , "&integerUnsaturation=",integerUnsaturation, "&minUnsaturation=",minUnsaturation, "&numberOfResultsOnly=",numberOfResultsOnly, "&typedResult=", typedResult  , sep="")
					resJSON <- fromJSON(paste(readLines(query), collapse=""))
					if (length(resJSON$results)!=0){
					  
					  nbMax = length(resJSON$results)	
					  if (nbMax >= max.nresults){nbMax=max.nresults}
					  
					  for (mf.index in 1:nbMax){	
					    result <- data.frame( "j"=j , "mz"=csv_active$mz[j] , "into"=csv_active$into[j] , "mf"=resJSON$results[[mf.index]]$mf , "em"=resJSON$results[[mf.index]]$em , "error"=resJSON$results[[mf.index]]$error , "unsaturation"=resJSON$results[[mf.index]]$unsat,stringsAsFactors = FALSE)
					    unsaturation	<-resJSON$results[[mf.index]]$unsat
					    my.mass = resJSON$results[[mf.index]]$em 
					    my.unsat = ((abs(2*unsaturation))%%2)
					    if (grepl('ESI_',input_files[[i]])==TRUE){
					      if (my.mass > 0) {
					        if (my.unsat==1){
					          ifelse(!exists('report'),report <- result,  report <- rbind(report,result,stringsAsFactors = FALSE))
					        }
					      }
					    }	
					    if (grepl('ESI_',input_files[[i]])==FALSE){
					      if (my.mass > 0) {
					        ifelse(!exists('report'),report <- result,  report <- rbind(report,result,stringsAsFactors = FALSE))
					      }	
					    }
					    
					  }
					}
					else
					{					
					  result <- data.frame( "j"=j , "mz"=csv_active$mz[j] , "into"=csv_active$into[j] , "mf"="-" , "em"="-" , "error"="-" , "unsaturation"="-",stringsAsFactors = FALSE)
					  ifelse(!exists('report'),report <- result,  report <- rbind(report,result,stringsAsFactors = FALSE))	
					}
				} 
				MFCSVfile <- paste(substr(input_files[[i]],1,nchar(input_files[[i]])-18), "_RawMF.CSV",sep ="")
				write.csv(report, file = MFCSVfile)
		}
}