library("xcms")
library("msdata")
library("MSnbase")
library("CAMERA")
library("readMzXmlData")

# Set directory
setwd("C:/Users/Simon Ollivier/Desktop/working directory")
input_dir <- 'C:/Users/Simon Ollivier/Desktop/working directory'
input_files <- list.files(input_dir, pattern = "[.]mzXML")

# Detection parameters
for (i in 1:length(input_files)){
  file <- readMzXmlFile(input_files[[i]],removeMetaData = FALSE)
  
# For Thermo F. Q-Exactive
  if(grepl('FTMS',file[[1]]$metaData$msInstrument$msMassAnalyzer)==TRUE){
p <- MassifquantParam(ppm = 15, peakwidth = c(15, 150), snthresh = 2,
                      prefilter = c(3, 100), mzCenterFun = "wMean", integrate = 1L,
                      mzdiff = -0.001, fitgauss = FALSE, noise = 0,
                      verboseColumns = FALSE, criticalValue = 1.125,
                      consecMissedLimit = 50, unions = 1, checkBack = 0,
                      withWave = FALSE)
  }
  
## For Bruker MaXis 4G
  if(grepl('Quadrupole TOF',file[[1]]$metaData$msInstrument$msMassAnalyzer)==TRUE){
p <- MassifquantParam(ppm = 15, peakwidth = c(15, 150), snthresh = 2,
                      prefilter = c(5, 2000), mzCenterFun = "wMean", integrate = 1L,
                      mzdiff = -0.001, fitgauss = FALSE, noise = 0,
                      verboseColumns = FALSE, criticalValue = 1.125,
                      consecMissedLimit = 50, unions = 1, checkBack = 0,
                      withWave = FALSE)  
  }
  
register(SerialParam())

# Generates Peak Lists
monfichier <- readMSData(input_files[[i]], msLevel. = 1, mode = "onDisk")
monfichier2 <- findChromPeaks(monfichier, param = p, return.type = "xcmsSet")
monfichier3 <- xsAnnotate(monfichier2)
#monfichier4 <- groupFWHM(monfichier3, perfwhm = 100) ## Chromatographic parameter
monfichier5 <- findIsotopes (monfichier3, mzabs = 0.01)
#monfichier6 <- groupCorr(monfichier5, cor_eic_th = 0.75) ## Chromatographic parameter
#monfichier7 <- findAdducts(monfichier6, polarity="positive") ## Only considers ESI ## Only considers LC data whereas DIMS includes more [M+Na]+
peaklist <- getPeaklist(monfichier5)

# Creates a deisotoped peak list with isotope cluster analysis
  
  ## First step : Generating an appropriate matrix
  basenames <- c("is_monoiso", "estim_nC", "estim_nN", "estim_nO", "estim_nS", "estim_nCl", "estim_nBr")

  isopeaklist_gen <- matrix(ncol=length(basenames))
  colnames(isopeaklist_gen) <- c(basenames)
  
  ## Second step : Pasting the data from the original peaklist in our matrix
  data.matrix(peaklist)
  mat.list <- list(peaklist,isopeaklist_gen)
  df.list <- lapply(mat.list, as.data.frame)
  cat.df <- function(d1,d2) {d1[names(d2)] <- d2; d1}
  as_one.df <- Reduce(cat.df, df.list)
  isopeaklist_mat <- data.matrix(as_one.df)
  isopeaklist <- isopeaklist_mat[do.call("order", as.data.frame(isopeaklist_mat)),]
  
  ## Third step : Isotope clustering
  print("Calculating isotopic contributions...")
  
  #Setting polarity for electron mass correction
  if(grepl('+',file[[1]]$metaData$polarity)==TRUE){
    correction <- +0.000549 
  }
  if(grepl('-',file[[1]]$metaData$polarity)==TRUE){
    correction <- -0.000549 
  }
  
  #Setting parameters
  for (p in (1:nrow(isopeaklist))){
    for (q in (1:nrow(isopeaklist))){

      mz_M <- as.numeric(isopeaklist[p,"mz"])
      mz_M1 <- as.numeric(isopeaklist[q,"mz"])
      into_M <- as.numeric(isopeaklist[p,"into"])
      into_M1 <- as.numeric(isopeaklist[q,"into"])
      
	# Set ppm cut-off
	  ppm_value <- 3
      deiso_ppm <- 1e-6 * mz_M * ppm_value
  
  # Set maximum intensity for isotope contribution depanding on m/z (see Table)
    ## Correction of electron mass
      
        target <- mz_M+correction 
     
    ## Setting upper limit for the estimated number of elements depending on m/z
    if((target <= 200)==TRUE){
      Cmax <- 15
      Nmax <- 8
      Omax <- 7
      Smax <- 6
      Clmax <- 4
      Brmax <- 2
    } 
    
    if(((target>200)==TRUE)&((target<=400)==TRUE)){
      Cmax <- 30
      Nmax <- 10
      Omax <- 14
      Smax <- 12
      Clmax <- 7
      Brmax <- 4      
    }

    if(((target>400)==TRUE)&((target<=600)==TRUE)){
      Cmax <- 42
      Nmax <- 13
      Omax <- 21
      Smax <- 12
      Clmax <- 8
      Brmax <- 6      
    }  
    
    if(((target>600)==TRUE)&((target<=800)==TRUE)){
      Cmax <- 56
      Nmax <- 16
      Omax <- 25
      Smax <- 20
      Clmax <- 10
      Brmax <- 8      
    }
    
    if(((target>800)==TRUE)&((target<=1000)==TRUE)){
      Cmax <- 66
      Nmax <- 25
      Omax <- 37
      Smax <- 20
      Clmax <- 11
      Brmax <- 8      
    }
    
    if(((target>1000)==TRUE)&((target<=1500)==TRUE)){
      Cmax <- 100
      Nmax <- 26
      Omax <- 44
      Smax <- 20
      Clmax <- 11
      Brmax <- 8      
    }
    
    ## Implementing isotopic abundance
    percent_M1 <- (into_M1/into_M)*100
    est_nC <- percent_M1/1.10800
    est_nO <- percent_M1/0.20004
    est_nN <- percent_M1/0.36630
    est_nS33 <- percent_M1/0.74869
    est_nS34 <- percent_M1/4.19599
    est_nCl <- percent_M1/24.23530
    est_nBr <- percent_M1/49.31400
    
	# Set relative intensity threshold for research of isotopic cluster
	
    bp_into <- max(as.numeric(isopeaklist[p,"into"]))
    q_thresh <- 0.01/100*bp_into
    
    ## If peak intensity is over detection limit assign value > detected max
    if(is.na(into_M)==TRUE){isopeaklist[p,"into"]=1e10}
    if(is.na(into_M1)==TRUE){isopeaklist[q,"into"]=1e10}
	
	# Isotope counter
    
########################################################################################################################
    # WARNING: DO NOT ATTEMPT TO IMPLEMENT AN 'IFELSE' LOOP!!! THE COMPARISON WITH THE OTHER PEAKS OF THE MATRIX WOULD
    # RETURN AN ESTIMATED VALUE OF 0!!
########################################################################################################################
    
  ## Carbon
     if (((mz_M + 1.003355 - deiso_ppm <= mz_M1)==TRUE)&((mz_M1 <= mz_M + 1.003355 + deiso_ppm)==TRUE)) {
       if(((into_M1 >= q_thresh)==TRUE)&((into_M >= q_thresh)==TRUE)&((est_nC <= Cmax)==TRUE)) {
       isopeaklist[q,"is_monoiso"] = "13C"
       isopeaklist[p,"estim_nC"] = round(est_nC)
       }
     }
  ## Oxygen
    if (((mz_M + 2.004245 - deiso_ppm <= mz_M1)==TRUE)&((mz_M1 <= mz_M + 2.004245 + deiso_ppm)==TRUE)) {
      if(((into_M1 >= q_thresh)==TRUE)&((into_M >= q_thresh)==TRUE)&((est_nO <= Omax)==TRUE)) {
        isopeaklist[q,"is_monoiso"] = "18O"
        isopeaklist[p,"estim_nO"] = round(est_nO)
      }
    }
	## Nitrogen
    if (((mz_M + 0.997035 - deiso_ppm <= mz_M1)==TRUE)&((mz_M1 <= mz_M + 0.997035 + deiso_ppm)==TRUE)) {
      if(((into_M1 >= q_thresh)==TRUE)&((into_M >= q_thresh)==TRUE)&((est_nN <= Nmax)==TRUE)) {
        isopeaklist[q,"is_monoiso"] = "15N"
        isopeaklist[p,"estim_nN"] = round(est_nN)
      }
    }
	## Sulphur
    ### 34S
    if (((mz_M + 1.995796 - deiso_ppm <= mz_M1)==TRUE)&((mz_M1 <= mz_M + 1.995796 + deiso_ppm)==TRUE)) {
      if(((into_M1 >= q_thresh)==TRUE)&((into_M >= q_thresh)==TRUE)&((est_nS34 <= Smax)==TRUE)) {
        isopeaklist[q,"is_monoiso"] = "34S"
        isopeaklist[p,"estim_nS"] = round(est_nS34)
      } 
    } 
    ### 33S (might be detected as well if the number of S atoms is high, but is least adapted for calculation of estim_nS)
    if (((mz_M + 0.999388 - deiso_ppm <= mz_M1)==TRUE)&((mz_M1 <= mz_M + 0.999388 + deiso_ppm)==TRUE)) {
      if(((into_M1 >= q_thresh)==TRUE)&((into_M >= q_thresh)==TRUE)&((est_nS33 <= Smax)==TRUE)) {
        isopeaklist[q,"is_monoiso"] = "33S"
      }
      }
	## Halogens
    #### Bromine (unlock code if necessary)
#      if (((mz_M + 1.997952 - deiso_ppm <= mz_M1)==TRUE)&((mz_M1 <= mz_M + 1.997952 + deiso_ppm)==TRUE)) {
#        if(((into_M1 >= q_thresh)==TRUE)&((into_M >= q_thresh)==TRUE)&((est_nBr <= Brmax)==TRUE)) {
#          isopeaklist[q,"is_monoiso"] = "81Br"
#          isopeaklist[p,"estim_nBr"] = round(est_nBr)
#        } 
#      } 
    ### Chlorine
      if (((mz_M + 1.997050 - deiso_ppm <= mz_M1)==TRUE)&((mz_M1 <= mz_M + 1.997050 + deiso_ppm)==TRUE)) {
        if(((into_M1 >= q_thresh)==TRUE)&((into_M >= q_thresh)==TRUE)&((est_nCl <= Clmax)==TRUE)) {
          isopeaklist[q,"is_monoiso"] = "37Cl"
          isopeaklist[p,"estim_nCl"] = round(est_nCl)
        } 
      }}}
  
  print("Calculated isotopic contributions!")
  
  for(x in (1:nrow(isopeaklist))){
    if(is.na(isopeaklist[x,"is_monoiso"]==TRUE)==TRUE){isopeaklist[x,"is_monoiso"]=1}
  }
      ## Creates a CSV export file (unlock to print peaklist w/ identified isotopic clusters)
          
     # CSVfile <- paste(substr(input_files[[i]],1,nchar(input_files[[i]])-6), "_isopeaklist.CSV",sep ="")
     # write.csv(isopeaklist, file = CSVfile)
    
	## Fourth step : Deisotoping (i.e. 'is_monoiso' stays NA) 

    deisopeaklist <- subset(isopeaklist, isopeaklist[,"is_monoiso"] == 1)

# Apply a lock-mass calibration on the resulting peak list w/ higher accuracy than constructor software 
# for Bruker maXis (lock on palmitic acid C16H32O2)       
    ## Instrument
    if(grepl('Quadrupole TOF',file[[1]]$metaData$msInstrument$msMassAnalyzer)==TRUE){
      
      ##Polarity
      if(grepl('+',file[[1]]$metaData$polarity)==TRUE){
        lock.ref=257.247507
      }
      if(grepl('-',file[[1]]$metaData$polarity)==TRUE){
        lock.ref=255.232955
      }
      a<-which(abs(as.numeric(deisopeaklist[,"mz"])-lock.ref)==min(abs(as.numeric(deisopeaklist[,"mz"])-lock.ref)))
      lock.mes<-as.numeric(deisopeaklist[a,"mz"])
      lock.into<-as.numeric(deisopeaklist[a,"into"])
      c <- lock.ref/lock.mes
      
      for(z in 1:nrow(deisopeaklist)){
        mz.mes <- as.numeric(deisopeaklist[z,"mz"])
        mz.into <- as.numeric(deisopeaklist[z,"into"])
        deisopeaklist[z,"mz"] = c*mz.mes
      }
    }
    
## Creates a CSV export file
CSVfile <- paste(substr(input_files[[i]],1,nchar(input_files[[i]])-6), "_deisopeaklist.CSV",sep ="")
write.csv(deisopeaklist, file = CSVfile)

  }
