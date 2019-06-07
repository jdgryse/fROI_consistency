
#Local working directories
scripts <- "PATH/TO/SCRIPTS/"
copestor <- "PATH/TO/DATA/"
threshstor <- "PATH/TO/THRESHOLDED/MAPS"
resultstor <- "PATH/TO/RESULTS"

## PRELIMINARY

#Libraries
library(neuRosim)
library(oro.nifti)
library(AnalyzeFMRI)

#Parameters
nrun <- 10
nscan <- 165 

#Source files
conditions_mLR <- read.csv(paste(scripts,"conditions_mLR.txt",sep=""))
conditions_abt <- read.csv(paste(scripts,"conditions_abt.txt",sep=""))

#Storage matrix
results.null <- array(NA, dim=c(1,2))
colnames(results.null) <- c("uncnhst","fdr")
results.abt60 <- array(NA, dim=c(1,6))
colnames(results.abt60) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
results.abt75 <- array(NA, dim=c(1,6))
colnames(results.abt75) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
results.abt90 <- array(NA, dim=c(1,6))
colnames(results.abt90) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
results.abt95 <- array(NA, dim=c(1,6))
colnames(results.abt95) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
results.mLR60 <- array(NA, dim=c(1,4))
colnames(results.mLR60) <- c("0.88","1.5","3.68","8")
results.mLR75 <- array(NA, dim=c(1,4))
colnames(results.mLR75) <- c("0.88","1.5","3.68","8")
results.mLR90 <- array(NA, dim=c(1,4))
colnames(results.mLR90) <- c("0.88","1.5","3.68","8")
results.mLR95 <- array(NA, dim=c(1,4))
colnames(results.mLR95) <- c("0.88","1.5","3.68","8")



## COMPUTE COHERENCE FOR THE THREE METHODS

## NHST
null <- readNIfTI(paste(resultstor, "NULL_coherencemap",sep=""))
null_start <- c(EMbinom(Y=na.omit(c(null)),N=nrun,iniL=0.1, iniPI1=0.6, iniPI2=0.05, max.iter=500, tolerance=0.0001))
results.null[1,1] <- CohenKappa(null_start$lambda, null_start$PI1,null_start$PI2)

fdr <- readNIfTI(paste(resultstor, "FDR_coherencemap",sep=""))
fdr_start <- c(EMbinom(Y=na.omit(c(fdr)),N=nrun,iniL=0.1, iniPI1=0.6, iniPI2=0.05, max.iter=500, tolerance=0.0001))
results.null[1,2] <- CohenKappa(fdr_start$lambda, fdr_start$PI1,fdr_start$PI2)


##ABT
for(j in 1:dim(conditions_abt)) {

    perc <- conditions_abt[j,1]
    alpha <- conditions_abt[j,2]
    beta <- conditions_abt[j,3]

    abt <- readNIfTI(paste(resultstor,"ABT_delta_",perc,"_alpha_",alpha,"_beta_",beta,"_coherencemap",sep=""))
	  abt_start <- c(EMbinom(Y=na.omit(c(abt)),N=nrun,iniL=0.1, iniPI1=0.6, iniPI2=0.05, max.iter=500, tolerance=0.0001))


   if(perc==60) {
      
      if(alpha == 0.05 && beta == 0.1)
        results.abt60[1,1] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.05 && beta == 0.2)
        results.abt60[1,2] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.05 && beta == 0.3)
        results.abt60[1,3] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.001 && beta == 0.1)
        results.abt60[1,4] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.001 && beta == 0.2)
        results.abt60[1,5] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.001 && beta == 0.3)
        results.abt60[1,6] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
    
    }

    if(perc==75) {
      
      if(alpha == 0.05 && beta == 0.1)
        results.abt75[1,1] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.05 && beta == 0.2)
        results.abt75[1,2] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.05 && beta == 0.3)
        results.abt75[1,3] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.001 && beta == 0.1)
        results.abt75[1,4] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.001 && beta == 0.2)
        results.abt75[1,5] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.001 && beta == 0.3)
        results.abt75[1,6] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      
    }

    if(perc==90) {
      
      if(alpha == 0.05 && beta == 0.1)
        results.abt90[1,1] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.05 && beta == 0.2)
        results.abt90[1,2] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.05 && beta == 0.3)
        results.abt90[1,3] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.001 && beta == 0.1)
        results.abt90[1,4] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.001 && beta == 0.2)
        results.abt90[1,5] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.001 && beta == 0.3)
        results.abt90[1,6] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      
    }

    if(perc==95) {
      
      if(alpha == 0.05 && beta == 0.1)
        results.abt95[1,1] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.05 && beta == 0.2)
        results.abt95[1,2] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.05 && beta == 0.3)
        results.abt95[1,3] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.001 && beta == 0.1)
        results.abt95[1,4] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.001 && beta == 0.2)
        results.abt95[1,5] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      if(alpha == 0.001 && beta == 0.3)
        results.abt95[1,6] <- CohenKappa(abt_start$lambda, abt_start$PI1,abt_start$PI2)
      
    }


}


##mLR

for(j in 1:dim(conditions_mLR)[1]) {
  
    perc <- conditions_mLR[j,1]
    k <- conditions_mLR[j,2]

	mLR <- readNIfTI(paste(resultstor, "mLR_delta_",perc,"_k_",k,"_coherencemap",sep=""))
	mLR_start <- c(EMbinom(Y=na.omit(c(mLR)),N=nrun,iniL=0.1, iniPI1=0.6, iniPI2=0.05, max.iter=500, tolerance=0.0001))

	if(perc==60 && k==0.88)
	  results.mLR60[1,1] <- CohenKappa(mLR_start$lambda, mLR_start$PI1,mLR_start$PI2)
	if(perc==60 && k==1.5)
	  results.mLR60[1,2] <- CohenKappa(mLR_start$lambda, mLR_start$PI1,mLR_start$PI2)
	if(perc==60 && k==3.68)
	  results.mLR60[1,3] <- CohenKappa(mLR_start$lambda, mLR_start$PI1,mLR_start$PI2)
	if(perc==60 && k==8)
	  results.mLR60[1,4] <- CohenKappa(mLR_start$lambda, mLR_start$PI1,mLR_start$PI2)     

	if(perc==75 && k==0.88)
	  results.mLR75[1,1] <- CohenKappa(mLR_start$lambda, mLR_start$PI1,mLR_start$PI2)
	if(perc==75 && k==1.5)
	  results.mLR75[1,2] <- CohenKappa(mLR_start$lambda, mLR_start$PI1,mLR_start$PI2)
	if(perc==75 && k==3.68)
	  results.mLR75[1,3] <- CohenKappa(mLR_start$lambda, mLR_start$PI1,mLR_start$PI2)
	if(perc==75 && k==8)
	  results.mLR75[1,4] <- CohenKappa(mLR_start$lambda, mLR_start$PI1,mLR_start$PI2)     

	if(perc==90 && k==0.88)
	  results.mLR90[1,1] <- CohenKappa(mLR_start$lambda, mLR_start$PI1,mLR_start$PI2)
	if(perc==90 && k==1.5)
	  results.mLR90[1,2] <- CohenKappa(mLR_start$lambda, mLR_start$PI1,mLR_start$PI2)
	if(perc==90 && k==3.68)
	  results.mLR90[1,3] <- CohenKappa(mLR_start$lambda, mLR_start$PI1,mLR_start$PI2)
	if(perc==90 && k==8)
	  results.mLR90[1,4] <- CohenKappa(mLR_start$lambda, mLR_start$PI1,mLR_start$PI2)     

	if(perc==95 && k==0.88)
	  results.mLR95[1,1] <- CohenKappa(mLR_start$lambda, mLR_start$PI1,mLR_start$PI2)
	if(perc==95 && k==1.5)
	  results.mLR95[1,2] <- CohenKappa(mLR_start$lambda, mLR_start$PI1,mLR_start$PI2)
	if(perc==95 && k==3.68)
	  results.mLR95[1,3] <- CohenKappa(mLR_start$lambda, mLR_start$PI1,mLR_start$PI2)
	if(perc==95 && k==8)
	  results.mLR95[1,4] <- CohenKappa(mLR_start$lambda, mLR_start$PI1,mLR_start$PI2)     

}

