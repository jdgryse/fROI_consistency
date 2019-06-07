
#Local working directories
scripts <- "PATH/TO/SCRIPTS/"
copestor <- "PATH/TO/DATA/"
threshstor <- "PATH/TO/THRESHOLDED/MAPS/"
results <- "PATH/TO/RESULTS/"


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

#Storage arrays
results.null <- array(NA, dim=c(100,2))
colnames(results.null) <- c("uncnhst","fdr")
results.abt60 <- array(NA, dim=c(100,6))
colnames(results.abt60) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
results.abt75 <- array(NA, dim=c(100,6))
colnames(results.abt75) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
results.abt90 <- array(NA, dim=c(100,6))
colnames(results.abt90) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
results.abt95 <- array(NA, dim=c(100,6))
colnames(results.abt95) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
results.mLR60 <- array(NA, dim=c(100,4))
colnames(results.mLR60) <- c("0.88","1.5","3.68","8")
results.mLR75 <- array(NA, dim=c(100,4))
colnames(results.mLR75) <- c("0.88","1.5","3.68","8")
results.mLR90 <- array(NA, dim=c(100,4))
colnames(results.mLR90) <- c("0.88","1.5","3.68","8")
results.mLR95 <- array(NA, dim=c(100,4))
colnames(results.mLR95) <- c("0.88","1.5","3.68","8")




## ANALYSIS OF 10 COPES AND VARCOPES

for(i in 1:nrun) {

  print(i)

  ## NHST: UNCORRECTED AND FDR CORRECTED
  null <- readNIfTI(paste(threshstor, "NULL_run_",i,sep=""))
  fdr <- readNIfTI(paste(threshstor, "FDR_run_",i,sep=""))

  results.null[i,1] <- sum(null, na.rm=TRUE)
  results.null[i,2] <- sum(fdr, na.rm=TRUE)


  ## ABT
  for(j in 1:dim(conditions_abt)[1]) {
  
    perc <- conditions_abt[j,1]
    alpha <- conditions_abt[j,2]
    beta <- conditions_abt[j,3]

    abt <- readNIfTI(paste(threshstor,"ABT_delta_",perc,"_alpha_",alpha,"_beta_",beta,"_run_",i,sep=""))

    if(perc==60) {
      
      if(alpha == 0.05 && beta == 0.1)
        results.abt60[i,1] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.05 && beta == 0.2)
        results.abt60[i,2] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.05 && beta == 0.3)
        results.abt60[i,3] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.001 && beta == 0.1)
        results.abt60[i,4] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.001 && beta == 0.2)
        results.abt60[i,5] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.001 && beta == 0.3)
        results.abt60[i,6] <- sum(abt, na.rm=TRUE)
    
    }

    if(perc==75) {
      
      if(alpha == 0.05 && beta == 0.1)
        results.abt75[i,1] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.05 && beta == 0.2)
        results.abt75[i,2] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.05 && beta == 0.3)
        results.abt75[i,3] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.001 && beta == 0.1)
        results.abt75[i,4] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.001 && beta == 0.2)
        results.abt75[i,5] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.001 && beta == 0.3)
        results.abt75[i,6] <- sum(abt, na.rm=TRUE)
    
    }

    if(perc==90) {
      
      if(alpha == 0.05 && beta == 0.1)
        results.abt90[i,1] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.05 && beta == 0.2)
        results.abt90[i,2] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.05 && beta == 0.3)
        results.abt90[i,3] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.001 && beta == 0.1)
        results.abt90[i,4] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.001 && beta == 0.2)
        results.abt90[i,5] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.001 && beta == 0.3)
        results.abt90[i,6] <- sum(abt, na.rm=TRUE)
    
    }

    if(perc==95) {
      
      if(alpha == 0.05 && beta == 0.1)
        results.abt95[i,1] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.05 && beta == 0.2)
        results.abt95[i,2] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.05 && beta == 0.3)
        results.abt95[i,3] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.001 && beta == 0.1)
        results.abt95[i,4] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.001 && beta == 0.2)
        results.abt95[i,5] <- sum(abt, na.rm=TRUE)
      if(alpha == 0.001 && beta == 0.3)
        results.abt95[i,6] <- sum(abt, na.rm=TRUE)
    
    }

  }


  ## mLR

  for(j in 1:dim(conditions_mLR)[1]) {
  
    perc <- conditions_mLR[j,1]
    k <- conditions_mLR[j,2]

    mLR <- readNIfTI(paste(threshstor,"mLR_delta_",perc,"_k_",k,"_run_",i,sep=""))

    if(perc==60 && k==0.88)
      results.mLR60[i,1] <- sum(mLR, na.rm=TRUE)
    if(perc==60 && k==1.5)
      results.mLR60[i,2] <- sum(mLR, na.rm=TRUE)
    if(perc==60 && k==3.68)
      results.mLR60[i,3] <- sum(mLR, na.rm=TRUE)
    if(perc==60 && k==8)
      results.mLR60[i,4] <- sum(mLR, na.rm=TRUE)      

    if(perc==75 && k==0.88)
      results.mLR75[i,1] <- sum(mLR, na.rm=TRUE)
    if(perc==75 && k==1.5)
      results.mLR75[i,2] <- sum(mLR, na.rm=TRUE)
    if(perc==75 && k==3.68)
      results.mLR75[i,3] <- sum(mLR, na.rm=TRUE)
    if(perc==75 && k==8)
      results.mLR75[i,4] <- sum(mLR, na.rm=TRUE)

    if(perc==90 && k==0.88)
      results.mLR90[i,1] <- sum(mLR, na.rm=TRUE)
    if(perc==90 && k==1.5)
      results.mLR90[i,2] <- sum(mLR, na.rm=TRUE)
    if(perc==90 && k==3.68)
      results.mLR90[i,3] <- sum(mLR, na.rm=TRUE)
    if(perc==90 && k==8)
      results.mLR90[i,4] <- sum(mLR, na.rm=TRUE)

    if(perc==95 && k==0.88)
      results.mLR95[i,1] <- sum(mLR, na.rm=TRUE)
    if(perc==95 && k==1.5)
      results.mLR95[i,2] <- sum(mLR, na.rm=TRUE)
    if(perc==95 && k==3.68)
      results.mLR95[i,3] <- sum(mLR, na.rm=TRUE)
    if(perc==95 && k==8)
      results.mLR95[i,4] <- sum(mLR, na.rm=TRUE)
    
  }
  
}


## WRITE AWAY RESULTS

setwd(results)
write.csv(results.null, "voxnum_null_cond",row.names=FALSE)
write.csv(results.abt60, "voxnum_abt60_cond",row.names=FALSE)
write.csv(results.abt75, "voxnum_abt75_cond",row.names=FALSE)
write.csv(results.abt90, "voxnum_abt90_cond",row.names=FALSE)
write.csv(results.abt95, "voxnum_abt95_cond",row.names=FALSE)
write.csv(results.mLR60, "voxnum_mLR60_cond",row.names=FALSE)
write.csv(results.mLR75, "voxnum_mLR75_cond",row.names=FALSE)
write.csv(results.mLR90, "voxnum_mLR90_cond",row.names=FALSE)
write.csv(results.mLR95, "voxnum_mLR95_cond",row.names=FALSE)


## CALCULATE STANDARD DEVIATION OF THE DETECTED NUMBER OF VOXELS


results.null <- as.matrix(read.csv(paste(results,"voxnum_null_cond",sep="")))
results.abt60 <- as.matrix(read.csv(paste(results,"voxnum_abt60_cond",sep="")))
results.abt75 <- as.matrix(read.csv(paste(results,"voxnum_abt75_cond",sep="")))
results.abt90 <- as.matrix(read.csv(paste(results,"voxnum_abt90_cond",sep="")))
results.abt95 <- as.matrix(read.csv(paste(results,"voxnum_abt95_cond",sep="")))
results.mLR60 <- as.matrix(read.csv(paste(results,"voxnum_mLR60_cond",sep="")))
results.mLR75 <- as.matrix(read.csv(paste(results,"voxnum_mLR75_cond",sep="")))
results.mLR90 <- as.matrix(read.csv(paste(results,"voxnum_mLR90_cond",sep="")))
results.mLR95 <- as.matrix(read.csv(paste(results,"voxnum_mLR95_cond",sep="")))


vecsd.null <- colSds(results.null)
vecsd.abt60 <- colSds(results.abt60)
vecsd.abt75 <- colSds(results.abt75)
vecsd.abt90 <- colSds(results.abt90)
vecsd.abt95 <- colSds(results.abt95)
vecsd.mLR60 <- colSds(results.mLR60)
vecsd.mLR75 <- colSds(results.mLR75)
vecsd.mLR90 <- colSds(results.mLR90)
vecsd.mLR95 <- colSds(results.mLR95)

