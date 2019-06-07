
#Local working directories
scripts <- "PATH/TO/SCRIPTS/"
copestor <- "PATH/TO/DATA/"
threshstor <- "PATH/TO/THRESHOLDED/MAPS"
results <- "PATH/TO/RESULTS"

## PRELIMINARY

#Libraries
library(neuRosim)
library(oro.nifti)
library(AnalyzeFMRI)
library(lattice)

#Parameters
nrun <- 10
nscan <- 165 

#All combinations of runs
pairs <- t(combn(c(1:100),2,simplify=TRUE))

#Source files
conditions_mLR <- read.csv(paste(scripts,"conditions_mLR.txt",sep=""))
conditions_abt <- read.csv(paste(scripts,"conditions_abt.txt",sep=""))

#Function for Maitra similarity index
maitra <- function(image1, image2) {

  nvox1 <- sum(image1==1,na.rm=TRUE)
  nvox2 <- sum(image2==1,na.rm=TRUE)
  overlap <- sum(image1 == image2 & image1 == 1,na.rm=TRUE)
  
  if(nvox1 == 0 & nvox2 == 0) {
    cons <- 0
  } else {
    cons <- overlap/(nvox1+nvox2-overlap)
  }

  return(c(cons,nvox1,nvox2))

}


#Storage arrays
results.null <- array(NA, dim=c(dim(pairs)[1],1))
colnames(results.null) <- c("cons alpha=0.001")
results.fdr <- array(NA, dim=c(dim(pairs)[1],1))
colnames(results.fdr) <- c("cons qval=0.05")
results.abt60 <- array(NA, dim=c(dim(pairs)[1],6))
colnames(results.abt60) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
results.abt75 <- array(NA, dim=c(dim(pairs)[1],6))
colnames(results.abt75) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
results.abt90 <- array(NA, dim=c(dim(pairs)[1],6))
colnames(results.abt90) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
results.abt95 <- array(NA, dim=c(dim(pairs)[1],6))
colnames(results.abt95) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
results.mLR60 <- array(NA, dim=c(dim(pairs)[1],4))
colnames(results.mLR60) <- c("0.88","1.5","3.68","8")
results.mLR75 <- array(NA, dim=c(dim(pairs)[1],4))
colnames(results.mLR75) <- c("0.88","1.5","3.68","8")
results.mLR90 <- array(NA, dim=c(dim(pairs)[1],4))
colnames(results.mLR90) <- c("0.88","1.5","3.68","8")
results.mLR95 <- array(NA, dim=c(dim(pairs)[1],4))
colnames(results.mLR95) <- c("0.88","1.5","3.68","8")


##EVALUATE MSI FOR EVERY POSSIBLE PAIR OF TWO RUNS


for(i in 1:dim(pairs)[1]) {

  print(i)

  run1 <- pairs[i,1]
  run2 <- pairs[i,2]

  ## UNCORRECTED NULL + FDR

  null1 <- readNIfTI(paste(threshstor, "NULL_run_",run1,sep=""))
  null2 <- readNIfTI(paste(threshstor, "NULL_run_",run2,sep=""))
  fdr1 <- readNIfTI(paste(threshstor, "FDR_run_",run1,sep=""))
  fdr2 <- readNIfTI(paste(threshstor, "FDR_run_",run2,sep=""))

  results.null[i,1] <- maitra(null1,null2)[1]
  results.fdr[i,1] <- maitra(fdr1,fdr2)[1]


  ## ABT
  for(j in 1:dim(conditions_abt)[1]) {
  
    perc <- conditions_abt[j,1]
    alpha <- conditions_abt[j,2]
    beta <- conditions_abt[j,3]

    abt1 <- readNIfTI(paste(threshstor,"ABT_delta_",perc,"_alpha_",alpha,"_beta_",beta,"_run_",run1,sep=""))
    abt2 <- readNIfTI(paste(threshstor,"ABT_delta_",perc,"_alpha_",alpha,"_beta_",beta,"_run_",run2,sep=""))

    if(perc==60) {
      
      if(alpha == 0.05 && beta == 0.1)
        results.abt60[i,1] <- maitra(abt1,abt2)[1]
      if(alpha == 0.05 && beta == 0.2)
        results.abt60[i,2] <- maitra(abt1,abt2)[1]
      if(alpha == 0.05 && beta == 0.3)
        results.abt60[i,3] <- maitra(abt1,abt2)[1]
      if(alpha == 0.001 && beta == 0.1)
        results.abt60[i,4] <- maitra(abt1,abt2)[1]
      if(alpha == 0.001 && beta == 0.2)
        results.abt60[i,5] <- maitra(abt1,abt2)[1]
      if(alpha == 0.001 && beta == 0.3)
        results.abt60[i,6] <- maitra(abt1,abt2)[1]
    
    }

    if(perc==75) {
      
      if(alpha == 0.05 && beta == 0.1)
        results.abt75[i,1] <- maitra(abt1,abt2)[1]
      if(alpha == 0.05 && beta == 0.2)
        results.abt75[i,2] <- maitra(abt1,abt2)[1]
      if(alpha == 0.05 && beta == 0.3)
        results.abt75[i,3] <- maitra(abt1,abt2)[1]
      if(alpha == 0.001 && beta == 0.1)
        results.abt75[i,4] <- maitra(abt1,abt2)[1]
      if(alpha == 0.001 && beta == 0.2)
        results.abt75[i,5] <- maitra(abt1,abt2)[1]
      if(alpha == 0.001 && beta == 0.3)
        results.abt75[i,6] <- maitra(abt1,abt2)[1]
    
    }

    if(perc==90) {
      
      if(alpha == 0.05 && beta == 0.1)
        results.abt90[i,1] <- maitra(abt1,abt2)[1]
      if(alpha == 0.05 && beta == 0.2)
        results.abt90[i,2] <- maitra(abt1,abt2)[1]
      if(alpha == 0.05 && beta == 0.3)
        results.abt90[i,3] <- maitra(abt1,abt2)[1]
      if(alpha == 0.001 && beta == 0.1)
        results.abt90[i,4] <- maitra(abt1,abt2)[1]
      if(alpha == 0.001 && beta == 0.2)
        results.abt90[i,5] <- maitra(abt1,abt2)[1]
      if(alpha == 0.001 && beta == 0.3)
        results.abt90[i,6] <- maitra(abt1,abt2)[1]
    
    }

    if(perc==95) {
      
      if(alpha == 0.05 && beta == 0.1)
        results.abt95[i,1] <- maitra(abt1,abt2)[1]
      if(alpha == 0.05 && beta == 0.2)
        results.abt95[i,2] <- maitra(abt1,abt2)[1]
      if(alpha == 0.05 && beta == 0.3)
        results.abt95[i,3] <- maitra(abt1,abt2)[1]
      if(alpha == 0.001 && beta == 0.1)
        results.abt95[i,4] <- maitra(abt1,abt2)[1]
      if(alpha == 0.001 && beta == 0.2)
        results.abt95[i,5] <- maitra(abt1,abt2)[1]
      if(alpha == 0.001 && beta == 0.3)
        results.abt95[i,6] <- maitra(abt1,abt2)[1]
    
    }

  }




  ## mLR

  for(j in 1:dim(conditions_mLR)[1]) {
  
    perc <- conditions_mLR[j,1]
    k <- conditions_mLR[j,2]

    mLR1 <- readNIfTI(paste(threshstor,"mLR_delta_",perc,"_k_",k,"_run_",run1,sep=""))
    mLR2 <- readNIfTI(paste(threshstor,"mLR_delta_",perc,"_k_",k,"_run_",run2,sep=""))

    if(perc==60 && k==0.88)
      results.mLR60[i,1] <- maitra(mLR1,mLR2)[1]
    if(perc==60 && k==1.5)
      results.mLR60[i,2] <- maitra(mLR1,mLR2)[1]
    if(perc==60 && k==3.68)
      results.mLR60[i,3] <- maitra(mLR1,mLR2)[1]
    if(perc==60 && k==8)
      results.mLR60[i,4] <- maitra(mLR1,mLR2)[1]      

    if(perc==75 && k==0.88)
      results.mLR75[i,1] <- maitra(mLR1,mLR2)[1]
    if(perc==75 && k==1.5)
      results.mLR75[i,2] <- maitra(mLR1,mLR2)[1]
    if(perc==75 && k==3.68)
      results.mLR75[i,3] <- maitra(mLR1,mLR2)[1]
    if(perc==75 && k==8)
      results.mLR75[i,4] <- maitra(mLR1,mLR2)[1]

    if(perc==90 && k==0.88)
      results.mLR90[i,1] <- maitra(mLR1,mLR2)[1]
    if(perc==90 && k==1.5)
      results.mLR90[i,2] <- maitra(mLR1,mLR2)[1]
    if(perc==90 && k==3.68)
      results.mLR90[i,3] <- maitra(mLR1,mLR2)[1]
    if(perc==90 && k==8)
      results.mLR90[i,4] <- maitra(mLR1,mLR2)[1]

    if(perc==95 && k==0.88)
      results.mLR95[i,1] <- maitra(mLR1,mLR2)[1]
    if(perc==95 && k==1.5)
      results.mLR95[i,2] <- maitra(mLR1,mLR2)[1]
    if(perc==95 && k==3.68)
      results.mLR95[i,3] <- maitra(mLR1,mLR2)[1]
    if(perc==95 && k==8)
      results.mLR95[i,4] <- maitra(mLR1,mLR2)[1]

  }
}

## WRITE AWAY RESULTS

setwd(results)
write.csv(results.null, "NULL_maitra_perc",row.names=FALSE)
write.csv(results.fdr, "FDR_maitra_perc",row.names=FALSE)
write.csv(results.abt60, "ABT_maitra_perc_60",row.names=FALSE)
write.csv(results.abt75, "ABT_maitra_perc_75",row.names=FALSE)
write.csv(results.abt90, "ABT_maitra_perc_90",row.names=FALSE)
write.csv(results.abt95, "ABT_maitra_perc_95",row.names=FALSE)
write.csv(results.mLR60, "mLR_maitra_perc_60",row.names=FALSE)
write.csv(results.mLR75, "mLR_maitra_perc_75",row.names=FALSE)
write.csv(results.mLR90, "mLR_maitra_perc_90",row.names=FALSE)
write.csv(results.mLR95, "mLR_maitra_perc_95",row.names=FALSE)
