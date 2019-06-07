
#Local working directories
scripts <- "PATH/TO/SCRIPTS/"
copestor <- "PATH/TO/DATA/"
threshstor <- "PATH/TO/THRESHOLDED/MAPS/"
resultstor <- "PATH/TO/RESULTS/"


## PRELIMINARY

#Libraries
library(neuRosim)
library(oro.nifti)
library(AnalyzeFMRI)
library(lattice)

#Source files
conditions_mLR <- read.csv(paste(scripts,"conditions_mLR.txt",sep=""))
conditions_abt <- read.csv(paste(scripts,"conditions_abt.txt",sep=""))

#Parameters
nrun <- 10
nscan <- 165 


## SUM OF ALL THRESHOLDED MAPS

null <- readNIfTI(paste(threshstor, "NULL_run_",1,sep=""))

for(i in 2:nrun) {

  null <- null + readNIfTI(paste(threshstor, "NULL_run_",i,sep=""))

}

writeNIfTI(null, paste(resultstor,"NULL_coherencemap",sep=""),gzipped=FALSE)


fdr <- readNIfTI(paste(threshstor, "FDR_run_",1,sep=""))

for(i in 2:nrun) {

  fdr <- fdr + readNIfTI(paste(threshstor, "FDR_run_",i,sep=""))

}

writeNIfTI(fdr, paste(resultstor,"FDR_coherencemap",sep=""),gzipped=FALSE)



for(j in 1:dim(conditions_abt)[1]) {

	perc <- conditions_abt[j,1]
    alpha <- conditions_abt[j,2]
    beta <- conditions_abt[j,3]

	abt <- readNIfTI(paste(threshstor,"ABT_delta_",perc,"_alpha_",alpha,"_beta_",beta,"_run_",1,sep=""))

	for(i in 2:nrun) {

  		abt <- abt + readNIfTI(paste(threshstor, "ABT_delta_",perc,"_alpha_",alpha,"_beta_",beta,"_run_",i,sep=""))

	}

	writeNIfTI(abt, paste(resultstor,"ABT_delta_",perc,"_alpha_",alpha,"_beta_",beta,"_coherencemap",sep=""),gzipped=FALSE)

}



for(j in 1:dim(conditions_mLR)[1]) {

	perc <- conditions_mLR[j,1]
    k <- conditions_mLR[j,2]

	mLR <- readNIfTI(paste(threshstor,"mLR_delta_",perc,"_k_",k,"_run_",1,sep=""))

	for(i in 2:nrun) {

  		mLR <- mLR + readNIfTI(paste(threshstor, "mLR_delta_",perc,"_k_",k,"_run_",i,sep=""))

	}

	writeNIfTI(mLR, paste(resultstor,"mLR_delta_",perc,"_k_",k,"_coherencemap",sep=""),gzipped=FALSE)

}

