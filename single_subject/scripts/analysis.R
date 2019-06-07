
#Local working directories
scripts <- "PATH/TO/SCRIPTS/"
copestor <- "PATH/TO/DATA/"
threshstor <- "PATH/TO/THRESHOLDED/MAPS/"


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
tau <- 0.5
qval <- 0.05
nrun <- 10
nscan <- 165


## ANALYSIS OF 10 COPES AND VARCOPES

for(i in 1:nrun) {

  print(i)

  #DELTA CALCULATION

  EStruth <- readNIfTI(paste(copestor,"truth_ES_",i,".nii.gz",sep=""))
  mask <- readNIfTI(paste(copestor,"real_OVERALLMASK.nii.gz",sep=""))
  mask <- ifelse(mask==0, NA, mask)

  truth <- EStruth * mask

  #Define the effect size of interest based on percentile of observed contrast estimates
  es1 <- quantile(truth, 0.60, na.rm=TRUE)
  es2 <- quantile(truth, 0.75,na.rm=TRUE)
  es3 <- quantile(truth, 0.90,na.rm=TRUE)
  es4 <- quantile(truth, 0.95,na.rm=TRUE)

  perc <- c(60,75,90,95)

  ## READING IN B1 and SB1 files

  cope <- readNIfTI(paste(copestor, "cope_",i,".nii.gz",sep=""),reorient=FALSE)
  varcope <- readNIfTI(paste(copestor, "varcope_",i,".nii.gz",sep=""))
  tmap <- cope/sqrt(varcope)

  ## NHST uncorrected
  alpha <- 0.001
  pmap.null <- pt(tmap, nscan-2,lower.tail=FALSE)
  active.null <- ifelse(pmap.null <= alpha, 1, 0)

  writeNIfTI(active.null, paste(threshstor,"NULL_run_",i,sep=""),gzipped=TRUE)delta


  ## NHST FDR CORRECTED

  pvalms <- sort(pmap.null)
  orderpvalms <- rank(pmap.null)
  FDRqval <- (orderpvalms/(length(pmap.null)))*qval
  pr <- ifelse(pvalms[orderpvalms] < FDRqval,1,0)
  pthres.FDR<-ifelse(sum(pr,na.rm=TRUE)==0,0,max(FDRqval[pr==1],na.rm=TRUE))
  active.fdr <- ifelse(pmap.null <= pthres.FDR, 1, 0)

  writeNIfTI(active.fdr, paste(threshstor,"FDR_run_",i,sep=""),gzipped=TRUE)


  ## ABT

  for(j in 1:dim(conditions_abt)[1]) {
  
    perc <- conditions_abt[j,1]
    alpha <- conditions_abt[j,2]
    beta <- conditions_abt[j,3]

    if(perc==60)
      mudelta <- es1
    if(perc==75)
      mudelta <- es2    
    if(perc==90)
      mudelta <- es3
    if(perc==95)
      mudelta <- es4

    pmap.alt <- pnorm(tmap, (mudelta/sqrt(varcope)), sqrt(((varcope + tau^2)/varcope)))
    active.abt <- ifelse(pmap.null <= alpha & pmap.alt >= beta, 1, 0)

    writeNIfTI(active.abt, paste(threshstor,"ABT_delta_",perc,"_alpha_",alpha,"_beta_",beta,"_run_",i,sep=""),gzipped=TRUE)


  }


  ## MAXIMIZED LIKELIHOOD RATIO

  for(m in 1:dim(conditions_mLR)[1]) {

    perc <- conditions_mLR[m,1]
    k <- conditions_mLR[m,2]

    if(perc==60)
      mudelta <- es1
    if(perc==75)
      mudelta <- es2    
    if(perc==90)
      mudelta <- es3
    if(perc==95)
      mudelta <- es4

    b1stats <- c(cope)
    sb1stats <- c(sqrt(varcope))

    teller <- c()
    noemer <- c()

    for(l in 1:length(b1stats)) {
  
      if (b1stats[l] > mudelta) {teller[l] <- dnorm(b1stats[l], mean=b1stats[l], sd=sb1stats[l])}
      if (b1stats[l] <= mudelta) {teller[l] <- dnorm(b1stats[l], mudelta, sd=sb1stats[l])}
      if (b1stats[l] < mudelta) {noemer[l]<-dnorm(b1stats[l], mean=b1stats[l], sd=sb1stats[l])}
      if (b1stats[l] >= mudelta){noemer[l]<-dnorm(b1stats[l], mudelta, sd=sb1stats[l])}

    }

    mLR <-teller/noemer
    mLR <- array(mLR, dim=dim(tmap))
    active.mLR <- ifelse(mLR >= k , 1, 0)

    writeNIfTI(active.mLR, paste(threshstor,"mLR_delta_",perc,"_k_",k,"_run_",i,sep=""),gzipped=TRUE)
  }

}





