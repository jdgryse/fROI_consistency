
## PRELIMINARY


#Local working directories
scripts <- "PATH/TO/SCRIPTS/"
copestor <- "PATH/TO/COPES/"
varcopestor <- "PATH/TO/VARCOPES/"
maskstor <- "PATH/TO/MASKS/"
resultstor <- "PATH/TO/STORE/RESULT/ARRAYS/"

#Libraries
library(neuRosim)
library(oro.nifti)
library(AnalyzeFMRI)
library(fmri)


#Source files
conditions_mLR <- read.csv(paste(scripts,"conditions_mLR_hcp.txt",sep=""))
conditions_abt <- read.csv(paste(scripts,"conditions_abt_hcp.txt",sep=""))
source(paste(scripts,"functions_hcp.R",sep=""))


#Parameters
qval <- 0.05
nsub <- 50
nscan <- 176
tau <- 0.5


#Hybrid mask right FFA
rightffa <- readNIfTI(paste(maskstor, "FFA_hybrid_right.nii.gz",sep=""))
rightffa <- ifelse(rightffa==0, NA, 1) #3700 voxels


#Storage matrices for euclidean distance between peaks as defined under the three methods within the general rFFA fROI
peakresults.null <- array(NA, dim=c(nsub,1))
colnames(peakresults.null) <- c("euclid alpha=0.001")
peakresults.fdr <- array(NA, dim=c(nsub,1))
colnames(peakresults.fdr) <- c("euclid qval=0.05")
peakresults.abt60 <- array(NA, dim=c(nsub,6))
colnames(peakresults.abt60) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
peakresults.abt75 <- array(NA, dim=c(nsub,6))
colnames(peakresults.abt75) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
peakresults.abt85 <- array(NA, dim=c(nsub,6))
colnames(peakresults.abt85) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
peakresults.abt95 <- array(NA, dim=c(nsub,6))
colnames(peakresults.abt95) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
peakresults.mLR60 <- array(NA, dim=c(nsub,4))
colnames(peakresults.mLR60) <- c("0.88","1.5","3.68","8")
peakresults.mLR75 <- array(NA, dim=c(nsub,4))
colnames(peakresults.mLR75) <- c("0.88","1.5","3.68","8")
peakresults.mLR85 <- array(NA, dim=c(nsub,4))
colnames(peakresults.mLR85) <- c("0.88","1.5","3.68","8")
peakresults.mLR95 <- array(NA, dim=c(nsub,4))
colnames(peakresults.mLR95) <- c("0.88","1.5","3.68","8")

#Storage matrices for Maitra similarity index for the general rFFA fROI
maitraresults.null <- array(NA, dim=c(nsub,1))
colnames(maitraresults.null) <- c("euclid alpha=0.001")
maitraresults.fdr <- array(NA, dim=c(nsub,1))
colnames(maitraresults.fdr) <- c("euclid qval=0.05")
maitraresults.abt60 <- array(NA, dim=c(nsub,6))
colnames(maitraresults.abt60) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
maitraresults.abt75 <- array(NA, dim=c(nsub,6))
colnames(maitraresults.abt75) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
maitraresults.abt85 <- array(NA, dim=c(nsub,6))
colnames(maitraresults.abt85) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
maitraresults.abt95 <- array(NA, dim=c(nsub,6))
colnames(maitraresults.abt95) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
maitraresults.mLR60 <- array(NA, dim=c(nsub,4))
colnames(maitraresults.mLR60) <- c("0.88","1.5","3.68","8")
maitraresults.mLR75 <- array(NA, dim=c(nsub,4))
colnames(maitraresults.mLR75) <- c("0.88","1.5","3.68","8")
maitraresults.mLR85 <- array(NA, dim=c(nsub,4))
colnames(maitraresults.mLR85) <- c("0.88","1.5","3.68","8")
maitraresults.mLR95 <- array(NA, dim=c(nsub,4))
colnames(maitraresults.mLR95) <- c("0.88","1.5","3.68","8")

#Storage matrices for Maitra similarity index for the summary rFFA fROI
summmaitraresults.null <- array(NA, dim=c(nsub,1))
colnames(summmaitraresults.null) <- c("euclid alpha=0.001")
summmaitraresults.fdr <- array(NA, dim=c(nsub,1))
colnames(summmaitraresults.fdr) <- c("euclid qval=0.05")
summmaitraresults.abt60 <- array(NA, dim=c(nsub,6))
colnames(summmaitraresults.abt60) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
summmaitraresults.abt75 <- array(NA, dim=c(nsub,6))
colnames(summmaitraresults.abt75) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
summmaitraresults.abt85 <- array(NA, dim=c(nsub,6))
colnames(summmaitraresults.abt85) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
summmaitraresults.abt95 <- array(NA, dim=c(nsub,6))
colnames(summmaitraresults.abt95) <- c("0.05_0.1","0.05_0.2","0.05_0.3","0.001_0.1","0.001_0.2","0.001_0.3")
summmaitraresults.mLR60 <- array(NA, dim=c(nsub,4))
colnames(summmaitraresults.mLR60) <- c("0.88","1.5","3.68","8")
summmaitraresults.mLR75 <- array(NA, dim=c(nsub,4))
colnames(summmaitraresults.mLR75) <- c("0.88","1.5","3.68","8")
summmaitraresults.mLR85 <- array(NA, dim=c(nsub,4))
colnames(summmaitraresults.mLR85) <- c("0.88","1.5","3.68","8")
summmaitraresults.mLR95 <- array(NA, dim=c(nsub,4))
colnames(summmaitraresults.mLR95) <- c("0.88","1.5","3.68","8")



## ANALYSIS OF 10 COPES AND VARCOPES + CALCULATING CONSISTENCY MEASURES

for(i in 1:10) {
  
  print(i)
  
  ## READING IN B1 and SB1 files
  cope_run1 <- readNIfTI(paste(copestor, "cope_",i,"_run_1.nii.gz",sep=""))
  cope_run2 <- readNIfTI(paste(copestor, "cope_",i,"_run_2.nii.gz",sep=""))
  varcope_run1 <- readNIfTI(paste(varcopestor, "varcope_",i,"_run_1",sep=""))
  varcope_run2 <- readNIfTI(paste(varcopestor, "varcope_",i,"_run_2.nii.gz",sep=""))
  tmap1 <- (cope_run1/sqrt(varcope_run1))
  tmap2 <- (cope_run2/sqrt(varcope_run2))

  
  ## NHST uncorrected 
  alpha <- 0.001
  pmap.null1 <- pt(tmap1, nscan-2,lower.tail=FALSE)
  active.null1 <- ifelse(pmap.null1 <= alpha, 1, 0)
  
  pmap.null2 <- pt(tmap2, nscan-2,lower.tail=FALSE)
  active.null2 <- ifelse(pmap.null2 <= alpha, 1, 0)

  #Create summary fROIs
  nullsumm1 <- contiguousroinull(cope_run1,sqrt(varcope_run1),active.null1,rightffa)
  nullsumm2 <- contiguousroinull(cope_run2,sqrt(varcope_run2),active.null2,rightffa)
  
  #Calculate consistency measures using functions in hcp_functions.R
  peakresults.null[i,1] <- eucliddistnull(active.null1,active.null2,tmap1,tmap2,rightffa)[1]
  maitraresults.null[i,1] <- maitra(active.null1,active.null2,rightffa)[1]
  summmaitraresults.null[i,1] <- maitra(nullsumm1,nullsumm2,rightffa)[1]
  



  ## NHST FDR CORRECTED
  
  #Run 1
  pvalms1 <- sort(pmap.null1)
  orderpvalms1 <- rank(pmap.null1)
  FDRqval1 <- (orderpvalms1/(length(pmap.null1)))*qval
  pr1 <- ifelse(pvalms1[orderpvalms1] < FDRqval1,1,0)
  pthres.FDR1<-ifelse(sum(pr1,na.rm=TRUE)==0,0,max(FDRqval1[pr1==1],na.rm=TRUE))
  active.fdr1 <- ifelse(pmap.null1 <= pthres.FDR1, 1, 0)
  
  #Run 2
  pvalms2 <- sort(pmap.null2)
  orderpvalms2 <- rank(pmap.null2)
  FDRqval2 <- (orderpvalms2/(length(pmap.null2)))*qval
  pr2 <- ifelse(pvalms2[orderpvalms2] < FDRqval2,1,0)
  pthres.FDR2<-ifelse(sum(pr2,na.rm=TRUE)==0,0,max(FDRqval2[pr2==1],na.rm=TRUE))
  active.fdr2 <- ifelse(pmap.null2 <= pthres.FDR2, 1, 0)
  
  #Create summary rFFA fROIs
  fdrsumm1 <- contiguousroinull(cope_run1,sqrt(varcope_run1),active.fdr1,rightffa)
  fdrsumm2 <- contiguousroinull(cope_run2,sqrt(varcope_run2),active.fdr2,rightffa)

  #Calculate consistency measures
  peakresults.fdr[i,1] <- eucliddistnull(active.fdr1,active.fdr2,tmap1,tmap2,rightffa)[1]
  maitraresults.fdr[i,1] <- maitra(active.fdr1,active.fdr2,rightffa)[1]
  summmaitraresults.fdr[i,1] <- maitra(fdrsumm1,fdrsumm2,rightffa)[1]

  

  ## ALTERNATIVE-BASED THRESHOLDING

  for(j in 1:dim(conditions_abt)[1]) {

    print(j)
    perc <- conditions_abt[j,1]
    alpha <- conditions_abt[j,2]
    beta <- conditions_abt[j,3]
  
    #Delta estimated on other run in hybrid mask
    rightdelta_run1 <- quantile(ifelse(cope_run2*rightffa==0, NA, cope_run2*rightffa),perc/100,na.rm=TRUE)
    rightdelta_run2 <- quantile(ifelse(cope_run1*rightffa==0, NA, cope_run1*rightffa),perc/100,na.rm=TRUE)

    rightpmap.alt1 <- pnorm(tmap1, (rightdelta_run1/sqrt(varcope_run1)), sqrt(((varcope_run1 + tau^2)/varcope_run1)))
    rightactive.abt1 <- ifelse(pmap.null1 <= alpha & rightpmap.alt1 >= beta, 1, 0)
  
    rightpmap.alt2 <- pnorm(tmap2, (rightdelta_run2/sqrt(varcope_run2)), sqrt(((varcope_run2 + tau^2)/varcope_run2)))
    rightactive.abt2 <- ifelse(pmap.null2 <= alpha & rightpmap.alt2 >= beta, 1, 0)

    #Create summary rFFA fROI
    abtsumm1 <- contiguousroiabt(cope_run1, sqrt(varcope_run1), rightdelta_run1, rightactive.abt1, rightffa)
    abtsumm2 <- contiguousroiabt(cope_run2, sqrt(varcope_run2), rightdelta_run2, rightactive.abt2, rightffa)

    #Calculate the consistency measures
    if(perc==60) {
      
      if(alpha == 0.05 && beta == 0.1) {
         peakresults.abt60[i,1] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
         maitraresults.abt60[i,1] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
         summmaitraresults.abt60[i,1] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.05 && beta == 0.2) {
         peakresults.abt60[i,2] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
         maitraresults.abt60[i,2] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
         summmaitraresults.abt60[i,2] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.05 && beta == 0.3) {
         peakresults.abt60[i,3] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
         maitraresults.abt60[i,3] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
         summmaitraresults.abt60[i,3] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.001 && beta == 0.1){
         peakresults.abt60[i,4] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
         maitraresults.abt60[i,4] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
         summmaitraresults.abt60[i,4] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.001 && beta == 0.2){
         peakresults.abt60[i,5] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
         maitraresults.abt60[i,5] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
         summmaitraresults.abt60[i,5] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.001 && beta == 0.3){
         peakresults.abt60[i,6] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
         maitraresults.abt60[i,6] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
         summmaitraresults.abt60[i,6] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
    }


    if(perc==75) {
      
      if(alpha == 0.05 && beta == 0.1) {
        peakresults.abt75[i,1] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt75[i,1] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt75[i,1] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.05 && beta == 0.2) {
        peakresults.abt75[i,2] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt75[i,2] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt75[i,2] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.05 && beta == 0.3) {
        peakresults.abt75[i,3] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt75[i,3] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt75[i,3] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.001 && beta == 0.1){
        peakresults.abt75[i,4] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt75[i,4] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt75[i,4] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.001 && beta == 0.2){
        peakresults.abt75[i,5] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt75[i,5] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt75[i,5] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.001 && beta == 0.3){
        peakresults.abt75[i,6] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt75[i,6] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt75[i,6] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
    }

    if(perc==85) {
      
      if(alpha == 0.05 && beta == 0.1) {
        peakresults.abt85[i,1] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt85[i,1] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt85[i,1] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.05 && beta == 0.2) {
        peakresults.abt85[i,2] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt85[i,2] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt85[i,2] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.05 && beta == 0.3) {
        peakresults.abt85[i,3] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt85[i,3] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt85[i,3] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.001 && beta == 0.1){
        peakresults.abt85[i,4] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt85[i,4] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt85[i,4] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.001 && beta == 0.2){
        peakresults.abt85[i,5] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt85[i,5] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt85[i,5] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.001 && beta == 0.3){
        peakresults.abt85[i,6] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt85[i,6] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt85[i,6] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
    }

    if(perc==95) {
      
      if(alpha == 0.05 && beta == 0.1) {
        peakresults.abt95[i,1] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt95[i,1] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt95[i,1] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.05 && beta == 0.2) {
        peakresults.abt95[i,2] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt95[i,2] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt95[i,2] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.05 && beta == 0.3) {
        peakresults.abt95[i,3] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt95[i,3] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt95[i,3] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.001 && beta == 0.1){
        peakresults.abt95[i,4] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt95[i,4] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt95[i,4] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.001 && beta == 0.2){
        peakresults.abt95[i,5] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt95[i,5] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt95[i,5] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
      if(alpha == 0.001 && beta == 0.3){
        peakresults.abt95[i,6] <- eucliddistabt(rightactive.abt1, rightactive.abt2, cope_run1, cope_run2, varcope_run1, varcope_run2, rightdelta_run1, rightdelta_run2,rightffa)
        maitraresults.abt95[i,6] <- maitra(rightactive.abt1, rightactive.abt2,rightffa)[1]
        summmaitraresults.abt95[i,6] <- maitra(abtsumm1, abtsumm2,rightffa)[1]
      }
    }
  }


  ## MAXIMIZED LIKELIHOOD RATIO

  for(m in 1:dim(conditions_mLR)[1]) {

    print(m)

    perc <- conditions_mLR[m,1]
    k <- conditions_mLR[m,2]
  
    #Delta estimated on other run in hybrid mask
    rightdelta_run1 <- quantile(ifelse(cope_run2*rightffa==0, NA, cope_run2*rightffa),perc/100,na.rm=TRUE)
    rightdelta_run2 <- quantile(ifelse(cope_run1*rightffa==0, NA, cope_run1*rightffa),perc/100,na.rm=TRUE)
  
    #mLR is computed in each voxel for run 1
    b1stats1 <- c(cope_run1)
    sb1stats1 <- c(sqrt(varcope_run1))
  
    rightteller1 <- c()
    rightnoemer1 <- c()
  
    for(l in 1:length(b1stats1)) {
    
      if (b1stats1[l] > rightdelta_run1) {rightteller1[l] <- dnorm(b1stats1[l], mean=b1stats1[l], sd=sb1stats1[l])}
      if (b1stats1[l] <= rightdelta_run1) {rightteller1[l] <- dnorm(b1stats1[l], rightdelta_run1, sd=sb1stats1[l])}
      if (b1stats1[l] < rightdelta_run1) {rightnoemer1[l]<-dnorm(b1stats1[l], mean=b1stats1[l], sd=sb1stats1[l])}
      if (b1stats1[l] >= rightdelta_run1){rightnoemer1[l]<-dnorm(b1stats1[l], rightdelta_run1, sd=sb1stats1[l])}
    
    }
  
    rightmLR1 <-rightteller1/rightnoemer1
    rightmLR1 <- array(rightmLR1, dim=dim(tmap1))
    rightactive.mLR1 <- ifelse(rightmLR1 >= k , 1, 0)

    #mLR is computed in each voxel for run 2  
    b1stats2 <- c(cope_run2)
    sb1stats2 <- c(sqrt(varcope_run2))
  
    rightteller2 <- c()
    rightnoemer2 <- c()
  
    for(l in 1:length(b1stats2)) {
    
      if (b1stats2[l] > rightdelta_run2) {rightteller2[l] <- dnorm(b1stats2[l], mean=b1stats2[l], sd=sb1stats2[l])}
      if (b1stats2[l] <= rightdelta_run2) {rightteller2[l] <- dnorm(b1stats2[l], rightdelta_run2, sd=sb1stats2[l])}
      if (b1stats2[l] < rightdelta_run2) {rightnoemer2[l]<-dnorm(b1stats2[l], mean=b1stats2[l], sd=sb1stats2[l])}
      if (b1stats2[l] >= rightdelta_run2){rightnoemer2[l]<-dnorm(b1stats2[l], rightdelta_run2, sd=sb1stats2[l])}
    
    } 
  
    rightmLR2 <-rightteller2/rightnoemer2
    rightmLR2 <- array(rightmLR2, dim=dim(tmap2))
    rightactive.mLR2 <- ifelse(rightmLR2 >= k , 1, 0)

    #Create summary rFFA fROI
    mLRsumm1 <- contiguousroimLR(rightteller1, rightmLR1, rightactive.mLR1,rightffa)
    mLRsumm2 <- contiguousroimLR(rightteller2, rightmLR2, rightactive.mLR2,rightffa)

    #Calculate the consistency measures
    if(perc==60 && k==0.88) {
      peakresults.mLR60[i,1] <- eucliddistmLR(rightactive.mLR1,rightactive.mLR2,rightmLR1,rightmLR2,rightteller1,rightteller2,rightffa)
      maitraresults.mLR60[i,1] <- maitra(rightactive.mLR1,rightactive.mLR2,rightffa)[1]
      summmaitraresults.mLR60[i,1] <- maitra(mLRsumm1,mLRsumm2,rightffa)[1]
    }
    if(perc==60 && k==1.5) {
      peakresults.mLR60[i,2] <- eucliddistmLR(rightactive.mLR1,rightactive.mLR2,rightmLR1,rightmLR2,rightteller1,rightteller2,rightffa)
      maitraresults.mLR60[i,2] <- maitra(rightactive.mLR1,rightactive.mLR2,rightffa)[1]
      summmaitraresults.mLR60[i,2] <- maitra(mLRsumm1,mLRsumm2,rightffa)[1]
    }
    if(perc==60 && k==3.68) {
      peakresults.mLR60[i,3] <- eucliddistmLR(rightactive.mLR1,rightactive.mLR2,rightmLR1,rightmLR2,rightteller1,rightteller2,rightffa)
      maitraresults.mLR60[i,3] <- maitra(rightactive.mLR1,rightactive.mLR2,rightffa)[1]
      summmaitraresults.mLR60[i,3] <- maitra(mLRsumm1,mLRsumm2,rightffa)[1]
    }
    if(perc==60 && k==8) {
      peakresults.mLR60[i,4] <- eucliddistmLR(rightactive.mLR1,rightactive.mLR2,rightmLR1,rightmLR2,rightteller1,rightteller2,rightffa)
      maitraresults.mLR60[i,4] <- maitra(rightactive.mLR1,rightactive.mLR2,rightffa)[1]
      summmaitraresults.mLR60[i,4] <- maitra(mLRsumm1,mLRsumm2,rightffa)[1]
    }


    if(perc==75 && k==0.88) {
      peakresults.mLR75[i,1] <- eucliddistmLR(rightactive.mLR1,rightactive.mLR2,rightmLR1,rightmLR2,rightteller1,rightteller2,rightffa)
      maitraresults.mLR75[i,1] <- maitra(rightactive.mLR1,rightactive.mLR2,rightffa)[1]
      summmaitraresults.mLR75[i,1] <- maitra(mLRsumm1,mLRsumm2,rightffa)[1]
    }
    if(perc==75 && k==1.5) {
      peakresults.mLR75[i,2] <- eucliddistmLR(rightactive.mLR1,rightactive.mLR2,rightmLR1,rightmLR2,rightteller1,rightteller2,rightffa)
      maitraresults.mLR75[i,2] <- maitra(rightactive.mLR1,rightactive.mLR2,rightffa)[1]
      summmaitraresults.mLR75[i,2] <- maitra(mLRsumm1,mLRsumm2,rightffa)[1]
    }
    if(perc==75 && k==3.68) {
      peakresults.mLR75[i,3] <- eucliddistmLR(rightactive.mLR1,rightactive.mLR2,rightmLR1,rightmLR2,rightteller1,rightteller2,rightffa)
      maitraresults.mLR75[i,3] <- maitra(rightactive.mLR1,rightactive.mLR2,rightffa)[1]
      summmaitraresults.mLR75[i,3] <- maitra(mLRsumm1,mLRsumm2,rightffa)[1]
    }
    if(perc==75 && k==8) {
      peakresults.mLR75[i,4] <- eucliddistmLR(rightactive.mLR1,rightactive.mLR2,rightmLR1,rightmLR2,rightteller1,rightteller2,rightffa)
      maitraresults.mLR75[i,4] <- maitra(rightactive.mLR1,rightactive.mLR2,rightffa)[1]
      summmaitraresults.mLR75[i,4] <- maitra(mLRsumm1,mLRsumm2,rightffa)[1]
    }


    if(perc==85 && k==0.88) {
      peakresults.mLR85[i,1] <- eucliddistmLR(rightactive.mLR1,rightactive.mLR2,rightmLR1,rightmLR2,rightteller1,rightteller2,rightffa)
      maitraresults.mLR85[i,1] <- maitra(rightactive.mLR1,rightactive.mLR2,rightffa)[1]
      summmaitraresults.mLR85[i,1] <- maitra(mLRsumm1,mLRsumm2,rightffa)[1]
    }
    if(perc==85 && k==1.5) {
      peakresults.mLR85[i,2] <- eucliddistmLR(rightactive.mLR1,rightactive.mLR2,rightmLR1,rightmLR2,rightteller1,rightteller2,rightffa)
      maitraresults.mLR85[i,2] <- maitra(rightactive.mLR1,rightactive.mLR2,rightffa)[1]
      summmaitraresults.mLR85[i,2] <- maitra(mLRsumm1,mLRsumm2,rightffa)[1]
    }
    if(perc==85 && k==3.68) {
      peakresults.mLR85[i,3] <- eucliddistmLR(rightactive.mLR1,rightactive.mLR2,rightmLR1,rightmLR2,rightteller1,rightteller2,rightffa)
      maitraresults.mLR85[i,3] <- maitra(rightactive.mLR1,rightactive.mLR2,rightffa)[1]
      summmaitraresults.mLR85[i,3] <- maitra(mLRsumm1,mLRsumm2,rightffa)[1]
    }
    if(perc==85 && k==8) {
      peakresults.mLR85[i,4] <- eucliddistmLR(rightactive.mLR1,rightactive.mLR2,rightmLR1,rightmLR2,rightteller1,rightteller2,rightffa)
      maitraresults.mLR85[i,4] <- maitra(rightactive.mLR1,rightactive.mLR2,rightffa)[1]
      summmaitraresults.mLR85[i,4] <- maitra(mLRsumm1,mLRsumm2,rightffa)[1]
    }


    if(perc==95 && k==0.88) {
      peakresults.mLR95[i,1] <- eucliddistmLR(rightactive.mLR1,rightactive.mLR2,rightmLR1,rightmLR2,rightteller1,rightteller2,rightffa)
      maitraresults.mLR95[i,1] <- maitra(rightactive.mLR1,rightactive.mLR2,rightffa)[1]
      summmaitraresults.mLR95[i,1] <- maitra(mLRsumm1,mLRsumm2,rightffa)[1]
    }
    if(perc==95 && k==1.5) {
      peakresults.mLR95[i,2] <- eucliddistmLR(rightactive.mLR1,rightactive.mLR2,rightmLR1,rightmLR2,rightteller1,rightteller2,rightffa)
      maitraresults.mLR95[i,2] <- maitra(rightactive.mLR1,rightactive.mLR2,rightffa)[1]
      summmaitraresults.mLR95[i,2] <- maitra(mLRsumm1,mLRsumm2,rightffa)[1]
    }
    if(perc==95 && k==3.68) {
      peakresults.mLR95[i,3] <- eucliddistmLR(rightactive.mLR1,rightactive.mLR2,rightmLR1,rightmLR2,rightteller1,rightteller2,rightffa)
      maitraresults.mLR95[i,3] <- maitra(rightactive.mLR1,rightactive.mLR2,rightffa)[1]
      summmaitraresults.mLR95[i,3] <- maitra(mLRsumm1,mLRsumm2,rightffa)[1]
    }
    if(perc==95 && k==8) {
      peakresults.mLR95[i,4] <- eucliddistmLR(rightactive.mLR1,rightactive.mLR2,rightmLR1,rightmLR2,rightteller1,rightteller2,rightffa)
      maitraresults.mLR95[i,4] <- maitra(rightactive.mLR1,rightactive.mLR2,rightffa)[1]
      summmaitraresults.mLR95[i,4] <- maitra(mLRsumm1,mLRsumm2,rightffa)[1]
    }
  }
}



## WRITE AWAY RESULT ARRAYS CONSISTENCY MEASURES


setwd(resultstor)
write.csv(peakresults.null, "NULL_peak_perc",row.names=FALSE)
write.csv(peakresults.fdr, "FDR_peak_perc",row.names=FALSE)
write.csv(peakresults.abt60, "ABT_peak_perc_60",row.names=FALSE)
write.csv(peakresults.abt75, "ABT_peak_perc_75",row.names=FALSE)
write.csv(peakresults.abt85, "ABT_peak_perc_85",row.names=FALSE)
write.csv(peakresults.abt95, "ABT_peak_perc_95",row.names=FALSE)
write.csv(peakresults.mLR60, "mLR_peak_perc_60",row.names=FALSE)
write.csv(peakresults.mLR75, "mLR_peak_perc_75",row.names=FALSE)
write.csv(peakresults.mLR85, "mLR_peak_perc_85",row.names=FALSE)
write.csv(peakresults.mLR95, "mLR_peak_perc_95",row.names=FALSE)

write.csv(maitraresults.null, "NULL_maitra_perc",row.names=FALSE)
write.csv(maitraresults.fdr, "FDR_maitra_perc",row.names=FALSE)
write.csv(maitraresults.abt60, "ABT_maitra_perc_60",row.names=FALSE)
write.csv(maitraresults.abt75, "ABT_maitra_perc_75",row.names=FALSE)
write.csv(maitraresults.abt85, "ABT_maitra_perc_85",row.names=FALSE)
write.csv(maitraresults.abt95, "ABT_maitra_perc_95",row.names=FALSE)
write.csv(maitraresults.mLR60, "mLR_maitra_perc_60",row.names=FALSE)
write.csv(maitraresults.mLR75, "mLR_maitra_perc_75",row.names=FALSE)
write.csv(maitraresults.mLR85, "mLR_maitra_perc_85",row.names=FALSE)
write.csv(maitraresults.mLR95, "mLR_maitra_perc_95",row.names=FALSE)

write.csv(summmaitraresults.null, "NULL_summmaitra_perc",row.names=FALSE)
write.csv(summmaitraresults.fdr, "FDR_summmaitra_perc",row.names=FALSE)
write.csv(summmaitraresults.abt60, "ABT_summmaitra_perc_60",row.names=FALSE)
write.csv(summmaitraresults.abt75, "ABT_summmaitra_perc_75",row.names=FALSE)
write.csv(summmaitraresults.abt85, "ABT_summmaitra_perc_85",row.names=FALSE)
write.csv(summmaitraresults.abt95, "ABT_summmaitra_perc_95",row.names=FALSE)
write.csv(summmaitraresults.mLR60, "mLR_summmaitra_perc_60",row.names=FALSE)
write.csv(summmaitraresults.mLR75, "mLR_summmaitra_perc_75",row.names=FALSE)
write.csv(summmaitraresults.mLR85, "mLR_summmaitra_perc_85",row.names=FALSE)
write.csv(summmaitraresults.mLR95, "mLR_summmaitra_perc_95",row.names=FALSE)




