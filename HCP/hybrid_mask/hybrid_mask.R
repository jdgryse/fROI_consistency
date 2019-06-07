
#Local working directories
masks <- "PATH/TO/MASKS"

#Libraries
library(oro.nifti)

#Anatomical masks of right and left hemisphere
lefthem <- readNIfTI(paste(masks,"left.nii.gz", sep=""))
lefthem <- ifelse(lefthem > 0,1,0)
righthem <- readNIfTI(paste(masks,"right.nii.gz", sep=""))
righthem <- ifelse(righthem > 0,1,0)

#Anatomical mask of bilateral fusiform gyrus
ffa1 <- readNIfTI(paste(masks,"occipital_fusiform.nii.gz", sep=""))
ffa2 <- readNIfTI(paste(masks,"occipital_temporal_fusiform.nii.gz", sep=""))
ffa3 <- readNIfTI(paste(masks,"posterior_fusiform.nii.gz", sep=""))
ffa4 <- readNIfTI(paste(masks,"anterior_fusiform.nii.gz", sep=""))

ffa1 <- ifelse(ffa1 > 0, 1, 0)
ffa2 <- ifelse(ffa2 > 0, 1, 0)
ffa3 <- ifelse(ffa3 > 0, 1, 0)
ffa4 <- ifelse(ffa4 > 0, 1, 0)

combined <- ffa1 + ffa2 + ffa3 + ffa4
bilffa <- ifelse(combined > 0, 1, 0)

#Functional mask bilateral FFA
func <- readNIfTI(paste(masks,"FFA_functional.nii.gz", sep=""))
func <- ifelse(func > 0, 1, 0)

#Creating bilateral hybrid mask
overlap <- bilffa + func
hybrid <- ifelse(overlap==2, 1, 0)

#Creating left hybrid mask (not used in paper)
lefthyb <- hybrid + lefthem
lefthyb <- ifelse(lefthyb==2, 1, 0)

#Creating right hybrid mask (used in paper)
righthyb <- hybrid + righthem
righthyb <- ifelse(righthyb==2, 1, 0)

#Write away masks
writeNIfTI(righthyb, paste(masks, "FFA_hybrid_right",sep=""),gzipped=TRUE)
writeNIfTI(lefthyb, paste(masks, "FFA_hybrid_left", sep=""),gzipped=TRUE)