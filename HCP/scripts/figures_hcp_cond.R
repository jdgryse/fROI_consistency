
#Local working directories
scripts <- "PATH/TO/SCRIPTS/"
results <- "PATH/TO/RESULTS"

## PRELIMINARY

#Libraries
library(neuRosim)
library(oro.nifti)
library(AnalyzeFMRI)
library(lattice)
library(ggplot2)
library(ggjoy)
library(gplots)
library(ggpubr)
library(ggridges)

## Read in result arrays with consistency measures created with analysis_hcp.R
null_maitra <- read.csv(paste(results, "NULL_maitra_perc",sep=""))
fdr_maitra <- read.csv(paste(results, "FDR_maitra_perc",sep=""))
abt60_maitra <- read.csv(paste(results,"ABT_maitra_perc_60",sep=""))
abt75_maitra <- read.csv(paste(results,"ABT_maitra_perc_75",sep=""))
abt85_maitra <- read.csv(paste(results,"ABT_maitra_perc_85",sep=""))
abt95_maitra <- read.csv(paste(results, "ABT_maitra_perc_95",sep=""))
mLR60_maitra <- read.csv(paste(results,"mLR_maitra_perc_60",sep=""))
mLR75_maitra <- read.csv(paste(results,"mLR_maitra_perc_75",sep=""))
mLR85_maitra <- read.csv(paste(results,"mLR_maitra_perc_85",sep=""))
mLR95_maitra <- read.csv(paste(results,"mLR_maitra_perc_95",sep=""))

#If you want to create the Figure with the Maitra similarity index for the summary rFFA fROI, please comment the arrays above and uncomment these below
#null_maitra <- read.csv(paste(results, "NULL_summmaitra_perc",sep=""))
#fdr_maitra <- read.csv(paste(results, "FDR_summmaitra_perc",sep=""))
#abt60_maitra <- read.csv(paste(results,"ABT_summmaitra_perc_60",sep=""))
#abt75_maitra <- read.csv(paste(results,"ABT_summmaitra_perc_75",sep=""))
#abt85_maitra <- read.csv(paste(results,"ABT_summmaitra_perc_85",sep=""))
#abt95_maitra <- read.csv(paste(results, "ABT_summmaitra_perc_95",sep=""))
#mLR60_maitra <- read.csv(paste(results,"mLR_summmaitra_perc_60",sep=""))
#mLR75_maitra <- read.csv(paste(results,"mLR_summmaitra_perc_75",sep=""))
#mLR85_maitra <- read.csv(paste(results,"mLR_summmaitra_perc_85",sep=""))
#mLR95_maitra <- read.csv(paste(results,"mLR_summmaitra_perc_95",sep=""))

null_peak <- read.csv(paste(results, "NULL_peak_perc",sep=""))*2
fdr_peak <- read.csv(paste(results, "FDR_peak_perc",sep=""))*2
abt60_peak <- read.csv(paste(results,"ABT_peak_perc_60",sep=""))*2
abt75_peak <- read.csv(paste(results,"ABT_peak_perc_75",sep=""))*2
abt85_peak <- read.csv(paste(results,"ABT_peak_perc_85",sep=""))*2
abt95_peak <- read.csv(paste(results, "ABT_peak_perc_95",sep=""))*2
mLR60_peak <- read.csv(paste(results,"mLR_peak_perc_60",sep=""))*2
mLR75_peak <- read.csv(paste(results,"mLR_peak_perc_75",sep=""))*2
mLR85_peak <- read.csv(paste(results,"mLR_peak_perc_85",sep=""))*2
mLR95_peak <- read.csv(paste(results,"mLR_peak_perc_95",sep=""))*2



##############################
## VIOLIN PLOT MAITRA       ##
##############################

#Calculating largest MSI value, so y-axis is comparable for all figures
maxy <- max(null_maitra, fdr_maitra, abt60_maitra,abt75_maitra,abt85_maitra,abt95_maitra,mLR60_maitra,mLR75_maitra,mLR85_maitra,mLR95_maitra)


## ALTERNATIVE-BASED THRESHOLDING

#NHST
nhst_indic <- c(rep("Uncorr.",times=length(null_maitra[,1])),rep("FDR",times=length(fdr_maitra[,1])))
nhst_measure <- c(null_maitra[,1],fdr_maitra[,1])
nhst.maitra<- data.frame(matrix(NA, nrow = length(nhst_measure), ncol = 2))
nhst.maitra[,1] <- nhst_measure
nhst.maitra[,2] <- nhst_indic
colnames(nhst.maitra) <- c("Maitra similarity index","Method")

null <- ggplot(nhst.maitra, aes(x=nhst.maitra[,2], y=nhst.maitra[,1])) +
  geom_violin(aes(color=nhst.maitra[,2]),trim=TRUE) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  xlab("Method") + ylab("Maitra similarity index")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("NHST") + theme(legend.position = "none")


#Violin plot ABT60 alpha = 0.05 Maitra
abt60_005_indic <- c(rep("0.1",times=length(abt60_maitra[,1])),rep("0.2",times=length(abt60_maitra[,2])),rep("0.3",times=length(abt60_maitra[,3])))
abt60_005_measure <- c(abt60_maitra[,1],abt60_maitra[,2],abt60_maitra[,3])
abt60_005.maitra <- data.frame(matrix(NA, nrow = length(abt60_005_measure), ncol = 2))
abt60_005.maitra[,1] <- abt60_005_measure
abt60_005.maitra[,2] <- abt60_005_indic
colnames(abt60_005.maitra) <- c("Maitra similarity index",expression(beta))

abt60_005 <- ggplot(abt60_005.maitra, aes(x=abt60_005.maitra[,2], y=abt60_005.maitra[,1])) +
  geom_violin(aes(color=abt60_005.maitra[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab(expression(beta)) + ylab("Maitra similarity index")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("ABT + 60th") + theme(legend.position = "none")

#Violin plot ABT75 alpha = 0.05 Maitra
abt75_005_indic <- c(rep("0.1",times=length(abt75_maitra[,1])),rep("0.2",times=length(abt75_maitra[,2])),rep("0.3",times=length(abt75_maitra[,3])))
abt75_005_measure <- c(abt75_maitra[,1],abt75_maitra[,2],abt75_maitra[,3])
abt75_005.maitra <- data.frame(matrix(NA, nrow = length(abt75_005_measure), ncol = 2))
abt75_005.maitra[,1] <- abt75_005_measure
abt75_005.maitra[,2] <- abt75_005_indic
colnames(abt75_005.maitra) <- c("Maitra similarity index",expression(beta))

abt75_005 <- ggplot(abt75_005.maitra, aes(x=abt75_005.maitra[,2], y=abt75_005.maitra[,1])) +
  geom_violin(aes(color=abt75_005.maitra[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab(expression(beta)) + ylab("Maitra similarity index")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("ABT + 75th") + theme(legend.position = "none")

#Violin plot ABT90 alpha = 0.05 Maitra
abt85_005_indic <- c(rep("0.1",times=length(abt85_maitra[,1])),rep("0.2",times=length(abt85_maitra[,2])),rep("0.3",times=length(abt85_maitra[,3])))
abt85_005_measure <- c(abt85_maitra[,1],abt85_maitra[,2],abt85_maitra[,3])
abt85_005.maitra <- data.frame(matrix(NA, nrow = length(abt85_005_measure), ncol = 2))
abt85_005.maitra[,1] <- abt85_005_measure
abt85_005.maitra[,2] <- abt85_005_indic
colnames(abt90_005.maitra) <- c("Maitra similarity index",expression(beta))

abt85_005 <- ggplot(abt85_005.maitra, aes(x=abt85_005.maitra[,2], y=abt85_005.maitra[,1])) +
  geom_violin(aes(color=abt85_005.maitra[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab(expression(beta)) + ylab("Maitra similarity index")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("ABT + 85th") + theme(legend.position = "none")

#Violin plot ABT95 alpha = 0.05 Maitra
abt95_005_indic <- c(rep("0.1",times=length(abt95_maitra[,1])),rep("0.2",times=length(abt95_maitra[,2])),rep("0.3",times=length(abt95_maitra[,3])))
abt95_005_measure <- c(abt95_maitra[,1],abt95_maitra[,2],abt95_maitra[,3])
abt95_005.maitra <- data.frame(matrix(NA, nrow = length(abt95_005_measure), ncol = 2))
abt95_005.maitra[,1] <- abt95_005_measure
abt95_005.maitra[,2] <- abt95_005_indic
colnames(abt95_005.maitra) <- c("Maitra similarity index",expression(beta))

abt95_005 <- ggplot(abt95_005.maitra, aes(x=abt95_005.maitra[,2], y=abt95_005.maitra[,1])) +
  geom_violin(aes(color=abt95_005.maitra[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab(expression(beta)) + ylab("Maitra similarity index")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("ABT + 95th") + theme(legend.position = "none")

#Violin plot ABT60 alpha = 0.001 Maitra
abt60_0001_indic <- c(rep("0.1",times=length(abt60_maitra[,1])),rep("0.2",times=length(abt60_maitra[,2])),rep("0.3",times=length(abt60_maitra[,3])))
abt60_0001_measure <- c(abt60_maitra[,4],abt60_maitra[,5],abt60_maitra[,6])
abt60_0001.maitra <- data.frame(matrix(NA, nrow = length(abt60_0001_measure), ncol = 2))
abt60_0001.maitra[,1] <- abt60_0001_measure
abt60_0001.maitra[,2] <- abt60_0001_indic
colnames(abt60_0001.maitra) <- c("Maitra similarity index",expression(beta))

abt60_0001 <- ggplot(abt60_0001.maitra, aes(x=abt60_0001.maitra[,2], y=abt60_0001.maitra[,1])) +
  geom_violin(aes(color=abt60_0001.maitra[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab(expression(beta)) + ylab("Maitra similarity index")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("ABT + 60th") + theme(legend.position = "none")

#Violin plot ABT75 alpha = 0.001 Maitra
abt75_0001_indic <- c(rep("0.1",times=length(abt75_maitra[,1])),rep("0.2",times=length(abt75_maitra[,2])),rep("0.3",times=length(abt75_maitra[,3])))
abt75_0001_measure <- c(abt75_maitra[,4],abt75_maitra[,5],abt75_maitra[,6])
abt75_0001.maitra <- data.frame(matrix(NA, nrow = length(abt75_0001_measure), ncol = 2))
abt75_0001.maitra[,1] <- abt75_0001_measure
abt75_0001.maitra[,2] <- abt75_0001_indic
colnames(abt75_0001.maitra) <- c("Maitra similarity index",expression(beta))

abt75_0001 <- ggplot(abt75_0001.maitra, aes(x=abt75_0001.maitra[,2], y=abt75_0001.maitra[,1])) +
  geom_violin(aes(color=abt75_0001.maitra[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab(expression(beta)) + ylab("Maitra similarity index")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("ABT + 75th") + theme(legend.position = "none")

#Violin plot ABT90 alpha = 0.001 Maitra
abt85_0001_indic <- c(rep("0.1",times=length(abt85_maitra[,1])),rep("0.2",times=length(abt85_maitra[,2])),rep("0.3",times=length(abt85_maitra[,3])))
abt85_0001_measure <- c(abt85_maitra[,4],abt85_maitra[,5],abt85_maitra[,6])
abt85_0001.maitra <- data.frame(matrix(NA, nrow = length(abt85_0001_measure), ncol = 2))
abt85_0001.maitra[,1] <- abt85_0001_measure
abt85_0001.maitra[,2] <- abt85_0001_indic
colnames(abt90_0001.maitra) <- c("Maitra similarity index",expression(beta))

abt85_0001 <- ggplot(abt85_0001.maitra, aes(x=abt85_0001.maitra[,2], y=abt85_0001.maitra[,1])) +
  geom_violin(aes(color=abt85_0001.maitra[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab(expression(beta)) + ylab("Maitra similarity index")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("ABT + 85th") + theme(legend.position = "none")

#Violin plot ABT95 alpha = 0.001 Maitra
abt95_0001_indic <- c(rep("0.1",times=length(abt95_maitra[,1])),rep("0.2",times=length(abt95_maitra[,2])),rep("0.3",times=length(abt95_maitra[,3])))
abt95_0001_measure <- c(abt95_maitra[,4],abt95_maitra[,5],abt95_maitra[,6])
abt95_0001.maitra <- data.frame(matrix(NA, nrow = length(abt95_0001_measure), ncol = 2))
abt95_0001.maitra[,1] <- abt95_0001_measure
abt95_0001.maitra[,2] <- abt95_0001_indic
colnames(abt95_0001.maitra) <- c("Maitra similarity index",expression(beta))

abt95_0001 <- ggplot(abt95_0001.maitra, aes(x=abt95_0001.maitra[,2], y=abt95_0001.maitra[,1])) +
  geom_violin(aes(color=abt95_0001.maitra[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab(expression(beta)) + ylab("Maitra similarity index")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("ABT + 95th") + theme(legend.position = "none")

#Combine violin plots for ABT with alpha = 0.05
figure_abt_005 <- ggarrange(null, abt60_005, abt75_005, abt85_005,abt95_005,ncol=5,nrow=1)
figure_abt_005

#Combine violin plots for ABT with alpha = 0.001
figure_abt_0001 <- ggarrange(null, abt60_0001, abt75_0001, abt85_0001,abt95_0001,ncol=5,nrow=1)
figure_abt_0001

## MAXIMIZED LIKELIHOOD RATIO

#mLR + 60
mLR60_indic <- c(rep("0.88",times=length(mLR60_maitra[,1])),rep("1.5",times=length(mLR60_maitra[,2])),rep("3.68",times=length(mLR60_maitra[,3])),rep("8",times=length(mLR60_maitra[,4])))
mLR60_measure <- c(mLR60_maitra[,1],mLR60_maitra[,2],mLR60_maitra[,3],mLR60_maitra[,4])
mLR60.maitra <- data.frame(matrix(NA, nrow = length(mLR60_measure), ncol = 2))
mLR60.maitra[,1] <- mLR60_measure
mLR60.maitra[,2] <- mLR60_indic
colnames(mLR60.maitra) <- c("Maitra similarity index","k")

mLR60<- ggplot(mLR60.maitra, aes(x=mLR60.maitra[,2], y=mLR60.maitra[,1])) +
  geom_violin(aes(color=mLR60.maitra[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab("k") + ylab("Maitra similarity index")  + ylim(0,maxy) + geom_boxplot(width=0.1)  + theme_bw() + ggtitle("mLR + 60th") + theme(legend.position = "none")


#mLR + 75
mLR75_indic <- c(rep("0.88",times=length(mLR75_maitra[,1])),rep("1.5",times=length(mLR75_maitra[,2])),rep("3.68",times=length(mLR75_maitra[,3])),rep("8",times=length(mLR75_maitra[,4])))
mLR75_measure <- c(mLR75_maitra[,1],mLR75_maitra[,2],mLR75_maitra[,3],mLR75_maitra[,4])
mLR75.maitra <- data.frame(matrix(NA, nrow = length(mLR75_measure), ncol = 2))
mLR75.maitra[,1] <- mLR75_measure
mLR75.maitra[,2] <- mLR75_indic
colnames(mLR75.maitra) <- c("Maitra similarity index","k")

mLR75<- ggplot(mLR75.maitra, aes(x=mLR75.maitra[,2], y=mLR75.maitra[,1])) +
  geom_violin(aes(color=mLR75.maitra[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab("k") + ylab("Maitra similarity index")  + ylim(0,maxy) + geom_boxplot(width=0.1)  + theme_bw() + ggtitle("mLR + 75th") + theme(legend.position = "none")

#mLR + 90
mLR85_indic <- c(rep("0.88",times=length(mLR85_maitra[,1])),rep("1.5",times=length(mLR85_maitra[,2])),rep("3.68",times=length(mLR85_maitra[,3])),rep("8",times=length(mLR85_maitra[,4])))
mLR85_measure <- c(mLR85_maitra[,1],mLR85_maitra[,2],mLR85_maitra[,3],mLR85_maitra[,4])
mLR85.maitra <- data.frame(matrix(NA, nrow = length(mLR85_measure), ncol = 2))
mLR85.maitra[,1] <- mLR85_measure
mLR85.maitra[,2] <- mLR85_indic
colnames(mLR85.maitra) <- c("Maitra similarity index","k")

mLR85<- ggplot(mLR85.maitra, aes(x=mLR85.maitra[,2], y=mLR85.maitra[,1])) +
  geom_violin(aes(color=mLR85.maitra[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab("k") + ylab("Maitra similarity index")  + ylim(0,maxy) + geom_boxplot(width=0.1)  + theme_bw() + ggtitle("mLR + 85th") + theme(legend.position = "none")

#mLR + 95
mLR95_indic <- c(rep("0.88",times=length(mLR95_maitra[,1])),rep("1.5",times=length(mLR95_maitra[,2])),rep("3.68",times=length(mLR95_maitra[,3])),rep("8",times=length(mLR95_maitra[,4])))
mLR95_measure <- c(mLR95_maitra[,1],mLR95_maitra[,2],mLR95_maitra[,3],mLR95_maitra[,4])
mLR95.maitra <- data.frame(matrix(NA, nrow = length(mLR95_measure), ncol = 2))
mLR95.maitra[,1] <- mLR95_measure
mLR95.maitra[,2] <- mLR95_indic
colnames(mLR95.maitra) <- c("Maitra similarity index","k")

mLR95<- ggplot(mLR95.maitra, aes(x=mLR95.maitra[,2], y=mLR95.maitra[,1])) +
  geom_violin(aes(color=mLR95.maitra[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab("k") + ylab("Maitra similarity index")  + ylim(0,maxy) + geom_boxplot(width=0.1)  + theme_bw() + ggtitle("mLR + 95th") + theme(legend.position = "none")

#Combine mLR plots
figure_mLR <- ggarrange(null, mLR60, mLR75, mLR85,mLR95,ncol=5,nrow=1)
figure_mLR




#################################################
## VIOLIN PLOT EUCLIDEAN DISTANCE PEAK         ##
#################################################

## ALTERNATIVE-BASED THRESHOLDING

maxy <- 50

#NHST
nhst_indic <- c(rep("Uncorr.",times=length(null_peak[,1])),rep("FDR",times=length(fdr_peak[,1])))
nhst_measure <- c(null_peak[,1],fdr_peak[,1])
nhst.peak<- data.frame(matrix(NA, nrow = length(nhst_measure), ncol = 2))
nhst.peak[,1] <- nhst_measure
nhst.peak[,2] <- nhst_indic
colnames(nhst.peak) <- c("Euclidian distance (voxels)","Method")

null <- ggplot(nhst.peak, aes(x=nhst.peak[,2], y=nhst.peak[,1])) +
  geom_violin(aes(color=nhst.peak[,2]),trim=TRUE) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  xlab("Method") + ylab("Euclidean distance (mm)")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("NHST") + theme(legend.position = "none")


#Violin plot ABT60 alpha = 0.05 peak
abt60_005_indic <- c(rep("0.1",times=length(abt60_peak[,1])),rep("0.2",times=length(abt60_peak[,2])),rep("0.3",times=length(abt60_peak[,3])))
abt60_005_measure <- c(abt60_peak[,1],abt60_peak[,2],abt60_peak[,3])
abt60_005.peak <- data.frame(matrix(NA, nrow = length(abt60_005_measure), ncol = 2))
abt60_005.peak[,1] <- abt60_005_measure
abt60_005.peak[,2] <- abt60_005_indic
colnames(abt60_005.peak) <- c("Euclidian distance (voxels)","beta")

abt60_005 <- ggplot(abt60_005.peak, aes(x=abt60_005.peak[,2], y=abt60_005.peak[,1])) +
  geom_violin(aes(color=abt60_005.peak[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab(expression(beta)) + ylab("Euclidean distance (mm)")  + ylim(0,maxy) + geom_boxplot(width=0.1)  + theme_bw() + ggtitle("ABT + 60th")+ theme(legend.position = "none")

#Violin plot ABT75 alpha = 0.05 peak
abt75_005_indic <- c(rep("0.1",times=length(abt75_peak[,1])),rep("0.2",times=length(abt75_peak[,2])),rep("0.3",times=length(abt75_peak[,3])))
abt75_005_measure <- c(abt75_peak[,1],abt75_peak[,2],abt75_peak[,3])
abt75_005.peak <- data.frame(matrix(NA, nrow = length(abt75_005_measure), ncol = 2))
abt75_005.peak[,1] <- abt75_005_measure
abt75_005.peak[,2] <- abt75_005_indic
colnames(abt75_005.peak) <- c("Euclidian distance (voxels)","beta")

abt75_005 <- ggplot(abt75_005.peak, aes(x=abt75_005.peak[,2], y=abt75_005.peak[,1])) +
  geom_violin(aes(color=abt75_005.peak[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab(expression(beta)) + ylab("Euclidean distance (mm)")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("ABT + 75th") + theme(legend.position = "none")

#Violin plot ABT85 alpha = 0.05 peak
abt85_005_indic <- c(rep("0.1",times=length(abt85_peak[,1])),rep("0.2",times=length(abt85_peak[,2])),rep("0.3",times=length(abt85_peak[,3])))
abt85_005_measure <- c(abt85_peak[,1],abt85_peak[,2],abt85_peak[,3])
abt85_005.peak <- data.frame(matrix(NA, nrow = length(abt85_005_measure), ncol = 2))
abt85_005.peak[,1] <- abt85_005_measure
abt85_005.peak[,2] <- abt85_005_indic
colnames(abt85_005.peak) <- c("Euclidian distance (voxels)","beta")

abt85_005 <- ggplot(abt85_005.peak, aes(x=abt85_005.peak[,2], y=abt85_005.peak[,1])) +
  geom_violin(aes(color=abt85_005.peak[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab(expression(beta)) + ylab("Euclidean distance (mm)")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("ABT + 85th") + theme(legend.position = "none")

#Violin plot ABT95 alpha = 0.05 peak
abt95_005_indic <- c(rep("0.1",times=length(abt95_peak[,1])),rep("0.2",times=length(abt95_peak[,2])),rep("0.3",times=length(abt95_peak[,3])))
abt95_005_measure <- c(abt95_peak[,1],abt95_peak[,2],abt95_peak[,3])
abt95_005.peak <- data.frame(matrix(NA, nrow = length(abt95_005_measure), ncol = 2))
abt95_005.peak[,1] <- abt95_005_measure
abt95_005.peak[,2] <- abt95_005_indic
colnames(abt95_005.peak) <- c("Euclidian distance (voxels)","beta")

abt95_005 <- ggplot(abt95_005.peak, aes(x=abt95_005.peak[,2], y=abt95_005.peak[,1])) +
  geom_violin(aes(color=abt95_005.peak[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab(expression(beta)) + ylab("Euclidean distance (mm)")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("ABT + 95th") + theme(legend.position = "none")

#Violin plot ABT60 alpha = 0.001 peak
abt60_0001_indic <- c(rep("0.1",times=length(abt60_peak[,1])),rep("0.2",times=length(abt60_peak[,2])),rep("0.3",times=length(abt60_peak[,3])))
abt60_0001_measure <- c(abt60_peak[,4],abt60_peak[,5],abt60_peak[,6])
abt60_0001.peak <- data.frame(matrix(NA, nrow = length(abt60_0001_measure), ncol = 2))
abt60_0001.peak[,1] <- abt60_0001_measure
abt60_0001.peak[,2] <- abt60_0001_indic
colnames(abt60_0001.peak) <- c("Euclidian distance (voxels)","beta")

abt60_0001 <- ggplot(abt60_0001.peak, aes(x=abt60_0001.peak[,2], y=abt60_0001.peak[,1])) +
  geom_violin(aes(color=abt60_0001.peak[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab(expression(beta)) + ylab("Euclidean distance (mm)")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("ABT + 60th") + theme(legend.position = "none")

#Violin plot ABT75 alpha = 0.001 peak
abt75_0001_indic <- c(rep("0.1",times=length(abt75_peak[,1])),rep("0.2",times=length(abt75_peak[,2])),rep("0.3",times=length(abt75_peak[,3])))
abt75_0001_measure <- c(abt75_peak[,4],abt75_peak[,5],abt75_peak[,6])
abt75_0001.peak <- data.frame(matrix(NA, nrow = length(abt75_0001_measure), ncol = 2))
abt75_0001.peak[,1] <- abt75_0001_measure
abt75_0001.peak[,2] <- abt75_0001_indic
colnames(abt75_0001.peak) <- c("Euclidian distance (voxels)","beta")

abt75_0001 <- ggplot(abt75_0001.peak, aes(x=abt75_0001.peak[,2], y=abt75_0001.peak[,1])) +
  geom_violin(aes(color=abt75_0001.peak[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab(expression(beta)) + ylab("Euclidean distance (mm)")  + ylim(0,maxy) + geom_boxplot(width=0.1)  + theme_bw() + ggtitle("ABT + 75th") + theme(legend.position = "none")

#Violin plot ABT85 alpha = 0.001 peak
abt85_0001_indic <- c(rep("0.1",times=length(abt85_peak[,1])),rep("0.2",times=length(abt85_peak[,2])),rep("0.3",times=length(abt85_peak[,3])))
abt85_0001_measure <- c(abt85_peak[,4],abt85_peak[,5],abt85_peak[,6])
abt85_0001.peak <- data.frame(matrix(NA, nrow = length(abt85_0001_measure), ncol = 2))
abt85_0001.peak[,1] <- abt85_0001_measure
abt85_0001.peak[,2] <- abt85_0001_indic
colnames(abt85_0001.peak) <- c("Euclidian distance (voxels)","beta")

abt85_0001 <- ggplot(abt85_0001.peak, aes(x=abt85_0001.peak[,2], y=abt85_0001.peak[,1])) +
  geom_violin(aes(color=abt85_0001.peak[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab(expression(beta)) +ylab("Euclidean distance (mm)")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("ABT + 85th") + theme(legend.position = "none")

#Violin plot ABT95 alpha = 0.001 peak
abt95_0001_indic <- c(rep("0.1",times=length(abt95_peak[,1])),rep("0.2",times=length(abt95_peak[,2])),rep("0.3",times=length(abt95_peak[,3])))
abt95_0001_measure <- c(abt95_peak[,4],abt95_peak[,5],abt95_peak[,6])
abt95_0001.peak <- data.frame(matrix(NA, nrow = length(abt95_0001_measure), ncol = 2))
abt95_0001.peak[,1] <- abt95_0001_measure
abt95_0001.peak[,2] <- abt95_0001_indic
colnames(abt95_0001.peak) <- c("Euclidian distance (voxels)","beta")

abt95_0001 <- ggplot(abt95_0001.peak, aes(x=abt95_0001.peak[,2], y=abt95_0001.peak[,1])) +
  geom_violin(aes(color=abt95_0001.peak[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab(expression(beta)) + ylab("Euclidean distance (mm)")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("ABT + 95th") + theme(legend.position = "none")

#Combine plots for ABT with alpha = 0.05
figure_abt_005 <- ggarrange(null, abt60_005, abt75_005, abt85_005,abt95_005,ncol=5,nrow=1)
figure_abt_005

#Combine plots for ABT with alpha = 0.001
figure_abt_0001 <- ggarrange(null, abt60_0001, abt75_0001, abt85_0001,abt95_0001,ncol=5,nrow=1)
figure_abt_0001




## MAXIMIZED LIKELIHOOD RATIO


#mLR + 60
mLR60_indic <- c(rep("0.88",times=length(mLR60_peak[,1])),rep("1.5",times=length(mLR60_peak[,2])),rep("3.68",times=length(mLR60_peak[,3])),rep("8",times=length(mLR60_peak[,4])))
mLR60_measure <- c(mLR60_peak[,1],mLR60_peak[,2],mLR60_peak[,3],mLR60_peak[,4])
mLR60.peak <- data.frame(matrix(NA, nrow = length(mLR60_measure), ncol = 2))
mLR60.peak[,1] <- mLR60_measure
mLR60.peak[,2] <- mLR60_indic
colnames(mLR60.peak) <- c("Euclidean distance (voxels)","k")

mLR60<- ggplot(mLR60.peak, aes(x=mLR60.peak[,2], y=mLR60.peak[,1])) +
  geom_violin(aes(color=mLR60.peak[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab("k") + ylab("Euclidean distance (mm)")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("mLR + 60th") + theme(legend.position = "none")


#mLR + 75
mLR75_indic <- c(rep("0.88",times=length(mLR75_peak[,1])),rep("1.5",times=length(mLR75_peak[,2])),rep("3.68",times=length(mLR75_peak[,3])),rep("8",times=length(mLR75_peak[,4])))
mLR75_measure <- c(mLR75_peak[,1],mLR75_peak[,2],mLR75_peak[,3],mLR75_peak[,4])
mLR75.peak <- data.frame(matrix(NA, nrow = length(mLR75_measure), ncol = 2))
mLR75.peak[,1] <- mLR75_measure
mLR75.peak[,2] <- mLR75_indic
colnames(mLR75.peak) <- c("Euclidian distance (voxels)","k")

mLR75<- ggplot(mLR75.peak, aes(x=mLR75.peak[,2], y=mLR75.peak[,1])) +
  geom_violin(aes(color=mLR75.peak[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab("k") + ylab("Euclidean distance (mm)")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("mLR + 75th") + theme(legend.position = "none")

#mLR + 85
mLR85_indic <- c(rep("0.88",times=length(mLR85_peak[,1])),rep("1.5",times=length(mLR85_peak[,2])),rep("3.68",times=length(mLR85_peak[,3])),rep("8",times=length(mLR85_peak[,4])))
mLR85_measure <- c(mLR85_peak[,1],mLR85_peak[,2],mLR85_peak[,3],mLR85_peak[,4])
mLR85.peak <- data.frame(matrix(NA, nrow = length(mLR85_measure), ncol = 2))
mLR85.peak[,1] <- mLR85_measure
mLR85.peak[,2] <- mLR85_indic
colnames(mLR85.peak) <- c("Euclidian distance (voxels)","k")

mLR85<- ggplot(mLR85.peak, aes(x=mLR85.peak[,2], y=mLR85.peak[,1])) +
  geom_violin(aes(color=mLR85.peak[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab("k") + ylab("Euclidean distance (mm)")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("mLR + 85th") + theme(legend.position = "none")

#mLR + 95
mLR95_indic <- c(rep("0.88",times=length(mLR95_peak[,1])),rep("1.5",times=length(mLR95_peak[,2])),rep("3.68",times=length(mLR95_peak[,3])),rep("8",times=length(mLR95_peak[,4])))
mLR95_measure <- c(mLR95_peak[,1],mLR95_peak[,2],mLR95_peak[,3],mLR95_peak[,4])
mLR95.peak <- data.frame(matrix(NA, nrow = length(mLR95_measure), ncol = 2))
mLR95.peak[,1] <- mLR95_measure
mLR95.peak[,2] <- mLR95_indic
colnames(mLR95.peak) <- c("Euclidian distance (voxels)","k")

mLR95<- ggplot(mLR95.peak, aes(x=mLR95.peak[,2], y=mLR95.peak[,1])) +
  geom_violin(aes(color=mLR95.peak[,2])) + scale_fill_manual(values=c("#E99999","#E69F00", "#56B4E9")) +
  xlab("k") + ylab("Euclidean distance (mm)")  + ylim(0,maxy) + geom_boxplot(width=0.1) + theme_bw() + ggtitle("mLR + 95th") + theme(legend.position = "none")

#Combine all plots with mLR
figure_mLR <- ggarrange(null, mLR60, mLR75, mLR85,mLR95,ncol=5,nrow=1)
figure_mLR


