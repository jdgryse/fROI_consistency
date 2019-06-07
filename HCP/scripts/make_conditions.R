##########################################
## PARAMETER CONDITIONS HCP
##########################################

###MAXIMIZED LR

percval <- c(60,75,85,95)
kval <- c(0.88,1.5,3.68,8)

perc <- rep(percval, each=length(kval))
k <- rep(kval, times=length(percval))

conditions <- as.data.frame(cbind(perc,k))
colnames(conditions) <- c("perc","k")

write.csv(conditions,file="PATH/TO/SCRIPTS/conditions_mLR_hcp.txt",row.names=FALSE)


###ABT

percval <- c(60,75,85,95)
alphaval <- c(0.05,0.001)
betaval <- c(0.1,0.2,0.3)

perc <- rep(percval, each=length(alphaval)*length(betaval))
alpha <- rep(alphaval, times=length(betaval)*length(percval))
beta <- rep(betaval,times=length(alphaval)*length(percval))

conditions <- as.data.frame(cbind(perc,alpha,beta))
colnames(conditions) <- c("perc","alpha","beta")

write.csv(conditions,file="PATH/TO/SCRIPTS/conditions_abt_hcp.txt",row.names=FALSE)
