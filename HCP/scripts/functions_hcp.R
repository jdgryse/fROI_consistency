#Function for maitra similarity index

maitra <- function(image1, image2, mask) {

  #only voxels that overlap with mask
  image1 <- ifelse(image1*mask==1, 1, NA)
  image2 <- ifelse(image2*mask==1, 1, NA)
  
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


#Function for euclidian distance between peaks with nhst

eucliddistnull <- function(thresh1, thresh2,t1, t2, mask) {

  #only voxels that overlap with hybrid mask
  thresh1 <- ifelse(thresh1*mask==1, 1, NA)
  thresh2 <- ifelse(thresh2*mask==1, 1, NA)

  tnull1 <- t1 * thresh1
  tnull2 <- t2 * thresh2

  #Coordinates of maximum t-value under null
  coord1 <- which(tnull1==max(tnull1,na.rm=TRUE), arr.ind=TRUE)
  coord2 <- which(tnull2==max(tnull2,na.rm=TRUE), arr.ind=TRUE)

  eucdist <- sqrt((coord1[1]-coord2[1])^2+(coord1[2]-coord2[2])^2+(coord1[3]-coord2[3])^2)

  return(eucdist)

}

#Function for euclidian distance between peaks with ABT

eucliddistabt <- function(thresh1, thresh2,c1, c2, v1, v2, es1, es2, mask) {

  #Find voxels overlapping with hybrid mask
  thresh1 <- ifelse(thresh1*mask==1, 1, NA)
  thresh2 <- ifelse(thresh2*mask==1, 1, NA)

  #Compute t-values under the alternative
  t1 <- (c1-es1)/sqrt(v1)
  t2 <- (c1-es2)/sqrt(v2)

  tabt1 <- t1 * thresh1
  tabt2 <- t2 * thresh2

  #Coordinates of maximum t-value under null
  coord1 <- which(tabt1==max(tabt1,na.rm=TRUE), arr.ind=TRUE)
  coord2 <- which(tabt2==max(tabt2,na.rm=TRUE), arr.ind=TRUE)

  eucdist <- sqrt((coord1[1]-coord2[1])^2+(coord1[2]-coord2[2])^2+(coord1[3]-coord2[3])^2)

  return(eucdist)

}

#Function for euclidian distance between peaks with mLR

eucliddistmLR <- function(thresh1, thresh2,meas1, meas2, tel1, tel2, mask) {


  #Find voxels overlapping with hybrid mask
  thresh1 <- ifelse(thresh1*mask==1, 1, NA)
  thresh2 <- ifelse(thresh2*mask==1, 1, NA)
  tel1 <- tel1*mask
  tel2 <- tel2*mask
  meas1 <- meas1*mask
  meas2 <- meas2*mask

  #Coordinates of voxel with maximum mLR (sometimes mLR is infinite and multiple infinite values are detected, in that case use numerator of mLR)
  if(max(meas1,na.rm=TRUE)==Inf) {
    tel1 <- ifelse(meas1==Inf, tel1, NA)
    coord1 <- which(tel1==max(tel1,na.rm=TRUE), arr.ind=TRUE)
  }
  if(max(meas1,na.rm=TRUE)!=Inf)
    coord1 <- which(meas1==max(meas1,na.rm=TRUE), arr.ind=TRUE) 
  if(max(meas2,na.rm=TRUE)==Inf) {
    tel2 <- ifelse(meas2==Inf, tel2, NA)
    coord2 <- which(tel2==max(tel2,na.rm=TRUE), arr.ind=TRUE)
  }
   if(max(meas2,na.rm=TRUE)!=Inf)
    coord2 <- which(meas2==max(meas2,na.rm=TRUE), arr.ind=TRUE) 

  eucdist <- sqrt((coord1[1]-coord2[1])^2+(coord1[2]-coord2[2])^2+(coord1[3]-coord2[3])^2)

  return(eucdist)

}


#Function to detect the top 50 voxels contiguous to the peak for summary fROI with NHST

contiguousroinull <- function(b1, sb1, thresh,mask) {

  tmap <- b1/sb1
  tmap <- ifelse(thresh*mask==1, tmap, 0)

  #Find the peak activated voxel in roi and its contiguous voxels and their t-values
  coord <- c(which(tmap==max(tmap, na.rm=TRUE), arr.ind=TRUE))
  
  x <- coord[1]
  y <- coord[2]
  z <- coord[3]

    surroundings <- c(tmap[x+1,y+1,z-1],
                    tmap[x,y+1,z-1],
                    tmap[x-1,y+1,z-1],
                    tmap[x+1,y,z-1],
                    tmap[x,y,z-1],
                    tmap[x-1,y,z-1],
                    tmap[x+1,y-1,z-1],
                    tmap[x,y-1,z-1],
                    tmap[x-1,y-1,z-1],

          tmap[x+1,y+1,z],
                    tmap[x,y+1,z],
                    tmap[x-1,y+1,z],
                    tmap[x+1,y,z],
                    tmap[x-1,y,z],
                    tmap[x+1,y-1,z],
                    tmap[x,y-1,z],
                    tmap[x-1,y-1,z],

          tmap[x+1,y+1,z+1],
                    tmap[x,y+1,z+1],
                    tmap[x-1,y+1,z+1],
                    tmap[x+1,y,z+1],
                    tmap[x,y,z+1],
                    tmap[x-1,y,z+1],
                    tmap[x+1,y-1,z+1],
                    tmap[x,y-1,z+1],
                    tmap[x-1,y-1,z+1])                    


    #Sort the surrounding t-values
    sort_surr <- sort(surroundings[which(surroundings!=0)], decreasing=TRUE)


    if(length(sort_surr)==0) {

      selectmap50 <- array(0, dim=dim(tmap))
      selectmap50[x,y,z] <- 1

    } else {



    contig1 <- array(NA, dim=c(3,length(surroundings)))

    for (j in 1:length(sort_surr)) {

      contig1[,j] <- c(which(tmap==sort_surr[j], arr.ind=TRUE))

    }

    #Make arrays to check whether the second, third etc, round already contains voxels from previous rounds
    tvalues1 <- c(max(tmap, na.rm = TRUE),sort_surr)
    tvalues2 <- c(max(tmap, na.rm = TRUE),sort_surr)
    tvalues3 <- c()

    arrsurr <- array(NA, dim=c(26,length(sort_surr)))

    #Find the surrounding voxels of the voxels surrounding the peak in the roi
    for(j in 1:length(sort_surr)) {

      coord2 <- c(which(tmap==sort_surr[j], arr.ind=TRUE))

      x <- coord2[1]
      y <- coord2[2]
      z <- coord2[3]

      #2nd round to get 50 top activated voxels around peak
      surroundings2 <- c(tmap[x+1,y+1,z-1],
                      tmap[x,y+1,z-1],
                      tmap[x-1,y+1,z-1],
                      tmap[x+1,y,z-1],
                      tmap[x,y,z-1],
                      tmap[x-1,y,z-1],
                      tmap[x+1,y-1,z-1],
                      tmap[x,y-1,z-1],
                      tmap[x-1,y-1,z-1],

            tmap[x+1,y+1,z],
                      tmap[x,y+1,z],
                      tmap[x-1,y+1,z],
                      tmap[x+1,y,z],
                      tmap[x-1,y,z],
                      tmap[x+1,y-1,z],
                      tmap[x,y-1,z],
                      tmap[x-1,y-1,z],

            tmap[x+1,y+1,z+1],
                      tmap[x,y+1,z+1],
                      tmap[x-1,y+1,z+1],
                      tmap[x+1,y,z+1],
                      tmap[x,y,z+1],
                      tmap[x-1,y,z+1],
                      tmap[x+1,y-1,z+1],
                      tmap[x,y-1,z+1],
                      tmap[x-1,y-1,z+1])    

      arrsurr[,j] <- surroundings2

      }

    #Get the coordinates of the t-values surrounding the peak
    alltvals <- sort(unique(c(tvalues1,c(arrsurr[which(arrsurr!=0)]))),decreasing=TRUE)

    if(length(alltvals) >= 50) {

      alltvals <- alltvals[1:50]
      allcoords <- array(NA, dim=c(3,length(alltvals)))

      #Find the coordinates of these t-values
      for (j in 1:length(alltvals)) {

        allcoords[,j] <- c(which(tmap==alltvals[j], arr.ind=TRUE))

      }

    }

    if(length(alltvals) < 50) {

      allmLRvals <- alltvals[1:length(alltvals)]
      allcoords <- array(NA, dim=c(3,length(alltvals)))

      #Find the coordinates of these t-values
      for (j in 1:length(alltvals)) {

        allcoords[,j] <- c(which(tmap==alltvals[j], arr.ind=TRUE))

      }


    }

    #Put all voxels that we want on 1, rest on 0
    selectmap50 <- array(0, dim=dim(tmap))

    #Create mask for b1 map depending on how many voxels we will compute the mean with
    for (j in 1:length(alltvals)) {

      selectmap50[allcoords[1,j], allcoords[2,j], allcoords[3,j]] <- 1

    }
  }
    return(selectmap50)

}


#Function to detect the top 50 voxels contiguous to the peak for summary fROI with ABT

contiguousroiabt <- function(b1, sb1, delta, thresh, mask) {

  tmap <- (b1-delta)/sb1
  tmap <- ifelse(thresh*mask==1, tmap, 0)

  #Find the peak activated voxel in roi and its contiguous voxels and their t-values
  coord <- c(which(tmap==max(tmap, na.rm=TRUE), arr.ind=TRUE))
  
  x <- coord[1]
  y <- coord[2]
  z <- coord[3]

    surroundings <- c(tmap[x+1,y+1,z-1],
                    tmap[x,y+1,z-1],
                    tmap[x-1,y+1,z-1],
                    tmap[x+1,y,z-1],
                    tmap[x,y,z-1],
                    tmap[x-1,y,z-1],
                    tmap[x+1,y-1,z-1],
                    tmap[x,y-1,z-1],
                    tmap[x-1,y-1,z-1],

          tmap[x+1,y+1,z],
                    tmap[x,y+1,z],
                    tmap[x-1,y+1,z],
                    tmap[x+1,y,z],
                    tmap[x-1,y,z],
                    tmap[x+1,y-1,z],
                    tmap[x,y-1,z],
                    tmap[x-1,y-1,z],

          tmap[x+1,y+1,z+1],
                    tmap[x,y+1,z+1],
                    tmap[x-1,y+1,z+1],
                    tmap[x+1,y,z+1],
                    tmap[x,y,z+1],
                    tmap[x-1,y,z+1],
                    tmap[x+1,y-1,z+1],
                    tmap[x,y-1,z+1],
                    tmap[x-1,y-1,z+1])                    


    #Sort the surrounding t-values
    sort_surr <- sort(surroundings[which(surroundings!=0)], decreasing=TRUE)


    if(length(sort_surr)==0) {

      selectmap50 <- array(0, dim=dim(tmap))
      selectmap50[x,y,z] <- 1

    } else {


    contig1 <- array(NA, dim=c(3,length(sort_surr)))

    for (j in 1:length(sort_surr)) {

      contig1[,j] <- c(which(tmap==sort_surr[j], arr.ind=TRUE))

    }

    #Make arrays to check whether the second, third etc, round already contains voxels from previous rounds
    tvalues1 <- c(max(tmap, na.rm = TRUE),sort_surr)
    tvalues2 <- c(max(tmap, na.rm = TRUE),sort_surr)
    tvalues3 <- c()

    arrsurr <- array(NA, dim=c(26,length(sort_surr)))

    #Find the surrounding voxels of the voxels surrounding the peak in the roi
    for(j in 1:length(sort_surr)) {

      coord2 <- c(which(tmap==sort_surr[j], arr.ind=TRUE))

      x <- coord2[1]
      y <- coord2[2]
      z <- coord2[3]

      #2nd round to get 50 top activated voxels around peak
      surroundings2 <- c(tmap[x+1,y+1,z-1],
                      tmap[x,y+1,z-1],
                      tmap[x-1,y+1,z-1],
                      tmap[x+1,y,z-1],
                      tmap[x,y,z-1],
                      tmap[x-1,y,z-1],
                      tmap[x+1,y-1,z-1],
                      tmap[x,y-1,z-1],
                      tmap[x-1,y-1,z-1],

            tmap[x+1,y+1,z],
                      tmap[x,y+1,z],
                      tmap[x-1,y+1,z],
                      tmap[x+1,y,z],
                      tmap[x-1,y,z],
                      tmap[x+1,y-1,z],
                      tmap[x,y-1,z],
                      tmap[x-1,y-1,z],

            tmap[x+1,y+1,z+1],
                      tmap[x,y+1,z+1],
                      tmap[x-1,y+1,z+1],
                      tmap[x+1,y,z+1],
                      tmap[x,y,z+1],
                      tmap[x-1,y,z+1],
                      tmap[x+1,y-1,z+1],
                      tmap[x,y-1,z+1],
                      tmap[x-1,y-1,z+1])    

      arrsurr[,j] <- surroundings2

      }

    #Get the coordinates of the t-values surrounding the peak
    alltvals <- sort(unique(c(tvalues1,c(arrsurr[which(arrsurr!=0)]))),decreasing=TRUE)

    if(length(alltvals) >= 50) {

      alltvals <- alltvals[1:50]
      allcoords <- array(NA, dim=c(3,length(alltvals)))

      #Find the coordinates of these t-values
      for (j in 1:length(alltvals)) {

        allcoords[,j] <- c(which(tmap==alltvals[j], arr.ind=TRUE))

      }

    }

    if(length(alltvals) < 50) {

      allmLRvals <- alltvals[1:length(alltvals)]
      allcoords <- array(NA, dim=c(3,length(alltvals)))

      #Find the coordinates of these t-values
      for (j in 1:length(alltvals)) {

        allcoords[,j] <- c(which(tmap==alltvals[j], arr.ind=TRUE))

      }


    }

    #Put all voxels that we want on 1, rest on 0
    selectmap50 <- array(0, dim=dim(tmap))

    for (j in 1:length(alltvals)) {

      selectmap50[allcoords[1,j], allcoords[2,j], allcoords[3,j]] <- 1

    }
  }
    return(selectmap50)

}


#Function to detect the top 50 voxels contiguous to the peak for summary fROI with mLR

contiguousroimLR <- function(tel, meas, thresh, mask) {

  meas <- ifelse(thresh*mask == 1, meas, 0)
  tel <- ifelse(thresh*mask == 1, tel, 0)

  #Find the peak activated voxel in roi and its contiguous voxels and their mLR values
  if(max(meas,na.rm=TRUE)==Inf) {
    tel <- ifelse(meas==Inf, tel, NA)
    coord <- which(tel==max(tel,na.rm=TRUE), arr.ind=TRUE)
    maxmLR <- max(tel,na.rm=TRUE)
  }
  if(max(meas,na.rm=TRUE)!=Inf) {
    coord <- which(meas==max(meas,na.rm=TRUE), arr.ind=TRUE) 
    maxmLR <- max(meas,na.rm=TRUE)
  }

  x <- coord[1]
  y <- coord[2]
  z <- coord[3]

  surroundings <- c(meas[x+1,y+1,z-1],
                    meas[x,y+1,z-1],
                    meas[x-1,y+1,z-1],
                    meas[x+1,y,z-1],
                    meas[x,y,z-1],
                    meas[x-1,y,z-1],
                    meas[x+1,y-1,z-1],
                    meas[x,y-1,z-1],
                    meas[x-1,y-1,z-1],
                    
                    meas[x+1,y+1,z],
                    meas[x,y+1,z],
                    meas[x-1,y+1,z],
                    meas[x+1,y,z],
                    meas[x-1,y,z],
                    meas[x+1,y-1,z],
                    meas[x,y-1,z],
                    meas[x-1,y-1,z],
                    
                    meas[x+1,y+1,z+1],
                    meas[x,y+1,z+1],
                    meas[x-1,y+1,z+1],
                    meas[x+1,y,z+1],
                    meas[x,y,z+1],
                    meas[x-1,y,z+1],
                    meas[x+1,y-1,z+1],
                    meas[x,y-1,z+1],
                    meas[x-1,y-1,z+1])                    
                 


    #Sort the surrounding t-values
    sort_surr <- sort(surroundings[which(surroundings!=0)], decreasing=TRUE)

    if(length(sort_surr)==0) {

      selectmap50 <- array(0, dim=dim(meas))
      selectmap50[x,y,z] <- 1

    } else {

    contig1 <- array(NA, dim=c(3,length(sort_surr)))

    for (j in 1:length(sort_surr)) {

      contig1[,j] <- c(which(meas==sort_surr[j], arr.ind=TRUE))

    }

    #Make arrays to check whether the second, third etc, round already contains voxels from previous rounds
    mLRvalues1 <- c(maxmLR,sort_surr)
    mLRvalues2 <- c(maxmLR,sort_surr)

    arrsurr <- array(NA, dim=c(26,length(sort_surr)))

    #Find the surrounding voxels of the voxels surrounding the peak in the roi
    for(j in 1:length(sort_surr)) {

      coord2 <- c(which(meas==sort_surr[j], arr.ind=TRUE))

      x <- coord2[1]
      y <- coord2[2]
      z <- coord2[3]

      #2nd round to get 50 top activated voxels around peak
      surroundings2 <- c(meas[x+1,y+1,z-1],
                        meas[x,y+1,z-1],
                        meas[x-1,y+1,z-1],
                        meas[x+1,y,z-1],
                        meas[x,y,z-1],
                        meas[x-1,y,z-1],
                        meas[x+1,y-1,z-1],
                        meas[x,y-1,z-1],
                        meas[x-1,y-1,z-1],
                        
                        meas[x+1,y+1,z],
                        meas[x,y+1,z],
                        meas[x-1,y+1,z],
                        meas[x+1,y,z],
                        meas[x-1,y,z],
                        meas[x+1,y-1,z],
                        meas[x,y-1,z],
                        meas[x-1,y-1,z],
                        
                        meas[x+1,y+1,z+1],
                        meas[x,y+1,z+1],
                        meas[x-1,y+1,z+1],
                        meas[x+1,y,z+1],
                        meas[x,y,z+1],
                        meas[x-1,y,z+1],
                        meas[x+1,y-1,z+1],
                        meas[x,y-1,z+1],
                        meas[x-1,y-1,z+1])        

      arrsurr[,j] <- surroundings2

      }

    #Get coordinates of contiguous t-values
    allmLRvals <- sort(unique(c(mLRvalues1,c(arrsurr[which(arrsurr!=0)]))),decreasing=TRUE)

    if(length(allmLRvals) >= 50) {

      allmLRvals <- allmLRvals[1:50]
      allcoords <- array(NA, dim=c(3,length(allmLRvals)))

      #Find the coordinates of these t-values
      for (j in 1:length(allmLRvals)) {

        allcoords[,j] <- c(which(meas==allmLRvals[j], arr.ind=TRUE))

      }

    }

    if(length(allmLRvals) < 50) {

      allmLRvals <- allmLRvals[1:length(allmLRvals)]
      allcoords <- array(NA, dim=c(3,length(allmLRvals)))

      #Find the coordinates of these t-values
      for (j in 1:length(allmLRvals)) {

        allcoords[,j] <- c(which(meas==allmLRvals[j], arr.ind=TRUE))

      }


    }

    #Put all voxels that we want on 1, rest on 0
    selectmap50 <- array(0, dim=dim(meas))

    for (j in 1:length(allmLRvals)) {

      selectmap50[allcoords[1,j], allcoords[2,j], allcoords[3,j]] <- 1

    }
  }

    return(selectmap50)

}
