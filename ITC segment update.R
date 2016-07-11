# Changes made:

# 160607: - Enabled the search window to size to be scalable.
#         - Tidied up some of the code.
#         - Made the maxima search tool work from tallest to smallest.

# 160711: - Changed the distance of pixels to the local maximum a raial distance rather than straight line distance.
#         - Scale the maximum distance according to the allometry.
#           The programme should work through the trees from tallest to shortest and identify the crown area on the fly; 
#           This would over come problems associated with climbing up the canopy on to adjacent crowns.

# !!! take the window scaling as linear for the moment. Later it would be worth converting to an appropriate allometry.
# !!! Make it possible to specify in the function an alternative height maximum to scale the search window;
#     - if not given this will revert to the maximum tree height in the image - Although that probably isn't much good
#     - either way a warning should be given if the tallest point in the image exceeds the max value provided by the user.




#rm(list=ls())
# Distance is for the radius (not diameter) in pixels (not meters):
#imagery = chm.plot.buffer; searchWinSize = 3; TRESHSeed = 0.45; TRESHCrown = 0.55; DIST = 10; specT = 0

#TRESHSeed=0.45 # Proportional height relative to tree maximum.
#TRESHCrown=0.55 # Proportional height relative to tree mean.
#DIST=8 # Distance in pixels from tree maximum.
#specT=2 # Minimum height in m for tree or crown.
#win.min = 3
#win.max = 19
#z.max = 45
#z.min = specT

# Fits a basic linear model to set up the relationship between tree height and window size:
win.lmfun<-function(win.min, win.max, z.min, z.max)
{
  y<-c(win.min, win.max)
  z<-c(z.min, z.max)
  fm.win<-lm(y~z)
  return(fm.win)
}

# Calculates the window size when given the relationship between tree height and window size as a fitted model (fm.win)
# and z the tree heights:
win.sizer<-function(z, fm.win)
{
  win.pred<-predict(fm.win, newdata=data.frame(z=z))
  win.size<-2*round((win.pred+1)/2)-1
  names(win.size)<-NULL
  return(win.size)
}

# Extracts just the locations of the tree maxima:
itIMG_ts<-function (imagery = NULL, win.min, win.max, z.min, z.max, specT = 0) 
{
  # !!! You need to make some appropriate checks of the arguments here !!!
  
  # Blurs the chm:
  imagery <- raster::focal(imagery, w = matrix(1, 3, 3), 
                           fun = function(x) {
                             mean(x, na.rm = T)
                           })
  # Extracts the image data as a matrix:
  Max <- matrix(dim(imagery)[2], dim(imagery)[1], data = imagery[, 
                                                                 ], byrow = FALSE)
  # Flips the image the right way again:
  Max <- Max[1:dim(imagery)[2], dim(imagery)[1]:1]
  Gnew <- Max     # Copies the max matrix.
  Max[, ] <- 0    # Sets max to 0s. This will be used later for storing ...
  Index <- Max    # Copies max again.
  Index[, ] <- 0  # Sets index to 0s.
  Gnew[is.na(Gnew)] <- 0 # Sets any nas to 0s.
  Gnew[Gnew < specT] <- 0 # Any values beneath the minimum height are set to 0s.
  index = 1       # Initiates the index
  II <- which(Gnew != 0, arr.ind = T) # Extracts the locations of pixels which are bigger than the min tree height.
  dim(II)
  
  fm.win<-win.lmfun(win.min, win.max, z.min, z.max) # generates the relationship between height and search window size.
  # Extracts only the pixels that are sufficiently far from the image edge for the search window. In each direction.
  # The search window is selected according to the height of the tree:
  z<-Gnew[II[,1]+nrow(Gnew)*(II[,2]-1)] # extracts the tree heights from the matrix
  WinSize<-win.sizer(z, fm.win) # Finds the window size for each pixel
  half.WinSize<-ceiling(WinSize/2) # Halfs the window size for subsetting by pixels far enough from the image edge
  
  # Extracts only the pixels far enough from the image edge.
  II.ind<-II[, 1] >= half.WinSize & 
    II[, 1] <= (nrow(Gnew) - half.WinSize) &
    II[,2] >= half.WinSize &
    II[,2] <= (ncol(Gnew) - half.WinSize)
  II<-II[II.ind,]
  WinSize<-WinSize[II.ind]
  z<-z[II.ind]
  dim(II)
  
  # reorder from greatest to least according to z value (i.e. big trees first):
  z.order<-order(z, decreasing = TRUE)
  II<-II[z.order,]
  WinSize<-WinSize[z.order]
  
  # Works through each pixel one by one:
  indexII<-1
  #indexII<-(1:nrow(II)[1])[1]
  for (indexII in 1:dim(II)[1]) 
  {
    r = as.numeric(II[indexII, 1]) # Extracts the row pos.
    k = as.numeric(II[indexII, 2]) # Extracts the column pos.
    searchWinSize<-WinSize[indexII]
    hc.sWS<-ceiling(searchWinSize/2) # half the search window size rounded to floor
    hf.sWS<-floor(searchWinSize/2) # half the search window size rounded to ceiling
    FIL <- matrix(searchWinSize, searchWinSize, data = NA) 
    # Extracts the search window for the pixel:
    FIL <- Gnew[(r - hf.sWS):(r + hf.sWS), 
                (k - hf.sWS):(k + hf.sWS)]
    # Extracts the window from Max indicating whether a tree has already being designated or not:
    Max.chk<-Max[(r - hf.sWS):(r + hf.sWS),
                 (k - hf.sWS):(k + hf.sWS)]
    
    # If the focal pixel has the greatest value in the window & there is no tree already assigned in the output matrix within the window & the max value is no 0...
    # because the order is from tallest to shortest, large trees will always suppress the designation of small trees.
    if (FIL[hc.sWS, hc.sWS] == max(FIL, na.rm = T) &
        max(Max.chk, na.rm = T) == 0 &
        max(FIL, na.rm = T) != 0) 
    {
      Max[r, k] <- 1 # A logical assignment of tallest tree
      Index[r, k] <- index # Assigns the sequence in which the trees were found.
      index <- index + 1 # increments
    }
  }
  Ntrees <- index # Number of trees encountered.
  
  # converts back to a spatial grid:
  Cb <- imagery
  Mb <- imagery
  Cb[] <- as.numeric(Index[1:dim(Index)[1], dim(Index)[2]:1], 
                     byrow = TRUE) # given the correct orientation again.
  Cb[Cb==0]<-NA # Excludes non maxima pixels.
  
  m2 <- methods::as(Cb, "SpatialGridDataFrame")
  m3 <- raster::raster(m2, layer = 1)
  m3.shp <- raster::rasterToPoints(m3, fun = , dissolve = TRUE)
  
  IT <- sp::SpatialPointsDataFrame(m3.shp, data = data.frame(tree=m3.shp[,'layer']), 
                                   match.ID = F)
  return(IT)
}

itcIMG_ts<-function (imagery = NULL, win.min, win.max, z.min, z.max, TRESHSeed = 0.45, 
                     TRESHCrown = 0.55, DIST = 10, specT = 0) 
{
  
  # !!! You need to make some appropriate checks of the arguments here !!!
  
  # Blurs the chm:
  imagery <- raster::focal(imagery, w = matrix(1, 3, 3), 
                           fun = function(x) {
                             mean(x, na.rm = T)
                           })
  # Extracts the image data as a matrix:
  Max <- matrix(dim(imagery)[2], dim(imagery)[1], data = imagery[, 
                                                                 ], byrow = FALSE)
  # Flips the image the right way again:
  Max <- Max[1:dim(imagery)[2], dim(imagery)[1]:1]
  Gnew <- Max     # Copies the max matrix.
  Max[, ] <- 0    # Sets max to 0s. This will be used later for storing ...
  Index <- Max    # Copies max again.
  Index[, ] <- 0  # Sets index to 0s.
  Gnew[is.na(Gnew)] <- 0 # Sets any nas to 0s.
  Gnew[Gnew < specT] <- 0 # Any values beneath the minimum height are set to 0s.
  index = 1       # Initiates the index
  II <- which(Gnew != 0, arr.ind = T) # Extracts the locations of pixels which are bigger than the min tree height.
  dim(II)
  
  fm.win<-win.lmfun(win.min, win.max, z.min, z.max) # generates the relationship between height and search window size.
  # Extracts only the pixels that are sufficiently far from the image edge for the search window. In each direction.
  # The search window is selected according to the height of the tree:
  z<-Gnew[II[,1]+nrow(Gnew)*(II[,2]-1)] # extracts the tree heights from the matrix
  WinSize<-win.sizer(z, fm.win) # Finds the window size for each pixel
  half.WinSize<-ceiling(WinSize/2) # Halfs the window size for subsetting by pixels far enough from the image edge
  
  # Extracts only the pixels far enough from the image edge.
  II.ind<-II[, 1] >= half.WinSize & 
    II[, 1] <= (nrow(Gnew) - half.WinSize) &
    II[,2] >= half.WinSize &
    II[,2] <= (ncol(Gnew) - half.WinSize)
  II<-II[II.ind,]
  WinSize<-WinSize[II.ind]
  z<-z[II.ind]
  dim(II)
  
  # reorder from greatest to least according to z value (i.e. big trees first):
  z.order<-order(z, decreasing = TRUE)
  II<-II[z.order,]
  WinSize<-WinSize[z.order]
  
  # Works through each pixel one by one:
  indexII<-1
  #indexII<-(1:nrow(II)[1])[1]
  for (indexII in 1:dim(II)[1]) 
  {
    r = as.numeric(II[indexII, 1]) # Extracts the row pos.
    k = as.numeric(II[indexII, 2]) # Extracts the column pos.
    searchWinSize<-WinSize[indexII]
    hc.sWS<-ceiling(searchWinSize/2) # half the search window size rounded to floor
    hf.sWS<-floor(searchWinSize/2) # half the search window size rounded to ceiling
    FIL <- matrix(searchWinSize, searchWinSize, data = NA) 
    # Extracts the search window for the pixel:
    FIL <- Gnew[(r - hf.sWS):(r + hf.sWS), 
                (k - hf.sWS):(k + hf.sWS)]
    # Extracts the window from Max indicating whether a tree has already being designated or not:
    Max.chk<-Max[(r - hf.sWS):(r + hf.sWS),
                 (k - hf.sWS):(k + hf.sWS)]
    
    # If the focal pixel has the greatest value in the window & there is no tree already assigned in the output matrix within the window & the max value is no 0...
    # because the order is from tallest to shortest, large trees will always suppress the designation of small trees.
    if (FIL[hc.sWS, hc.sWS] == max(FIL, na.rm = T) &
        max(Max.chk, na.rm = T) == 0 &
        max(FIL, na.rm = T) != 0) 
    {
      Max[r, k] <- 1 # A logical assignment of tallest tree
      Index[r, k] <- index # Assigns the sequence in which the trees were found.
      index <- index + 1 # increments
    }
  }
  Ntrees <- index # Number of trees encountered.
  
  if (Ntrees > 0) 
  {
    # Extracts the chm values:
    Cb <- imagery
    Mb <- imagery
    Cb[] <- as.numeric(Gnew[1:dim(Gnew)[1], dim(Gnew)[2]:1], 
                       byrow = TRUE)
    Mb[] <- as.numeric(Max[1:dim(Max)[1], dim(Max)[2]:1], 
                       byrow = TRUE)
    
    Crowns <- Index # Assigns the crown sequence within the raster to Crowns
    OldCrowns <- Crowns # Copies crowns
    Check <- OldCrowns # Copies again.
    Check[, ] <- 0 # Sets check to 0.
    #filsize <- 3 # ????
    #Niter <- 100 # ????
    it = 1
    while (it == 1) 
    {
      it = 0
      II <- which(Crowns != 0 & Check == 0, arr.ind = T) # Extracts the crown pixels that have not been done yet; seems a bit inefficient.
      if (length(II) > 0)  # should be nrow and not length.
      {
        #indexII<-(1:nrow(II))[1]
        for (indexII in 1:dim(II)[1])  # Works through all the crown pixels simultaneously
        {
          r = as.numeric(II[indexII, 1])
          k = as.numeric(II[indexII, 2]) # Extracts the tree location.
          if (r != 1 & r != dim(Gnew)[1] & k != 1 & k != dim(Gnew)[2]) 
          { # So longs as the pixel is not right on the boundary; !!! this might need to be changed depending on window size.
            ind <- Crowns[r, k]
            coordSeed <- which(Index == ind, arr.ind = TRUE) # This finds the tree seed location.
            coordCrown <- which(Crowns == ind, arr.ind = TRUE) # ... this finds the current areas occupied by the crown.
            rvSeed <- Gnew[coordSeed] # Extracts the tree height.
            rvCrown <- mean(Gnew[coordCrown], na.rm = T) # Extracts the mean tree height.
            # Makes a matrix of the coordinates and chm values for the adjacent cells (surrounding the focal pixel).
            filData <- matrix(4, 3, data = 0)
            filData[1, 1] <- r - 1
            filData[1, 2] <- k
            filData[1, 3] <- Gnew[r - 1, k]
            filData[2, 1] <- r
            filData[2, 2] <- k - 1
            filData[2, 3] <- Gnew[r, k - 1]
            filData[3, 1] <- r
            filData[3, 2] <- k + 1
            filData[3, 3] <- Gnew[r, k + 1]
            filData[4, 1] <- r + 1
            filData[4, 2] <- k
            filData[4, 3] <- Gnew[r + 1, k]
            
            # Calculates distance of pixels from the focal pixel
            fil.dists<-as.matrix(dist(rbind(coordSeed, filData[,1:2])))[1,-1]
            
            # Checks which (if any) of the values are greater than the max or mean tree heights adjusted by the thresholds.
            # and less than 5% taller than the max tree height.
            # and the euclidean distance from the maximum crown radius is less than the specified DIST.
            GFIL <- (filData[, 3] > (rvSeed * TRESHSeed) & 
                       (filData[, 3] > (rvCrown * TRESHCrown)) & 
                       (filData[, 3] <= (rvSeed + (rvSeed * 0.05))) &
                       (fil.dists<DIST))
            #(abs(coordSeed[1] - filData[,1]) < DIST) &
            #(abs(coordSeed[2] - filData[, 2]) < DIST))
            filData <- filData[GFIL, ] # Subsets filData by the decision tree output.
            if (length(filData) > 3) 
            {
              # pp<-nrow(filData)[1]
              # Assigns each remaining pixel to the crown:
              for (pp in 1:dim(filData)[1]) 
              {
                rr <- filData[pp, 1]
                kk <- filData[pp, 2]
                if (Crowns[rr, kk] == 0 & Gnew[rr, kk] != 0) 
                {
                  Crowns[rr, kk] <- Crowns[r, k]
                  it <- 1
                }
              }
            }
          }
        }
      }
      Check <- OldCrowns # Marks the pixels that have been completed:
      OldCrowns <- Crowns # Recreates the oldcrowns matrix.
    }
    
    Cb <- imagery
    Mb <- imagery
    Cb[] <- as.numeric(Crowns[1:dim(Crowns)[1], dim(Crowns)[2]:1], 
                       byrow = TRUE)
    Mb[] <- as.numeric(Max[1:dim(Max)[1], dim(Max)[2]:1], 
                       byrow = TRUE)
    m2 <- methods::as(Cb, "SpatialGridDataFrame")
    m3 <- raster::raster(m2, layer = 1)
    m3.shp <- raster::rasterToPolygons(m3, fun = , dissolve = TRUE)
    
    # Extracts the mean and max tree heights:      
    names(m3.shp@data) <- "tree"
    HyperCrowns <- m3.shp[m3.shp@data[, 1] != 0, ]
    HCbuf <- rgeos::gBuffer(HyperCrowns, width = -res(imagery)[1]/2, byid = T)
    HyperCrowns<-HyperCrowns[HyperCrowns$tree %in% HCbuf$tree,] # excludes any trees that have been buffered out.
    ITCcv <- rgeos::gConvexHull(HCbuf, byid = T)
    ITCcvSD <- sp::SpatialPolygonsDataFrame(ITCcv, data = HyperCrowns@data, match.ID = F)
    #ITCcvSD <- sp::SpatialPolygonsDataFrame(HyperCrowns, data = HyperCrowns@data, match.ID = F)
    
    ITCcvSD$CA_m2 <- unlist(lapply(ITCcvSD@polygons, function(x) methods::slot(x, "area")))
    
    # Adds the tree height info:
    u.trees<-ITCcvSD$tree
    CH<-lapply(u.trees, function(indexIII)
    {
      CH<-Gnew[Crowns==indexIII]
      data.frame(tree=indexIII, CH_mean=mean(CH, na.rm=TRUE), CH_max=max(CH, na.rm=TRUE))
    }
    )
    CH<-do.call(rbind, CH)
    
    ITCcvSD$CH_mean <- CH$CH_mean
    ITCcvSD$CH_max <- CH$CH_max
    ITCcvSD <- ITCcvSD[ITCcvSD$CA_m2 > 1, ] # excludes the trees smaller than 1m^2 
    return<-ITCcvSD
    #      if (exists("ITCcvSD")) {
    #        return <- ITCcvSD[, -1] # excludes the tree label.
    #      }
  }
}