#############################################################################################################################################
###  Method for optimizing pedogenon classes as soil districts for the Basque Country
###  in the context of the European soil Monitoring Law
###  In this script: Modeling and optimization the number of clusters and covariate selection
###  Helper functions

### Desired extent: Basque Country
### Resolution: 25m
### CRS: EPSG=25830

###  Author: Mercedes Roman Dobarco
###  Date: 27/02/2024


# 1. Calculate centroids of soil properties by pedogenon and depth --------

### Function to calculate the inter-profile distance
### Calculated as average between inter-horizon distance
### Distance here is euclidean, we could also use Mahalanobis

### Short function to compute the centroids by pedogenon, and depth layer
centroids.SoilVars.Pedogenon.fun <- function(df.soil, depth.var, target.vars) {
  
  ### Rename variables
  colnames(df.soil)[colnames(df.soil) == depth.var] <- "Layer_depth"
  
  centroids.targetvars <- df.soil %>% 
    group_by(., PdGn, Layer_depth) %>%
    summarise_at(., .vars=target.vars, mean, na.rm = TRUE) %>%
    as.data.frame()
  
  return(centroids.targetvars)
  
}


# 2. Intra cluster distance - distance of each soil profile to its --------

### This function creates a dataframe with the soil profile distance (unique ID) 
### to the centroid of its respective pedogenon class
### Required input
### The input dataframe "df.soil" neeeds to have the columns:
### "uniqueID": ID that designates the unique combination of location (coordinates) and date ("newID"),
###             or just unique coordinates
### "Dataset": Basonet, LUCAS, Carbosol,etc.
### "depth.var": column that describes the depth intervals produced from splines (0-10, 10-20, 20-40, 40-60, 60-100 cm) 
###              or the two depth intervals from Basonet (0-20 and 20-40 cm)
### "PdGn": Pedogenon class, the column needs to be changed before applying this function

### target.vars refers to the soil properties considered in this exercise.
### depth.var: column that describes the depth interval

### The output, dataframe, needs to be averaged (through the third column) to calculate the Dintra average distance

Dintra.function <- function(df.soil, uniqueID, depth.var, target.vars) {
  
  ### Calculate distance between each observation to their centroid. (Euclidean).
  
  ### Variables for test
  # uniqueID <- "newID"
  # depth.var <-"Layer_depth" 
  # target.vars <- c("Silt","Clay","pH", "TOC")
  
  ### Rename variables
  colnames(df.soil)[colnames(df.soil) == depth.var] <- "Layer_depth"
  colnames(df.soil)[colnames(df.soil) == uniqueID] <- "uniqueID"
  
  ### Subset the variables of interest for this case
  df.soil <- df.soil[, c("uniqueID","Dataset","Layer_depth","PdGn", target.vars)]
 
  ### Subset only complete observations
  df.soil <- df.soil[complete.cases(df.soil),]
  
  ### Check number of depth layers available
  depths_here <- unique(df.soil[,"Layer_depth"]) 
  
  ### Unique Pedogenon classes
  PdGn.classes.here <- unique(df.soil[,"PdGn"]) 
  
  ### Rearrange by pedogenon, uniqueID and depth layer
  df.soil <- df.soil %>%
    arrange(., PdGn, uniqueID, Layer_depth) %>% as.data.frame()
  
  ### Calculate centroid by Pedogenon and depth interval
  centroids.targetvars <- df.soil %>% 
    group_by(., PdGn, Layer_depth) %>%
    summarise_at(., .vars=target.vars, mean, na.rm = TRUE) %>%
    as.data.frame()
  
  # ### Calculate all possible combination between pedogenon classes and depth layers.
  # cmbn.PdGn.Depth <- expand.grid(PdGn.classes.here,depths_here)
  # colnames(cmbn.PdGn.Depth) <- c("PdGn","Layer_depth")
  # ### Left join and attach centroids (can be NA)
  # centroids.expanded <- left_join(cmbn.PdGn.Depth, centroids.targetvars) %>%
  #   arrange(., PdGn, Layer_depth) %>% as.data.frame()
  
  ### All uniqueIDs in this dataframe
  all.uniqueIDs <- sort(unique(df.soil[,"uniqueID"]))
  length(all.uniqueIDs)
  
  ### Output dataframe will have uniqueID, PdGn, and dist_to_centroid
  Dist.Centroid.DF <- data.frame(uniqueID=numeric(), ### Soil profile ID
                                 PdGn=numeric(),     ### Pedogenon it belongs to
                                 dist_to_centroid=numeric()) ### Profile distance to pedogenon centroid
  
  ### for each uniqueID - compute soil profile distance to centroid
  ### For each depth interval and each pedogenon class:
  ### Create empty dataframes for the N depths and fill, 
  ### First row is the centroid
  ### Subsequent rows are all the observations for that pedogenon and depth including ALL unique IDs
  
  for(soil.profile in 1:length(all.uniqueIDs)) {
    
    ### subset soil profile
    soil.prof.i <- df.soil[df.soil$uniqueID ==all.uniqueIDs[[soil.profile]],]
    
    ### depth intervals with observations?
    soil.prof.depths <- unique(soil.prof.i$Layer_depth)
    
    ### vector with depths
    hor.distances.i <- c(rep(NA,length(soil.prof.depths)))
    
    for(depth in 1:length(soil.prof.depths)) {
      
      df.hor.j <- rbind(soil.prof.i[soil.prof.i$Layer_depth == soil.prof.depths[[depth]],c("PdGn","Layer_depth",target.vars)],
                        centroids.targetvars[centroids.targetvars$PdGn == unique(soil.prof.i$PdGn) &
                                               centroids.targetvars$Layer_depth == soil.prof.depths[[depth]], ])
      
      ### Euclidean horizon distance to horizon centroid
      dist.hor.j <- dist(df.hor.j[,target.vars])
      
      ### Assign to vector
      hor.distances.i[[depth]] <- dist.hor.j
      
    }
    
    ### Average of the horizon distances
    soil.prof.dist_to_centroid <- mean(hor.distances.i, na.rm=TRUE)
    
    ### Output vector
    soil.prof.Dintra <- data.frame("uniqueID" = unique(soil.prof.i$uniqueID), 
                                     "PdGn" =  unique(soil.prof.i$PdGn),
                                   "dist_to_centroid" = soil.prof.dist_to_centroid)
    
    ### Rbind to output dataframe
    Dist.Centroid.DF <- rbind(Dist.Centroid.DF,soil.prof.Dintra)
    
  }
  
    
  Dist.Centroid.DF <- Dist.Centroid.DF %>% arrange(., PdGn,uniqueID) %>% as.data.frame()
  
  return(Dist.Centroid.DF)
  
}



# 3. Inter-centroid distances, calculated from the soil profile di --------

### This function creates a dist object with the inter-centroid (from average soil profile distances)
### distances - Soil-profile Centroids are treated as individuals
### Required input
### The input dataframe "df.soil" neeeds to have the columns:
### "uniqueID": ID that designates the unique combination of location (coordinates) and date ("newID"),
###             or just unique coordinates
### "Dataset": Basonet, LUCAS, Carbosol,etc.
### "depth.var": column that describes the depth intervals produced from splines (0-10, 10-20, 20-40, 40-60, 60-100 cm) 
###              or the two depth intervals from Basonet (0-20 and 20-40 cm)
### "PdGn": Pedogenon class, the column needs to be changed before applying this function

### target.vars refers to the soil properties considered in this exercise.
### depth.var: column that describes the depth interval

### The output, a dist object, needs to be averaged to calculate the Dinter distance


Dinter.function <- function(df.soil, uniqueID, depth.var, target.vars) {
  
  
  ### Calculate distance between pedogenon centroids, from soil data
  
  ### Variables for test
  # df.soil <- Soil.df.subA
  # uniqueID <- "newID"
  # depth.var <-"Layer_depth" 
  # target.vars <- c("Silt","Clay","pH", "TOC")
  
  ### Rename variables
  colnames(df.soil)[colnames(df.soil) == depth.var] <- "Layer_depth"
  colnames(df.soil)[colnames(df.soil) == uniqueID] <- "uniqueID"
  
  ### Subset the variables of interest for this case
  df.soil <- df.soil[, c("uniqueID","Dataset","Layer_depth","PdGn", target.vars)]
  
  ### Subset only complete observations
  df.soil <- df.soil[complete.cases(df.soil),]
  
  ### Check number of depth layers available
  depths_here <- unique(df.soil[,"Layer_depth"]) 
  
  ### Unique Pedogenon classes
  PdGn.classes.here <- unique(df.soil[,"PdGn"]) 
  
  ### Rearrange by pedogenon,and depth layer
  df.soil <- df.soil %>%
    arrange(., PdGn, Layer_depth, uniqueID) %>% as.data.frame()
  
  ### Calculate centroid by Pedogenon and depth interval
  centroids.targetvars <- df.soil %>% 
    group_by(., PdGn, Layer_depth) %>%
    summarise_at(., .vars=target.vars, mean, na.rm = TRUE) %>%
    as.data.frame()
  
  ### Calculate all possible combination between pedogenon classes and depth layers.
  cmbn.PdGn.Depth <- expand.grid(PdGn.classes.here,depths_here)
  colnames(cmbn.PdGn.Depth) <- c("PdGn","Layer_depth")
  ### Left join and attach centroids (can be NA)
  centroids.expanded <- left_join(cmbn.PdGn.Depth, centroids.targetvars) %>%
    arrange(., PdGn, Layer_depth) %>% as.data.frame()

  ### subset by depth layers in depths_here
  centroids.depths <- split(centroids.expanded, centroids.expanded$Layer_depth)
  
  ### Within each element of the list, this is, horizon, calculate distance between centroids
  centroids.dists <- lapply(centroids.depths, 
                            FUN = function(x, target.vars) {as.matrix(dist(x[,target.vars]))} )
  
  ### Create an array from the list
  centroids.dist.array <- array( data = unlist(centroids.dists), 
                                 dim = c(nrow(centroids.dists[[1]]),
                                         ncol(centroids.dists[[1]]),
                                         length(centroids.dists)))
  
  
  ### Average inter-cluster centroid profile-distances, from the array
  centroids.ave.prfl.dist <- apply(centroids.dist.array, 
                                   MARGIN = c(1,2), 
                                   function(x) mean(x,na.rm=TRUE))
  
  ### Transform into dist object
  centroids.ave.prfl.d <- as.dist(centroids.ave.prfl.dist)
  
  return(centroids.ave.prfl.d)
  
}


### End of the script