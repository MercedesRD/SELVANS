#############################################################################################################################################
###  Method for optimizing pedogenon classes as soil districts for the Basque Country
###  in the context of the European soil Monitoring Law
###  In this script: Modeling and optimization the number of clusters and covariate selection

### Desired extent: Basque Country
### Resolution: 25m
### CRS: EPSG=25830

###  Author: Mercedes Roman Dobarco
###  Date: 15/02/2024

####### Load packages
### Spatial
library(sf)
library(terra)
library(gdalUtilities)
library(XML)
library(reproducible)

### Visualization
library(lattice)
library(ggplot2)
library(viridis) # color palettes
library(scales)
library(rasterVis)
library(gridExtra)
library(rasterVis)
library(RColorBrewer)
library(tmap)    # for static and interactive maps
library(leaflet) # for interactive maps
library(mapview) # for interactive maps
library(shiny)   # for web applications

###Data carpentry
library(dplyr)
library(tidyverse)

### Parallel computing
library(foreach)
library(parallel)
library(doParallel)

### clustering
#install.packages("remotes")
#remotes::install_github("andrewthomasjones/tkmeans")
library(ClusterR)
library(clusterSim)

### Steps in pedogenon mapping:

### Load the helper functions
source("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_scripts/Euskadi/5.PedogenonModeling_helper.R")

# ### 1. Tentative number of pedogenons -----------------------------------

### 1. Tentative number of pedogenons determined by:
### 1.a combination of soil order, environmental zone (1), geology and major vegetation group.
### 1.b number of soil map units from traditional soil map.
### 1.c Guo et al (2003) equations applied to the Basque Country
### 1.d Minimum number of observations per pedogenon class required for validation (~10)
### something between 2 to 140 classes
search_space <- c(2:140)

# ### 2. Selection of SCORPAN variables --------------------------------------

### Number of variables per soil-forming factor:
### Soil variables: clay/silt, clay/sand, silt/sand, CEC - between 2-3
### the spatial patterns of the SoilGrids mapsdiffered considerably from the 
### texture maps elaborated by the Basque Government in 2019, at a smaller resolution and produced
### with a regional dataset. 
### The particle size fractions and SOC stock maps for the 0-30 cm depth interval were produced with
### around 12,000 observationsand covariates of parent material, land use, climate, and relief, following 
### a digital soil mapping approach and scorpan modeling. 
### (https://www.euskadi.eus/mapa-de-existencias-de-carbono-y-mapa-de-textura-para-los-suelos-de-la-capv/web01-a2inglur/es/)

### Climate = 5 variables selected from correlation plots 
### and variable importance for predicting soil properties
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/4.CovariateSelection/scaled/")
clim.vars <- c("clim_bio1.tif","clim_bio4.tif","clim_bio5.tif","clim_bio12.tif","clim_bio15.tif")
clim.r <- rast(clim.vars)
plot(clim.r)

### Relief: between 5 to 10 - These will be selected with iterations
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/4.CovariateSelection/scaled/")
relief.fixed <- c("relief_dem.tif")
#dem <- rast(relief.fixed)

### Variables that logically limit soil CAPACITY for performing functions 
### and their resilience to disturbance:
relief.potential <-list.files(pattern="relief")
relief.potential <- c("relief_slope.tif","relief_slope_5.tif","relief_slope_10.tif", # slope
                      "relief_northerness.tif","relief_easterness.tif",  ### aspect
                      "relief_north_slope.tif", ### aspect x slope. This variable is correlated with northness so we may want to keep one of both
                      "relief_twi.tif","relief_mrrtf.tif","relief_mrvbf.tif", # my preferred is TWI 
                      "relief_valley_depth.tif","relief_st_height.tif", # Preferred hydrological variables
                      "relief_mid_slope.tif","relief_norm_height.tif","relief_slope_height.tif", ### Other hydrological variables
                      "relief_p_curv.tif", "relief_pr_curv_5.tif", "relief_pr_curv_10.tif", # profile curvature (3,5,10)
                      "relief_pl_curv_5.tif", "relief_pl_curv_10.tif", # planar curvature (5,10)
                      "relief_lg_curv_5.tif", "relief_lg_curv_10.tif", # longitudinal curvature (5,10)
                      "relief_t_curv.tif", "relief_cs_curv_5.tif", "relief_cs_curv_10.tif", # tangential curvature (3) cross-sectional curvature (5,10) (same)
                      "relief_tpi_8_3.tif") ## I eliminate "relief_tpi_20_5.tif",
relief.r <- rast(c(relief.fixed,relief.potential))
names(relief.r) <- paste0("r_",names(relief.r))

### Vegetation: 4 or 5 PCs of MVG PCA scores (54-66 % variance)
setwd("C:/Covariates/Euskadi/Organisms/")
vegetation.vars <- c("MVG_PC1.tif","MVG_PC2.tif","MVG_PC3.tif","MVG_PC4.tif","MVG_PC5.tif")
vegetation.r <- rast(vegetation.vars)
names(vegetation.r) <- paste0("o_",names(vegetation.r))

### Parent material: 8 or 50 % of variance
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/1.CovariatesEus/")
parentmaterial.vars <- c("pm_mca_1.tif","pm_mca_2.tif","pm_mca_3.tif",
                         "pm_mca_4.tif","pm_mca_5.tif","pm_mca_6.tif",
                         "pm_mca_7.tif","pm_mca_8.tif")
parentmaterial.r <- rast(parentmaterial.vars)

### All variables
scorpan <- c(clim.r, vegetation.r, relief.r, parentmaterial.r)

### 2.a Take sample
### Regular sample
set.seed(2233)
scorpanSample <- terra::spatSample(x = scorpan, 
                                   size=50000,
                                   method="regular", 
                                   as.df=TRUE, 
                                   xy=TRUE)
### Only complete cases
scorpanSample <- scorpanSample[complete.cases(scorpanSample),]
dim(scorpanSample)
coords.scorpanSample <- scorpanSample[,1:2] ### Keep coordinates separately


# ### 2.b Extract covariates at soil observations -------------------------

### 2.b Extract covariates at soil observations, that will be also used to

### Table with soil data and covariates, with same names as in the dataframe used for clustering
### Load the soil datasets
load("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/1.SoilDatasets/SoilDataInputPedM.RData")
### I keep Soil.df.harmonised and Soil.df.harmonised_sf - splined data
### soil_datasets_3 (before splines and with correct coordinates for BASONET) is also an option

### Think if I want to use newID or just myID - all years? just 2021?
### Same for LUCAS - subset one year?

### Transform soil_datasets_3 to spatial, extracting only BASONET data for either 
### for 2001 (but coordinates may not be entirely correct)
### or 2021 (coordinates are correct but I don´t know if I put the texture data too)
BASONET.soil <- soil_datasets_3[soil_datasets_3$Dataset== "BASONET",]
BASONET.soil.sf <-st_as_sf(BASONET.soil, coords = c("UTM_X","UTM_Y"), crs = 25830)
#BASONET.soil.WGS84 <-  st_transform(BASONET.soil.sf, 4326)

### Extract
soil.covariates <- terra::extract(x=scorpan, y=Soil.df.harmonised_sf, method="simple")
BASONET.covariates <- terra::extract(x=scorpan, y=BASONET.soil.sf, method="simple")

### Bind to the soil data
Soil.df.harmonised.scorpan <- cbind(Soil.df.harmonised,soil.covariates)
BASONET.scorpan <- cbind(BASONET.soil,BASONET.covariates)
summary(Soil.df.harmonised.scorpan)
summary(BASONET.scorpan)

### We can work with these two dataframes:
### Soil.df.harmonised.scorpan  Splined dataframe, all datasets
### BASONET.scorpan Basonet dataframe, before splines (0-20 cm and 20-40 cm)

### Should we consider all years from BASONET and LUCAS?

### TEXTURE from 2009 for LUCAS
unique(Soil.df.harmonised.scorpan[Soil.df.harmonised.scorpan$Dataset =="LUCAS" & !is.na(Soil.df.harmonised.scorpan$Silt),]$Date)

### For BASONET, let´s use 2021 in terms of coordinates, they are more accurate
unique(BASONET.scorpan[!is.na(BASONET.scorpan$Silt),]$Date) 
### Texture has been duplicated in 2021, no problem

### I exclude BASONET 2001 from both datasets to eliminate duplicity for texture on Pedogenon centroids
Soil.df.harmonised.scorpan <- Soil.df.harmonised.scorpan[!(Soil.df.harmonised.scorpan$Dataset == "BASONET" &
                                                             Soil.df.harmonised.scorpan$Date == 2001),]
BASONET.scorpan <- BASONET.scorpan[BASONET.scorpan$Date != 2001,]
### With LUCAS there is no problem because TEXTURE is only available for 2009

### Eliminate rows without covariates
Soil.df.harmonised.scorpan <- Soil.df.harmonised.scorpan[complete.cases(Soil.df.harmonised.scorpan[,colnames(soil.covariates)]),]
BASONET.scorpan <- BASONET.scorpan[complete.cases(BASONET.scorpan[,colnames(soil.covariates)]),]

### Which climate variables are correlated?
par(mfrow=c(1,1))
library(Hmisc)
library(corrplot)
Soil.df.harmonised.scorpan[,c("Sand","Silt","Clay","pH","TOC")] %>%
  cor(., use = "pairwise.complete.obs") %>%
  corrplot.mixed(.,upper = "ellipse", 
                 lower = "number", 
                 number.cex=0.7, tl.cex=0.6, tl.col = "black")

BASONET.scorpan[,c("Sand","Silt","Clay","CECef",
                   "TOC","Carbonates","pH","EC",               
                   "TN_ppm","P_ppm","K_ppm",
                   "Ca_ppm","Mg_ppm","Na_ppm")] %>%
  cor(., use = "pairwise.complete.obs") %>%
  corrplot.mixed(.,upper = "ellipse", 
                 order = 'hclust',
                 lower = "number", 
                 number.cex=0.7, 
                 tl.cex=0.6, 
                 tl.col = "black")

### Do not use P_ppm, Na_ppm, Mg_ppm because there are a lot of missing observations
### I use CECef and exclude Ca_ppm because they are highly correlated
### I use TOC and exclude TN_ppm because they are correlated
### Then I use Silt, Clay, pH, CECef, EC, TOC.



# ### 2.c Create all combinations of relief covariates --------------------


### 2.c Create all combinations of relief covariates

### Taking n=5, ... n =10
### function combn from package utils
library(utils)
relief.potential <- names(relief.r)[-1] ### Take out DEM, which will be fixed

### My choice, "expernt knowledge" or just my personal preference
relief.preferred <- c("r_slope","r_northerness","r_easterness","r_north_slope",
                      "r_twi","r_valley_depth","r_st_height")

### Split by 3,5,10 window
relief.potential.3 <-  c("r_slope",
                         "r_northerness","r_easterness","r_north_slope",
                         "r_twi","r_mrrtf","r_mrvbf",
                         "r_valley_depth","r_st_height","r_mid_slope",
                         "r_norm_height","r_slope_height",
                         "r_p_curv","r_t_curv","r_tpi_8_3")

relief.potential.5 <-  c("r_slope_5",
                         "r_northerness","r_easterness","r_north_slope",
                         "r_twi","r_mrrtf","r_mrvbf",
                         "r_valley_depth","r_st_height","r_mid_slope",
                         "r_norm_height","r_slope_height",
                         "r_pr_curv_5","r_pl_curv_5","r_lg_curv_5","r_cs_curv_5",
                         "r_tpi_8_3")

relief.potential.10 <-  c("r_slope_10",
                         "r_northerness","r_easterness","r_north_slope",
                         "r_twi","r_mrrtf","r_mrvbf",
                         "r_valley_depth","r_st_height","r_mid_slope",
                         "r_norm_height","r_slope_height",
                         "r_pr_curv_10","r_pl_curv_10","r_lg_curv_10","r_cs_curv_10",
                         "r_tpi_8_3")
  
### n=4
### I start from my subset of preferred variables
combn.4.relief <- utils::combn(x=relief.preferred, m = 4)
### this gives 35 combinations
combn.4.relief <- t(combn.4.relief)

### eliminate those where "r_north_slope" and  "r_northerness" are together
subv <- c("r_northerness","r_north_slope")
conds <- apply(X = combn.4.relief, MARGIN = 1,  FUN = function(x){sum(subv %in% x)!=2})
combn.4.relief <- combn.4.relief[conds,]
### This reduces that to 25 possible combinations
combn.4.relief <- cbind(combn.4.relief, matrix(data=NA, nrow=nrow(combn.4.relief), ncol = 2))

combn.5.relief <- utils::combn(x=relief.preferred, m = 5)
### this gives 21 combinations
combn.5.relief <- t(combn.5.relief)

### eliminate those where "r_north_slope" and  "r_northerness" are together
subv <- c("r_northerness","r_north_slope")
conds <- apply(X = combn.5.relief, MARGIN = 1,  FUN = function(x){sum(subv %in% x)!=2})
combn.5.relief <- combn.5.relief[conds,]
combn.5.relief <- cbind(combn.5.relief, matrix(data=NA, nrow=nrow(combn.5.relief), ncol = 1))

combn.6.relief <- utils::combn(x=relief.preferred, m = 6)
### this gives 21 combinations
combn.6.relief <- t(combn.6.relief)

### eliminate those where "r_north_slope" and  "r_northerness" are together
subv <- c("r_northerness","r_north_slope")
conds <- apply(X = combn.6.relief, MARGIN = 1,  FUN = function(x){sum(subv %in% x)!=2})
combn.6.relief <- combn.6.relief[conds,]

### bind matrices
relief.combi <- rbind(combn.4.relief,combn.5.relief,combn.6.relief)
### (expand later with more possible combinations of curvatures, etc.)
### Collapse to vector
relief.combi.v <- apply(X = relief.combi, MARGIN = 1,  FUN = function(x){paste(x,collapse=",")})

### Keep the combination in a cata.frame with a code
relief.combi.df <- data.frame (ID = 1:nrow(relief.combi), relief.combi = relief.combi.v)

### fixed variables
fixed.columns <- c(names(clim.r), names(vegetation.r), "r_dem", names(parentmaterial.r))

### directory to store the models
OutDir <- "C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModeling/"

### Delete temporal files
tmpFiles(current = FALSE,orphan = TRUE,old = TRUE,remove = TRUE)

# ### 3. For each combination of relief covariates: -----------------------


### 3. For each combination of relief covariates:

kquality_combi_list <- list() ### List with the dataframes of internal quality and soil-profile derived indices

#icq_df_combi_list <- list() ### this list will store the internal clustering quality indices
### for each combination of environmental covariates

#BASONET_df_combi_list <- list() ### this list will store the soil profile indices - 
### Only BASONET 2021 data
### for each combination of environmental covariates

#SoilHarmonised_df_combi_list <- list() ### this list will store the soil profile indices 
### All datasets
### for each combination of environmental covariates

#cluster_assignment_combi_list <- list() ### this list will store the dataframe for combination j,
### with all the cluster assignments from the search space

### List to store the dataframes with BASONET and HARMONISED soil datasets cluster assignments for different j combinations
BASONET.dfs <- list()
HARMONISED.dfs <- list()

### Create a for loop to test the different metrics with different soil datasets
#j <- 5

for(j in 1:nrow(relief.combi)) { 
  setwd(OutDir)
  dir.create(paste0("scorpan_combi_",j))
  }

# 
# cl <- makePSOCKcluster(20) # I defined cl by this commend
# registerDoParallel(cl)
# clusterExport(cl, varlist = c("terra","ClusterR", "clusterSim","lowmemtkmeans", "tidyverse","dplyr"))
# clusterEvalQ(cl, .libPaths("~/C:/Users/mercedes.roman/AppData/Local/R/win-library/4.3")) # pass libpath
# clusterEvalQ(cl, library(ClusterR)) # pass My package which includes Rcpp functions
# stopCluster(cl)
# print("FinishCL")# for testing that cl works
library(doParallel)
library(foreach)

tic <- Sys.time()
detectCores()
cl <- makeCluster(13)   ###
registerDoParallel(cl)
getDoParWorkers()

kquality_scorpan_combis <- foreach(j = 1:nrow(relief.combi),
                                   .packages=c("terra","ClusterR", "clusterSim","lowmemtkmeans", "tidyverse","dplyr"),
                                   .export = c("OutDir", "fixed.columns","relief.combi",
                                               "relief.combi.df", "relief.combi.v",
                                               "scorpanSample", "BASONET.scorpan", "search_space",
                                               "Soil.df.harmonised.scorpan","Dintra.function", "Dinter.function")) %dopar% { 
  
#for(j in 1:nrow(relief.combi)) {
  
# ### a. SCORPAN sample, input for clustering for this combination --------

  ### a. SCORPAN sample, input for clustering for this combination of covariates
  setwd(paste0(OutDir,"scorpan_combi_",j,"/"))
  print(paste0("Working on scorpan combination ",j))
  
  ### Subset variables in scorpan dataset
  #fixed.columns <- c(names(clim.r), names(vegetation.r), "r_dem", names(parentmaterial.r))
  var.columns <- c(t(relief.combi)[,j])
  var.columns <- var.columns[!is.na(var.columns)] ### eliminate NA if any
  
  ### scorpan dataframe for the combination j
  scorpanSample.j <- scorpanSample[,c(fixed.columns,var.columns)]
  
  ### a.1 Perform Cholesky transformation to decorrelate the data
  
  # The basic Euclidean distance treats each variable as equally important in calculating the distance.
  # An alternative approach is to scale the contribution of individual variables to the distance value according
  # to the variability of each variable. This approach is illustrated by the Mahalanobis distance, 
  # which is a measure of the distance between each observation in a multidimensional cloud of points and
  # the centroid of the cloud.
  ### Calculate the Mahalanobis distance, as Euclidean distance after applying the Cholesky decomposition
  
  # # Rescale the data
  C.j <- chol(var(as.matrix(scorpanSample.j)))
  scorpanSample.j.rs <- as.matrix(scorpanSample.j) %*% solve(C.j)
  
  ### Output of transformation of the sampled scorpan dataset for variable combination j: scorpanSample.j.rs
  # scorpanSample.j.rs
  ### Variance-Covariance matrix of the scorpan dataset for variable combination j: C.j
  # C.j
  
  ### I may be better working with arrays (ncombiJ x nIndices x search_space_K). 
  ### Dataframe to store the results of the clustering - internal quality indices
  ### These will be the results for this variable combination
  out.clustering.indices <- as.data.frame(matrix(data=NA, ncol=length(search_space)+1, nrow= 5))
  colnames(out.clustering.indices) <- c("Index", paste0("K.",search_space))
  #out.clustering.indices$Index <- c("Total SSE", "Sum WCSE","BetweenSS to TotalSS","Calinski-Harbasz","Silhouette","BIC")
  out.clustering.indices$Index <- c("Total SSE", "Sum WCSE","BetweenSS to TotalSS","Calinski-Harbasz","BIC")

#   ### b. JUST BASONET - 430 sites ---------------------------------------

  
  ### b. JUST BASONET - 430 sites
  ### N = 430 soil profiles with 1 or 2 horizons
  
  ### copy the soil dataframe
  BASONET.scorpan.j <- as.data.frame(BASONET.scorpan)
  
  ### Columns to store cluster assignment
  kassignments <- as.data.frame(matrix(data=NA,
                                       nrow=nrow(BASONET.scorpan.j),
                                       ncol = length(search_space)))
  colnames(kassignments) <- paste0("K.",search_space)
  knames <- paste0("K.",search_space)
  BASONET.scorpan.j <- cbind(BASONET.scorpan.j, kassignments)
  rm(kassignments)
  
  ### Now decorrelate the SCORPAN variables at the locations of the soil observations
  ### scorpan variables corresponding to combination j
  B.scorpan.j <- BASONET.scorpan.j[,c(fixed.columns,var.columns)]
  ### Rescale the data
  B.scorpan.j.rs <- as.matrix(B.scorpan.j) %*% solve(C.j)
  
  ### Dataframe to store output quality indices
  
  ### I may be better working with arrays (ncombiJ x nIndices x search_space_K). 
  ### Dataframe to store the results of the clustering - soil data profile indicators - BASONET
  ### These will be the results for this variable combination
  out.indices.BASONET <- as.data.frame(matrix(data=NA, ncol=length(search_space)+1, nrow= 6))
  colnames(out.indices.BASONET) <- c("Index",paste0("K.",search_space))
  out.indices.BASONET$Index <- c("B_N_PdGn_sites", # Number of pedogenons with any soil observation
                                 "B_Perc_PdGn_sites", ## What percentage of the number of classes do have observations?
                                 "B_Min_sites", # Min number of observations per pedogenon (of those with any)
                                 "B_Median_sites", # Median number of observations per pedogenon (of those with any)
                                 "B_Max_sites", # Min number of observations per pedogenon (of those with any)
                                 "B_Din_Dex" # Ratio of SOIL PROFILE intra cluster distance to between cluster distance
                                 ### Din is average distance from each observation to the centroid of the cluster
                                 ### Dex is the average distance between cluster centroids
                                 )
  
#   ###  c. ALL DATASETS -  -----------------------------------------------

  ###  c. ALL DATASETS - 
  
  ### copy the soil dataframe
  Soil.df.harmonised.scorpan.j <- as.data.frame(Soil.df.harmonised.scorpan)
  
  ### Columns to store cluster assignment
  kassignments <- as.data.frame(matrix(data=NA,
                                       nrow=nrow(Soil.df.harmonised.scorpan.j),
                                       ncol = length(search_space)))
  colnames(kassignments) <- paste0("K.",search_space)
  knames <- paste0("K.",search_space)
  Soil.df.harmonised.scorpan.j <- cbind(Soil.df.harmonised.scorpan.j, kassignments)
  rm(kassignments)
  
  ### Now decorrelate the SCORPAN variables at the locations of the soil observations
  ### scorpan variables corresponding to combination j
  SH.scorpan.j <- Soil.df.harmonised.scorpan.j[,c(fixed.columns,var.columns)]
  ### Rescale the data
  SH.scorpan.j.rs <- as.matrix(SH.scorpan.j) %*% solve(C.j)
  
  ### I may be better working with arrays (ncombiJ x nIndices x search_space_K). 
  ### Dataframe to store the results of the clustering - soil data profile indicators - 
  ### subset of selected properties within pedogenon: sand, silt, clay, TOC, pH, Mg, Ca, K, N
  ### These will be the results for this variable combination
  out.indices.SoilHarmonised <- as.data.frame(matrix(data=NA, ncol=length(search_space)+1, nrow= 6))
  colnames(out.indices.SoilHarmonised) <- c("Index",paste0("K.",search_space))
  out.indices.SoilHarmonised$Index <- c("SH_N_PdGn_sites", # Number of pedogenons with any soil observation
                                        "SH_Perc_PdGn_sites", ## What percentage of the number of classes do have observations?
                                        "SH_Min_sites", # Min number of observations per pedogenon (of those with any)
                                        "SH_Median_sites", # Median number of observations per pedogenon (of those with any)
                                        "SH_Max_sites", # Min number of observations per pedogenon (of those with any)
                                        "SH_Din_Dex") # Ratio of SOIL PROFILE intra cluster distance to between cluster distance
  ### Din is average distance from each observation to the centroid of the cluster
  ### Dex is the average distance between cluster centroids
 
  setwd(paste0(OutDir,"scorpan_combi_",j,"/"))
  gc()

#   ### 4. for each k in the search space, search_space: ------------------

  ### 4. for each k in the search space, search_space:
  #k <- 30
 for( k in 1:length(search_space)){
    
   print(paste0("Calculating indices for k=",search_space[[k]]))
    
    ### Note, I can repeat this step 30 or 100 times (or as many as I want)
    ### by changing the seed for clustering, to obtain median estimates of the clustering indices
    
    ### 4.a Run k-means clustering and calculate
    set.seed(1990)
    kmeans_clorpt.jk <- ClusterR::KMeans_rcpp(scorpanSample.j.rs, 
                                              clusters = search_space[[k]], 
                                              num_init = 20, 
                                              max_iters = 10000,
                                              fuzzy = FALSE,
                                              initializer = 'kmeans++', 
                                              verbose = F)
    
    ### 4.b total SSE, sum of within cluster SE, between-cluster SSE / total SSE
    ### Output from function ClusterR::KMeans_rcpp
    Tot_SSE_jk <- kmeans_clorpt.jk$total_SSE
    sumWCSE_jk <- sum(kmeans_clorpt.jk$WCSS_per_cluster,na.rm=TRUE)   
    Btwcse_jk <- kmeans_clorpt.jk$between.SS_DIV_total.SS
    
    ### 4.c Internal cluster quality indices like:
    
    require(clusterSim)
    ### - Calinski-Harbasz 
    icqG1.jk <- clusterSim::index.G1(x=scorpanSample.j.rs,
                                     cl=kmeans_clorpt.jk$clusters, 
                                     centrotypes="centroids")
    
    ### - Silhouette - This is very memory demanding, so I use a subset of the data of 20,000 observations
    # icqS.jk <- clusterSim::index.S(d=dist(scorpanSample.j.rs), 
    #                                cl=kmeans_clorpt.jk$clusters)
    
    ### - Bayesian information criterion to penalize larger number of clusters (cluster_BIC {lowmemtkmeans})
    require(lowmemtkmeans)
    BIC.jk <- lowmemtkmeans::cluster_BIC(data=as.matrix(scorpanSample.j.rs),
                                         centres=as.matrix(kmeans_clorpt.jk$centroids))
    
    ### 4.d Store the results of the clustering indices
    #out.clustering.indices[,k+1] <- c(Tot_SSE_jk,sumWCSE_jk,Btwcse_jk,icqG1.jk,icqS.jk,BIC.jk)
    out.clustering.indices[,k+1] <- c(Tot_SSE_jk,sumWCSE_jk,Btwcse_jk,icqG1.jk,BIC.jk)
    
    # Save the k-means model and centroids
    save(kmeans_clorpt.jk, file=paste0(OutDir,"scorpan_combi_",j,"/kmeans_scorpanID",j,".k",search_space[[k]],".RData"))
    gc()
    
# ### 5. Predict cluster assignment to soil profiles BASONET--------

    ### 5. Predict cluster assignment to soil properties observations - BASONET dataset:
    
    ### Extract the index of the dataframe rows that are na/nan/Inf
    df.na <- which(apply(B.scorpan.j.rs, 
                         MARGIN = 1, 
                         FUN = function(x) {any(is.na(x))}))
    
    if(length(df.na) ==0) {
      
    ### Predict cluster assignment - BASONET
    cluster  <- predict_KMeans(data = B.scorpan.j.rs, CENTROIDS = kmeans_clorpt.jk$centroids)
    ### Assign to the dataframe with soil observations
    BASONET.scorpan.j[,knames[[k]]] <- cluster
    
    } else if (length(df.na) > 0) {
      
      cluster  <- predict_KMeans(data = B.scorpan.j.rs[-df.na,], CENTROIDS = kmeans_clorpt.jk$centroids)
      BASONET.scorpan.j[-df.na, knames[[k]]] <- cluster
      
    }
    
    ### 6. Summarise number of observations per cluster: average, min and max.
    ### Here I only use observations from 2021, so there is no problem of double counting same coordinates and different years
    
    ### 6.a Subset of soil properties, those "more stable" - for BASONET, I decided (seeing also the correlation plots)
    target.vars.BASONET <- c("Silt","Clay","CECef")
    
    ### Create the variable "Layer_depth"
    BASONET.scorpan.j$Layer_depth <- ifelse(BASONET.scorpan.j$Lower_limit == 19, "000_020_cm", "020_040_cm" )
    
    ### subset data for the Dintra and Dinter calculations
    Soil.df.BASONET <- BASONET.scorpan.j[,c("newID","Dataset","Layer_depth","Date", ### In this case either newID or myID design unique location
                                            knames[[k]],
                                            target.vars.BASONET)]
    
    ### Change name of pedogenon column
    colnames(Soil.df.BASONET)[colnames(Soil.df.BASONET) ==knames[[k]]] <- "PdGn"
    
    ### Subset only complete observations
    Soil.df.BASONET <- Soil.df.BASONET[complete.cases(Soil.df.BASONET),]
    
    ### Number of observations per pedogenon? Individual locations (unique coordinates + date)
    ### subset only one year for BASONET and LUCAS - Done
    summary.PdGn.BASONET <- Soil.df.BASONET[,c("newID","PdGn")] %>% distinct(.,newID,PdGn) %>% count(., PdGn)
    
    #out.indices.BASONET$Index 
    #hist(summary.PdGn.BASONET$n, breaks=20)
    ### How many of the pedogenons have any observation?
    B_N_PdGn_sites <- length(unique(summary.PdGn.BASONET$PdGn))
    ### What percentage does this represent from all the classes?
    B_Perc_PdGn_sites <- round(B_N_PdGn_sites/search_space[[k]]*100, digits=1)
    B_Min_sites <- min(summary.PdGn.BASONET$n)
    B_Median_sites <- median(summary.PdGn.BASONET$n)
    B_Max_sites <- max(summary.PdGn.BASONET$n)
    

# 6.b Dintra/Dinter BASONET -----------------------------------------------

    
    ### 6.b Din/Dex
    ### Din is average distance from each observation to the centroid of its cluster
    ### Dex is the average distance between cluster centroids
    ### Calculate with functions from "5.PedogenonModeling_helper.R"
    
    BASONET.DF.intra <- Dintra.function(df.soil = Soil.df.BASONET,
                                      uniqueID = "newID",
                                      depth.var = "Layer_depth",
                                      target.vars = target.vars.BASONET)
    
    ### Average Distance between each observation to their centroid.
    Dintra.BASONET <- mean(BASONET.DF.intra$dist_to_centroid, na.rm=TRUE) 
  
    
    BASONET.DF.inter <- Dinter.function(df.soil = Soil.df.BASONET,
                                        uniqueID = "newID",
                                        depth.var = "Layer_depth",
                                        target.vars = target.vars.BASONET)
    
    ### Calculate average distance between centroids.
    ### These are the distances between centroids.
    Dinter.BASONET <- mean(BASONET.DF.inter, na.rm=TRUE)
    
    ### Ratio Din to Dex
    Din_Dex_BASONET <- Dintra.BASONET/Dinter.BASONET
    
    ### Store the results of the soil profile distances indices
    out.indices.BASONET[,k+1] <- c(round(B_N_PdGn_sites, digits=0),
                                   round(B_Perc_PdGn_sites, digits=1),
                                   round(B_Min_sites, digits=0),
                                   round(B_Median_sites, digits=0),
                                   round(B_Max_sites, digits=0),
                                   round(Din_Dex_BASONET, digits=3))
    
    gc()
                                  
    
# ### 7. Predict cluster assignment to soil profiles HARMONISED DATASET--------
    
    ### 7. Predict cluster assignment to soil properties observations - HARMONISED dataset:
    
    ### Extract the index of the dataframe rows that are na/nan/Inf
    df.na <- which(apply(SH.scorpan.j.rs, 
                         MARGIN = 1, 
                         FUN = function(x) {any(is.na(x))}))
    
    if(length(df.na) ==0) {
      
      ### Predict cluster assignment - HARMONISED
      cluster  <- predict_KMeans(data = SH.scorpan.j.rs, CENTROIDS = kmeans_clorpt.jk$centroids)
      
      ### Assign to the dataframe with soil observations
      Soil.df.harmonised.scorpan.j[,knames[[k]]] <- cluster
      
    } else if (length(df.na) > 0) {
      
      cluster  <- predict_KMeans(data = SH.scorpan.j.rs[-df.na,], CENTROIDS = kmeans_clorpt.jk$centroids)
      Soil.df.harmonised.scorpan.j[-df.na, knames[[k]]] <- cluster
      
    }
    
    ### 6. Summarise number of observations per cluster: average, min and max.
    ### Here I only use observations from 2021, so there is no problem of double counting same coordinates and different years
    ### and because I require soil texture, only LUCAS 2009 is taken into account
    
    ### 6.a Subset of soil properties, those with more observations 
    target.vars.HARMONISED <- c("Silt","Clay","TOC","pH")
  
    ### subset data for the Dintra and Dinter calculations
    Soil.df.HARMONISED <- Soil.df.harmonised.scorpan.j[,c("newID","Dataset","Layer_depth","Date",
                                                          ### In this case either newID or myID design unique location
                                                          knames[[k]],
                                                          target.vars.HARMONISED)]
    
    ### Change name of pedogenon column
    colnames(Soil.df.HARMONISED)[colnames(Soil.df.HARMONISED) ==knames[[k]]] <- "PdGn"
    
    ### Subset only complete observations
    Soil.df.HARMONISED <- Soil.df.HARMONISED[complete.cases(Soil.df.HARMONISED),]
    
    ### Number of observations per pedogenon? Individual locations (unique coordinates + date)
    ### subset only one year for BASONET and LUCAS - Done
    summary.PdGn.HARMONISED <- Soil.df.HARMONISED[,c("newID","PdGn")] %>% distinct(.,newID,PdGn) %>% count(., PdGn)
    
    #out.indices.SoilHarmonised$Index 
    #hist(summary.PdGn.BASONET$n, breaks=20)
    ### How many of the pedogenons have any observation?
    SH_N_PdGn_sites <- length(unique(summary.PdGn.HARMONISED$PdGn))
    ### What percentage does this represent from all the classes?
    SH_Perc_PdGn_sites <- round(SH_N_PdGn_sites/search_space[[k]]*100, digits=1)
    SH_Min_sites <- min(summary.PdGn.HARMONISED$n)
    SH_Median_sites <- median(summary.PdGn.HARMONISED$n)
    SH_Max_sites <- max(summary.PdGn.HARMONISED$n)
    
    
    # 6.b Dintra/Dinter HARMONISED -----------------------------------------------
    
    
    ### 6.b Din/Dex
    ### Din is average distance from each observation to the centroid of its cluster
    ### Dex is the average distance between cluster centroids
    ### Calculate with functions from "5.PedogenonModeling_helper.R"
    
    # centroids.test <- centroids.SoilVars.Pedogenon.fun(df.soil =Soil.df.HARMONISED,
    #                                                    depth.var = "Layer_depth",
    #                                                    target.vars = target.vars.HARMONISED  )
    
    ### DEBUG THIS FUNCTION
    HARMONISED.DF.intra <- Dintra.function(df.soil = Soil.df.HARMONISED,
                                        uniqueID = "newID",
                                        depth.var = "Layer_depth",
                                        target.vars = target.vars.HARMONISED)
    
    ### Average Distance between each observation to their centroid.
    Dintra.HARMONISED <- mean(HARMONISED.DF.intra$dist_to_centroid, na.rm=TRUE) 
    
    
    HARMONISED.DF.inter <- Dinter.function(df.soil = Soil.df.HARMONISED,
                                        uniqueID = "newID",
                                        depth.var = "Layer_depth",
                                        target.vars = target.vars.HARMONISED)
    
    ### Calculate average distance between centroids.
    ### These are the distances between centroids.
    Dinter.HARMONISED <- mean(HARMONISED.DF.inter, na.rm=TRUE)
    
    ### Ratio Din to Dex
    Din_Dex_HARMONISED <- Dintra.HARMONISED/Dinter.HARMONISED
    
    ### Store the results of the soil profile distances indices
    out.indices.SoilHarmonised[,k+1] <- c(round(SH_N_PdGn_sites, digits=0),
                                   round(SH_Perc_PdGn_sites, digits=1),
                                   round(SH_Min_sites, digits=0),
                                   round(SH_Median_sites, digits=0),
                                   round(SH_Max_sites, digits=0),
                                   round(Din_Dex_HARMONISED, digits=3))  
    gc()
    
  }
  
  ### Now we run through all the K
  ### Store the results in the same dataframe
  indices.j <- rbind(out.clustering.indices,out.indices.BASONET,out.indices.SoilHarmonised)
  
  ### Write csv file in case foreach crashes
  write.csv(indices.j, file = paste0(OutDir,"scorpan_combi_",j,"/cl_indices_comb",j,".csv"))
  
  ### Keep in a list outside the loop
  kquality_combi_list[[j]] <- indices.j
  
  ### Export BASONET dataframe with cluster assignments
  BASONET.dfs[[j]] <- BASONET.scorpan.j
  write.csv(BASONET.scorpan.j, file = paste0(OutDir,"scorpan_combi_",j,"/BASONET_assign_",j,".csv"))
  
  ### Export harmonised soil dataframe with cluster assignments
  HARMONISED.dfs[[j]] <- Soil.df.harmonised.scorpan.j
  write.csv(Soil.df.harmonised.scorpan.j, file = paste0(OutDir,"scorpan_combi_",j,"/HARMONISED_assign_",j,".csv"))
  
  tmpFiles(current = FALSE, orphan = TRUE, old = TRUE, remove = TRUE)
  gc()
  
  indices.j ### We return this
 

}

stopCluster(cl)
tac <- Sys.time()
tac-tic

save.image("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModeling/relief_combi1.RData")


#########################################################################################################################

### 8. Create HEATMAPS of cluster quality for each k and each covariate combination and choose optimal k and covariate combination



### 9. Map pedogenon classes



### 10. Assign soil health indicators (soil condition indications) to each pedogenon class and plot values
### TOC, TOC:clay, EC, pH, bulk density, N, Mg, K, Ca, etc.
### Boxplots by pedogenon
### Boxplots by vegetation class (conifers vs broadleaved, native vs exotic, etc.)
### Diagram phases by pedogenon class? e.g., bulk density vs TOC or TOC/clay


####################################################################################################################
