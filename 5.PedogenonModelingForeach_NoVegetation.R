#############################################################################################################################################
###  Method for optimizing pedogenon classes as soil districts for the Basque Country
###  in the context of the European soil Monitoring Law
###  In this script: Modeling and optimization the number of clusters and covariate selection
###                  EXCLUDING POTENTIAL VEGETATION

### Desired extent: Basque Country
### Resolution: 25m
### CRS: EPSG=25830

###  Author: Mercedes Roman Dobarco
###  Date: 30/04/2024

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
library(Hmisc)
library(corrplot)

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
search_space <- c(2:50)

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
# setwd("C:/Covariates/Euskadi/Organisms/")
# vegetation.vars <- c("MVG_PC1.tif","MVG_PC2.tif","MVG_PC3.tif","MVG_PC4.tif","MVG_PC5.tif")
# vegetation.r <- rast(vegetation.vars)
# names(vegetation.r) <- paste0("o_",names(vegetation.r))

### Parent material: 8 or 50 % of variance
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/1.CovariatesEus/")
parentmaterial.vars <- c("pm_mca_1.tif","pm_mca_2.tif","pm_mca_3.tif",
                         "pm_mca_4.tif","pm_mca_5.tif","pm_mca_6.tif",
                         "pm_mca_7.tif","pm_mca_8.tif")
parentmaterial.r <- rast(parentmaterial.vars)

### All variables
scorpan <- c(clim.r, relief.r, parentmaterial.r)

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

# ### Split by 3,5,10 window
# relief.potential.3 <-  c("r_slope",
#                          "r_northerness","r_easterness","r_north_slope",
#                          "r_twi","r_mrrtf","r_mrvbf",
#                          "r_valley_depth","r_st_height","r_mid_slope",
#                          "r_norm_height","r_slope_height",
#                          "r_p_curv","r_t_curv","r_tpi_8_3")
# 
# relief.potential.5 <-  c("r_slope_5",
#                          "r_northerness","r_easterness","r_north_slope",
#                          "r_twi","r_mrrtf","r_mrvbf",
#                          "r_valley_depth","r_st_height","r_mid_slope",
#                          "r_norm_height","r_slope_height",
#                          "r_pr_curv_5","r_pl_curv_5","r_lg_curv_5","r_cs_curv_5",
#                          "r_tpi_8_3")
# 
# relief.potential.10 <-  c("r_slope_10",
#                          "r_northerness","r_easterness","r_north_slope",
#                          "r_twi","r_mrrtf","r_mrvbf",
#                          "r_valley_depth","r_st_height","r_mid_slope",
#                          "r_norm_height","r_slope_height",
#                          "r_pr_curv_10","r_pl_curv_10","r_lg_curv_10","r_cs_curv_10",
#                          "r_tpi_8_3")
  
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

### fixed variables  ### I eliminate vegetation
fixed.columns <- c(names(clim.r),"r_dem",names(parentmaterial.r))

### directory to store the models
OutDir <- "C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModelingNoVegetation/"

### Delete temporal files
tmpFiles(current = FALSE,orphan = TRUE,old = TRUE,remove = TRUE)

# ### 3. For each combination of relief covariates: -----------------------

### some did not finished running...
#relief.combi <- relief.combi[c(27:38),]


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
cl <- makeCluster(20)   ###
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
    

# 5.b Dintra/Dinter BASONET -----------------------------------------------

    
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
                                  
    
# ### 6. Predict cluster assignment to soil profiles HARMONISED DATASET--------
    
    ### 6. Predict cluster assignment to soil properties observations - HARMONISED dataset:
    
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

### Does not exist
### save.image("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModelingNoVegetation/relief_noVeg_combi27032024.RData")


# 8. Create HEATMAPS of cluster quality for each k and each covari --------

### 8. Create HEATMAPS of cluster quality for each k and each covariate combination and choose optimal k and covariate combination

### 26/04/2024
### load("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModelingNoVegetation/relief_noVeg_combi27032024.RData")

OutDir <- "C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModelingNoVegetation/"
setwd(OutDir)

### Load the tables with clustering indices results in a for loop
for(j in 1:nrow(relief.combi)) { 
  
  if(j == 1) {
    
  setwd(paste0(OutDir,"scorpan_combi_",j))
  
  table.1 <- read.csv(paste0("cl_indices_comb",j,".csv"))
  
  table.1$scorpan_comb <- j
  
  } else if (j >1) {
    
    setwd(paste0(OutDir,"scorpan_combi_",j))
    
    table.j <- read.csv(paste0("cl_indices_comb",j,".csv"))
    
    table.j$scorpan_comb <- j
    
    table.1 <- rbind(table.1,table.j)
    
  }
  
  
}

dim(table.1)

### Rename table
table_indices <- table.1
rm(table.1, table.j,j)

### Make heatmaps in for loop
indices <- unique(table_indices$Index)

dir.create("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModelingNoVegetation/Plots")
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModelingNoVegetation/Plots/")


for(index in 1:length(indices)){
  
  ### Select index
  index.table <- table_indices[table_indices$Index == indices[[index]], ]
  
  ### Columns to rows
  index.table.long <- tidyr::pivot_longer(index.table, 
                                          cols =  starts_with("K."), 
                                          names_to = "clusters",
                                          values_to = "index")
  
  index.table.long$clusters <- factor(index.table.long$clusters,
                                      levels=unique(index.table.long$clusters))
  
    plot.index <- ggplot(index.table.long,
                       aes(x =scorpan_comb, y = clusters, fill = index)) +
    geom_tile() +
    scale_fill_viridis() +
    labs(title = paste0(indices[[index]]))+
    labs(x="SCORPAN combination")+
    labs(y="Number of clusters")+
    theme(
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 16),
      axis.text.y = element_text(size = 8),
      legend.text=element_text(size=16), #change font size of legend text
      legend.title=element_text(size=16),
        plot.title=element_text(size=20))
  
  file_name = paste("index_plot_",indices[[index]] , ".jpeg", sep="")
  jpeg(file_name, width = 1200, height = 900, units = "px")
  print(plot.index)
  dev.off()
  
}


### 9. Preferred options ---------------------------------------------------

### I am searching for the combination that:
### - has at least 3 observations per pedogenon class
### - has the smallest Din/Dex with BASONET and the complete dataset
### - has a large Calinski-Harbasz index

index.table <- table_indices[table_indices$Index %in% 
                               c("Calinski-Harbasz","B_Min_sites",
                                 "SH_Min_sites","B_Din_Dex","SH_Din_Dex"), ]
index.table.long <- tidyr::pivot_longer(index.table, 
                                        cols =  starts_with("K."), 
                                        names_to = "clusters",
                                        values_to = "index")

index.table.wide <- tidyr::pivot_wider(index.table.long, 
                                       id_cols=c(scorpan_comb ,clusters ),
                                       names_from = "Index",
                                       values_from = "index")
colnames(index.table.wide)<- c("scorpan_comb","clusters","Calinski_Harbasz",
                               "B_Min_sites","B_Din_Dex","SH_Min_sites","SH_Din_Dex")

### At least 3 onservations per pedogenon class
index.table.wide <- index.table.wide[index.table.wide$B_Min_sites >= 3,]
#index.table.wide <- index.table.wide[index.table.wide$scorpan_comb != 38,] ### This one has a weird result 
### The Din/Dex smaller than 1 (at least the same dispersion, no more!)
index.table.wide <- index.table.wide[index.table.wide$B_Din_Dex <=1,]
### Arrange
index.table.wide <- arrange(index.table.wide, B_Din_Dex, desc(Calinski_Harbasz))
write.csv(index.table.wide, file="CH_BDinDex_Min3.csv")
index.table.wide[1:10,]

index.table.wide$clusters <- factor(index.table.wide$clusters,
                                    levels=paste0("K.",c(4:10)))
index.table.wide <- arrange(index.table.wide, B_Din_Dex, desc(Calinski_Harbasz))

ggplot(index.table.wide,
       aes(x =scorpan_comb, y = clusters, fill = B_Din_Dex)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Din/Dex BASONET")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))

ggplot(index.table.wide,
       aes(x =scorpan_comb, y = clusters, fill = Calinski_Harbasz)) +
  geom_tile() +
  scale_fill_viridis(direction=-1) +
  labs(title = "Calinski-Harbasz BASONET")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))


#sort(unique(index.table.wide$scorpan_comb))
index.table.wide$scorpan_comb_f <- factor(index.table.wide$scorpan_comb,
                                          levels = sort(unique(index.table.wide$scorpan_comb)))
ggplot(index.table.wide,
       aes(x = B_Din_Dex, y = Calinski_Harbasz, 
           color=scorpan_comb_f, shape=clusters)) +
  geom_point(size=4) +
  scale_shape_manual(values = c('K.4'= 0,'K.5'= 8,'K.6'=15, 'K.7'=18, 'K.8'=16,
                                'K.9'=17,'K.10'=1)) +
  # # scale_color_viridis(discrete=TRUE)+
  labs(title = "Calinski-Harbasz index and Din/Dex")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=14))









top <- data.s[data.s$bi_class %in% c("4-4"),]
index.table.wide$DinDexNeg <- index.table.wide$B_Din_Dex * (-1)

data <- bi_class(top, x = "Calinski_Harbasz", y = "DinDexNeg", 
                 style = "equal", dim = 3)
data <- as.data.frame(data)
write.csv(data, file="biscale_Basonet_equal.csv")

ggplot(data,
       aes(x =scorpan_comb, y = clusters, fill = bi_class)) +
  geom_tile() +
  bi_scale_fill(pal = "DkBlue2", dim = 3) +
  labs(title = "Calinski-Harbasz index and Din/Dex")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=14))



# sample 4x4 legend
legend <- bi_legend(pal = "DkBlue2",
                    dim = 3,
                    xlab = "Calinski-Harbasz",
                    ylab = "- Din/Dex",
                    size = 16)
bi_pal(pal="DkBlue2", dim = 4, preview = TRUE, flip_axes = FALSE, rotate_pal = FALSE)







# create classes
index.table.wide$DinDexNeg <- index.table.wide$B_Din_Dex * (-1)
data <- bi_class(index.table.wide, x = "Calinski_Harbasz", y = "DinDexNeg", 
                 style = "quantile", dim = 4)

hist(data$Calinski_Harbasz,breaks=30)
hist(data$B_Din_Dex,breaks=30)

ggplot(data,
       aes(x =scorpan_comb, y = clusters, fill = bi_class)) +
  geom_tile() +
  bi_scale_fill(pal = "DkBlue2", dim = 4) +
  labs(title = "Calinski-Harbasz index and Din/Dex")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=14))

data[data$bi_class =="4-4",]
# scorpan_comb clusters Calinski_Harbasz B_Min_sites B_Din_Dex SH_Min_sites SH_Din_Dex DinDexNeg bi_class
#         <int> <fct>               <dbl>        <dbl>     <dbl>        <dbl>      <dbl>    <dbl> <chr>   
# 1           16 K.9                 1344.           4     0.795            6       1.10    -0.795 4-4     
# 2           11 K.10                1350.           4     0.83             6       1.13    -0.83  4-4     
# 3            5 K.10                1358.           4     0.836            6       1.19    -0.836 4-4     
# 4           20 K.11                1322.           3     0.843            6       1.08    -0.843 4-4     
# 5           22 K.11                1377.           4     0.861            6       1.20    -0.861 4-4     
# 6            7 K.12                1334.           4     0.864            6       1.18    -0.864 4-4   



























































































































### Zoom to optimal number of classes
index.table <- table_indices[table_indices$Index == "Calinski-Harbasz", ]
index.table.long <- tidyr::pivot_longer(index.table, 
                                        cols =  starts_with("K."), 
                                        names_to = "clusters",
                                        values_to = "index")

index.table.long <- index.table.long[index.table.long$clusters %in% paste0("K.",2:40),]
index.table.long$clusters <- factor(index.table.long$clusters,
                                    levels=paste0("K.",2:40))
ggplot(index.table.long,
       aes(x =scorpan_comb, y = clusters, fill = index)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Calinski-Harbasz")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))

### what combination has the maximum index C-H?
index.table.long <- arrange(index.table.long, desc(index))
index.table.long[1:15,]$scorpan_comb

### what combinations?
relief.combi.df[unique(index.table.long[1:15,]$scorpan_comb),]

### BIC
index.table <- table_indices[table_indices$Index == "BIC", ]
index.table.long <- tidyr::pivot_longer(index.table, 
                                        cols =  starts_with("K."), 
                                        names_to = "clusters",
                                        values_to = "index")

index.table.long <- index.table.long[index.table.long$clusters %in% paste0("K.",2:140),]
index.table.long$clusters <- factor(index.table.long$clusters,
                                    levels=paste0("K.",2:140))
ggplot(index.table.long,
       aes(x =scorpan_comb, y = clusters, fill = index)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Calinski-Harbasz")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))

### what combination has the minimum BIC?
index.table.long <- arrange(index.table.long, index)
index.table.long[1:20,]

### Ratio between to total SS
index.table <- table_indices[table_indices$Index == "BetweenSS to TotalSS", ]
index.table.long <- tidyr::pivot_longer(index.table, 
                                        cols =  starts_with("K."), 
                                        names_to = "clusters",
                                        values_to = "index")

index.table.long <- index.table.long[index.table.long$clusters %in% paste0("K.",2:140),]
index.table.long$clusters <- factor(index.table.long$clusters,
                                    levels=paste0("K.",2:140))

### what combination has the greater ratio?
index.table.long <- arrange(index.table.long, desc(index))
table(index.table.long[1:20,]$scorpan_comb)
# 6 13 15 16 24 25 
# 1  1  4 10  1  3

ggplot(index.table.long,
       aes(x =scorpan_comb, y = clusters, fill = index)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Calinski-Harbasz")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))

### what combination has the minimum BIC?
index.table.long <- arrange(index.table.long, index)
index.table.long[1:20,]

### Minimum number of observations

### Zoom to optimal number of classes
index.table <- table_indices[table_indices$Index == "B_Min_sites", ]
index.table.long <- tidyr::pivot_longer(index.table, 
                                        cols =  starts_with("K."), 
                                        names_to = "clusters",
                                        values_to = "index")

index.table.long <- index.table.long[index.table.long$clusters %in% paste0("K.",2:140),]
## Minimum 3 observations per soil district 
index.table.long <- index.table.long[index.table.long$index >= 3,]

index.table.long$clusters <- factor(index.table.long$clusters,
                                    levels=paste0("K.",2:140))
ggplot(index.table.long,
       aes(x =scorpan_comb, y = clusters, fill = index)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Minimum number of soil observations per district")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))


### Zoom to optimal number of classes
index.table <- table_indices[table_indices$Index == "SH_Min_sites", ]
index.table.long <- tidyr::pivot_longer(index.table, 
                                        cols =  starts_with("K."), 
                                        names_to = "clusters",
                                        values_to = "index")

index.table.long <- index.table.long[index.table.long$clusters %in% paste0("K.",2:140),]
## Minimum 3 observations per soil district 
index.table.long <- index.table.long[index.table.long$index >= 3,]

index.table.long$clusters <- factor(index.table.long$clusters,
                                    levels=paste0("K.",2:140))

ggplot(index.table.long,
       aes(x =scorpan_comb, y = clusters, fill = index)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Minimum number of soil observations per district")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))

### what combination has max number of minimum observations for a greater number of clusters?
index.table <- table_indices[table_indices$Index %in% c("B_Min_sites","Calinski-Harbasz"), ]
index.table.long <- tidyr::pivot_longer(index.table, 
                                        cols =  starts_with("K."), 
                                        names_to = "clusters",
                                        values_to = "index")

index.table.wide <- tidyr::pivot_wider(index.table.long, 
                                       id_cols=c(scorpan_comb ,clusters ),
                                        names_from = "Index",
                                        values_from = "index")

#index.table.wide <- index.table.wide[index.table.wide$B_Min_sites >= 3,]
index.table.wide <- arrange(index.table.wide, desc(`Calinski-Harbasz`))
write.csv(index.table.wide, file="CH_MinB.csv")

index.table.wide$clusters <- factor(index.table.wide$clusters,
                                    levels=paste0("K.",2:140))

ggplot(index.table.wide[index.table.wide$B_Min_sites >= 3,],
       aes(x =scorpan_comb, y = clusters, fill = `Calinski-Harbasz`)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Calinski-Harbasz index and minimum number of soil observations per")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))

ggplot(index.table.wide[index.table.wide$B_Min_sites >= 3,],
       aes(x =scorpan_comb, y = clusters, fill = B_Min_sites)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Calinski-Harbasz index of classes with at least 3soil profiles")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))

### I use a biscale ramp to plot at the same time minimum number of classes and C-H index
install.packages("biscale")
library(biscale)
# create classes
colnames(index.table.wide) <- c("scorpan_comb","clusters","Calinski-Harbasz", "B_Min_sites")
index.table.wide <- index.table.wide[index.table.wide$clusters %in% c(paste0("K.",2:40)), ]
index.table.wide$clusters <- factor(index.table.wide$clusters,
                                    levels=paste0("K.",2:40))

data <- bi_class(index.table.wide, x = "Calinski-Harbasz", y = "B_Min_sites", style = "equal", dim = 4)
data <- as.data.frame(data)

ggplot(data,
       aes(x =scorpan_comb, y = clusters, fill = bi_class)) +
  geom_tile() +
  bi_scale_fill(pal = "DkBlue2", dim = 4) +
  labs(title = "Calinski-Harbasz index and minimum number of soil observations per class")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=14))


# sample 4x4 legend
legend <- bi_legend(pal = "DkBlue2",
                    dim = 4,
                    xlab = "Calinski-Harbasz ",
                    ylab = "Min number soil observations per class",
                    size = 16)
bi_pal(pal, dim = 3, preview = TRUE, flip_axes = FALSE, rotate_pal = FALSE)

#### Din/Dex BASONET
index.table <- table_indices[table_indices$Index == "B_Din_Dex", ]
index.table.long <- tidyr::pivot_longer(index.table, 
                                        cols =  starts_with("K."), 
                                        names_to = "clusters",
                                        values_to = "index")
### Exclude one combination that has too high index
index.table.long <- index.table.long[index.table.long$scorpan_comb != 38,]
index.table.long$clusters <- factor(index.table.long$clusters,
                                    levels=paste0("K.",2:140))
ggplot(index.table.long,
       aes(x =scorpan_comb, y = clusters, fill = index)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Din / Dex")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))

### what combination has the maximum index C-H?
index.table.long <- arrange(index.table.long, index)
table(index.table.long[1:20,]$scorpan_comb)

#### Din/Dex Harmonised dataset
index.table <- table_indices[table_indices$Index == "SH_Din_Dex", ]
index.table.long <- tidyr::pivot_longer(index.table, 
                                        cols =  starts_with("K."), 
                                        names_to = "clusters",
                                        values_to = "index")
### Exclude one combination that has too high index
index.table.long <- index.table.long[index.table.long$scorpan_comb != 38,]
index.table.long$clusters <- factor(index.table.long$clusters,
                                    levels=paste0("K.",2:140))
ggplot(index.table.long,
       aes(x =scorpan_comb, y = clusters, fill = index)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Din / Dex")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))

### what combination has the maximum index C-H?
index.table.long <- arrange(index.table.long, index)
table(index.table.long[1:20,]$scorpan_comb)

### For those around 20 clusters

#### Din/Dex BASONET
index.table <- table_indices[table_indices$Index == "B_Din_Dex", ]
index.table.long <- tidyr::pivot_longer(index.table, 
                                        cols =  starts_with("K."), 
                                        names_to = "clusters",
                                        values_to = "index")
### Exclude one combination that has too high index
index.table.long <- index.table.long[index.table.long$scorpan_comb != 38,]
index.table.long <- index.table.long[index.table.long$clusters %in% paste0("K.", 2:20),]
index.table.long$clusters <- factor(index.table.long$clusters,
                                    levels=paste0("K.",2:20))
ggplot(index.table.long,
       aes(x =scorpan_comb, y = clusters, fill = index)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Din / Dex")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))

### what combination has the maximum index C-H?
index.table.long <- arrange(index.table.long, index)
table(index.table.long[1:20,]$scorpan_comb)

index.table <- table_indices[table_indices$Index %in% 
                               c("B_Min_sites","Calinski-Harbasz","B_Din_Dex"), ]
index.table.long <- tidyr::pivot_longer(index.table, 
                                        cols =  starts_with("K."), 
                                        names_to = "clusters",
                                        values_to = "index")

index.table.wide <- tidyr::pivot_wider(index.table.long, 
                                       id_cols=c(scorpan_comb ,clusters ),
                                       names_from = "Index",
                                       values_from = "index")
colnames(index.table.wide)<- c("scorpan_comb","clusters","Calinski_Harbasz",
                               "B_Min_sites","B_Din_Dex")

index.table.wide <- index.table.wide[index.table.wide$B_Min_sites >= 3,]
index.table.wide <- index.table.wide[index.table.wide$scorpan_comb != 38,]
index.table.wide <- arrange(index.table.wide, desc(Calinski_Harbasz),B_Din_Dex)
write.csv(index.table.wide, file="Basonet_CH_MinO_DinDex.csv")

index.table.wide$clusters <- factor(index.table.wide$clusters,
                                    levels=paste0("K.",2:140))

ggplot(index.table.wide,
       aes(x =scorpan_comb, y = clusters, fill = Calinski_Harbasz)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Calinski-Harbasz index")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))

ggplot(index.table.wide,
       aes(x =scorpan_comb, y = clusters, fill = B_Din_Dex)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Din / Dex")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))

ggplot(index.table.wide,
       aes(x =scorpan_comb, y = clusters, fill = B_Min_sites)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Calinski-Harbasz index of classes with at least 3soil profiles")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))

### I use a biscale ramp to plot at the same time minimum number of classes and C-H index
install.packages("biscale")
library(biscale)
# create classes
index.table.wide$DinDexNeg <- index.table.wide$B_Din_Dex * (-1)
data <- bi_class(index.table.wide, x = "Calinski_Harbasz", y = "DinDexNeg", 
                 style = "equal", dim = 4)
data <- as.data.frame(data)
write.csv(data, file="biscale_Basonet_equal.csv")

ggplot(data,
       aes(x =scorpan_comb, y = clusters, fill = bi_class)) +
  geom_tile() +
  bi_scale_fill(pal = "DkBlue2", dim = 4) +
  labs(title = "Calinski-Harbasz index and Din/Dex")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=14))

# sample 4x4 legend
legend <- bi_legend(pal = "DkBlue2",
                    dim = 4,
                    xlab = "Calinski-Harbasz",
                    ylab = "- Din/Dex",
                    size = 16)
bi_pal(pal, dim = 3, preview = TRUE, flip_axes = FALSE, rotate_pal = FALSE)

### the combination that:
### - has at least 3 observations per pedogenon class
### - has the smallest Din/Dex
### - has a large Calinski-Harbasz index
data <- arrange(data, B_Din_Dex)
data[1:20,]

### Now with the harmonised dataset
index.table <- table_indices[table_indices$Index %in% 
                               c("SH_Min_sites","Calinski-Harbasz","SH_Din_Dex"), ]
index.table.long <- tidyr::pivot_longer(index.table, 
                                        cols =  starts_with("K."), 
                                        names_to = "clusters",
                                        values_to = "index")

index.table.wide <- tidyr::pivot_wider(index.table.long, 
                                       id_cols=c(scorpan_comb ,clusters ),
                                       names_from = "Index",
                                       values_from = "index")
colnames(index.table.wide)<- c("scorpan_comb","clusters","Calinski_Harbasz",
                               "SH_Min_sites","SH_Din_Dex")

index.table.wide <- index.table.wide[index.table.wide$SH_Min_sites >= 3,]
index.table.wide <- index.table.wide[index.table.wide$scorpan_comb != 38,]
index.table.wide <- arrange(index.table.wide, desc(Calinski_Harbasz),SH_Din_Dex)
write.csv(index.table.wide, file="SH_CH_MinObs_DinDex.csv")

index.table.wide$clusters <- factor(index.table.wide$clusters,
                                    levels=paste0("K.",2:140))

ggplot(index.table.wide,
       aes(x =scorpan_comb, y = clusters, fill = Calinski_Harbasz)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Calinski-Harbasz index")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))

ggplot(index.table.wide,
       aes(x =scorpan_comb, y = clusters, fill = SH_Din_Dex)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Din / Dex")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))

### I use a biscale ramp to plot at the same time minimum number of classes and C-H index
install.packages("biscale")
library(biscale)
# create classes
index.table.wide$DinDexNeg <- index.table.wide$SH_Din_Dex * (-1)
data <- bi_class(index.table.wide, x = "Calinski_Harbasz", y = "DinDexNeg", 
                 style = "quantile", dim = 4)
data <- arrange(data, SH_Din_Dex)
data <- as.data.frame(data)
data[data$bi_class =="3-3",]
write.csv(data, file="biscale_SH_equal.csv")

data.s <- data[data$bi_class %in% c("3-4","4-3","4-4","3-3"),]
data.s <- arrange(data.s, SH_Din_Dex)
data.s[data.s$bi_class =="4-3",]
data.s[data.s$bi_class =="3-3",]
data.s[1:20,]
data.s <- arrange(data.s, desc(Calinski_Harbasz))
data.s[1:20,]


hist(data$Calinski_Harbasz,breaks=30)
hist(data$SH_Din_Dex,breaks=30)

ggplot(data,
       aes(x =scorpan_comb, y = clusters, fill = bi_class)) +
  geom_tile() +
  bi_scale_fill(pal = "DkBlue2", dim = 4) +
  labs(title = "Calinski-Harbasz index and Din/Dex")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=14))

# sample 4x4 legend
legend <- bi_legend(pal = "DkBlue2",
                    dim = 4,
                    xlab = "Calinski-Harbasz",
                    ylab = "- Din/Dex",
                    size = 16)
a <- bi_pal("DkBlue2", dim = 4)
data <- arrange(data, B_Din_Dex)
data[1:20,]


### 9. Preferred options ---------------------------------------------------

### I am searching for the combination that:
### - has at least 3 observations per pedogenon class
### - has the smallest Din/Dex with BASONET and the complete dataset
### - has a large Calinski-Harbasz index

### Now both datasets
index.table <- table_indices[table_indices$Index %in% 
                               c("Calinski-Harbasz","B_Min_sites",
                                 "SH_Min_sites","B_Din_Dex","SH_Din_Dex"), ]
index.table.long <- tidyr::pivot_longer(index.table, 
                                        cols =  starts_with("K."), 
                                        names_to = "clusters",
                                        values_to = "index")

index.table.wide <- tidyr::pivot_wider(index.table.long, 
                                       id_cols=c(scorpan_comb ,clusters ),
                                       names_from = "Index",
                                       values_from = "index")
colnames(index.table.wide)<- c("scorpan_comb","clusters","Calinski_Harbasz",
                               "B_Min_sites","B_Din_Dex","SH_Min_sites","SH_Din_Dex")

### At least 3 onservations per pedogenon class
index.table.wide <- index.table.wide[index.table.wide$B_Min_sites >= 3,]
index.table.wide <- index.table.wide[index.table.wide$scorpan_comb != 38,] ### This one has a weird result 
### The Din/Dex smaller than 1 (at least the same dispersion, no more!)
index.table.wide <- index.table.wide[index.table.wide$B_Din_Dex <=1,]
### Arrange
index.table.wide <- arrange(index.table.wide, B_Din_Dex, desc(Calinski_Harbasz))
### write.csv(index.table.wide, file="SH_CH_MinObs_DinDex.csv")

index.table.wide$clusters <- factor(index.table.wide$clusters,
                                    levels=paste0("K.",2:140))
index.table.wide <- arrange(index.table.wide, B_Din_Dex, desc(Calinski_Harbasz))

relief.combi.df[c(16,32),]

ggplot(index.table.wide,
       aes(x =scorpan_comb, y = clusters, fill = SH_Din_Dex)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(title = "Din/Dex harmonised dataset")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=20))


# create classes
index.table.wide$DinDexNeg <- index.table.wide$B_Din_Dex * (-1)
data <- bi_class(index.table.wide, x = "Calinski_Harbasz", y = "DinDexNeg", 
                 style = "quantile", dim = 4)

hist(data$Calinski_Harbasz,breaks=30)
hist(data$B_Din_Dex,breaks=30)

ggplot(data,
       aes(x =scorpan_comb, y = clusters, fill = bi_class)) +
  geom_tile() +
  bi_scale_fill(pal = "DkBlue2", dim = 4) +
  labs(title = "Calinski-Harbasz index and Din/Dex")+
  labs(x="SCORPAN combination")+
  labs(y="Number of clusters")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    legend.text=element_text(size=16), 
    legend.title=element_text(size=16),
    plot.title=element_text(size=14))

data[data$bi_class =="4-4",]
# scorpan_comb clusters Calinski_Harbasz B_Min_sites B_Din_Dex SH_Min_sites SH_Din_Dex DinDexNeg bi_class
#         <int> <fct>               <dbl>        <dbl>     <dbl>        <dbl>      <dbl>    <dbl> <chr>   
# 1           16 K.9                 1344.           4     0.795            6       1.10    -0.795 4-4     
# 2           11 K.10                1350.           4     0.83             6       1.13    -0.83  4-4     
# 3            5 K.10                1358.           4     0.836            6       1.19    -0.836 4-4     
# 4           20 K.11                1322.           3     0.843            6       1.08    -0.843 4-4     
# 5           22 K.11                1377.           4     0.861            6       1.20    -0.861 4-4     
# 6            7 K.12                1334.           4     0.864            6       1.18    -0.864 4-4   


# ### 10. Map pedogenon classes -------------------------------------------

### Number of variables per soil-forming factor:

### Climate = 5 variables selected from correlation plots 
### and variable importance for predicting soil properties
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/4.CovariateSelection/scaled/")
clim.vars <- c("clim_bio1.tif","clim_bio4.tif","clim_bio5.tif","clim_bio12.tif","clim_bio15.tif")
clim.r <- rast(clim.vars)
plot(clim.r)

### Relief: between 5 to 10 - These will be selected with iterations
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/4.CovariateSelection/scaled/")
relief.fixed <- c("relief_dem.tif")

### Option combi 16
relief.combi16 <- c("relief_slope.tif", # slope
                    "relief_twi.tif", # TWI 
                    "relief_valley_depth.tif",
                    "relief_st_height.tif") # Preferred hydrological variables
relief.r <- rast(c(relief.fixed, relief.combi16))
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

### Take sample
### Regular sample
set.seed(2233)
scorpanSample <- terra::spatSample(x = scorpan, 
                                   size=200000,
                                   method="regular", 
                                   as.df=TRUE, 
                                   xy=TRUE)

### complete cases
scorpanSample <- scorpanSample[complete.cases(scorpanSample),]

## Take out the coordinates
CLORPT.df.coords <- scorpanSample[,1:2]
CLORPT.df <- scorpanSample[,3:ncol(scorpanSample)]

### Perform Cholesky transformation to decorrelate the data
# The basic Euclidean distance treats each variable as equally important in calculating the distance.
# An alternative approach is to scale the contribution of individual variables to the distance value according
# to the variability of each variable. This approach is illustrated by the Mahalanobis distance, 
# which is a measure of the distance between each observation in a multidimensional cloud of points and
# the centroid of the cloud.
### Calculate the Mahalanobis distance, as Euclidean distance after applying the Cholesky decomposition

# # Rescale the data
C <- chol(var(as.matrix(CLORPT.df)))
CLORPT.rs <- as.matrix(CLORPT.df) %*% solve(C)

# perform KMeans_rcpp clustering with the selected number of clusters
library(ClusterR)

my_seed <- 4587 ### 
myk <- 9
pedogenons.eus.comb16.k9 <-KMeans_rcpp(data=CLORPT.rs, 
                                  clusters=myk,
                                  num_init = 20, 
                                  max_iters = 20000,
                                  initializer = "kmeans++",
                                  fuzzy = FALSE, 
                                  verbose = TRUE,
                                  seed = my_seed)

#  Mapping the clusters and write layers
### Previously, the stack with the scaled covariates were at
scorpan
plot(scorpan)

### Extract the index of the centroids that are na/nan/Inf
Kcent.nan <- which(apply(pedogenons.eus.comb16.k9[["centroids"]], MARGIN = 1, FUN = function(y){any(is.na(y))}))
### None are NA

### Define the size of the blocks --- At each raster row do we start and finish each crop?
bs <- blocks(scorpan)
### I crop all raster files (across variables stacks) in a parallel process
### with a %dopar% from the foreach package
### this is thought for a large study area, but of course for the example is not needed

## My desired cluster number
K <- 9

k_rast_list <- list()

for(i in 1:bs$n){
  
  ### Get one tile of the raster stack
  # tile <- crop(covariates.stack, 
  #                ext(covariates.stack, bs$row[[i]], bs$row[[i]]+bs$nrows[[i]], 1, ncol(covariates.stack)))
  ### Crop with this syntax from https://stackoverflow.com/questions/69965770/subset-a-raster-using-row-column-index-in-terra 
  if(i < bs$n)  {
    tile <- scorpan[ bs$row[[i]]:(bs$row[[i]]+bs$nrows[[i]]), 1:ncol(scorpan), drop=FALSE]
  }  else {
    tile <- scorpan[ bs$row[[i]]:(bs$row[[i]]+bs$nrows[[i]]-5), 1:ncol(scorpan), drop=FALSE]
  }
  
  ### Transform into a dataframe
  tile.df <- as.data.frame(tile, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=FALSE)
  
  ### For each new pixel, I first have to rescale its values
  ## Take out the coordinates
  tile.df.coords <- tile.df[,1:2]
  tile.df <- tile.df[,3:ncol(tile.df)]
  
  # Rescale the data with CLORPT.df (sample from the stack covariates that was used to calibrate the kmeans in the first place)
  tile.df.rs <- as.matrix(tile.df) %*% solve(C)
  tile.df.rs <- as.data.frame(tile.df.rs)
  
  ### Predict cluster assignment
  
  ### Extract the index of the dataframe rows that are na/nan/Inf
  df.na <- which(apply(tile.df.rs, MARGIN = 1, FUN = function(x) {any(is.na(x))}))
  
  ### Create empty prediction column
  tile.df.rs$cluster <- NA
  
  ### If K is more than one instead of km.pedogenon.rcpp being a kmeans model, it would be a list of models
  ### km.pedogenon.rcpp[[m]]$centroids
  ### predict in those rows where there are not na
  tile.df.rs[-df.na, ]$cluster  <- predict_KMeans(data = tile.df.rs[-df.na,1:(ncol(tile.df.rs)-1)], 
                                                  CENTROIDS = pedogenons.eus.comb16.k9$centroids)
  ### Assign the values to a new raster
  k.pred <- setValues(tile[[1]], tile.df.rs$cluster)
  names(k.pred) <- "PdGn"
  k_rast_list[[i]] <- k.pred # Return this
}

#stopCluster(cl)
## Assign function to mosaic
k_sprc <- sprc(k_rast_list)
## Create mosaic for whole France
k.raster <- terra::mosaic(k_sprc, 
                          fun="min", 
                          overwrite=TRUE, 
                          filename="C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModeling/PdGnCmb16K9.tif")
plot(k.raster)
k.raster <- terra::rast("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModeling/PdGnCmb16K9.tif")
plot(k.raster)

setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModeling/")


### Play with map output

### Area of pedogenon classes
pedogenon.area.func <- function(kmap, fname) {
  
  areaPixels <-terra::cellSize(kmap)
  s <- c(kmap, areaPixels)
  k.A <-values(s)
  k.A <- as.data.frame(k.A)
  colnames(k.A) <- c("Pedogenon", "Area_m2")
  k.A <- k.A[!is.na(k.A$Pedogenon),]
  
  k.area.df <-  k.A %>% 
    group_by(.,as.factor(Pedogenon), .drop=TRUE ) %>% ## Group by Pedogenon
    summarise(Pedogenon_area = sum(Area_m2, na.rm=TRUE)) ### sum area by Pedogenon class
  k.area.df <- as.data.frame(k.area.df)
  colnames(k.area.df) <- c("Pedogenon", "Area_m2")
  k.area.df$Area_Km2 <- round(k.area.df$Area_m2 /1000000, digits=1)
  k.area.df <- k.area.df[,c("Pedogenon", "Area_Km2")]
  write.csv(k.area.df, file=paste0(fname,".csv")) ### Write table to csv file
  return(k.area.df) ## and return
}

PdGn.area <- pedogenon.area.func(kmap=k.raster, fname="PdGnAreaCmb16K9")

### Returns a table with a row per Pedogenon indicating the closer Pedogenon class,
### the Mahalanobis distance between Pedogenons calculated with CLORPT covariates,
### and the areas that they occupy in NSW, or the study area.
### Note: the distance is calculated with the Euclidean method, but since the data of the
### training dataset was rescaled with the inverse Cholesky transformation, 
### the resulting distance is the same as the Mahalanobis distance calculated on the original data
### Inputs:
### kmodel - kmeans model from the package ClusterR
### k.area.df - is the output of the pedogenon.area.func function
### fname - name to export the table to csv
centroid.dist.func <- function(kmodel, k.area.df, fname){
  
  ### kmodel is a kmeans model
  ### k.area.df is the output of the Pedogenon.area.func function
  
  # extract the centroids
  K.centroids <- kmodel$centroids
  K.centroids <- as.data.frame(K.centroids)
  K.centroids$Pedogenon <- c(1:nrow(K.centroids))
  #rownames(K.centroids) <- c(1:nrow(K.centroids))
  
  ## Is any centroid NA?
  ### Extract the index of the centroids that are na/nan/Inf
  Kcent.nan <- which(apply(K.centroids, MARGIN = 1, FUN = function(x) {any(is.na(x))}))
  
  ### Calculate distance between all centroids
  dist.centroids <- dist(x=K.centroids[,!names(K.centroids) %in% c("Pedogenon")],
                         method = "euclidean")
  
  ### Create empty dataframe to store output
  outs <- data.frame(Pedogenon=rep(as.integer(NA),nrow(K.centroids)),
                     ClosestP=rep(as.integer(NA),nrow(K.centroids)),
                     Distance=rep(as.double(NA), nrow(K.centroids)))
  
  outs$Pedogenon <- K.centroids$Pedogenon ### Assign Pedogenon
  Gs <- as.numeric(as.character(outs$Pedogenon))
  dist.centroids <- as.matrix(dist.centroids)
  ### Calculate distance to the closest Pedogenon
  for(i in 1:nrow(outs)){
    min.dist <- sort(dist.centroids[rownames(dist.centroids)[Gs[[i]]],])[2]
    outs[i,"Distance"] <- min.dist
    outs[i,"ClosestP"] <- names(min.dist)
  }
  ### Remember that those Pedogenons that don't exist are NA
  outs$ClosestP <- ifelse(outs$Pedogenon %in% Kcent.nan, NA, outs$ClosestP )
  outs$Distance <- ifelse(outs$Pedogenon %in% Kcent.nan, NA, outs$Distance )
  colnames(outs) <- c("Pedogenon", "Closest Pedogenon", "Distance")
  outs$Distance <-round(outs$Distance, digits = 3)
  
  ### Join with the Pedogenon area
  outs$Pedogenon <- as.character(outs$Pedogenon)
  k.area.df$Pedogenon <- as.character(k.area.df$Pedogenon)
  outs <- left_join(outs, k.area.df, by ="Pedogenon")
  
  ### Create column with area of the closest Pedogenon
  outs$Pedo2.Area <- NA
  if(length(Kcent.nan) >0) {
    G.exists <- c(1:nrow(outs))[-Kcent.nan]
  } else if(length(Kcent.nan) == 0) {
    G.exists <- c(1:nrow(outs))
  }
  
  for(i in 1:length(G.exists)){
    target.G <- outs[outs$Pedogenon == G.exists[[i]], ]$`Closest Pedogenon`
    target.A <- outs[outs$Pedogenon == target.G, ]$Area_Km2
    outs[outs$Pedogenon == G.exists[[i]], ]$Pedo2.Area <- target.A
  }
  
  colnames(outs) <- c("Pedogenon", "Closest.Pedogenon", "MahabDist", "Area_Km2", "Closests.Pedo.Area_Km2")
  
  write.csv(outs, file=paste0(fname,".csv")) ### Write table to csv file
  return(as.data.frame(outs)) ## and return
}

closetPdGn_k9Cmb16 <- centroid.dist.func(kmodel = pedogenons.eus.comb16.k9,
                   k.area.df = PdGn.area,
                   fname = "closetPdGn_k9Cmb16")


### First, perform the hierarchical clustering and save it to plot
### Input: kmodel - kmeans model from the package ClusterR
### Output: Hierarchical cluster (ward.D2 distance) of pedogenon centroids, hclust object
viz.map.legend.hclust <- function(kmodel) {
  
  ### Extract centroids from model
  centroids <- kmodel$centroids
  ### Extract the index of the centroids that are na/nan/Inf
  Kcent <- as.data.frame(centroids)
  Kcent.nan <- which(apply(Kcent, MARGIN = 1, FUN = function(x) {any(is.na(x))}))
  ### Exclude these clusters from everywhere
  if(length(Kcent.nan) >0) {
    Kcent.exist <- Kcent[-Kcent.nan,]
  } else if(length(Kcent.nan) == 0) {
    Kcent.exist <- Kcent
  }
  # Kcent.exist <- Kcent[-Kcent.nan,]
  ### Hierarchical clustering
  hc <- hclust(dist(Kcent.exist), method="ward.D2")
  plot(dendsort(hc), main="Hierarchical clustering of kmeans centroids", sub="", xlab="")
  return(hc)
  
}

hclust_PdGnk9_Cmb16 <- viz.map.legend.hclust(pedogenons.eus.comb16.k9)

### function to choose the number of branches for color ramps
### Input:
### hc.object - hclust object, hierarchical cluster, output from viz.map.legend.hclust function
### branchN - number of branches that we are considering for this kmeans model
### Output: a plot with the dendrogram and colored branches
viz.branches <- function(hc.object, branchN) {
  hc.object %>% as.dendrogram(.) %>% color_branches(., k = branchN) %>%
    plot(., main = paste0("Colored ",branchN," branches"))
}

viz.branches(hclust_PdGnk9_Cmb16,9)


### Choose my color palette manually:
hcl_palettes(plot = TRUE)

mypalette <- c(sequential_hcl("PurpOr", n = 5),
               sequential_hcl("TealGrn", n = 5),
               sequential_hcl("BurgYl", n = 5),
               sequential_hcl("RdPu", n = 5),
               sequential_hcl("GnBu", n = 5),
               sequential_hcl("OrYel", n = 5)
               )

plot(1:30, 1:30, pch=19, cex=3, col=mypalette)
mypalette <- mypalette[c(2,4,6,8,12,18,22,26,28)]
mypalette <- c("#9D50A6","#177F97", "#2EC6AF","#0090BA","#A74F5A","#EF4868", "#F39B4C","cornsilk3", "#EDAA7D")


viz.map.legend.pal <- function(kmodel, ### k-means model from ClusterR
                               manual, mypalette, ### If manual, provide the color palette (will not know to which class they are assigned)
                               branchN =NULL, pal.names, ### for higher number of classes (e.g., > 15) provide palette names
                               legend.name, ###  
                               kmap, need.proj){
  
  # kmodel = pedogenons.eus.comb16.k9
  # branchN = 8
  # pal.names =  c("PurpOr","TealGrn","OrYel","BurgYl","RdPu",
  #                             "GnBu","YlOrRd","Peach","Turku","Lajolla",
  #                             "OrRd", "Greens", "Burg", "Heat 2", "Dark Mint",
  #                             "Blues", "SunsetDark", "PuBuGn", "Viridis", "Heat")
  # 
  # legend.name = "SoilDistricts"
  # kmap =k.raster
  # need.proj =TRUE
  
  ### Extract centroids from model
  centroids <- kmodel$centroids
  ### Extract the index of the centroids that are na/nan/Inf
  Kcent <- as.data.frame(centroids)
  Kcent$Pedogenon <- c(1:nrow(Kcent))
  Kcent.nan <- which(apply(Kcent, MARGIN = 1, FUN = function(x) {any(is.na(x))}))
  ### Exclude these clusters from everywhere
  if(length(Kcent.nan) >0) {
    Kcent.exist <- Kcent[-Kcent.nan,]
  } else if(length(Kcent.nan) == 0) {
    Kcent.exist <- Kcent
  }
  
  ### Perform hierarchical clustering on centroids, with Ward.D2 method
  hc <- hclust(dist(Kcent.exist[,!names(Kcent.exist) %in% c("Pedogenon")]), method="ward.D2")
  
  ### Extract labels
  hc.labels <- hc %>% as.dendrogram(.) %>% labels %>% as.numeric()
  
  ### Extract the membership from the tree
  dend <- hc %>% as.dendrogram(.)
  
  if(is.null(branchN)){
    branchN = nrow(Kcent.exist) ### We ignore the grouping of classes with the dendrogram
  }
  
  Kcent.exist$branch <- dend %>% dendextend:::cutree.dendrogram(., k = branchN)
  
  # branch.centroids <- as.data.frame(cbind(c(as.numeric(as.character(rownames(Kcent.exist)))),
  #                                         as.numeric(as.character(dendextend:::cutree.dendrogram(dend,k = branchN)))))
  
  branch.centroids <- Kcent.exist[,c("Pedogenon", "branch")]
  branch.centroids$Pedogenon <- as.numeric(branch.centroids$Pedogenon)
  branch.centroids$branch <- as.numeric(branch.centroids$branch)
  colnames(branch.centroids) <- c("Centroid", "Branch")
  
  ### sort the dataframe of branch and Pedogenon by the dendrogram labels
  # This line, using functions from dplyr or tidyverse does not work anymore
  # branch.centroids.ord <- branch.centroids %>% left_join(tibble(Centroid = hc.labels), by = "Centroid")
  
  reorder_idx <- match(hc.labels,branch.centroids$Centroid) # Saving indices for how to reorder `branch.centroids$Centroid` to match `hc.labels`
  branch.centroids.ord <- branch.centroids[reorder_idx,]
  
  if (manual == "no") {
    
    numbs.pal <- c((table(Kcent.exist$branch)))
    branch.count <- as.data.frame(cbind(c(1:branchN), numbs.pal))
    colnames(branch.count) <- c("Branch", "Count")
    #branch.count <- branch.count[order(- branch.count$Count),]
    branch.count <- branch.count %>% arrange(., -Count)
    
    ###Assign color to each 
    branch.centroids.ord$colors <- NA
    
    for(i in 1:length(numbs.pal)){
      ## Generate as many colors for each pallete as centroids in the branch
      branch.centroids.ord[branch.centroids.ord$Branch == branch.count[i,]$Branch,]$colors <-
        sequential_hcl(pal.names[[i]], n = branch.count[i,]$Count)
    }
  } else if (manual == "yes") {
    
    branch.centroids.ord$colors <- mypalette
    
  }
  
  ### Create legend
  ### Reorder the colors depending on the labels
  legend.plot <- dend %>%  set("labels_col", branch.centroids.ord$colors) %>% 
    set("branches_k_color", branch.centroids.ord$colors) %>% 
    set("branches_lwd", 4) %>% 
    set("labels_cex", 2)
  
  pdf(file = paste0("Map_legend",legend.name,".pdf"), width = 10, height = 5 )
  plot(legend.plot,
      # main = "Soil districts",
       horiz = FALSE) # change color 
  dev.off()
  
  jpeg(file = paste0("Map_legend",legend.name,".jpeg"),
       width = 900,
       height = 500 )
  plot(legend.plot,
       # main = "Soil districts",
       horiz = FALSE) # change color 
  dev.off()
  
  
  ### Now, reorder by Pedogenon class
  branch.centroids.ord <- branch.centroids.ord %>% arrange(., Centroid)
  #branch.centroids.ord <- branch.centroids.ord[order(branch.centroids.ord$Centroid),]
  
  ### Create palette for leaflet
  #pal <- branch.centroids.ord$colors
  
  # binpal <- colorBin(palette = branch.centroids.ord$colors,
  #                    bins = c(as.numeric(as.character(rownames(Kcent.exist))),
  #                             tail(as.numeric(as.character(rownames(Kcent.exist))),1)+1),
  #                    na.color = "transparent")
  
  ### Project the map into the leaflet projection
  if (need.proj == TRUE) {
    kmap <- terra::project(kmap, "epsg:3857", method="near")
    
  } else if (need.proj == FALSE) { 
    kmap <- kmap
  }
  
  binpal <- colorBin(palette = branch.centroids.ord$colors,
                     bins = c(branch.centroids.ord$Centroid,
                              tail(branch.centroids.ord$Centroid,1)+1),
                     na.color = "transparent")
  
  map.out <- leaflet() %>%
    # Base groups
    addTiles(group="OSM (default)") %>%
    ##addProviderTiles("Esri.WorldImagery", group = "World Imagery") %>% # , group = "World Imagery"
    ##addProviderTiles("OpenTopoMap", group = "Topo Map") %>%
    addRasterImage(kmap, opacity = 1, colors=binpal, project=FALSE, 
                   maxBytes = 300000000, group = "Pedogenons") %>%
    addScaleBar("bottomright",
                options = scaleBarOptions(metric=TRUE,
                                          maxWidth=200))
  # %>%
    #fitBounds(lng1=140, lat1=-38, lng2=154, lat2=-28) %>%
    #leafem::addMouseCoordinates() %>%
    # addLayersControl(
    #   baseGroups = c("OSM (default)"),
    #   overlayGroups = c("soil districts"),
    #   options = layersControlOptions(collapsed = FALSE)
    # )
    # 
  #mapshot(map.out, file = paste0(OutDir,"/Map_",legend.name,".pdf"), remove_url = FALSE)
  output <- list("hc"=hc, "branch.centroids.ord"=branch.centroids.ord,
                 "legend.plot"=legend.plot, "map.out"=map.out)
  return(output)
  
}


map.PdGn.Comb26.k9 <- viz.map.legend.pal(kmodel = pedogenons.eus.comb16.k9,
                                         manual = "yes",
                                         mypalette = mypalette,
                                         legend.name = "SoilDistrictsCmb16K9",
                                         kmap = k.raster,need.proj = TRUE)


### 10. Assign soil health indicators (soil condition indications) to each pedogenon class and plot values
### TOC, TOC:clay, EC, pH, bulk density, N, Mg, K, Ca, etc.
### Boxplots by pedogenon
### Boxplots by vegetation class (conifers vs broadleaved, native vs exotic, etc.)
### Diagram phases by pedogenon class? e.g., bulk density vs TOC or TOC/clay

PdGn.Cmb16.K9 <- terra::rast("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModeling/PdGnCmb16K9.tif")
plot(PdGn.Cmb16.K9)

BASONET.soil.sf <-st_as_sf(BASONET.soil, coords = c("UTM_X","UTM_Y"), crs = 25830)
#BASONET.soil.WGS84 <-  st_transform(BASONET.soil.sf, 4326)

### Extract
harmonisedsoil.pedogenon <- terra::extract(x=PdGn.Cmb16.K9, y=Soil.df.harmonised_sf, method="simple")
BASONET.pedogenon <- terra::extract(x=PdGn.Cmb16.K9, y=BASONET.soil.sf, method="simple")

### Bind to the soil data
Soil.df.harmonised.PdGn <- cbind(Soil.df.harmonised,harmonisedsoil.pedogenon)
BASONET.PdGn <- cbind(BASONET.soil,BASONET.pedogenon)
summary(Soil.df.harmonised.PdGn)
summary(BASONET.PdGn)
BASONET.PdGn$Depth_interval <- ifelse(BASONET.PdGn$Upper_limit==0, "0-20 cm", "20-40 cm")


# 10.1. Electric conductivity ---------------------------------------------


### Order factor levels as per dendrogram
#BASONET.PdGn$PdGn.f <- factor(as.character(BASONET.PdGn$PdGn),
 #                             levels= myorder.f)

ggplot(data=BASONET.PdGn[!is.na(BASONET.PdGn$EC) & !is.na(BASONET.PdGn$PdGn), ]) +  
  geom_boxplot(aes(x = as.factor(PdGn), y= EC, fill=as.factor(PdGn))) +
  ylab(label="Electric Conductivity") +
  xlab(label="Soil district") +
  theme_bw() +
  scale_fill_manual(values= map.PdGn.Comb26.k9$branch.centroids.ord$colors,
                    name="Soil district")+ 
  facet_wrap(~ Depth_interval) 

BASONET.EC <- BASONET.PdGn[!is.na(BASONET.PdGn$EC) & !is.na(BASONET.PdGn$PdGn), ] %>%
  group_by(., PdGn, Depth_interval ) %>%
  summarize(.,
            EC.min = round(min(EC),1),
            EC.q25 = round(quantile(EC,0.25 ),1),
            EC.mean = round(mean(EC),1),
            EC.q75 = round(quantile(EC,0.75 ),1),
            EC.max = round(max(EC),1),
            EC.sd = round(sd(EC),1),
            count=n())




# 10.1. TOC:Clay ---------------------------------------------

### Order factor levels as per dendrogram
#BASONET.PdGn$PdGn.f <- factor(as.character(BASONET.PdGn$PdGn),
#                             levels= myorder.f)

summary(BASONET.PdGn$Clay)

BASONET.PdGn$TOC_g_kg <- BASONET.PdGn$TOC *10
BASONET.PdGn$Clay_g_kg <- BASONET.PdGn$Clay *10
BASONET.PdGn$TOC_Clay <- BASONET.PdGn$TOC_g_kg / BASONET.PdGn$Clay_g_kg


ggplot(data=BASONET.PdGn[!is.na(BASONET.PdGn$TOC) & 
                           !is.na(BASONET.PdGn$Clay) &
                           !is.na(BASONET.PdGn$PdGn), ]) +  
  geom_boxplot(aes(x = as.factor(PdGn), y= TOC_Clay, fill=as.factor(PdGn))) +
  ylab(label="TOC / Clay") +
  xlab(label="Soil district") +
  theme_bw() +
  scale_fill_manual(values= map.PdGn.Comb26.k9$branch.centroids.ord$colors,
                    name="Soil district")+ 
  facet_wrap(~ Depth_interval) +
  geom_hline(yintercept =1/13,
             color = "red",
             lwd = 1,  
             linetype = "dashed")+
  ylim(0,0.6)

BASONET.PdGn$TOC_Clay_Healthy <- ifelse(BASONET.PdGn$TOC_Clay > 1/13, "healthy", "unhealthy")


BASONET.TOC_Clay <- BASONET.PdGn[!is.na(BASONET.PdGn$TOC_Clay) & !is.na(BASONET.PdGn$PdGn), ] %>%
  group_by(., PdGn, Depth_interval ) %>%
  summarize(.,
            TOC_Clay.min = round(min(TOC_Clay, na.rm=TRUE),2),
            TOC_Clay.q25 = round(quantile(TOC_Clay,0.25, na.rm=TRUE ),2),
            TOC_Clay.mean = round(mean(TOC_Clay,na.rm=TRUE),2),
            TOC_Clay.q75 = round(quantile(TOC_Clay,0.75, na.rm=TRUE),2),
            TOC_Clay.max = round(max(TOC_Clay, na.rm=TRUE),2),
            TOC_Clay.sd = round(sd(TOC_Clay, na.rm=TRUE),2),
            count=n())


BASONET.TOC_Clay.Health <- BASONET.PdGn[!is.na(BASONET.PdGn$TOC_Clay) & !is.na(BASONET.PdGn$PdGn), ] %>%
  group_by(., PdGn,Depth_interval, TOC_Clay_Healthy ) %>%
  summarize(.,
            TOC_Clay.mean = round(mean(TOC_Clay,na.rm=TRUE),2),
            count=n())

#BASONET.TOC_Clay.Health$count <- ifelse(is.na(BASONET.TOC_Clay.Health$count),0, BASONET.TOC_Clay.Health$count)

BASONET.TOC_Clay.Health.wide  <- tidyr::pivot_wider(BASONET.TOC_Clay.Health, 
                                                    id_cols=c(PdGn ,Depth_interval ),
                                                    names_from = "TOC_Clay_Healthy",
                                                    values_from = "count")
BASONET.TOC_Clay.Health.wide$healthy <- 
  ifelse(is.na(BASONET.TOC_Clay.Health.wide$healthy), 0, BASONET.TOC_Clay.Health.wide$healthy)


BASONET.TOC_Clay.Health.wide$count_tot <- BASONET.TOC_Clay.Health.wide$healthy + BASONET.TOC_Clay.Health.wide$unhealthy


BASONET.TOC_Clay.Health.wide$unhealthy_perc <- round(BASONET.TOC_Clay.Health.wide$unhealthy / BASONET.TOC_Clay.Health.wide$count_tot * 100, digits=0)
BASONET.TOC_Clay.Health.wide$unhealthy_perc
hist(BASONET.TOC_Clay.Health.wide$unhealthy_perc, breaks=10)
BASONET.TOC_Clay.Health.wide <- arrange(BASONET.TOC_Clay.Health.wide, unhealthy_perc)


# 10.3 Subsoil compaction -------------------------------------------------

library(soiltexture)
TT.plot( class.sys = "USDA.TT" )


### No missing data allowed
BASONET.PdGn.Texture <- BASONET.PdGn[!is.na(BASONET.PdGn$Clay) & !is.na(BASONET.PdGn$PdGn),]

BASONET.PdGn.Texture$CLAY <- BASONET.PdGn.Texture$Clay
BASONET.PdGn.Texture$SILT <- BASONET.PdGn.Texture$Silt
BASONET.PdGn.Texture$SAND <- BASONET.PdGn.Texture$Sand
BASONET.PdGn.Text <- TT.normalise.sum( as.data.frame(BASONET.PdGn.Texture[,c("CLAY", "SILT", "SAND")]),
                                       css.names = c("CLAY", "SILT", "SAND"))

BASONET.PdGn.Texture$texture.class <- TT.points.in.classes(BASONET.PdGn.Text, class.sys = "USDA.TT", PiC.type = "t") 

ggplot(data=BASONET.PdGn.Texture[!is.na(BASONET.PdGn.Texture$Bulk_Density) & 
                                   BASONET.PdGn.Texture$Depth_interval == "20-40 cm" &
                                   !is.na(BASONET.PdGn.Texture$PdGn), ]) +  
  geom_boxplot(aes(x = as.factor(PdGn), y= Bulk_Density, fill=as.factor(PdGn))) +
  ylab(label="Bulk density") +
  xlab(label="Soil district") +
  theme_bw() +
  scale_fill_manual(values= map.PdGn.Comb26.k9$branch.centroids.ord$colors,
                    name="Soil district")+ 
  facet_wrap(~ Depth_interval) 

### Criteria soil health for subsoil bulk density
# 1 Cl clay
# 2 SiCl silty clay
# 3 SaCl sandy clay
# 4 ClLo clay loam
# 5 SiClLo silty clay loam
# 6 SaClLo sandy clay loam
# 7 Lo loam
# 8 SiLo silty loam
# 9 SaLo sandy loam
# 10 Si silt
# 11 LoSa loamy sand
# 12 Sa sand

criteriaSubsoilcompaction <- function(texture, bulk_density, clay) {
  
  cond1 <-  texture %in% c("Sa","LoSa", "SaLo", "Lo", "Lo, SiLo",  "SaClLo, SaLo") & 
    bulk_density < 1.8
  
  cond2 <- texture %in% c("SaClLo", "ClLo", "Si", "SiLo") &
    bulk_density < 1.75
  
  cond3 <- texture %in% c("SiClLo", "SiCl, SiClLo") &
    bulk_density < 1.65
  
  cond4 <- texture %in% c("SiCl","SaCl") &
    bulk_density < 1.58
  
  cond5 <- texture == "ClLo" & clay <= 45 & clay >=35
  
  cond6 <- texture == "Cl" & bulk_density < 1.47
  
  condition <- ifelse(cond1, "healthy",
                    ifelse(cond2, "healthy",
                           ifelse(cond3, "healthy",
                                  ifelse(cond4, "healthy",
                                         ifelse(cond5, "healthy",
                                                ifelse(cond6, "healthy", "unhealthy"))))))
  return(condition)
  
}

BASONET.SubsoilCompaction.PdGn <- BASONET.PdGn.Texture[!is.na(BASONET.PdGn.Texture$Bulk_Density) & 
                                                         BASONET.PdGn.Texture$Depth_interval == "20-40 cm" &
                                                         !is.na(BASONET.PdGn.Texture$PdGn), ]


BASONET.SubsoilCompaction.PdGn$Condition <- criteriaSubsoilcompaction(texture = BASONET.SubsoilCompaction.PdGn$texture.class, 
                                                                      bulk_density = BASONET.SubsoilCompaction.PdGn$Bulk_Density,
                                                                      clay = BASONET.SubsoilCompaction.PdGn$Clay)

length(BASONET.SubsoilCompaction.PdGn$Condition)
perc.unhealthy <- length(BASONET.SubsoilCompaction.PdGn$Condition[BASONET.SubsoilCompaction.PdGn$Condition=="unhealthy"])/length(BASONET.SubsoilCompaction.PdGn$Condition) *100

ggplot(data=BASONET.SubsoilCompaction.PdGn) +  
  geom_boxplot(aes(x = as.factor(Condition), y= Bulk_Density, fill=as.factor(PdGn))) +
  ylab(label="Bulk density") +
  xlab(label="Soil district") +
  theme_bw() +
  scale_fill_manual(values= map.PdGn.Comb26.k9$branch.centroids.ord$colors,
                    name="Soil district") # + 
  #facet_wrap(~ as.factor(PdGn))
BASONET.SubsoilCompaction.PdGn.summary <- BASONET.SubsoilCompaction.PdGn %>%
  group_by(., PdGn, Condition ) %>%
  summarize(.,
            BulkDensity.mean = round(mean(Bulk_Density,na.rm=TRUE),2),
            count=n())


# 11. Other indicators ----------------------------------------------------

### Soil Phosphorus
BASONET.PdGn$P_mg_g <- BASONET.PdGn$P_ppm / 1000

ggplot(data=BASONET.PdGn[!is.na(BASONET.PdGn$P_mg_g) & 
                           !is.na(BASONET.PdGn$PdGn), ]) +  
  geom_boxplot(aes(x = as.factor(PdGn), y= P_mg_g, fill=as.factor(PdGn))) +
  ylab(label="Phosphorus (mg/g)") +
  xlab(label="Soil district") +
  theme_bw() +
  scale_fill_manual(values= map.PdGn.Comb26.k9$branch.centroids.ord$colors,
                    name="Soil district")+ 
  facet_wrap(~ Depth_interval) 


### Soil Nitrogen
BASONET.PdGn$N_mg_g <- BASONET.PdGn$TN_ppm / 1000

ggplot(data=BASONET.PdGn[!is.na(BASONET.PdGn$N_mg_g) & 
                           !is.na(BASONET.PdGn$PdGn), ]) +  
  geom_boxplot(aes(x = as.factor(PdGn), y= N_mg_g, fill=as.factor(PdGn))) +
  ylab(label="Nitrogen (mg/g)") +
  xlab(label="Soil district") +
  theme_bw() +
  scale_fill_manual(values= map.PdGn.Comb26.k9$branch.centroids.ord$colors,
                    name="Soil district")+ 
  facet_wrap(~ Depth_interval) +
  ylim(0,10)

### pH - without criteria
ggplot(data=BASONET.PdGn[!is.na(BASONET.PdGn$pH) & 
                           !is.na(BASONET.PdGn$PdGn), ]) +  
  geom_boxplot(aes(x = as.factor(PdGn), y= pH, fill=as.factor(PdGn))) +
  ylab(label="pH") +
  xlab(label="Soil district") +
  theme_bw() +
  scale_fill_manual(values= map.PdGn.Comb26.k9$branch.centroids.ord$colors,
                    name="Soil district")+ 
  facet_wrap(~ Depth_interval) 

### Topsoil compaction - without criteria
ggplot(data=BASONET.PdGn[!is.na(BASONET.PdGn$Bulk_Density) & 
                           BASONET.PdGn$Depth_interval == "0-20 cm" &
                                   !is.na(BASONET.PdGn$PdGn), ]) +  
  geom_boxplot(aes(x = as.factor(PdGn), y= Bulk_Density, fill=as.factor(PdGn))) +
  ylab(label="Bulk density") +
  xlab(label="Soil district") +
  theme_bw() +
  scale_fill_manual(values= map.PdGn.Comb26.k9$branch.centroids.ord$colors,
                    name="Soil district")+ 
  facet_wrap(~ Depth_interval) 



# 12. PERMANOVA -----------------------------------------------------------

BASONET.PdGn.subset <- BASONET.PdGn[complete.cases(BASONET.PdGn[,c("PdGn","Depth_interval",
                                                                    "Silt","Clay","pH","TOC","TN_ppm","K_ppm",
                                                                    "Ca_ppm", "Mg_ppm","CECef")]) &
                                      BASONET.PdGn$Date == 2021, ]

Response <- BASONET.PdGn.subset[,c("Silt","Clay","pH","TOC",
                                   "TN_ppm","K_ppm","Ca_ppm", "Mg_ppm","CECef")]

### PERMANOVA
library(vegan)
set.seed(1984)
permanova.BASONET <- adonis2(Response ~ PdGn+Depth_interval,
                             data=BASONET.PdGn.subset,
                             method = "mahalanobis", by="terms",
                             permutations=9999)
permanova.BASONET
densityplot(permustats(permanova.BASONET))

BASONET.dist <- vegdist(x = Response, method = "mahalanobis")

set.seed(1984)
BASONET_disper <- betadisper(d =BASONET.dist, group=BASONET.PdGn.subset$PdGn)
plot(BASONET_disper)

# PCoA.BASONET <- cmdscale(BASONET.dist, k=3)
# plot(PCoA.BASONET,pch=19,
#      col= map.PdGn.Comb26.k9$branch.centroids.ord$colors[BASONET.PdGn.subset$PdGn])

save.image("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/6.SoilCondition/6.SoilCondition22032024.RData")
