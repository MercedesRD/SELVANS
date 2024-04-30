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

### Load the R session from the previous script
load("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/6.SoilCondition/6.SoilCondition22032024.RData")

### Load the k-means model
load("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModeling/scorpan_combi_16/kmeans_scorpanID16.k9.RData")

# 1. SCORPAN variables - GIS layers ------------------------------------------

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

### This is combination 16
var.columns <- c(t(relief.combi)[,16])
var.columns <- var.columns[!is.na(var.columns)] ### eliminate NA if any


# 2. Plot SCORPAN vs Soil Profile distances - individual soil profiles --------

### this plot is independent of cluster assignment, just depends on scorpan selection

### I load the result from the selected combination (for starters)
BASONET_assignmentCmb16K9 <- read.csv("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModeling/scorpan_combi_16/BASONET_assign_16.csv")

### Predict pedogenon with the model fit with more observations (the pedogenon number class may have changed)
#BASONET_assignmentCmb16K9$PdGn <- predict(, BASONET_assignmentCmb16K9)
kmeans_clorpt.jk

### We have 804 observations for 428 unique soil profile locations
length(unique(BASONET_assignmentCmb16K9$newID))

### Arrange by soil profile ID and depth.
### Soil properties of interest:
### Silt, Clay, pH, CECef, EC, TOC.
BASONET_subset <- BASONET_assignmentCmb16K9[,c("Dataset","newID","Layer_depth","K.9",
                                               fixed.columns,var.columns,
                                               "Silt","Clay","CECef",
                                               "TOC","pH","EC")]

### Rearrange by uniqueID and depth layer
BASONET_subset <- BASONET_subset %>%
  arrange(., newID, Layer_depth) %>% as.data.frame()

### Subset only complete observations
BASONET_subset <- BASONET_subset[complete.cases(BASONET_subset),]
### 613 observations at 347 soil profile locations
length(unique(BASONET_subset$newID))

### Calculate SCORPAN distance between soil profile locations

### scorpan dataframe for the combination j
scorpan.j <- dplyr::distinct(BASONET_subset[,c("newID",fixed.columns,var.columns)])
length(unique(scorpan.j$newID))
dim(scorpan.j)
basonet_combi16_SoilProfileIds <- unique(scorpan.j$newID)
### Rename rownames by "newID" and use later as id
rownames(scorpan.j) <- scorpan.j$newID 

### Rescale the SCORPAN variables
all.equal(colnames(scorpan.j),dimnames(C)[[1]])

### Reorder columns
scorpan.j <- scorpan.j[,dimnames(C)[[1]]]
all.equal(colnames(scorpan.j),dimnames(C)[[1]])

### Decorrelate the SCORPAN variables at the locations of the soil observations
B.scorpan.j.rs <- as.matrix(scorpan.j) %*% solve(C) 

### Calculate distance (Mahalanobis distance between scorpan variables)
dist.B.scorpan <- dist(B.scorpan.j.rs, diag=FALSE, upper =FALSE)
### As matrix
dist.B.scorpan.M <- as.matrix(dist.B.scorpan)
### As vector
dist.B.scorpan.v <- as.vector(dist.B.scorpan)
hist(dist.B.scorpan.v)
length(dist.B.scorpan.v)

#### Distance between soil properties
BASONET_subset.soil <- BASONET_subset[,c("newID","Layer_depth","Silt","Clay","CECef","TOC","pH","EC")]

### subset by depth layers
basonet.soildf.hors <- split(BASONET_subset.soil, BASONET_subset$Layer_depth)

### Rownames for identifying the soil profiles
rownames(basonet.soildf.hors$`000_020_cm`) <- basonet.soildf.hors$`000_020_cm`$newID
rownames(basonet.soildf.hors$`020_040_cm`) <- basonet.soildf.hors$`020_040_cm`$newID

### Within each element of the list, this is, horizon, calculate distance between soil profiles
horizon.dists <- lapply(basonet.soildf.hors, FUN = function(x) {as.matrix(dist(x[,c("Silt","Clay","CECef","TOC","pH","EC")]))} )

### Do dimnames inform of X and Y soil profiles?
setdiff(dimnames(horizon.dists$`000_020_cm`)[[1]], dimnames(horizon.dists$`020_040_cm`)[[1]])
setdiff(dimnames(horizon.dists$`020_040_cm`)[[1]],dimnames(horizon.dists$`000_020_cm`)[[1]])
setdiff(basonet_combi16_SoilProfileIds, dimnames(horizon.dists$`000_020_cm`)[[1]])

### From Chat GPT (even if I complain, it helped me finding the solution to my problem because my brain does not work today
# Create a function to reshape matrices
reshape_and_fill <- function(mat, dimnames) {
  # Create a new matrix with desired dimnames filled with NA
  new_mat <- matrix(NA, nrow = length(dimnames[[1]]), ncol = length(dimnames[[2]]), dimnames = dimnames)
  # Assign values from original matrix to corresponding positions
  new_mat[rownames(mat), colnames(mat)] <- mat
  return(new_mat)
}

# Create a new matrix with desired dimnames filled with NA
new_mat_000_020_cm <- reshape_and_fill(horizon.dists$`000_020_cm`, 
                                       dimnames =list(as.character(basonet_combi16_SoilProfileIds),
                                                      as.character(basonet_combi16_SoilProfileIds)))
  
new_mat_020_040_cm <- reshape_and_fill(horizon.dists$`020_040_cm`, 
                                       dimnames =list(as.character(basonet_combi16_SoilProfileIds),
                                                      as.character(basonet_combi16_SoilProfileIds)))

### Create an array 
basonetSoilProfiless.dist.array <- array(c(new_mat_000_020_cm,new_mat_020_040_cm),
                                         dim = c(dim(new_mat_000_020_cm), 2))

### Average soil profile distances, from the array
basonet.ave.prfl.dist <- apply(basonetSoilProfiless.dist.array, 
                                 MARGIN = c(1,2), 
                                 function(x) mean(x, na.rm=TRUE))

### Transform into dist object
basonet.ave.prfl.dist.d <- as.dist(basonet.ave.prfl.dist, upper = FALSE, diag = FALSE)
### As vector
basonet.ave.prfl.dist.v <- as.vector(basonet.ave.prfl.dist.d)
hist(basonet.ave.prfl.dist.v)
length(basonet.ave.prfl.dist.v)

### With ggplot
relationship.scorpan.SoilProfile.df <- data.frame(scorpan.dist = dist.B.scorpan.v,
                                                 soil.profile.dist = basonet.ave.prfl.dist.v)

ggplot(relationship.scorpan.SoilProfile.df,
       aes(x=scorpan.dist, y=soil.profile.dist)) +
  geom_point(alpha=0.1, shape=1, )+
  geom_smooth(fill="#95c8ff", color="#2b90ff", lwd=2)+
  theme_bw()+
  labs(title = "BASONET dataset")+
  labs(x="Mahalanobis distance in the SCORPAN space")+
  labs(y="Euclidean soil profile distance")


# 3. Distance between pedogenon centroids ---------------------------------

centroids.k9 <- kmeans_clorpt.jk$centroids

### Calculate distance (Mahalanobis distance between scorpan variables)
dist.B.scorpan <- dist(centroids.k9, diag=FALSE, upper =FALSE)
### As matrix
dist.B.scorpan.M <- as.matrix(dist.B.scorpan)
### As vector
dist.B.scorpan.v <- as.vector(dist.B.scorpan)
hist(dist.B.scorpan.v, breaks=10)
length(dist.B.scorpan.v)

### Change name of pedogenon column
colnames(BASONET_subset)[colnames(BASONET_subset) =="K.9" ] <- "PdGn"

### Calculate inter-centroid distance
BASONET.DF.inter <- Dinter.function(df.soil = BASONET_subset,
                                    uniqueID = "newID",
                                    depth.var = "Layer_depth",
                                    target.vars = c("Silt","Clay","CECef","TOC","pH","EC"))


### As vector
BASONET.dist.inter.v <- as.vector(BASONET.DF.inter)
hist(BASONET.dist.inter.v, breaks=10)
length(BASONET.dist.inter.v)

### With ggplot
relationship.centroids.distances.k9 <- data.frame(scorpan.dist = dist.B.scorpan.v,
                                                  soil.profile.dist = BASONET.dist.inter.v)

ggplot(relationship.centroids.distances.k9,
       aes(x=scorpan.dist, y=soil.profile.dist)) +
  geom_point(alpha=1, shape=19, size=2)+
  geom_smooth(fill="#95c8ff", color="#2b90ff", lwd=1)+
  theme_bw()+
  labs(title = "Pedogenon centroids")+
  labs(x="Mahalanobis distance in the SCORPAN space")+
  labs(y="BASONET soil profiles inter-centroid distances")



# ### Repeat with a different combinations --------------------------------

### This is combination 38
var.columns <- c(t(relief.combi)[,38])
var.columns <- var.columns[!is.na(var.columns)] ### eliminate NA if any


### I load the result from the selected combination (for starters)
BASONET_assignmentCmb38 <- read.csv("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModeling/scorpan_combi_38/BASONET_assign_38.csv")

cl_indices_comb38 <- read_csv("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModeling/scorpan_combi_38/cl_indices_comb38.csv")
View(cl_indices_comb38)

### Transpose dataframe, columns to rows etc

### what combination has max number of minimum observations for a greater number of clusters?
index.table.long <- tidyr::pivot_longer(cl_indices_comb38, 
                                        cols =  starts_with("K."), 
                                        names_to = "clusters",
                                        values_to = "index")

index.table.wide <- tidyr::pivot_wider(index.table.long, 
                                       id_cols=c(clusters ),
                                       names_from = "Index",
                                       values_from = "index")

index.table.wide$K <- as.numeric(gsub(x=index.table.wide$clusters, pattern="K.",replacement = ""))

index.table.wide$clusters <- factor(index.table.wide$clusters,
                                    levels=paste0("K.",2:140))

ggplot(index.table.wide,
       aes(x=K, y = B_Din_Dex)) +
  geom_point() +
  scale_fill_viridis() +
  ylim(0,1.5)

ggplot(index.table.wide,
       aes(x=K, y = `Calinski-Harbasz`)) +
  geom_point() +
  scale_fill_viridis() 

ggplot(index.table.wide,
       aes(x=K, y = `BetweenSS to TotalSS`)) +
  geom_point() +
  scale_fill_viridis() 

ggplot(index.table.wide,
       aes(x=K, y = B_Median_sites)) +
  geom_point() +
  scale_fill_viridis()

colnames(ggplot(index.table.wide,
       aes(x=K, y = `BetweenSS to TotalSS`)) +
  geom_point() +
  scale_fill_viridis())

index.table.wide[index.table.wide$`Calinski-Harbasz`==(max(index.table.wide$`Calinski-Harbasz`)),]
index.table.wide[index.table.wide$`Calinski-Harbasz`==(max(index.table.wide$`Calinski-Harbasz`)),]
sort(index.table.wide[index.table.wide$B_Din_Dex <= 0.6,]$K)

### We have 804 observations for 428 unique soil profile locations
length(unique(BASONET_assignmentCmb38$newID))

### Arrange by soil profile ID and depth.
### Soil properties of interest:
### Silt, Clay, pH, CECef, EC, TOC.
BASONET_subset <- BASONET_assignmentCmb38[,c("Dataset","newID","Layer_depth","K.30",
                                               fixed.columns,var.columns,
                                               "Silt","Clay","CECef",
                                               "TOC","pH","EC")]

### Rearrange by uniqueID and depth layer
BASONET_subset <- BASONET_subset %>%
  arrange(., newID, Layer_depth) %>% as.data.frame()

### Subset only complete observations
BASONET_subset <- BASONET_subset[complete.cases(BASONET_subset),]
### 613 observations at 347 soil profile locations
length(unique(BASONET_subset$newID))

### Calculate SCORPAN distance between soil profile locations

### scorpan dataframe for the combination j
scorpan.j <- dplyr::distinct(BASONET_subset[,c("newID",fixed.columns,var.columns)])
length(unique(scorpan.j$newID))
dim(scorpan.j)
basonet_combi38_SoilProfileIds <- unique(scorpan.j$newID)
### Rename rownames by "newID" and use later as id
rownames(scorpan.j) <- scorpan.j$newID 

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
CLORPT.df <- scorpanSample[,c(fixed.columns,var.columns)]

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

### Rescale the SCORPAN variables
all.equal(colnames(scorpan.j),dimnames(C)[[1]])

### Reorder columns
scorpan.j <- scorpan.j[,dimnames(C)[[1]]]
all.equal(colnames(scorpan.j),dimnames(C)[[1]])

### Decorrelate the SCORPAN variables at the locations of the soil observations
B.scorpan.j.rs <- as.matrix(scorpan.j) %*% solve(C) 

### Calculate distance (Mahalanobis distance between scorpan variables)
dist.B.scorpan <- dist(B.scorpan.j.rs, diag=FALSE, upper =FALSE)
### As matrix
dist.B.scorpan.M <- as.matrix(dist.B.scorpan)
### As vector
dist.B.scorpan.v <- as.vector(dist.B.scorpan)
hist(dist.B.scorpan.v)
length(dist.B.scorpan.v)

#### Distance between soil properties
BASONET_subset.soil <- BASONET_subset[,c("newID","Layer_depth","Silt","Clay","CECef","TOC","pH","EC")]

### subset by depth layers
basonet.soildf.hors <- split(BASONET_subset.soil, BASONET_subset$Layer_depth)

### Rownames for identifying the soil profiles
rownames(basonet.soildf.hors$`000_020_cm`) <- basonet.soildf.hors$`000_020_cm`$newID
rownames(basonet.soildf.hors$`020_040_cm`) <- basonet.soildf.hors$`020_040_cm`$newID

### Within each element of the list, this is, horizon, calculate distance between soil profiles
horizon.dists <- lapply(basonet.soildf.hors, FUN = function(x) {as.matrix(dist(x[,c("Silt","Clay","CECef","TOC","pH","EC")]))} )

### Do dimnames inform of X and Y soil profiles?
setdiff(dimnames(horizon.dists$`000_020_cm`)[[1]], dimnames(horizon.dists$`020_040_cm`)[[1]])
setdiff(dimnames(horizon.dists$`020_040_cm`)[[1]],dimnames(horizon.dists$`000_020_cm`)[[1]])
setdiff(basonet_combi16_SoilProfileIds, dimnames(horizon.dists$`000_020_cm`)[[1]])

### From Chat GPT (even if I complain, it helped me finding the solution to my problem because my brain does not work today
# Create a function to reshape matrices
reshape_and_fill <- function(mat, dimnames) {
  # Create a new matrix with desired dimnames filled with NA
  new_mat <- matrix(NA, nrow = length(dimnames[[1]]), ncol = length(dimnames[[2]]), dimnames = dimnames)
  # Assign values from original matrix to corresponding positions
  new_mat[rownames(mat), colnames(mat)] <- mat
  return(new_mat)
}

# Create a new matrix with desired dimnames filled with NA
new_mat_000_020_cm <- reshape_and_fill(horizon.dists$`000_020_cm`, 
                                       dimnames =list(as.character(basonet_combi16_SoilProfileIds),
                                                      as.character(basonet_combi16_SoilProfileIds)))

new_mat_020_040_cm <- reshape_and_fill(horizon.dists$`020_040_cm`, 
                                       dimnames =list(as.character(basonet_combi16_SoilProfileIds),
                                                      as.character(basonet_combi16_SoilProfileIds)))

### Create an array 
basonetSoilProfiless.dist.array <- array(c(new_mat_000_020_cm,new_mat_020_040_cm),
                                         dim = c(dim(new_mat_000_020_cm), 2))

### Average soil profile distances, from the array
basonet.ave.prfl.dist <- apply(basonetSoilProfiless.dist.array, 
                               MARGIN = c(1,2), 
                               function(x) mean(x, na.rm=TRUE))

### Transform into dist object
basonet.ave.prfl.dist.d <- as.dist(basonet.ave.prfl.dist, upper = FALSE, diag = FALSE)
### As vector
basonet.ave.prfl.dist.v <- as.vector(basonet.ave.prfl.dist.d)
hist(basonet.ave.prfl.dist.v)
length(basonet.ave.prfl.dist.v)

### With ggplot
relationship.scorpan.SoilProfile.df <- data.frame(scorpan.dist = dist.B.scorpan.v,
                                                  soil.profile.dist = basonet.ave.prfl.dist.v)

ggplot(relationship.scorpan.SoilProfile.df,
       aes(x=scorpan.dist, y=soil.profile.dist)) +
  geom_point(alpha=0.1, shape=1, )+
  geom_smooth(fill="#95c8ff", color="#2b90ff", lwd=2)+
  theme_bw()+
  labs(title = "BASONET dataset")+
  labs(x="Mahalanobis distance in the SCORPAN space")+
  labs(y="Euclidean soil profile distance")

### It does not change much
### I try with another k-means model, k=30
load("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/5.PedogenonModeling/scorpan_combi_38/kmeans_scorpanID38.k30.RData")
centroids.k30 <- kmeans_clorpt.jk$centroids

### Calculate distance (Mahalanobis distance between scorpan variables)
dist.B.scorpan <- dist(centroids.k30, diag=FALSE, upper =FALSE)
### As matrix
dist.B.scorpan.M <- as.matrix(dist.B.scorpan)
### As vector
dist.B.scorpan.v <- as.vector(dist.B.scorpan.M)
hist(dist.B.scorpan.v, breaks=10)
length(dist.B.scorpan.v)

### Change name of pedogenon column
colnames(BASONET_subset)[colnames(BASONET_subset) =="K.30" ] <- "PdGn"

### Calculate inter-centroid distance
BASONET.DF.inter <- Dinter.function(df.soil = BASONET_subset,
                                    uniqueID = "newID",
                                    depth.var = "Layer_depth",
                                    target.vars = c("Silt","Clay","CECef","TOC","pH","EC"))
### Transform to matrix
BASONET.DF.inter.M <- as.matrix(BASONET.DF.inter)

### Rename dimanmes to match the pedogenons present (with observations)
dimnames(BASONET.DF.inter.M)[[1]] <- as.character(sort(unique(BASONET_subset$PdGn)))
dimnames(BASONET.DF.inter.M)[[2]] <- as.character(sort(unique(BASONET_subset$PdGn)))

# Create a new matrix with desired dimnames filled with NA
new_mat_BASONET.DF.inter<- reshape_and_fill(BASONET.DF.inter.M, 
                                       dimnames =list(as.character(1:30),
                                                      as.character(1:30)))

### As vector
BASONET.dist.inter.v <- as.vector(new_mat_BASONET.DF.inter)
hist(BASONET.dist.inter.v, breaks=10)
length(BASONET.dist.inter.v)

### With ggplot
relationship.centroids.distances.k30 <- data.frame(scorpan.dist = dist.B.scorpan.v,
                                                  soil.profile.dist = BASONET.dist.inter.v)

### Eliminate 0 distances
relationship.centroids.distances.k30 <- relationship.centroids.distances.k30[relationship.centroids.distances.k30$scorpan.dist !=0,]
### and repeated distances (from upper diagonal)
relationship.centroids.distances.k30 <- distinct(relationship.centroids.distances.k30)

ggplot(relationship.centroids.distances.k30,
       aes(x=scorpan.dist, y=soil.profile.dist)) +
  geom_point(alpha=1, shape=19, size=2)+
  geom_smooth(fill="#95c8ff", color="#2b90ff", lwd=1)+
  theme_bw()+
  labs(title = "Pedogenon centroids K= 30")+
  labs(x="Mahalanobis distance in the SCORPAN space")+
  labs(y="BASONET soil profiles inter-centroid distances")+
  ylim(0,50)

### End of the script