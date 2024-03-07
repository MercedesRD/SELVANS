#### Pedogenon maps for France

### Date 05/12/2023
### Author: Mercedes Roman Dobarco
### Objective: Loas covariates for pedogenon mapping

# analysis packages
library(dplyr)
library(tidyverse)
library(sf)
library(terra)
library(Hmisc)
library(FactoMineR)
library(factoextra)

# visualization packages
library(ggplot2)
library(gganimate)
library(patchwork)
library(viridis)
library(scales)
library(corrplot)
#library(GGally)


# Study area ------------------------------------------------------

### France buffer
france_buffer_WGS84 <- st_read("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Administrative/france_bufferWGS84.shp")
### Transform to spatvector
france_v <- vect(france_buffer_WGS84)
myext <- ext(st_bbox(france_v))

# 1. Method for variable selection: Climate and relief variables selected from correlation plots. ------------

### SCORPAN variables

### Soil
### At the moment, no soil information will be included

### Climate
### Here we select them based on correlation plots (see power point "GenosoilPhenosoilFrance")
### Subset of climate variables were: bio1, bio3, bio4, bio8, bio9,bio2 and bio15
climDir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Climate/"
setwd(paste0(climDir,"TimeID_-9/"))
bioclim <- list.files(pattern="Fr.tif")
bioclimrast <- terra::rast(paste0(climDir,"TimeID_-9/",bioclim))
### subset of variables selected (correlation < 0.7)
BioclimSub <- terra::subset(bioclimrast, c(1,3,4,8,9,12,15))
names(BioclimSub) <- c("Bio1", "Bio3", "Bio4", "Bio8","Bio9", "Bio12", "Bio15")

### Organisms
### Vegetation - Forest cover 1000 BCE (~3000 BP)
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/vegetation/Zanon/")
forest_cover_3000BP <- terra::rast("forest_cover_3000.tif")
plot(forest_cover_3000BP)

### HYDE 3.2
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/HYDE3.2/h1000BCE")
conv_rangeland1000BCE <- terra::rast("conv_rangeland1000BC_Fr.tif")
cropland1000BCE <- terra::rast("cropland1000BC_Fr.tif")
grazing1000BCE <- terra::rast("grazing1000BC_Fr.tif")

### Create binary variable Bare soil from Biome4 data
# setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/vegetation/biome05deg")
# biome4_2900BP <- terra::rast("biome4_2900BP.tif")
# biome4_3000BP <- terra::rast("biome4_3000BP.tif")
# unique(biome4_3000BP)
# unique(biome4_2900BP)
# 
# ### Crop to France
# biome4_2900BP_Fr <- terra::crop(biome4_2900BP,france_v)
# biome4_3000BP_Fr <- terra::crop(biome4_3000BP,france_v)
# plot(biome4_3000BP_Fr) ### --> 1050 BCE
# plot(biome4_2900BP_Fr)
### Classes as 5, 6, 8, 10, 15, 23
# 5  Temperate conifer forest
# 6  Warm mixed forest
# 8  Cool conifer forest
# 10 Evegreen taiga/montane forest
# 15 Temperate sclerophyll woodland
# 23 Shrub tundra

# list of the 28 biomes assigned, including land ice
# !     1  Tropical evergreen forest
# !     2  Tropical semi-deciduous forest
# !     3  Tropical deciduous forest/woodland
# !     4  Temperate deciduous forest
# !     5  Temperate conifer forest
# !     6  Warm mixed forest
# !     7  Cool mixed forest
# !     8  Cool conifer forest
# !     9  Cold mixed forest
# !     10 Evegreen taiga/montane forest
# !     11 Deciduous taiga/montane forest
# !     12 Tropical savanna
# !     13 Tropical xerophytic shrubland
# !     14 Temperate xerophytic shrubland
# !     15 Temperate sclerophyll woodland
# !     16 Temperate broadleaved savanna
# !     17 Open conifer woodland
# !     18 Boreal parkland
# !     19 Tropical grassland
# !     20 Temperate grassland
# !     21 Desert
# !     22 Steppe tundra
# !     23 Shrub tundra
# !     24 Dwarf shrub tundra
# !     25 Prostrate shrub tundra
# !     26 Cushion forb lichen moss tundra
# !     27 Barren
# !     28 Land ice

### Create binary variable alpine grasslands from Biome4 data
#shrub_tundra <- app(biome4_3000BP_Fr, )

### Relief

### Covariates selected from correlation plot, 
### from those that I prefer (excluding curvature, basically)

# Selected variables:
# CTI, easterness, HLI, mrrtf, mrvbf, north_slope, northerness, scale_pos, slope,  SRR, srtm 
ReliefDir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Relief/"
setwd(ReliefDir)
list.files(pattern=".tif")
relief.covar <- c("cti.tif","easterness.tif","hli.tif","mrrtf.tif",
                  "mrvbf.tif","north_slope.tif","northerness.tif",
                  "scale_pos.tif","slope.tif","srr.tif","srtm.tif")
reliefrast <- terra::rast(paste0(ReliefDir,relief.covar))
names(reliefrast) <- gsub(relief.covar,replacement = "", pattern=".tif")

### Project all to EPSG:4326 because that is the CRS I am working at
reliefrast <- project(reliefrast, "EPSG:4326", method="bilinear")

# 2. Scale the numerical variables ----------------------------------------

library(foreach)
library(parallel)
library(doParallel)

ext(BioclimSub)
ext(forest_cover_3000BP)
ext(reliefrast)

### It all should have the same extent, but it fails. I resample to the extent of reliefrast
template.r <- reliefrast[[1]]
names(template.r) <- template.r
values(template.r) <- NA

BioclimSub.e <- terra::resample(BioclimSub, y = template.r, method="bilinear")
forest_cover_3000BP.e <- terra::resample(forest_cover_3000BP, y = template.r, method="bilinear")
conv_rangeland1000BCE.e <- terra::resample(conv_rangeland1000BCE, y = template.r, method="bilinear")
cropland1000BCE.e <- terra::resample(cropland1000BCE, y = template.r, method="bilinear")
grazing1000BCE.e <- terra::resample(grazing1000BCE, y = template.r, method="bilinear")

###change names
names(forest_cover_3000BP.e) <- "forest_cover"
names(conv_rangeland1000BCE.e) <- "conv_rangelands"
names(cropland1000BCE.e) <- "cropland"
names(grazing1000BCE.e) <- "grazing"

### Create mask with dem
mymask <- reliefrast[[11]]
f <- function(x) ifelse(!is.na(x), 1, NA)
mymask <- app(mymask, f)
plot(mymask)

covariates.selection.cont <- c(BioclimSub.e, ### Climate
                               forest_cover_3000BP.e,conv_rangeland1000BCE.e,cropland1000BCE.e,grazing1000BCE.e,    # organisms
                               reliefrast)   ### relief
plot(covariates.selection.cont)

### Create new directory to store raster files, scaled
#dir.create("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/scaled/")
#setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/scaled/")
raster.names <- names(covariates.selection.cont)

#detectCores()
#cl <- makeCluster(8)   ### Create cluster
#registerDoParallel(cl)
#getDoParWorkers()
# covariates.cont.out <- foreach(i=1:length(covariates.selection.cont), .packages=c("terra", "sf"),
#                                            .export = c("covariates.selection.cont","raster.names", "mymask")) %dopar% {

for(i in 1:length(raster.names)){
  print(covariates.selection.cont[[i]])
  m <- terra::mask(x=covariates.selection.cont[[i]], mask=mymask) # Mask pixels outside France
  s <- terra::scale(x=m,center=TRUE, scale=TRUE) # Scale because it is a continuous variable
  setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/scaled/") # change wd
  writeRaster(s, filename = paste0(raster.names[[i]],".tif"),overwrite = TRUE)
  # format = "GTiff",na.rm=T, inf.rm=T, ) # Write to file
  #s # Return scaled raster
}
#   }
#stopCluster(cl)
rm(i,s,m,r)

### Summary
### 7 Climate variables
### 4 organisms variables
### 11 relief variables
### Consider reducing number of relief covariates to balance the number of variables by soil-forming factor
### Parent material is categorical - 9 classes
### Time
### In this case we don´t incorporate a time proxy

# 3. Categorical variables -----------------------------------------------

### Parent material
### Until we have the map from BRGM we use the 1:1M parent material map
mat11 <- rast("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/ParentMaterial/mat11.tif")
mat11 <- project(mat11, "EPSG:4326", method="near")
### Reclass 0 as NA
f <- function(x) ifelse(x == 0, NA, x)
mat11 <- app(mat11, f)
writeRaster(mat11, filename="C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/ParentMaterial/mat11_4326.tif")
mat11 <- rast("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/ParentMaterial/mat11_4326.tif")

### Transform to continuous varible with a PCA
### Create dummy variables
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/ParentMaterial/")
PM.classes <- sort(unique(values(mat11)))
dummy.list <- list()
for(i in 1:length(PM.classes)){
  print(i)
  dummy.pm <- app(mat11, 
                  fun = function(x) {ifelse(is.na(x), NA, ifelse(x==PM.classes[[i]],1,0))},
                  filename=paste0("MAT11.dummy",PM.classes[[i]],".tif"),
                  overwrite=TRUE)
  dummy.list[[i]] <- dummy.pm
}

dummy.files <- list.files(pattern="dummy")
dummy.rast <- rast(dummy.files)
rm(dummy.files)
plot(dummy.rast)
### change names
names(dummy.rast) <- paste0("mat_",PM.classes)

### Perform PCA
### First, let's get a hint of how many components we need
Npixels <- values(dummy.rast[[1]])
Npixels <- Npixels[!is.na(Npixels)]
N.sample <- length(Npixels)*0.1 ### 10% of pixels aprox

set.seed(1946)
# Regular sampling
sampleMAT11<- spatSample(dummy.rast, size = 6000000 , method="regular",replace=FALSE, as.df=TRUE, xy=TRUE, na.rm = TRUE)
### select only complete cases
sampleMAT11 <-sampleMAT11[complete.cases(sampleMAT11),]
sampleMAT11 <- as.data.frame(sampleMAT11)
summary(sampleMAT11)

### Apply scaled PCA
MAT11.pca <- prcomp(sampleMAT11[,3:ncol(sampleMAT11)], scale=TRUE)
library(factoextra)
fviz_eig(MAT11.pca)
fviz_pca_var(MAT11.pca,axes = c(1,2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )

# Eigenvalues
eig.val <- get_eigenvalue(MAT11.pca)
eig.val ### Let's retain ~ 100 % variability with 8 PCs

### Try the rasterPCA function from the RStool package
pca.pred <- predict(dummy.rast, MAT11.pca, index=1:8)
plot(pca.pred)
### Write them to file 
for (i in 1:nlyr(pca.pred)){
  writeRaster(pca.pred[[i]], filename = paste0("mat11_PC",i,".tif"), overwrite=TRUE)
}

mat11.files <- list.files(pattern="_PC")
mat11.rast <- rast(mat11.files)
plot(mat11.rast)

### Clean environment

# 4. Regular sample of all covariates --------------------------------------

### Start the script from here for modeling pedogenon classes. Covariates processed already

### We load scaled continuous variables
ScaledContDir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/scaled/"
setwd(ScaledContDir)

### My preferred order for covariates
covariates.sel <- c("Bio1.tif","Bio3.tif","Bio4.tif","Bio8.tif","Bio9.tif","Bio12.tif","Bio15.tif", ### climate
                    "forest_cover.tif","conv_rangelands.tif","grazing.tif","cropland.tif",          ### organisms
                    "cti.tif","easterness.tif","hli.tif","mrrtf.tif","mrvbf.tif","north_slope.tif",
                    "northerness.tif","scale_pos.tif","slope.tif","srr.tif","srtm.tif")             ### relief  
cont.stack <- rast(covariates.sel)

### Categorical variables
ParentMatDir <-"C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/ParentMaterial/"
setwd(ParentMatDir)
cat.stack <- rast(list.files(pattern="mat11_PC"))

### Make stack
covariates.stack <- c(cont.stack,cat.stack);plot(covariates.stack)

### For searching the optimal number of clusters with Silhouette index 
### I sample an even smaller number of pixels.

# Set the random number generator to reproduce the results
set.seed(1946)
# Regular sampling
sampleCLORPT <- spatSample(covariates.stack, size = 40000 , method="regular", replace=FALSE, as.df=TRUE, xy=TRUE, na.rm = TRUE)
  
### Returns a dataframe
# transform to sf
# sampleCLORPT <- st_as_sf(sampleCLORPT)
# ### Transform into a dataframe
# CLORPT.df <- st_drop_geometry(sampleCLORPT)

### select only complete cases
CLORPT.df <-sampleCLORPT[complete.cases(sampleCLORPT),]
dim(CLORPT.df)
## Clean 
rm(sampleCLORPT)
### The 40,000 sample for running silhouette index, etc.

### For searching the optimal number of clusters I sample a smaller number of pixels.
### I repeat with more pixels, now that I can perform it in parallel in a day, and I change the seed
# Set the random number generator to reproduce the results
set.seed(1946)
# Regular sampling
sampleCLORPT <- spatSample(covariates.stack, size = 400000 , method="regular", replace=FALSE, as.df=TRUE, xy=TRUE, na.rm = TRUE)

### Returns a dataframe
# transform to sf
# sampleCLORPT <- st_as_sf(sampleCLORPT)
# ### Transform into a dataframe
# CLORPT.df <- st_drop_geometry(sampleCLORPT)

### select only complete cases
CLORPT.df <-sampleCLORPT[complete.cases(sampleCLORPT),]
dim(CLORPT.df)
## Clean 
rm(sampleCLORPT)


# 5. Apply Cholesky decomposition to sample data --------------------------

# The basic Euclidean distance treats each variable as equally important in calculating the distance.
# An alternative approach is to scale the contribution of individual variables to the distance value according
# to the variability of each variable. This approach is illustrated by the Mahalanobis distance, 
# which is a measure of the distance between each observation in a multidimensional cloud of points and
# the centroid of the cloud.

### Calculate the Mahalanobis distance, as Euclidean distance after applying the Cholesky decomposition
## Take out the coordinates
CLORPT.df.coords <- CLORPT.df[,1:2]
CLORPT.df.subset <- CLORPT.df[,3:ncol(CLORPT.df)]

# Rescale the data 
C <- chol(var(as.matrix(CLORPT.df.subset)))
CLORPT.rs <- as.matrix(CLORPT.df.subset) %*% solve(C)

### euclidean distance of the first 10
b <- dist(CLORPT.rs[1:10,],method = "euclidean")
library(biotools)
a <- D2.dist(CLORPT.df.subset[1:10,], cov=var(CLORPT.df.subset))
sqrt(a);b ### distances are the same
### Clean
rm(a,b)


# ### 5. Optimal number of clusters ---------------------------------------

### 5. Optimal number of clusters

### Skip this step if you already have an optimal number

### how can we calculate the optimal number of clusters?
### there are several methods.
### A first one is to check the sum of within-cluster distances across all clusters
### Also, the package NbClust can calculate several indices to choose the optimal k

### With the package ClusterR
library(ClusterR)
search_space <- c(seq(from=10, to=170, by=10)) ### Let's say that we want to check between 10 to 170 clusters
set.seed(1991)
system.time(opt_kmeans <- Optimal_Clusters_KMeans(data = CLORPT.rs, 
                                                  max_clusters = search_space,
                                                  criterion = "WCSSE", 
                                                  num_init = 10,
                                                  max_iters = 10000, 
                                                  initializer = "kmeans++", 
                                                  plot_clusters = TRUE,
                                                  verbose = TRUE))
plot(search_space, opt_kmeans,pch=20, col="blue")
lines(search_space, opt_kmeans)
### As always, the elbow method is not very informative


### An index that combines the within cluster similarity and the between cluster dissimilarity is the Calinski-Harbasz index
### Calculate the optimal number of k-means clusters with NbClust package
library(NbClust)
set.seed(1984)
system.time(opt.Clusters.CH <- NbClust(data = CLORPT.rs, 
                                       diss=NULL, 
                                       distance = "euclidean",
                                       min.nc = 10, 
                                       max.nc = 170,
                                       method = "kmeans", 
                                       index = "ch"))
# user   system  elapsed 
# 3697.02  1832.89 62257.79 

opt.Clusters.CH
summary(opt.Clusters.CH)
opt.Clusters.CH$Best.nc
plot(10:170, opt.Clusters.CH$All.index,pch=20, col="blue")
lines(10:170, opt.Clusters.CH$All.index)
abline(v=10, col="red", lty=2)

### the silhouette method is also very popular (check here)
### this is taking too long so I will do it in parallel
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/NumberPedogenons")
library(foreach)
library(parallel)
library(iterators)
library(doParallel)
library(NbClust)

detectCores()
cl <- makeCluster(8)   ### Create cluster
registerDoParallel(cl)
getDoParWorkers()
search_spaces <- list(c(2,20),c(21,30),c(31,40),c(41,50),
                      c(51,60), c(61,70), c(71,80),
                       c(81,90), c(91,100), c(101,110), 
                      c(111,120),c(121,130),c(131,140),
                       c(141,150),c(151,160),c(161,170))
system.time(
  silhhouette_list <- foreach (i=1:length(search_spaces),.packages=c("NbClust"),
                               .export = c("CLORPT.rs","search_spaces")) %dopar% {
  set.seed(1991)
  #set.seed(1984)
  opt.Clusters.silhouette <- NbClust(data = CLORPT.rs, 
                                     diss=NULL, 
                                     distance = "euclidean",
                                     min.nc = search_spaces[[i]][1], 
                                     max.nc = search_spaces[[i]][2],
                                     method = "kmeans",
                                     index = "silhouette")
  gc()
  # Save the object as we go
  save(opt.Clusters.silhouette, file=paste0("opt.Clusters.silhouette50s_",i,".RData"))
  opt.Clusters.silhouette # We return this
}
)
stopCluster(cl)

### Extract the Silhouette index for each k of the search space
all.index <- c(silhhouette_list[[1]]$All.index,silhhouette_list[[2]]$All.index,silhhouette_list[[3]]$All.index,
               silhhouette_list[[4]]$All.index,silhhouette_list[[5]]$All.index,silhhouette_list[[6]]$All.index,
               silhhouette_list[[7]]$All.index,silhhouette_list[[8]]$All.index,silhhouette_list[[9]]$All.index,
               silhhouette_list[[10]]$All.index,silhhouette_list[[11]]$All.index,silhhouette_list[[12]]$All.index,
               silhhouette_list[[13]]$All.index,silhhouette_list[[14]]$All.index,silhhouette_list[[15]]$All.index,
               silhhouette_list[[16]]$All.index)
plot(2:170, all.index,pch=20, col="blue")
lines(2:170, all.index)
abline(v=31, col="red", lty=2)


### My parameters for running the clustering

### Apply k-means with foreach (if you want to create several maps)
# setwd(OutDir)
# detectCores()
# cl <- makeCluster(3)   ### Create cluster
# registerDoParallel(cl)
# getDoParWorkers()
#  K <- c(5,6,12)
# system.time(kmeansCLORPT <- foreach (i=1:length(K),.packages=c("ClusterR"),.export = c("CLORPT.rs","K")) %dopar% {
#   set.seed(1991+i)
#   kmeans_clorpt <- KMeans_rcpp(CLORPT.rs, clusters = K[[i]], num_init = 10, max_iters = 5000,
#                                  fuzzy = TRUE, initializer = 'kmeans++', verbose = T)
#   
#   # Save the object as we go
#   save(kmeans_clorpt, file=paste0(OutDir,"kmeans_clorpt.k",K[[i]],".RData"))
#   kmeans_clorpt # We return this
# })
# 
# stopCluster(cl)

### In this case, we choose 31 classes as optimum, identified by the Silhouette index

# perform KMeans_rcpp clustering
library(ClusterR)
my_seed <- 4587 ### Your set.seed() number
FR1000BCE.pedogenon <-KMeans_rcpp(data=CLORPT.rs, 
                                  clusters=31,
                                  num_init = 10, 
                                  max_iters = 20000,
                                  initializer = "kmeans++",
                                  fuzzy = FALSE, 
                                  verbose = TRUE,
                                  seed = my_seed)

### save the cluster number in the original dataframe
# CLORPT.df$Cluster.N <- as.factor(FR1000BCE.pedogenon$clusters)


# 6. Create pedogenon map ----------------------------------------------------

# 1. Mapping the clusters and write layers -------------------------------
### Previously, the stack with the scaled covariates were at
covariates.stack
plot(covariates.stack)

### Can terra apply a matricial multiplication directly?

tile.df.rs <- covariates.stack %*% solve(C)


### Predict the k class with a nested foreach loop 

### Extract the index of the centroids that are na/nan/Inf
Kcent.nan <- which(apply(FR1000BCE.pedogenon[["centroids"]], MARGIN = 1, FUN = function(y){any(is.na(y))}))
### None are NA

### Define the size of the blocks --- At each raster row do we start and finish each crop?
bs <- blocks(covariates.stack)
### I crop all raster files (across variables stacks) in a parallel process
### with a %dopar% from the foreach package
### this is thought for a large study area, but of course for the example is not needed

## My desired cluster number
K <- 31
#K <- c(6,12)
### If K is more than one instead of km.pedogenon.rcpp being a kmeans model, it would be a list of models

system.time(
  
  for(m in 1:length(K)){
    
    cl <- makeCluster(6)   
    registerDoParallel(cl)
    
    k_rast_list <- foreach(i=1:bs$n, .packages=c("terra", "ClusterR"), .export = c("C", "FR1000BCE.pedogenon")) %dopar% {
      
      ### Get one tile of the raster stack
      tile <- crop(covariates.stack, ext(covariates.stack, bs$row[[i]], bs$row[[i]]+bs$nrows[[i]], 1, ncol(covariates.stack)))
      ### Transform into a dataframe
      tile.df <- as.data.frame(tile, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=FALSE, long=FALSE)
      
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
      tile.df.rs[-df.na, ]$cluster  <- predict_KMeans(data = tile.df.rs[-df.na,1:(ncol(tile.df.rs)-1)], CENTROIDS = km.pedogenon.rcpp$centroids)
      ### Assign the values to a new raster
      k.pred <- setValues(tile[[1]], tile.df.rs$cluster)
      names(k.pred) <- paste0("K",K[[m]])
      k.pred # Return this
    }
    
    stopCluster(cl)
    
    ## Assign function to mosaic
    k_rast_list$na.rm <- TRUE
    k_rast_list$fun <- min
    ## Create mosaic for whole NSW
    k.raster <- do.call(mosaic, k_rast_list)
    names(k.raster) <- paste0("K",K[[m]])
    
    ## Write to file
    #writeRaster(k.raster, filename= paste0("K",K[[m]],".tif"), na.rm=T,inf.rm=T, format="GTiff", overwrite=TRUE )
    gc()
    
  }
)




























#### end of this script