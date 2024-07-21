#### Pedogenon maps for France

### Date 10/11/2023
### Author: Mercedes Roman Dobarco
### Objective: 1. Check the correlation and 
###            2. evolution of Relief for France,
###            3. Perform PCA 

# analysis packages
library(dplyr)
library(tidyverse)
library(sf)
library(terra)
library(Hmisc)
library(FactoMineR)
library(factoextra)
library(ClusterR)
library(dendextend)
library(dendsort)

# visualization packages
library(ggplot2)
library(gganimate)
library(patchwork)
library(viridis)
library(scales)
library(corrplot)
library(GGally)

### 1. Study area and-----------------------------------------------------------

### I bring shapefile
france_regions <- read_sf("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Administrative/regions-20180101-shp/regions-20180101.shp")

### Only metropolitan France
france_regions <- france_regions[france_regions$code_insee %in% c(11:93),]
plot(france_regions["nom"])
france <- st_union(france_regions)
### Project to RGF93 v2b / Lambert-93
franceL93 <- st_transform(france,9794)
france_buffer_WGS84 <- st_read("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Administrative/france_bufferWGS84.shp")

### Load Soil type map
soil <- terra::rast("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Soil/soil1.tif")

# ### Load DEM from France (same as GSM AWC products) to get extent
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Relief")
relief <- list.files(pattern="tif$")
relief <- c("cti.tif","curv_long.tif","curv_trans.tif","curvature.tif",
            "easterness.tif","exposition.tif",
            "hli.tif","linear_aspect.tif","mrrtf.tif",
            "mrvbf.tif","north_slope.tif","northerness.tif" ,
            "roughness.tif","sar.tif","scale_pos.tif",
            "slope.tif","slopeascos.tif","slopeassin.tif","slopeastrasp.tif",
            "srr.tif","srtm.tif")
reliefrast <- terra::rast(relief)
plot(reliefrast)

relief.names <- gsub(x=relief, pattern=".tif", replacement="")
names(reliefrast) <- relief.names

# dem <- rast("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Relief/srtm.tif")
# st_crs(dem)
# dem <- terra::project(dem, "EPSG:4326")
# plot(dem)
# ext(dem)
# 
# ### Create a mask from DEM
# demmask <- function(x) {ifelse(!is.na(x),1,NA)}
# frmask <- terra::app(dem, demmask, cores=6)
# 
# ### Crop because we have too much space in the south
# frmask <- terra::trim(frmask)
# plot(frmask)
# gc()
# ext(frmask)


# 1. Prediction ability of Relief covariates for soil type ----------------

SoilReliefS <- c(soil,reliefrast)
### Take regular sample
set.seed(2233)
SoilReliefSample <- terra::spatSample(x = SoilReliefS, 
                                  size=400000,
                                  method="regular", 
                                  as.df=TRUE, 
                                  xy=TRUE)
### Only complete cases
SoilReliefSample <- SoilReliefSample[complete.cases(SoilReliefSample),]
dim(SoilReliefSample)

### eliminate observations with Soil type = 0
SoilReliefSample <- SoilReliefSample[SoilReliefSample$soil1 != 0,]

### Trnasform into factor
SoilReliefSample$soil1 <- as.factor(SoilReliefSample$soil1)

### Exclude coordinates
SoilReliefSample <- SoilReliefSample[,3:ncol(SoilReliefSample)]
str(SoilReliefSample)

### Run random forest model
library(caret)
library(mlbench)
library(Hmisc)
library(randomForest)
library(ranger)

set.seed(71)
soil.rf <- ranger(soil1  ~ ., 
                  data=SoilReliefSample,
                  num.trees = 3000,
                  importance='impurity',
                  write.forest=FALSE)
print(soil.rf)
## Look at variable importance:
importanceRF <- sort(round(importance(soil.rf), 2))

impRF <- data.frame(variables=names(importanceRF) , importance=importanceRF)

ggplot(impRF) +
  geom_col(aes(y=fct_reorder(variables, importance),
               x=importance))

### Perform recursive feature elimination
set.seed(10)
ctrl <- rfeControl(functions = 'ranger',
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

rfProfile <- rfe(x=SoilReliefSample[,-1],
                 y =SoilReliefSample$soil1,
                 sizes = subsets,
                 rfeControl = ctrl)

lmProfile


### 2. Correlation plot -----------------

FrBf_WGS84 <- rast("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Administrative/FrBf_WGS84Mask.tif")

### Take regular sample
set.seed(2233)
ReliefSample <- terra::spatSample(x = reliefrast, 
                                   size=400000,
                                   method="regular", 
                                   as.df=TRUE, 
                                   xy=TRUE)
### Only complete cases
ReliefSample <- ReliefSample[complete.cases(ReliefSample),]
dim(ReliefSample)

### Which climate variables are correlated?
par(mfrow=c(1,1))
ReliefSample[,3:ncol(ReliefSample)] %>%
  cor(., use = "pairwise.complete.obs") %>%
  corrplot.mixed(.,upper = "ellipse", order="hclust",
                 lower = "number", 
                 number.cex=0.7, tl.cex=0.6, tl.col = "black")

### checking the correlations are significant
testRes = cor.mtest(ReliefSample[,3:ncol(ReliefSample)], conf.level = 0.95)
corr_bioclim = rcorr(as.matrix(ReliefSample[,3:ncol(ReliefSample)]))
M <- corr_bioclim$r
p_mat <- corr_bioclim$P

corrplot(M, 
         diag=FALSE,
         type = "upper", 
         order = "hclust", 
         method = "number",
         number.cex=0.6, 
         tl.cex=0.6,
         p.mat = p_mat, 
         sig.level = 0.05)

ReliefSample[,3:ncol(ReliefSample)] %>%
  cor(., use = "pairwise.complete.obs") %>%
  corrplot.mixed(.,upper = "ellipse", 
                 lower = "number", order="hclust",
                 number.cex=0.7, 
                 tl.cex=0.6, tl.col = "black")


### Highlight only those with correlation above a certain values
M0.70 <- M # Copy matrix
M0.70[ M0.70 < 0.70 & M > -0.70 ] = 0
corrplot(M0.70)
corrplot(M0.70, 
         diag=FALSE,
         type = "upper", 
         order = "hclust", 
         method = "number",
         number.cex=0.6, 
         tl.cex=0.6,
         p.mat = p_mat, 
         sig.level = 0.05)

rm(M,M0.70,p_mat,corr_bioclim)

# ### Create paired plot using GGally
# ### Subset candidate variables
# BioclimSubsetSample <- ReliefSample[,c("x","y","bio01","bio03","bio05", "bio08","bio09", "bio12", "bio15", "bio17")]
# corrplotsPairs <- GGally::ggpairs(data=BioclimSubsetSample)
# 
# ### Plot subset of variables based on correlation threshold 0.7
# BioclimSub <- terra::subset(bioclimrast, c(1,3,5,8,9,12,15,17))
# names(BioclimSub) <- c("Bio1", "Bio3", "Bio5", "Bio8",
#                        "Bio9", "Bio12", "Bio15", "Bio17")
# plot(BioclimSub)
# 
# plot(BioclimSub["Bio9"])
# plot(france, add=TRUE)


### 2.2 PCA on Bioclim variables ------------------------------------------------

### Perform a PCA to check variability of the data and check loadings
library("FactoMineR")
pca_relief <- FactoMineR::PCA(X = ReliefSample[,3:ncol(ReliefSample)], 
                               ncp=10,
                               scale.unit=TRUE,
                               graph=TRUE)

### The first 5 components retain almost 95% of the variance
eigenvalues <- pca_relief$eig
eigenvalues[1:10, 1:3]

### Check variance explained by each dimension
barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], 
      type="b", pch=19, col = "black")

### I can make the plot also with package factoextra
library(factoextra)
fviz_screeplot(pca_relief, ncp=10)

### Biplot with variable contribution to the pricipal components
plot(pca_relief, choix="var")
fviz_pca_var(pca_relief,col.var="contrib", axes = c(1,2))

# Contributions of variables to PC1
fviz_contrib(pca_relief, choice = "var", axes = 1, top = 10)

contrib.df <- pca_relief$var$contrib[,1:7]
#setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/Climate")
#write.csv(contrib.df, file="pca_climate.csv")

contrib.df %>%
  as.data.frame(.,) %>%
  rownames_to_column("Relief") %>%
  pivot_longer(-c(Relief), names_to = "PCs", values_to = "contrib") %>%
  mutate(PCs= fct_relevel(PCs,colnames(contrib.df ))) %>%
  ggplot(aes(y=PCs, x=Relief, fill=contrib)) + 
  geom_raster() + 
  scale_fill_viridis(direction=-1, option = "A",
                    trans = scales::pseudo_log_trans(sigma = 1))+
  theme(axis.text.x = element_text(angle = 45))

### Attach coordinates to ReliefSample
PCsScores <- as.data.frame(pca_relief$ind$coord)
ReliefSample <- cbind(ReliefSample,PCsScores)
### copy dataframe
PCsMapDF <- ReliefSample

### Separate maps
PC1map <- ggplot() +
  geom_sf(data = franceL93) +
  geom_point(aes(y = y, x = x, color=Dim.1), data=PCsMapDF) +
  scale_color_viridis(discrete = FALSE, option="D", direction = -1)
PC2map <- ggplot() +
  geom_sf(data = franceL93) +
  geom_point(aes(y = y, x = x, color=Dim.2), data=PCsMapDF) +
  scale_color_viridis(discrete = FALSE, option="D", direction = -1)
PC3map <- ggplot() +
  geom_sf(data = franceL93) +
  geom_point(aes(y = y, x = x, color=Dim.3), data=PCsMapDF) +
  scale_color_viridis(discrete = FALSE, option="D", direction = -1)


library(gridExtra)
grid.arrange(PC1map, PC2map, PC3map, PC4map, PC5map, PC6map, PC7map, ncol=4)
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/Climate")
PC1map; ggsave("PC1map.jpg")
PC2map; ggsave("PC2map.jpg")
PC3map; ggsave("PC3map.jpg")
PC4map; ggsave("PC4map.jpg")
PC5map; ggsave("PC5map.jpg")
PC6map; ggsave("PC6map.jpg")
PC7map; ggsave("PC7map.jpg")

# perform KMeans_rcpp clustering
library(ClusterR)
my_seed <- 11 ### Your set.seed() number
climatePCs_cl <-KMeans_rcpp(data=PCsMapDF[,c( "Dim.1","Dim.2","Dim.3","Dim.4")],
                            clusters=8,
                            num_init = 20,
                            max_iters = 5000,
                            initializer = "kmeans++",
                            fuzzy = FALSE,
                            verbose = FALSE,
                            seed = my_seed)

### save the cluster number in the original dataframe
PCsMapDF$PCs8k <- as.factor(climatePCs_cl$clusters)
ggmap(france_background) +
  geom_point(aes(y = y, x = x, color=PCs8k),
             cex=0.1, pch=19,
             data=PCsMapDF) +
  scale_color_viridis(discrete=TRUE,direction=-1, option="A")

library(dendsort) 
library(dendextend) # visualize dendrograms
library(colorspace) # colors

### Hierarchical clustering
hc <- hclust(dist(climatePCs_cl$centroids), method="ward.D2")
plot(dendsort(hc), main="Hierarchical clustering of 8000 BCE climate PCs", sub="", xlab="")

### Extract labels
hc.labels <- hc %>% as.dendrogram(.) %>% labels %>% as.numeric()

### Extract the membership from the tree
dend <- hc %>% as.dendrogram(.)
plot(dend)

### Colours
mypalette <- c(sequential_hcl("OrYel", n = 4),
               sequential_hcl("TealGrn", n = 2), 
               sequential_hcl("PurpOr", n = 2))
legend.plot <- dend %>%  
  set("labels_cex", 2) %>% 
  set("branches_lwd", 2)%>% 
  set("labels_col", mypalette) %>% 
  set("branches_k_color", mypalette)
plot(legend.plot)

### Order for map
cluster_colors <- data.frame(orderID = hc.labels, colors = mypalette)
palette_order <- cluster_colors %>% arrange(cluster_colors, orderID)

### Plot with hierarchical clustering
ggmap(france_background) +
  geom_point(aes(y = y, x = x, color=PCs8k),
             cex=0.1, pch=19,
             data=PCsMapDF) +
  scale_colour_manual(values=palette_order$colors)
  
centroids <- climatePCs_cl$centroids
dend %>% dendextend:::cutree.dendrogram(., k = 3)


### What id a different number of clusters?
library(ClusterR)
set.seed(1991)
search_space <- c(5:15) ### Let's say that we want to check between 2 to 10 clusters, because there are 8 climatic types
system.time(opt_kmeans <- Optimal_Clusters_KMeans(data = PCsMapDF[,c( "Dim.1","Dim.2","Dim.3","Dim.4")], 
                                                  max_clusters = search_space,
                                                  criterion = "WCSSE", num_init = 10,
                                                  max_iters = 5000, initializer = "kmeans++", plot_clusters = TRUE,
                                                  verbose = TRUE))
plot(search_space, opt_kmeans,pch=20, col="blue")
lines(search_space, opt_kmeans)

# perform KMeans_rcpp clustering
library(ClusterR)
my_seed <- 13 ### Your set.seed() number
climatePCs_cl <-KMeans_rcpp(data=PCsMapDF[,c( "Dim.1","Dim.2","Dim.3","Dim.4")],
                            clusters=9,
                            num_init = 20,
                            max_iters = 5000,
                            initializer = "kmeans++",
                            fuzzy = FALSE,
                            verbose = FALSE,
                            seed = my_seed)

### save the cluster number in the original dataframe
PCsMapDF$PCs9k <- as.factor(climatePCs_cl$clusters)
ggmap(france_background) +
  geom_point(aes(y = y, x = x, color=PCs8k),
             cex=0.1, pch=19,
             data=PCsMapDF) +
  scale_color_viridis(discrete=TRUE,direction=-1, option="A")

### 2.3. Perform clustering and check if other time steps fall in thes --------

### Scale variables with the cholesky transformation
### although the correlation is already acceptable
### But I don´t want to scale (SD = 1) so I can predict on the new climate variables

# Rescale the data - Only first 22 (numerical variables)
C <- chol(var(as.matrix(BioclimSubsetSample[,3:10])))
Bioclim.Cholesky <- as.matrix(BioclimSubsetSample[,3:10]) %*% solve(C)

### Optimal number of clusters
### how can we calculate the optimal number of clusters?
### there are several methods.
### A first one is to check the sum of within-cluster distances across all clusters
### Also, the package NbClust can calculate several indices to choose the optimal k

### With the package ClusterR
library(ClusterR)
set.seed(1991)
search_space <- c(2:20) ### Let's say that we want to check between 2 to 10 clusters, because there are 8 climatic types
system.time(opt_kmeans <- Optimal_Clusters_KMeans(data = Bioclim.Cholesky, 
                                                  max_clusters = search_space,
                                                  criterion = "WCSSE", num_init = 10,
                                                  max_iters = 5000, initializer = "kmeans++", plot_clusters = TRUE,
                                                  verbose = TRUE))
plot(search_space, opt_kmeans,pch=20, col="blue")
lines(search_space, opt_kmeans)

### Take smaller sample for this exercise
### Take regular sample
set.seed(2233)
gc()
ReliefSample2 <- terra::spatSample(x = bioclimrast, 
                                   size= 100000,
                                   method="regular", 
                                   as.df=TRUE, 
                                   xy=TRUE)
### Only complete cases
ReliefSample2 <- ReliefSample2[complete.cases(ReliefSample2),]
dim(ReliefSample2)

### change column names
colnames(ReliefSample2) <- c("x","y",biovars,"mask")
### Subset candidate variables
ReliefSample2 <- ReliefSample2[,c("x","y","bio01","bio03","bio05", "bio08","bio09", "bio12", "bio15", "bio17")]

## Scale
Bioclim2.Cholesky <- as.matrix(ReliefSample2[,3:10]) %*% solve(C)
gc()

### Calculate the optimal number of k-means clusters with NbClust package
library(NbClust)
set.seed(1984)
system.time(opt.Clusters.dunn <- NbClust(data = Bioclim2.Cholesky, 
                                         diss=NULL, 
                                         distance = "euclidean",
                                         min.nc = 2, 
                                         max.nc = 20,
                                         method = "kmeans",
                                         index = "dunn"))
gc()
opt.Clusters.dunn
summary(opt.Clusters.dunn)
opt.Clusters.dunn$Best.nc
plot(search_space, opt.Clusters.dunn$All.index,pch=20, col="blue")
lines(search_space, opt.Clusters.dunn$All.index)
abline(v=opt.Clusters.dunn$Best.nc, col="red", lty=2)
#### Last session was lost but I do remember 11 was the optimal number

### In this case, we choose 11 classes as optimum
# perform KMeans_rcpp clustering
my_seed <- 4587 ### Your set.seed() number
climate.10BP.rcpp <-KMeans_rcpp(data=Bioclim.Cholesky, 
                                clusters=11,
                                num_init = 10, 
                                max_iters = 5000,
                                initializer = "kmeans++", 
                                fuzzy = FALSE, 
                                verbose = FALSE,
                                seed = my_seed)

### save the cluster number in the original dataframe
BioclimSubsetSample$K10BP <- as.factor(climate.10BP.rcpp$clusters)

### Plot clusters
ggmap(france_background) +
  geom_point(aes(y = y, x = x, color=K10BP),
             cex=0.1, pch=19,
             data=BioclimSubsetSample) +
  scale_color_viridis(discrete=TRUE,direction=-1, option="A")

### In this case, we choose 11 classes as optimum
# perform KMeans_rcpp clustering
my_seed <- 46637 ### Your set.seed() number
climate.10BP.rcpp <-KMeans_rcpp(data=Bioclim.Cholesky, 
                                clusters=8,
                                num_init = 10, 
                                max_iters = 5000,
                                initializer = "kmeans++", 
                                fuzzy = FALSE, 
                                verbose = FALSE,
                                seed = my_seed)

### save the cluster number in the original dataframe
BioclimSubsetSample$K10BP <- as.factor(climate.10BP.rcpp$clusters)

### Plot clusters
ggmap(france_background) +
  geom_point(aes(y = y, x = x, color=K10BP),
             cex=0.1, pch=19,
             data=BioclimSubsetSample) +
  scale_color_viridis(discrete=TRUE,direction=-1, option="A")


### 2.4 Download and crop climate covariates for other time steps --------------

get_chelsa_paleo_france_crop <- function(climdir, timeIDs, vars, ExtD){
  timeout_old <- getOption('timeout')
  options(timeout=1000000)
  
  for(timeID in 1:length(timeIDs)){
    
    for(var in 1:length(vars)){
      
      ### Download file
      name <- paste0("CHELSA_TraCE21k_",vars[[var]],"_",timeIDs[[timeID]], "_V1.0.tif")
      source_url <- file.path(paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_trace/bio/",name))
      destination <- file.path(paste0(climdir,"TimeID_",timeIDs[[timeID]],"/",name))
      if(!file.exists(destination)){
        download.file(source_url, destination, method="wget")
        } else { 
        message(paste0(destination, " already downloaded, skipping to next"))
      }
    }
    options(timeout = timeout_old)
    
    ### Process raster files for this time step
    setwd(paste0(climdir,"TimeID_",timeIDs[[timeID]],"/"))
    print(paste0("Cropping and storing rasters for France in ",paste0(climdir,"TimeID_",timeIDs[[timeID]],"/")))
    ### List tif files
    bioclim <- list.files(pattern=".tif")
    
    ### Crop to the extent of France (keep same alignment and resolution) 
    for(r in 1:length(bioclim)){
      
      ### Load raster
      bioclimrast <- terra::rast(paste0(climdir,"TimeID_",timeIDs[[timeID]],"/",bioclim[[r]]))
      ### Crop to the extent of France
      bioclimrastc <- terra::crop(bioclimrast, ExtD)
      ### Write Raster
      writeRaster(bioclimrastc, filename = paste0(names(bioclimrastc),"_Fr.tif"))
      ### Remove raster global extent
      file.remove(paste0(climdir,"TimeID_",timeIDs[[timeID]],"/",bioclim[[r]]))
      
    } 
  }
}

# Bring new layers using the function
setwd(climDir)
### List of time steps
timez <- as.character(c(1:20))
### List of bioclim variables
### Only for the subset I selected
biovars <- paste0("bio", c("01","03","05","08","09","12","15","17"))
### Extent to crop
myext <- ext(st_bbox(france_buffer_WGS84))

get_chelsa_paleo_france_crop(climdir = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Climate/",
                             timeIDs = timez, 
                             vars = biovars,
                             ExtD = myext)

get_chelsa_paleo_france_crop(climdir = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Climate/",
                             timeIDs = "-76", 
                             vars = biovars,
                             ExtD = myext)

### 2.5 Scaled variables -----------------------------------------------------

### Mask oceans (in this case, use buffer around France)
### Subset bioclim variables + mask
BioclimSub <- terra::subset(bioclimrast, c(1,3,5,8,9,12,15,17,20))
names(BioclimSub) <- c("Bio1", "Bio3", "Bio5", "Bio8",
                       "Bio9", "Bio12", "Bio15", "Bio17", "mask")
target.vars <- c("Bio1", "Bio3", "Bio5", "Bio8",
                 "Bio9", "Bio12", "Bio15", "Bio17")

### Mask and replace in the same raster stack
for(i in 1:length(target.vars)){
  BioclimSub[[i]] <- terra::mask(BioclimSub[[i]], BioclimSub[["mask"]])
}

plot(BioclimSub)

### Calculate mean and sd of 8k BCE bioclimatic variables
rs8kBCE <- data.frame(mean=c(rep(NA,8)), sd=c(rep(NA,8)))

for(i in 1:length(target.vars)){
  rs8kBCE[i,1] <- mean(values(BioclimSub[[i]]),na.rm=TRUE)
  rs8kBCE[i,2] <- sd(values(BioclimSub[[i]]),na.rm=TRUE)
  ### Learnt that this function exists
  ### global(x, c("sum", "mean", "sd"), na.rm=TRUE)
}

### Scale 8k BCE bioclimatic variables
Bioclim.sc <- list()
for(i in 1:length(target.vars)){
  Bioclim.sc[[i]] <- terra::scale(BioclimSub[[i]], center=TRUE, scale=TRUE)
}
Bioclim.sc <- rast(Bioclim.sc)

### Take smaller sample for this exercise

### Take regular sample
set.seed(2233)
Bioclim.sc.Samp <- terra::spatSample(x = Bioclim.sc, 
                                    size= 40000,
                                    method="regular", 
                                    as.df=TRUE, 
                                    xy=TRUE)
### Only complete cases
Bioclim.sc.Samp <- Bioclim.sc.Samp[complete.cases(Bioclim.sc.Samp),]
dim(Bioclim.sc.Samp)

### Calculate the optimal number of k-means clusters with NbClust package
library(NbClust)
set.seed(1984)
system.time(opt.Clusters.dunn <- NbClust(data = Bioclim.sc.Samp[,3:10], 
                                         diss=NULL, 
                                         distance = "euclidean",
                                         min.nc = 2, 
                                         max.nc = 20,
                                         method = "kmeans",
                                         index = "dunn"))
gc()
opt.Clusters.dunn
summary(opt.Clusters.dunn)
opt.Clusters.dunn$Best.nc
plot(2:20, opt.Clusters.dunn$All.index,pch=20, col="blue")
lines(2:20, opt.Clusters.dunn$All.index)
abline(v=opt.Clusters.dunn$Best.nc, col="red", lty=2)
#### Optimum number is 9, but results differ so much...

# perform KMeans_rcpp clustering
my_seed <- 2245 ### Your set.seed() number
climate.10BP.rcpp <-KMeans_rcpp(data=Bioclim.sc.Samp[,3:10], 
                                clusters=9,
                                num_init = 10, 
                                max_iters = 5000,
                                initializer = "kmeans++", 
                                fuzzy = FALSE, 
                                verbose = FALSE,
                                seed = my_seed)

### save the cluster number in the original dataframe
Bioclim.sc.Samp$K10BP <- climate.10BP.rcpp$clusters

### Plot clusters
library(ggmap)
ggmap(france_background) +
  geom_point(aes(y = y, x = x, color=K10BP),
             cex=1, pch=19,
             data=Bioclim.sc.Samp) +
  scale_color_viridis(direction=-1, option="A")


### 2.6 Predict clusters from 10k BP in every time step ------------------

### function to extract coordinates, scale, 
### and predict cluster from 8k BCE in new data
### @climdir the directory whwre the paleoclimate variables for that time step are stored
### @timeID time step of interest
### @vars Bioclim variables of interest, name of columns in the dataframe
### @dfScale dataframe with mean and sd of the Bioclimatic variables from the first time step
### @SampleDF the dataframe with the sample that we use for clustering 
### the first time step (8000 BCE in this case)
### @kmeansM the kmeans model (package ClusterR) that we use for predicting

### The function returns the predictions in the new time step. 
### We can attach to the dataframe SampleDF as output

# save.image("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/Rsessions/22112023.RData")
# load("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/Rsessions/22112023.RData")

assign_8kBCE_cluster_timeSteps <- function(climdir, timeID, 
                                           vars, dfScale,
                                           SampleDF, kmeansM){
  
    ### Process raster files for this time step
    setwd(paste0(climdir,"TimeID_",timeID,"/"))
   
    print(paste0("Extracting bioclim variables for ",paste0("TimeID = ",timeID)))
    
    ### List tif files
    bioclim <- list.files(pattern="Fr.tif")
    
    ### Load raster stack
    bioclimrast <- terra::rast(paste0(climdir,"TimeID_",timeID,"/",bioclim))
    
    ### Scale with the mean and sd from 8k BCE
    Bioclim.sc <- terra::scale(bioclimrast, 
                               center = dfScale[,"mean"],
                               scale = dfScale[,"sd"])
    
    ### Take regular sample in same coordinates as for 8000 BCE
    extractCoords <- as.matrix(SampleDF[,c("x", "y")])
    extraction <- terra::extract(Bioclim.sc, extractCoords)
    
    ### Rename columns with vars
    colnames(extraction) <- vars
    
    ## Create empty prediction column
    SampleDF$cluster <- NA
   
    ### Predict assignment in new raster
    SampleDF$cluster <- predict_KMeans(data = extraction, CENTROIDS = kmeansM$centroids)
    
    ### Rename cluster for that time step
    colnames(SampleDF)[colnames(SampleDF)=="cluster"] <- paste0("cl",timeID)
    
    ### Return dataframe, extraction df, and
    return(SampleDF)
    
    } 


### Apply function in for loop
timez <- as.character(c(-79:20))

for(i in 1:length(timez)) {
  
  if(i == 1) {
    print(timez[[i]])
    ### Apply function for each time step
    out <- assign_8kBCE_cluster_timeSteps(climdir = climDir, 
                                          timeID = timez[[i]], 
                                          vars = c("Bio1", "Bio3", "Bio5", "Bio8", "Bio9", "Bio12", "Bio15", "Bio17"),
                                          dfScale = rs8kBCE,
                                          SampleDF = Bioclim.sc.Samp,
                                          kmeansM = climate.10BP.rcpp)
    
  } else if (i > 1) {
    
    print(timez[[i]])
    
    ### Apply function for each time step
    out <- assign_8kBCE_cluster_timeSteps(climdir = climDir, 
                                          timeID = timez[[i]], 
                                          vars = c("Bio1", "Bio3", "Bio5", "Bio8", "Bio9", "Bio12", "Bio15", "Bio17"),
                                          dfScale = rs8kBCE,
                                          SampleDF = out,
                                          kmeansM = climate.10BP.rcpp)
  }
}


### Plot clusters
ggmap(france_background) +
  geom_point(aes(y = y, x = x, color=`cl-79`),
             cex=1, pch=19,
             data=out) +
  scale_color_viridis(direction=-1, option="A")

### Plot evolution over several time steps (not all of them)

### change column names
colnames(out) <- c("x","y", target.vars, paste0("timeID_",-80:20))

p1 <- ggmap(france_background) +
  geom_point(aes(y = y, x = x, color=as.factor(`timeID_-80`)),
             cex=1, pch=19,
             data=out) +
  scale_color_viridis(discrete=TRUE,direction=-1, option="A")

p2 <- ggmap(france_background) +
  geom_point(aes(y = y, x = x, color=as.factor(`timeID_0`)),
             cex=1, pch=19,
             data=out) +
  scale_color_viridis(discrete=TRUE,direction=-1, option="A")

p3 <- ggmap(france_background) +
  geom_point(aes(y = y, x = x, color=as.factor(`timeID_7`)),
             cex=1, pch=19,
             data=out) +
  scale_color_viridis(discrete=TRUE,direction=-1, option="A")

library(gridExtra)
grid.arrange(p1,p2,p3, ncol=3)


### Hierarchical clustering
hc <- hclust(dist(climate.10BP.rcpp$centroids), method="ward.D2")
par(mfrow=c(1,1))
plot(dendsort(hc), main="Hierarchical clustering of 8000 BCE bioclim", sub="", xlab="")

### Extract labels
hc.labels <- hc %>% as.dendrogram(.) %>% labels %>% as.numeric()

### Extract the membership from the tree
dend <- hc %>% as.dendrogram(.)
plot(dend)

### Colours
library(colorspace)
mypalette <- c(sequential_hcl("OrYel", n = 5),
               sequential_hcl("TealGrn", n = 4))
legend.plot <- dend %>%  
  set("labels_cex", 2) %>% 
  set("branches_lwd", 2)%>% 
  set("labels_col", mypalette) %>% 
  set("branches_k_color", mypalette)
plot(legend.plot)

### Order for map
cluster_colors <- data.frame(orderID = hc.labels, colors = mypalette)
palette_order <- cluster_colors %>% arrange(cluster_colors, orderID)

### Plot with hierarchical clustering
ggmap(france_background) +
  geom_point(aes(y = y, x = x, color=as.factor(`timeID_-80`)),
             cex=1, pch=19,
             data=out) +
  scale_colour_manual(values=palette_order$colors)

ggmap(france_background) +
  geom_point(aes(y = y, x = x, color=as.factor(`timeID_0`)),
             cex=1, pch=19,
             data=out) +
  scale_colour_manual(values=palette_order$colors)

ggmap(france_background) +
  geom_point(aes(y = y, x = x, color=as.factor(`timeID_20`)),
             cex=1, pch=19,
             data=out) +
  scale_colour_manual(values=palette_order$colors)


### Pivor to long
cluster.long <- out %>%
  pivot_longer(., cols=11:111, names_to = "timeID", values_to = "cluster"  )

cluster.summary <- cluster.long %>%
  group_by(., timeID, cluster) %>%
  summarise(., N=n()) 

cluster.summary$timeStep <- as.numeric(gsub(x=cluster.summary$timeID, 
                                 pattern= "timeID_",
                                 replacement = ""))

ggplot()+
  geom_point(aes(x=timeStep, y=N, color=as.factor(cluster)), 
             data=cluster.summary) +
  geom_smooth(aes(x=timeStep, y=N, color=as.factor(cluster)), 
             data=cluster.summary, se=FALSE) +
  scale_color_manual(values=palette_order$colors)

save.image("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/Rsessions/22112023.RData")

# ### 3.1 Analyses on 900-1000 BCE ---------------------

### Apply to time step -9
#biovars <- paste0("bio", c("02","04","06","07","10","11","13","14","16","18","19"))
### get extent
myext <- ext(st_bbox(france_buffer_WGS84))
biovars <- paste0("bio", c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19"))

### I had selected the bioclim variables based on 8000 BCE correlations.
### I repeat the analysis with all bioclimatic variables
get_chelsa_paleo_france_crop(climdir = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Climate/",
                             timeIDs = "-9", 
                             vars = biovars,
                             ExtD = myext)

### Make raster stack with climate 1000 BCE
setwd(paste0(climDir,"TimeID_-9/"))
bioclim <- list.files(pattern="Fr.tif")
bioclimrast <- terra::rast(paste0(climDir,"TimeID_-9/",bioclim))
plot(bioclimrast)
FrBf_WGS84 <- rast("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Administrative/FrBf_WGS84Mask.tif")

### I attach the mask
bioclimrast <- c(bioclimrast,FrBf_WGS84)

### Take regular sample
set.seed(2233)
ReliefSample <- terra::spatSample(x = bioclimrast, 
                                   size=700000,
                                   method="regular", 
                                   as.df=TRUE, 
                                   xy=TRUE)
### Only complete cases
ReliefSample <- ReliefSample[complete.cases(ReliefSample),]
dim(ReliefSample)

### change column names
colnames(ReliefSample) <- c("x","y",biovars,"mask")

### What is the distance between pixels?
# sampledists <- distance(lonlat=TRUE, as.matrix(ReliefSample[400000:401000,c("x","y")]))
# sampledists <- as.matrix(sampledists)
# m2 <- subset(melt(sampledists), value!=0)
# m2 %>% arrange(.,value) %>% slice_head(., n=10) ~ 930m

### Which climate variables are correlated?
par(mfrow=c(1,1))
ReliefSample[,3:21] %>%
  cor(., use = "pairwise.complete.obs") %>%
  corrplot.mixed(.,upper = "ellipse", 
                 lower = "number",  order="hclust",
                 number.cex=0.7, tl.cex=0.6, tl.col = "black")

### checking the correlations are significant
testRes = cor.mtest(ReliefSample[,3:21], conf.level = 0.95)
corr_bioclim = rcorr(as.matrix(ReliefSample[,3:21]))
M <- corr_bioclim$r
p_mat <- corr_bioclim$P

corrplot(M, 
         diag=FALSE,
         type = "upper", 
         order = "hclust", 
         method = "number",
         number.cex=0.6, 
         tl.cex=0.6,
         p.mat = p_mat, 
         sig.level = 0.05)

### Highlight only those with correlation above a certain values
M0.70 <- M # Copy matrix
M0.70[ M0.70 < 0.7 & M > -0.7 ] = 0
corrplot(M0.70)
corrplot(M0.70, 
         diag=FALSE,
         type = "upper", 
         order = "hclust", 
         method = "number",
         number.cex=0.6, 
         tl.cex=0.6,
         p.mat = p_mat, 
         sig.level = 0.05)
rm(M,M0.70,p_mat,corr_bioclim)

### Create paired plot using GGally
### Subset candidate variables
BioclimSubsetSample <- ReliefSample[,c("x","y","bio01","bio03","bio04", "bio08","bio09", "bio12", "bio15")]
corrplotsPairs <- GGally::ggpairs(data=BioclimSubsetSample)

### Plot subset of variables based on correlation threshold 0.7
BioclimSub <- terra::subset(bioclimrast, c(1,3,4,8,9,12,15))
names(BioclimSub) <- c("Bio1", "Bio3", "Bio4", "Bio8",
                       "Bio9", "Bio12", "Bio15")
plot(BioclimSub)
plot(BioclimSub["Bio9"])
plot(france, add=TRUE)

### Perform PCA

### Perform a PCA to check variability of the data and check loadings
library("FactoMineR")
pca_relief <- FactoMineR::PCA(X = ReliefSample[,3:21], 
                               ncp=10,
                               scale.unit=TRUE,
                               graph=TRUE)

save.image("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/Rsessions/01122023.RData")

### The first 5 components retain almost 95% of the variance
eigenvalues <- pca_relief$eig
eigenvalues[1:10, 1:3]

### Check variance explained by each dimension
barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], 
      type="b", pch=19, col = "black")

### I can make the plot also with package factoextra
library(factoextra)
fviz_screeplot(pca_relief, ncp=10)

### Biplot with variable contribution to the pricipal components
plot(pca_relief, choix="var")
fviz_pca_var(pca_relief,col.var="contrib")
plot(pca_relief, choix="var", axes=c(3,4))

# Contributions of variables to PC1
fviz_contrib(pca_relief, choice = "var", axes = 4, top = 10)

contrib.df <- pca_relief$var$contrib[,1:5]
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/Climate/1000BCE/")
write.csv(contrib.df, file="pca_climate1000BCE.csv")

contrib.df %>%
  as.data.frame(.,) %>%
  rownames_to_column("Bioclim") %>%
  pivot_longer(-c(Bioclim), names_to = "PCs", values_to = "contrib") %>%
  mutate(PCs= fct_relevel(PCs,colnames(contrib.df ))) %>%
  ggplot(aes(y=PCs, x=Bioclim, fill=contrib)) + 
  geom_raster() + 
  scale_fill_viridis(direction=-1, option = "A",
                     trans = scales::pseudo_log_trans(sigma = 1)) +
  theme(axis.text.x = element_text(angle = 45))

### Can I plot the principal components over France?
remotes::install_github("https://github.com/dkahle/ggmap")
library(ggmap)
register_stadiamaps("562bf9a7-9967-49e8-b6bd-eb15699dfdb9", write = TRUE)

### Check bbox
# myext
# france_background <- get_stadiamap(bbox=c(left=-5.2, bottom=42.2, right=8.4, top=51.2),  
#                                    maptype="stamen_terrain_background", zoom=6)
# ggmap(france_background)

### Attach coordinates to ReliefSample
PCsScores <- as.data.frame(pca_relief$ind$coord)
ReliefSample <- cbind(ReliefSample,PCsScores)
### copy dataframe
PCsMapDF <- ReliefSample

ggplot() +
  geom_sf(data = france) +
  geom_point(aes(y = y, x = x, color=Dim.4), data=PCsMapDF) +
  scale_color_viridis(discrete = FALSE, option="D", direction = -1)

ggplot() +
  geom_sf(data = france) +
  geom_point(aes(y = y, x = x, color=Dim.5), data=PCsMapDF) +
  scale_color_viridis(discrete = FALSE, option="D", direction = -1)


# ### Separate maps
# PC1map <- ggmap(france_background) +
#   geom_point(aes(y = y, x = x, color=Dim.1), data=PCsMapDF) +
#   scale_color_viridis(discrete = FALSE, option="D", direction = -1)
# PC2map <- ggmap(france_background) +
#   geom_point(aes(y = y, x = x, color=Dim.2), data=PCsMapDF) +
#   scale_color_viridis(discrete = FALSE, option="D", direction = -1)
# PC3map <- ggmap(france_background) +
#   geom_point(aes(y = y, x = x, color=Dim.3), data=PCsMapDF) +
#   scale_color_viridis(discrete = FALSE, option="D", direction = -1)
# PC4map <- ggmap(france_background) +
#   geom_point(aes(y = y, x = x, color=Dim.4), data=PCsMapDF) +
#   scale_color_viridis(discrete = FALSE, option="D", direction = -1)
# PC5map <- ggmap(france_background) +
#   geom_point(aes(y = y, x = x, color=Dim.5), data=PCsMapDF) +
#   scale_color_viridis(discrete = FALSE, option="D", direction = -1)
# PC6map <- ggmap(france_background) +
#   geom_point(aes(y = y, x = x, color=Dim.6), data=PCsMapDF) +
#   scale_color_viridis(discrete = FALSE, option="D", direction = -1)
# PC7map <- ggmap(france_background) +
#   geom_point(aes(y = y, x = x, color=Dim.7), data=PCsMapDF) +
#   scale_color_viridis(discrete = FALSE, option="D", direction = -1)
# 
# library(gridExtra)
# grid.arrange(PC1map, PC2map, PC3map, PC4map, PC5map, PC6map, PC7map, ncol=4)
# setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/Climate")
# PC1map; ggsave("PC1map.jpg")
# PC2map; ggsave("PC2map.jpg")
# PC3map; ggsave("PC3map.jpg")
# PC4map; ggsave("PC4map.jpg")
# PC5map; ggsave("PC5map.jpg")
# PC6map; ggsave("PC6map.jpg")
# PC7map; ggsave("PC7map.jpg")

# perform KMeans_rcpp clustering
library(ClusterR)
my_seed <- 1993 ### Your set.seed() number
climatePCs_cl <-KMeans_rcpp(data=PCsMapDF[,c( "Dim.1","Dim.2","Dim.3","Dim.4","Dim.5")],
                            clusters=8,
                            num_init = 5,
                            max_iters = 5000,
                            initializer = "kmeans++",
                            fuzzy = FALSE,
                            verbose = FALSE,
                            seed = my_seed)


library(ClusterR)
my_seed <- 1993 ### Your set.seed() number
climatePCs_cl2 <-KMeans_rcpp(data=PCsMapDF[,c( "Dim.1","Dim.2","Dim.3","Dim.4")],
                            clusters=8,
                            num_init = 5,
                            max_iters = 5000,
                            initializer = "kmeans++",
                            fuzzy = FALSE,
                            verbose = FALSE,
                            seed = my_seed)

### save the cluster number in the original dataframe
PCsMapDF$PC5_8k <- as.factor(climatePCs_cl$clusters)
PCsMapDF$PC4_8k <- as.factor(climatePCs_cl2$clusters)

ggplot() +
  geom_sf(data = france) +
  geom_point(aes(y = y, x = x, color=PC5_8k),
             cex=0.1, pch=19,
             data=PCsMapDF) +
  scale_color_viridis(discrete=TRUE,direction=-1, option="A")

ggplot() +
  geom_sf(data = france) +
  geom_point(aes(y = y, x = x, color=PC4_8k),
             cex=0.1, pch=19,
             data=PCsMapDF) +
  scale_color_viridis(discrete=TRUE,direction=-1, option="A")

library(dendsort) 
library(dendextend) # visualize dendrograms
library(colorspace) # colors

### Hierarchical clustering
hc <- hclust(dist(climatePCs_cl$centroids), method="ward.D2")
plot(dendsort(hc), main="Hierarchical clustering of 1000 BCE climate PCs", sub="", xlab="")

### Extract labels
hc.labels <- hc %>% as.dendrogram(.) %>% labels %>% as.numeric()

### Extract the membership from the tree
dend <- hc %>% as.dendrogram(.)
plot(dend)

### Colours
mypalette <- c(sequential_hcl("TealGrn", n = 2),
               sequential_hcl("OrYel", n = 4),
               sequential_hcl("PurpOr", n = 2))
legend.plot <- dend %>%  
  set("labels_cex", 2) %>% 
  set("branches_lwd", 2)%>% 
  set("labels_col", mypalette) %>% 
  set("branches_k_color", mypalette)
plot(legend.plot)

### Order for map
cluster_colors <- data.frame(orderID = hc.labels, colors = mypalette)
palette_order <- cluster_colors %>% arrange(cluster_colors, orderID)

### Plot with hierarchical clustering
ggplot() +
  geom_sf(data = france) +
  geom_point(aes(y = y, x = x, color=PC5_8k),
             cex=0.1, pch=19,
             data=PCsMapDF) +
  scale_colour_manual(values=palette_order$colors)

centroids <- climatePCs_cl$centroids
dend %>% dendextend:::cutree.dendrogram(., k = 3)
save.image("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/Rsessions/04122023.RData")


###################################################################

### What id a different number of clusters?
library(ClusterR)
set.seed(1991)
search_space <- c(5:12) ### Let's say that we want to check between 2 to 10 clusters, because there are 8 climatic types
system.time(opt_kmeans <- Optimal_Clusters_KMeans(data = PCsMapDF[,c( "Dim.1","Dim.2","Dim.3","Dim.4", "Dim.5")], 
                                                  max_clusters = search_space,
                                                  criterion = "WCSSE", num_init = 5,
                                                  max_iters = 5000, initializer = "kmeans++", plot_clusters = TRUE,
                                                  verbose = TRUE))
plot(search_space, opt_kmeans,pch=20, col="blue")
lines(search_space, opt_kmeans)

##### Here today

### Run silhouette better

# # perform KMeans_rcpp clustering
# library(ClusterR)
# my_seed <- 13 ### Your set.seed() number
# climatePCs_cl <-KMeans_rcpp(data=PCsMapDF[,c( "Dim.1","Dim.2","Dim.3","Dim.4","Dim.5")],
#                             clusters=9,
#                             num_init = 20,
#                             max_iters = 5000,
#                             initializer = "kmeans++",
#                             fuzzy = FALSE,
#                             verbose = FALSE,
#                             seed = my_seed)
# 
# ### save the cluster number in the original dataframe
# PCsMapDF$PCs9k <- as.factor(climatePCs_cl$clusters)
# ggmap(france_background) +
#   geom_point(aes(y = y, x = x, color=PCs8k),
#              cex=0.1, pch=19,
#              data=PCsMapDF) +
#   scale_color_viridis(discrete=TRUE,direction=-1, option="A")


################################################################################

### Mask oceans (in this case, use buffer around France)
### Subset bioclim variables + mask
BioclimSub <- terra::subset(bioclimrast, c(1,3,4,8,9,12,15,20))
names(BioclimSub) <- c("Bio1", "Bio3", "Bio4", "Bio8",
                       "Bio9", "Bio12", "Bio15", "mask")
target.vars <- c("Bio1", "Bio3", "Bio4", "Bio8",
                 "Bio9", "Bio12", "Bio15")

### Mask and replace in the same raster stack
for(i in 1:length(target.vars)){
  BioclimSub[[i]] <- terra::mask(BioclimSub[[i]], BioclimSub[["mask"]])
}

plot(BioclimSub)

### Calculate mean and sd of 8k BCE bioclimatic variables
rs1000BCE <- data.frame(mean=c(rep(NA,7)), sd=c(rep(NA,7)))

for(i in 1:length(target.vars)){
  # rs8kBCE[i,1] <- mean(values(BioclimSub[[i]]),na.rm=TRUE)
  # rs8kBCE[i,2] <- sd(values(BioclimSub[[i]]),na.rm=TRUE)
  ### Learnt that this function exists
  ### global(x, c("sum", "mean", "sd"), na.rm=TRUE)
  rs1000BCE[i,1] <- global(BioclimSub[[i]],"mean", na.rm=TRUE)
  rs1000BCE[i,2] <- global(BioclimSub[[i]],"sd",na.rm=TRUE)
}

### Scale 8k BCE bioclimatic variables
Bioclim.sc <- list()
for(i in 1:length(target.vars)){
  Bioclim.sc[[i]] <- terra::scale(BioclimSub[[i]], center=TRUE, scale=TRUE)
}
Bioclim.sc <- rast(Bioclim.sc)

### Take smaller sample for this exercise
set.seed(2233)
Bioclim.sc.Samp <- terra::spatSample(x = Bioclim.sc, 
                                     size= 40000,
                                     method="regular", 
                                     as.df=TRUE, 
                                     xy=TRUE)
### Only complete cases
Bioclim.sc.Samp <- Bioclim.sc.Samp[complete.cases(Bioclim.sc.Samp),]
dim(Bioclim.sc.Samp)

### Take larger sample for mapping and final clustering
set.seed(2233)
Bioclim.sc.Samp.l <- terra::spatSample(x = Bioclim.sc, 
                                     size= 700000,
                                     method="regular", 
                                     as.df=TRUE, 
                                     xy=TRUE)
### Only complete cases
Bioclim.sc.Samp.l <- Bioclim.sc.Samp.l[complete.cases(Bioclim.sc.Samp.l),]
dim(Bioclim.sc.Samp.l)

### Calculate the optimal number of k-means clusters with NbClust package
library(NbClust)
set.seed(1984)
system.time(opt.Clusters.dunn <- NbClust(data = Bioclim.sc.Samp[,3:ncol(Bioclim.sc.Samp)], 
                                         diss=NULL, 
                                         distance = "euclidean",
                                         min.nc = 2, 
                                         max.nc = 15,
                                         method = "kmeans",
                                         index = "dunn"))
gc()
opt.Clusters.dunn
summary(opt.Clusters.dunn)
opt.Clusters.dunn$Best.nc
plot(2:15, opt.Clusters.dunn$All.index,pch=20, col="blue")
lines(2:15, opt.Clusters.dunn$All.index)
abline(v=opt.Clusters.dunn$Best.nc, col="red", lty=2)
#### Optimum number is 15, after 6 and 7.

### the silhouette method is also very popular (check here)

library(NbClust)
set.seed(1984)
system.time(opt.Clusters.silhouette <- NbClust(data = Bioclim.sc.Samp[,3:ncol(Bioclim.sc.Samp)], 
                                               diss=NULL, 
                                               distance = "euclidean",
                                               min.nc = 2, 
                                               max.nc = 15,
                                               method = "kmeans", 
                                               index = "silhouette"))
opt.Clusters.silhouette
summary(opt.Clusters.silhouette)
opt.Clusters.silhouette$Best.nc
plot(2:15, opt.Clusters.silhouette$All.index,pch=20, col="blue")
lines(2:15, opt.Clusters.silhouette$All.index)
abline(v=opt.Clusters.silhouette$Best.nc, col="red", lty=2)
### 14 and then 6

# perform KMeans_rcpp clustering
my_seed <- 2245 ### Your set.seed() number
climate.1000BCE.rcpp <-KMeans_rcpp(data=Bioclim.sc.Samp[,3:9], 
                                clusters=14,
                                num_init = 10, 
                                max_iters = 5000,
                                initializer = "kmeans++", 
                                fuzzy = FALSE, 
                                verbose = FALSE,
                                seed = my_seed)

### Repeat with more samples
my_seed <- 2245 ### Your set.seed() number
climate.1000BCE.rcpp <-KMeans_rcpp(data=Bioclim.sc.Samp.l[,target.vars], 
                                   clusters=14,
                                   num_init = 10, 
                                   max_iters = 5000,
                                   initializer = "kmeans++", 
                                   fuzzy = FALSE, 
                                   verbose = FALSE,
                                   seed = my_seed)

### save the cluster number in the original dataframe
Bioclim.sc.Samp.l$K1000BCE <- climate.1000BCE.rcpp$clusters

### Plot clusters
ggplot() +
  geom_sf(data = france) +
  geom_point(aes(y = y, x = x, color=as.factor(K1000BCE)),
             cex=0.1, pch=19,
             data=Bioclim.sc.Samp.l) +
  scale_color_viridis(discrete=TRUE,direction=-1, option="A")

library(dendsort) 
library(dendextend) # visualize dendrograms
library(colorspace) # colors

### Hierarchical clustering
hc <- hclust(dist(climate.1000BCE.rcpp$centroids), method="ward.D2")
plot(dendsort(hc), main="Hierarchical clustering of 1000 BCE bioclim variables", sub="", xlab="")

### Extract labels
hc.labels <- hc %>% as.dendrogram(.) %>% labels %>% as.numeric()

### Extract the membership from the tree
dend <- hc %>% as.dendrogram(.)
plot(dend)

### Colours
mypalette <- c(sequential_hcl("TealGrn", n = 2),
               sequential_hcl("Blues", n = 2),
               sequential_hcl("PurpOr", n = 3),
               sequential_hcl("BurgYl", n = 2),
               "yellow",
               sequential_hcl("OrYel", n = 4))
legend.plot <- dend %>%  
  set("labels_cex", 2) %>% 
  set("branches_lwd", 2)%>% 
  set("labels_col", mypalette) %>% 
  set("branches_k_color", mypalette)
plot(legend.plot)

### Order for map
cluster_colors <- data.frame(orderID = hc.labels, colors = mypalette)
palette_order <- cluster_colors %>% arrange(cluster_colors, orderID)

### Plot with hierarchical clustering
ggplot() +
  geom_sf(data = france) +
  geom_point(aes(y = y, x = x, color=as.factor(K1000BCE)),
             cex=0.5, pch=19,
             data=Bioclim.sc.Samp.l) +
  scale_colour_manual(values=palette_order$colors)


### 3.2 Predict clusters from 10k BP in every time step ------------------

## Need to download and process Bio4

# Bring new layers using the function
setwd(climDir)
### List of time steps
timez <- as.character(c(-7:20))
### Extent to crop
myext <- ext(st_bbox(france_buffer_WGS84))

get_chelsa_paleo_france_crop <- function(climdir, timeIDs, vars, ExtD){
  #timeout_old <- getOption('timeout')
  #options(timeout=1000000)
  
  for(timeID in 1:length(timeIDs)){
    
    for(var in 1:length(vars)){
      
      ### Download file
      name <- paste0("CHELSA_TraCE21k_",vars[[var]],"_",timeIDs[[timeID]], "_V1.0.tif")
      source_url <- file.path(paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_trace/bio/",name))
      destination <- file.path(paste0(climdir,"TimeID_",timeIDs[[timeID]],"/",name))
      if(!file.exists(destination)){
        download.file(source_url, destination, method="wget")
      } else { 
        message(paste0(destination, " already downloaded, skipping to next"))
      }
    }
    #options(timeout = timeout_old)
    
    ### Process raster files for this time step
    setwd(paste0(climdir,"TimeID_",timeIDs[[timeID]],"/"))
    print(paste0("Cropping and storing rasters for France in ",paste0(climdir,"TimeID_",timeIDs[[timeID]],"/")))
    
    ### List tif files
    bioclim <- list.files(pattern=".tif")
    ### subset the variable in question
    bioclim <- bioclim[grep(pattern = vars, x = bioclim)]
    
    ### Crop to the extent of France (keep same alignment and resolution) 
    for(r in 1:length(bioclim)){
      
      ### Load raster
      bioclimrast <- terra::rast(paste0(climdir,"TimeID_",timeIDs[[timeID]],"/",bioclim[[r]]))
      ### Crop to the extent of France
      bioclimrastc <- terra::crop(bioclimrast, ExtD)
      ### Write Raster
      writeRaster(bioclimrastc, filename = paste0(names(bioclimrastc),"_Fr.tif"))
      ### Remove raster global extent
      file.remove(paste0(climdir,"TimeID_",timeIDs[[timeID]],"/",bioclim[[r]]))
      
    } 
  }
}

save.image("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/Rsessions/05122023.RData")

### continue here
### Download Bio4 for all the time steps
get_chelsa_paleo_france_crop(climdir = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Climate/",
                             timeIDs = timez, 
                             vars = "bio04",
                             ExtD = myext)

### function to extract coordinates, scale, 
### and predict cluster from 8k BCE in new data
### @climdir the directory whwre the paleoclimate variables for that time step are stored
### @timeID time step of interest
### @vars Bioclim variables of interest, name of columns in the dataframe
### @dfScale dataframe with mean and sd of the Bioclimatic variables from the first time step
### @SampleDF the dataframe with the sample that we use for clustering 
### the first time step (1000 BCE in this case)
### @kmeansM the kmeans model (package ClusterR) that we use for predicting

### The function returns the predictions in the new time step. 
### We can attach to the dataframe SampleDF as output

### Function to assing cluster
assign_1000BCE_cluster_timeSteps <- function(climdir,
                                             timeID, 
                                             vars,
                                             vars.names,
                                             dfScale,
                                             SampleDF,
                                             kmeansM){
  
  ### Process raster files for this time step
  setwd(paste0(climdir,"TimeID_",timeID,"/"))
  
  print(paste0("Extracting bioclim variables for ",paste0("TimeID = ",timeID)))
  
  ### List tif files
  #bioclim <- list.files(pattern="Fr.tif")
  ### subset the variable in question
  bioclim <- paste0("CHELSA_TraCE21k_",vars.names,"_",timeID,"_V1.0_Fr.tif")
  
  ### Load raster stack
  bioclimrast <- terra::rast(paste0(climdir,"TimeID_",timeID,"/",bioclim))
  
  ### Scale with the mean and sd from 8k BCE
  Bioclim.sc <- terra::scale(bioclimrast, 
                             center = dfScale[,"mean"],
                             scale = dfScale[,"sd"])
  
  ### Take regular sample in same coordinates as for 8000 BCE
  extractCoords <- as.matrix(SampleDF[,c("x", "y")])
  extraction <- terra::extract(Bioclim.sc, extractCoords)
  
  ### Rename columns with vars
  colnames(extraction) <- vars
  
  ## Create empty prediction column
  SampleDF$cluster <- NA
  
  ### Predict assignment in new raster
  SampleDF$cluster <- predict_KMeans(data = extraction, CENTROIDS = kmeansM$centroids)
  
  ### Rename cluster for that time step
  colnames(SampleDF)[colnames(SampleDF)=="cluster"] <- paste0("cl",timeID)
  
  ### Return dataframe, extraction df, and
  return(SampleDF)
  
} 


### Apply function in for loop
timez <- as.character(c(-8:20))

for(i in 1:length(timez)) {
  
  if(i == 1) {
    print(timez[[i]])
    ### Apply function for each time step
    out <- assign_1000BCE_cluster_timeSteps(climdir = climDir, 
                                            timeID = timez[[i]], 
                                            vars = c("Bio1", "Bio3", "Bio4", "Bio8", "Bio9", "Bio12", "Bio15"),
                                            vars.names = c("bio01", "bio03", "bio04", "bio08", "bio09", "bio12", "bio15"),
                                            dfScale = rs1000BCE,
                                            SampleDF = Bioclim.sc.Samp.l,
                                            kmeansM = climate.1000BCE.rcpp)
    
  } else if (i > 1) {
    
    print(timez[[i]])
    
    ### Apply function for each time step
    out <- assign_1000BCE_cluster_timeSteps(climdir = climDir, 
                                          timeID = timez[[i]], 
                                          vars = c("Bio1", "Bio3", "Bio4", "Bio8", "Bio9", "Bio12", "Bio15"),
                                          vars.names = c("bio01", "bio03", "bio04", "bio08", "bio09", "bio12", "bio15"),
                                          dfScale = rs1000BCE,
                                          SampleDF = out,
                                          kmeansM = climate.1000BCE.rcpp)
  }
}


### Plot clusters
ggplot() +
  geom_sf(data = france) +
  geom_point(aes(y = y, x = x, color=as.factor(`cl-8`)),
             cex=1, pch=19,
             data=out) +
  scale_colour_manual(values=palette_order$colors)

### Plot evolution over several time steps (not all of them)

### change column names
colnames(out) <- c("x","y", target.vars, paste0("timeID_",-9:20))

p1 <- ggplot() +
  geom_sf(data = france) +
  geom_point(aes(y = y, x = x, color=as.factor(`timeID_-8`)),
             cex=1, pch=19,
             data=out) +
  scale_colour_manual(values=palette_order$colors)

p2 <- ggplot() +
  geom_sf(data = france) +
  geom_point(aes(y = y, x = x, color=as.factor(`timeID_0`)),
             cex=1, pch=19,
             data=out) +
  scale_colour_manual(values=palette_order$colors)

p3 <- ggplot() +
  geom_sf(data = france) +
  geom_point(aes(y = y, x = x, color=as.factor(`timeID_7`)),
             cex=1, pch=19,
             data=out) +
  scale_colour_manual(values=palette_order$colors)

p4 <- ggplot() +
  geom_sf(data = france) +
  geom_point(aes(y = y, x = x, color=as.factor(`timeID_20`)),
             cex=1, pch=19,
             data=out) +
  scale_colour_manual(values=palette_order$colors)

# library(gridExtra)
# grid.arrange(p1,p2,p3, ncol=3)

### Pivor to long
cluster.long <- out %>%
  pivot_longer(., cols=10:39,
               names_to = "timeID",
               values_to = "cluster"  )

cluster.summary <- cluster.long %>%
  group_by(., timeID, cluster) %>%
  summarise(., N=n()) 

cluster.summary$timeStep <- as.numeric(gsub(x=cluster.summary$timeID, 
                                            pattern= "timeID_",
                                            replacement = ""))

ggplot()+
  geom_point(aes(x=timeStep, y=N, color=as.factor(cluster)), 
             data=cluster.summary) +
  geom_smooth(aes(x=timeStep, y=N, color=as.factor(cluster)), 
              data=cluster.summary, se=FALSE) +
  scale_color_manual(values=palette_order$colors)


############################################################################