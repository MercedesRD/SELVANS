#############################################################################################################################################
###  Selection of scorpan variables for generating pedogenon classes for the Basque Country
###
###  In this script: Covariate selection based on predictive ability for several soil properties

### Desired extent: Basque Country
### Resolution: 25m
### CRS: EPSG=25830
###  Author: Mercedes Roman Dobarco
###  Date: 09/02/2024

####### Load packages
### Spatial
library(sf)
library(terra)

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

### Machine learning
library(Cubist)
require(ranger)

# 1 Load all scorpan variables --------------------------------------

HomeDir <- "C:/Covariates/Euskadi/"
setwd(HomeDir)
eus_buffer <- rast("eus_buffer_125m.tif")

### Soil
setwd("C:/Covariates/Euskadi/Soil/SoilGrids/")
soil.files <- list.files(pattern=".tif$")
soil.rast <- rast(c("cec_30-60cm_mean_eus.tif","cec_60-100cm_mean_eus.tif","cec_100-200cm_mean_eus.tif",  
                    "clay_30-60cm_mean_eus.tif","clay_60-100cm_mean_eus.tif","clay_100-200cm_mean_eus.tif", 
                    "sand_30-60cm_mean_eus.tif","sand_60-100cm_mean_eus.tif","sand_100-200cm_mean_eus.tif",
                    "silt_30-60cm_mean_eus.tif","silt_60-100cm_mean_eus.tif","silt_100-200cm_mean_eus.tif"))
names(soil.rast) <- c("cec_30_60","cec_60_100","cec_100_200",
                      "clay_30_60","clay_60_100","clay_100_200",
                      "sand_30_60","sand_60_100","sand_100_200",
                      "silt_30_60","silt_60_100","silt_100_200")
# raster.names <- c("cec_30_60","cec_60_100","cec_100_200",
#                  "clay_30_60","clay_60_100","clay_100_200",
#                  "sand_30_60","sand_60_100","sand_100_200",
#                  "silt_30_60","silt_60_100","silt_100_200")
# plot(soil.rast)
# 
# for(i in 1:nlyr(soil.rast)){
#   print(soil.rast[[i]])
#   m <- terra::mask(x=soil.rast[[i]], mask=eus_buffer) # Mask pixels outside Euskadi
#   s <- terra::scale(x=m,center=TRUE, scale=TRUE) # Scale because it is a continuous variable
#   setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/4.CovariateSelection/scaled/") # change wd
#   writeRaster(s, filename = paste0(raster.names[[i]],".tif"),overwrite = TRUE)
#   # format = "GTiff",na.rm=T, inf.rm=T, ) # Write to file
#   #s # Return scaled raster
# }

### Climate
setwd(paste0(HomeDir,"Climate/"))
climate.files <- sort(list.files(pattern=".tif"))
clim.rast <- rast(climate.files)

### Aridity index
clim_ari <- clim.rast$`CHELSA_pet_penman_mean_1981-2010_V.2.1_eus`/ clim.rast$`CHELSA_bio12_1981-2010_V.2.1_eus`
names(clim_ari) <- "clim_ari"
writeRaster(clim_ari, "clim_ari.tif")

climate.files <- sort(list.files(pattern=".tif"))
clim.rast <- rast(climate.files)

names(clim.rast) <- c("clim_bio1","clim_bio10","clim_bio11",
                      "clim_bio12","clim_bio13","clim_bio14",
                      "clim_bio15","clim_bio16","clim_bio17",
                      "clim_bio18","clim_bio19","clim_bio2",
                      "clim_bio3","clim_bio4","clim_bio5",
                      "clim_bio6","clim_bio7","clim_bio8",
                      "clim_bio9","clim_pet","clim_ari") #Mean monthly potential evapotranspiration

### Relief
setwd("C:/Covariates/Euskadi/Relief/")

### Covariates subset
relief.files <- c("mdt_lidar_2017_25m_etrs89.tif", # DEM
                          "slope.tif", "slope_5.tif", "slope_10.tif", # Slope calculated from 3,5,10 cells
                          "easterness.tif", # Easterness 
                          "northerness.tif", # northerness
                          "nor_slope.tif" , # northness x slope
                          "twi_saga.tif", # SAGA Topographic wetness index
                          "p_curv.tif", "pr_curv_5.tif", "pr_curv_10.tif", # profile curvature (3,5,10)
                          "pl_curv_5.tif", "pl_curv_10.tif", # planar curvature (5,10)
                          "lg_curv_5.tif", "lg_curv_10.tif", # longitudinal curvature (5,10)
                          "t_curv.tif", "cs_curv_5.tif", "cs_curv_10.tif" , # tangential curvature (3) cross-sectional curvature (5,10) (same)
                          "mrvbf.tif", "mrrtf.tif", ### MRVBF, MRRTF
                          "slope_height.tif", # Slope Height
                          "norm_height.tif", # Normalized height
                          "st_height.tif", # Standardized Height
                          "valley_depth.tif", # Valley depth
                          "mid_slope.tif",  # Mid slope
                          "tpi_8_3.tif", "tpi_20_5.tif")  # Multiscale topographic position index, calculated with a search radius of 8,20 cells, and 3,5 scales respectively
                         # "geomorphon.tif", "geomorphon2.tif")  ## Geomorphons (search radius of 11, skip 1 & 3, flat 1 & 1.5)

relief.rast <- terra::rast(relief.files)
plot(relief.rast)
plot(relief.rast[[c(16:19)]]) 
plot(relief.rast[[c(16:17)]])  
names(relief.rast)[names(relief.rast)=="mdt_lidar_2017_25m_etrs89"] <- "dem"
names(relief.rast)[names(relief.rast)=="twi_saga"] <- "twi"
names(relief.rast)[names(relief.rast)=="slope@relief"] <- "slope"
names(relief.rast)[names(relief.rast)=="elevation_easterness@relief"] <- "easterness"
names(relief.rast)[names(relief.rast)=="elevation_northerness@relief"] <- "northerness"
names(relief.rast)[names(relief.rast)=="elevation_northerness_slope@relief"] <- "north_slope"
names(relief.rast)[names(relief.rast)=="t_curv@relief"] <- "t_curv"
names(relief.rast)[names(relief.rast)=="p_curv@relief"] <- "p_curv"
names(relief.rast)[names(relief.rast)=="mrvbf@relief"] <- "mrvbf"
names(relief.rast)[names(relief.rast)=="mrrtf@relief"] <- "mrrtf"

### Parent material and organisms from MCA dimensions
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/1.CovariatesEus/")
mca.files <- paste0("MCA_",1:32,".tif")
mca.rast <- rast(mca.files)
names(mca.rast) <- paste0("MCA_",1:32)

### Merge all
covariates.rast <- c(soil.rast,clim.rast,relief.rast,mca.rast)
raster.names <- names(covariates.rast)

#### scale all variables
for(i in 1:nlyr(covariates.rast)){
  print(covariates.rast[[i]])
  m <- terra::mask(x=covariates.rast[[i]], mask=eus_buffer) # Mask pixels outside Euskadi
  s <- terra::scale(x=m,center=TRUE, scale=TRUE) # Scale because it is a continuous variable
  setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/4.CovariateSelection/scaled/") # change wd
  writeRaster(s, filename = paste0(raster.names[[i]],".tif"),overwrite = TRUE)
  # format = "GTiff",na.rm=T, inf.rm=T, ) # Write to file
  #s # Return scaled raster
}

### Note: I manually added "relief_" before all relief variables

### load scaled rasters
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/4.CovariateSelection/scaled/")
scaled_files <- sort(list.files(pattern=".tif"))
scaled_rast <- rast(scaled_files)
names(scaled_rast)[names(scaled_rast)=="clim_bio21"] <- "clim_bio2"

# ### 2. Extract covariates at soil observations --------------------------

load("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/SoilDatasets/Hamonised_Dataset.RData")
### clean space

### I keep Soil.df.harmonised and Soil.df.harmonised_sf - splined
### soil_datasets_3 (before splines and with correct coordinates for BASONET)

### Extract
soil.covariates <- terra::extract(x=scaled_rast, y=Soil.df.harmonised_sf, method="simple")

### bind
soil.scorpan <- cbind(Soil.df.harmonised, soil.covariates)

### How many observations per soil property?

soil.scorpan.summary.dataset <- soil.scorpan %>% 
  dplyr::select(c("newID","Dataset","Layer_depth",
              "Sand","Silt","Clay","CEC","CECef","CoarseFrag",
              "CoarseFrag_Grav","Bulk_Density","pH","EC","TOC","Carbonates",
              "TN_ppm","P_ppm","K_ppm","Ca_ppm","Mg_ppm","Na_ppm")) %>%
  pivot_longer(cols =c("Sand","Silt","Clay","CEC","CECef","CoarseFrag",
                       "CoarseFrag_Grav","Bulk_Density","pH","EC","TOC","Carbonates",
                       "TN_ppm","P_ppm","K_ppm","Ca_ppm","Mg_ppm","Na_ppm"),
               names_to = "soil_prop",values_to ="value")  %>%
  filter(., !is.na(value))%>%
  count(Dataset,Layer_depth, soil_prop, sort = TRUE) %>%
  pivot_wider(names_from = soil_prop, values_from=n) %>% as.data.frame()
  

soil.scorpan.summary <- soil.scorpan %>% 
  dplyr::select(c("newID","Dataset","Layer_depth",
                  "Sand","Silt","Clay","CEC","CECef","CoarseFrag",
                  "CoarseFrag_Grav","Bulk_Density","pH","EC","TOC","Carbonates",
                  "TN_ppm","P_ppm","K_ppm","Ca_ppm","Mg_ppm","Na_ppm")) %>%
  pivot_longer(cols =c("Sand","Silt","Clay","CEC","CECef","CoarseFrag",
                       "CoarseFrag_Grav","Bulk_Density","pH","EC","TOC","Carbonates",
                       "TN_ppm","P_ppm","K_ppm","Ca_ppm","Mg_ppm","Na_ppm"),
               names_to = "soil_prop",values_to ="value")  %>%
  filter(., !is.na(value))%>%
  count(Layer_depth, soil_prop, sort = TRUE) %>%
  pivot_wider(names_from = soil_prop, values_from=n) %>% as.data.frame()

soil.scorpan.summary <- soil.scorpan.summary[,c("Layer_depth",names(sort(colSums(soil.scorpan.summary[,2:ncol(soil.scorpan.summary)],na.rm=TRUE), decreasing = TRUE)))]
soil.scorpan.summary <- soil.scorpan.summary %>% bind_rows(summarise(.,
                                                                     across(where(is.numeric), sum, na.rm=TRUE),
                                                                     across(where(is.character), ~"Total")))
soil.scorpan.summary <- soil.scorpan.summary[1:6,]

# 3. Fit randomForest models ----------------------------------------------

### Change colnames - problem with "-" in soil variables
# names(soil.scorpan)[names(soil.scorpan)=="cec_100-200cm"] <- "cec_100_200"
# names(soil.scorpan)[names(soil.scorpan)=="cec_30-60cm"] <- "cec_30_60"
# names(soil.scorpan)[names(soil.scorpan)=="cec_60-100cm"] <- "cec_60_100"
# names(soil.scorpan)[names(soil.scorpan)=="clay_100-200cm"] <- "clay_100_200"
# names(soil.scorpan)[names(soil.scorpan)=="clay_30-60cm"] <- "clay_30_60"
# names(soil.scorpan)[names(soil.scorpan)=="clay_60-100cm"] <- "clay_60_100"
# names(soil.scorpan)[names(soil.scorpan)=="sand_100-200cm"] <- "sand_100_200"
# names(soil.scorpan)[names(soil.scorpan)=="sand_30-60cm"] <- "sand_30_60"
# names(soil.scorpan)[names(soil.scorpan)=="sand_60-100cm_"] <- "sand_60_100"
# names(soil.scorpan)[names(soil.scorpan)=="silt_100-200cm"] <- "silt_100_200"
# names(soil.scorpan)[names(soil.scorpan)=="silt_30-60cm"] <- "silt_30_60"
# names(soil.scorpan)[names(soil.scorpan)=="silt_60-100cm"] <- "silt_60_100"

## Create list of variables by soil forming factor

### soil vars
soil_vars <- c("cec_30_60","cec_60_100","cec_100_200",
               "clay_30_60","clay_60_100","clay_100_200",
               "sand_30_60","sand_60_100","sand_100_200",
               "silt_30_60","silt_60_100","silt_100_200")

### climate vars
clim_vars <- c("clim_bio1","clim_bio2",
               "clim_bio3","clim_bio4","clim_bio5",
               "clim_bio6","clim_bio7","clim_bio8",
               "clim_bio9","clim_bio10","clim_bio11",
              "clim_bio12","clim_bio13","clim_bio14",
              "clim_bio15","clim_bio16","clim_bio17",
              "clim_bio18","clim_bio19","clim_pet","clim_ari")

### relief vars
relief_vars <- c("dem", # DEM
                  "slope", "slope_5", "slope_10", # Slope calculated from 3,5,10 cells
                  "easterness", # Easterness 
                  "northerness", # northerness
                  "north_slope" , # northness x slope
                  "twi", # SAGA Topographic wetness index
                  "p_curv", "pr_curv_5", "pr_curv_10", # profile curvature (3,5,10)
                  "pl_curv_5", "pl_curv_10", # planar curvature (5,10)
                  "lg_curv_5", "lg_curv_10", # longitudinal curvature (5,10)
                  "t_curv", "cs_curv_5", "cs_curv_10" , # tangential curvature (3) cross-sectional curvature (5,10) (same)
                  "mrvbf", "mrrtf", ### MRVBF, MRRTF
                  "slope_height", # Slope Height
                  "norm_height", # Normalized height
                  "st_height", # Standardized Height
                  "valley_depth", # Valley depth
                  "mid_slope",  # Mid slope
                  "tpi_8_3", "tpi_20_5") 

### mca.dim vars
mca_vars <- paste0("MCA_",1:32)

### Write for foreach loop to fit models for all depth intervals and variable subset

### Depths
depths <- c("000_010_cm","010_020_cm","020_040_cm","040_060_cm","060_100_cm")

### Variables
target_soil <- c("Sand","Silt","Clay","CEC","CECef","CoarseFrag",
                 "CoarseFrag_Grav","Bulk_Density","pH","EC","TOC","Carbonates",
                 "TN_ppm","P_ppm","K_ppm","Ca_ppm","Mg_ppm","Na_ppm")

### All scorpan variables
scorpan_vars <- c(soil_vars,clim_vars,relief_vars,mca_vars) ## 92 variables

### Short function to plot variable importance
### Plot varImp
myplotvarImp <- function(ranger.model, var.target) { 
  require(ggplot2)
  ### Dataframe with variable importance
  DF <- data.frame(covariates=names(ranger.model$variable.importance),
                   Importance=as.vector(ranger.model$variable.importance))
  ### Arrange
  DF <- DF %>% dplyr::arrange(., Importance) %>% as.data.frame()
  DF$covariates <- factor(DF$covariates, levels = DF$covariates)
  plotVImp <- ggplot(DF, aes(x=Importance, y= covariates, fill=Importance))+ 
    geom_bar(stat="identity", position="dodge") + 
    xlab("Covariate Importance") +
    ylab("") +
    ggtitle(var.target)+
    guides(fill=F)+
    scale_fill_gradient(low="gray80", high="blue")
  
  ### Arrange again, with most important variables first
  DF <- DF %>% dplyr::arrange(., desc(Importance)) %>% as.data.frame()
  return(list(DF,plotVImp))
  }

### Fit random model by variable and depth interval for clay, silt, sand, CEC, coarse fragments, carbonates, CECef

detectCores()
cl <- makeCluster(20)   ### Create cluster
registerDoParallel(cl)
getDoParWorkers()

tic <- Sys.time()

# out.rfmodels <- foreach(i=1:length(depths),.packages=c("ranger", "Boruta", "ggplot2"),
#                         .export = c("target_soil","soil.scorpan", "myplotvarImp")) %:% 
  
out.rfmodels <- foreach (j=1:length(target_soil), .packages=c("ranger","ggplot2","dplyr","tidyverse"), 
                         .export = c("target_soil","soil.scorpan","myplotvarImp",
                                     "soil_vars","clim_vars","relief_vars","mca_vars")) %dopar% {
                                       
 ### All scorpan variables
 scorpan_vars <- c(soil_vars,clim_vars,relief_vars,mca_vars)
                           
 ### Transform depth to factor
 soil.scorpan$Layer_depth <- factor(soil.scorpan$Layer_depth, levels=depths)
    
 ### subset depth
 # scorpan.df <- soil.scorpan[soil.scorpan$Layer_depth == depths[[i]],]
 # # make the formula for target variable and all scorpan
 # form1 <- as.formula(paste0(target_soil[[j]],"~",paste(scorpan_vars, collapse="+")))
 # scorpan.df <- scorpan.df[,c(target_soil[[j]],scorpan_vars)]
 # scorpan.df <- scorpan.df[complete.cases(scorpan.df),]
    
 ### without subsetting depth
 ### all variables
 # make the formula for target variable and all scorpan
 form1 <- as.formula(paste0(target_soil[[j]],"~ Layer_depth +",paste(scorpan_vars, collapse="+")))
 ### subset dataframe with complete observations for target variable
 scorpan.df <- soil.scorpan[,c(target_soil[[j]],"Layer_depth",scorpan_vars)]
 scorpan.df <- scorpan.df[complete.cases(scorpan.df),]
 
 set.seed(3899)
 rf.j.impurity <- ranger(form1, 
                 data = scorpan.df,
                 num.trees = 5000,
                 importance = 'impurity_corrected')
 
 df.varImp.impurity <- myplotvarImp(ranger.model = rf.j.impurity, var.target = target_soil[[j]])
 
 set.seed(3899)
 rf.j.permutation <- ranger(form1, 
                    data = scorpan.df,
                    num.trees = 5000,
                    importance = 'permutation',
                    scale.permutation.importance = TRUE)
    
 df.varImp.permutation <- myplotvarImp(ranger.model = rf.j.permutation, var.target = target_soil[[j]])

 ### climate variables
 # make the formula for target variable and climate variables
 form1 <- as.formula(paste0(target_soil[[j]],"~ Layer_depth +",paste(clim_vars, collapse="+")))
 ### subset dataframe with complete observations
 climate.df <- soil.scorpan[,c(target_soil[[j]],"Layer_depth",clim_vars)]
 climate.df <- climate.df[complete.cases(climate.df),]
 
 set.seed(3899)
 rf.j.impurity <- ranger(form1, 
                         data = climate.df,
                         num.trees = 5000,
                         importance = 'impurity_corrected')
 
 df.varImp.impurity.cl <- myplotvarImp(ranger.model = rf.j.impurity, var.target = target_soil[[j]])
 
 set.seed(3899)
 rf.j.permutation <- ranger(form1, 
                            data = climate.df,
                            num.trees = 5000,
                            importance = 'permutation',
                            scale.permutation.importance = TRUE)
 
 df.varImp.permutation.cl <- myplotvarImp(ranger.model = rf.j.permutation, var.target = target_soil[[j]])
 
 ### relief variables
 # make the formula for target variable and climate variables
 form1 <- as.formula(paste0(target_soil[[j]],"~ Layer_depth +",paste(relief_vars, collapse="+")))
 ### subset dataframe with complete observations
 relief.df <- soil.scorpan[,c(target_soil[[j]],"Layer_depth",relief_vars)]
 relief.df <- relief.df[complete.cases(relief.df),]
 
 set.seed(3899)
 rf.j.impurity <- ranger(form1, 
                         data = relief.df,
                         num.trees = 5000,
                         importance = 'impurity_corrected')
 
 df.varImp.impurity.relf <- myplotvarImp(ranger.model = rf.j.impurity, var.target = target_soil[[j]])
 
 set.seed(3899)
 rf.j.permutation <- ranger(form1, 
                            data = relief.df,
                            num.trees = 5000,
                            importance = 'permutation',
                            scale.permutation.importance = TRUE)
 
 df.varImp.permutation.relf <- myplotvarImp(ranger.model = rf.j.permutation, var.target = target_soil[[j]])
 
 ### MCA variables
 # make the formula for target variable and climate variables
 form1 <- as.formula(paste0(target_soil[[j]],"~ Layer_depth +",paste(mca_vars, collapse="+")))
 ### subset dataframe with complete observations
 mca.df <- soil.scorpan[,c(target_soil[[j]],"Layer_depth",mca_vars)]
 mca.df <- mca.df[complete.cases(mca.df),]
 
 set.seed(3899)
 rf.j.impurity <- ranger(form1, 
                         data = mca.df,
                         num.trees = 5000,
                         importance = 'impurity_corrected')
 
 df.varImp.impurity.mca <- myplotvarImp(ranger.model = rf.j.impurity, var.target = target_soil[[j]])
 
 set.seed(3899)
 rf.j.permutation <- ranger(form1, 
                            data = mca.df,
                            num.trees = 5000,
                            importance = 'permutation',
                            scale.permutation.importance = TRUE)
 
 df.varImp.permutation.mca <- myplotvarImp(ranger.model = rf.j.permutation, var.target = target_soil[[j]])
 
  
 ### Soil variables
 # make the formula for target variable and climate variables
 form1 <- as.formula(paste0(target_soil[[j]],"~ Layer_depth +",paste(soil_vars, collapse="+")))
 ### subset dataframe with complete observations
 soil.df <- soil.scorpan[,c(target_soil[[j]],"Layer_depth",soil_vars)]
 soil.df <- soil.df[complete.cases(soil.df),]
 
 set.seed(3899)
 rf.j.impurity <- ranger(form1, 
                         data = soil.df,
                         num.trees = 5000,
                         importance = 'impurity_corrected')
 
 df.varImp.impurity.soil <- myplotvarImp(ranger.model = rf.j.impurity, var.target = target_soil[[j]])
 
 set.seed(3899)
 rf.j.permutation <- ranger(form1, 
                            data = soil.df,
                            num.trees = 5000,
                            importance = 'permutation',
                            scale.permutation.importance = TRUE)
 
 df.varImp.permutation.soil <- myplotvarImp(ranger.model = rf.j.permutation, var.target = target_soil[[j]])
 
 gc()
 ### Save in a list
 all.dfs.imp <- list(df.varImp.impurity, df.varImp.permutation,
                     df.varImp.impurity.cl,df.varImp.permutation.cl,
                     df.varImp.impurity.relf,df.varImp.permutation.relf,
                     df.varImp.impurity.mca,df.varImp.permutation.mca,
                     df.varImp.impurity.soil,df.varImp.permutation.soil)
 
 all.dfs.imp # We return this
 
 }
gc()

tac <- Sys.time()
stopCluster(cl)
tac-tic  
# Time difference of 10.70619 mins

### Rank variables - top 10 or top 5 by soil-forming factor
# climate.dfs.permutation <- data.frame(covariates=c(clim_vars,"Layer_depth"))
# #climate.dfs.permutation$covariates <- as.factor(climate.dfs.permutation$covariates)
# relief.dfs.permutation <- data.frame(covariates=c(relief_vars,"Layer_depth"))
# #relief.dfs.permutation$covariates <- as.factor(relief.dfs.permutation$covariates)
# mca.dfs.permutation <- data.frame(covariates=c(mca_vars,"Layer_depth"))
# #mca.dfs.permutation$covariates <- as.factor(mca.dfs.permutation$covariates)
# soil.dfs.permutation <- data.frame(covariates=c(soil_vars,"Layer_depth"))
# #soil.dfs.permutation$covariates <- as.factor(soil.dfs.permutation$covariates)

### Merge all dataframes into a single one, by soil-forming factor
for (i in 1:length(target_soil)) {
  
  if(i==1) {
    
    climate.dfs.permutation <- out.rfmodels[[i]][[4]][[1]]
    climate.dfs.permutation$covariates <- as.character(climate.dfs.permutation$covariates )
    colnames(climate.dfs.permutation) <- c("covariates",target_soil[[i]])
    
    relief.dfs.permutation <- out.rfmodels[[i]][[6]][[1]]
    relief.dfs.permutation$covariates <- as.character(relief.dfs.permutation$covariates )
    colnames(relief.dfs.permutation) <- c("covariates",target_soil[[i]])
    
    mca.dfs.permutation <- out.rfmodels[[i]][[8]][[1]]
    mca.dfs.permutation$covariates <- as.character(mca.dfs.permutation$covariates)
    colnames(mca.dfs.permutation) <- c("covariates",target_soil[[i]])
    
    soil.dfs.permutation <- out.rfmodels[[i]][[10]][[1]]
    soil.dfs.permutation$covariates <- as.character(soil.dfs.permutation$covariates)
    colnames(soil.dfs.permutation) <- c("covariates",target_soil[[i]])
    
  } else if (i > 1) {
    
    out.rfmodels[[i]][[4]][[1]]$covariates <- as.character(out.rfmodels[[i]][[4]][[1]]$covariates)
    colnames(out.rfmodels[[i]][[4]][[1]]) <- c("covariates",target_soil[[i]])
    climate.dfs.permutation <- dplyr::left_join(climate.dfs.permutation, out.rfmodels[[i]][[4]][[1]], by="covariates")
    
    out.rfmodels[[i]][[6]][[1]]$covariates <- as.character(out.rfmodels[[i]][[6]][[1]]$covariates)
    colnames(out.rfmodels[[i]][[6]][[1]]) <- c("covariates",target_soil[[i]])
    relief.dfs.permutation <- dplyr::left_join(relief.dfs.permutation,out.rfmodels[[i]][[6]][[1]], by="covariates")
    
    out.rfmodels[[i]][[8]][[1]]$covariates <- as.character(out.rfmodels[[i]][[8]][[1]]$covariates)
    colnames(out.rfmodels[[i]][[8]][[1]]) <- c("covariates",target_soil[[i]])
    mca.dfs.permutation <- dplyr::left_join(mca.dfs.permutation,out.rfmodels[[i]][[8]][[1]], by="covariates")
    
    out.rfmodels[[i]][[10]][[1]]$covariates <- as.character(out.rfmodels[[i]][[10]][[1]]$covariates)
    colnames(out.rfmodels[[i]][[10]][[1]]) <- c("covariates",target_soil[[i]])
    soil.dfs.permutation <- dplyr::left_join(soil.dfs.permutation,out.rfmodels[[i]][[10]][[1]], by="covariates")

  }
  
  }

### Are numbers comparable across soil-forming factors?


# 4. Climate subset -----------------------------------------------------------------

setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/4.CovariateSelection")
#### Create Heatmap figure
### Reorder covariates by soil-forming factor
library(tidyverse)

### Climate
clim_vars <- c("clim_bio1","clim_bio2",
               "clim_bio3","clim_bio4","clim_bio5",
               "clim_bio6","clim_bio7","clim_bio8",
               "clim_bio9","clim_bio10","clim_bio11",
               "clim_bio12","clim_bio13","clim_bio14",
               "clim_bio15","clim_bio16","clim_bio17",
               "clim_bio18","clim_bio19","clim_pet","clim_ari")

rownames(climate.dfs.permutation) <- climate.dfs.permutation$covariates
### eliminate soil depth
climate.dfs.permutation <- climate.dfs.permutation[!rownames(climate.dfs.permutation)=="Layer_depth",]
### Arrange by variable name
climate.dfs.permutation <- climate.dfs.permutation %>% arrange(factor(covariates, levels = c(clim_vars)))
varImp <-climate.dfs.permutation

### Order the rows based on stable/dynamic proerties and number of observations
colnames(soil.scorpan.summary)

soil.order <- c("Clay","Sand","Silt","Carbonates","CECef","CoarseFrag","CEC","CoarseFrag_Grav",
  "TOC","pH","Ca_ppm","K_ppm","TN_ppm","Mg_ppm","EC","P_ppm","Bulk_Density","Na_ppm")          
#varImp <- dplyr::left_join(data.frame(covariates=c(clim_vars,"Layer_depth")),climate.dfs.permutation, by="covariates")

#library(ggplot2)
library(gplots)
library(tidyverse)
heatmap(t(as.matrix(varImp[,soil.order])),
        cexRow = 0.6, 
        cexCol = 0.8,
        margins = c(10, 6))

mypal <-viridis_pal(option = "A",direction=-1)(20)

heatmap.2(t(as.matrix(varImp[,soil.order])),
          dendrogram = "none",
          Colv = FALSE,Rowv = FALSE,
          cexRow = 0.8, notecol="transparent", tracecol = "transparent", 
          cexCol = 0.8,scale="none",col=mypal,
          margins = c(10, 6),
          key.title="Variable importance")

heatmap.2(t(as.matrix(varImp[,c("Clay","Sand","Silt", "TOC","pH",
                               "Ca_ppm","K_ppm","Mg_ppm",
                               "Carbonates","CECef")])), 
          dendrogram = "none",
          Colv = FALSE, Rowv = FALSE,
          cexRow = 0.8, notecol="transparent", tracecol = "transparent", 
          cexCol = 0.8, scale="none",col=mypal,
          margins = c(10, 6),
          key.title="Variable importance")

### Sort all variables by importance
### Can I extract ranks across serveral dataframes?
### Extract top 10 variables by soil-forming factor and for clay, silt, sand, CEC,...
library(RobustRankAggreg)
### Rank variables - top 10 or top 5 by soil-forming factor
climate.list.perm <- list()

target_soil
### Merge all dataframes into a single one, by soil-forming factor
### But only of a subset of variables , and top 10 from each
climate.list.perm <- list(out.rfmodels[[1]][[4]][[1]][1:10,1], # sand
                          out.rfmodels[[2]][[4]][[1]][1:10,1], # silt
                          out.rfmodels[[3]][[4]][[1]][1:10,1], # clay
                          out.rfmodels[[9]][[4]][[1]][1:10,1], # pH
                          out.rfmodels[[11]][[4]][[1]][1:10,1], #TOC
                          out.rfmodels[[16]][[4]][[1]][1:10,1], #Ca
                          out.rfmodels[[15]][[4]][[1]][1:10,1], # K
                          out.rfmodels[[17]][[4]][[1]][1:10,1], # Mg
                          out.rfmodels[[12]][[4]][[1]][1:10,1], # Carbonates
                          out.rfmodels[[5]][[4]][[1]][1:10,1]) # CECef
aggregateRanks(glist = climate.list.perm,method='stuart')

### But only of a subset of variables, all of them
climate.list.perm <- list(out.rfmodels[[1]][[4]][[1]][,1], # sand
                          out.rfmodels[[2]][[4]][[1]][,1], # silt
                          out.rfmodels[[3]][[4]][[1]][,1], # clay
                          #out.rfmodels[[9]][[4]][[1]][,1], # pH
                          #out.rfmodels[[11]][[4]][[1]][,1], #TOC
                          out.rfmodels[[16]][[4]][[1]][,1], # Ca
                          #out.rfmodels[[15]][[4]][[1]][,1], # K
                          #out.rfmodels[[17]][[4]][[1]][,1], # Mg
                          out.rfmodels[[12]][[4]][[1]][,1], # Carbonates
                          out.rfmodels[[5]][[4]][[1]][,1]) # CECef
aggregateRanks(glist = climate.list.perm)[1:10,]$Name

### Climate
setwd(paste0(HomeDir,"Climate/"))
climate.files <- sort(list.files(pattern=".tif"))
clim.rast <- rast(climate.files)

names(clim.rast) <- c("clim_bio1","clim_bio10","clim_bio11",
                      "clim_bio12","clim_bio13","clim_bio14",
                      "clim_bio15","clim_bio16","clim_bio17",
                      "clim_bio18","clim_bio19","clim_bio2",
                      "clim_bio3","clim_bio4","clim_bio5",
                      "clim_bio6","clim_bio7","clim_bio8",
                      "clim_bio9","clim_pet","clim_ari") #Mean monthly potential evapotranspiration
### Mask
clim.rast <- terra::mask(x=clim.rast, mask=eus_buffer) 
plot(clim.rast$clim_ari)
plot(clim.rast$clim_pet)

### Take regular sample
set.seed(2233)
bioclimSample <- terra::spatSample(x = clim.rast, 
                                   size=700000,
                                   method="regular", 
                                   as.df=TRUE, 
                                   xy=TRUE)
### Only complete cases
bioclimSample <- bioclimSample[complete.cases(bioclimSample),]
dim(bioclimSample)

library(Hmisc)
library(dendextend)
library(dendsort)
library(corrplot)

### Which climate variables are correlated?
par(mfrow=c(1,1))
bioclimSample[,3:ncol(bioclimSample)] %>%
  cor(., use = "pairwise.complete.obs") %>%
  corrplot.mixed(.,upper = "ellipse", 
                 lower = "number", 
                 number.cex=0.7, tl.cex=0.6, tl.col = "black")

### checking the correlations are significant
testRes = cor.mtest(bioclimSample[,3:ncol(bioclimSample)], conf.level = 0.95)
corr_bioclim = rcorr(as.matrix(bioclimSample[,3:ncol(bioclimSample)]))
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

bioclimSample[,3:ncol(bioclimSample)] %>%
  cor(., use = "pairwise.complete.obs") %>%
  corrplot.mixed(.,upper = "ellipse", 
                 lower = "number", order="hclust",
                 number.cex=0.7, 
                 tl.cex=0.6, tl.col = "black")

### Highlight only those with correlation above a certain values
M0.85 <- M # Copy matrix
M0.85[ M0.85 < 0.85 & M > -0.85 ] = 0
corrplot(M0.85)
corrplot(M0.85, 
         diag=FALSE,
         type = "upper", 
         order = "hclust", 
         method = "number",
         number.cex=0.6, 
         tl.cex=0.6,
         p.mat = p_mat, 
         sig.level = 0.05)

### Highlight only those with correlation above a certain values
M0.80 <- M # Copy matrix
M0.80[ M0.80 < 0.8 & M > -0.8 ] = 0
corrplot(M0.80)
corrplot(M0.80, 
         diag=FALSE,
         type = "upper", 
         order = "hclust", 
         method = "number",
         number.cex=0.6, 
         tl.cex=0.6,
         p.mat = p_mat, 
         sig.level = 0.05)

### Highlight only those with correlation above a certain values
M0.75 <- M # Copy matrix
M0.75[ M0.75 < 0.75 & M > -0.75 ] = 0
corrplot(M0.75)
corrplot(M0.75, 
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
BioclimSubsetSample <- bioclimSample[,c("clim_bio1","clim_bio4","clim_bio5", 
                                         "clim_bio12", "clim_bio15")]
corrplotsPairs <- GGally::ggpairs(data=BioclimSubsetSample)

### Plot subset of variables based on correlation threshold 0.7 and ranking of RF models
BioclimSub <- terra::subset(clim.rast, c(1,4,7,14,15))
plot(BioclimSub)

save.image("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/4.CovariateSelection/12022024.RData")

# 5. Relief subset -----------------------------------------------------------------

#### Create Heatmap figure
### Reorder covariates by soil-forming factor
library(tidyverse)

### Relief vars
relief_vars
rownames(relief.dfs.permutation) <- relief.dfs.permutation$covariates
### eliminate soil depth
relief.dfs.permutation <- relief.dfs.permutation[!rownames(relief.dfs.permutation)=="Layer_depth",]
### Arrange by variable name
relief.dfs.permutation <- relief.dfs.permutation %>% arrange(factor(covariates, levels = c(relief_vars)))
varImp <-relief.dfs.permutation

### Order the rows based on stable/dynamic properties and number of observations
colnames(soil.scorpan.summary)

soil.order <- c("Clay","Sand","Silt","Carbonates","CECef","CoarseFrag","CEC","CoarseFrag_Grav",
                "TOC","pH","Ca_ppm","K_ppm","TN_ppm","Mg_ppm","EC","P_ppm","Bulk_Density","Na_ppm")          
#varImp <- dplyr::left_join(data.frame(covariates=c(clim_vars,"Layer_depth")),climate.dfs.permutation, by="covariates")

#library(ggplot2)
library(gplots)
library(tidyverse)
heatmap(t(as.matrix(varImp[,soil.order])),
        cexRow = 0.6, 
        cexCol = 0.8,
        margins = c(10, 6))

mypal <-viridis_pal(option = "A",direction=-1)(20)

heatmap.2(t(as.matrix(varImp[,soil.order])),
          dendrogram = "none",
          Colv = FALSE,Rowv = FALSE,
          cexRow = 0.8, notecol="transparent", tracecol = "transparent", 
          cexCol = 0.8,scale="none",col=mypal,
          margins = c(10, 6),
          key.title="Variable importance")

heatmap.2(t(as.matrix(varImp[,c("Clay","Sand","Silt", "TOC","pH",
                                "Ca_ppm","K_ppm","Mg_ppm",
                                "Carbonates","CECef")])), 
          dendrogram = "none",
          Colv = FALSE, Rowv = FALSE,
          cexRow = 0.8, notecol="transparent", tracecol = "transparent", 
          cexCol = 0.8, scale="none",col=mypal,
          margins = c(10, 6),
          key.title="Variable importance")

### Sort all variables by importance
### Can I extract ranks across serveral dataframes?
### Extract top 10 variables by soil-forming factor and for clay, silt, sand, CEC,...

library(RobustRankAggreg)
### Rank variables - top 10 or top 5 by soil-forming factor
relief.list.perm <- list()

### Merge all dataframes into a single one, by soil-forming factor
### But only of a subset of variables , and top 10 from each
relief.list.perm <- list(out.rfmodels[[1]][[6]][[1]][1:10,1], # sand
                          out.rfmodels[[2]][[6]][[1]][1:10,1], # silt
                          out.rfmodels[[3]][[6]][[1]][1:10,1], # clay
                          out.rfmodels[[9]][[6]][[1]][1:10,1], # pH
                          out.rfmodels[[11]][[6]][[1]][1:10,1], # TOC
                          out.rfmodels[[16]][[6]][[1]][1:10,1], # Ca
                          out.rfmodels[[15]][[6]][[1]][1:10,1], # K
                          out.rfmodels[[17]][[6]][[1]][1:10,1], # Mg
                          out.rfmodels[[12]][[6]][[1]][1:10,1], # Carbonates
                          out.rfmodels[[5]][[6]][[1]][1:10,1]) # CECef
relief.list.perm <- list(out.rfmodels[[1]][[6]][[1]][,1], # sand
                         out.rfmodels[[2]][[6]][[1]][,1], # silt
                         out.rfmodels[[3]][[6]][[1]][,1], # clay
                         out.rfmodels[[9]][[6]][[1]][,1], # pH
                         out.rfmodels[[11]][[6]][[1]][,1], # TOC
                         out.rfmodels[[16]][[6]][[1]][,1], # Ca
                         out.rfmodels[[15]][[6]][[1]][,1], # K
                         out.rfmodels[[17]][[6]][[1]][,1])#, # Mg
                         #out.rfmodels[[12]][[6]][[1]][,1], # Carbonates
                         #out.rfmodels[[5]][[6]][[1]][,1]) # CECef

### Check the graphs with ranking of variable importance
aggregateRanks(glist = relief.list.perm)[1:11,]$Name

### But only of a subset of variables, all of them
relief.list.perm <- list()

relief.list.perm <- list(out.rfmodels[[1]][[6]][[1]][,1], # sand
                          out.rfmodels[[2]][[6]][[1]][,1], # silt
                          out.rfmodels[[3]][[6]][[1]][,1], # clay
                          #out.rfmodels[[9]][[4]][[1]][,1], #pH
                          #out.rfmodels[[11]][[4]][[1]][,1], #TOC
                          out.rfmodels[[16]][[6]][[1]][,1], #Ca
                          #out.rfmodels[[15]][[4]][[1]][,1], # K
                          #out.rfmodels[[17]][[4]][[1]][,1], # Mg
                         # out.rfmodels[[12]][[6]][[1]][,1], # Carbonates
                          out.rfmodels[[5]][[6]][[1]][,1]) # CECef
aggregateRanks(glist = relief.list.perm)[1:10,]$Name

### Plot
out.rfmodels[[1]][[6]][[2]] # sand
out.rfmodels[[2]][[6]][[2]] # silt
out.rfmodels[[3]][[6]][[2]] # clay
out.rfmodels[[5]][[6]][[2]] # CECef
out.rfmodels[[16]][[6]][[2]] # Ca
out.rfmodels[[12]][[6]][[2]] # Carbonates

par(mfrow=c(2,2))
out.rfmodels[[9]][[6]][[2]] # pH
out.rfmodels[[11]][[6]][[2]] # TOC
out.rfmodels[[15]][[6]][[2]] # K
out.rfmodels[[17]][[6]][[2]] # Mg

### Relief
setwd(paste0(HomeDir,"Relief/"))
### Covariates subset
relief.files <- c("mdt_lidar_2017_25m_etrs89.tif", # DEM
                  "slope.tif", "slope_5.tif", "slope_10.tif", # Slope calculated from 3,5,10 cells
                  "easterness.tif", # Easterness 
                  "northerness.tif", # northerness
                  "nor_slope.tif" , # northness x slope
                  "twi_saga.tif", # SAGA Topographic wetness index
                  "p_curv.tif", "pr_curv_5.tif", "pr_curv_10.tif", # profile curvature (3,5,10)
                  "pl_curv_5.tif", "pl_curv_10.tif", # planar curvature (5,10)
                  "lg_curv_5.tif", "lg_curv_10.tif", # longitudinal curvature (5,10)
                  "t_curv.tif", "cs_curv_5.tif", "cs_curv_10.tif" , # tangential curvature (3) cross-sectional curvature (5,10) (same)
                  "mrvbf.tif", "mrrtf.tif", ### MRVBF, MRRTF
                  "slope_height.tif", # Slope Height
                  "norm_height.tif", # Normalized height
                  "st_height.tif", # Standardized Height
                  "valley_depth.tif", # Valley depth
                  "mid_slope.tif",  # Mid slope
                  "tpi_8_3.tif", "tpi_20_5.tif")  # Multiscale topographic position index, calculated with a search radius of 8,20 cells, and 3,5 scales respectively
# "geomorphon.tif", "geomorphon2.tif")  ## Geomorphons (search radius of 11, skip 1 & 3, flat 1 & 1.5)

relief.rast <- terra::rast(relief.files)
plot(relief.rast)
names(relief.rast)[names(relief.rast)=="mdt_lidar_2017_25m_etrs89"] <- "dem"
names(relief.rast)[names(relief.rast)=="twi_saga"] <- "twi"
names(relief.rast)[names(relief.rast)=="slope@relief"] <- "slope"
names(relief.rast)[names(relief.rast)=="elevation_easterness@relief"] <- "easterness"
names(relief.rast)[names(relief.rast)=="elevation_northerness@relief"] <- "northerness"
names(relief.rast)[names(relief.rast)=="elevation_northerness_slope@relief"] <- "north_slope"
names(relief.rast)[names(relief.rast)=="t_curv@relief"] <- "t_curv"
names(relief.rast)[names(relief.rast)=="p_curv@relief"] <- "p_curv"
names(relief.rast)[names(relief.rast)=="mrvbf@relief"] <- "mrvbf"
names(relief.rast)[names(relief.rast)=="mrrtf@relief"] <- "mrrtf"
# ### Mask
# relief.rast <- terra::mask(x=relief.rast, mask=eus_buffer) 

### Take regular sample
set.seed(2233)
ReliefSample <- terra::spatSample(x = relief.rast, 
                                   size=700000,
                                   method="regular", 
                                   as.df=TRUE, 
                                   xy=TRUE)
### Only complete cases
ReliefSample <- ReliefSample[complete.cases(ReliefSample),]
dim(ReliefSample)

library(Hmisc)
library(dendextend)
library(dendsort)
library(corrplot)

### Which climate variables are correlated?
par(mfrow=c(1,1))
ReliefSample[,3:ncol(ReliefSample)] %>%
  cor(., use = "pairwise.complete.obs") %>%
  corrplot.mixed(.,upper = "ellipse", 
                 lower = "number", 
                 number.cex=0.7, tl.cex=0.6, tl.col = "black")

### checking the correlations are significant
testRes = cor.mtest(ReliefSample[,3:ncol(ReliefSample)], conf.level = 0.95)
corr_Relief = rcorr(as.matrix(ReliefSample[,3:ncol(ReliefSample)]))
M <- corr_Relief$r
p_mat <- corr_Relief$P

corrplot(M, 
         diag=FALSE,
         type = "upper", 
         order = "hclust", 
         method = "number",
         number.cex=0.6, 
         tl.cex=0.6,
         # p.mat = p_mat, 
         sig.level = 0.05)

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
M0.85 <- M # Copy matrix
M0.85[ M0.85 < 0.85 & M > -0.85 ] = 0
corrplot(M0.85)
corrplot(M0.85, 
         diag=FALSE,
         type = "upper", 
         order = "hclust", 
         method = "number",
         number.cex=0.6, 
         tl.cex=0.6,
        # p.mat = p_mat, 
         sig.level = 0.05)

### Highlight only those with correlation above a certain values
M0.80 <- M # Copy matrix
M0.80[ M0.80 < 0.8 & M > -0.8 ] = 0
corrplot(M0.80)
corrplot(M0.80, 
         diag=FALSE,
         type = "upper", 
         order = "hclust", 
         method = "number",
         number.cex=0.6, 
         tl.cex=0.6,
         #p.mat = p_mat, 
         sig.level = 0.05)

### Highlight only those with correlation above a certain values
M0.75 <- M # Copy matrix
M0.75[ M0.75 < 0.75 & M > -0.75 ] = 0
corrplot(M0.75)
corrplot(M0.75, 
         diag=FALSE,
         type = "upper", 
         order = "hclust", 
         method = "number",
         number.cex=0.6, 
         tl.cex=0.6,
        # p.mat = p_mat, 
         sig.level = 0.05)

rm(M,M0.75,p_mat,corr_Relief)

### Hard choice
### Perform and RDA and check for variables explaining variance
SoilRelief <- soil.scorpan[, c("Clay","Sand","Silt","TOC","pH",
                               relief_vars, "Layer_depth")]
library(vegan)
SoilRelief <- SoilRelief[complete.cases(SoilRelief),]
dim(SoilRelief)

### condition with depth
SoilRelief$Layer_depth <- as.factor(SoilRelief$Layer_depth)

### Search the optimal variables with RDA
X = SoilRelief[,c("Clay","Sand","Silt","TOC","pH")] #,"Ca_ppm","K_ppm","Mg_ppm")]
mod0 <- rda(X ~ 1, SoilRelief[,c(relief_vars, "Layer_depth")])  # Model with intercept only
mod1 <- rda(X ~ ., SoilRelief[,c(relief_vars, "Layer_depth")])  # Model with all explanatory variables

## With scope present, the default direction is "both"
set.seed(3344)
step.res <- ordiR2step(mod0, scope = formula(mod1), perm.max = 5000, direction="both")
step.res$anova
plot(step.res)

SoilRelief.rds <- rda(Y = SoilRelief[,c("Clay","Sand","Silt","TOC","pH")],
                      X = SoilRelief[,c("dem","twi","slope",
                                        "mid_slope","north_slope","valley_depth",
                                        "easterness","northerness")],
                      Z = SoilRelief$Layer_depth)
plot(SoilRelief.rds)


save.image("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/4.CovariateSelection/14022024.RData")
# load("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/4.CovariateSelection/14022024.RData")
### end of the script