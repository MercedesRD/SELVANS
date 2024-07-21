### Number of pedogenon classes for France.
### Date: 30/11/2023
### Author: Mercedes Roman Dobarco

### Objective: How many unique combinations between 
### soil type, geological class and climate type can we find in France?
### This will give us an initial number of Pedogenon classes for 1000 BCE

### load libraries
library(terra)
#library(raster)

setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/")

### Load soil map 1:1000000
soil <- terra::rast("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Soil/soil1.tif")
plot(soil)
ext(soil)
unique(soil)
soilL93 <- project(soil, "EPSG:2154", method="near")
### Reclass 0 as NA
f <- function(x) ifelse(x == 0, NA, x)
soilL93.r <- app(soilL93, f)
plot(soilL93.r)

### Parent material 1:1000000
pmat <- terra::rast("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/ParentMaterial/mat11.tif")
plot(pmat)
### Transform
pmatL93 <- project(pmat, "EPSG:2154", method="near")
### Reclass 0 as NA
f <- function(x) ifelse(x == 0, NA, x)
pmatL93.r <- app(pmatL93, f)
plot(pmatL93.r)

## create raster stack
combi <- c(soilL93.r, pmatL93.r)
plot(combi)

### Get number of unique combinations
combi_cat <- terra::unique(combi, incomparables =FALSE,
                           na.rm=TRUE, as.raster=TRUE, 
                           digits=0)
plot(combi_cat)
my_levels <- levels(combi_cat)[[1]]
length(unique(my_levels$label)) ### 45 combinations

### climate type
climate <- terra::rast("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Climate/ClimateType/TYPO_RGF93.tif")
plot(climate)
ext(climate)
crs("epsg:2154", describe=TRUE)
### Resample to the resolution and extent of the other two raster files
climate90 <- terra::resample(x=climate, y=soilL93, method="near")
writeRaster(climate90, 
            filename = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Climate/ClimateType/climtype90.tif")

### Reclass 9 as NA
f <- function(x) ifelse(x == 9, NA, x)
climate90.r <- app(climate90, f)
plot(climate90.r)

plot(combi_cat)
my_levels <- levels(combi_cat)[[1]]
length(unique(my_levels$label)) 

combi2 <- c(climate90.r, combi_cat)
combi2_cat <- terra::unique(combi2, incomparables = FALSE,
                            na.rm=TRUE, as.raster=TRUE, 
                            digits=0)

plot(combi2_cat)
my_levels2 <- levels(combi2_cat)[[1]]
length(unique(my_levels2$label)) ### 296 combinations

levels2.na <- grep(pattern="NA", x=my_levels2$label)
whole_levels2 <- my_levels2$label[-levels2.na]
length(whole_levels2) ### We get 264 levels

### Combine soil, pmat, and climate into a unique identifier
### simply...
climate1 <- climate90.r *1000
pmat1 <- pmatL93.r * 100
combis <- climate1 + pmat1 + soilL93.r
length(unique(values(combis)[!is.nan(values(combis))]))### We get 264 levels

# ## create raster stack
# combi <- c(soilL93, pmatL93)
# 
# ### Get number of unique combinations
# combi_cat <- terra::unique(combi, incomparables =TRUE,
#                            na.rm=TRUE, as.raster=TRUE, 
#                            digits=0)
# plot(combi_cat)
# my_levels <- levels(combi_cat)[[1]]
# length(unique(my_levels$label)) 
# 
# combi2 <- c(climate90, combi_cat)
# combi2_cat <- terra::unique(combi2, incomparables =TRUE,
#                             na.rm=TRUE, as.raster=TRUE, 
#                             digits=0)
# 
# plot(combi2_cat)
# my_levels2 <- levels(combi2_cat)[[1]]
# length(unique(my_levels2$label)) 
# 
# levels2.na <- grep(pattern="NA", x=my_levels2$label)
# whole_levels2 <- my_levels2$label[-levels2.na]
# length(whole_levels2)
# n0.levels <- grep(pattern="0_0", x=whole_levels)
# whole_levels <- whole_levels[-n0.levels]
# ### 285 combinations
# 
# ### Get number of unique combinations
# combi_cat3 <- terra::unique(combi, incomparables =FALSE,
#                            na.rm=TRUE, as.raster=TRUE, 
#                            digits=0)
# plot(combi_cat3)
# my_levels3 <- levels(combi_cat3)[[1]]
# length(unique(my_levels3$label)) ### 46 combinations
# 
# combi4 <- c(climate90, combi_cat3)
# combi4_cat <- terra::unique(combi4, incomparables =FALSE,
#                             na.rm=TRUE, as.raster=TRUE, 
#                             digits=0)
# 
# plot(combi4_cat)
# my_levels4 <- levels(combi4_cat)[[1]]
# length(unique(my_levels4$label)) 
# 
# levels2.na <- grep(pattern="NA", x=my_levels4$label)
# whole_levels4 <- my_levels4$label[-levels2.na]
# length(whole_levels4)
# n0.levels <- grep(pattern="_0_0", x=whole_levels4)
# whole_levels4 <- whole_levels4[-n0.levels]
# ### We get 276 levels
# 
# ### Eliminate level 9 climate class
# #n9.levels <- grep(pattern="9_", x=whole_levels4)
# whole_levels4 <- whole_levels4[-c(265:276)]
# length(whole_levels4)
# ### We get 264 levels
# 
# save.image("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/Rsessions/combis.RData")
# load("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/Rsessions/combis.RData")

# Environmental Zones Europe Metzger et al. 2005 --------------------------

library(sf)
EnSv8 <- st_read("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/EnvSEv8/EnSv8/EnSv8/ens_v8.shp")
plot(EnSv8["EnZ"])

### Project
EnSv8 <- st_transform(EnSv8,2154)
### Crop
EnSv8 <- st_crop(x=EnSv8, y=st_bbox(soilL93))
### Rasterize
EnSv8.Fr <- terra::rasterize(x=EnS_name, 
                             y=soilL93, 
                             field="EnZ", 
                             fun="min")
### Less categories
EnZ.Fr <- terra::rasterize(x=EnSv8, 
                           y=soilL93, 
                           field="EnZ", 
                           fun="min")

### Mask to the same surface as soil and pmat
EnZ.Fr <- mask(EnZ.Fr, mask = soilL93.r)
plot(EnZ.Fr)

### Now combine
combi3 <- c(EnZ.Fr, combi_cat)
combi3_cat <- terra::unique(combi3, incomparables = FALSE,
                            na.rm=TRUE, as.raster=TRUE, 
                            digits=0)

plot(combi3_cat)
my_levels3 <- levels(combi3_cat)[[1]]
length(unique(my_levels3$label)) ### 180 combinations

levels3.na <- grep(pattern="NA", x=my_levels3$label)
whole_levels3 <- my_levels3$label[-levels3.na]
length(whole_levels3) ### We get 167 levels

save.image("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/Rsessions/combis.RData")

### end of the script