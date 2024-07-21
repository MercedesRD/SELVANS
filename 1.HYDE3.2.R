#### Pedogenon maps for France

### Date 29/11/2023
### Author: Mercedes Roman Dobarco
### Objective: Crop HYDE 3.2 data to the extent of France

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
library(GGally)

### 1. Study area -----------------------------------------------------------

### I bring shapefile
france_regions <- read_sf("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Administrative/regions-20180101-shp/regions-20180101.shp")

### Only metropolitan France
france_regions <- france_regions[france_regions$code_insee %in% c(11:93),]
plot(france_regions["nom"])
france <- st_union(france_regions)
plot(france)
st_crs(france)
st_bbox(france)
rm(france_regions)

# ### Project to RGF93 v2b / Lambert-93
# franceL93 <- st_transform(france,9794)
# 
# ### Create a buffer around france
# sf_use_s2(TRUE)
# france_buffer <- sf::st_buffer(franceL93, dist = 10000)
# par(mfrow=c(1,1))
# plot(france_buffer)
# plot(franceL93,col='blue',add=TRUE)
# st_write(france_buffer, "C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Administrative/france_bufferL93.shp")
# 
# ### Project back to WGS84
# france_buffer_WGS84 <- st_transform(france_buffer,4326)
# st_write(france_buffer_WGS84, "C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Administrative/france_bufferWGS84.shp")
france_buffer_WGS84 <- st_read("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Administrative/france_bufferWGS84.shp")


# Import HYDE 3.2 general area --------------------------------------------

### Import Maximum land area 
maxln <- terra::rast("D:/HYDE/HYDE3.2/general_files/general_files/maxln_cr.asc")
st_crs(maxln)
st_crs(france_buffer_WGS84)
### Transform to spatvector
france_v <- vect(france_buffer_WGS84)
### Assign correct CRS to HYDE 3.2
crs(maxln)  <- "epsg:4326"

#### Crop area of France
maxln_Fr <- terra::crop(maxln,france_v)
plot(maxln_Fr)
plot(france_v,add=TRUE)
### Write to file
writeRaster(maxln_Fr,filename="D:/HYDE/HYDE3.2/France/maxln_Fr.tif")
maxln_Fr <- rast("D:/HYDE/HYDE3.2/France/maxln_Fr.tif")

# Import HYDE 3.2 cropping and grazing for 8000 BCE until 700 CE --------------------------------------------

### Import HYDE 3.2 data
baseline.dir <- "D:/HYDE/HYDE3.2/baseline/baseline/zip/"
setwd(baseline.dir)

### Names of the directories of interest
varz <- c("8000BC_lu","7000BC_lu","6000BC_lu","5000BC_lu",
          "4000BC_lu","3000BC_lu","2000BC_lu","1000BC_lu",
          "0AD_lu","100AD_lu","200AD_lu","300AD_lu","400AD_lu",
          "500AD_lu","600AD_lu","700AD_lu","800AD_lu")

### List ALL files for land use, I will select dates later
varz <- gsub(x=list.files(pattern="_lu"), pattern=".zip", replacement = "")

### Create directories
unzip.dir <- "D:/HYDE/HYDE3.2/unzip/"
dir.create(unzip.dir)
setwd(unzip.dir)
for(i in 1:length(varz)) {
  dir.create(paste0("h", varz[[i]]))
}

### Unzip to each file
for(i in 1:length(varz)) {
  unzip(zipfile = paste0(baseline.dir,varz[[i]],".zip"),
        exdir = paste0(unzip.dir,"h", varz[[i]]))
}

### Select HYDE 3.2 variables
### cropping
### grazing
### total irrigation
HydeVars <- c("cropland","grazing","pasture","rangeland",
              "conv_rangeland","tot_irri","tot_rainfed","tot_rice",
              "ir_norice","rf_norice","ir_rice","rf_rice")
#HydeVars <- c("cropland","grazing","tot_irri")  
              ###"working on time slice 1890AD_lu"
### Crop to Fr###ance and save as tif file
for(i in 31:length(varz)) {
  print(paste0("working on time slice ",varz[[i]]))
  setwd(paste0(unzip.dir,"h", varz[[i]]))
  # list files
  rast.files <- list.files(pattern=".asc")
  ### stack
  lustack <- rast(rast.files)
  ### Crop to the extent of France
  lustack_Fr <- terra::crop(lustack,france_v)
  ### Divide by area of land
  percfun=function(x,y){return(x/y *100)}
  #lu_relj.s <- list()
  for(j in 1:12){
    lu_relj <- percfun(lustack_Fr[[j]], maxln_Fr)
    writeRaster(lu_relj,filename=paste0(names(lu_relj),"_Fr.tif"),
                overwrite=TRUE)
    #lu_relj.s[[j]] <- lu_relj
  }
  #lu_relj.s <- rast(lu_relj.s)
}

### end of the script