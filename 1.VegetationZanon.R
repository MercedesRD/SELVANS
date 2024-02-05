#### Pedogenon maps for France

### Date 29/11/2023
### Author: Mercedes Roman Dobarco
### Objective: Crop Zanon et al. (2018) data to the extent of France

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
france_buffer_WGS84 <- st_read("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Administrative/france_bufferWGS84.shp")

# Import Zanon et al. (2018) data from 1000BCE to latest available --------------------------------------------

### The fossil dataset was divided into 49 timeslices ranging from 12,000 to 0 
### calibrated years BP (years before AD 1950; hereafter BP). 
### Each timeslice covers a 250-year window with the exception 
### of the most recent one, which is asymmetric as it cannot project 
### into the future, and therefore covers an interval of 185 years 
### centered on 0 BP.

### therefore 950 BCE --> 2900 BP it is comprised by the layer 3000 BP (3125 BP - 2875 BP)

### Import forest cover for 3000 BP
setwd("C:/Covariates/Europe/Zanon_et_al_2018_maps/forest_cover/")
forest_cover_3000 <- terra::rast("forest_cover_3000.grd")
plot(forest_cover_3000)
plot(france, add=TRUE)
st_crs(forest_cover_3000)
st_crs(france_buffer_WGS84)

### Project France buffer to same CRS
france_buffer_9122 <- transform(france_buffer_WGS84,9122)
### Transform to spatvector
france_v <- vect(france_buffer_9122)
#### Crop area of France
forest_cover_3000_Fr <- terra::crop(forest_cover_3000,france_v)
plot(forest_cover_3000_Fr)
plot(france_v,add=TRUE)
### Write to file
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/vegetation/Zanon/")
writeRaster(forest_cover_3000_Fr,filename="forest_cover_3000.tif")

### Repeat for other time steps
setwd("C:/Covariates/Europe/Zanon_et_al_2018_maps/forest_cover/")
list.files()
### Names of the directories of interest
varz <- c("forest_cover_0.grd","forest_cover_250.grd","forest_cover_500.grd",
          "forest_cover_750.grd","forest_cover_1000.grd","forest_cover_1250.grd",
          "forest_cover_1500.grd","forest_cover_1750.grd","forest_cover_2000.grd",
          "forest_cover_2250.grd","forest_cover_2500.grd","forest_cover_2750.grd")
names.out <- gsub(pattern=".grd", replacement = ".tif", x=varz)

### Crop to France and save as tif file
for(i in 1:length(varz)) {
  print(paste0("working on layer ",varz[[i]]))
  setwd("C:/Covariates/Europe/Zanon_et_al_2018_maps/forest_cover/")
  ### Load
  forest_cover <- terra::rast(varz[[i]])
  ### Crop to the extent of France
  forest_cover_Fr <- terra::crop(forest_cover,france_v)
    ### Write Raster
  setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/vegetation/Zanon/")
  writeRaster(forest_cover_Fr,
              filename=names.out[[i]],
              overwrite=TRUE)
}

### end of the script