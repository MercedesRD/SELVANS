### Timeslices - preparation of files for each timeID

### Date 19/03/2024
### Author: Mercedes Roman Dobarco
### Objective: Prepare the HYDE3.2 data layer for the period 1000 BCE - 1 CE at thedesired temporal resolution

# analysis packages
library(dplyr)
library(tidyverse)
library(sf)
library(terra)

# visualization packages
library(ggplot2)
library(viridis)
library(scales)
library(viridisLite)
library(rasterVis)


# 1. HYDE - 1000 BCE - 0CE ---------------------------------------------------

### We have HYDE for 1000 BCE - 1 BCE (millenial stime step), e.g.  grazing1000BC_Fr.tif -- timeID_0
### And we have HYDE 0 CE -	99 CE (centennial time step), e.g., grazing0AD_Fr.tif -- -- timeID_10

### set output wd
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices")
timeSteps <- c(0:10) 
for(i in 1:length(timeSteps)){dir.create(paste0("TimeID_",timeSteps[[i]]))}

increment.hyde.1000BCE_0CE.func <- function(TimeID_10,TimeID_0,type_increment, name) {
  
  ### Make a SpatRaster
  input <- c(TimeID_0,TimeID_10)
  
  ### Prepare input raster stacks
  total_change <-  terra::lapp(input, fun=function(x,y){return(y-x)})
  b  <-  terra::lapp(input, fun=function(x,y){return((y/x)^0.1)})
  
  input2 <- c(TimeID_0,total_change)
  input3 <- c(TimeID_0, b)
  
  if (type_increment=="linear") {
    
    TimeID_1 = terra::lapp(input2, fun=function(x,y){return(x+(0.1*y))},
                              filename=paste0("TimeID_1/",name,"_timeID_1.tif"), 
                              overwrite=TRUE)
    TimeID_2 = terra::lapp(input2, fun=function(x,y){return(x+(0.2*y))},
                              filename=paste0("TimeID_2/",name,"_timeID_2.tif"), 
                              overwrite=TRUE)
    TimeID_3 = terra::lapp(input2, fun=function(x,y){return(x+(0.3*y))},
                              filename=paste0("TimeID_3/",name,"_timeID_3.tif"), 
                              overwrite=TRUE)
    TimeID_4 = terra::lapp(input2, fun=function(x,y){return(x+(0.4*y))},
                              filename=paste0("TimeID_4/",name,"_timeID_4.tif"), 
                              overwrite=TRUE)
    TimeID_5 = terra::lapp(input2, fun=function(x,y){return(x+(0.5*y))},
                              filename=paste0("TimeID_5/",name,"_timeID_5.tif"), 
                              overwrite=TRUE)
    TimeID_6 = terra::lapp(input2, fun=function(x,y){return(x+(0.6*y))},
                              filename=paste0("TimeID_6/",name,"_timeID_6.tif"), 
                              overwrite=TRUE)
    TimeID_7 = terra::lapp(input2, fun=function(x,y){return(x+(0.7*y))},
                              filename=paste0("TimeID_7/",name,"_timeID_7.tif"), 
                              overwrite=TRUE)
    TimeID_8 = terra::lapp(input2, fun=function(x,y){return(x+(0.8*y))},
                              filename=paste0("TimeID_8/",name,"_timeID_8.tif"), 
                              overwrite=TRUE)
    TimeID_9 = terra::lapp(input2, fun=function(x,y){return(x+(0.9*y))},
                              filename=paste0("TimeID_9/",name,"_timeID_9.tif"), 
                              overwrite=TRUE)
    
  } else if (type_increment=="exponential"){
    
    TimeID_1 = terra::lapp(input3, fun=function(x,y){return(x*(y^1))},
                                filename=paste0("TimeID_1/",name,"_timeID_1.tif"), 
                                overwrite=TRUE)
    TimeID_2 = terra::lapp(input3, fun=function(x,y){return(x*(y^2))},
                              filename=paste0("TimeID_2/",name,"_timeID_2.tif"), 
                              overwrite=TRUE)
    TimeID_3 = terra::lapp(input3, fun=function(x,y){return(x*(y^3))},
                              filename=paste0("TimeID_3/",name,"_timeID_3.tif"), 
                              overwrite=TRUE)
    TimeID_4 = terra::lapp(input3, fun=function(x,y){return(x*(y^4))},
                              filename=paste0("TimeID_4/",name,"_timeID_4.tif"), 
                              overwrite=TRUE)
    TimeID_5 = terra::lapp(input3, fun=function(x,y){return(x*(y^5))},
                              filename=paste0("TimeID_5/",name,"_timeID_5.tif"), 
                              overwrite=TRUE)
    TimeID_6 = terra::lapp(input3, fun=function(x,y){return(x*(y^6))},
                              filename=paste0("TimeID_6/",name,"_timeID_6.tif"), 
                              overwrite=TRUE)
    TimeID_7 = terra::lapp(input3, fun=function(x,y){return(x*(y^7))},
                              filename=paste0("TimeID_7/",name,"_timeID_7.tif"), 
                              overwrite=TRUE)
    TimeID_8 = terra::lapp(input3, fun=function(x,y){return(x*(y^8))},
                              filename=paste0("TimeID_8/",name,"_timeID_8.tif"), 
                              overwrite=TRUE)
    TimeID_9 = terra::lapp(input3, fun=function(x,y){return(x*(y^9))},
                              filename=paste0("TimeID_9/",name,"_timeID_9.tif"), 
                              overwrite=TRUE)

  }
  
  ### return raster stack
  hyde.gradual <- c(TimeID_1,TimeID_2,TimeID_3,TimeID_4,TimeID_5,TimeID_6,TimeID_7,TimeID_8,TimeID_9)
  plot(hyde.gradual)
  return(hyde.gradual)
  
}

# TimeID10 <- 25
# TimeID0 <- 0.6
# #TimeID10= a*(b^10)
# #TimeID0= a*(b^0)
# #TimeID10/TimeID0 = (b^10)
# b=(TimeID10/TimeID0)^(1/10)
# a=TimeID0
# exponential.inc <- function(x) {a*(b^x)}
# y <- exponential.inc(x)
# x <- 0:10
# plot(x,y)

### Time ID 0
### Copy and rename files
setwd("D:/HYDE/HYDE3.2/unzip/h1000BC_lu/")
files_1000BCE <- list.files(pattern=".tif")

file.copy(from = paste0("D:/HYDE/HYDE3.2/unzip/h1000BC_lu/",files_1000BCE),
          to = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_0/", files_1000BCE))
file.rename(from = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_0/", files_1000BCE),
            to = gsub(x = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_0/", files_1000BCE),
                      pattern="1000BC_Fr",
                      replacement="_timeID_0"))

### Time ID 10
### Copy and rename files
setwd("D:/HYDE/HYDE3.2/unzip/h0AD_lu/")
files_0CE <- list.files(pattern=".tif")

file.copy(from = paste0("D:/HYDE/HYDE3.2/unzip/h0AD_lu/",files_0CE),
          to = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_10/", files_0CE))
file.rename(from = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_10/", files_0CE),
            to = gsub(x = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_10/", files_0CE),
                      pattern="0AD_Fr",
                      replacement="_timeID_10"))

### Process HYDE variables for time steps in between
### Select HYDE 3.2 variables
### cropping
### grazing
### total irrigation
HydeVars <- c("cropland","grazing","pasture","rangeland",
              "conv_rangeland","tot_irri","tot_rainfed","tot_rice",
              "ir_norice","rf_norice","ir_rice","rf_rice")



### Cropland
cropland_timeID_0 <- rast("./TimeID_0/cropland_timeID_0.tif")
cropland_timeID_10 <- rast("./TimeID_10/cropland_timeID_10.tif")
plot(cropland_timeID_0)
plot(cropland_timeID_10)

cropland_1000BCE_0CE <- increment.hyde.1000BCE_0CE.func(TimeID_10 = cropland_timeID_10, 
                                                        TimeID_0 = cropland_timeID_0,
                                                        type_increment = "linear", 
                                                        name = "cropland_l")

cropland_1000BCE_0CE <- increment.hyde.1000BCE_0CE.func(TimeID_10 = cropland_timeID_10, 
                                                        TimeID_0 = cropland_timeID_0,
                                                        type_increment = "exponential", 
                                                        name = "cropland_e")

cropland <- c(cropland_timeID_0,cropland_1000BCE_0CE,cropland_timeID_10)
names(cropland) <- c("crop2900BP", "crop2800BP", "crop2700BP", "crop2600BP", "crop2500BP", "crop2400BP",
                     "crop2300BP", "crop2200BP", "crop2100BP", "crop2000BP", "crop1900BP")
plot(cropland)
rasterVis::levelplot(cropland, par.settings = viridisTheme)

### grazing
grazing_timeID_0 <- rast("./TimeID_0/grazing_timeID_0.tif")
grazing_timeID_10 <- rast("./TimeID_10/grazing_timeID_10.tif")
plot(grazing_timeID_0)
plot(grazing_timeID_10)

grazing_1000BCE_0CE <- increment.hyde.1000BCE_0CE.func(TimeID_10 = grazing_timeID_10, 
                                                        TimeID_0 = grazing_timeID_0,
                                                        type_increment = "linear", 
                                                        name = "grazing_l")

grazing_1000BCE_0CE <- increment.hyde.1000BCE_0CE.func(TimeID_10 = grazing_timeID_10, 
                                                        TimeID_0 = grazing_timeID_0,
                                                        type_increment = "exponential", 
                                                        name = "grazing_e")

grazing <- c(grazing_timeID_0,grazing_1000BCE_0CE,grazing_timeID_10)
names(grazing) <- c("graz2900BP", "graz2800BP", "graz2700BP", "graz2600BP", "graz2500BP", "graz2400BP",
                     "graz2300BP", "graz2200BP", "graz2100BP", "graz2000BP", "graz1900BP")
plot(grazing)
rasterVis::levelplot(grazing, par.settings = viridisTheme)


### pasture
pasture_timeID_0 <- rast("./TimeID_0/pasture_timeID_0.tif")
pasture_timeID_10 <- rast("./TimeID_10/pasture_timeID_10.tif")
plot(pasture_timeID_0)
plot(pasture_timeID_10)

pasture_1000BCE_0CE <- increment.hyde.1000BCE_0CE.func(TimeID_10 = pasture_timeID_10, 
                                                       TimeID_0 = pasture_timeID_0,
                                                       type_increment = "linear", 
                                                       name = "pasture_l")
pasture <- c(pasture_timeID_0,pasture_1000BCE_0CE,pasture_timeID_10)
names(pasture) <- c("pasture2900BP", "pasture2800BP", "pasture2700BP", "pasture2600BP", "pasture2500BP", "pasture2400BP",
                    "pasture2300BP", "pasture2200BP", "pasture2100BP", "pasture2000BP", "pasture1900BP")
plot(pasture)
rasterVis::levelplot(pasture, par.settings = viridisTheme)

# pasture_1000BCE_0CE <- increment.hyde.1000BCE_0CE.func(TimeID_10 = pasture_timeID_10, 
#                                                        TimeID_0 = pasture_timeID_0,
#                                                        type_increment = "exponential", 
#                                                        name = "pasture_e")
# 
# pasture <- c(pasture_timeID_0,pasture_1000BCE_0CE,pasture_timeID_10)
# names(pasture) <- c("2900BP", "2800BP", "2700BP", "2600BP", "2500BP", "2400BP",
#                     "2300BP", "2200BP", "2100BP", "2000BP", "1900BP")
# plot(pasture)
### Exponential change does not work well because the term "a" ia 0 and everything is predicted as 0

### rangeland
# rangeland_timeID_0 <- rast("./TimeID_0/rangeland_timeID_0.tif")
# rangeland_timeID_10 <- rast("./TimeID_10/rangeland_timeID_10.tif")
# plot(rangeland_timeID_0)
# plot(rangeland_timeID_10)
# 
# rangeland_1000BCE_0CE <- increment.hyde.1000BCE_0CE.func(TimeID_10 = rangeland_timeID_10, 
#                                                        TimeID_0 = rangeland_timeID_0,
#                                                        type_increment = "linear", 
#                                                        name = "rangeland_l")
# rangeland <- c(rangeland_timeID_0,rangeland_1000BCE_0CE,rangeland_timeID_10)
# names(rangeland) <- c("2900BP", "2800BP", "2700BP", "2600BP", "2500BP", "2400BP",
#                     "2300BP", "2200BP", "2100BP", "2000BP", "1900BP")
# plot(rangeland)

### conv_rangeland
conv_rangeland_timeID_0 <- rast("./TimeID_0/conv_rangeland_timeID_0.tif")
conv_rangeland_timeID_10 <- rast("./TimeID_10/conv_rangeland_timeID_10.tif")
plot(conv_rangeland_timeID_0)
plot(conv_rangeland_timeID_10)

conv_rangeland_1000BCE_0CE <- increment.hyde.1000BCE_0CE.func(TimeID_10 = conv_rangeland_timeID_10, 
                                                         TimeID_0 = conv_rangeland_timeID_0,
                                                         type_increment = "linear", 
                                                         name = "conv_rangeland_l")
conv_rangeland <- c(conv_rangeland_timeID_0,conv_rangeland_1000BCE_0CE,conv_rangeland_timeID_10)
names(conv_rangeland) <- c("conv_2900BP", "conv_2800BP", "conv_2700BP", "conv_2600BP", "conv_2500BP",
                           "conv_2400BP", "conv_2300BP", "conv_2200BP", "conv_2100BP", "conv_2000BP", 
                           "conv_1900BP")
plot(conv_rangeland)
rasterVis::levelplot(conv_rangeland, par.settings = viridisTheme)

### Other Hyde variables
otherHydeVars <- c("tot_irri","tot_rainfed","tot_rice",
              "ir_norice","rf_norice","ir_rice","rf_rice")

for(i in 1:length(otherHydeVars)){
  print(i) 
  ### X variable
  timeID_0 <- rast(paste0("./TimeID_0/",otherHydeVars[[i]],"_timeID_0.tif"))
  timeID_10 <- rast(paste0("./TimeID_10/",otherHydeVars[[i]],"_timeID_10.tif"))
  plot(timeID_0)
  plot(timeID_10)
  
  r_1000BCE_0CE <- increment.hyde.1000BCE_0CE.func(TimeID_10 = timeID_10, 
                                                   TimeID_0 = timeID_0,
                                                   type_increment = "linear", 
                                                   name = paste0(otherHydeVars[[i]],"_l"))
  r_stack <- c(timeID_0,r_1000BCE_0CE,timeID_10)
  names(r_stack) <- paste0(otherHydeVars[[i]],c("_2900BP", "_2800BP", "_2700BP", "_2600BP", "_2500BP",
                             "_2400BP", "_2300BP", "_2200BP", "_2100BP", "_2000BP", 
                             "_1900BP"))
  plot(r_stack)
  print(rasterVis::levelplot(r_stack, par.settings = viridisTheme))
  
}

#########################################################################################################

#  2. Climate  - 1000 BCE - 0CE ----------------------------------------------------------------

### Copy and rename files
setwd("D:/FRANCE/Covariates/Climate/TimeID_-9/")
files_1000BCE <- list.files(pattern="Fr.tif")

file.copy(from = paste0("./",files_1000BCE),
          to = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_0/", files_1000BCE))

file.rename(from = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_0/", files_1000BCE),
            to = gsub(x = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_0/", files_1000BCE),
                      pattern="CHELSA_TraCE21k_",
                      replacement="CHELSA_"))
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_0/")
files_1000BCE <- list.files(pattern="CHELSA")
file.rename(from = paste0("./", files_1000BCE),
            to = gsub(x = paste0("./", files_1000BCE),
                      pattern="_-9_V1.0_Fr",
                      replacement="_timeID_0"))

setwd("D:/FRANCE/Covariates/Climate/TimeID_-8/")
files_900BCE <- list.files(pattern="Fr.tif")

file.copy(from = paste0("./",files_900BCE),
          to = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_1/", files_900BCE))

file.rename(from = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_1/", files_900BCE),
            to = gsub(x = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_1/", files_900BCE),
                      pattern="CHELSA_TraCE21k_",
                      replacement="CHELSA_"))
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_1/")
files_900BCE <- list.files(pattern="CHELSA")
file.rename(from = paste0("./", files_900BCE),
            to = gsub(x = paste0("./", files_900BCE),
                      pattern="_-8_V1.0_Fr",
                      replacement="_timeID_1"))

setwd("D:/FRANCE/Covariates/Climate/TimeID_-7/")
files_800BCE <- list.files(pattern="Fr.tif")

file.copy(from = paste0("./",files_800BCE),
          to = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_2/", files_800BCE))

file.rename(from = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_2/", files_800BCE),
            to = gsub(x = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_2/", files_800BCE),
                      pattern="CHELSA_TraCE21k_",
                      replacement="CHELSA_"))
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_2/")
files_800BCE <- list.files(pattern="CHELSA")
file.rename(from = paste0("./", files_800BCE),
            to = gsub(x = paste0("./", files_800BCE),
                      pattern="_-7_V1.0_Fr",
                      replacement="_timeID_2"))

setwd("D:/FRANCE/Covariates/Climate/TimeID_-6/")
files_700BCE <- list.files(pattern="Fr.tif")

file.copy(from = paste0("./",files_700BCE),
          to = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_3/", files_700BCE))

file.rename(from = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_3/", files_700BCE),
            to = gsub(x = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_3/", files_700BCE),
                      pattern="CHELSA_TraCE21k_",
                      replacement="CHELSA_"))
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_3/")
files_700BCE <- list.files(pattern="CHELSA")
file.rename(from = paste0("./", files_700BCE),
            to = gsub(x = paste0("./", files_700BCE),
                      pattern="_-6_V1.0_Fr",
                      replacement="_timeID_3"))

setwd("D:/FRANCE/Covariates/Climate/TimeID_-5/")
files_600BCE <- list.files(pattern="Fr.tif")

file.copy(from = paste0("./",files_600BCE),
          to = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_4/", files_600BCE))

file.rename(from = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_4/", files_600BCE),
            to = gsub(x = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_4/", files_600BCE),
                      pattern="CHELSA_TraCE21k_",
                      replacement="CHELSA_"))
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_4/")
files_600BCE <- list.files(pattern="CHELSA")
file.rename(from = paste0("./", files_600BCE),
            to = gsub(x = paste0("./", files_600BCE),
                      pattern="_-5_V1.0_Fr",
                      replacement="_timeID_4"))

setwd("D:/FRANCE/Covariates/Climate/TimeID_-4/")
files_500BCE <- list.files(pattern="Fr.tif")

file.copy(from = paste0("./",files_500BCE),
          to = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_5/", files_500BCE))

file.rename(from = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_5/", files_500BCE),
            to = gsub(x = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_5/", files_500BCE),
                      pattern="CHELSA_TraCE21k_",
                      replacement="CHELSA_"))
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_5/")
files_500BCE <- list.files(pattern="CHELSA")
file.rename(from = paste0("./", files_500BCE),
            to = gsub(x = paste0("./", files_500BCE),
                      pattern="_-4_V1.0_Fr",
                      replacement="_timeID_5"))

setwd("D:/FRANCE/Covariates/Climate/TimeID_-3/")
files_400BCE <- list.files(pattern="Fr.tif")

file.copy(from = paste0("./",files_400BCE),
          to = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_6/", files_400BCE))

file.rename(from = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_6/", files_400BCE),
            to = gsub(x = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_6/", files_400BCE),
                      pattern="CHELSA_TraCE21k_",
                      replacement="CHELSA_"))
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_6/")
files_400BCE <- list.files(pattern="CHELSA")
file.rename(from = paste0("./", files_400BCE),
            to = gsub(x = paste0("./", files_400BCE),
                      pattern="_-3_V1.0_Fr",
                      replacement="_timeID_6"))


setwd("D:/FRANCE/Covariates/Climate/TimeID_-2/")
files_300BCE <- list.files(pattern="Fr.tif")

file.copy(from = paste0("./",files_300BCE),
          to = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_7/", files_300BCE))

file.rename(from = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_7/", files_300BCE),
            to = gsub(x = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_7/", files_300BCE),
                      pattern="CHELSA_TraCE21k_",
                      replacement="CHELSA_"))
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_7/")
files_300BCE <- list.files(pattern="CHELSA")
file.rename(from = paste0("./", files_300BCE),
            to = gsub(x = paste0("./", files_300BCE),
                      pattern="_-2_V1.0_Fr",
                      replacement="_timeID_7"))


setwd("D:/FRANCE/Covariates/Climate/TimeID_-1/")
files_200BCE <- list.files(pattern="Fr.tif")

file.copy(from = paste0("./",files_200BCE),
          to = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_8/", files_200BCE))

file.rename(from = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_8/", files_200BCE),
            to = gsub(x = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_8/", files_200BCE),
                      pattern="CHELSA_TraCE21k_",
                      replacement="CHELSA_"))
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_8/")
files_200BCE <- list.files(pattern="CHELSA")
file.rename(from = paste0("./", files_200BCE),
            to = gsub(x = paste0("./", files_200BCE),
                      pattern="_-1_V1.0_Fr",
                      replacement="_timeID_8"))

setwd("D:/FRANCE/Covariates/Climate/TimeID_0/")
files_100BCE <- list.files(pattern="Fr.tif")

file.copy(from = paste0("./",files_100BCE),
          to = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_9/", files_100BCE))

file.rename(from = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_9/", files_100BCE),
            to = gsub(x = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_9/", files_100BCE),
                      pattern="CHELSA_TraCE21k_",
                      replacement="CHELSA_"))
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_9/")
files_100BCE <- list.files(pattern="CHELSA")
file.rename(from = paste0("./", files_100BCE),
            to = gsub(x = paste0("./", files_100BCE),
                      pattern="_0_V1.0_Fr",
                      replacement="_timeID_9"))


setwd("D:/FRANCE/Covariates/Climate/TimeID_1/")
files_0CE <- list.files(pattern="Fr.tif")

file.copy(from = paste0("./",files_0CE),
          to = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_10/", files_0CE))

file.rename(from = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_10/", files_0CE),
            to = gsub(x = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_10/", files_0CE),
                      pattern="CHELSA_TraCE21k_",
                      replacement="CHELSA_"))
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_10/")
files_0CE <- list.files(pattern="CHELSA")
file.rename(from = paste0("./", files_0CE),
            to = gsub(x = paste0("./", files_0CE),
                      pattern="_1_V1.0_Fr",
                      replacement="_timeID_10"))


#  2. Forest cover (Zanon et al., 2018)  - 1000 BCE - 0CE ----------------------------------------------------------------

### The fossil dataset was divided into 49 timeslices ranging from 12,000 to 0 
### calibrated years BP (years before AD 1950; hereafter BP). 
### Each timeslice covers a 250-year window with the exception 
### of the most recent one, which is asymmetric as it cannot project 
### into the future, and therefore covers an interval of 185 years 
### centered on 0 BP.

### For my new Time IDs between TimeID_0 to TimeID_10 I will need the files named
### 1750 BP, 2000 BP, 2250 BP, 2500 BP, 2750 BP and 3000 BP that correspond to 
### 200 CE,  50 BCE, 300 BCE, 550 BCE,  800 BCE and 1050 BCE

setwd("D:/FRANCE/Covariates/vegetation/Zanon/")

zanon_files <- list.files(pattern=".tif")

### some Zanon maps correspond with my TimeIDs (midpoint, not whole interval)

file.copy(from = "D:/FRANCE/Covariates/vegetation/Zanon/forest_cover_2500.tif",
          to = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_4/forest_cover_2500.tif")
file.rename(from = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_4/forest_cover_2500.tif",
            to = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_4/forest_cover_timeID_4.tif")

file.copy(from = "D:/FRANCE/Covariates/vegetation/Zanon/forest_cover_2000.tif",
          to = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_9/forest_cover_2000.tif")
file.rename(from = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_9/forest_cover_2000.tif",
            to = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_9/forest_cover_timeID_9.tif")


### The time IDs correspond to the following interval midpoint

yearsBP <- c(-3000,-2900,-2800,-2700,-2600,-2500,-2400,-2300,-2200,-2100,-2000,-1900)
timeSteps <- (paste0("TimeID_",-1:10))

### function to make linear regression model

lmZanon.1000BCE_0CE.func <- function(Zanon0, Zanon1, timeZ0, timeZ1, TimeBP,outdir, name) {
  
  ### Make a SpatRaster
  input <- c(Zanon0,Zanon1)

  ### Calculate slope and intercept
  slope  <-  terra::lapp(input, fun=function(x,y){return((y-x)/(timeZ1-timeZ0))})
  input <- c(slope,Zanon1)
  intercept <- terra::lapp(input, fun=function(x,y){return(y - (x*timeZ1))})
  
  ### Estimate forest cover for target time
  input <- c(intercept,slope)
  forest_cover <- terra::lapp(input, fun=function(x,y){return(x+(y*TimeBP))},
                              filename=paste0(outdir, name,".tif"), 
                              overwrite=TRUE)
  return(forest_cover)
  
  }

### TimeID_0
setwd("D:/FRANCE/Covariates/vegetation/Zanon/")

Zanon_3000BP <- rast("forest_cover_3000.tif")
Zanon_2750BP <- rast("forest_cover_2750.tif")
Zanon_2500BP <- rast("forest_cover_2500.tif")
Zanon_2250BP <- rast("forest_cover_2250.tif")
Zanon_2000BP <- rast("forest_cover_2000.tif")
Zanon_1750BP <- rast("forest_cover_1750.tif")

putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_0/"
forest_cover_timeID_0 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_3000BP, 
                                                  Zanon1 = Zanon_2750BP, 
                                                  timeZ0 = -3000, 
                                                  timeZ1 = -2750, 
                                                  TimeBP = -2900,
                                                  outdir = putdir,
                                                  name="forest_cover_timeID_0")
plot(forest_cover_timeID_0)

### timeID_1
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_1/"
forest_cover_timeID_1 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_3000BP, 
                                                  Zanon1 = Zanon_2750BP, 
                                                  timeZ0 = -3000, 
                                                  timeZ1 = -2750, 
                                                  TimeBP = -2800,
                                                  outdir = putdir,
                                                  name="forest_cover_timeID_1")
plot(forest_cover_timeID_1)

### timeID_2
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_2/"
forest_cover_timeID_2 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_2750BP, 
                                                  Zanon1 = Zanon_2500BP, 
                                                  timeZ0 = -2750, 
                                                  timeZ1 = -2500, 
                                                  TimeBP = -2700,
                                                  outdir = putdir,
                                                  name="forest_cover_timeID_2")
plot(forest_cover_timeID_2)

### timeID_3
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_3/"
forest_cover_timeID_3 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_2750BP, 
                                                  Zanon1 = Zanon_2500BP, 
                                                  timeZ0 = -2750, 
                                                  timeZ1 = -2500, 
                                                  TimeBP = -2600,
                                                  outdir = putdir,
                                                  name="forest_cover_timeID_3")
plot(forest_cover_timeID_3)

### timeID_5
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_5/"
forest_cover_timeID_5 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_2500BP, 
                                                  Zanon1 = Zanon_2250BP, 
                                                  timeZ0 = -2500, 
                                                  timeZ1 = -2250, 
                                                  TimeBP = -2400,
                                                  outdir = putdir,
                                                  name="forest_cover_timeID_5")
plot(forest_cover_timeID_5)

### timeID_6
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_6/"
forest_cover_timeID_6 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_2500BP, 
                                                  Zanon1 = Zanon_2250BP, 
                                                  timeZ0 = -2500, 
                                                  timeZ1 = -2250, 
                                                  TimeBP = -2300,
                                                  outdir = putdir,
                                                  name="forest_cover_timeID_6")
plot(forest_cover_timeID_6)

### timeID_7
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_7/"
forest_cover_timeID_7 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_2250BP, 
                                                  Zanon1 = Zanon_2000BP, 
                                                  timeZ0 = -2250, 
                                                  timeZ1 = -2000, 
                                                  TimeBP = -2200,
                                                  outdir = putdir,
                                                  name="forest_cover_timeID_7")
plot(forest_cover_timeID_7)

### timeID_8
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_8/"
forest_cover_timeID_8 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_2250BP, 
                                                  Zanon1 = Zanon_2000BP, 
                                                  timeZ0 = -2250, 
                                                  timeZ1 = -2000, 
                                                  TimeBP = -2100,
                                                  outdir = putdir,
                                                  name="forest_cover_timeID_8")
plot(forest_cover_timeID_8)

### timeID_10
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_10/"
forest_cover_timeID_10 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_2000BP, 
                                                  Zanon1 = Zanon_1750BP, 
                                                  timeZ0 = -2000, 
                                                  timeZ1 = -1750, 
                                                  TimeBP = -1900,
                                                  outdir = putdir,
                                                  name="forest_cover_timeID_10")
plot(forest_cover_timeID_10)

### end of the script
