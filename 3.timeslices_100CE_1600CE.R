### Timeslices - preparation of files for each timeID
### range: 100 CE - 1600 CE (start time)

### Date 08/04/2024
### Author: Mercedes Roman Dobarco
### Objective: Prepare forest cover variables and copy files for bioclimate and land use.

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

### Prepare output directories

### set output wd
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices")
timeSteps <- c(11:26) 
for(i in 1:length(timeSteps)){dir.create(paste0("TimeID_",timeSteps[[i]]))}


# 1. HYDE - 100 CE - 1600 CE ---------------------------------------------------

### We have HYDE (centennial time step), e.g., grazing0AD_Fr.tif -- -- timeID_10

inputHyde <- c("h100AD_lu","h200AD_lu","h300AD_lu","h400AD_lu","h500AD_lu","h600AD_lu",
               "h700AD_lu","h800AD_lu","h900AD_lu","h1000AD_lu","h1100AD_lu","h1200AD_lu",
               "h1300AD_lu","h1400AD_lu","h1500AD_lu","h1600AD_lu")
inputnames <- c("100AD_Fr","200AD_Fr","300AD_Fr","400AD_Fr","500AD_Fr","600AD_Fr",
                "700AD_Fr","800AD_Fr","900AD_Fr","1000AD_Fr","1100AD_Fr","1200AD_Fr",
                "1300AD_Fr","1400AD_Fr","1500AD_Fr","1600AD_Fr")

for(i in 1:length(timeSteps)){
  print(i)
  setwd(paste0("D:/HYDE/HYDE3.2/unzip/",inputHyde[[i]],"/"))
  files_hyde <- list.files(pattern=".tif")
  
  file.copy(from = paste0("D:/HYDE/HYDE3.2/unzip/",inputHyde[[i]],"/",files_hyde),
            to = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_",timeSteps[[i]], "/", files_hyde))
  
  file.rename(from = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_",timeSteps[[i]], "/", files_hyde),
              to = gsub(x = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_",timeSteps[[i]], "/", files_hyde),
                        pattern=inputnames[[i]],
                        replacement=paste0("_timeID_",timeSteps[[i]])))
  }



#########################################################################################################

#  2. Climate  - 1000 BCE - 0CE ----------------------------------------------------------------

### Copy and rename files

### We have CHELSA (centennial time step), 

inputChelsa <- paste0("TimeID_",2:17)
original_timeIDs <- 2:17
my

for(i in 1:length(timeSteps)){
  print(i)
  setwd(paste0("D:/FRANCE/Covariates/Climate/",inputChelsa[[i]],"/"))
  files_chelsa <- list.files(pattern="Fr.tif")
  
  file.copy(from = paste0("./",files_chelsa),
            to = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_",timeSteps[[i]], "/", files_chelsa))
  
  file.rename(from = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_",timeSteps[[i]], "/", files_chelsa),
              to = gsub(x = paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_",timeSteps[[i]], "/", files_chelsa),
                        pattern="CHELSA_TraCE21k_",
                        replacement="CHELSA_"))
  
  setwd(paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_",timeSteps[[i]], "/"))
  files_chelsa <- list.files(pattern="CHELSA")
  file.rename(from = paste0("./", files_chelsa),
              to = gsub(x = paste0("./", files_chelsa),
                        pattern=paste0(original_timeIDs[[i]],"_V1.0_Fr"),
                        replacement=paste0("timeID_",timeSteps[[i]])))
  
  
}



#  2. Forest cover (Zanon et al., 2018)  - 100 CE - 1600 CE ----------------------------------------------------------------

### The fossil dataset was divided into 49 timeslices ranging from 12,000 to 0 
### calibrated years BP (years before AD 1950; hereafter BP). 
### Each timeslice covers a 250-year window with the exception 
### of the most recent one, which is asymmetric as it cannot project 
### into the future, and therefore covers an interval of 185 years 
### centered on 0 BP.

### For my new Time IDs between TimeID_11 to TimeID_26 I will need the files named
### 250 BP, 500 BP, 750 BP, 1000 BP, 1250 BP, 1500 BP, 1750 BP, 2000 BP that correspond to 
### 1700 CE, 1450 CE, 1200 CE, 950 CE, 700 CE, 450 CE, 200 CE, 50 BCE

setwd("D:/FRANCE/Covariates/vegetation/Zanon/")

zanon_files <- list.files(pattern=".tif")

### some Zanon maps correspond with my TimeIDs (midpoint, not whole interval)

file.copy(from = "D:/FRANCE/Covariates/vegetation/Zanon/forest_cover_1500.tif",
          to = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_14/forest_cover_1500.tif")
file.rename(from = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_14/forest_cover_1500.tif",
            to = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_14/forest_cover_timeID_14.tif")

file.copy(from = "D:/FRANCE/Covariates/vegetation/Zanon/forest_cover_1000.tif",
          to = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_19/forest_cover_1000.tif")
file.rename(from = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_19/forest_cover_1000.tif",
            to = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_19/forest_cover_timeID_19.tif")

file.copy(from = "D:/FRANCE/Covariates/vegetation/Zanon/forest_cover_500.tif",
          to = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_24/forest_cover_500.tif")
file.rename(from = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_24/forest_cover_500.tif",
            to = "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_24/forest_cover_timeID_24.tif")


### The time IDs correspond to the following interval midpoint

yearsBP <- seq(from=-1900, to=-300,by=100)
timeSteps <- (paste0("TimeID_",11:26))

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

### TimeID_11

setwd("D:/FRANCE/Covariates/vegetation/Zanon/")
Zanon_3000BP <- rast("forest_cover_3000.tif")
Zanon_2750BP <- rast("forest_cover_2750.tif")
Zanon_2500BP <- rast("forest_cover_2500.tif")
Zanon_2250BP <- rast("forest_cover_2250.tif")
Zanon_2000BP <- rast("forest_cover_2000.tif")
Zanon_1750BP <- rast("forest_cover_1750.tif")
Zanon_1500BP <- rast("forest_cover_1500.tif")
Zanon_1250BP <- rast("forest_cover_1250.tif")
Zanon_1000BP <- rast("forest_cover_1000.tif")
Zanon_750BP <- rast("forest_cover_750.tif")
Zanon_500BP <- rast("forest_cover_500.tif")
Zanon_250BP <- rast("forest_cover_250.tif")

putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_11/"
forest_cover_timeID_11 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_2000BP,
                                                   Zanon1 = Zanon_1750BP, 
                                                   timeZ0 = -2000, 
                                                   timeZ1 = -1750, 
                                                   TimeBP = -1800,
                                                   outdir = putdir,
                                                   name="forest_cover_timeID_11")
plot(forest_cover_timeID_11)

### timeID_12
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_12/"
forest_cover_timeID_12 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_1750BP, 
                                                   Zanon1 = Zanon_1500BP,
                                                   timeZ0 = -1750, 
                                                   timeZ1 = -1500,
                                                   TimeBP = -1700,
                                                   outdir = putdir,
                                                   name="forest_cover_timeID_12")
plot(forest_cover_timeID_12)

### timeID_13
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_13/"
forest_cover_timeID_13 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_1750BP, 
                                                   Zanon1 = Zanon_1500BP,
                                                   timeZ0 = -1750, 
                                                   timeZ1 = -1500,
                                                   TimeBP = -1600,
                                                   outdir = putdir,
                                                   name="forest_cover_timeID_13")
plot(forest_cover_timeID_13)

### timeID_15
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_15/"
forest_cover_timeID_15 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_1500BP, 
                                                   Zanon1 = Zanon_1250BP,
                                                   timeZ0 = -1500, 
                                                   timeZ1 = -1250,
                                                   TimeBP = -1400,
                                                   outdir = putdir,
                                                   name="forest_cover_timeID_15")
plot(forest_cover_timeID_15)

### timeID_16
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_16/"
forest_cover_timeID_16 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_1500BP, 
                                                   Zanon1 = Zanon_1250BP,
                                                   timeZ0 = -1500, 
                                                   timeZ1 = -1250,
                                                   TimeBP = -1300,
                                                   outdir = putdir,
                                                   name="forest_cover_timeID_16")
plot(forest_cover_timeID_16)

### timeID_17
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_17/"
forest_cover_timeID_17 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_1500BP, 
                                                   Zanon1 = Zanon_1250BP,
                                                   timeZ0 = -1500, 
                                                   timeZ1 = -1250,
                                                   TimeBP = -1200,
                                                   outdir = putdir,
                                                   name="forest_cover_timeID_17")
plot(forest_cover_timeID_17)

### timeID_18
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_18/"
forest_cover_timeID_18 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_1250BP, 
                                                   Zanon1 = Zanon_1000BP,
                                                   timeZ0 = -1250, 
                                                   timeZ1 = -1000,
                                                   TimeBP = -1100,
                                                   outdir = putdir,
                                                   name="forest_cover_timeID_18")
plot(forest_cover_timeID_18)


### timeID_20
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_20/"
forest_cover_timeID_20 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_1000BP, 
                                                   Zanon1 = Zanon_750BP,
                                                   timeZ0 = -1000, 
                                                   timeZ1 = -750,
                                                   TimeBP = -900,
                                                   outdir = putdir,
                                                   name = "forest_cover_timeID_20")
plot(forest_cover_timeID_20)

### timeID_21
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_21/"
forest_cover_timeID_21 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_1000BP, 
                                                   Zanon1 = Zanon_750BP,
                                                   timeZ0 = -1000, 
                                                   timeZ1 = -750,
                                                   TimeBP = -800,
                                                   outdir = putdir,
                                                   name = "forest_cover_timeID_21")
plot(forest_cover_timeID_21)

### timeID_22
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_22/"
forest_cover_timeID_22 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_750BP, 
                                                   Zanon1 = Zanon_500BP,
                                                   timeZ0 = -750, 
                                                   timeZ1 = -500,
                                                   TimeBP = -700,
                                                   outdir = putdir,
                                                   name = "forest_cover_timeID_22")
plot(forest_cover_timeID_22)


### timeID_23
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_23/"
forest_cover_timeID_23 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_750BP, 
                                                   Zanon1 = Zanon_500BP,
                                                   timeZ0 = -750, 
                                                   timeZ1 = -500,
                                                   TimeBP = -600,
                                                   outdir = putdir,
                                                   name = "forest_cover_timeID_23")
plot(forest_cover_timeID_23)

### timeID_25
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_25/"
forest_cover_timeID_25 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_500BP, 
                                                   Zanon1 = Zanon_250BP,
                                                   timeZ0 = -500, 
                                                   timeZ1 = -250,
                                                   TimeBP = -400,
                                                   outdir = putdir,
                                                   name = "forest_cover_timeID_25")
plot(forest_cover_timeID_25)

### timeID_26
putdir <- "C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_26/"
forest_cover_timeID_26 <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_500BP, 
                                                   Zanon1 = Zanon_250BP,
                                                   timeZ0 = -500, 
                                                   timeZ1 = -250,
                                                   TimeBP = -300,
                                                   outdir = putdir,
                                                   name = "forest_cover_timeID_26")
plot(forest_cover_timeID_26)


### end of the script