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
timeSteps <- c(27:46) 
for(i in 1:length(timeSteps)){dir.create(paste0("TimeID_",timeSteps[[i]]))}


# 1. HYDE - 1700 CE - 1890 CE ---------------------------------------------------

### We have HYDE (decadal time step), e.g., grazing0AD_Fr.tif -- -- timeID_10
### We copy and rename the files

inputHyde <- paste0("h",seq(from=1700, to=1890, by =10),"AD_lu")
inputnames <- paste0(seq(from=1700, to=1890, by =10),"AD_Fr")

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

#  2. Climate  - 1700 BCE - 1890CE ----------------------------------------------------------------

### Calculate gradual change, from four 100-year intervals to decadal estimates for bioclimate variables

### We have CHELSA (decadal time steps), 

### function to make linear regression model

lmCHELSA.1700CE_1890CE.func <- function(Y0, Y1, X0, X1, TimeCE, outdir, name) {
  
  ### Make a SpatRaster
  input <- c(Y0,Y1)
  
  ### Calculate slope and intercept
  slope  <-  terra::lapp(input, fun=function(x,y){return((y-x)/(X1-X0))})
  input <- c(slope,Y1)
  intercept <- terra::lapp(input, fun=function(x,y){return(y - (x*X1))})
  
  ### Estimate forest cover for target time
  input <- c(intercept,slope)
  bioclim_var <- terra::lapp(input, fun=function(x,y){return(x+(y*TimeCE))},
                              filename=paste0(outdir, name,".tif"), 
                              overwrite=TRUE)
  return(bioclim_var)
  
}

inputdirs <- c(paste0("D:/FRANCE/Covariates/Climate/","TimeID_",17:19,"/"),
               "D:/FRANCE/Covariates/Climate/CHELSAcruts/TimeID_20/")

inputFiles_17 <- paste0("CHELSA_TraCE21k_",
                        c("bio01","bio02","bio03","bio04","bio05","bio06","bio07","bio08","bio09",
                          "bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19"),
                        "_17_V1.0_Fr")

inputFiles_18 <- paste0("CHELSA_TraCE21k_",
                        c("bio01","bio02","bio03","bio04","bio05","bio06","bio07","bio08","bio09",
                          "bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19"),
                        "_18_V1.0_Fr")

inputFiles_19 <- paste0("CHELSA_TraCE21k_",
                        c("bio01","bio02","bio03","bio04","bio05","bio06","bio07","bio08","bio09",
                          "bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19"),
                        "_19_V1.0_Fr")

inputFiles_20 <- paste0(c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9",
                          "bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19"),
                        "_1901_1909")

inputFiles <- list(inputFiles_17,inputFiles_18,inputFiles_19,inputFiles_20)


timeSteps <- c(27:46) 
inputTimeCE <- seq(from=1705, to=1895, by=10)
           # 1705  1715  1725  1735  1745  1755  1765  1775  1785  1795  1805  1815  1825  1835  1845  1855  1865  1875  1885  1895
inputX0 <- c(1650, 1650, 1650, 1650, 1650, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1850, 1850, 1850, 1850, 1850)
inputX1 <- c(1750, 1750, 1750, 1750, 1750, 1850, 1850, 1850, 1850, 1850, 1850, 1850, 1850, 1850, 1850, 1905, 1905, 1905, 1905, 1905)
inputY1 <- c(rep(2,5),rep(3,10),rep(4,5))
inputY0 <- inputY1-1


namesout <- paste0("CHELSA_",
                   c("bio01","bio02","bio03","bio04","bio05","bio06","bio07","bio08","bio09",
                     "bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19"),
                   "_timeID_")


for(i in 1: length(timeSteps)){
  
  ### working in one time step
  
  for(j in 1:length(namesout)) {
    
    ### processing one bioclimate variable
    
    ### Load input raster files
    d0 <- inputY0[[i]]
    d1 <- inputY1[[i]]
    Y0raster <- terra::rast(paste0(inputdirs[[d0]],inputFiles[[d0]][[j]],".tif"))
    Y1raster <- rast(paste0(inputdirs[[d1]],inputFiles[[d1]][[j]],".tif"))
    
    putdir <- paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_", timeSteps[[i]], "/")
    bioclim_ij <- lmCHELSA.1700CE_1890CE.func(Y0 = Y0raster,
                                              Y1 = Y1raster, 
                                              X0 = inputX0[[i]],
                                              X1 = inputX1[[i]],
                                              TimeCE = inputTimeCE[[i]],
                                              outdir = putdir,
                                              name=paste0(namesout[[j]],timeSteps[[i]]))
    plot(bioclim_ij)
    print(bioclim_ij)
    
  }
  
}


#  2. Forest cover (Zanon et al., 2018)  - 1700 CE - 1900 CE ----------------------------------------------------------------

### The fossil dataset was divided into 49 timeslices ranging from 12,000 to 0 
### calibrated years BP (years before AD 1950; hereafter BP). 
### Each timeslice covers a 250-year window with the exception 
### of the most recent one, which is asymmetric as it cannot project 
### into the future, and therefore covers an interval of 185 years 
### centered on 0 BP.

### For my new Time IDs between TimeID_27 to TimeID_46 I will need the files named
### 250 BP, 0 BP that correspond to 
### 1700 CE, 1950 CE

setwd("D:/FRANCE/Covariates/vegetation/Zanon/")

zanon_files <- list.files(pattern=".tif")

### The time IDs correspond to the following interval midpoint

yearsCE <- seq(from=1705, to=1895, by=10)
timeSteps <- 27:46

### function to make linear regression model

lmZanon.1700CE_1900CE.func <- function(Zanon0, Zanon1, timeZ0, timeZ1, TimeCE, outdir, name) {
  
  ### Make a SpatRaster
  input <- c(Zanon0, Zanon1)

  ### Calculate slope and intercept
  slope  <-  terra::lapp(input, fun=function(x,y){return((y-x)/(timeZ1-timeZ0))})
  input <- c(slope,Zanon1)
  intercept <- terra::lapp(input, fun=function(x,y){return(y-(x*timeZ1))})
  
  ### Estimate forest cover for target time
  input <- c(intercept,slope)
  forest_cover <- terra::lapp(input, fun=function(x,y){return(x+(y*TimeCE))},
                              filename=paste0(outdir, name,".tif"), 
                              overwrite=TRUE)
  return(forest_cover)
  
  }

### I can do it in a loop since they all use the same input

setwd("D:/FRANCE/Covariates/vegetation/Zanon/")
Zanon_0BP <- rast("forest_cover_0.tif")
Zanon_250BP <- rast("forest_cover_250.tif")


timeSteps <- 27:46

for(i in 1: length(timeSteps)){
putdir <- paste0("C:/Users/mercedes.roman/Desktop/SELVANS/France/Output/mytimeslices/TimeID_", timeSteps[[i]], "/")
forest_cover_timei <- lmZanon.1000BCE_0CE.func(Zanon0 = Zanon_250BP,
                                               Zanon1 = Zanon_0BP,
                                               timeZ0 = 1700, 
                                               timeZ1 = 1950, 
                                               TimeCE = yearsCE[[i]],
                                                   outdir = putdir,
                                                   name=paste0("forest_cover_timeID_",timeSteps[[i]]))
plot(forest_cover_timei)
}


### end of the script
