#### Genosoil and phenosoil mapping for France
### Date 21/02/2024
### Author: Mercedes Roman Dobarco
### Objective: Prepare climate variables for the period 1980-2018

### Packages
library(terra)
library(devtools)
#install_github("JoshOBrien/gdalUtilities")
#library(gdalUtilities)
library(sf)
tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)

### Rasterized france_buffer_WGS84 to create mask, created in the script 1.Paleoclimate.R
FrBf_WGS84 <- rast("D:/FRANCE/Covariates/Administrative/FrBf_WGS84Mask.tif")
plot(FrBf_WGS84)


### Information on CHELSA climatologies

### All variables are saved as integers with a given offset and scale embedded in the geotiff file to arrive at e.g. Celsius or kg m^-2 (mm) (only for climatologies).  

### Variables: pr_01, …, pr_12 
### Defnition: Monthly precipitation amount kg m-2 month-1 
### Scale: 0.1 
### Offset: 0 
### Explanation: Precipitation amount for each month; "Amount" means mass per unit area. "Precipitation" in the Earth's atmosphere means precipitation of water in all phases. 




# ### Data from CHELSA timeseries (1980-2018) -----------------------------------

### Create folders for the different decades: 1901-1980
setwd("D:/FRANCE/Covariates/Climate/")
dir.create("CHELSApresent")
climDir <- "D:/FRANCE/Covariates/Climate/CHELSApresent/"
setwd(climDir)

### My time intervals are:
# TimeID_28: 1980 - 1989 (1985)
# TimeID_29: 1985 - 1995 (1990)
# TimeID_30: 1995	- 2005 (2000)
# TimeID_31: 2001 - 2011 (2006)
# TimeID_32: 2007 - 2017 (2012)
# TimeID_33: 2018 - 2018 (2018)


timeSteps <- c(28:33) 
for(i in 1:length(timeSteps)){
  dir.create(paste0("TimeID_",timeSteps[[i]]))
}

timeIDs <- paste0("TimeID_",c(28:33))
vars <- c("pr", "tasmax", "tasmin")
starts <- c(1980, 1985, 1995, 2001, 2007, 2018)
ends <- c(1989, 1995, 2005, 2011, 2017, 2018)
months <- c("01","02","03","04","05","06","07","08","09","10","11","12")

### function modified from 
# https://rdrr.io/github/kapitzas/WorldClimTiles/

# https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/monthly/pr/CHELSA_pr_01_1979_V.2.1.tif
# https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/monthly/pr/CHELSA_pr_01_1980_V.2.1.tif

### this function downloads, resamples and crops (gdalwarp) to the extent of France
  
get_chelsa_timeseries_france <- function(climdir, mystart, myend, thistimeID, months, vars){
  timeout_old <- getOption('timeout')
  options(timeout=1000000)
  for(var in 1:length(vars)){
    
    for(month in 1:length(months)){
    
      target.years <- mystart: myend
      
      for(year in 1:length(target.years)){
        
        ### Name of the file to download
        name <- paste0("CHELSA_",vars[[var]],"_",months[[month]],"_",target.years[[year]], "_V.2.1.tif")
        source_url <- file.path(paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/monthly/",vars[[var]],"/",name))
        destination <- file.path(paste0(climdir, thistimeID,"/",name))  
        
        if(!file.exists(destination)){
          download.file(source_url, destination, method="wget")
        } else { 
          message(paste0(destination, " already downloaded, skipping to next"))
        }
        
        ### Extent for cropping
        FrBf_WGS84 <- rast("D:/FRANCE/Covariates/Administrative/FrBf_WGS84Mask.tif")
        
        ### don´t know why gdalwarp does not work...
        ### I do it with terra instead
        sraster <- terra::rast(destination) ### Same CRS so no need to project, just resample
        
        ### I believe the NA is given by -32678
        NAflag(sraster) <- -32768
        f <- function(x) ifelse(x == -32768, NA, x)
        sraster <- app(sraster,f)

        ### crop
        traster <- terra::crop(sraster, 
                               FrBf_WGS84,
                               filename= gsub(x=destination, pattern= ".tif", replacement = "_Fr.tif"),
                               overwrite=TRUE)
        
        ### Delete the large file
        file.remove(destination)
      
        options(timeout = timeout_old)
        
      }
    }
    }
  }



### within each Time_ID, average for the decade
average_chelsa_timeseries_france <- function(climdir, mystarts, myends, thistimeID, months, vars, timename){
  for(var in 1:length(vars)){
    for(month in 1:length(months)){
     # for(timeID in 1:length(timeIDs)){
        
      target.years <- mystarts:myends
      ### List all the files for that Time_ID, variable, and month
      inputdir <- paste0(climdir,thistimeID,"/")
      setwd(inputdir)
      #my_pattern  <- paste0("CHELSA_",vars[[var]],"_",months[[month]],"_",target.years, "_V.2.1.tif")
      input_files <- paste0("CHELSA_",vars[[var]],"_",months[[month]],"_",target.years, "_V.2.1_Fr.tif")
      input_rasters <- rast(input_files)
      ### average for this month
      raster_mean <- terra::app(input_rasters,
                                fun=mean, 
                                filename=paste0("CHELSA_",vars[[var]],"_",months[[month]],"_",timename,"_V.2.1_Fr.tif"), 
                                overwrite=TRUE)
      ### Delete all yearly files
      file.remove(input_files)
    }
  }
}

### Can I do this in parallel?
setwd("D:/FRANCE/Covariates/Climate/CHELSApresent/")

### Define my parameters
timeIDs <- paste0("TimeID_",c(28:33))
vars <- c("pr", "tasmax", "tasmin")
starts <- c(1980, 1985, 1995, 2001, 2007, 2018)
ends <- c(1989, 1995, 2005, 2011, 2017, 2018)
months <- c("01","02","03","04","05","06","07","08","09","10","11","12")

timenames <- list(TimeID_28="1980_1989",
                  TimeID_29="1985_1995",
                  TimeID_30="1995_2005 ",
                  TimeID_31="2001_2011",
                  TimeID_32="2007_2017",
                  TimeID_33="2018_2018")

library(foreach)
library(doParallel)

tic <- Sys.time()
detectCores()
cl <- makeCluster(3)   ###
registerDoParallel(cl)
getDoParWorkers()

foreach(timeID = 1:length(timeIDs),
        .packages=c("terra","predicts"),
        .export = c("get_chelsa_timeseries_france",
                    "average_chelsa_timeseries_france",
                    "timeIDs","starts","ends","vars","timenames")) %dopar% {

          
          ### download and crop files
          get_chelsa_timeseries_france(climdir = "D:/FRANCE/Covariates/Climate/CHELSApresent/",
                                                   mystart =  starts[[timeID]], 
                                                   myend = ends[[timeID]],
                                                   thistimeID = timeIDs[[timeID]],
                                                   months = months,
                                                   vars = c("pr", "tasmax", "tasmin"))
          tmpFiles(current=FALSE, orphan=TRUE, old=TRUE, remove=TRUE)
          gc()
          
          ### within each Time_ID, average for the time interval
          average_chelsa_timeseries_france(climdir = "D:/FRANCE/Covariates/Climate/CHELSApresent/",
                                           mystart =  starts[[timeID]], 
                                           myend = ends[[timeID]],
                                           thistimeID = timeIDs[[timeID]],
                                           months = months,
                                           vars = c("pr", "tasmax", "tasmin"),
                                           timename=timenames[[timeID]])
          tmpFiles(current=FALSE, orphan=TRUE, old=TRUE, remove=TRUE)
          gc()
          
          ### Transform UNITS
          ### pr are in kg m-2 month-1 /100 and needs to be divided by 100 to pass to kg m-2 month-1
          ### Tasmax and Tasmin is in K/10 so first I need to divide by 10 and to pass to celcius, subtract 273,15
          fun_K_to_C <- function(x) {(x/10)-273.15}
          fun_pr <- function(x) {x/100}
          
          inputdir <- paste0("D:/FRANCE/Covariates/Climate/CHELSAcruts/",timeIDs[[timeID]],"/")
          setwd(inputdir)
          
          prec_files  <- paste0("CHELSA_pr_",months,"_",timenames[[timeID]],"_V.1.0_Fr.tif")
          tmax_files  <- paste0("CHELSA_tasmax_",months,"_",timenames[[timeID]],"_V.1.0_Fr.tif")
          tmin_files  <- paste0("CHELSA_tasmin_",months,"_",timenames[[timeID]],"_V.1.0_Fr.tif")
          
          prec_r <- rast(prec_files); plot(prec_r)
          tmax_r <- rast(tmax_files) ; plot(tmax_r)
          tmin_r <- rast(tmin_files); plot(tmin_r)
          
          ### change names
          names(prec_r) <- paste0("prec_",months,"_",timenames[[timeID]])
          names(tmax_r) <- paste0("tmax_",months,"_",timenames[[timeID]])
          names(tmin_r) <- paste0("tmin_",months,"_",timenames[[timeID]])
          
          
          prec_r <- terra::app(prec_r,
                               fun=fun_pr, 
                               filename=paste0("CHELSA_prec_",months[[month]],"_",timename,"_V.2.1_Fr_u.tif"), 
                               overwrite=TRUE)
          plot(prec_r)
          
          tmax_r <- terra::app(tmax_r,
                               fun=fun_K_to_C, 
                               filename=paste0("CHELSA_tmax_",months[[month]],"_",timename,"_V.2.1_Fr_u.tif"), 
                               overwrite=TRUE)
          plot(tmax_r)
          
          tmin_r <- terra::app(tmin_r,
                               fun=fun_K_to_C, 
                               filename=paste0("CHELSA_tmin_",months[[month]],"_",timename,"_V.2.1_Fr_u.tif"), 
                               overwrite=TRUE)
          plot(tmin_r)
          
          ### Calculate bioclimatic variables
          library("predicts")
          bioclim_t <-bcvars(prec=prec_r,
                             tmin=tmin_r, 
                             tmax=tmax_r)
          file.remove(tmax_files)
          file.remove(prec_files)
          file.remove(tmin_files)
          
          writeRaster(bioclim_t,
                      filename=paste0("bioclim",c(1:19),"_",timenames[[timeID]],".tif"), 
                      overwrite=TRUE)
          
          tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)
          
          gc()
}

stopCluster(cl)

tac <- Sys.time()
tac-tic

### end of this script