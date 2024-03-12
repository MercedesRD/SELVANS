#### Genosoil and phenosoil mapping for France
### Date 21/02/2024
### Author: Mercedes Roman Dobarco
### Objective: Prepare climate variables for the period 1901-1980

### Packages
library(terra)
library(devtools)
#install_github("JoshOBrien/gdalUtilities")
#library(gdalUtilities)
library(sf)

### Input data 
### Extent, CRS and resolution from the DEM file
### CRS: EPSG:4326
### Extent: -5.52182637678798, 8.685779331278, 41.3037378828989, 51.0902184585588 (xmin, xmax, ymin, ymax)
### Resolution: 0.001003433, 0.001003433 
# dem <- rast("D:/FRANCE/Covariates/scaled/srtm.tif")
# textent <- ext(dem)
# ### for gdalwarp
# t_ext <- c(textent[1], textent[3], textent[2], textent[4])
# t_res <- res(dem)
# rm(textent)
# template.r <- dem
# values(template.r) <- NA

### These settings were too memory demanding. I prefer to keep the original resolution

### Rasterized france_buffer_WGS84 to create mask, created in the script 1.Paleoclimate.R
FrBf_WGS84 <- rast("D:/FRANCE/Covariates/Administrative/FrBf_WGS84Mask.tif")
plot(FrBf_WGS84)

# ### Data from CHELSAcruts (1901-1980) -----------------------------------

### Create folders for the different decades: 1901-1980
setwd("D:/FRANCE/Covariates/Climate/")
dir.create("CHELSAcruts")
climDir <- "D:/FRANCE/Covariates/Climate/CHELSAcruts/"
setwd(climDir)

### My time intervals are:
# TimeID_20: 1901 - 1909
# TimeID_21: 1910 - 1919
# TimeID_22: 1920 - 1929
# TimeID_23: 1930 - 1939
# TimeID_24: 1940 - 1949
# TimeID_25: 1950 - 1959
# TimeID_26: 1960 - 1969
# TimeID_27: 1970 - 1979

timeSteps <- c(20:27) 
for(i in 1:length(timeSteps)){
  dir.create(paste0("TimeID_",timeSteps[[i]]))
}

timeIDs <- paste0("TimeID_",c(20:27))
vars <- c("prec", "tmax", "tmin")
years <- c(1901:1979)
months <- c(1:12)

### function modified from 
# https://rdrr.io/github/kapitzas/WorldClimTiles/

### this function downloads, resamples and crops (gdalwarp) to the extent of France
  
get_chelsa_cruts_france <- function(climdir, years, months, vars, timeID){
  timeout_old <- getOption('timeout')
  options(timeout=1000000)
  for(var in 1:length(vars)){
    for(year in 1:length(years)){
      for(month in 1:length(months)){
        
        ### Name of the file to download
        name <- paste0("CHELSAcruts_",vars[[var]],"_",months[[month]],"_",years[[year]], "_V.1.0.tif")
        source_url <- file.path(paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_cruts/",vars[[var]],"/",name))
        destination <- file.path(paste0(climdir, timeID,"/",name))  
        
        if(!file.exists(destination)){
          download.file(source_url, destination, method="wget")
        } else { 
          message(paste0(destination, " already downloaded, skipping to next"))
        }
        
        ### Extent for cropping
        FrBf_WGS84 <- rast("D:/FRANCE/Covariates/Administrative/FrBf_WGS84Mask.tif")
        
        ### Crop to the extent of France using the dem as template
        # gdalwarp(
        #   srcfile = destination, ### Destination from the previous download. It may be confusing 
        #   dstfile = gsub(x=destination, pattern= ".tif", replacement = "_Fr.tif"), ### Add a FR at the end
        #   t_srs = "EPSG:4326", ### target CRS
        #   te = t_ext, ### Target extent, from DEM
        #   tr = t_res, ### Same with resolution
        #   r = "bilinear", ## bilinear resampling because it is a continuous variable
        #   dryrun = TRUE
        #   )
        
        ### don´t know why gdalwarp does not work...
        ### I do it with terra instead
        sraster <- terra::rast(destination) ### Same CRS so no need to project, just resample
        
        ### I believe the NA is given by -32678
        NAflag(sraster) <- -32768
        f <- function(x) ifelse(x == -32768, NA, x)
        sraster <- app(sraster,f)
    
        ### I take the DEM as template
        # traster <- terra::resample(sraster, 
        #                            template, 
        #                            method="bilinear",
        #                            filename= gsub(x=destination, pattern= ".tif", replacement = "_Fr.tif"),
        #                            overwrite=TRUE
        #                            )
        
        ### crop
        traster <- terra::crop(sraster, 
                               FrBf_WGS84,
                               filename= gsub(x=destination, pattern= ".tif", replacement = "_Fr.tif"),
                               overwrite=TRUE)
        
        ### Delete the large file
        file.remove(destination)
        
      }
      options(timeout = timeout_old)
    }
  }
}

### Apply to time step 20
get_chelsa_cruts_france(climdir = "D:/FRANCE/Covariates/Climate/CHELSAcruts/",
                        timeID = "TimeID_20", 
                        years = c(1901:1909),
                        months = c(1:12),
                        vars = c("prec","tmax","tmin"))


### within each Time_ID, average for the decade
average_chelsa_cruts_france <- function(climdir, years, months, vars, timeID, timename){
  for(var in 1:length(vars)){
    for(month in 1:length(months)){
      ### List all the files for that Time_ID, variable, and month
      inputdir <- paste0(climdir,timeID,"/")
      setwd(inputdir)
      my_pattern  <- paste0("CHELSAcruts_",vars[[var]],"_",months[[month]],"_")
      input_files <- list.files(pattern=my_pattern)
      input_rasters <- rast(input_files)
      ### average for this month
      terra::app(input_rasters,
                 fun=mean, 
                 filename=paste0("CHELSAcruts_",vars[[var]],"_",months[[month]],"_",timename,"_V.1.0_Fr.tif"), 
                 overwrite=TRUE)
      ### Delete all yearly files
      file.remove(input_files)
    }
  }
}

average_chelsa_cruts_france(climdir = "D:/FRANCE/Covariates/Climate/CHELSAcruts/",
                            timeID = "TimeID_20", 
                            years = c(1901:1909),
                            months = c(1:12),
                            vars = c("prec","tmax","tmin"),
                            timename="1901_1909")

### check the average is ok, and delete rest of the files
setwd("D:/FRANCE/Covariates/Climate/CHELSAcruts/TimeID_20/")
### Create raster stack with precipitation data (12 layers)
### Load the files in order (although it should not make any difference for the bioclim variables)

vars <- c("prec","tmax","tmin")
years <- c(1901:1979)
months <- c(1:12)

prec_files  <- paste0("CHELSAcruts_prec_",c(1:12),"_1901_1909_V.1.0_Fr.tif")
tmax_files  <- paste0("CHELSAcruts_tmax_",c(1:12),"_1901_1909_V.1.0_Fr.tif")
tmin_files  <- paste0("CHELSAcruts_tmin_",c(1:12),"_1901_1909_V.1.0_Fr.tif")

prec_r <- rast(prec_files); plot(prec_r)
tmax_r <- rast(tmax_files); plot(tmax_r)
tmin_r <- rast(tmin_files); plot(tmin_r)

### change names
names(prec_r) <- paste0("prec_",c(1:12),"_1901_1909")
names(tmax_r) <- paste0("tmax_",c(1:12),"_1901_1909")
names(tmin_r) <- paste0("tmin_",c(1:12),"_1901_1909")

### Calculate bioclimatic variables
#install.packages('predicts', repos='https://rspatial.r-universe.dev')
library("predicts")
bioclim_1901_1909 <-bcvars(prec=prec_r, 
                           tmin=tmin_r, 
                           tmax=tmax_r)#,
                           #datatype="INT4S",
                          # filename=paste0("bioclim_",c(1:19),"_1901_1909.tif"), 
                          # overwrite=TRUE)
plot(bioclim_1901_1909)
bioclim_1901_1909 <- rast("bioclim_1_1901_1909.tif") ### Rename to "bioclim_1901_1909_all.tif"

writeRaster(bioclim_1901_1909,
            filename=paste0("bioclim",c(1:19),"_1901_1909.tif"), 
            overwrite=TRUE)
file.remove("bioclim_1_1901_1909.tif")
 ### Delete and keep single file with 19 layers

### Can I do this in parallel?
setwd("D:/FRANCE/Covariates/Climate/CHELSAcruts/")

### Define my parameters
timeIDs <- paste0("TimeID_",c(20:27))

year_intervals <- list(TimeID_20=c(1901:1909),
                       TimeID_21=c(1910:1919),
                       TimeID_22=c(1920:1929),
                       TimeID_23=c(1930:1939),
                       TimeID_24=c(1940:1949),
                       TimeID_25=c(1950:1959),
                       TimeID_26=c(1960:1969),
                       TimeID_27=c(1970:1979))

timenames <- list(TimeID_20="1901_1909",
                  TimeID_21="1910_1919",
                  TimeID_22="1920_1929",
                  TimeID_23="1930_1939",
                  TimeID_24="1940_1949",
                  TimeID_25="1950_1959",
                  TimeID_26="1960_1969",
                  TimeID_27="1970_1979")

library(foreach)
library(doParallel)

# tic <- Sys.time()
# detectCores()
# cl <- makeCluster(2)   ###
# registerDoParallel(cl)
# getDoParWorkers()
# 
# foreach(timeID = 2:length(timeIDs),
#         .packages=c("terra","predicts"),
#         .export = c("get_chelsa_cruts_france","average_chelsa_cruts_france",
#                     "timeIDs","year_intervals","timenames")) %dopar% {
          #timeID = 7
  
for(timeID in 6:length(timeIDs)) {
 
#for(timeID in 2:5) { 
print(timeID)
  
  timeID=8
      
          tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)
          #timeID = 1          
          ### download and crop files
          get_chelsa_cruts_france(climdir = "D:/FRANCE/Covariates/Climate/CHELSAcruts/",
                                  timeID = timeIDs[[timeID]],
                                  years = year_intervals[[timeID]],
                                  months = c(1:12),
                                  vars = c("prec","tmax","tmin"))
          tmpFiles(current=FALSE, orphan=TRUE, old=TRUE, remove=TRUE)
          
          ### within each Time_ID, average for the decade
          average_chelsa_cruts_france(climdir = "D:/FRANCE/Covariates/Climate/CHELSAcruts/",
                                      timeID = timeIDs[[timeID]], 
                                      years = year_intervals[[timeID]],
                                      months = c(1:12),
                                      vars = c("prec","tmax","tmin"),
                                      timename=timenames[[timeID]])
          tmpFiles(current=FALSE, orphan=TRUE, old=TRUE, remove=TRUE)
          
          inputdir <- paste0("D:/FRANCE/Covariates/Climate/CHELSAcruts/",timeIDs[[timeID]],"/")
          setwd(inputdir)
          
          prec_files  <- paste0("CHELSAcruts_prec_",c(1:12),"_",timenames[[timeID]],"_V.1.0_Fr.tif")
          tmax_files  <- paste0("CHELSAcruts_tmax_",c(1:12),"_",timenames[[timeID]],"_V.1.0_Fr.tif")
          tmin_files  <- paste0("CHELSAcruts_tmin_",c(1:12),"_",timenames[[timeID]],"_V.1.0_Fr.tif")
          
          prec_r <- rast(prec_files)
          tmax_r <- rast(tmax_files)
          tmin_r <- rast(tmin_files)
          
          ### change names
          names(prec_r) <- paste0("prec_",c(1:12),"_",timenames[[timeID]])
          plot(prec_r) ; prec_r
          names(tmax_r) <- paste0("tmax_",c(1:12),"_",timenames[[timeID]])
          names(tmin_r) <- paste0("tmin_",c(1:12),"_",timenames[[timeID]])
          
          ### Before calculating bioclim variables, pass the temperature to Celsius
          tmax_r <- terra::app(tmax_r,
                               fun=function(x){x/10}, 
                               filename=paste0("CHELSAcruts_tmax_",timenames[[timeID]],"_V.1.0_Fr_u.tif"), 
                               overwrite=TRUE)
          plot(tmax_r); tmax_r
          
          tmin_r <- terra::app(tmin_r,
                               fun=function(x){x/10}, 
                               filename=paste0("CHELSAcruts_tmin_",timenames[[timeID]],"_V.1.0_Fr_u.tif"), 
                               overwrite=TRUE)
          plot(tmin_r) ; tmin_r
          
          tmax_r <- rast(paste0("CHELSAcruts_tmax_",timenames[[timeID]],"_V.1.0_Fr_u.tif"))
          tmin_r <- rast(paste0("CHELSAcruts_tmin_",timenames[[timeID]],"_V.1.0_Fr_u.tif"))
         
          ### Calculate bioclimatic variables
          library("predicts")
          bioclim_t <-bcvars(prec=prec_r,
                             tmin=tmin_r, 
                             tmax=tmax_r)
          
          plot(bioclim_t); bioclim_t
          
          writeRaster(bioclim_t,
                     filename=paste0("bio",c(1:19),"_",timenames[[timeID]],".tif"), 
                     overwrite=TRUE)
          
          file.remove(tmax_files)
          file.remove(prec_files)
          file.remove(tmin_files)
          file.remove(paste0("CHELSAcruts_tmax_",timenames[[timeID]],"_V.1.0_Fr_u.tif"))
          file.remove(paste0("CHELSAcruts_tmin_",timenames[[timeID]],"_V.1.0_Fr_u.tif"))
          
          tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)
          
          gc()
         
}

# stopCluster(cl)
# 
# tac <- Sys.time()
# tac-tic

### end of this script