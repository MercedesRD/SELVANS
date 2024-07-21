#############################################################################################################################################
###  Preprocessing of environmental covariates for generating pedogenon classes for the Basque Country
###
###  In this script: Align and resample raster layers

######## 1. Align, resample, crop

### Desired extent: Basque Country
### Resolution: 25m
### CRS: EPSG=25830

###  Author: Mercedes Roman Dobarco
###  Date: 11/07/2023

####### Load packages
### Spatial
#library(rgdal)
#library(sp)
library(sf)
#library(gstat)
#library(raster)
library(terra)
library(gdalUtilities)
library(XML)
library(reproducible)

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


# ### 1. Boundaries Basque Country and RELIEF covariates --------------------------------------------

setwd("C:/Covariates/Euskadi/")
### Load a buffer of 125 m from the boundaries of the Basque Country
### This layer was derived in QGIS, from the shapefile with the administrative boundaries of
### the CAPV ("U:/Covariates/Administrative/CB_CAPV_5000_ETRS89/U11.CB_CAPV_5000_ETRS89.shp")
eus_buffer <- rast("eus_buffer_125m.tif")

### project to WGS84
eus_buffer_WGS84 <- terra::project(eus_buffer,"EPSG:4326", method="near")
plot(eus_buffer_WGS84)
writeRaster(eus_buffer_WGS84,
            "C:/Users/mercedes.roman/Desktop/SELVANS/WP1/TimeSync/Input/eus_buffer_WGS84.tif")
eus_buffer_WGS84 <- rast("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/TimeSync/Input/eus_buffer_WGS84.tif")

### what is the original extent of the DEM for the CAPV?
setwd("C:/Covariates/Euskadi/Relief/")
#setwd("U:/Covariates/Euskadi/Relief/")

### Covariates subset
covariates.selection <- c("mdt_lidar_2017_25m_etrs89.tif", # DEM
                          "slope.tif", "slope_5.tif", "slope_10.tif", # Slope calculated from 3,5,10 cells
                          "easterness.tif", # Easterness 
                          "northerness.tif", # northerness
                          "nor_slope.tif" , # northness x slope
                          "twi_saga.tif", # SAGA Topographic wetness index
                          # "p_curv.tif", "pr_curv_5.tif", "pr_curv_10.tif", # profile curvature (3,5,10)
                          # "pl_curv_5.tif", "pl_curv_10.tif", # planar curvature (5,10)
                          # "lg_curv_5.tif", "lg_curv_10.tif", # longitudinal curvature (5,10)
                          # "t_curv.tif", "cs_curv_5.tif", "cs_curv_10.tif" , # tangential curvature (3) cross-sectional curvature (5,10) (same)
                          "mrvbf.tif", "mrrtf.tif", ### MRVBF, MRRTF
                          "slope_height.tif", # Slope Height
                          "norm_height.tif", # Normalized height
                          "st_height.tif", # Standardized Height
                          "valley_depth.tif", # Valley depth
                          "mid_slope.tif",  # Mid slope
                          "tpi_8_3.tif", "tpi_20_5.tif" , # Multiscale topographic position index, calculated with a search radius of 8,20 cells, and 3,5 scales respectively
                          "geomorphon.tif", "geomorphon2.tif")  ## Geomorphons (search radius of 11, skip 1 & 3, flat 1 & 1.5)

relief.rast <- terra::rast(covariates.selection)
plot(relief.rast)
plot(relief.rast[[c(16:19)]]) 
plot(relief.rast[[c(16:17)]])  

### Check boundaries
ext(relief.rast)
# SpatExtent : 461050, 606500, 4700750, 4811750 (xmin, xmax, ymin, ymax)
crs(relief.rast)
crs("epsg:25830", describe=TRUE)

# ETRS89 / UTM zone 30N
# <25830> +proj=utm +zone=30 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 
#    CRS "+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"
###I use this as reference for the subsequent analysis



# ### 2. SOIL variables. SoilGrids for Basque Country ------------------------

library(terra)
library(gdalUtilities)
library(sf)
### Downloaded SoilGrids for Spain in another script
# setwd("U:/Covariates/Soil/SoilGrids") ### Here stored in the original projection and at 250 m resolution
# setwd("C:/Covariates/Spain/SoilGrids") ### Here for whole spain at 25 m resolution and ETRS 1989 / UTM 30N

### Take a reference raster, the DEM as many other times
C <- rast("C:/Covariates/Euskadi/relief/mdt_lidar_2017_25m_etrs89.tif")
ext(dem);res(dem)

### input Soil Grids (downloaded)
# "U:/Covariates/Spain/Soil/SoilGrids"
# U:\Covariates\Spain\Soil\SoilGrids
HomeDir <- "U:/Covariates/"
  
### I use gdalwarp
setwd("U:/Covariates/Spain/Soil/SoilGrids/")
soil.files <- list.files(pattern=".tif$")
### what is the CRS of origin?
t <- rast(paste0(HomeDir,"Spain/Soil/SoilGrids/",soil.files[[1]]))
library(reproducible)
assessDataType(t) #[1] "INT2U"

for(i in 1:length(soil.files)){
  setwd(paste0(HomeDir,"Spain/Soil/SoilGrids/"))
  print(i)
  gdalUtilities::gdalwarp(
    srcfile = paste0(HomeDir,"Spain/Soil/SoilGrids/",soil.files[[i]]),
    dstfile = gsub(".tif","_eus.tif", paste0("C:/Covariates/Euskadi/Soil/SoilGrids/",soil.files[[i]])),
    s_srs = '+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs', ### Homolosine
    t_srs = "EPSG:25830", #'+proj=utm +zone=30 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs', ### EPSG:25830
    te = c(461050,4700750, 606500, 4811750), ### <xmin ymin xmax ymax>
    tr=c(25, 25), ## 25 m - Output file resolution in target georeferenced units
    ot= "Int16",
    r="bilinear",  ## Continuous variable, but resampled to the nearest eighbour, to preserve original predictions (just downscaling)
    overwrite=TRUE)
}
gc()
rm(i, t,soil.files)


### Code for terra::project
# for(i in 1:length(soil.files)){
#   setwd(paste0(HomeDir,"Spain/Soil/SoilGrids/"))
#   print(i)
#   t <- rast(paste0(HomeDir,"Spain/Soil/SoilGrids/",soil.files[[i]]))
#   tp <- terra::project(x=t, y=dem, method="bilinear",
#                        filename=gsub(".tif","_eus.tif", paste0("C:/Covariates/Euskadi/Soil/SoilGrids/",soil.files[[i]])),
#                        overwrite=TRUE)
# }
# rm(i,soil.files, t, tp)



# ### 3. Anthromes ------------------------------------------------------------

# analysis packages
library(remotes)
remotes::install_github("nick-gauthier/anthromes")
library(anthromes)
library(stars)
library(dplyr)
library(ggplot2)
# visualization packages
library(gganimate)
library(patchwork)
library(tidyverse)
library(sf)
##devtools::load_all()

### mys study area is the Basque Country
### I bring shapefile
eus_buffer <- read_sf('C:/Covariates/Euskadi/boundaries/eus_buffer_125m.shp')

### Project to EPSG:4326 for Anthromes
demWGS84 <- terra::project(dem, "epsg:4326")
eus_bufferWGS84 <- st_transform(eus_buffer, 4326)

# Import the precomputed Anthromes-12k-DGG dataset. These include fixed inputs like land area, region, and biome as well as the anthrome classifications for HYDE 3.2 baseline, upper, and lower scenarios.

# HYDE/Anthromes fixed inputs and baseline scenario
an12_dgg_inputs <- read_sf('U:/Covariates/World/Anthromes/Anthromes-12k-DGG/an12_dgg_inputs.shp')
an12_dgg_baseline <- read_csv('U:/Covariates/World/Anthromes/Anthromes-12k-DGG/an12_dgg_baseline.csv')

## Merge
an12_dgg <- an12_dgg_inputs %>%
   left_join(an12_dgg_baseline, by = 'id')

### remove unneeded files to save RAM
rm(an12_dgg_baseline,  an12_dgg_inputs)

### subset for Basque Country
bounding_box <- st_bbox(demWGS84) %>% st_as_sfc()
an12_dgg_Eus <- st_intersection(an12_dgg, bounding_box)
an12_dgg_Eus_buffer <- st_intersection(an12_dgg, eus_bufferWGS84)

### Rasterise, taking as reference the DEM 25 m
colnames(an12_dgg_Eus)
rm(an12_dgg); gc()

## Project to ETRS89 / UTM zone 30N (EPSG:25830) 
an12_dgg_Eus_25830 <- st_transform(an12_dgg_Eus, 25830)
an12_dgg_Eus_buffer_25830 <- st_transform(an12_dgg_Eus_buffer, 25830)
plot(an12_dgg_Eus_buffer_25830["X1000AD"])

### from sf to SpatVector
#an12_dgg_Eus_25830 <- as(an12_dgg_Eus_25830, "Spatial")
an12_dgg_Eus_25830_Vect <- vect(an12_dgg_Eus_25830)
an12_dgg_Eus_buffer_25830 <- vect(an12_dgg_Eus_buffer_25830)
plot(an12_dgg_Eus_25830_Vect["1000AD"])

### what are the unique values present for anthromes?
names(an12_dgg_Eus_25830_Vect)
unique_anthromes <- apply(X = as.data.frame(an12_dgg_Eus_25830_Vect[,6:80]), MARGIN = 2, FUN=function(x) unique(x))
sort(unique(unlist(unique_anthromes)))
#  11 12 22 23 24 31 32 33 34 41 42 43 51 52 53 54
# 11 Urban
# 12 Mixed settlements
# 22 Irrigated villages
# 23 Rainfed villages
# 24 Pastoral villages
# 31 Residential irrigated croplands
# 32 Residential rainfed croplands
# 33 Populated croplands
# 34 Remote croplands
# 41 Residential rangelands
# 42 Populated rangelands
# 43 Remote rangelands
# 51 Residential woodlands
# 52 Populated woodlands
# 53 Remote woodlands
# 54 Inhabited drylands


OutDir <- "C:/Covariates/Euskadi/Organisms/anthromes/"
setwd(OutDir)
years <- names(an12_dgg_Eus_25830_Vect)[6:80]

### Seems parallel cannot work with many terra functions. so I do the rasterize in a foor loop.
# detectCores()
# cl <- makeCluster(14)   ### Create cluster
# registerDoParallel(cl)
# getDoParWorkers()

anthromes_rasts <- foreach (i=1:length(years)) %do% {

  myyear <- gsub(pattern = "X", replacement = "", x = years[[i]])
  myfilename <- paste0(OutDir,"anthrome_",myyear,".tif")
  
  rast <- terra::rasterize(x = an12_dgg_Eus_25830_Vect, 
                           y = dem, 
                           field =years[[i]], 
                           fun = "min", 
                           touches=TRUE,
                           overwrite=TRUE,
                           filename= myfilename)
  rast # We return this
  
}

anthromes_rasts_stack <- rast(anthromes_rasts)
plot(anthromes_rasts_stack)
#stopCluster(cl)

dem1km <- terra::aggregate(dem, fact=40, fun="mean")
years <- names(an12_dgg_Eus_buffer_25830)[6:80]
anthromes_rasts <- foreach (i=1:length(years)) %do% {
  
  myyear <- gsub(pattern = "X", replacement = "", x = years[[i]])
  #myfilename <- paste0(OutDir,"anthrome_",myyear,".tif")
  
  rast <- terra::rasterize(x = an12_dgg_Eus_buffer_25830, 
                           y = dem1km, 
                           field =years[[i]], 
                           fun = "min", 
                           touches=TRUE)
  # overwrite=TRUE, 
  # filename= myfilename)
  rast # We return this
  
}

anthromes_rasts_stack <- rast(anthromes_rasts)
plot(anthromes_rasts_stack[["X10000BC"]], col=q16 )
plot(anthromes_rasts_stack[["X2000AD"]], col=q16 )
plot(anthromes_rasts_stack[["X0AD"]], col=q16 )

### Transform to dataframe
anthromes_rasts_df <- as.data.frame(anthromes_rasts_stack)
my_cols <- colnames(anthromes_rasts_df)
colnames(anthromes_rasts_df) <- gsub(pattern = "X", replacement = "", x = my_cols)
my_cols <- colnames(anthromes_rasts_df)

library(dplyr)
#Only complete cases
anthromes_rasts_df <- anthromes_rasts_df[complete.cases(anthromes_rasts_df),]

## long format
anthromes_rasts_df_long <- anthromes_rasts_df %>% 
  pivot_longer(., names_to = "year", cols= my_cols, values_to="Anthromes")

# Transform into factor
anthromes_rasts_df_long$Anthromes <- factor(anthromes_rasts_df_long$Anthromes, levels=c("11", "12", "22", "23", "24",
                                                       "31", "32", "33", "34",
                                                       "41", "42", "43",
                                                       "51", "52", "53", "54"))

### Summarise, count number of cells by anthrome and year
anthromes_summary <- anthromes_rasts_df_long %>%
  group_by(., year, Anthromes) %>%
 summarise(Area = n())
anthromes_summary <- anthromes_summary %>% 
  mutate(time_step = ordered(year, levels = my_cols)) ### New variable, ordered levels of time step
area_total <- sum(anthromes_summary[anthromes_summary$year=="0AD",]$Area)
anthromes_summary$Area_perc <- round(anthromes_summary$Area/area_total*100, digits=1)


length(unique(anthromes_summary$Anthromes)) ### 16 colors I need
library(colorspace)
q16 <- qualitative_hcl(16, palette = "Dark 3")

q16 <-c("11"="brown","12"="brown1",
        "22"="cornflowerblue","23"="darkorchid","24"="#F161AE", 
        "31"="#72F6D2","32"="#CBC937","33"="#F2F081","34"="#F0EFA9",
        "41"="#E09F2A","42"="#ECB776","43"="#F6DCB5",
        "51"="#19AA34","52"="#5FCE2C","53"="#B6E284","54"="#F6F8E3")

# 11 Urban "brown"
# 12 Mixed settlements "brown1"
# 22 Irrigated villages "cornflowerblue"
# 23 Rainfed villages "darkorchid"
# 24 Pastoral villages "#F161AE"
# 31 Residential irrigated croplands "#72F6D2"
# 32 Residential rainfed croplands "#CBC937"
# 33 Populated croplands "#F2F081"
# 34 Remote croplands "#F0EFA9"
# 41 Residential rangelands "#E09F2A"
# 42 Populated rangelands "#ECB776"
# 43 Remote rangelands "#F6DCB5"
# 51 Residential woodlands "#19AA34"
# 52 Populated woodlands "#5FCE2C"
# 53 Remote woodlands "#B6E284"
# 54 Inhabited drylands "#F6F8E3"

anthromes_plot <- anthromes_summary %>% 
  ggplot(aes(time_step, Area_perc)) +
  geom_col(aes(fill =  Anthromes, color =  Anthromes), width = .8, size = .1)+
  geom_vline(xintercept = c('0AD', '1700AD'), color = 'black') +
  geom_hline(yintercept = .5, linetype = 2, color = 'black') +
  scale_fill_manual(values = q16,
                    labels= c("Urban"," Mixed settlements","Irrigated villages","Rainfed villages",
                              "Pastoral villages", 
                              "Residential irrigated croplands","Residential rainfed croplands",
                              "Populated croplands", "Remote croplands",
                              "Residential rangelands","Populated rangelands","Remote rangelands",
                              "Residential woodlands","Populated woodlands","Remote woodlands",
                              "Inhabited drylands"))+
  scale_color_manual(values = q16,
                     labels = c("Urban"," Mixed settlements","Irrigated villages","Rainfed villages",
                                "Pastoral villages",
                                "Residential irrigated croplands","Residential rainfed croplands",
                                "Populated croplands", "Remote croplands",
                                "Residential rangelands","Populated rangelands","Remote rangelands",
                                "Residential woodlands","Populated woodlands","Remote woodlands",
                                "Inhabited drylands")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=1))+
 # scale_y_continuous(expand = c(0.01,0.01), labels = scales::percent_format(accuracy = 1), n.breaks = 3) +
  scale_x_discrete(breaks = c('8000BC','200BC','0AD','800AD','1700AD','1750AD',
                              '1800AD','1850AD','1900AD','1950AD',
                              '2000AD','2005AD','2010AD','2015AD'), 
                   labels = c('8000ʙᴄᴇ','200ʙᴄᴇ','1ᴄᴇ','800ᴄᴇ','1700ᴄᴇ','1750ᴄᴇ',
                              '1800ᴄᴇ','1850ᴄᴇ','1900ᴄᴇ','1950ᴄᴇ',
                              '2000ᴄᴇ','2005ᴄᴇ','2010ᴄᴇ','2015ᴄᴇ')) +
  labs(x = 'Year', y = 'Area (%)') +
  theme(legend.position = 'bottom',
        strip.text.x = element_text(hjust = 0), panel.spacing.x=unit(.7, "lines"),
        strip.background = element_blank(), axis.line.y = element_blank(), 
        axis.text=element_text(size=rel(1)),
        panel.border = element_rect(colour = "grey20", fill=NA, size = 1))

save.image("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/1.Anthromes_Eus.RData")



# ### 4. Climate covariates at 25m resolution -----------------------------

### Take a reference raster, the DEM as many other times
HomeDir <- "C:/Covariates/"
dem <- rast("C:/Covariates/Euskadi/relief/mdt_lidar_2017_25m_etrs89.tif")
st_bbox(dem)
ext(dem);res(dem)

### crop climate covariates to Euskadi
### I use gdalwarp
setwd("C:/Covariates/World/bio/")
climate.files <- list.files(pattern=".tif")
t <- rast("CHELSA_bio1_1981-2010_V.2.1.tif" )

### what is the CRS of origin?
#t <- rast("CHELSA_bio1_1981-2010_V.2.1.tif")
crs("epsg:4326", describe=TRUE)
crs("epsg:25830", describe=TRUE)

system.time(for(i in 1:length(climate.files)){
  setwd(paste0(HomeDir,"Euskadi/Climate/"))
  print(i)
  gdalUtilities::gdalwarp(
    srcfile = paste0(HomeDir,"/World/bio/",climate.files[[i]]),
    dstfile = gsub(".tif","_eus.tif",
                   paste0(HomeDir,"Euskadi/Climate/",climate.files[[i]])),
           #s_srs = '+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs', ### Homolosine
    s_srs = "EPSG:4326", #"+proj=longlat +datum=WGS84 +no_defs +type=crs", ### EPSG:4326
    t_srs = "EPSG:25830", #'+proj=utm +zone=30 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs', ### EPSG:25830
    te = c(461050,4700750, 606500, 4811750), ###  <xmin ymin xmax ymax>
    tr=c(25, 25), ## 25 m - Output file resolution in target georeferenced units
    r="bilinear", ## Continuous variable
    overwrite=TRUE)
  
})
rm(i, t,climate.files)

setwd(paste0(HomeDir,"Euskadi/Climate/"))
climate.files <- list.files(pattern=".tif")
climate.rast <- terra::rast(climate.files)
plot(climate.rast)



# ### 5. Suelos Euskadi ----------------------------------------------------------

HomeDir <- "C:/Covariates/"  ### Change this with your Home directory 
setwd(HomeDir)
library(sf)
edafo_eus <- st_read("C:/Covariates/Euskadi/Soil/Edafologia-Suelos y capacidad de uso/edafologia_biz_gip_ETRS89.shp")
edafo_eus <- st_set_crs(edafo_eus,25830)
plot(edafo_eus[!is.na(edafo_eus$UNIDAD),"UNIDAD"])
plot(edafo_eus["UNIDAD"])
plot(edafo_eus["SUBUNIDAD"])

unique(edafo_eus$UNIDAD)
# [1] NA   "A - Acrisol"   "B - Cambisol"  "0 - Cauce"     "E - Rendsina"  "0 - Embalse"   "G - Gleysol"  
# [8] "H - Histosol"  "I - Litosol"   "J - Fluvisol"  "L - Luvisol"   "Q - Arenosol"  "R - Regosol"   "0 - Sin suelo"
# [15] "T - Andosol"   "U - Ranker"    "Z - Solonchak"

unique(edafo_eus$SUBUNIDAD)
# [1] NA                       "Ag - Acrisol gleico"    "Ah - Acrisol humico"    "Ao - Acrisol ortico"   
# [5] "Bc - Cambisol cromico"  "Bd - Cambisol districo" "Be - Cambisol eutrico"  "Bg - Cambisol gleico"  
# [9] "Bh - Cambisol humico"   "Bk - Cambisol calcico"  "00 - Cauce"             "E - Rendsina"          
# [13] "00 - Embalse"           "Ge - Gleysol eutrico"   "Hd - Histosol districo" "I - Litosol"           
# [17] "Jc - Fluvisol cromico"  "Je - Fluvisol eutrico"  "Lc - Luvisol cromico"   "Lg - Luvisol gleico"   
# [21] "Lk - Luvisol calcico"   "Lo - Luvisol ortico"    "Q - Arenosol"           "Rc - Regosol calcico"  
# [25] "Rd - Regosol districo"  "Re - Regosol eutrico"   "00 - Sin suelo"         "Th - Andosol humico"   
# [29] "U - Ranker"             "Z - Solonchak"          "Zg - Solonchak gleico" 

unique(edafo_eus$SUB_ASOCIA)
# [1] NA                       "- - -"                  "Bd - Cambisol districo" "Bg - Cambisol gleico"  
# [5] "Bh - Cambisol humico"   "Bk - Cambisol calcico"  "Ag - Acrisol gleico"    "Ao - Acrisol ortico"   
# [9] "Be - Cambisol eutrico"  "Lc - Luvisol cromico"   "Lg - Luvisol gleico"    "Lo - Luvisol ortico"   
# [13] "Th - Andosol humico"    "Bc - Cambisol cromico"  "U - Ranker"             "E - Rendsina"          
# [17] "I - Litosol"            "Re - Regosol eutrico"   "Rd - Regosol districo"  "Gm - Gleysol mollico"  
# [21] "Je - Fluvisol eutrico"  "Ah - Acrisol humico"    "Lk - Luvisol calcico"   "Ge - Gleysol eutrico"

SMU <- unique(edafo_eus$CODIGO)
SMU <- SMU[-1]
length(SMU)

### Rewrite code SMU
edafo_eus$SMU <- NA
end <- str_locate(edafo_eus[!is.na(edafo_eus$CODIGO),]$CODIGO, pattern = "/")[,1]-1
edafo_eus[!is.na(edafo_eus$CODIGO),]$SMU <- str_sub(edafo_eus[!is.na(edafo_eus$CODIGO),]$CODIGO, start = 1, end = end)

### Unique SMU
SMU <- unique(edafo_eus$SMU)
SMU <- SMU[-1]
length(SMU)
### 147 SMU

## Area of each polygon in m2
edafo_eus$area <- st_area(edafo_eus)
Map_area <- sum(edafo_eus$area)
Map_area_km2 <- Map_area/1000000

### The soil map for Bizkaia and Guipuzkoa
A <- 4.187307 # thousand km2
Ssuborder <- 5.7984 *(A^0.2265)
Ssuborder ### 8 suborders
Sfam <- 15.2540 *(A^0.6315)
Sfam ### we would be talking about 38 soil classes at the family level
Sser <- 16.342 *(A^0.7106)
Sser ### and around 45 at the level of series
### the soil map recognises 147 different codes

### Because we don´t have the right report of the map, if we look at unique combinations of association and sub association
edafo_eus$SMU2 <- NA
edafo_eus$SMU2 <- paste(edafo_eus$SUBUNIDAD, edafo_eus$SUB_ASOCIA, sep=" , ")
unique(edafo_eus$SMU2)
edafo_eus$SMU2 <-  gsub( replacement = "", x = edafo_eus$SMU2,pattern = " , - - -")
edafo_eus$SMU2 <- ifelse(edafo_eus$SMU2  %in% c("00 - Embalse","NA , NA","00 - Cauce","00 - Sin suelo"),
                         yes= NA, no=edafo_eus$SMU2)
SMU2 <- unique(edafo_eus$SMU2)
edafo_eus$SMU2 <- factor(edafo_eus$SMU2, levels=SMU2)
### Transform to numeric
edafo_eus$SMU2num <- as.numeric(edafo_eus$SMU2)
edafo_eus[927,]
### In the table of equivalences subtract one
CODE_SMU2 <- data.frame(ID=1:117,SMU2 = SMU2[-1])

### This gives 117 SMU
### a lineal relationship would give around 200 classes for the Basque Country
### but a power relationship around 140 classes (a lot!)
### I want to rasterize this value
template_rast <- dem
values(template_rast) <- NA
setwd("C:/Covariates/Euskadi/Soil/")
edafo_eus_rast <- terra::rasterize(edafo_eus, template_rast, field="SMU2num", filename="edafo_eus.tif")
plot(edafo_eus_rast)
levels(edafo_eus_rast) <- CODE_SMU2 ### as factor


# ### 6. Potential vegetation ---------------------------------------------

potential_veg <- st_read("U:/Covariates/Euskadi/Organisms/VegPotencial/CT_VEGETACION_POTENCIAL_100000_ETRS89.shp")
#copy file to C:/
potential_veg <- st_set_crs(potential_veg,25830)
plot(potential_veg["VEGETACION"])
unique(potential_veg$VEGETACION)

area.tot <- sum(potential_veg$SHAPE_AREA, na.rm=TRUE)
potential_veg$area_rel <- round(potential_veg$SHAPE_AREA / area.tot *100, digits=2)
hist(potential_veg$area_rel, breaks=30, xlim=c(0,10))

potential_veg <- potential_veg %>% arrange(., SHAPE_AREA) 
potential_veg$VEGETACION

# [1] "Complejo de comunidades ligadas a las rocas silíceas"   "Vegetación de cubetas endorreicas"                     
# [3] "Alcornocal"                                             "Vegetación herbácea ligada al agua"                    
# [5] "Pinar de pino carrasco"                                 "Pinar de pino albar (Pinus sylvestris)"                
# [7] "Vegetación de arenales costeros"                        "Encinar del interior (carrascal estellés)"             
# [9] "Vegetación de acantilados litorales"                    "Robledal de Quercus petrea"                            
# [11] "Vegetación de marismas"                                 "Quejigal-robledal calcícola (con Quercus pubescens)"   
# [13] "Quejigal con boj"                                       "Carrascal montano con boj"                             
# [15] "Hayedo con boj"                                         "Complejo de comunidades ligadas a las rocas calcáreas" 
# [17] "Alameda-aliseda mediterránea y/o de transición"         "Quejigal atlántico (con Smilax aspera y Quercus robur)"
# [19] "Carrascal mediterráneo seco"                            "Encinar cantábrico"                                    
# [21] "Aliseda cantábrica"                                     "Carrascal montano seco"                                
# [23] "Robledal eutrofo subatlántico"                          "Quejigal submediterráneo"                              
# [25] "Marojal"                                                "Hayedo calcícola o eutrofo"                            
# [27] "Quejigal subcantábrico"                                 "Hayedo acidófilo"                                      
# [29] "Robledal acidófilo y robledal-bosque mixto atlántico"  

### I want to rasterize this variable, prior to grouping
template_rast <- dem
values(template_rast) <- NA
setwd("C:/Covariates/Euskadi/Organisms/")
potential_veg_rast <- terra::rasterize(potential_veg, template_rast, field="VEGETACION",
                                       overwrite=TRUE, filename="potential_veg.tif")
plot(potential_veg_rast)
potential_veg_rast <-terra::rast("potential_veg.tif")
### There is no need to transform into numeric before rasterizing

### Dataframe with levels
veg.levels <- levels(potential_veg_rast)[[1]]

### Marojal - Quercus pyrenaica - agrupar con robles/quejigos?
### Group by broader categories (major vegetation group)
### Some categories were grouped because they occupied a very small area, 
### and were assigned to thei major group, closest to their category. E.g., "Alcornocal" --> robledal
### "Complejo de comunidades ligadas a las rocas silíceas" --> hayedo
### "Complejo de comunidades ligadas a las rocas calcáreas" --> hayedo
### "Vegetación de cubetas endorreicas" --> "Carrascal"
### "Vegetación herbácea ligada al agua" --> "Carrascal"
### Marismas (marshes) y vegetación costera y de acantilados (coastal and coastal cliffs), are grouped together

potential_veg$MVG <- NA
potential_veg$MVG <- ifelse(potential_veg$VEGETACION %in% c("Alameda-aliseda mediterránea y/o de transición","Aliseda cantábrica"), "Aliseda",
                            ifelse(potential_veg$VEGETACION %in% c("Carrascal mediterráneo seco","Carrascal montano con boj","Carrascal montano seco"), "Carrascal",
                                   # ifelse(potential_veg$VEGETACION == "Complejo de comunidades ligadas a las rocas calcáreas", "Formaciones calcareas",
                                          #       ifelse(potential_veg$VEGETACION == "Complejo de comunidades ligadas a las rocas silíceas", "Formaciones siliceas",
                                   ifelse(potential_veg$VEGETACION %in% c("Encinar cantábrico","Encinar del interior (carrascal estellés)"), "Encinar",
                                          ifelse(potential_veg$VEGETACION %in% c("Hayedo acidófilo","Hayedo calcícola o eutrofo","Hayedo con boj",
                                                                                 "Complejo de comunidades ligadas a las rocas calcáreas",
                                                                                 "Complejo de comunidades ligadas a las rocas silíceas"), "Hayedo",
                                                 ifelse(potential_veg$VEGETACION == "Marojal", "Marojal",
                                                        ifelse(potential_veg$VEGETACION %in% c("Pinar de pino albar (Pinus sylvestris)","Pinar de pino carrasco"), "Pinar",
                                                               #ifelse(potential_veg$VEGETACION == "Pinar de pino carrasco", "Pino carrasco",
                                                               ifelse(potential_veg$VEGETACION %in% c("Quejigal-robledal calcícola (con Quercus pubescens)",
                                                                                                      "Quejigal atlántico (con Smilax aspera y Quercus robur)",
                                                                                                      "Quejigal con boj","Quejigal subcantábrico","Quejigal submediterráneo"), "Quejigal",
                                                                      ifelse(potential_veg$VEGETACION %in% 
                                                                               c("Alcornocal", "Robledal acidófilo y robledal-bosque mixto atlántico",
                                                                                 "Robledal de Quercus petrea","Robledal eutrofo subatlántico"), "Robledal",
                                                                             ifelse(potential_veg$VEGETACION %in% c("Vegetación de acantilados litorales",
                                                                                                                    "Vegetación de arenales costeros",
                                                                                                                    "Vegetación de marismas"),"Costera",
                                                                                    ifelse(potential_veg$VEGETACION == "Vegetación de cubetas endorreicas", "Carrascal",
                                                                                           # ifelse(potential_veg$VEGETACION == "Vegetación de marismas", "Marismas",
                                                                                           ifelse(potential_veg$VEGETACION == "Vegetación herbácea ligada al agua", "Carrascal",
                                                                                                  NA)))))))))))
length(unique(potential_veg$MVG)) ### Reduced from 29 to 9 classes
plot(potential_veg["MVG"])

### Rasterize new aggregated categories
mvg_rast <- terra::rasterize(potential_veg, template_rast, field="MVG",overwrite=TRUE, filename="mvg.tif")
plot(mvg_rast)
mvg_rast <- terra::rast("mvg.tif")

### Dataframe with levels
mvg.levels <- levels(mvg_rast)[[1]]

### Transform both variables to continuous with a PCA
### Create dummy variables
setwd("C:/Covariates/Euskadi/Organisms/")
veg.classes <- sort(unique(values(potential_veg_rast)))
dummy.list <- list()
for(i in 1:length(veg.classes)){
  print(i)
  dummy.pm <- app(potential_veg_rast, 
                  fun = function(x) {ifelse(is.na(x), NA, ifelse(x==veg.classes[[i]],1,0))},
                  filename=paste0("veg.dummy",veg.classes[[i]],".tif"),
                  overwrite=TRUE)
  dummy.list[[i]] <- dummy.pm
}

dummy.files <- list.files(pattern="dummy")
dummy.rast <- rast(dummy.files)
rm(dummy.files)
plot(dummy.rast)
### change names
names(dummy.rast) <- paste0("veg_",veg.classes)

### Perform PCA
### First, let's get a hint of how many components we need
Npixels <- values(dummy.rast[[1]])
Npixels <- Npixels[!is.na(Npixels)]
N.sample <- length(Npixels)*0.1 ### 10% of pixels aprox

set.seed(1946)
# Regular sampling
sampleVEG<- spatSample(dummy.rast, size = 1000000 , method="regular",replace=FALSE, as.df=TRUE, xy=TRUE, na.rm = TRUE)
### select only complete cases
sampleVEG <-sampleVEG[complete.cases(sampleVEG),]
sampleVEG <- as.data.frame(sampleVEG)
summary(sampleVEG)

### Apply scaled PCA
VEG.pca <- prcomp(sampleVEG[,3:ncol(sampleVEG)], scale=TRUE)
library(factoextra)
fviz_eig(VEG.pca)
fviz_pca_var(VEG.pca,axes = c(1,2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# Eigenvalues
eig.val <- get_eigenvalue(VEG.pca)
eig.val ### Let's retain ~ 50 % variability with 14 PCs

### Try the rasterPCA function from the RStool package
pca.pred <- predict(dummy.rast, VEG.pca, index=1:14)
plot(pca.pred)
### Write them to file 
for (i in 1:nlyr(pca.pred)){
  writeRaster(pca.pred[[i]], filename = paste0("VEG_PC",i,".tif"), overwrite=TRUE)
}

VEG.files <- list.files(pattern="MVG_")
VEG.rast <- rast(VEG.files)
plot(VEG.rast)

### now MVG PCA
### Transform to continuous varible with a PCA
### Create dummy variables
setwd("C:/Covariates/Euskadi/Organisms/")
mvg.classes <- mvg.levels$value
dummy.list <- list()
for(i in 1:length(mvg.classes)){
  print(i)
  dummy.pm <- app(mvg_rast, 
                  fun = function(x) {ifelse(is.na(x), NA, ifelse(x==mvg.classes[[i]],1,0))},
                  filename=paste0("mvg.classes",mvg.classes[[i]],".tif"),
                  overwrite=TRUE)
  dummy.list[[i]] <- dummy.pm
}

dummy.files <- list.files(pattern="mvg.classes")
dummy.rast <- rast(dummy.files)
rm(dummy.files)
plot(dummy.rast)
### change names
names(dummy.rast) <- paste0("mvg_",mvg.classes)

### Perform PCA
### First, let's get a hint of how many components we need
Npixels <- values(dummy.rast[[1]])
Npixels <- Npixels[!is.na(Npixels)]
N.sample <- length(Npixels)*0.1 ### 10% of pixels aprox

set.seed(1946)
# Regular sampling
sampleVEG<- spatSample(dummy.rast, size = 3000000,
                       method="regular",replace=FALSE, 
                       as.df=TRUE, xy=TRUE, na.rm = TRUE)
### select only complete cases
sampleVEG <-sampleVEG[complete.cases(sampleVEG),]
sampleVEG <- as.data.frame(sampleVEG)
summary(sampleVEG)

### Apply scaled PCA
VEG.pca <- prcomp(sampleVEG[,3:ncol(sampleVEG)], scale=TRUE)
library(factoextra)
fviz_eig(VEG.pca)
fviz_pca_var(VEG.pca,axes = c(1,2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# Eigenvalues
eig.val <- get_eigenvalue(VEG.pca)
eig.val ### Let's retain ~ 100 % variability with 8 PCs

### Try the rasterPCA function from the RStool package
pca.pred <- predict(dummy.rast, VEG.pca, index=1:8)
plot(pca.pred)
### Write them to file 
for (i in 1:nlyr(pca.pred)){
  writeRaster(pca.pred[[i]], filename = paste0("MVG_PC",i,".tif"), overwrite=TRUE)
}
VEG.files <- list.files(pattern="MVG_PC")
VEG.rast <- rast(VEG.files)
plot(VEG.rast)
### the first 5 PCs explain the spatial variability well enough

### Clean environment


# ### 7. Lithology and regolith depth ------------------------------------------------------------------

### I have to reclass and rasterize the lithology map
library(sf)
lithology <-  st_read("U:/Covariates/Euskadi/Parent_Material/CT_LITOLOGICO_25000_ETRS89.shp")
lithology <- st_set_crs(lithology,25830)
plot(lithology["LITOLOGIA"])

### Reclassify lithology
sort(unique(lithology$LITOLOGIA))

# [1] "00 - Embalses, ríos y ?"  --> transform into "01 - Depósitos superficiales"                                        
# [2] "01 - Depósitos superficiales"                                   
# [3] "02 - Rocas detríticas de grano grueso (Areniscas). Dominante"   
# [4] "03 - Rocas detríticas de grano medio (Limolitas). Dominante"    
# [5] "04 - Rocas detríticas de grano fino (Lutitas). Dominante"       
# [6] "08 - Detríticos alternantes"                                    
# [7] "09 - Margas descarbonatadas"                                    
# [8] "10 - Margas"                                                    
# [9] "11 - Calizas impuras y calcarenitas"                            
# [10] "12 - Calizas"                                                   
# [11] "13 - Rocas volcánicas piroclásticas"                            
# [12] "14 - Rocas volcánicas en coladas"                               
# [13] "15 - Ofitas"                                                    
# [14] "16 - Arcillas con yesos y otras sales"                          
# [15] "17 - Alternancia de margocalizas, margas calizas y calcarenitas"
# [16] "18 - Dolomías"                                                  
# [17] "19 - Pizarras"                                                  
# [18] "20 - Rocas ígneas"                                              
# [19] "22 - Granitos de grano grueso"                                  
# [20] "23 - Granodioritas"                                             
# [21] "24 - Rocas filonianas" 

lithology$LITOLOGIA <- ifelse(lithology$LITOLOGIA  == "00 - Embalses, ríos y ?",
                         yes= "01 - Depósitos superficiales",
                         no=lithology$LITOLOGIA )

### I want to rasterize this value
template_rast <- dem
values(template_rast) <- NA
setwd("C:/Covariates/Euskadi/ParentMaterial/")
lithology_rast <- terra::rasterize(lithology, template_rast, field="LITOLOGIA",
                                   overwrite=TRUE, filename="lithology.tif")
lithology_rast <- terra::rast("lithology.tif")
plot(lithology_rast)

lito.levels <- levels(lithology_rast)[[1]]

### Pass to numeric
lithology_rast.num <- as.numeric(lithology_rast, 1)


#### Reclass into level 1 Lithology
rcl_level1 <- matrix(ncol=2, nrow=20, 
                     c(1:20, 
                       1, # Superficial deposits
                       2,2,2,2, # Sedimentary detrital rocks
                       3,3,3,3, # Limestone, dolostone, marl, calcareous
                       6,6, #Volcanic and subvolcanic
                       7, # Granitoids
                       4, # Clays with gypsum and other salts
                       3,3, # Limestone, dolostone, marl, calcareous
                       5, # Metamorphic
                       6, #Volcanic and subvolcanic
                       7,7, # Granitoids
                       6)) #Volcanic and subvolcanic

# [2] "01 - Depósitos superficiales"                                      ---> 1 Superficial deposits                                  
# [3] "02 - Rocas detríticas de grano grueso (Areniscas). Dominante"      ---> 2 Sedimentary detrital rocks 
# [4] "03 - Rocas detríticas de grano medio (Limolitas). Dominante"       ---> 2 Sedimentary detrital rocks    
# [5] "04 - Rocas detríticas de grano fino (Lutitas). Dominante"          ---> 2 Sedimentary detrital rocks       
# [6] "08 - Detríticos alternantes"                                       ---> 2 Sedimentary detrital rocks                                 
# [7] "09 - Margas descarbonatadas"                                       ---> 3 Limestone, dolostone, marl, calcareous                                   
# [8] "10 - Margas"                                                       ---> 3 Limestone, dolostone, marl, calcareous  
# [9] "11 - Calizas impuras y calcarenitas"                               ---> 3 Limestone, dolostone, marl, calcareous   
# [10] "12 - Calizas"                                                     ---> 3 Limestone, dolostone, marl, calcareous  
# [11] "13 - Rocas volcánicas piroclásticas"                              ---> 6 Volcanic and subvolcanic  
# [12] "14 - Rocas volcánicas en coladas"                                 ---> 6 Volcanic and subvolcanic 
# [13] "15 - Ofitas"                                                      ---> 7 Granitoids
# [14] "16 - Arcillas con yesos y otras sales"                            ---> 4 Clays with gypsum and other salts  
# [15] "17 - Alternancia de margocalizas, margas calizas y calcarenitas"  ---> 3 Limestone, dolostone, marl, calcareous
# [16] "18 - Dolomías"                                                    ---> 3 Limestone, dolostone, marl, calcareous
# [17] "19 - Pizarras"                                                    ---> 5 Metamorphic                                                  
# [18] "20 - Rocas ígneas"                                                ---> 6 Volcanic and subvolcanic                                              
# [19] "22 - Granitos de grano grueso"                                    ---> 7 Granitoids
# [20] "23 - Granodioritas"                                               ---> 7 Granitoids
# [21] "24 - Rocas filonianas"                                            ---> 6 Volcanic and subvolcanic

lithology_lv1 <- classify(lithology_rast.num,
                          rcl_level1,
                          filename="lithology_lv1.tif",
                          overwrite=TRUE)
lithology_lv1 <- rast("lithology_lv1.tif")
plot(lithology_lv1)

#### Reclass into level 2 lithology
rcl_level2 <- matrix(ncol=2, nrow=20, 
                     c(1:20,
                       1,1,1,1,1,1,1,1,1, # Sedimentary
                       3,3,3, # Igneous
                       1,1,1, # Sedimentary
                       2,     # Metamorphic
                       3,3,3,3)) # Igneous

lithology_lv2 <- classify(lithology_rast.num,
                          rcl_level2,
                          filename="lithology_lv2.tif",
                          overwrite=TRUE)
lithology_lv2 <- rast("lithology_lv2.tif")
plot(lithology_lv2)

#### Regolith depth
library(sf)
regolith <-  st_read("U:/Covariates/Euskadi/Parent_Material/CT_ESPESOR_REGOLITO_25000_ETRS89.shp")
regolith <- st_set_crs(regolith,25830)
names(regolith)
plot(regolith["CODIGO"])
plot(regolith["ESPESOR_RE"])
sort(unique(regolith$CODIGO))

regolith$CODIGO <- ifelse(regolith$CODIGO  == "EMBALSE", yes= "3", no=regolith$CODIGO)
regolith$ESPESOR_RE <- ifelse(regolith$ESPESOR_RE  == "Embalse",
                              yes= "(3) - Espesor de 1 a 2 m", 
                              no=regolith$ESPESOR_RE)
plot(regolith["ESPESOR_RE"])
sort(unique(regolith$CODIGO))
sort(unique(regolith$ESPESOR_RE))
regolith.lvls <- c("(1) - Espesor de 0 a 0.5 m",
                   "(2) - Espesor de 0.5 a 1 m",
                   "(3) - Espesor de 1 a 2 m",
                   "(4) - Espesor de 2 a 4 m",
                   "(5) - Espesor mayor de 4 m")

### I want to rasterize this value
template_rast <- dem
values(template_rast) <- NA
setwd("C:/Covariates/Euskadi/ParentMaterial/")
regolith_rast <- terra::rasterize(regolith, template_rast, field="ESPESOR_RE",
                                  overwrite=TRUE, filename="regolith.tif")
regolith_rast <- terra::rast("regolith.tif")
plot(regolith_rast, col=c(viridis(n=3,direction = 1, option = "D"),
                          viridis(n=2,direction = 1, option = "A")))


# ### 8. Magnetic field intensity --------------------------------------------------------------

# ### Load and create mosaic
# setwd("U:/Covariates/Euskadi/Parent_Material/SIGEOF_MALLA_RASTER_Gravimetria/MFI")
# mfi.files <- list.files(pattern=glob2rx("H_VMG00*.png"))
# t <- raster::stack(mfi.files[[1]])
# plot(t)
# plotRGB(t, r=1, g=2, b=3, stretch="hist")
# 
# ### make a list of rasters
# mylist <- lapply(mfi.files, FUN=raster::stack)
# myfun = function(x) { 
#   crs(x) <- "epsg:23030"
#   return(x)}
# mylist <- lapply(mylist, FUN=myfun)
# 
# plotRGB(mylist[[7]], r=1, g=2, b=3, stretch="hist")
# plot(mylist[[7]],band=3)
# 
# #mfi.rats <- sprc(mylist)
# ### mosaic
# # MagField <- terra::mosaic(mfi.rats, fun="mean", overwrite=TRUE,
# #                           filename="C:/Covariates/Euskadi/ParentMaterial/icm.tif")
# # plotRGB(MagField, r=1, g=2, b=3, stretch="hist")
# 
# ###Change wd
# setwd("C:/Covariates/Euskadi/ParentMaterial/")
# 
# ### Mosaic band by band
# for(j in 1:3){
#     print(j)
#     ### Make list with one band from each raster stack
#   myfun = function(x) {
#     rsub <- subset(x, subset=j, drop=TRUE)
#     return(rsub)}
#   list.band <- lapply(mylist, FUN=myfun)
#     ## Assign function to mosaic
#   list.band$fun <- mean
#   ## Create mosaic 
#   band.mosaic <- do.call(mosaic,list.band)
#     writeRaster(band.mosaic,
#               na.rm=T, inf.rm=T, format="GTiff", overwrite=TRUE,
#               filename=paste0("icm_band",j,".tif"))
#   plot(band.mosaic)
# }
# 
# setwd("C:/Covariates/Euskadi/ParentMaterial/")
# mfi.files <- list.files(pattern="icm_band")
# mfi <- terra::rast(mfi.files)
# plotRGB(mfi, r=1, g=2, b=3, stretch="hist")
# plot(mfi)
# 
# ### Use focal mean to fill voids, several times
# mfi1   <- terra::focal(mfi[[1]], w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi1.2 <- terra::focal(mfi1, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi1.3 <- terra::focal(mfi1.2, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi1.4 <- terra::focal(mfi1.3, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi1.5 <- terra::focal(mfi1.4, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi1.6 <- terra::focal(mfi1.5, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi1.7 <- terra::focal(mfi1.6, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# plot(mfi1.7)
# writeRaster(mfi1.7, overwrite=TRUE, filename="icm_band1_filled.tif")
# 
# mfi2   <- terra::focal(mfi[[2]], w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi2.2 <- terra::focal(mfi2, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi2.3 <- terra::focal(mfi2.2, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi2.4 <- terra::focal(mfi2.3, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi2.5 <- terra::focal(mfi2.4, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi2.6 <- terra::focal(mfi2.5, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi2.7 <- terra::focal(mfi2.6, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# writeRaster(mfi2.7, overwrite=TRUE, filename="icm_band2_filled.tif")
# plot(mfi2.7)
# 
# mfi3   <- terra::focal(mfi[[3]], w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi3.2 <- terra::focal(mfi3, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi3.3 <- terra::focal(mfi3.2, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi3.4 <- terra::focal(mfi3.3, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi3.5 <- terra::focal(mfi3.4, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi3.6 <- terra::focal(mfi3.5, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# mfi3.7 <- terra::focal(mfi3.6, w=3, fun=mean, na.policy="only", na.rm=TRUE)
# writeRaster(mfi3.7, overwrite=TRUE, filename="icm_band3_filled.tif")
# plot(mfi3.7)
# 
# mfi.files <- list.files(pattern="filled")
# mfi <- terra::rast(mfi.files)
# plotRGB(mfi, r=1, g=2, b=3, stretch="hist")
# plot(mfi)

### But because it is hard to match the RGB colors to MFI anomaly, it is hard to use this layer... 
### how to extract the numerical value?
## But also.. it looks really bad


#############################################################################################

# 9. MCA on categorical variables ---------------------------------------------------------------------

### Vegetation MVG - 9 classes

### Load and stack
setwd("C:/Covariates/Euskadi/Organisms/")
mvg <- rast("mvg.tif")
plot(mvg)
### Transform into categorical raster
levels(mvg) [[1]]

### all categories
veg_pot <- rast("potential_veg.tif")
levels(veg_pot)[[1]]

### Geology, three layers
setwd("C:/Covariates/Euskadi/ParentMaterial/")

lithology <- terra::rast("lithology.tif")
plot(lithology)
lithology_lv1 <- rast("lithology_lv1.tif")
plot(lithology_lv1, col=viridis(n=7,direction = 1, option = "D"))
CODE_lito_v1 <- data.frame(ID=1:7, 
                           Litho_v1 = c("Superficial deposits",
                                        "Sedimentary detrital rocks",
                                        "Limestone, dolostone, marl, calcareous",
                                        "Clays with gypsum and other salts",
                                        "Metamorphic",
                                        "Volcanic and subvolcanic",
                                        "Granitoids"))
levels(lithology_lv1) <- CODE_lito_v1 ### as factor
is.factor(lithology_lv1)

lithology_lv2 <- rast("lithology_lv2.tif")
CODE_lito_v2 <- data.frame(ID=1:3, 
                           Litho_v2 = c("Sedimentary",
                                        "Metamorphic",
                                        "Igneous"))
levels(lithology_lv2) <- CODE_lito_v2 ### as factor
is.factor(lithology_lv2)
plot(lithology_lv2)

### Regolith
regolith <- rast("regolith.tif")
plot(regolith)
# regolith.lvls <- data.frame(ID=1:5,
#                             regolith =c("(1) - Espesor de 0 a 0.5 m",
#                                         "(2) - Espesor de 0.5 a 1 m",
#                                         "(3) - Espesor de 1 a 2 m",
#                                         "(4) - Espesor de 2 a 4 m",
#                                         "(5) - Espesor mayor de 4 m"))
# levels(regolith) <- regolith.lvls

## Create stack
cat.s <- c(mvg, lithology, lithology_lv1, regolith)
plot(cat.s)
### Rename variables
names(cat.s) <- c("MVG", "Lithology", "Lithology_v1", "Regolith")

### After exploring MCA with Potential Vegetation (29 variables) I decide to eliminate it
# cat.s <- c(mvg, lithology, lithology_lv1, regolith)
# plot(cat.s)
# ### Rename variables
# names(cat.s) <- c("PotVeg", "Lithology", "Lithology_v1", "Regolith")

### Generate a regular sample
# Set the random number generator to reproduce the results
set.seed(1946)
# Regular sampling
sampleCAT <- spatSample(cat.s, size = 5000000 , method="regular", 
                        replace=FALSE, as.df=TRUE, xy=TRUE, na.rm = TRUE)
dim(sampleCAT)

### It is a dataframe
### select only complete cases 
sampleCAT <-sampleCAT[complete.cases(sampleCAT),]
str(sampleCAT)
length(unique(sampleCAT$MVG))
length(unique(sampleCAT$Lithology))

### Perform MCA with MASS package/or FactoMineR package
library(MASS)
library(FactoMineR)
library(factoextra)

mca.cat2 <- MCA(sampleCAT[,3:ncol(sampleCAT)], ncp = 45, graph = TRUE)
eig.val <- get_eigenvalue(mca.cat2)
fviz_screeplot(mca.cat2, addlabels = TRUE, ylim = c(0, 45), ncp=31,"variance")
mca.eig.df <- as.data.frame(get_eig(mca.cat2))
mca.eig.df$dim <- 1:nrow(mca.eig.df)

# plot
ggplot(data = mca.eig.df, 
       aes(x = dim, y = cumulative.variance.percent)) +
    ggtitle("Cummulative variance explained by MCA dimensions")+
   geom_line()+
  geom_point(color="blue", size=2)

ggplot(data = mca.eig.df, 
       aes(x = dim, y = variance.percent)) +
  ggtitle("Variance explained by MCA dimensions")+
  geom_point()
### Keep the first 6 axes explaining 35% of variance. NOT VERY GOOD

ggplot(data = mca.eig.df, 
       aes(x = dim, y = eigenvalue)) +
  ggtitle("eigenvalue by MCA dimensions")+
  geom_point()

fviz_mca_var(mca.cat2,col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
fviz_mca_var(mca.cat2,col.var = "contrib",axes = c(3,4),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
fviz_mca_var(mca.cat2,col.var = "contrib",axes = c(5,6),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# Hide individuals
fviz_mca_biplot(mca.cat2, invisible ="ind")

# Select the top 10 variables
fviz_mca_biplot(mca.cat2, select.var = list(contrib = 10), invisible ="ind")

### Perform MCA with MASS package
mca.cat <- MASS::mca(df=sampleCAT[,3:ncol(sampleCAT)], nf=31) ### Keep ~100% variability with 32 variables
plot(mca.cat)
# eigenvalues
plot(1:31,mca.cat$d^2)

# number of categories per variable
cats = apply(sampleCAT[,3:ncol(sampleCAT)], 2, function(x) nlevels(as.factor(x)))

# data frame for ggplot
mca_vars_df = data.frame(mca.cat$cs, Variable = rep(names(cats), cats))

# plot
ggplot(data = mca_vars_df, 
       aes(x = X1, y = X2, label = rownames(mca_vars_df))) +
  geom_hline(yintercept = 0, colour = "gray70") +
  geom_vline(xintercept = 0, colour = "gray70") +
  geom_text(aes(colour = Variable)) +
  ggtitle("MCA plot of lithology and potential vegetation")

# Color by cos2 values: quality on the factor map
fviz_mca_var(mca.cat2, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # Avoid text overlapping
             ggtheme = theme_minimal())

### Predict the coordinates for all the pixels. In other words, map the MCA components
cat.s.df <- as.data.frame(cat.s) ### This worked with raster package. now with terra I have to adapt the code

### Can I do it directly with terra? The stack is small (25831920 pixels)
MCA.pred <- terra::predict(object=cat.s,
                           model=mca.cat, 
                           fun=predict,
                           na.rm=TRUE,
                           type = "row")
plot(MCA.pred)

### Write Raster
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/1.CovariatesEus/") # set wd
writeRaster(MCA.pred, 
            filename=paste0("MCA_",1:nlyr(MCA.pred),".tif"), 
            overwrite=TRUE)


# 10. MCA on parent material variables only -----------------------------------

## Create stack
cat.s <- c(lithology, lithology_lv1, regolith)
plot(cat.s)
### Rename variables
names(cat.s) <- c("Lithology", "Lithology_v1", "Regolith")

### Generate a regular sample
# Set the random number generator to reproduce the results
set.seed(1946)
# Regular sampling
sampleCAT <- spatSample(cat.s, size = 6000000 , method="regular", 
                        replace=FALSE, as.df=TRUE, xy=TRUE, na.rm = TRUE)
dim(sampleCAT)

### It is a dataframe
### select only complete cases 
sampleCAT <-sampleCAT[complete.cases(sampleCAT),]
str(sampleCAT)
length(unique(sampleCAT$Lithology))

### Perform MCA with MASS package/or FactoMineR package
library(MASS)
library(FactoMineR)
library(factoextra)

mca.cat2 <- MCA(sampleCAT[,3:ncol(sampleCAT)], ncp = 45, graph = TRUE)
eig.val <- get_eigenvalue(mca.cat2)
fviz_screeplot(mca.cat2, addlabels = TRUE, ylim = c(0, 45), ncp=31,"variance")
mca.eig.df <- as.data.frame(get_eig(mca.cat2))
mca.eig.df$dim <- 1:nrow(mca.eig.df)
### The first 8 dimensions explain 50 % of variability. That is a bit better

# plot
ggplot(data = mca.eig.df, 
       aes(x = dim, y = cumulative.variance.percent)) +
  ggtitle("Cummulative variance explained by MCA dimensions")+
  geom_line()+
  geom_point(color="blue", size=2)

ggplot(data = mca.eig.df, 
       aes(x = dim, y = variance.percent)) +
  ggtitle("Variance explained by MCA dimensions")+
  geom_point()
### Keep the first 8 axes

ggplot(data = mca.eig.df, 
       aes(x = dim, y = eigenvalue)) +
  ggtitle("eigenvalue by MCA dimensions")+
  geom_point()

fviz_mca_var(mca.cat2,col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
fviz_mca_var(mca.cat2,col.var = "contrib",axes = c(3,5),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
fviz_mca_var(mca.cat2,col.var = "contrib",axes = c(7,8),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# Hide individuals
fviz_mca_biplot(mca.cat2, invisible ="ind")


### Perform MCA with MASS package
mca.cat <- MASS::mca(df=sampleCAT[,3:ncol(sampleCAT)], nf=15) ### Keep ~100% variability with 32 variables
plot(mca.cat)
# eigenvalues
plot(1:15,mca.cat$d^2)

# number of categories per variable
cats = apply(sampleCAT[,3:ncol(sampleCAT)], 2, function(x) nlevels(as.factor(x)))

# data frame for ggplot
mca_vars_df = data.frame(mca.cat$cs, Variable = rep(names(cats), cats))

# plot
ggplot(data = mca_vars_df, 
       aes(x = X1, y = X2, label = rownames(mca_vars_df))) +
  geom_hline(yintercept = 0, colour = "gray70") +
  geom_vline(xintercept = 0, colour = "gray70") +
  geom_text(aes(colour = Variable)) +
  ggtitle("MCA plot of lithology and potential vegetation")

# Color by cos2 values: quality on the factor map
fviz_mca_var(mca.cat2, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # Avoid text overlapping
             ggtheme = theme_minimal())

### Predict the coordinates for all the pixels. In other words, map the MCA components
#cat.s.df <- as.data.frame(cat.s) ### This worked with raster package. now with terra I have to adapt the code

### Can I do it directly with terra? The stack is small (25831920 pixels)
MCA.pred <- terra::predict(object=cat.s,
                           model=mca.cat, 
                           fun=predict,
                           na.rm=TRUE,
                           type = "row")
plot(MCA.pred[[1:8]])
names(MCA.pred) <- paste0("pm_mca_",1:15)

### Write Raster
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/1.CovariatesEus/") # set wd
writeRaster(MCA.pred, 
            filename=paste0("pm_mca_",1:nlyr(MCA.pred),".tif"), 
            overwrite=TRUE)

### end of the script