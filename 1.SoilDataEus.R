#############################################################################################################################################
### author: Mercedes Roman Dobarco
### email: mercedes.roman@bc3research.org
### Objective:  Compile datasets of soil data for the Basque Country

### 1.Plot with leaflet

### Desired extent: Basque Country
### Resolution: Not decided yet. probably 30m
### CRS: EPSG=4326 (???)

###  Author: Mercedes Roman Dobarco
###  Date: 18/05/2023

####### Load packages

### Read and manage data
library(readxl)
library(xlsx)
library(dplyr)
library(tidyverse)

### Spatial
library(sf)
library(terra)
library(gdalUtilities)

### Visualization
library(lattice)
library(ggplot2)
library(viridis) # color palettes
library(scales)
library(rasterVis)
library(gridExtra)
library(rasterVis)
library(RColorBrewer)
#remotes::install_github('r-tmap/tmap')
library(tmap)    # for static and interactive maps
library(leaflet) # for interactive maps
library(mapview) # for interactive maps
library(shiny)   # for web applications

# 1. INES data ------------------------------------------------------------

### Read INES data
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/INES datuak (2018)")
ines_a <- read_excel("Araba/DatParcelasA.xlsx")
ines_b <- read_excel("Bizkaia/DatParcelas.xlsx")
ines_c <- read_excel("Gipuzkoa/DatParcelasG.xlsx")

setdiff(colnames(ines_a), colnames(ines_b))
setdiff(colnames(ines_a), colnames(ines_c))
setdiff(colnames(ines_c), colnames(ines_a))
setdiff(colnames(ines_b), colnames(ines_c))
setdiff(colnames(ines_c), colnames(ines_b))

### change colnames and create missing column in ines_b
colnames(ines_a)[colnames(ines_a) == "Probintzia"] <- "Provincia"
ines_b$Provincia <- "Bizkaia"

### ines dataset
ines <- full_join(ines_a,ines_b)
ines <- full_join(ines,ines_c)
rm(ines_a, ines_b, ines_c)
names(ines)

### Transform into spatial
ines_sf <- st_as_sf(ines, coords = c("CoordXpm","CoordYpm"), crs = 25830)

### Project
ines_WGS84 <- st_transform(ines_sf, 4326)
coords_wgs84 <- st_coordinates(ines_WGS84)
ines$Latitude <- coords_wgs84[,2]
ines$Longitude <- coords_wgs84[,1]
rm(ines_WGS84, coords_wgs84)

### INES has some field observations on soil physical and hydraulic properties,
### like erosion, permeability, soil structure

table(ines$Estructura)
# 0   1   2   3   4 
# 7 119  90  65  46 
table(ines$CodLitologiaCampo)
table(ines$CodPermeabilidad)
# 0   2   3   4   5   6 
# 7   5   1  76  40 198 

### Transform some variables
### Probably I will just keep clay, sand, silt, OM%, coarse fragments (%), bulk density (g/cm3)
### From an email from TRAGSATEC I know that
### "Las muestras de suelo del INES que se analizan en laboratorio son de los 10 cm superiores del suelo (0-10cm)."
ines$Upper_limit <- 0
ines$Lower_limit <- 10

### Change colnames or copy columns to have same name as in Carbosol
colnames(ines)[colnames(ines)== "Provincia"] <- "Province"
colnames(ines)[colnames(ines)== "CoordXpm"] <- "UTM_X"
colnames(ines)[colnames(ines)== "CoordYpm"] <- "UTM_Y"
colnames(ines)[colnames(ines)== "Altpm"] <- "Elevation"
colnames(ines)[colnames(ines)== "Arena"] <- "Sand"
colnames(ines)[colnames(ines)== "Limo"] <- "Silt"
colnames(ines)[colnames(ines)== "Arcilla"] <- "Clay"
colnames(ines)[colnames(ines)== "MO"] <- "OM"

### We don´t know for certain that the coarse fragments are gravimetric 
### but since they are determined in the laboratory this is the most likely
ines$CoarseFrag_Grav <- ines$EG 

### Transformation of bulk density to the correct units
ines$Bulk_Density <- ines$DA/1000
### Indicate information on Bulk Density method
ines$BD_method <- "Disturbed sample with coarse fragments"
ines$BD_quality <- "low"

### Transformation of OM to TOC % with conversion factor
ines$TOC <- ines$OM / 1.724

### I keep active carbonate because there is data from BASONET on this
ines$Act_carbonates <- ines$CalizaActiva

### I keep Thickness O horizon with new name
ines$Thickness_O <- ines$ProfMO

### Indicate dataset
ines$Dataset <- "INES"

### Indicate year of sampling. The exact date of fieldwork is unknown, 
### but the closest year is 2018.
ines$Date <- 2018

### Create my ID
ines$myID <- paste0(ines$Dataset,"_",ines$CodParcela )

### Subset of columns to merge
ines_sub <- ines[,c("myID","Dataset","Province","Date",
                    "UTM_X","UTM_Y","Longitude","Latitude",
                    "Upper_limit","Lower_limit",
                    "Sand","Silt","Clay",
                    "OM","TOC",
                    "Bulk_Density","BD_method", "BD_quality",
                    "CoarseFrag_Grav","Act_carbonates")]
summary(ines_sub$Bulk_Density)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.4500  0.8300  0.9300  0.9267  1.0300  1.2000       8
### Is bulk density of whole soil or fine earth?

### It is in between. It is calculated with particles > 2mm 
### but excluding big stones and roots.
### from a DISTURBED sample, a cylinder of known volume (Kopecky cylinder)
### Is filled with dry soil (before sieving) avoiding big stones and roots.
### Weight, and calculate bulk density after subtracting the weight of the cylinder.

### Procedimiento: Se rellena hasta el borde el cilindro de Kopecky, 
### de volumen conocido, con suelo seco y sin tamizar, 
### evitC!ndose elementos que no forman parte de ninguna de las fracciones del suelo,
### como piedras grandes o raC-ces. Se pesa y se anota el valor. 
### La densidad aparente serC! el cociente entre el peso medido y 
### el volumen fijo del cilindro de Kopecky.

### I will exclude the data of bulk density from data analysis

### INES 2018
colnames(ines_sub)
sel.columns.ines <- c("Dataset","myID","UTM_X","UTM_Y","Longitude","Latitude",
                      "Upper_limit","Lower_limit","Date",
                       # "Bulk_Density","BD_method","BD_quality",
                      "Sand","Silt","Clay","OM","TOC", ### Minimum dataset
                      "CoarseFrag_Grav","Act_carbonates")
ines_m <- ines_sub[,sel.columns.ines]


# 2. CARBOSOL -------------------------------------------------------------

setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/CARBOSOL/")
carbosol_horizons <- read_csv("carbosol_horizons_mrd.csv")
carbosol_profiles <- read_csv("carbosol_profiles_mrd.csv")
### Join
carbosol <- right_join(carbosol_profiles,carbosol_horizons)
#, by=c("Id_Profile","Province", "Location", "Latitude", "Longitude", "UTM_X", "UTM_Y"))
### Transform to spatial
carbosol_sf <- st_as_sf(carbosol, coords = c("UTM_X","UTM_Y"), crs = 25830)
rm(carbosol_horizons,carbosol_profiles)

### bulk density?
summary(carbosol$Bulk_Density)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.620   1.200   1.370   1.353   1.500   2.650   14208 
### These values agree more with bulk density of whole soil
summary(carbosol$EC)

### Sampling date?
summary(carbosol$Date)

### Load administrative boundaries of Euskadi
### the CAPV ("U:/Covariates/Administrative/CB_CAPV_5000_ETRS89/U11.CB_CAPV_5000_ETRS89.shp")
euskadi <- st_read("U:/Covariates/Administrative/CB_CAPV_5000_ETRS89/U11.CB_CAPV_5000_ETRS89.shp")

### and the rest of Spain
HomeDir <- "U:/Covariates/" ### Change this with your Home directory 

setwd(paste0(HomeDir,"Administrative/SHP_ETRS89/recintos_autonomicas_inspire_peninbal_etrs89/"))
autonomias <- st_read("recintos_autonomicas_inspire_peninbal_etrs89.shp")
autonomias25830 <- st_transform(autonomias, 25830)
st_bbox(autonomias25830)

### Countries 
world <- st_read("U:/Covariates/Europe/ref-countries-2020-01m.shp/CNTR_RG_01M_2020_4326.shp/CNTR_RG_01M_2020_4326.shp")
world25830 <- st_transform(world, 25830)

### subset neighbours around Euskadi
buff_Eus <- st_buffer(euskadi, dist=1000000)
plot(buff_Eus["CCAA"])
world_eus_int <- st_intersects(world25830, buff_Eus, sparse = FALSE)
world_eus <- world25830[world_eus_int,]
plot(world_eus["CNTR_ID"])
rm(buff_Eus,world,world25830,world_eus_int)

### Grey fill
tm_shape(world_eus, bbox=st_bbox(euskadi), crs=25830) +
  tm_borders() +
  tm_shape(autonomias25830) +
  tm_fill(col="gray95") +
  tm_borders() +
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_scalebar(breaks = c(0,20,40), text.size = 0.8, position = c("right", "bottom")) +
  tm_grid(lines=FALSE) +
  tm_shape(carbosol_sf) +
  tm_symbols(shape=21, 
             size=0.3,
             legend.shape.show = TRUE,
             col="indianred3", 
             border.col="black", 
             title.shape = "CARBOSOL") +
  tm_shape(ines_sf) +
  tm_symbols(shape=24, 
             size=0.2, 
             legend.shape.show = TRUE,
             col="cornflowerblue", 
             border.col="black",
             title.shape = "INES") +
  tm_layout(legend.outside=TRUE,
            legend.outside.position="right", 
            legend.title.size=1.5,
            legend.text.size=1.1)

eus_buffer <- read_sf('C:/Covariates/Euskadi/boundaries/eus_buffer_125m.shp')
eus_buffer_WGS84 <- st_transform(eus_buffer, 4326)
tm_shape(eus_buffer, bbox=st_bbox(euskadi), crs=25830) + tm_borders()

### I don't use the buffer to subset Carbosol locations, 
### but the  administrative limits 
carbosol_Eus <- st_intersection(carbosol_sf, euskadi)
### This leaves 907 observations (I lose 36 compared if I use a buffer to subset)

tm_shape(world_eus, bbox=st_bbox(euskadi), crs=25830) +
  tm_borders() +
  tm_shape(autonomias25830) +
 # tm_fill(col="gray95") +
  tm_borders() +
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_scalebar(breaks = c(0,20,40), text.size = 0.8, position = c("right", "bottom")) +
  tm_grid(lines=FALSE) +
  tm_shape(carbosol_Eus) +
  tm_symbols(shape=21, 
             size=0.3,
             legend.shape.show = TRUE,
             col="indianred3", 
             border.col="black", 
             title.shape = "CARBOSOL") +
  tm_shape(ines_sf) +
  tm_symbols(shape=24, 
             size=0.2, 
             legend.shape.show = TRUE,
             col="cornflowerblue", 
             border.col="black",
             title.shape = "INES") +
  tm_layout(legend.outside=TRUE,
            legend.outside.position="right", 
            legend.title.size=1.5,
            legend.text.size=1.1)

### Drop the geometry. Back to dataframe
carbosol_Eus_df <- carbosol_Eus
st_geometry(carbosol_Eus_df) <- NULL

coords_carbosol_utm <- st_coordinates(carbosol_Eus)
carbosol_Eus_df$UTM_X <- coords_carbosol_utm[,1]
carbosol_Eus_df$UTM_Y <- coords_carbosol_utm[,2]
carbosol_Eus_df <- as.data.frame(carbosol_Eus_df)

### bbox for plot larger than just Euskadi
bbox_plot_Eus <- st_buffer(euskadi, dist=20000)

tm_shape(world_eus, bbox=st_bbox(bbox_plot_Eus), crs=25830) +
  tm_borders() +
  tm_shape(autonomias25830) +
 # tm_fill(col="gray95") +
  tm_borders() +
  tm_compass(type = "arrow", 
             position = c("left", "top")) +
  tm_scalebar(breaks = c(0,20,40),
               text.size = 0.6, 
               position = c("left", "bottom")) +
  tm_grid(lines=FALSE) +
  tm_shape(carbosol_Eus) +
  tm_symbols(col="Date",
             palette="viridis",
             size=0.4,
             title.shape = "CARBOSOL") +
  tm_layout(inner.margins = c(0,0,0,0),
            legend.outside=TRUE,
            legend.outside.position="right", 
            legend.title.size=1.2,
            legend.text.size=0.8)

### How many locations?
length(unique(carbosol_Eus_df$Id_Profile)) ### 252 locations
length(unique(carbosol_Eus_df$Id_Horiz)) ### 907 horizons

library(ggplot2)
ggplot()+
  geom_histogram(aes(x=Upper_limit_m), data=carbosol_Eus_df, 
                 binwidth = 0.05,
                 fill="cadetblue3", col="gray30") +
  labs(x="Upper horizon depth (m)")+
  theme_bw()

ggplot()+
  geom_histogram(aes(x=Lower_limit_m), data=carbosol_Eus_df, 
                 binwidth = 0.05,
                 fill="cadetblue3", col="gray30") +
  labs(x="Lower horizon depth (m)")+
  theme_bw()

ggplot()+
  geom_histogram(aes(x=Date), data=carbosol_Eus_df, 
                 fill="cadetblue3", col="gray30") +
  labs(x="Date")+
  theme_bw()

### Indicate dataset
carbosol_Eus_df$Dataset <- "CARBOSOL"

### Create my ID
carbosol_Eus_df$myID <- paste0(carbosol_Eus_df$Dataset,"_",carbosol_Eus_df$Id_Profile)

### Horizon limits to cm
carbosol_Eus_df$Upper_limit <- carbosol_Eus_df$Upper_limit_m * 100
carbosol_Eus_df$Lower_limit <- carbosol_Eus_df$Lower_limit_m * 100

### Colnames
colnames(carbosol_Eus_df)

### CARBOSOL metadata (from "Metadata.xlxs" https://store.pangaea.de/Publications/LlorenteM_2017/Metadata.xlsx)

# "Id_Profile"  Unique identification number of profile
# "Id_ref"      Unique identification number of reference
# "Date"        Year soil sampling
# "Province"    Province
# "Location"    Municipality or local toponym
# "Latitude"    Latitude in WGS84
# "Longitude"   Longitude in WGS84
# "N_hor"       Number of horizons in the profile
# "N_hor_o"     Number of horizons in the profile with organic matter measurement
# "depth_m"     Depth (m) ### Note that in the metadata it says cm but this is not correct       
# "LCC"         Land cover class      
# "LCC_code"    Land cover class code      
# "CORINE_l1"          CORINE Land Cover Code Level 1
# "CORINE_l23"         CORINE Land Cover Code Level 2-3  
# "Vegetation_ref"     Vegetation details provided by original reference  
# "Elevation"          Altitude (meters above sea level)   
# "Aspect"             Orientation  
# "Slope_perc"         Slope (%)
# "ParentMaterial_ref" Parent material by original reference
# "PM_consistency"     Parent material consistency
# "PM_silica"          Parent material silica content
# "SoilClass_WRB"      Soil classification into WRB system
# "SoilClass_USDA"     Soil classification into USDA system 
# "Id_Horiz"           Unique identification number of horizon
# "Hor_descrip"        Horizon description/designation in the original reference
# "Hor_Pos"            Horizon position in the soil profile
# "RelPosit_profile"   Horizon relative position in the soil profile (qualitative)
# "Upper_limit_m"      Upper limit (m) ### Note that in the metadata it says cm but it is in m
# "Lower_limit_m"      Lower limit (m) ### Note that in the metadata it says cm but it is in m
# "Mid_Depth_m"        Midpoint of horizon depth 
# "Thickness_m"        Horizon thickness (m)
# "UTM_X"              UTM X coordinates
# "UTM_Y"              UTM Y coordinates
# "Dataset"            Dataset
# "myID"               Unique ID created for this data analysis
# "Upper_limit"        Upper limit in cm
# "Lower_limit"        Lower limit in cm


### Soil properties - CARBOSOL METADATA
# "Color_HLS"    Munsell Code
# "Bulk_Density" Bulk density g/cm3 Wide variety of methods, including application of pedotransfer functions. Block Method is the most widely used
# "CoarseFrag"   Coarse material (% VOLUMETRIC) > 2 mm; % of total volume
# "Sand"         Sand (%)  USDA System  (2mm - 50 microm) 
# "Silt"         Silt (%)  USDA System  (50-2 microm)           
# "Clay"         Clay (%)  USDA System  (< 2 microm)
# "OM"           Organic matter (%) - Conversion factor of 1.724 applied when needed
# "TOC"          TOC (%) Walkley-Black Method or Dry Combustion Methods (elemental analyzer)
# "pH"           pH in water, Soil:water ratio 1:1, 1:2.5 or 1:5
# "Carbonates"   Carbonates (%) Bernard’s Calcimeter or Titrimetry
# "C_N"          Carbon and nitrogen ratio - Direct calculation, from TOC and N data   
# "TN_ppm"       Nitrogen ppm - Kjeldahl Method or Dry Combustion Methods (elemental analyzer)
# "P_ppm"        Phosphorus ppm - Olsen or Mehlich Methods 
# "K_ppm"        Potasium ppm - Extraction with ammonium acetate and quantification by ICP spectroscopy or Flame Photometry Methods
# "Ca_ppm"       Calcium ppm - Extraction with ammonium acetate and quantification by ICP spectroscopy or Atomic Absorption Spectrometry Methods 
# "Mg_ppm"       Magnesium ppm - Extraction with ammonium acetate and quantification by ICP spectroscopy or Atomic Absorption Spectrometry Methods 
# "Na_ppm"       Sodium ppm - Extraction with ammonium acetate and quantification by ICP spectroscopy or Flame Photometry Methods  
# "CEC"          Cation-exchange capacity (cmol/Kg) - Cation Summation Method
# "EC"           Electric Conductivity (dS/m at 25C) - Soil:water ratio 1:5
# "Gypsum"       Gypsum content (%) - Conductometry in the water extract

### Indicate bulk density quality and method
carbosol_Eus_df$BD_method  <- "Core method and PTF"
carbosol_Eus_df$BD_quality <- "medium"

### Bulk density in Euskadi
summary(carbosol_Eus_df$Bulk_Density)
### quite high the maximum value!
hist(carbosol_Eus_df$Bulk_Density, breaks=30)

### We have soil order information of 106 sites
unique(carbosol_Eus_df[!is.na(carbosol_Eus_df$SoilClass_WRB),]$myID)

### Here for CARBOSOL unless I find I need to correct something else
### CARBOSOL
sel.columns.carbosol <- c("Dataset","myID","UTM_X","UTM_Y","Longitude","Latitude",
                          "Upper_limit","Lower_limit","Date",
                          "Bulk_Density","BD_method","BD_quality",
                          "Sand","Silt","Clay","OM","TOC", ### Minimum dataset
                          "Color_HLS","Carbonates",
                          "pH","EC","TN_ppm","P_ppm","K_ppm","Ca_ppm","Mg_ppm","Na_ppm","C_N",
                          "CEC","Gypsum","CoarseFrag", "SoilClass_WRB","SoilClass_USDA",
                          "ParentMaterial_ref","PM_consistency","PM_silica")
carbosol_m <- carbosol_Eus_df[,sel.columns.carbosol]


# 3.1 BASONET 2021 --------------------------------------------------------

library(readxl)
Basonet_2021 <- read_excel("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/Basonet/BasonetSuelos_mrd.xlsx",
                           sheet = "medicion2021-2022")
Basonet_2011 <- read_excel("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/Basonet/BasonetSuelos_mrd.xlsx",
                           sheet = "medicion2011")
Basonet_2001 <- read_excel("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/Basonet/BasonetSuelos_mrd.xlsx",
                           sheet = "medicion2001")
Basonet_Texture_2001 <- read_excel("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/Basonet/BasonetSuelos_mrd.xlsx",
                                   sheet = "texturas2001")
colnames(Basonet_2021)
# [1] "UTM_X"                  "UTM_Y"                  "th"                     "Plot"                   "Horizon_Depth_Interval"
# [6] "Unique_ID"              "Lab_ID"                 "pH"                     "Bulk_Density"           "Ex_Acidity_vol"        
# [11] "CECef_vol"              "Al_CECef"               "Carbonates"             "P_HCl_mg-l"             "OM"                    
# [16] "N_perc"                 "C_N"                    "EC_µS_cm"               "P_Olsen_mg-l"           "Ca_mg_l"               
# [21] "K_mg_l"                 "Mg_mg_l"                "Na_mg_l" 

setdiff(colnames(Basonet_2001), colnames(Basonet_2021)) ### In 2001 they had mesurements in ppm
setdiff(colnames(Basonet_2021), colnames(Basonet_2001)) ### In 2021 they measured also,
# "Ex_Acidity_vol", "CECef_vol" , "Al_CECef" , "EC_µS_cm","Na_mg_l"

### BASONET 2021 variables of interest
# "UTM_X" "UTM_Y" "Plot" "Horizon_Depth_Interval" 
# "pH"               pH - Water extraction (ratio 1:2.5 v/v) with a pH-meter
# "Bulk_Density"     to transform from mg/l to ppm - disturbed, sieved sample. it is not soil bulk density in field conditions
# "Ex_Acidity_vol"   Exchangeable acidity - meq/100 ml	Extraction with 0.6N barium chloride. Titration with 0.01N NaOH. (extraction ratio of 1:10 v/v)
# "CECef_vol"        Effective CEC (meq/100 ml) - 	Calculation. Effective CEC from exchangeable cations K+, Mg2+, Ca2+, Na+, H+ and Al3+
# "Al_CECef"         Al in Effective CEC (%) - Calculation.% Al in effective CEC
# "Carbonates"       % CaCO3 - Bernards calcimeter. Volume of CO2 released when adding 50% HCl. (Bernard's calcimeter)
###                  It is only analysed in samples with pH>7.5 
# "P_HCl_mg-l"       Phosphorus extracted in 3% HCl (mg/l P). Extraction in 3% HCl (extraction ratio 1:20 v/v) UV-VIS spectrophotometer
# "OM"               Oxidised Organic Matter (%). Oxidation of organic matter with 1N potassium dichromate 
###                  and concentrated sulfuric acid. Titration with Mohr's Salt 0.5N. (Walkley-Black Method)
# "N_perc"           Nitrogen Kjeldahl (%). Digestion in sulfuric acid and subsequent distillation in basic medium. 
###                  The distillate is collected in ammonia fixing solution and titrated with 0.1N HCl.
# "C_N"              Ratio C to N
# "EC_µS_cm"         Electric conductivity - Extraction in saturated calcium sulphate (extraction ratio 1:2.5 v/v) - Conductivity meter.
# "P_Olsen_mg-l"     Olsen extractable Phosphorus. Extraction in 0.5N sodium bicarbonate (extraction ratio 1:20 v/v)-UV-VIS spectrophotometer
# "Ca_mg_l"
# "K_mg_l" 
# "Mg_mg_l" 
# "Na_mg_l"         Extractable Calcium, Magnesium, Pottasium, and Sodium in Ammonium acetate		(mg/l)
###                 Extraction in 1N ammonium acetate (extraction ratio 1:20 v/v) OPTICAL ICP

### Harmonise variables to same units as CARBOSOL
### Indicate dataset
Basonet_2021$Dataset <- "BASONET"
### Create my ID
Basonet_2021$myID <- paste0(Basonet_2021$Dataset,"_",Basonet_2021$th,"_",Basonet_2021$Plot)

### Horizon limits
table(Basonet_2021$Horizon_Depth_Interval)
Basonet_2021$Upper_limit <- NA
Basonet_2021$Lower_limit <- NA
Basonet_2021$Upper_limit <- ifelse(Basonet_2021$Horizon_Depth_Interval == "0-20", 0,
                                   ifelse(Basonet_2021$Horizon_Depth_Interval == "20-40", 20, NA))
Basonet_2021$Lower_limit <- ifelse(Basonet_2021$Horizon_Depth_Interval == "0-20", 20,
                                   ifelse(Basonet_2021$Horizon_Depth_Interval == "20-40", 40, NA))

### Date
Basonet_2021$Date <- 2021

### Exchangeable acidity from meq/100 cm3 soil to meq/100 g soil
Basonet_2021$Ex_Acidity <- NA
Basonet_2021$Ex_Acidity <- Basonet_2021$Ex_Acidity_vol / Basonet_2021$Bulk_Density

### CECef from meq/100 cm3 soil to meq/ 100 g soil (equivalent to cmol/kg)
Basonet_2021$CECef <- NA
Basonet_2021$CECef <- Basonet_2021$CECef_vol / Basonet_2021$Bulk_Density

# "Carbonates"  Only analysed in samples with pH>7.5. 
table(Basonet_2021$Carbonates)
Basonet_2021[Basonet_2021$pH <= 7.5,]$Carbonates
Basonet_2021[Basonet_2021$pH > 7.5,]$Carbonates
summary(Basonet_2021[is.na(Basonet_2021$Carbonates),]$pH)
Basonet_2021[Basonet_2021$Carbonates == "-", ]$pH
### Assign carbonates == 0 % for pH < 7.5
Basonet_2021[Basonet_2021$pH <= 7.5,]$Carbonates <- "0"
### when carbonates "<1" I assign 0 % (In the previous version I assigned 1%)
Basonet_2021[Basonet_2021$Carbonates == "<1",]$Carbonates <- "0"
### Transform to numeric
Basonet_2021$Carbonates <- gsub(x= Basonet_2021$Carbonates, pattern=",", replacement = "."  )
Basonet_2021$Carbonates <- as.numeric(Basonet_2021$Carbonates)

# OM is in % as in Carbosol
# Conversion into TOC in %
Basonet_2021$TOC <- NA
Basonet_2021$TOC <- Basonet_2021$OM / 1.724

# N_perc"
### TN from % to ppm
Basonet_2021$TN_ppm <- NA
Basonet_2021$TN_ppm <- Basonet_2021$N_perc * 10000

# "EC_µS_cm" ### EC in Basonet in µS/cm transform to dS/m as in Carbosol
Basonet_2021$EC <-NA
Basonet_2021$EC <- Basonet_2021$EC_µS_cm / 1000

### Transform extractable cations into ppm of a dry soil weight basis dividing by bulk density
### what do I do with the ones <20? In this case I assigned NA
summary(Basonet_2021$Ca_mg_l) ### just 20 NA
Basonet_2021$Ca_ppm <- NA
Basonet_2021[!is.na(Basonet_2021$Ca_mg_l),]$Ca_ppm <-
  Basonet_2021[!is.na(Basonet_2021$Ca_mg_l),]$Ca_mg_l / 
  Basonet_2021[!is.na(Basonet_2021$Ca_mg_l),]$Bulk_Density

table(Basonet_2021$Mg_mg_l) ### 129 values <20. But there are some readings with 10...
Basonet_2021$Mg_ppm <- NA
Basonet_2021[Basonet_2021$Mg_mg_l == "< 20",]$Mg_mg_l <- "<20"
Basonet_2021[!is.na(Basonet_2021$Mg_mg_l) & Basonet_2021$Mg_mg_l != "<20", ]$Mg_ppm <-
  as.numeric(Basonet_2021[!is.na(Basonet_2021$Mg_mg_l) & Basonet_2021$Mg_mg_l != "<20", ]$Mg_mg_l) /
  Basonet_2021[!is.na(Basonet_2021$Mg_mg_l) & Basonet_2021$Mg_mg_l != "<20", ]$Bulk_Density

summary(Basonet_2021$K_mg_l) ### 48 NAs but no characters
Basonet_2021$K_ppm <- NA
Basonet_2021[!is.na(Basonet_2021$K_mg_l),]$K_ppm <-
  Basonet_2021[!is.na(Basonet_2021$K_mg_l),]$K_mg_l /
  Basonet_2021[!is.na(Basonet_2021$K_mg_l),]$Bulk_Density

table(Basonet_2021$Na_mg_l) ## 565 as "<20" I wonder if set that to 0 or leave as NA
Basonet_2021$Na_ppm <- NA
Basonet_2021[Basonet_2021$Na_mg_l == "< 20", ]$Na_mg_l <- "<20"
Basonet_2021[!is.na(Basonet_2021$Na_mg_l) & Basonet_2021$Na_mg_l != "<20", ]$Na_ppm <-
  as.numeric(Basonet_2021[!is.na(Basonet_2021$Na_mg_l) & Basonet_2021$Na_mg_l != "<20", ]$Na_mg_l) /
  Basonet_2021[!is.na(Basonet_2021$Na_mg_l) & Basonet_2021$Na_mg_l != "<20", ]$Bulk_Density

summary(as.numeric(Basonet_2021[!is.na(Basonet_2021$Na_mg_l) & Basonet_2021$Na_mg_l != "<20", ]$Na_mg_l))
hist(as.numeric(Basonet_2021[!is.na(Basonet_2021$Na_mg_l) &
                               Basonet_2021$Na_mg_l != "<20", ]$Na_mg_l),
     breaks=100,
     xlab="Na mg/l", xlim=c(0,100))
hist(Basonet_2021$Na_ppm, breaks=100)
summary(Basonet_2021$Na_ppm)

# "P_HCl_mg-l"
table(Basonet_2021$`P_HCl_mg-l`) ### I assign NA to those below limit detection and the max value to those >120
### Or should I assign 0?
### I assign 0 to those below detection limit (again, different from previous version)
Basonet_2021$P_HCl_ppm <- NA
Basonet_2021[Basonet_2021$`P_HCl_mg-l` == "<4,8",]$`P_HCl_mg-l` <- "<4.8"
Basonet_2021[Basonet_2021$`P_HCl_mg-l` == ">120",]$`P_HCl_mg-l` <- "120"
Basonet_2021[Basonet_2021$`P_HCl_mg-l` != "<4.8",]$P_HCl_ppm <-
  as.numeric(Basonet_2021[Basonet_2021$`P_HCl_mg-l` != "<4.8",]$`P_HCl_mg-l`)/
  Basonet_2021[Basonet_2021$`P_HCl_mg-l` != "<4.8",]$Bulk_Density
Basonet_2021[Basonet_2021$`P_HCl_mg-l` == "<4.8",]$P_HCl_ppm <- 0
summary(Basonet_2021$P_HCl_ppm)

# "P_Olsen_mg-l"
### Phosphorus from mg/l to ... again, assign NA to those values below detection limit
### I assign 0 to those below detection limit (again, different from previous version)
table(Basonet_2021$`P_Olsen_mg-l`)
Basonet_2021$P_ppm <- NA
Basonet_2021[Basonet_2021$`P_Olsen_mg-l` == "<4,8",]$`P_Olsen_mg-l` <- "<4.8"
Basonet_2021[Basonet_2021$`P_Olsen_mg-l` != "<4.8",]$P_ppm <-
  as.numeric(Basonet_2021[Basonet_2021$`P_Olsen_mg-l` != "<4.8",]$`P_Olsen_mg-l`)/
  Basonet_2021[Basonet_2021$`P_Olsen_mg-l` != "<4.8",]$Bulk_Density
Basonet_2021[Basonet_2021$`P_Olsen_mg-l` == "<4.8",]$P_ppm <- 0
summary(Basonet_2021$P_ppm)

### Join texture from year 2001
colnames(Basonet_2021)[colnames(Basonet_2021) %in% colnames(Basonet_Texture_2001)]
table(Basonet_Texture_2001$Lower_limit)

#Basonet_2021 <- left_join(Basonet_2021, Basonet_Texture_2001, by= c("th","Plot","Lower_limit"))
#
# Warning message:
#   In left_join(Basonet_2021, Basonet_Texture_2001, by = c("th", "Plot",  :
#                                                             Detected an unexpected many-to-many relationship between `x` and `y`.
#                                                           ℹ Row 198 of `x` matches multiple rows in `y`.
#                                                           ℹ Row 214 of `y` matches multiple rows in `x`.
#                                                           ℹ If a many-to-many relationship is expected, set `relationship = "many-to-many"` to silence this warning.
### There are a couple of mistatches
as.data.frame(Basonet_2021[198,])
Basonet_Texture_2001[Basonet_Texture_2001$th==1 & Basonet_Texture_2001$Plot==1198 & Basonet_Texture_2001$Lower_limit== 20,]
Basonet_Texture_2001[214,]
Basonet_2021[Basonet_2021$th==1 & Basonet_2021$Plot==1198 & Basonet_2021$Lower_limit== 20,]
### Only one record for Basonet_Texture_2001[214,]

### Calculate average value and substitute
Basonet_Texture_2001[Basonet_Texture_2001$th==1 & Basonet_Texture_2001$Plot==1198 & Basonet_Texture_2001$Lower_limit== 20,] %>%
  group_by(., th, Plot, Lower_limit) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

### I use a new dataframe
Basonet_Texture_2001_av <- Basonet_Texture_2001 %>%
  group_by(., th, Plot, Lower_limit) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
Basonet_Texture_2001[Basonet_Texture_2001$th==1 &
                       Basonet_Texture_2001$Plot==1198 &
                       Basonet_Texture_2001$Lower_limit== 20,]
Basonet_Texture_2001_av[Basonet_Texture_2001_av$th==1 &
                          Basonet_Texture_2001_av$Plot==1198 &
                          Basonet_Texture_2001_av$Lower_limit== 20,]
### Join with average texture
Basonet_2021 <- left_join(Basonet_2021, Basonet_Texture_2001_av, by= c("th","Plot","Lower_limit"))

### Transform to spatial
### Seems that BASONET 2021 was in ETRS89 UTM30 --> epsg:25830
basonet_2021_sf <- st_as_sf(Basonet_2021, coords = c("UTM_X","UTM_Y"), crs = 25830)
#basonet_2021_sf1 <- st_as_sf(Basonet_2021, coords = c("UTM_X","UTM_Y"), crs = 25830)
#basonet_2021_sf <- st_as_sf(Basonet_2021, coords = c("UTM_X","UTM_Y"), crs = 23030)
plot(basonet_2021_sf["TOC"])

### Project
#basonet_2021_WGS84_1 <- st_transform(basonet_2021_sf1, 4326)
basonet_2021_WGS84 <- st_transform(basonet_2021_sf, 4326)

tm_shape(eus_buffer_WGS84, bbox=st_bbox(eus_buffer_WGS84), crs=4326) +
  tm_borders() 
  tm_grid(lines=FALSE) +
  tm_shape(basonet_2021_WGS84) +
  tm_symbols(shape=19,
             size=0.3,
             legend.shape.show = TRUE,
             border.col="black",
             title.shape = "ED50") +
    # tm_shape(basonet_2021_WGS84_1) +
    # tm_symbols(shape=1,
    #            size=0.3,
    #            legend.shape.show = TRUE,
    #            border.col="red",
    #            title.shape = "ETRS89")  +
  tm_layout(inner.margins = c(0,0,0,0),
            legend.outside=TRUE,
            legend.outside.position="right",
            legend.title.size=1.2,
            legend.text.size=0.8)


coords_wgs84 <- st_coordinates(basonet_2021_WGS84)
Basonet_2021$Latitude <- coords_wgs84[,2]
Basonet_2021$Longitude <- coords_wgs84[,1]
rm(basonet_2021_WGS84, coords_wgs84)

tm_shape(world_eus, bbox=st_bbox(bbox_plot_Eus), crs=25830) +
  tm_borders() +
  tm_shape(autonomias25830) +
  #tm_fill(col="gray95") +
  tm_borders() +
   tm_grid(lines=FALSE) +
  tm_shape(basonet_2021_sf) +
  tm_symbols(shape=21,
             size=0.3,
             legend.shape.show = TRUE,
             col="Mg_ppm",
             palette="viridis",
             border.col="black",
             title.shape = "BASONET")  +
  tm_layout(inner.margins = c(0,0,0,0),
            legend.outside=TRUE,
            legend.outside.position="right",
            legend.title.size=1.2,
            legend.text.size=0.8)

### "Bulk_Density" - Lab method defined as "weight of 20 cm3 of soil"
### Include information about bulk density method
### Indicate information on Bulk Density method
Basonet_2021$BD_method <- "Disturbed sample fine fraction"
Basonet_2021$BD_quality <- "Unacceptable"

### Here for BASONET 2021 unless I find I need to correct something else

### BASONET 2021
sel.columns.basonet2021 <- c("Dataset","myID","UTM_X","UTM_Y","Longitude","Latitude",
                             "Upper_limit","Lower_limit","Date",
                             #"Bulk_Density","BD_method","BD_quality",
                             "Sand","Silt","Clay","OM","TOC", ### Minimum dataset
                             "pH","EC","TN_ppm","C_N",
                             "P_ppm","K_ppm","Ca_ppm","Mg_ppm","Na_ppm",
                             "Ex_Acidity","CECef","Al_CECef",
                             "Carbonates","P_HCl_ppm")
Basonet_2021_m <- Basonet_2021[,sel.columns.basonet2021]


# 3.2 BASONET 2001 --------------------------------------------------------

Basonet_2001 <- read_excel("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/Basonet/BasonetSuelos_mrd.xlsx",
                           sheet = "medicion2001")
Basonet_Texture_2001 <- read_excel("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/Basonet/BasonetSuelos_mrd.xlsx",
                                   sheet = "texturas2001")
colnames(Basonet_2001)
### Variables of interest
# "UTM_X"              
# "UTM_Y" 
# "Plot"              
# "Lower_limit"         Lower horizon limit (0 cm or 20 cm)
# "pH"                  pH (1:2.5 V/V)	- Water extraction (ratio 1:2.5 v/v) with a pH-m - Leave as is
# "Carbonates"          % CaCO3 (Bernards calcimeter). It is only analyzed in samples with pH>7.5 Volume of CO2 released when adding 50% HCl. (Bernard's calcimeter)      
# "Act_Carbonate_perc"  Active carbonates (%)
# "OM"                  Oxidised Organic Matter (%)- 	Oxidation of organic matter with 1N potassium dichromate and concentrated sulfuric acid. Titration with Mohr's Salt 0.5N. (Walkley-Black Method)
# "N_perc"              Nitrogen Kjeldahl (%)	-	Digestion in sulfuric acid and subsequent distillation in basic medium. The distillate is collected in ammonia fixing solution and titrated with 0.1N HCl.
# "C_N"                 Ratio C to N
# "P_Olsen_ppm"         Bicarbonate extractable phosphorus (mg/l) - Extraction in 0.5N sodium bicarbonate (extraction ratio 1:20 v/v)-UV-VIS spectrophotometer
# "P_HCl_ppm"           Phosphorus extracted in 3% HCl	(mg/l) - Extraction in 3% HCl (extraction ratio 1:20 v/v) UV-VIS spectrophotometer
# "K_ppm"               Extractable Calcium, Magnesium, Pottasium, and Sodium in Ammonium acetate		(ppm) Extraction in 1N ammonium acetate (extraction ratio 1:20 v/v) OPTICAL ICP  
# "Ca_ppm"                       
# "Mg_ppm"             

### Harmonise variables to same units as CARBOSOL
### Indicate dataset
Basonet_2001$Dataset <- "BASONET"
# ### Create my ID
# Basonet_2001$myID <- paste0(Basonet_2001$Dataset,"_",Basonet_2001$Plot)
### Create my ID
Basonet_2001$myID <- paste0(Basonet_2001$Dataset,"_",Basonet_2001$th,"_",Basonet_2001$Plot)
### Indicate date
Basonet_2001$Date <- 2001
### Horizon limits
table(Basonet_2001$Lower_limit)
Basonet_2001$Upper_limit <- NA
Basonet_2001$Upper_limit <- ifelse(Basonet_2001$Lower_limit == 20, 0,
                                   ifelse(Basonet_2001$Lower_limit == 40, 20, NA))

### Variables in BASONET 2001 - same as CARBOSOL

### pH
summary(Basonet_2001$pH)

### "Bulk_Density" - Lab method defined as "weight of 10 cm3 of soil"
### Indicate information on Bulk Density method
Basonet_2001$BD_method <- "Disturbed sample fine fraction"
Basonet_2001$BD_quality <- "Unacceptable"

# "Carbonates"  Only analysed in samples with pH>7.5. CO2 volume released when adding HCl 50% (Bernard calcimeter)
table(Basonet_2001$Carbonates)
Basonet_2001[Basonet_2001$pH <= 7.5,]$Carbonates
Basonet_2001[Basonet_2001$pH > 7.5,]$Carbonates
summary(Basonet_2001[is.na(Basonet_2001$Carbonates),]$pH)
Basonet_2001[Basonet_2001$Carbonates == "-", ]$pH
### Assign carbonates == 0 % for pH < 7.5
Basonet_2001[Basonet_2001$pH <= 7.5 & is.na(Basonet_2001$Carbonates),]$Carbonates <- 0
### Assign carbonates == 0 % for pH > 7.5 and NA
Basonet_2001[Basonet_2001$pH > 7.5 & is.na(Basonet_2001$Carbonates),]$Carbonates <- 0
# ### when carbonates "<1" I assign 0 %
# Basonet_2001[Basonet_2001$Carbonates == "<1",]$Carbonates <- "0"

# OM is in % as in Carbosol
# Conversion into TOC in %
Basonet_2001$TOC <- NA
Basonet_2001$TOC <- Basonet_2001$OM / 1.724

# [16] "N_perc"
### TN from % to ppm
Basonet_2001$TN_ppm <-NA
Basonet_2001$TN_ppm <- Basonet_2001$N_perc * 10000

# "C_N" as it is
# "Ca_ppm" exists
# "K_ppm" exists
# "Mg_ppm" exists
# "Na_ppm" not measured
### P_Olsen_ppm
### P_HCl_ppm

### Join with average texture from year 2001
Basonet_2001 <- left_join(Basonet_2001, Basonet_Texture_2001_av, by= c("th","Plot","Lower_limit"))

### Rename Active carbonates column
summary(Basonet_2001$Act_Carbonate_perc)
colnames(Basonet_2001)[colnames(Basonet_2001) == "Act_Carbonate_perc"] <- "Act_carbonates"
colnames(Basonet_2001)[colnames(Basonet_2001) =="P_Olsen_ppm"] <- "P_ppm"

### Transform to spatial
### It seems in 2001 the datum was ED50
basonet_2001_sf <- st_as_sf(Basonet_2001, coords = c("UTM_X","UTM_Y"), crs = 23030)
#basonet_2001_sf <- st_as_sf(Basonet_2001, coords = c("UTM_X","UTM_Y"), crs = 25830)
### Latitude and longitude
basonet_2001_WGS84 <- st_transform(basonet_2001_sf, 4326)
coords_wgs84 <- st_coordinates(basonet_2001_WGS84)
Basonet_2001$Latitude <- coords_wgs84[,2]
Basonet_2001$Longitude <- coords_wgs84[,1]

### Transform to ETRS89
st_can_transform(23030,25830) # TRUE
basonet_2001_sf <- st_transform(basonet_2001_sf, 25830)

plot(basonet_2001_sf["TOC"])
tm_shape(world_eus, bbox=st_bbox(euskadi), crs=25830) +
  tm_borders() +
  tm_shape(autonomias25830) +
  tm_fill(col="gray95") +
  tm_borders() +
  tm_grid(lines=FALSE) +
  tm_shape(basonet_2001_sf) +
  tm_symbols(shape=21,
             size=0.3,
             legend.shape.show = TRUE,
             col="C_N",
             palette="viridis",
             border.col="black",
             title.shape = "BASONET")  +
  tm_layout(inner.margins = c(0,0,0,0),
            legend.outside=TRUE,
            legend.outside.position="right",
            legend.title.size=1.2,
            legend.text.size=0.8)

rm(basonet_2001_WGS84, coords_wgs84, basonet_2001_sf)

### Here for BASONET 2001 unless I find I need to correct something else


# 3.3 BASONET Bulk Density 2001 -------------------------------------------

setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/Basonet/")
### Join with Bulk Density from cylinder...
BasonetDens01 <- read_excel("BasonetDens01.xlsx", sheet = "DensAparCilindros2001")
colnames(BasonetDens01)[colnames(BasonetDens01)== "TH"] <- "th"
colnames(BasonetDens01)[colnames(BasonetDens01)== "Horiz"] <- "Lower_limit"
colnames(BasonetDens01)[colnames(BasonetDens01)== "Parcela"] <- "Plot"

### Average by horizon - three cores were sampled by horizon
BasonetDensAve <- BasonetDens01 %>%
  group_by(., th, Plot, Lower_limit) %>%
  summarise(., Bulk_Density_Cylinder = mean(DensApar, na.rm = TRUE))

### Join with all data from 2001
colnames(Basonet_2001)[colnames(Basonet_2001) %in% colnames(BasonetDensAve)]
Basonet_2001 <- left_join(Basonet_2001, BasonetDensAve, by= c("th","Plot","Lower_limit"))
BD_df <- Basonet_2001[,c("UTM_X","UTM_Y","th","Plot","Lower_limit",
                         "Bulk_Density_Cylinder","TOC","Clay","Silt","Sand")]

### In fact, I am going to replace the bulk density with the right bulk density values
plot(Basonet_2001$Bulk_Density, Basonet_2001$Bulk_Density_Cylinder)
abline(0,1)
Basonet_2001$Bulk_Density <- NA
Basonet_2001$Bulk_Density <- Basonet_2001$Bulk_Density_Cylinder
Basonet_2001$BD_method <- "Core method"
Basonet_2001$BD_quality <- "medium"

### BASONET 2001
sel.columns.basonet2001 <- c("Dataset","myID","UTM_X","UTM_Y","Longitude","Latitude",
                             "Upper_limit","Lower_limit","Date",
                             "Bulk_Density","BD_method","BD_quality",
                             "Sand","Silt","Clay","OM","TOC", ### Minimum dataset
                             "pH","TN_ppm","C_N",
                             "P_ppm","K_ppm","Ca_ppm","Mg_ppm",
                             "Carbonates","Act_carbonates", 
                             "P_HCl_ppm")
Basonet_2001_m <- Basonet_2001[,sel.columns.basonet2001]


# 4.1 LUCAS 2018 ----------------------------------------------------------

setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/LUCAS_2018")
library(readr)
BulkDensity_2018_final_2 <- read_csv("LUCAS-SOIL-2018-v2/BulkDensity_2018_final-2.csv")
LUCAS_SOIL_2018 <- read_csv("LUCAS-SOIL-2018-v2/LUCAS-SOIL-2018.csv")

### Attach the locations (real) from microdata
microdata_lucas2018 <- read_csv("LUCAS-SOIL-2018-v2/ES_2018_20200213.csv")

### Transform Bulk Density to long
colnames(BulkDensity_2018_final_2) <- c("POINTID","0-10 cm","10-20 cm","20-30 cm","0-20 cm")
BulkDensity_2018_long <- BulkDensity_2018_final_2 %>%
  as.data.frame(.,) %>%
  pivot_longer(-c(POINTID), names_to = "Depth", values_to = "Bulk_Density")

### Separate 20-30 cm soil data
lucas_20_30 <- LUCAS_SOIL_2018 %>%
  select(., c("POINTID", "OC (20-30 cm)", "CaCO3 (20-30 cm)")) %>%
  filter(., !is.na(`OC (20-30 cm)`)) %>%
  mutate(., "Depth" = "20-30 cm")
colnames(lucas_20_30) <- c("POINTID","OC","CaCO3","Depth")
lucas_20_30$OC <- as.character(lucas_20_30$OC)
lucas_20_30$CaCO3 <- as.character(lucas_20_30$CaCO3)
sel_cols <- colnames(LUCAS_SOIL_2018)[!colnames(LUCAS_SOIL_2018) %in% c("OC (20-30 cm)","CaCO3 (20-30 cm)")]
LUCAS_SOIL_2018 <- LUCAS_SOIL_2018[,sel_cols]

### I subset only data from Spain
LUCAS_SOIL_2018 <- LUCAS_SOIL_2018[LUCAS_SOIL_2018$NUTS_0 == "ES", ]
### subset 20-30 cm observations in the locations from Spain
lucas_20_30 <- lucas_20_30[lucas_20_30$POINTID %in% LUCAS_SOIL_2018$POINTID,]
### None present

### Join soil data and bulk density
LUCAS_2018 <- left_join(LUCAS_SOIL_2018, BulkDensity_2018_long)

### Change name
colnames(microdata_lucas2018)[colnames(microdata_lucas2018)=="POINT_ID"] <- "POINTID"
### What columns are the same?
colnames(microdata_lucas2018)[colnames(microdata_lucas2018) %in% colnames(LUCAS_SOIL_2018)]
### SURVEY_DATE has a different format in both dataframes so I don´t keep it
### Only columns I want to use
sel_microdata <- c("POINTID","GPS_LONG","GPS_LAT")
### Add real coordinates
microdata_lucas2018[,sel_microdata] 
### Transform POINTID to numeric
microdata_lucas2018$POINTID <- as.numeric(microdata_lucas2018$POINTID)
### Add real coordinates
LUCAS_2018 <- left_join(LUCAS_2018, microdata_lucas2018[,sel_microdata])

### Exclude points without real coordinates
LUCAS_2018 <- LUCAS_2018[LUCAS_2018$GPS_LAT != 88.888888,]
LUCAS_2018 <- LUCAS_2018[!is.na(LUCAS_2018$GPS_LAT),]

### Correct the Latitude to negative where it corresponds (comparison theoretical with recorded GPS coordinates)
LUCAS_2018$GPS_LONG <- ifelse(LUCAS_2018$TH_LONG >= 0, yes=LUCAS_2018$GPS_LONG, no = - LUCAS_2018$GPS_LONG)
LUCAS_2018$Longitude <- LUCAS_2018$GPS_LONG
LUCAS_2018$Latitude <- LUCAS_2018$GPS_LAT

### Assign CRS to LUCAS (theoretical coordinates. Real ones are in EUROSTAT microdata)
### Transform to spatial
LUCAS_2018_sf <- st_as_sf(LUCAS_2018, coords = c("GPS_LONG","GPS_LAT"), crs = 4326)
plot(LUCAS_2018_sf["Bulk_Density"])
LUCAS_2018_sf <- st_transform(LUCAS_2018_sf,25830)

### Intersection with Euskadi
LUCAS_2018_Eus <- st_intersection(LUCAS_2018_sf, euskadi)
### 38 observations

tm_shape(world_eus, bbox=st_bbox(euskadi), crs=25830) +
  tm_borders() +
  tm_shape(autonomias25830) +
  tm_fill(col="gray95") +
  tm_borders() +
  tm_grid(lines=FALSE) +
  tm_shape(LUCAS_2018_Eus) +
  tm_symbols(shape=21,
             size=0.3,
             legend.shape.show = TRUE,
             col="Bulk_Density",
             palette="viridis",
             border.col="black",
             title.shape = "BASONET")  +
  tm_layout(inner.margins = c(0,0,0,0),
            legend.outside=TRUE,
            legend.outside.position="right",
            legend.title.size=1.2,
            legend.text.size=0.8)

tm_shape(autonomias25830) +
  tm_fill(col="gray95") +
  tm_borders() +
   tm_grid(lines=FALSE) +
  tm_shape(LUCAS_2018_sf) +
  tm_symbols(shape=21,
             size=0.3,
             legend.shape.show = TRUE,
             col="Bulk_Density",
             palette="viridis",
             border.col="black",
             title.shape = "BASONET") +
  tm_layout(inner.margins = c(0,0,0,0),
            legend.outside=TRUE,
            legend.outside.position="right",
            legend.title.size=1.2,
            legend.text.size=0.8)

### Drop the geometry. Back to dataframe
LUCAS_2018_Eus_df <- LUCAS_2018_Eus
coords_lucas_utm <- st_coordinates(LUCAS_2018_Eus)
LUCAS_2018_Eus_df$UTM_X <- coords_lucas_utm[,1]
LUCAS_2018_Eus_df$UTM_Y <- coords_lucas_utm[,2]
LUCAS_2018_Eus_df <- as.data.frame(LUCAS_2018_Eus_df)

### Add info on sampling year in a new column "Date"
LUCAS_2018_Eus_df$SURVEY_DATE
LUCAS_2018_Eus_df$Date <- 2018

### LUCAS 2018 variables of interest
# "Bulk_Density" "pH_CaCl2"     "pH_H2O"       "EC"     "OC"          
# "CaCO3"        "P"            "N"            "K"      "Ox_Al"        "Ox_Fe"

### bulk density method and quality is the highest because the protocol is very detailed
### Units g cm-3
LUCAS_2018_Eus_df$BD_method <- "Core method"
LUCAS_2018_Eus_df$BD_quality <- "high"

### Change variable units when needed in LUCAS
# pH_CaCl2 pH – measured in calcium chloride
# pH_H2O pH – measured in water - same as in CARBOSOL and BASONET
# Glass electrode in a 1:5 (V/V) suspension of soil in H20 and CaCl2
colnames(LUCAS_2018_Eus_df)[colnames(LUCAS_2018_Eus_df)=="pH_H2O"] <- "pH"

### EC (mS m-1) ISO 11265:1994 Metal electrodes in aqueous extract of soil
### EC Electrical conductivity (milliSiemens per meter – mS m-1)
summary(carbosol_Eus_df$EC)
summary(LUCAS_2018_Eus_df$EC)
# In Carbosol dS/m so I divide by 100
LUCAS_2018_Eus_df$EC <- LUCAS_2018_Eus_df$EC / 100

# "OC" Organic carbon content (g kg-1) ISO 10694:1995 Dry combustion (elementary analysis)
# and in Carbosol "TOC" it is in %
LUCAS_2018_Eus_df$TOC <- as.numeric(LUCAS_2018_Eus_df$OC) / 10

# "CaCO3" Calcium carbonate content (g kg-1)
# ISO 10693:1995 Volumetric method
# LOD is 1
# Transform to % as in Carbosol
LUCAS_2018_Eus_df$Carbonates <- as.numeric(LUCAS_2018_Eus_df$CaCO3)/ 10
LUCAS_2018_Eus_df[is.na(LUCAS_2018_Eus_df$Carbonates),]$pH
### At that pH carbonates are likely 0
LUCAS_2018_Eus_df[is.na(LUCAS_2018_Eus_df$Carbonates),]$Carbonates <- 0

### LUCAS "N" Modified Kjeldahl method Total nitrogen (g kg-1)
##  "TN_ppm" N ppm Kjeldahl Method or Dry Combustion Methods (elemental analyzer)
LUCAS_2018_Eus_df$TN_ppm <- as.numeric(LUCAS_2018_Eus_df$N) * 1000

# "P" Spectrometric determination of P soluble in sodium hydrogen CaCO3 solution
# ISO 11263:1994  Spectrometric determination of P soluble in sodium hydrogen CaCO3 solution
# Total (available) phosphorus (mg kg-1)
LUCAS_2018_Eus_df$P_ppm <- as.numeric(LUCAS_2018_Eus_df$P)
summary(LUCAS_2018_Eus_df$P_ppm)
summary(carbosol_Eus_df$P_ppm)

# "K" Atomic absorption spectrometry after extraction with NH4OAc
# Extractable potassium (mg kg-1)
# Carbosol "K_ppm" Potasium ppm	Extraction with ammonium acetate and
# quantification by ICP spectroscopy or Flame Photometry Methods
LUCAS_2018_Eus_df$K_ppm <- as.numeric(LUCAS_2018_Eus_df$K)

# "C_N"  Carbon to nitrogen ratio
LUCAS_2018_Eus_df$C_N <- as.numeric(LUCAS_2018_Eus_df$OC) / as.numeric(LUCAS_2018_Eus_df$N)
LUCAS_2018_Eus_df$C_N

# Oxalate extractable Fe and Al - Ross and Wang, (1993) - Acid ammonium oxalate method
LUCAS_2018_Eus_df$Ox_Al
LUCAS_2018_Eus_df$Ox_Fe

### Are measurements a bit better now?
summary(carbosol_Eus_df$TN_ppm)
summary(LUCAS_2018_Eus_df$TN_ppm)
summary(Basonet_2001$TN_ppm)
hist(LUCAS_2018_Eus_df$TN_ppm)
hist(carbosol_Eus_df$TN_ppm)
hist(Basonet_2001$TN_ppm)
summary(carbosol_Eus_df$P_ppm)
summary(LUCAS_2018_Eus_df$P_ppm)
summary(LUCAS_2018_Eus_df$K_ppm)
summary(carbosol_Eus_df$K_ppm)
summary(Basonet_2001$K_ppm)

### LUCAS 2018
colnames(LUCAS_2018_Eus_df)

### Add missing variables
LUCAS_2018_Eus_df$Dataset <- "LUCAS"
LUCAS_2018_Eus_df$myID <- paste0("LUCAS_",LUCAS_2018_Eus_df$POINTID)
table(LUCAS_2018_Eus_df$Depth)
LUCAS_2018_Eus_df$Upper_limit <- 0
LUCAS_2018_Eus_df$Lower_limit <- 20

sel.columns.lucas <- c("Dataset","myID","UTM_X","UTM_Y","Longitude","Latitude",
                       "Upper_limit","Lower_limit","Date", #"  for LUCAS
                       "Bulk_Density","BD_method","BD_quality",
                       "pH","EC","TOC","Carbonates","TN_ppm", "C_N",
                       "P_ppm", "K_ppm", ## If they exist
                       "pH_CaCl2", "Ox_Al","Ox_Fe") # for LUCAS 2018
LUCAS_2018_m <- LUCAS_2018_Eus_df[,sel.columns.lucas]


# 4.2 LUCAS 2015 ----------------------------------------------------------

setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/")
library(readr)
library(readxl)
LUCAS_2015 <- read_excel("LUCAS_2015/LUCAS_Topsoil_2015_20200323.xlsx")
LUCAS_2015_microdata <- read_csv("LUCAS_2015/ES_2015_20200225.CSV")
LUCAS_2015_microdata$GPS_LONG
summary(LUCAS_2015_microdata$GPS_LONG)
summary(LUCAS_2015_microdata$TH_LONG)

### Correct the Latitude to negative where it corresponds (comparison theoretical with recorded GPS coordinates)
LUCAS_2015_microdata <- LUCAS_2015_microdata[LUCAS_2015_microdata$GPS_LONG < 88,]
### Eliminate sites without real coordinates
LUCAS_2015_microdata <- LUCAS_2015_microdata[!is.na(LUCAS_2015_microdata$GPS_LAT),]
summary(LUCAS_2015_microdata$TH_LONG)
summary(LUCAS_2015_microdata$GPS_LONG)

### Correct the Latitude to negative where it corresponds (comparison theoretical with recorded GPS coordinates)
LUCAS_2015_microdata$GPS_LONG <- ifelse(LUCAS_2015_microdata$TH_LONG >= 0,
                                        yes=LUCAS_2015_microdata$GPS_LONG,
                                        no = - LUCAS_2015_microdata$GPS_LONG)
### Change name
colnames(LUCAS_2015_microdata)[colnames(LUCAS_2015_microdata)=="POINT_ID"] <- "POINTID"
colnames(LUCAS_2015)[colnames(LUCAS_2015)=="Point_ID"] <- "POINTID"
### What columns are the same?
colnames(LUCAS_2015_microdata)[colnames(LUCAS_2015_microdata) %in% colnames(LUCAS_2015)]
### Only columns I want to use
sel_microdata <- c("POINTID","GPS_LONG","GPS_LAT")
LUCAS_2015[LUCAS_2015$NUTS_0 == "ES", ]
### Transform POINTID to numeric
LUCAS_2015_microdata$POINTID <- as.numeric(LUCAS_2015_microdata$POINTID)
### Add coordinates
LUCAS_2015 <- inner_join(LUCAS_2015, LUCAS_2015_microdata[,sel_microdata], by=c("POINTID"))
summary(LUCAS_2015$GPS_LAT)
dim(LUCAS_2015[!is.na(LUCAS_2015$GPS_LAT),])

### Latitude and Longitude
LUCAS_2015$Latitude <- LUCAS_2015$GPS_LAT
LUCAS_2015$Longitude <- LUCAS_2015$GPS_LONG

### Transform to spatial and project
LUCAS_2015_sf <- st_as_sf(LUCAS_2015, coords = c("GPS_LONG","GPS_LAT"), crs = 4326)
colnames(LUCAS_2015)
plot(LUCAS_2015_sf["pH(H2O)"])
st_crs(LUCAS_2015_sf)
st_crs(euskadi)

### Project
LUCAS_2015_sf <- st_transform(LUCAS_2015_sf,25830)
### Intersection with Euskadi
LUCAS_2015_Eus <- st_intersection(LUCAS_2015_sf, euskadi)
### Only 43 points

tm_shape(world_eus, bbox=st_bbox(euskadi), crs=25830) +
  tm_borders() +
  tm_shape(autonomias25830) +
  tm_fill(col="gray95") +
  tm_borders() +
  tm_grid(lines=FALSE) +
  tm_shape(LUCAS_2015_Eus) +
  tm_symbols(shape=21,
             size=0.3,
             legend.shape.show = TRUE,
             col="pH.H2O.",
             palette="viridis",
             border.col="black",
             title.shape = "LUCAS 2015")  +
  tm_layout(inner.margins = c(0,0,0,0),
                          legend.outside=TRUE,
                          legend.outside.position="right",
                          legend.title.size=1.2,
                          legend.text.size=0.8)

### Now that I have the subset of data for Euskadi
### Drop the geometry. Back to dataframe
LUCAS_2015_Eus_df <- LUCAS_2015_Eus
st_geometry(LUCAS_2015_Eus_df) <- NULL
coords_lucas_utm <- st_coordinates(LUCAS_2015_Eus)
LUCAS_2015_Eus_df$UTM_X <- coords_lucas_utm[,1]
LUCAS_2015_Eus_df$UTM_Y <- coords_lucas_utm[,2]
LUCAS_2015_Eus_df <- as.data.frame(LUCAS_2015_Eus_df)

### Add info on sampling year in a new column "Date"
LUCAS_2015_Eus_df$Date <- 2015 ### All in 2015

### soil variables
# pH(CaCl2)	pH measured in a CaCl2 solution	─
# pH(H2O)	pH measured in a suspension of soil in water	─
# EC	Electrical conductivity	mS/m
# OC	Organic carbon content	g/kg
# N	Total nitrogen content	g/kg
# CaCO3	Carbonates content	g/kg
# P	Phosphorus content	mg/kg (ppm)
# K	Extractable potassium content	mg/kg (ppm)

### Change variable units or variable names when needed in LUCAS
# pH_CaCl2 pH – measured in calcium chloride
# pH_H2O pH – measured in water - same as in CARBOSOL and BASONET
### the ratio is not 1:2.5 but
# Glass electrode in a 1:5 (V/V) suspension of soil in H20 and CaCl2
colnames(LUCAS_2015_Eus_df)[colnames(LUCAS_2015_Eus_df) == "pH.H2O."] <- "pH"
colnames(LUCAS_2015_Eus_df)[colnames(LUCAS_2015_Eus_df) == "pH.CaCl2."] <- "pH_CaCl2"
colnames(LUCAS_2015_Eus_df)[colnames(LUCAS_2015_Eus_df) == "Coarse"] <- "CoarseFrag_Grav"

# EC Electrical conductivity (milliSiemens per meter – mS m-1)
# measured Metal electrodes in aqueous extract of soil (ISO 11265:1994)
# In Carbosol dS/m so I divide by 100
LUCAS_2015_Eus_df$EC <- LUCAS_2015_Eus_df$EC / 100

# "OC" Organic carbon content (g kg-1)
# and in Carbosol "TOC" it is in %
LUCAS_2015_Eus_df$TOC <- LUCAS_2015_Eus_df$OC / 10

# "CaCO3" Calcium carbonate content (g kg-1)
# "Carbonates"  Carbonates (%)
LUCAS_2015_Eus_df$Carbonates <- LUCAS_2015_Eus_df$CaCO3/ 10

### LUCAS "N" Modified Kjeldahl method Total nitrogen (g kg-1)
##  "TN_ppm" N ppm Kjeldahl Method or Dry Combustion Methods (elemental analyzer)
LUCAS_2015_Eus_df$TN_ppm <- as.numeric(LUCAS_2015_Eus_df$N) * 1000

# "P" Spectrometric determination of P soluble in sodium hydrogen CaCO3 solution
#  Total phosphorus mg/kg - this is already in ppm - ISO 11263:1194
# "P_ppm" Phosphorus ppm	Olsen or Mehlich Methods
LUCAS_2015_Eus_df$P_ppm <- LUCAS_2015_Eus_df$P

# "K" Atomic absorption spectrometry after extraction with NH4OAc
# Extractable potassium (mg/kg)
# Carbosol "K_ppm" Potasium ppm	Extraction with ammonium acetate and
# quantification by ICP spectroscopy or Flame Photometry Methods
LUCAS_2015_Eus_df$K_ppm <- LUCAS_2015_Eus_df$K

# "C_N"  Carbon to nitrogen ratio
LUCAS_2015_Eus_df$C_N <- LUCAS_2015_Eus_df$OC / LUCAS_2015_Eus_df$N

LUCAS_2015_Eus_df$Dataset <- "LUCAS"
LUCAS_2015_Eus_df$myID <- paste0("LUCAS_",LUCAS_2015_Eus_df$POINTID)
### Sample depth in this campaign was 0-20 cm
### See
### Jones, A, Fernandez-Ugalde, O., Scarpa, S. LUCAS 2015 Topsoil Survey.
### Presentation of dataset and results, EUR 30332 EN, Publications Office of the European Union:
### Luxembourg. 2020, ISBN 978-92-76-21080-1, doi:10.2760/616084, JRC121325
LUCAS_2015_Eus_df$Upper_limit <- 0
LUCAS_2015_Eus_df$Lower_limit <- 20

### where these points revisited from 2009?
LUCAS_2015_Eus_df$Revisited_point
### Extract site ID from the revisited points
lucas_2015_revisited <- LUCAS_2015_Eus_df[LUCAS_2015_Eus_df$Revisited_point =="Yes",]
lucas_2015_revisited <- lucas_2015_revisited[,c("POINTID","Revisited_point","UTM_X","UTM_Y","Longitude","Latitude")]
dim(lucas_2015_revisited) ### 40 points

### Transform to spatial and save
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/LUCAS_2015/")
LUCAS_2015_Eus <- st_as_sf(LUCAS_2015_Eus_df, coords = c("UTM_X","UTM_Y"), crs = 25830)
st_write(obj=LUCAS_2015_Eus, dsn="LUCAS_2015_Eus.shp", append=FALSE)
pointID_lucas2015 <- LUCAS_2015_Eus_df$POINTID


## LUCAS 2015
colnames(LUCAS_2015_Eus_df)
sel.columns.lucas <- c("Dataset","myID","UTM_X","UTM_Y","Longitude","Latitude",
                       "Upper_limit","Lower_limit","Date",
                       "CoarseFrag_Grav", #"Clay","Sand","Silt", I exclude clay, silt, and sand because there were only three measurmeents and they were obtained by laser diffraction
                       "pH","EC","TOC","Carbonates","TN_ppm","P_ppm","K_ppm","C_N", ## If they exist
                       "pH_CaCl2") # for LUCAS to compare with ppm values from Carbosol...
setdiff(sel.columns.lucas,colnames(LUCAS_2015_Eus_df))
LUCAS_2015_m <- LUCAS_2015_Eus_df[,sel.columns.lucas]


# 4.3 LUCAS 2009 ----------------------------------------------------------

library(readr)
LUCAS_2009_locs <- st_read("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/LUCAS_2009/LUCAS_21681_points/LUCAS_Points.shp")
NoSensitive <- read_delim("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/LUCAS_2009/LUCAS_21681_points/NoSensitive.csv",
                          delim = "|", escape_double = FALSE, trim_ws = TRUE)

### Keep this
Lucas2009_revisited_Eus <- NoSensitive[NoSensitive$POINT_ID %in% lucas_2015_revisited$POINTID,]
Lucas2009_pointID2015_Eus <- NoSensitive[NoSensitive$POINT_ID %in% pointID_lucas2015,]

### Change names
colnames(LUCAS_2009_locs)[colnames(LUCAS_2009_locs)=="POINT_ID"] <- "POINTID"
colnames(NoSensitive)[colnames(NoSensitive)=="POINT_ID"] <- "POINTID"

### What columns are the same?
colnames(NoSensitive)[colnames(NoSensitive) %in% colnames(LUCAS_2009_locs)]

### Join
LUCAS_2009 <- left_join(LUCAS_2009_locs, NoSensitive)

### Project
LUCAS_2009_sf <- st_transform(LUCAS_2009,25830)
### Intersection with Euskadi
LUCAS_2009_Eus <- st_intersection(LUCAS_2009_sf, euskadi)
colnames(LUCAS_2009_Eus)
### Only 43 points fall within the Basque Country, but of these, of only 40 we know the real coordinates.
### (I compared in ArcGIS LUCAS_2015_Eus with LUCAS_2009_sf and they are close but not in the exact coordinates)

### Keep the theoretical coordinates of the three that are not revisited
LUCAS_2009_Eus_nR <- LUCAS_2009_Eus[!LUCAS_2009_Eus$POINTID %in% pointID_lucas2015, ]
LUCAS_2009_Eus_nR$Revisited_point <- "No"

### Drop the geometry. Back to dataframe
LUCAS_2009_Eus_nR_df <- LUCAS_2009_Eus_nR

st_geometry(LUCAS_2009_Eus_nR_df) <- NULL
coords_lucas2009nR_utm <- st_coordinates(LUCAS_2009_Eus_nR)
LUCAS_2009_Eus_nR_df$UTM_X <- coords_lucas2009nR_utm[,1]
LUCAS_2009_Eus_nR_df$UTM_Y <- coords_lucas2009nR_utm[,2]
LUCAS_2009_Eus_nR_df <- as.data.frame(LUCAS_2009_Eus_nR_df)
LUCAS_2009_Eus_nR_df$Latitude <- LUCAS_2009_Eus_nR_df$GPS_LAT
LUCAS_2009_Eus_nR_df$Longitude <- LUCAS_2009_Eus_nR_df$GPS_LONG

### Add real coordinates to the other ones
colnames(Lucas2009_revisited_Eus)[colnames(Lucas2009_revisited_Eus)=="POINT_ID"] <- "POINTID"
Lucas2009_revisited_Eus <- left_join(Lucas2009_revisited_Eus,lucas_2015_revisited)
Lucas2009_revisited_Eus$Latitude <- Lucas2009_revisited_Eus$Latitude
Lucas2009_revisited_Eus$Longitude <- Lucas2009_revisited_Eus$Longitude

setdiff(colnames(LUCAS_2009_Eus_nR_df), colnames(Lucas2009_revisited_Eus))
setdiff(colnames(Lucas2009_revisited_Eus),colnames(LUCAS_2009_Eus_nR_df))
colnames(LUCAS_2009_Eus_nR_df)[colnames(LUCAS_2009_Eus_nR_df)=="LC.parent"] <- "LC-parent"
colnames(LUCAS_2009_Eus_nR_df)[colnames(LUCAS_2009_Eus_nR_df)=="LC.Group"] <- "LC-Group"
colnames(LUCAS_2009_Eus_nR_df)[colnames(LUCAS_2009_Eus_nR_df)=="Elevation.m."] <- "Elevation(m)"
colnames(LUCAS_2009_Eus_nR_df)[colnames(LUCAS_2009_Eus_nR_df)=="slope..Degrees."] <- "slope (Degrees)"

LUCAS_2009_Eus_df <- bind_rows(Lucas2009_revisited_Eus,LUCAS_2009_Eus_nR_df)
View(LUCAS_2009_Eus_df)

### clean
rm(LUCAS_2009_locs,NoSensitive,Lucas2009_revisited_Eus,LUCAS_2009_Eus_nR_df,LUCAS_2009_Eus_nR,
   LUCAS_2009,LUCAS_2009_sf,LUCAS_2009_Eus)

### Add info on sampling year in a new column "Date"
LUCAS_2009_Eus_df$SURV_DATE ### All in 2009
LUCAS_2009_Eus_df$Date <- 2009

### which P is extractable P?
summary(LUCAS_2009_Eus_df$P_x)
summary(LUCAS_2009_Eus_df$PTotal)

SoilAttr_LUCAS2009 <- st_read("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/LUCAS_2009/package-for-ESDAC-20190927/SoilAttr_LUCAS2009/SoilAttr_LUCAS2009.shp")
### Project
SoilAttr_LUCAS2009 <- st_transform(SoilAttr_LUCAS2009,25830)
### Intersection with Euskadi
SoilAttr_LUCAS2009 <- st_intersection(SoilAttr_LUCAS2009, euskadi)
SoilAttr_LUCAS2009
summary(as.numeric(SoilAttr_LUCAS2009$P))
summary(LUCAS_2009_Eus_df$P_x)### This is it
summary(LUCAS_2009_Eus_df$PTotal)
summary(as.numeric(SoilAttr_LUCAS2009$N))
summary(LUCAS_2009_Eus_df$N) ### This is it

colnames(LUCAS_2009_Eus_df)[colnames(LUCAS_2009_Eus_df) == "SURV_DATE"] <- "SURVEY_DATE"

### Eliminate some columns
LUCAS_2009_Eus_df <- LUCAS_2009_Eus_df[,c("POINTID","coarse","clay","silt","sand",
                                          "pH_in_H2O","pH_in_CaCl","OC","CaCO3",
                                          "N","P_x","K","CEC",#"Revisited_point","SURVEY_DATE",
                                          "UTM_X","UTM_Y","Date","Longitude","Latitude")]

### Sample depth in this campaign was 0-20 cm
### See
### Jones, A, Fernandez-Ugalde, O., Scarpa, S. LUCAS 2009 Topsoil Survey.
### Presentation of dataset and results, EUR 30332 EN, Publications Office of the European Union:
### Luxembourg. 2020, ISBN 978-92-76-21080-1, doi:10.2760/616084, JRC121325
LUCAS_2009_Eus_df$Upper_limit <- 0
LUCAS_2009_Eus_df$Lower_limit <- 20
### Indicate dataset
LUCAS_2009_Eus_df$Dataset <- "LUCAS"
LUCAS_2009_Eus_df$myID <- paste0("LUCAS_",LUCAS_2009_Eus_df$POINTID)

### Soil variables - rename those that don´t need transformation
# coarse -  coarse fragments in % - ISO 11464. 2006
# clay -  clay content (%) - ISO 11277. 1998 - sieving and sedimentation
# silt -  silt content (%) - ISO 11277. 1998 - sieving and sedimentation
# sand -  sand content (%) - ISO 11277. 1998 - sieving and sedimentation
# pH_in_H2O -  pH measured from water solution - ISO 10390. 1994 
# pH_in_CaCl -  pH measured from CaCl solution - ISO 10390. 1994
# routine determination of pH using a glass electrode in a 1:5 (V/V) suspension
# of soil in water (pH-H2O), in a solution of 0,01 mol/l calcium chloride (pH-CaCl2).

### Rename some variables
colnames(LUCAS_2009_Eus_df)[colnames(LUCAS_2009_Eus_df) == "coarse"] <- "CoarseFrag_Grav" 
### I assume it is gravimetric, but I am not sure
colnames(LUCAS_2009_Eus_df)[colnames(LUCAS_2009_Eus_df) == "clay"] <- "Clay"
colnames(LUCAS_2009_Eus_df)[colnames(LUCAS_2009_Eus_df) == "silt"] <- "Silt"
colnames(LUCAS_2009_Eus_df)[colnames(LUCAS_2009_Eus_df) == "sand"] <- "Sand"
colnames(LUCAS_2009_Eus_df)[colnames(LUCAS_2009_Eus_df) == "pH_in_H2O"] <- "pH"
colnames(LUCAS_2009_Eus_df)[colnames(LUCAS_2009_Eus_df) == "pH_in_CaCl"] <- "pH_CaCl2"

# OC -  organic carbon content (g/kg) - SO 10694. 1995
LUCAS_2009_Eus_df$TOC <- LUCAS_2009_Eus_df$OC / 10

# CaCO3 -  CaCO3 content (g/kg) - ISO 10693. 1994
LUCAS_2009_Eus_df$Carbonates <- LUCAS_2009_Eus_df$CaCO3/ 10

# N -  Nitrogen content (g/kg) - ISO 11261. 1995
# Determination of total nitrogen - Modified Kjeldahl method
LUCAS_2009_Eus_df$TN_ppm <- LUCAS_2009_Eus_df$N * 1000

# P -  Phosphorus content (mg/kg) - ISO 11263. 1994
### Spectrometric determination of phosphorus soluble in sodium hydrogen carbonate solution
colnames(LUCAS_2009_Eus_df)[colnames(LUCAS_2009_Eus_df) == "P_x"] <- "P_ppm"

# K -  Extractable Potassium content (mg/kg) - USDA, 2004
colnames(LUCAS_2009_Eus_df)[colnames(LUCAS_2009_Eus_df) == "K"] <- "K_ppm"

# CEC -  Cation Exchange Capacity (cmol(+)/kg)
LUCAS_2009_Eus_df$CEC ### as it is

# "C_N"  Carbon to nitrogen ratio
LUCAS_2009_Eus_df$C_N <- LUCAS_2009_Eus_df$OC / LUCAS_2009_Eus_df$N

## LUCAS 2009 
colnames(LUCAS_2009_Eus_df)
sel.columns.lucas <- c("Dataset","myID","UTM_X","UTM_Y","Longitude","Latitude",
                       "Upper_limit","Lower_limit","Date",
                       "CoarseFrag_Grav","Clay","Sand","Silt",
                       "pH","TOC","Carbonates","CEC","TN_ppm","P_ppm","K_ppm","C_N", ## If they exist
                       "pH_CaCl2") 
LUCAS_2009_m <- LUCAS_2009_Eus_df[,sel.columns.lucas]

## Clay, silt and sand for LUCAS follow the USDA limits
## Sand: 50 – 2000 µm
## Silt: 2 –50 µm
## Clay: < 2 µm

# # preparing LUCAS soil texture for TT.text.transform function
# library(soiltexture)
# lucas.t <- LUCAS_2009_Eus_df[,c("Sand","Silt","Clay")] %>%
#   mutate(SAND = Sand,
#          SILT = Silt,
#          CLAY = Clay)
# lucas.t <- lucas.t %>%
#   rowwise() %>%
#   drop_na(c(SAND, SILT, CLAY)) ## no NA allowed in the transformation
# 
# lucas.txt <- TT.normalise.sum(as.data.frame(select(lucas.t, c("SAND","SILT","CLAY")))) %>%
#   bind_cols(LUCAS_2009_Eus_df) %>%
#   drop_na(c("SAND","SILT","CLAY"))
# lucas.txt <- TT.normalise.sum(as.data.frame(select(lucas.t, c("SAND","SILT","CLAY")))) %>%
#   bind_cols(LUCAS_2009_Eus_df) %>%
#   drop_na(c("SAND","SILT","CLAY"))
# 
# lucas.trans <- TT.text.transf(
#   tri.data = lucas.txt,
#   base.css.ps.lim = c(0,2,50,2000),
#   dat.css.ps.lim = c(0,2,63,2000),
#   text.tol = 1/100
#   ) ## tolerance for the sum of clay, silt and sand of 1%
# lucas.trans
# ### Replace in the original dataframe
# LUCAS_2009_Eus_df$Clay <- lucas.trans$CLAY
# LUCAS_2009_Eus_df$Silt <- lucas.trans$SILT
# LUCAS_2009_Eus_df$Sand <- lucas.trans$SAND

# 5. Merge datasets ----------------------------------------------------------

### Merge
soil_datasets <- bind_rows(carbosol_m, ines_m)
soil_datasets <- bind_rows(soil_datasets,Basonet_2001_m)
soil_datasets <- bind_rows(soil_datasets,Basonet_2021_m)
soil_datasets <- bind_rows(soil_datasets,LUCAS_2018_m)
soil_datasets <- bind_rows(soil_datasets,LUCAS_2015_m)
soil_datasets <- bind_rows(soil_datasets,LUCAS_2009_m)

### 2965 observations at 1109 unique locations
dim(soil_datasets)
length(unique(soil_datasets$myID))

### After talking with Alejandro Cantero (HAZI) and checking differences in BASONET plot locations afgter transfromation from epsg 23030 yp epsg 25830
### It seems better to assign coordinates from 2021 to those plots class A1, and keep coordinates from 2001 to those that could not be completely revisited.
BasonetParcelas2001 <- read_excel("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Stand_data_Basonet/BasonetDasom.xlsx",
                                  sheet = "parcelas2001")
BasonetParcelas2021 <- read_excel("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Stand_data_Basonet/BasonetDasom.xlsx",
                                  sheet = "parcelas2021222")
head(BasonetParcelas2001)

### Create my ID
BasonetParcelas2001$myID <- paste0("BASONET_",BasonetParcelas2001$th,"_",BasonetParcelas2001$parc)
BasonetParcelas2021$myID <- paste0("BASONET_",BasonetParcelas2021$th,"_",BasonetParcelas2021$parc)

### subset A1 plots 
basonet2001_a1 <- BasonetParcelas2001[BasonetParcelas2001$Clase == "A1",]
basonet2021_a1 <- BasonetParcelas2021[BasonetParcelas2021$Clase == "A1",]

susbtitute.coords.ids <- basonet2021_a1$myID

all.equal(basonet2001_a1$myID,basonet2021_a1$myID )
setdiff(basonet2021_a1$myID,basonet2001_a1$myID)
setdiff(basonet2001_a1$myID,basonet2021_a1$myID)

### How are plots that were A1 in 2001 but not in 2021?
what.about <- setdiff(basonet2001_a1$myID, basonet2021_a1$myID)
BasonetParcelas2021[BasonetParcelas2021$myID %in% what.about, ]$Clase ### were A4

### And the opposite?
BasonetParcelas2001[BasonetParcelas2001$myID %in% 
                      setdiff(basonet2021_a1$myID, basonet2001_a1$myID), ]$Clase ### were A2 and A4
### what does that mean? I don´t know what A2 or A4 plots mean

soil_datasets[soil_datasets$Dataset == "BASONET" & soil_datasets$myID == susbtitute.coords.ids[[1]], ]


### anyway, I go with those A1 for 2021
for(i in 1:length(susbtitute.coords.ids)){
  
  print(i)
  
  correct_UTM_X <- unique(soil_datasets[soil_datasets$Dataset == "BASONET" &
                                   soil_datasets$Date == 2021 &
                                   soil_datasets$myID == susbtitute.coords.ids[[i]], ]$UTM_X)
  
  correct_UTM_Y <- unique(soil_datasets[soil_datasets$Dataset == "BASONET" &
                                   soil_datasets$Date == 2021 &
                                   soil_datasets$myID == susbtitute.coords.ids[[i]], ]$UTM_Y)
  
  target.length <- length(soil_datasets[soil_datasets$Dataset == "BASONET" &
                                          soil_datasets$Date == 2001 &
                                          soil_datasets$myID == susbtitute.coords.ids[[i]], ]$UTM_X)
  
  if (length(correct_UTM_X) == 1 & length(target.length) > 0) {
    
    print("Substitute 2001 coordinates with 2021 coordinates")

    ### Correct values
    soil_datasets[soil_datasets$Dataset == "BASONET" &
                    soil_datasets$Date == 2001 &
                    soil_datasets$myID == susbtitute.coords.ids[[i]], ]$UTM_X <- c(rep(correct_UTM_X, target.length))
    soil_datasets[soil_datasets$Dataset == "BASONET" &
                    soil_datasets$Date == 2001 &
                    soil_datasets$myID == susbtitute.coords.ids[[i]], ]$UTM_Y <- c(rep(correct_UTM_Y,target.length))
    
  } else if (length(correct_UTM_X) == 1 & length(target.length) == 0) {
    
    print("2021 was not sampled in 2001")
    
    } else {
      
      print("there is no soil data for this plot")
    
    }
}

### There are some exceptions, cases 23, 52 and 146
soil_datasets[soil_datasets$Dataset == "BASONET" &
                soil_datasets$Date == 2021 &
                soil_datasets$myID == susbtitute.coords.ids[[23]], ]
soil_datasets[soil_datasets$Dataset == "BASONET" &
                soil_datasets$Date == 2021 &
                soil_datasets$myID == susbtitute.coords.ids[[52]], ]
soil_datasets[soil_datasets$Dataset == "BASONET" &
                soil_datasets$Date == 2021 &
                soil_datasets$myID == susbtitute.coords.ids[[146]], ]

### Make spatial
soil_datasets_sf <- st_as_sf(soil_datasets, coords = c("UTM_X","UTM_Y"), crs = 25830)

tm_shape(world_eus, bbox=st_bbox(euskadi), crs=25830) +
  tm_borders() +
  tm_shape(autonomias25830) +
 # tm_fill(col="gray95") +
  tm_borders() +
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_scalebar(breaks = c(0,20,40), text.size = 0.8, position = c("right", "bottom")) +
  tm_grid(lines=FALSE) +
  tm_shape(soil_datasets_sf[soil_datasets_sf$TOC <=5,]) +
  tm_symbols(shape=21,
             size=0.3,
             legend.shape.show = TRUE,
             col="Dataset",
             palette="magma",
             border.col="black",
             title.shape = "Soil dataset")  +
  tm_layout(inner.margins = c(0,0,0,0),
            legend.outside=TRUE,
            legend.outside.position="right",
            legend.title.size=1.2,
            legend.text.size=0.8)

### how many unique locations?
length(unique(soil_datasets_sf$myID)) #1109

### How many unique IDs based on coordinates?
soil_datasets_sf$ID_coords <- paste0(soil_datasets$UTM_X, "_", soil_datasets$UTM_Y)
soil_datasets$ID_coords <- paste0(soil_datasets$UTM_X, "_", soil_datasets$UTM_Y)

length(unique(soil_datasets_sf$ID_coords)) ## 1306, even more places, but some BASONET plots have been corrected
length(unique(soil_datasets_sf$myID))
length(unique(soil_datasets$ID_coords))

BASONET <- soil_datasets_sf[soil_datasets_sf$Dataset =="BASONET",]
length(unique(BASONET$ID_coords)) # 647 locations
length(unique(BASONET$myID)) # 482 plots

### I will start analyse the change in soil properties with forest management in those revisited plots (198)
BASONET$Revisited <- ifelse(BASONET$myID %in% susbtitute.coords.ids,  "yes",  "no" )
susbtitute.coords.ids

### Export to RData
save(world_eus,soil_datasets,soil_datasets_sf,BASONET,
     euskadi,eus_buffer,eus_buffer_WGS84,autonomias,autonomias25830,
     file="C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/R_Output/soildatasets_correctBASONETcoords.RData")
save(world_eus,soil_datasets,soil_datasets_sf,BASONET,
     euskadi,eus_buffer,eus_buffer_WGS84,autonomias,autonomias25830,
     file="C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/RMarkdown/soildatasets_correctBASONETcoords.RData")

### subset one point per unique location
# uniqueIDs <- unique(soil_datasets_sf$myID)
# length(uniqueIDs) # 1109 sites
# Euskadi_sites <- soil_datasets_sf[,c("Dataset","myID","Longitude","Latitude","Date")]
# Euskadi_sites <- Euskadi_sites[!duplicated(Euskadi_sites),] ### 1537 sites
# Euskadi_sites[Euskadi_sites$Dataset=="BASONET",]

### Write to shapefile
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/RMarkdown")
st_write(obj=soil_datasets_sf, dsn="soil_datasets.shp", append=FALSE)
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/R_Output")
st_write(obj=soil_datasets_sf, dsn="soil_datasets.shp", append=FALSE)

### I want to use BASONET for TymeSync and for SEPAL CCDC tools
### For TimeSync:
# Make a point or polygon shapefile for areas of interest. 
# You will need to add at least two fields to your shapefile called "uniqID" and "yod", 
# where each feature has a unique ID (string), and a year of detection which can be null (int).
# The unique IDs will become the display names for each feature in the TimeSyncPlus application. 
# The year of detection will also be display for each feature in the aplication. 
# This shapefile will also need to be zipped in order to upload it to Google Earth Engine.

BASONET_sites <- BASONET[,c("Dataset","myID","Date","Revisited", "ID_coords")]
BASONET_sites <- BASONET_sites[!duplicated(BASONET_sites),] ### 833 sites
BASONET_sites[BASONET_sites$Dataset=="BASONET",]

BASONET_sites <- BASONET_sites %>%
  arrange(., Dataset, ID_coords, -Date) %>%
  group_by(., ID_coords) %>% 
  mutate(., newID =  cur_group_id())

### Eliminate 2001 plots in those revisited
BASONET_sites <- BASONET_sites[!duplicated(BASONET_sites$newID),]

### VARIABLES for TimeSync Plus
BASONET_sites$uniqID <- paste0(BASONET_sites$myID, "_", BASONET_sites$Date)
BASONET_sites$yod <- NA

dim(BASONET_sites[BASONET_sites$Date ==2021,]) ## 430 plots
dim(BASONET_sites[BASONET_sites$Date ==2001,]) ## 217 plots

### Write to shapefile
st_write(obj=BASONET_sites, dsn="BASONET_sites.shp", append=FALSE)
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/RMarkdown")
st_write(obj=BASONET_sites, dsn="BASONET_sites.shp", append=FALSE)


# 6.0 Pre-processing - average measurements by horizon layer -------------------

load("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Soil_data/R_Output/soildatasets_correctBASONETcoords.RData")

### Check one indicator
soil_datasets$TOC_Clay <- NA
soil_datasets[!is.na(soil_datasets$TOC) & !is.na(soil_datasets$Clay),]$TOC_Clay <-
  soil_datasets[!is.na(soil_datasets$TOC) & !is.na(soil_datasets$Clay),]$TOC /
  soil_datasets[!is.na(soil_datasets$TOC) & !is.na(soil_datasets$Clay),]$Clay

ggplot(soil_datasets, aes(y=TOC_Clay, Dataset)) +
  geom_boxplot() +
  geom_jitter(width = 0.3,alpha = 0.2) +
  ylim(0,1) +
  geom_hline(yintercept = 1/13,  color = "red")+
geom_hline(yintercept = 1/10, color = "blue")+
  geom_hline(yintercept = 1/8, color = "chartreuse")

# save.image("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/SoilDatasets/MergedDatasets_05012023.RData")
# load("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/SoilDatasets/MergedDatasets_05012023.RData")


## BASONET has same plot name assigned to differnt coordinates. I create a new ID for the whole dataset
### They also have several analysis by horizon (three samples per horizon at some plots)

## Create a column newID grouping by coordinates and date (just in case the same site was sampled in different dates)
## for example, LUCAS or BASONET
## this ID will be used for standarising by depth
soil_datasets$ID_coords <- paste0(soil_datasets$UTM_X, "_", soil_datasets$UTM_Y)
length(unique(soil_datasets$ID_coords)) # 1306

soil_datasets <- soil_datasets %>%
  arrange(., Dataset, -Date, myID, ID_coords, Upper_limit, Lower_limit) %>%
  group_by(., ID_coords,-Date) %>% 
  mutate(., newID =  cur_group_id())

length(unique(soil_datasets$newID)) # 1537

###  Arrange observations by ID, profle and depth
### For each sampling location, e.g., myID
soil_datasets <- arrange(soil_datasets, Dataset, -Date, newID, Upper_limit, Lower_limit) %>% as.data.frame()
view(soil_datasets)

## how many are individual samples that need to be averaged?
soil_datasets.summary.layers <- soil_datasets %>% 
  group_by(., newID, Upper_limit, Lower_limit) %>%
  summarise(., N = n())
soil_datasets.summary.layers %>% filter(., N >1) 

which.ones <- unique(soil_datasets.summary.layers[soil_datasets.summary.layers$N >1, ]$newID)
soil_datasets[soil_datasets$newID%in% which.ones ,]

kkeps <- unique(setdiff(soil_datasets$newID, which.ones))
length(unique(soil_datasets$newID))
length(which.ones) + length(kkeps)

### Create a "myLayer" variable
soil_datasets <- arrange(soil_datasets, newID, Upper_limit, Lower_limit)
soil_datasets <- soil_datasets %>% group_by(.,newID) %>% arrange(., newID, Upper_limit, Lower_limit) %>%
  mutate(., myLayer =  row_number()) %>% as.data.frame()

### Split into two dataframes
soil_datasets.Single <- soil_datasets[soil_datasets$newID %in% kkeps,]
soil_datasets.Multiple <- soil_datasets[soil_datasets$newID %in% which.ones,]
soil_datasets.Multiple <- arrange(soil_datasets.Multiple, newID, Upper_limit, Lower_limit )
soil_datasets.Multiple <- as.data.frame(soil_datasets.Multiple)

### Average by combination newID, Upper_limit, Lower_limit
soil_datasets.Multiple.Ave <- soil_datasets.Multiple[0,]
col.order <- colnames(soil_datasets.Multiple.Ave)

for (i in 1:length(which.ones)){
  print(i)
  df.i <- soil_datasets.Multiple[soil_datasets.Multiple$newID == which.ones[[i]],]
  ### Summarize numeric columns
  df.numeric <- df.i %>% group_by(., newID, Upper_limit, Lower_limit) %>% summarise_if(is.numeric, mean, na.rm = TRUE)
  df.numeric <- arrange(df.numeric,newID, Upper_limit, Lower_limit)
  df.numeric <- as.data.frame(df.numeric)
  
  ### Recalculate myLayer
  df.numeric$myLayer <- 1:nrow(df.numeric)
  
  ### Character variables and keep only the first rows
  df.char <- df.i %>% select_if(is.character)
  df.char <- df.char[1:nrow(df.numeric),]
  
  ### join both dataframes
  df.o <- cbind(df.numeric,df.char)
  ### Reorder columns
  df.o <- df.o[,col.order]
  ### Add to empty table
  soil_datasets.Multiple.Ave <- rbind(soil_datasets.Multiple.Ave, df.o )
}
rm(col.order,i,df.i,df.numeric,df.char,df.o, which.ones, kkeps)

### Merge together
soil_datasets_2 <- rbind(soil_datasets.Single, soil_datasets.Multiple.Ave)


# 6. Check for outliers ---------------------------------------------------

### Create index column
soil_datasets_2$index <- 1:nrow(soil_datasets_2)

## Remember I can change < LOD to NA for P (Basonet)
ggplot(soil_datasets_2, aes(y=TOC, Dataset)) + geom_boxplot()+ geom_jitter(width = 0.3,alpha = 0.2) 

### TOC
### Who are these horizons with TOC > 12%? I use limit from Soil Taxonomy to separate mineral from organic soils
TOC_organic_index <- soil_datasets_2[soil_datasets_2$TOC > 12 & !is.na(soil_datasets_2$TOC),]$index
soil_datasets_2[soil_datasets_2$index %in% TOC_organic_index,]$TOC <- NA

### Other properties

# pH is in the right scale
ggplot(soil_datasets_2, aes(y=pH, Dataset)) +
  geom_boxplot()+ geom_jitter(width = 0.3,alpha = 0.2)

# EC
ggplot(soil_datasets_2, aes(y=EC, Dataset)) + 
  geom_boxplot()+ geom_jitter(width = 0.3,alpha = 0.2)

### Very strange values
summary(soil_datasets_2[soil_datasets_2$Dataset =="CARBOSOL" & !is.na(soil_datasets_2$EC),]$EC)
hist(soil_datasets_2[soil_datasets_2$Dataset =="CARBOSOL" & !is.na(soil_datasets_2$EC),]$EC)
ggplot(soil_datasets_2[soil_datasets_2$Dataset =="CARBOSOL" & !is.na(soil_datasets_2$EC),],
       aes(y=EC, Dataset)) + geom_boxplot()+ geom_jitter(width = 0.3,alpha = 0.2) + ylim(0,2)
ggplot(soil_datasets_2[soil_datasets_2$Dataset !="CARBOSOL" & !is.na(soil_datasets_2$EC),],
       aes(y=EC, Dataset)) + geom_boxplot()+ geom_jitter(width = 0.3,alpha = 0.2)

### I will exclude EC - CARBOSOL from the analysis. 
### There may be a problem of units? Maybe they were nor dS / m?
#soil_datasets_3 <- soil_datasets_2[,colnames(soil_datasets_2)!="EC"]

## C_N
ggplot(soil_datasets_2, aes(y=C_N, Dataset)) + 
  geom_boxplot() + geom_jitter(width = 0.3,alpha = 0.2)
### I exclude those observations from very high TOC content
soil_datasets_2[soil_datasets_2$index %in% TOC_organic_index, ]$C_N <- NA


### HERE today 11/01/2024
# save.image("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/SoilDatasets/Cleaning_11012023.RData")

### HERE today 15/01/2024 - correcting EC data in BASONET 2021
save.image("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/SoilDatasets/Cleaning_15012023.RData")

### TN ppm
ggplot(soil_datasets_2, aes(y=TN_ppm, Dataset)) + 
  geom_boxplot() + geom_jitter(width = 0.3,alpha = 0.2)
soil_datasets_2[soil_datasets_2$TN_ppm > 10000 & !is.na(soil_datasets_2$TN_ppm), ]

### Extractable Phosphorus
ggplot(soil_datasets_2, aes(y=P_ppm, Dataset)) + 
  geom_boxplot() + geom_jitter(width = 0.3,alpha = 0.2)
### and cations
ggplot(soil_datasets_2, aes(y=K_ppm, Dataset)) + 
  geom_boxplot() + geom_jitter(width = 0.3,alpha = 0.2)
ggplot(soil_datasets_2, aes(y=Ca_ppm, Dataset)) + 
  geom_boxplot() + geom_jitter(width = 0.3,alpha = 0.2)
ggplot(soil_datasets_2, aes(y=Mg_ppm, Dataset)) + 
  geom_boxplot() + geom_jitter(width = 0.3,alpha = 0.2)

### Na ppm - these values are very clearly much higher than the rest
ggplot(soil_datasets_2, aes(y=Na_ppm, Dataset)) + geom_boxplot() +
  geom_jitter(width = 0.3,alpha = 0.2) +
  geom_hline(yintercept = 2500,  color = "red")
soil_datasets_2[soil_datasets_2$Na_ppm > 2500 &
                  !is.na(soil_datasets_2$Na_ppm),]$Na_ppm <- NA

### CEC
ggplot(soil_datasets_2, aes(y=CEC, Dataset)) + geom_boxplot() + geom_jitter(width = 0.3,alpha = 0.2)
ggplot(soil_datasets_2, aes(y=CECef, Dataset)) + geom_boxplot() + geom_jitter(width = 0.3,alpha = 0.2)

### Coarse fragments
ggplot(soil_datasets_2, aes(y=CoarseFrag, Dataset)) + geom_boxplot() + geom_jitter(width = 0.3,alpha = 0.2)
ggplot(soil_datasets_2, aes(y=CoarseFrag_Grav, Dataset)) + geom_boxplot() + geom_jitter(width = 0.3,alpha = 0.2)

### Carbonates
ggplot(soil_datasets_2, aes(y=Carbonates, Dataset)) + geom_boxplot() + geom_jitter(width = 0.3,alpha = 0.2)

### Bulk density
ggplot(soil_datasets_2, aes(y=Bulk_Density, Dataset)) + geom_boxplot() + geom_jitter(width = 0.3,alpha = 0.2)


# 7. Standarise to common depth intervals ---------------------------------

### the minimum dataset is:
### Particle size distribution
### TOC
### Bulk density
### TOC:clay
### Extractable phosphorus
### Total nitrogen
### pH
### EC
### Potassium
### Calcium
### Magnesium
### Sodium









### Fit mass-preserving splines for the depths...
plot(tern_SOCfractions.4$lower, tern_SOCfractions.4$upper)
abline(0,1, col="blue")
sort(unique(tern_SOCfractions.4$LDepth))
#  [1]  2  3  5  6  7  8  9 10 11 12 15 16 17 20 22 25 27 28 30
sort(unique(tern_SOCfractions.4$UDepth))
# [1]  0  2  3  5  6 10 15 20
length(unique(tern_SOCfractions.4$coord_ID))
length(unique(paste(tern_SOCfractions.4$coord_ID, tern_SOCfractions.4$UDepth)))

### How many of these have only one measurement per location?
## =1 --> aqp
## >1 --> Splines

### Now, separate those that have one observations per location from those with several per location
detach("package:aqp", unload=TRUE)
tern.summary.site <- tern_SOCfractions.4 %>% 
  group_by(., myID) %>%
  summarise(., N = n())%>% as.data.frame()
tern.summary.site %>% filter(., N ==1) %>% as.data.frame()

too.many <- unique(tern.summary.site[tern.summary.site$N > 1,]$myID)
kkeps <- unique(setdiff(tern.summary.site$myID, too.many))
length(unique(tern_SOCfractions.4$myID))
length(too.many) + length(kkeps)

### Split into two dataframes
tern_SOCfractions.4.aqp <- tern_SOCfractions.4[tern_SOCfractions.4$myID %in% kkeps,]
tern_SOCfractions.4.spl <- tern_SOCfractions.4[tern_SOCfractions.4$myID %in% too.many,]
tern_SOCfractions.4.aqp <- as.data.frame(tern_SOCfractions.4.aqp)
str(tern_SOCfractions.4.aqp)
summary(tern_SOCfractions.4.aqp$LDepth)
summary(tern_SOCfractions.4.aqp$UDepth)
tern_SOCfractions.4.aqp.bck <- tern_SOCfractions.4.aqp
#tern_SOCfractions.4.aqp <- tern_SOCfractions.4.aqp.bck

library(aqp)
# upgrade to SoilProfileCollection
# 'myID' is the name of the column containing the profile ID
# 'UDepth' is the name of the column containing horizon upper boundaries
# 'LDepth' is the name of the column containing horizon lower boundaries
#depths(tern_SOCfractions.4.aqp) <- myID ~ UDepth + LDepth
#print(tern_SOCfractions.4.aqp)
#checkSPC(tern_SOCfractions.4.aqp)
#spc_in_sync(tern_SOCfractions.4.aqp)
#checkHzDepthLogic(tern_SOCfractions.4.aqp)

### change of support according to GSM depths
### https://ncss-tech.github.io/AQP/aqp/aqp-intro.html#14_Aggregating_Soil_Profile_Collections_Along_Regular_%E2%80%9CSlabs%E2%80%9D

### Doing all at the same time
gsm.depths <- c(0, 5, 15, 30)
desired.properties <- c("POC_p","ROC_p","HOC_p","Vp")
forms <- list(paste("myID ~ POC_p"), paste("myID ~ ROC_p"),paste("myID ~ HOC_p"), paste("myID ~ Vp")) 

### Remember columns in the original dataset
# [1]  "mylabs"         "state_c"        "nir_lab"        "zone"           "easting"        "northing"       "matchedbarcode" "upper"          "lower"         
# [10] "type"           "pH"             "ec"             "longitude"      "latitude"       "OC"             "POC"            "HOC"            "ROC"           
# [19] "Lat_WGS84"      "Long_WGS84"     "coord_ID"       "UDepth"         "LDepth"         "HOC.st"         "POC.st"         "ROC.st"         "myLayer"       
# [28] "myID"           "TOC"            "SumFracOC.st"

## Let's keep "myID" "Lat_WGS84" "Long_WGS84" "coord_ID"  
out.aqp <- list()
for(j in 1:length(desired.properties)) {
  
  ### Eliminate NA for each property
  nona.idx <- !is.na(tern_SOCfractions.4.aqp[,desired.properties[[j]]])
  tern_SOCfractions.4.aqp.j <- tern_SOCfractions.4.aqp[nona.idx,]
  tern_SOCfractions.4.aqp.j <- tern_SOCfractions.4.aqp.j[, c("myID","Lat_WGS84","Long_WGS84","coord_ID","UDepth","LDepth","POC_p","ROC_p","HOC_p","Vp")]
  tern_SOCfractions.4.aqp.j.bck <- tern_SOCfractions.4.aqp.j
  depths(tern_SOCfractions.4.aqp.j) <- myID ~ UDepth + LDepth
  
  # print(tern_SOCfractions.4.aqp.j)
  # checkSPC(tern_SOCfractions.4.aqp.j)
  # spc_in_sync(tern_SOCfractions.4.aqp.j)
  # checkHzDepthLogic(tern_SOCfractions.4.aqp.j)
  
  ### For one property
  tern.gsm.j <- aqp::slab(tern_SOCfractions.4.aqp.j, fm= as.formula(forms[[j]]), slab.structure = gsm.depths, slab.fun = mean, na.rm=TRUE)
  head(tern.gsm.j)
  #ids <- aqp::profile_id(tern_SOCfractions.4.aqp.j)
  
  ### when contribution is less than 25% I set it to NA, only when another layer in that location already has the value assigned
  # tern.gsm.j.rev <- tern.gsm.j[0,]
  
  # for(i in 1:length(ids)){
  #   
  #   df.i <- tern.gsm.j[tern.gsm.j$myID ==ids[[i]], ]
  #   N.profile <- sum(!is.na(df.i$value))
  #   
  #   if(N.profile ==1){
  #     df.i$value <- ifelse(is.nan(df.i$value), NA, df.i$value)
  #     tern.gsm.j.rev <- rbind(tern.gsm.j.rev,df.i)
  #   }
  #   
  #   if(N.profile > 1 & any(df.i$contributing_fraction >= 0.25)){
  #     df.i$value <- ifelse(df.i$contributing_fraction < 0.25, NA, df.i$value)
  #     tern.gsm.j.rev <- rbind(tern.gsm.j.rev,df.i)
  #   } else {
  #     max.contrib <- max(df.i$contributing_fraction)
  #     df.i$value <- ifelse(df.i$contributing_fraction < max.contrib, NA, df.i$value)
  #     tern.gsm.j.rev <- rbind(tern.gsm.j.rev,df.i)
  #   }
  # }
  tern.gsm.j$value <- ifelse(is.nan(tern.gsm.j$value), NA, tern.gsm.j$value)
  # reshape to wide format
  tern.gsm.j$GMS_layer <- paste0(tern.gsm.j$top,"-",tern.gsm.j$bottom, " cm")
  tern.gsm.j.wide <- dcast(tern.gsm.j, myID + variable ~ GMS_layer , value.var = 'value', fun.aggregate=mean)
  ## Add some columns
  tern.gsm.j.wide <- merge(tern.gsm.j.wide, 
                           tern_SOCfractions.4.aqp.bck[, c("myID","coord_ID","Lat_WGS84","Long_WGS84")],
                           by="myID", all.x=TRUE)
  out.aqp[[j]] <- tern.gsm.j.wide
}

tern.aqp <- out.aqp

rm("tern.gsm.j",  "tern.gsm.j.rev","tern.gsm.j.rev.wide", "tern_SOCfractions.4.aqp.j",
   "tern_SOCfractions.4.aqp.j.bck", df.i,i,j,nona.idx, max.contrib,N.profile, forms, ids, out.aqp)

### Now perform splines and keep GSM intervals (we should just keep between 0-30 cm)
colnames(tern_SOCfractions.4.spl)
tern_SOCfractions.4.spl <- arrange(tern_SOCfractions.4.spl, myID, UDepth, LDepth) %>% as.data.frame()

### Unique SiteID
sampling.profiles <- unique(tern_SOCfractions.4.spl$myID)
desired.properties <- c("POC_p","ROC_p","HOC_p","Vp")

out.splines <- list()
for(j in 1:length(desired.properties)) {
  df.j <- as.data.frame(tern_SOCfractions.4.spl[, c("myID","UDepth","LDepth",desired.properties[[j]])])
  ### Fit mass-preserving spline for j property
  eaFit <- ea_spline(obj=df.j, var.name = desired.properties[[j]], d = t(c(0,5,15,30)), lam = 0.1, vlow = 0, show.progress = TRUE)
  eaFit$obs.preds[,2] <- as.numeric(eaFit$obs.preds[,2])
  eaFit$obs.preds[,3] <- as.numeric(eaFit$obs.preds[,3])
  eaFit$obs.preds[,4] <- as.numeric(eaFit$obs.preds[,4])
  out.splines[[j]] <- eaFit
  
}

TERN.POC_p <- out.splines[[1]]
TERN.ROC_p <- out.splines[[2]]
TERN.HOC_p <- out.splines[[3]]
TERN.Vp <- out.splines[[4]]

problematicIDs <- setdiff(sampling.profiles, TERN.POC_p$harmonised$id)
rm(out.splines, desired.properties, sampling.profiles,eaFit, df.j,j, problematicIDs, kkeps, too.many,gsm.depths, j,
   tern.summary.site, tern_SOCfractions.4.aqp, tern_SOCfractions.4.aqp.bck, tern_SOCfractions.4.spl, tern.gsm.j.wide)


# 6.1 Sand, silt, clay ----------------------------------------------------

### Add to 100% ?


