############################################################################################################################################
###  Preprocessing of environmental covariates for generating pedogenon classes for the Basque Country
###
###  In this script: Make maps to visualize covariates

### Desired extent: Basque Country
### Resolution: 25m
### CRS: EPSG=25830

###  Author: Mercedes Roman Dobarco
###  Date: 29/09/2023

####### Load packages
### Spatial
library(sf)
library(terra)

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


# Auxiliary layers --------------------------------------------------------

### Load administrative boundaries of Euskadi
### the CAPV ("U:/Covariates/Administrative/CB_CAPV_5000_ETRS89/U11.CB_CAPV_5000_ETRS89.shp")
euskadi <- st_read(("U:/Covariates/Administrative/CB_CAPV_5000_ETRS89/U11.CB_CAPV_5000_ETRS89.shp"))

### and the rest of Spain
HomeDir <- "U:/Covariates/"  ### Change this with your Home directory 

setwd(paste0(HomeDir,"Administrative/SHP_ETRS89/recintos_autonomicas_inspire_peninbal_etrs89/"))
autonomias <- st_read("recintos_autonomicas_inspire_peninbal_etrs89.shp")
autonomias25830 <- st_transform(autonomias, 25830)
st_bbox(autonomias25830)

### Test countries from EC
world <- st_read("U:/Covariates/Europe/ref-countries-2020-01m.shp/CNTR_RG_01M_2020_4326.shp/CNTR_RG_01M_2020_4326.shp")
world25830 <- st_transform(world, 25830)

### subset neighbours around Euskadi
buff_Eus <- st_buffer(euskadi, dist=1000000)
plot(buff_Eus["CCAA"])
world_eus_int <- st_intersects(world25830, buff_Eus, sparse = FALSE)
world_eus <- world25830[world_eus_int,]
plot(world_eus["CNTR_ID"])

### Grey fill
tmap_mode("plot")
tmap_options()

tm_shape(world_eus, bbox=c(-300000, 3700000, 1300000, 5100000)) +
  tm_borders() +
    tm_shape(autonomias25830) +
  tm_fill(col="gray90")+
  tm_borders() +
  tm_shape(euskadi) +
  tm_fill(col = "darkslategray3") +
  tm_borders(lwd=1.5) +
  tm_compass(type = "4star", position = c("left", "top")) +
  tm_scale_bar(breaks = c(0, 100, 200, 300, 400), text.size = 0.8, position = c("right", "bottom")) +
  tm_layout(main.title="Basque Country")+
  tm_grid(lines=FALSE)
 


# Relief ------------------------------------------------------------------

#### Plot some covariates 
setwd("C:/Covariates/Euskadi/Relief/")

outer <- tm_compass(type = "arrow", position = c("left", "top")) +
  tm_scale_bar(breaks = c(0, 20, 40), text.size = 0.8, position = c("right", "bottom")) +
  tm_shape(world_eus) +
  tm_borders() +
  tm_shape(autonomias25830) +
  tm_borders() +
  tm_shape(euskadi) +
  tm_borders(lwd=2, col="black") +
  tm_grid(lines=FALSE)


outer2 <- tm_compass(type = "arrow", position = c("left", "top")) +
  tm_scale_bar(breaks = c(0, 20, 40), text.size = 0.8, position = c("right", "bottom")) +
  tm_shape(world_eus) +
  tm_borders() +
  tm_shape(autonomias25830) +
  tm_borders() +
  tm_shape(euskadi) +
  tm_borders(lwd=1, col="black") 

library(stars)
dem <- terra::rast("mdt_lidar_2017_25m_etrs89.tif")
hillshade <- terra::rast("hillshade.tif")
bboxplot <- c(455000, 4699000, 609000, 4815000)

tm_shape(hillshade, bbox=bboxplot) +
  tm_raster(palette = gray(0:100 / 100), n = 100, legend.show = FALSE) +
  tm_shape(dem)+
  tm_raster(style = "cont",
            alpha = 0.6,
            title = "Elevation (m)", 
            palette = terrain.colors(50)) +
  tm_layout(legend.outside=TRUE, legend.outside.position="right", 
            legend.title.size=1.5,
            legend.text.size=1.1) +
  outer

### Load different layers

### Covariates subset
covariates.selection <- c("mdt_lidar_2017_25m_etrs89.tif", # DEM
                          "slope.tif", "slope_5.tif", "slope_10.tif", # Slope calculated from 3,5,10 cells
                          "easterness.tif", # Easterness 
                          "northerness.tif", # northerness
                          "nor_slope.tif" , # northness x slope
                          "twi_saga.tif", # SAGA Topographic wetness index
                          "p_curv.tif", "pr_curv_5.tif", "pr_curv_10.tif", # profile curvature (3,5,10)
                          "pl_curv_5.tif", "pl_curv_10.tif", # planar curvature (5,10)
                          "lg_curv_5.tif", "lg_curv_10.tif", # longitudinal curvature (5,10)
                          "t_curv.tif", "cs_curv_5.tif", "cs_curv_10.tif" , # tangential curvature (3) cross-sectional curvature (5,10) (same)
                          "mrvbf.tif", "mrrtf.tif", ### MRVBF, MRRTF
                          "slope_height.tif", # Slope Height
                          "norm_height.tif", # Normalized height
                          "st_height.tif", # Standardized Height
                          "valley_depth.tif", # Valley depth
                          "mid_slope.tif",  # Mid slope
                          "tpi_8_3.tif", "tpi_20_5.tif")# , # Multiscale topographic position index, calculated with a search radius of 8,20 cells, and 3,5 scales respectively
                         # "geomorphon.tif", "geomorphon2.tif")  ## Geomorphons (search radius of 11, skip 1 & 3, flat 1 & 1.5)

relief.rast <- terra::rast(covariates.selection)
names(relief.rast) <- c("dem","slope","slope_5","slope_10",
                        "easterness","northerness","nor_slope",
                        "twi","p_curv","pr_curv_5","pr_curv_10",
                        "pl_curv_5","pl_curv_10","lg_curv_5","lg_curv_10",
                        "t_curv","cs_curv_5","cs_curv_10",
                        "mrvbf","mrrtf",
                        "slope_height","norm_height","st_height","valley_depth","mid_slope",
                        "tpi_8_3","tpi_20_5")

library(colorspace)
colorspace::hcl_palettes(plot = TRUE)
rm(covariates.selection)

### Slopes
tm_shape(relief.rast[[2:4]], bbox=bboxplot) +
  tm_facets(ncol = 3) +
  tm_layout(panel.labels = c("slope (w=3)", "slope (w=5)","slope (w=10)")) +
  tm_raster(style= "quantile", n=5, palette = rev(hcl.colors(5, "OrYel"))) +
  # tm_layout(legend.outside=TRUE, 
  #           legend.outside.position="right", 
  #           legend.title.size=1.5,
  #           legend.text.size=1.1) +
  outer2

slope_map <- tm_shape(relief.rast$slope, bbox=bboxplot) +
  # tm_facets(ncol = 3) +
  #tm_layout(panel.labels = c("slope (w=3)", "slope (w=5)","slope (w=10)"))+
  tm_raster(style= "quantile", n=5, palette = rev(hcl.colors(5, "OrYel")),
            title="Slope (%)",
            labels= c("0-8 - 8.5", "8.5 - 19.2", "19.2 - 31.0", "31.0 - 45.2", "> 45.2"))+
  tm_layout(legend.title.size=1.5,
            legend.text.size=1.2,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.position = c(0.02,0.02)) +
  outer2

slopew5_map <- tm_shape(relief.rast$slope_5, bbox=bboxplot) +
  tm_raster(style= "quantile", n=5, palette = rev(hcl.colors(5, "OrYel")),
            title="Slope w=5 (%)")+ #,
            #labels= c("0-8 - 8.5", "8.5 - 19.2", "19.2 - 31.0", "31.0 - 45.2", "> 45.2"))+
  tm_layout(legend.title.size=1.5,
            legend.text.size=1.2,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.position = c(0.02,0.02)) +
  outer2

slopew10_map <- tm_shape(relief.rast$slope_10, bbox=bboxplot) +
  tm_raster(style= "quantile", n=5, palette = rev(hcl.colors(5, "OrYel")),
            title="Slope w=10 (%)")+ #,
  #labels= c("0-8 - 8.5", "8.5 - 19.2", "19.2 - 31.0", "31.0 - 45.2", "> 45.2"))+
  tm_layout(legend.title.size=1.5,
            legend.text.size=1.2,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.position = c(0.02,0.02)) +
  outer2
tmap_arrange(slope_map, slopew5_map, slopew10_map, ncol = 3)

### Easterness, northness, northnessxslope
easterness_map <- tm_shape(relief.rast$easterness, bbox=bboxplot) +
   tm_raster(style= "cont", palette = "RdBu",
            title="Easterness")+
  tm_layout(legend.title.size=1.5,
            legend.text.size=1.2,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.position = c(0.02,0.02)) +
  outer2

northness_map <- tm_shape(relief.rast$northerness, bbox=bboxplot) +
  tm_raster(style= "cont", palette = "RdBu",
            title="Northness")+
  tm_layout(legend.title.size=1.5,
            legend.text.size=1.2,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.position = c(0.02,0.02)) +
  outer2

norslope_map <- tm_shape(relief.rast$nor_slope, bbox=bboxplot) +
  tm_raster(style= "cont", palette = "RdBu",
            title="Northness x slope")+
  tm_layout(legend.title.size=1.5,
            legend.text.size=1.2,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.position = c(0.02,0.02)) +
  outer2

### TWI
twi <- tm_shape(relief.rast$twi, bbox=bboxplot) +
  tm_raster(style= "cont", palette = "-viridis",
            title="TWI")+
  tm_layout(legend.title.size=1.5,
            legend.text.size=1.2,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.position = c(0.02,0.02)) +
  outer2

### MRVBF
mrvbf <- tm_shape(relief.rast$mrvbf, bbox=bboxplot) +
  tm_raster(style= "cont", palette = "-viridis",
            title="MRVBF")+
  tm_layout(legend.title.size=1.5,
            legend.text.size=1.2,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.position = c(0.02,0.02)) +
  outer2

### MRRTF
mrrtf <- tm_shape(relief.rast$mrrtf, bbox=bboxplot) +
  tm_raster(style= "cont", palette = "-viridis",
            title="MRRTF")+
  tm_layout(legend.title.size=1.5,
            legend.text.size=1.2,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.position = c(0.02,0.02)) +
  outer2

### Curvatures
# tm_shape(relief.rast[[9:18]], bbox=bboxplot) +
#   tm_facets(ncol = 4) +
#  # tm_layout(panel.labels = c("slope (w=3)", "slope (w=5)","slope (w=10)")) +
#   tm_raster(style= "cont") +
#   outer2
plot(relief.rast[[9:18]])

### Channel attributes
get_brewer_pal("Blues", n = 7)

slope_height <- tm_shape(relief.rast$slope_height, bbox=bboxplot) +
  tm_raster(style= "cont", palette = "-viridis",
            title="Slope height")+
  tm_layout(legend.title.size=1.5,
            legend.text.size=1.2,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.position = c(0.02,0.02)) +
  outer2
norm_height <- tm_shape(relief.rast$norm_height, bbox=bboxplot) +
  tm_raster(style= "cont", palette = "-viridis",
            title="Normalized height")+
  tm_layout(legend.title.size=1.5,
            legend.text.size=1.2,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.position = c(0.02,0.02)) +
  outer2
st_height <- tm_shape(relief.rast$st_height, bbox=bboxplot) +
  tm_raster(style= "cont", palette = "-viridis",
            title="Standarized height")+
  tm_layout(legend.title.size=1.5,
            legend.text.size=1.2,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.position = c(0.02,0.02)) +
  outer2
valley_depth <- tm_shape(relief.rast$valley_depth, bbox=bboxplot) +
  tm_raster(style= "cont", palette = "-viridis",
            title="Valley depth")+
  tm_layout(legend.title.size=1.5,
            legend.text.size=1.2,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.position = c(0.02,0.02)) +
  outer2
mid_slope <- tm_shape(relief.rast$mid_slope, bbox=bboxplot) +
  tm_raster(style= "cont", palette = "-viridis",
            title="Mid slope position")+
  tm_layout(legend.title.size=1.5,
            legend.text.size=1.2,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.position = c(0.02,0.02)) +
  outer2

tpi_8_3 <- tm_shape(relief.rast$tpi_8_3, bbox=bboxplot) +
  tm_raster(style= "cont", palette = "-viridis",
            title="TPI 1")+
  tm_layout(legend.title.size=1.5,
            legend.text.size=1.2,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.position = c(0.02,0.02)) +
  outer2
tpi_20_5 <- tm_shape(relief.rast$tpi_20_5, bbox=bboxplot) +
  tm_raster(style= "cont", palette = "BrBG",
            title="TPI 2")+
  tm_layout(legend.title.size=1.5,
            legend.text.size=1.2,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.position = c(0.02,0.02)) +
  outer2

# Soil covariates -----------------------------------------------------

setwd("C:/Covariates/Euskadi/Soil/SoilGrids/")
soil.files <- list.files(pattern=".tif$")
soil.rast <- terra::rast(soil.files)
rm(soil.files)

names(soil.rast) <- c( "cec_100-200cm","cec_30-60cm","cec_60-100cm",
                       "clay_100-200cm","clay_30-60cm","clay_60-100",
                       "sand_100-200cm","sand_30-60cm","sand_60-100cm_",
                       "silt_100-200cm","silt_30-60cm","silt_60-100cm")
texture  <- c(soil.rast$`clay_30-60cm`, soil.rast$`clay_60-100`, soil.rast$`clay_100-200cm`,
              soil.rast$`sand_30-60cm`, soil.rast$`silt_60-100cm`, soil.rast$`silt_100-200cm`,
              soil.rast$`sand_30-60cm`, soil.rast$`sand_60-100cm_`, soil.rast$`sand_100-200cm`)
         
texture_maps <- tm_shape(texture, bbox=bboxplot) +
  tm_raster(style= "fixed", breaks = c(0,100,200,300,400,500,600,800), 
            n=7, 
            palette = rev(hcl.colors(7, "YlOrBr")),
            title="Particle content (mg)") +
  tm_facets(ncol = 3) +
  tm_layout(panel.labels = c("Clay (30-60cm)","Clay (60-100cm)", "Clay (100-200cm)",
                             "Silt (30-60cm)","Silt (60-100cm)", "Silt (100-200cm)",
                             "Sand (30-60cm)","Sand (60-100cm)", "Sand (100-200cm)"))+

  outer2

cec_map <- tm_shape(c(soil.rast$`cec_30-60cm`, soil.rast$`cec_60-100cm`, soil.rast$`cec_100-200cm`),
                    bbox=bboxplot, projection=25830) +
  tm_raster(style= "cont",
            palette = rev(hcl.colors("ag_Sunset")),
            title="CEC") +
  tm_facets(ncol = 3) +
  tm_layout(panel.labels = c("CEC (30-60cm)","CEC (60-100cm)", "CEC (100-200cm)")) +
  outer2

# Anthromes ---------------------------------------------------------------

AnthDir <- "C:/Covariates/Euskadi/Organisms/anthromes/"
setwd(AnthDir)
#anthro.files <- list.files(pattern=".tif$")
my_order_anthromes <- c("anthrome_10000BC.tif","anthrome_9000BC.tif","anthrome_8000BC.tif","anthrome_7000BC.tif",
                        "anthrome_6000BC.tif", "anthrome_5000BC.tif","anthrome_4000BC.tif","anthrome_3000BC.tif",
                        "anthrome_2000BC.tif", "anthrome_1000BC.tif","anthrome_0AD.tif",
                        "anthrome_100AD.tif",  "anthrome_200AD.tif", "anthrome_300AD.tif", "anthrome_400AD.tif",
                        "anthrome_500AD.tif",  "anthrome_600AD.tif", "anthrome_700AD.tif", "anthrome_800AD.tif",
                        "anthrome_900AD.tif",  "anthrome_1000AD.tif","anthrome_1100AD.tif","anthrome_1200AD.tif",
                        "anthrome_1300AD.tif", "anthrome_1400AD.tif","anthrome_1500AD.tif","anthrome_1600AD.tif",
                        "anthrome_1700AD.tif", "anthrome_1710AD.tif","anthrome_1720AD.tif","anthrome_1730AD.tif", 
                        "anthrome_1740AD.tif", "anthrome_1750AD.tif","anthrome_1760AD.tif","anthrome_1770AD.tif",
                        "anthrome_1780AD.tif", "anthrome_1790AD.tif","anthrome_1800AD.tif","anthrome_1810AD.tif",
                        "anthrome_1820AD.tif", "anthrome_1830AD.tif","anthrome_1840AD.tif","anthrome_1850AD.tif",
                        "anthrome_1860AD.tif", "anthrome_1870AD.tif","anthrome_1880AD.tif","anthrome_1890AD.tif",
                        "anthrome_1900AD.tif", "anthrome_1910AD.tif","anthrome_1920AD.tif","anthrome_1930AD.tif", 
                        "anthrome_1940AD.tif", "anthrome_1950AD.tif","anthrome_1960AD.tif","anthrome_1970AD.tif",
                        "anthrome_1980AD.tif", "anthrome_1990AD.tif","anthrome_2000AD.tif","anthrome_2001AD.tif",
                        "anthrome_2002AD.tif", "anthrome_2003AD.tif","anthrome_2004AD.tif","anthrome_2005AD.tif",
                        "anthrome_2006AD.tif", "anthrome_2007AD.tif","anthrome_2008AD.tif","anthrome_2009AD.tif",
                        "anthrome_2010AD.tif", "anthrome_2011AD.tif","anthrome_2012AD.tif","anthrome_2013AD.tif",
                        "anthrome_2014AD.tif", "anthrome_2015AD.tif","anthrome_2016AD.tif","anthrome_2017AD.tif")

names_anthromes <- gsub(my_order_anthromes,replacement = "", pattern=".tif" )
names_anthromes <- gsub(names_anthromes,replacement = "", pattern="anthrome_" )

### Create raster stack
anthro.rast <- terra::rast(my_order_anthromes)

### Colors by class
q17 <-c("11"="brown","12"="brown1",
        "21"= "lightskyblue","22"="cornflowerblue","23"="darkorchid","24"="#F161AE", 
        "31"="#72F6D2","32"="#CBC937","33"="#F2F081","34"="#F0EFA9",
        "41"="#E09F2A","42"="#ECB776","43"="#F6DCB5",
        "51"="#19AA34","52"="#5FCE2C","53"="#B6E284","54"="#F6F8E3")

categories <- data.frame(id=c(11,12,21,22,23,24,31,32,33,34,41,42,43,51,52,53,54),
                         class=c("Urban","Mixed settlements",
                                 "Rice villages","Irrigated villages","Rainfed villages",
                                 "Pastoral villages",
                                 "Residential irrigated croplands","Residential rainfed croplands",
                                 "Populated croplands", "Remote croplands",
                                 "Residential rangelands","Populated rangelands","Remote rangelands",
                                 "Residential woodlands","Populated woodlands","Remote woodlands",
                                 "Inhabited drylands"))
anthro.rast.cat <- anthro.rast

library(foreach)
foreach (r = 1:nlyr(anthro.rast.cat)) %do% {
  ### Assign levels
  levels(anthro.rast.cat[[r]]) <- categories
  ### Assign color table
  q17 <-c("brown","brown1",
          "lightskyblue","cornflowerblue","darkorchid","#F161AE", 
          "#72F6D2","#CBC937","#F2F081","#F0EFA9",
          "#E09F2A","#ECB776","#F6DCB5",
          "#19AA34","#5FCE2C","#B6E284","#F6F8E3")
  
    coltab(anthro.rast.cat[[r]]) <- data.frame(id=c(11,12,21,22,23,24,31,32,33,34,41,42,43,51,52,53,54), col=q17)
    }

names(anthro.rast.cat) <- names_anthromes

antht_1750AD <- tm_shape(anthro.rast.cat[["1750AD"]], bbox=bboxplot, projection=25830) +
  tm_raster(style= "cat",
            labels =c("Residential rainfed croplands","Remote croplands",
                      "Residential rangelands","Remote rangelands",
                      "Residential woodlands","Populated woodlands","Remote woodlands",
                      "Inhabited drylands") ,
            palette = c("#CBC937","#F0EFA9","#E09F2A","#F6DCB5",
                        "#19AA34","#5FCE2C","#B6E284","#F6F8E3"),
            title="1750 AD") +
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_shape(autonomias25830) +
  tm_borders() +
  tm_layout(legend.outside=TRUE,
            legend.outside.size=0.5,
            legend.outside.position="right",
            legend.title.size = 1.5,
            legend.text.size = 1.2,) +
  tm_shape(euskadi) +
  tm_borders(lwd=2, col="black") +
  tm_grid(lines=FALSE)

antht_1850AD <- tm_shape(anthro.rast.cat[["1850AD"]], bbox=bboxplot, projection=25830) +
  tm_raster(style= "cat",
            labels =c("Mixed settlements","Rainfed villages",
                      "Residential rainfed croplands","Populated croplands","Remote croplands",
                      "Residential rangelands",
                      "Residential woodlands","Populated woodlands","Remote woodlands",
                      "Inhabited drylands") ,
            palette = c("brown1","darkorchid",
                        "#CBC937","#F2F081","#F0EFA9",
                        "#E09F2A",
                        "#19AA34","#5FCE2C","#B6E284",
                        "#F6F8E3"),
            title="1850 AD") +
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_shape(autonomias25830) +
  tm_borders() +
  tm_layout(legend.outside=TRUE,
            legend.outside.size=0.5,
            legend.outside.position="right",
            legend.title.size = 1.5,
            legend.text.size = 1.2,) +
  tm_shape(euskadi) +
  tm_borders(lwd=2, col="black") +
  tm_grid(lines=FALSE)


# Climate -----------------------------------------------------------------

HomeDir <- "C:/Covariates/"
setwd(paste0(HomeDir,"Euskadi/Climate/"))
climate.files <- list.files(pattern=".tif")
climate.rast <- terra::rast(climate.files)
plot(climate.rast)

Bio1 <- tm_shape(climate.rast[[1]], bbox=bboxplot, projection=25830) +
  tm_raster(style= "cont", palette="-RdYlBu",
            title="MAT (°C)") +
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_shape(autonomias25830) +
  tm_borders() +
  tm_layout(legend.outside=TRUE,
            legend.outside.position="right",
            legend.title.size = 1.5,
            legend.text.size = 1.2,) +
  tm_shape(euskadi) +
  tm_borders(lwd=2, col="black") +
  tm_grid(lines=FALSE)

Bio12 <- tm_shape(climate.rast[[1]], bbox=bboxplot, projection=25830) +
  tm_raster(style= "cont", palette="PuBu",
            title="MAP (mm)") +
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_shape(autonomias25830) +
  tm_borders() +
  tm_layout(legend.outside=TRUE,
            legend.outside.position="right",
            legend.title.size = 1.5,
            legend.text.size = 1.2,) +
  tm_shape(euskadi) +
  tm_borders(lwd=2, col="black") +
  tm_grid(lines=FALSE)

ARI <- climate.rast$`CHELSA_pet_penman_mean_1981-2010_V.2.1_eus`/ climate.rast$`CHELSA_bio12_1981-2010_V.2.1_eus`

ARI_map <- tm_shape(ARI, bbox=bboxplot, projection=25830) +
  tm_raster(style= "quantile", palette="YlOrRd", n=10,
            title="Aridity Index") +
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_shape(autonomias25830) +
  tm_borders() +
  tm_layout(legend.outside=TRUE,
            legend.outside.position="right",
            legend.title.size = 1.5,
            legend.text.size = 1.2,) +
  tm_shape(euskadi) +
  tm_borders(lwd=2, col="black") +
  tm_grid(lines=FALSE)

Bio4 <- tm_shape(climate.rast$`CHELSA_bio4_1981-2010_V.2.1_eus`/100, bbox=bboxplot, projection=25830) +
  tm_raster(style= "cont", palette="Oranges", 
            title="Temp Seasonality (SD)") +
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_shape(autonomias25830) +
  tm_borders() +
  tm_layout(legend.outside=TRUE,
            legend.outside.position="right",
            legend.title.size = 1.5,
            legend.text.size = 1.2,) +
  tm_shape(euskadi) +
  tm_borders(lwd=2, col="black") +
  tm_grid(lines=FALSE)

# Soil map ----------------------------------------------------------------

### Soil map
setwd("C:/Covariates/Euskadi/Soil/")

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

edafo_eus$Order <- edafo_eus$UNIDAD
edafo_eus$Order <- ifelse(edafo_eus$Order %in% c("0 - Cauce","0 - Embalse","0 - Sin suelo"), NA, edafo_eus$Order)

### Plot
tm_shape(edafo_eus, bbox=bboxplot)+
  tm_fill(col="Order", 
          palette="Spectral")+
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_shape(autonomias25830) +
  tm_borders() +
  tm_layout(legend.outside=TRUE,
            legend.outside.position="right",
            legend.title.size = 1.5,
            legend.text.size = 1.2,) +
  tm_shape(euskadi) +
  tm_borders(lwd=2, col="black") +
  tm_grid(lines=FALSE)

### Plot
setwd("C:/Covariates/Euskadi/Soil/")
edafo_eus_r <- terra::rast("edafo_eus.tif")
tm_shape(edafo_eus_r, bbox=bboxplot)+
  tm_raster(style= "pretty", 
            palette="Spectral", 
            title="Soil Map Unit")+
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_shape(autonomias25830) +
  tm_borders() +
  tm_layout(legend.outside=TRUE,
            legend.outside.position="right",
            legend.title.size = 1.5,
            legend.text.size = 1.2,) +
  tm_shape(euskadi) +
  tm_borders(lwd=2, col="black") +
  tm_grid(lines=FALSE)


### The soil map for Bizkaia and Guipuzkoa
A <- 7.234 # thousand km2
Ssuborder <- 5.7984 *(A^0.2265)
Ssuborder ### 9 suborders
Sfam <- 15.2540 *(A^0.6315)
Sfam ### we would be talking about 53 soil classes at the family level
Sser <- 16.342 *(A^0.7106)
Sser ### and around 67 at the level of series

# Potential Vegetation ----------------------------------------------------------

potential_veg <- st_read("U:/Covariates/Euskadi/Organisms/VegPotencial/CT_VEGETACION_POTENCIAL_100000_ETRS89.shp")
potential_veg <- st_set_crs(potential_veg,25830)
plot(potential_veg["VEGETACION"])
sort(unique(potential_veg$VEGETACION))
veg.levels <- sort(unique(potential_veg$VEGETACION))

potential_veg$VEGETACION <- factor(potential_veg$VEGETACION, levels=veg.levels)

tm_shape(potential_veg, bbox=bboxplot)+
  tm_fill(col="VEGETACION",  
          title = "Potential Vegetation")+
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_shape(autonomias25830) +
  tm_borders() +
  tm_layout(legend.outside=TRUE,
            legend.width=2,
            legend.outside.position="right",
            legend.title.size = 1.5,
            legend.text.size = 1.1) +
  tm_shape(euskadi) +
  tm_borders(lwd=2, col="black") +
  tm_grid(lines=FALSE)


# Parent material ----------------------------------------------------

setwd("C:/Covariates/Euskadi/ParentMaterial/")

lithology_lv1 <- terra::rast("lithology_lv1.tif")
lithology_lv2 <- terra::rast("lithology_lv2.tif")
regolith_d <- terra::rast("regolith.tif")

### Plot
tm_shape(lithology_lv1, bbox=bboxplot)+
  tm_raster(style= "cat", 
            palette="Accent",n=6,
            labels=c("Superficial deposits",
                     "Sedimentary detrital rocks",
                     "Limestone, dolostone, marl, calcareous",
                      "Igneous rocks (volcanic rocks)",
                    "Igneous rocks",
                     "Igneous rocks (Granitoids)"),
            title="Lithology")+
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_shape(autonomias25830) +
  tm_borders() +
  tm_layout(legend.outside=TRUE,
            legend.outside.position="right") +
  tm_shape(euskadi) +
  tm_borders(lwd=2, col="black") +
  tm_grid(lines=FALSE)

tm_shape(lithology_lv1, bbox=bboxplot)+
  tm_raster(style= "cat", 
            palette="Set3",n=6,
            labels=c("Superficial deposits",
                     "Sedimentary detrital rocks",
                     "Limestone, dolostone, marl, calcareous",
                     "Igneous rocks (volcanic rocks)",
                     "Igneous rocks",
                     "Igneous rocks (Granitoids)"),
            title="Lithology")+
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_shape(autonomias25830) +
  tm_borders() +
  tm_layout(legend.outside=TRUE,
            legend.outside.position="right") +
  tm_shape(euskadi) +
  tm_borders(lwd=2, col="black") +
  tm_grid(lines=FALSE)

tm_shape(lithology_lv2, bbox=bboxplot)+
  tm_raster(style= "cat", 
            palette="Spectral",n=3,
            labels=c("Sedimentary",
                     "Metamorphic",
                     "Igneous"),
            title="Lithology")+
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_shape(autonomias25830) +
  tm_borders() +
  tm_layout(legend.outside=TRUE,
            legend.outside.position="right") +
  tm_shape(euskadi) +
  tm_borders(lwd=2, col="black") +
  tm_grid(lines=FALSE)


tm_shape(regolith_d, bbox=bboxplot)+
  tm_raster(style= "cat", 
            palette="YlOrBr",n=5,
            labels=c("0 - 0.5 m",
                     "0.5 - 1 m",
                     "1 - 2 m",
                     "2 - 4 m",
                     "> 4 m"),
            title="Regolith depth")+
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_shape(autonomias25830) +
  tm_borders() +
  tm_layout(legend.outside=TRUE,
            legend.outside.position="right") +
  tm_shape(euskadi) +
  tm_borders(lwd=2, col="black") +
  tm_grid(lines=FALSE)

### End of the script

### End of the script