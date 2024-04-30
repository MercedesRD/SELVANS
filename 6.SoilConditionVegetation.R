############################################################################################################################################
###  Method for optimizing pedogenon classes as soil districts for the Basque Country
###  in the context of the European soil Monitoring Law
###  In this script: Differences between vegetation groups

### Desired extent: Basque Country
### Resolution: 25m
### CRS: EPSG=25830

###  Author: Mercedes Roman Dobarco
###  Date: 15/02/2024

####### Load packages
### Spatial
library(sf)
library(terra)
library(gdalUtilities)
library(XML)
library(reproducible)
library(readxl)

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
library(Hmisc)
library(corrplot)

### Parallel computing
library(foreach)
library(parallel)
library(doParallel)

### clustering
library(ClusterR)

### Load information from previous script 5.PedogenonModelingForeach.R
setwd("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/6.SoilCondition")
load("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/6.SoilCondition/6.SoilCondition22032024.RData")

### Load stand information from BASONET plots

BASONET_spp_2001 <- read_excel("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Stand_data_Basonet/BASONET_spp.xlsx", 
                               sheet = "Plots2001")
View(BASONET_spp_2001)

BASONET_spp_2021 <- read_excel("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Stand_data_Basonet/BASONET_spp.xlsx", 
                               sheet = "Plots2021")
View(BASONET_spp_2021)

### Create unique ID including date
BASONET.PdGn$Basonet_Plot_Date <- paste0(BASONET.PdGn$myID,"_",BASONET.PdGn$Date)

### Same for stand data
BASONET_spp_2001$Basonet_Plot_Date <- paste0("BASONET_", BASONET_spp_2001$th, "_", BASONET_spp_2001$parc, "_2001")
BASONET_spp_2021$Basonet_Plot_Date <- paste0("BASONET_", BASONET_spp_2021$th, "_", BASONET_spp_2021$parc, "_2021")

### Subset columns
BASONET_spp_2001 <- BASONET_spp_2001[,c("Basonet_Plot_Date", "Clase","Sp1","Sp2","Sp3")]
BASONET_spp_2001$Date <- 2001
BASONET_spp_2021 <- BASONET_spp_2021[,c("Basonet_Plot_Date", "Clase","Sp1","Sp2","Sp3")]
BASONET_spp_2021$Date <- 2021

BASONET_spp <- rbind(BASONET_spp_2021,BASONET_spp_2001)

### What species are there?
sort(unique(c((BASONET_spp$Sp1),(BASONET_spp$Sp2),(BASONET_spp$Sp3))))

### how many plots the main species is Pinus radiata?
sort(unique(BASONET_spp[BASONET_spp$Sp1 == 28 & !is.na(BASONET_spp$Sp1),]$Basonet_Plot_Date)) ### 253 ()
sort(unique(BASONET_spp[BASONET_spp$Sp2 == 28 & !is.na(BASONET_spp$Sp2),]$Basonet_Plot_Date)) ### 17 ()

invasive_aloctonous <- c(7,11,92,207,292,307,392)
conifers_plantations <- c(28,18,30,33,34,35,235,335,435,319,229,230,17)
eucaliptus_plantations <- c(60,61,62,63,64,264,364,464)
conifers_native <- c(14,21,22,23,24,25,26,27,31,32,37,237,38,39,219,238)
broadleaved_native <- c(13,41,243,42,43,44,244,45,46,47,48,50,51,53,54,57,59,
                        257,357,457,557,657,757,857,957,55,255,56,256,356,245,
                        58,65,66,67,68,70,75,77,275,277,377,71,72,73,273,373,
                        74,1,80,82,84,86,87,88,89,253,268,293,294,858,81,83,283,
                        93,94,2,3,4,5,6,8,9,12,15,16,40,49,52,76,78,90,95,96,
                        97,98,99,215,276,278,291,295,297,299,315,355,376,378,
                        395,399,415,476,478,499,576,578,676,69,469,495,
                        279,281,282,289,389,455,456,489,515,569,595,599,776,778,258,678,79, 901,902)
#broadleaved_plantation <- c(,)


BASONET_spp$VegetationSp1 <- ifelse(BASONET_spp$Sp1 %in% invasive_aloctonous, "Invasive",
                                  ifelse(BASONET_spp$Sp1 %in% conifers_plantations, "Conifer plantation",
                                         ifelse(BASONET_spp$Sp1 %in% eucaliptus_plantations, "Eucalyptus",
                                                ifelse(BASONET_spp$Sp1 %in% conifers_native, "Conifers",
                                                       ifelse(BASONET_spp$Sp1 %in% broadleaved_native, "Broadleaved", BASONET_spp$Sp1)))))
BASONET_spp$VegetationSp2 <- ifelse(BASONET_spp$Sp2 %in% invasive_aloctonous, "Invasive",
                                    ifelse(BASONET_spp$Sp2 %in% conifers_plantations, "Conifer plantation",
                                           ifelse(BASONET_spp$Sp2 %in% eucaliptus_plantations, "Eucalyptus",
                                                  ifelse(BASONET_spp$Sp2 %in% conifers_native, "Conifers",
                                                         ifelse(BASONET_spp$Sp2 %in% broadleaved_native, "Broadleaved", NA)))))

table(BASONET_spp$VegetationSp1,BASONET_spp$VegetationSp2)
table(BASONET_spp$VegetationSp2)
BASONET_spp$Vegetation <- NA

BASONET_spp[((!is.na(BASONET_spp$VegetationSp1) &
                !is.na(BASONET_spp$VegetationSp2) &
               BASONET_spp$VegetationSp1 == "Broadleaved" & 
               BASONET_spp$VegetationSp2 %in% c("Broadleaved","Eucalyptus","Invasive")))|
              (!is.na(BASONET_spp$VegetationSp1) &
                 BASONET_spp$VegetationSp1 == "Broadleaved" & 
                is.na(BASONET_spp$VegetationSp2)),]$Vegetation <- "Broadleaved"

BASONET_spp[((!is.na(BASONET_spp$VegetationSp1) &
                !is.na(BASONET_spp$VegetationSp2) &
                BASONET_spp$VegetationSp1 == "Broadleaved" & 
                BASONET_spp$VegetationSp2 %in% c("Conifers"))|
               (!is.na(BASONET_spp$VegetationSp1) &
                  !is.na(BASONET_spp$VegetationSp2) &
                  BASONET_spp$VegetationSp2 %in% c("Broadleaved") & 
                  BASONET_spp$VegetationSp1 %in% c("Conifers"))|
               (!is.na(BASONET_spp$VegetationSp1) &
                  !is.na(BASONET_spp$VegetationSp2) &
                  BASONET_spp$VegetationSp2 %in% c("Invasive", "Eucalyptus") & 
                  BASONET_spp$VegetationSp1 %in% c("Conifers"))),]$Vegetation <- "Mixed"

BASONET_spp[((!is.na(BASONET_spp$VegetationSp1) &
                !is.na(BASONET_spp$VegetationSp2) &
                BASONET_spp$VegetationSp1 == "Broadleaved" & 
                BASONET_spp$VegetationSp2 %in% c( "Conifer plantation"))|
               (!is.na(BASONET_spp$VegetationSp1) &
                  !is.na(BASONET_spp$VegetationSp2) &
                  BASONET_spp$VegetationSp2 %in% c("Broadleaved") & 
                  BASONET_spp$VegetationSp1 %in% c("Conifer plantation"))),]$Vegetation <- "Mixed plantation"

BASONET_spp[(!is.na(BASONET_spp$VegetationSp1) &
               !is.na(BASONET_spp$VegetationSp2) &
               BASONET_spp$VegetationSp1 == "Conifers" & 
               BASONET_spp$VegetationSp2 == "Conifers")|
              (!is.na(BASONET_spp$VegetationSp1) &
                 BASONET_spp$VegetationSp1 == "Conifers" & 
                 is.na(BASONET_spp$VegetationSp2)), ]$Vegetation <- "Conifers"

BASONET_spp[(!is.na(BASONET_spp$VegetationSp1) &
               !is.na(BASONET_spp$VegetationSp2) &
               BASONET_spp$VegetationSp1 == "Conifer plantation" & 
               BASONET_spp$VegetationSp2 == "Conifer plantation")|
              (!is.na(BASONET_spp$VegetationSp1) &
                 BASONET_spp$VegetationSp1 == "Conifer plantation" & 
                 is.na(BASONET_spp$VegetationSp2))|
              (!is.na(BASONET_spp$VegetationSp1) &
                 !is.na(BASONET_spp$VegetationSp2) &
                 BASONET_spp$VegetationSp1 == "Conifer plantation" & 
                 BASONET_spp$VegetationSp2 == "Conifers") |
              (!is.na(BASONET_spp$VegetationSp1) &
                 !is.na(BASONET_spp$VegetationSp2) &
                 BASONET_spp$VegetationSp2 == "Conifer plantation" & 
                 BASONET_spp$VegetationSp1 == "Conifers")|
              (!is.na(BASONET_spp$VegetationSp1) &
                 !is.na(BASONET_spp$VegetationSp2) &
                 BASONET_spp$VegetationSp1 %in% c("Conifer plantation") & 
                 BASONET_spp$VegetationSp2 %in% c("Eucalyptus", "Invasive")),]$Vegetation <- "Conifer plantation"

BASONET_spp[(!is.na(BASONET_spp$VegetationSp1) &
               !is.na(BASONET_spp$VegetationSp2) &
               BASONET_spp$VegetationSp1 == "Eucalyptus" & 
               BASONET_spp$VegetationSp2 == "Eucalyptus")|
              (!is.na(BASONET_spp$VegetationSp1) &
              !is.na(BASONET_spp$VegetationSp2) &
              BASONET_spp$VegetationSp1 == "Eucalyptus" & 
              BASONET_spp$VegetationSp2 == "Broadleaved")|
              (!is.na(BASONET_spp$VegetationSp1) &
                 !is.na(BASONET_spp$VegetationSp2) &
                 BASONET_spp$VegetationSp1 == "Eucalyptus" & 
                 BASONET_spp$VegetationSp2 == "Invasive")|
              (!is.na(BASONET_spp$VegetationSp1) &
                 !is.na(BASONET_spp$VegetationSp2) &
                 BASONET_spp$VegetationSp1 == "Eucalyptus" & 
                 BASONET_spp$VegetationSp2 == "Conifer plantation")|
              (!is.na(BASONET_spp$VegetationSp1) &
                 !is.na(BASONET_spp$VegetationSp2) &
                 BASONET_spp$VegetationSp1 == "Eucalyptus" & 
                 BASONET_spp$VegetationSp2 == "Conifers")|
              (!is.na(BASONET_spp$VegetationSp1) &
                 BASONET_spp$VegetationSp1 == "Eucalyptus" & 
                 is.na(BASONET_spp$VegetationSp2)),]$Vegetation <- "Eucalyptus"

BASONET_spp[(!is.na(BASONET_spp$VegetationSp1) &
               BASONET_spp$VegetationSp1 == "Invasive" & 
               BASONET_spp$VegetationSp2 %in% c("Invasive","Conifer plantation", "Broadleaved"))|
              (!is.na(BASONET_spp$VegetationSp1) &
                 BASONET_spp$VegetationSp1 == "Invasive" & 
                 is.na(BASONET_spp$VegetationSp2)),]$Vegetation <- "Invasive"

table(BASONET_spp$Vegetation)

BASONET_spp$Vegetation <- factor(BASONET_spp$Vegetation,
                                 levels=c("Broadleaved", "Mixed", "Conifers",
                                          "Mixed plantation","Conifer plantation",
                                          "Eucalyptus","Invasive"))

### Join soil and stand data
BASONET.PdGn <- left_join(BASONET.PdGn,BASONET_spp)


ggplot(data=BASONET.PdGn[!is.na(BASONET.PdGn$TOC) & 
                           !is.na(BASONET.PdGn$Clay) &
                           !is.na(BASONET.PdGn$PdGn) &
                           !is.na(BASONET.PdGn$Vegetation)&
                           BASONET.PdGn$Date ==2021,]) +  
  geom_boxplot(aes(x = Vegetation, y = TOC, fill = Vegetation)) +
  ylab(label="TOC") +
  #xlab(label="Soil district") +
  theme_bw() +
  # scale_fill_manual(values= map.PdGn.Comb26.k9$branch.centroids.ord$colors,
  #                   name="Soil district")+ 
  facet_wrap(~ Depth_interval) +
  geom_hline(yintercept =1/13,
             color = "red",
             lwd = 1,  
             linetype = "dashed") #+
 # ylim(0,0.6)


ggplot(data=BASONET.PdGn[!is.na(BASONET.PdGn$TOC) & 
                           !is.na(BASONET.PdGn$Clay) &
                           !is.na(BASONET.PdGn$PdGn) &
                           !is.na(BASONET.PdGn$Vegetation)&
                           BASONET.PdGn$Date ==2021,]) +  
  geom_boxplot(aes(x = Vegetation,
                   y = TOC_Clay, 
                   fill = as.factor(PdGn))) +
  ylab(label="TOC / Clay") +
  xlab(label="Vegetation") +
  theme_bw() +
  scale_fill_manual(values= map.PdGn.Comb26.k9$branch.centroids.ord$colors,
                    name="Soil district") +
  geom_hline(yintercept =1/13,
             color = "red",
             lwd = 1,  
             linetype = "dashed")+
  ylim(0,0.4)

ggplot(data=BASONET.PdGn[!is.na(BASONET.PdGn$pH) & 
                           
                           !is.na(BASONET.PdGn$PdGn) &
                           !is.na(BASONET.PdGn$Vegetation),]) +  
  geom_boxplot(aes(x = as.factor(Date), y= pH, fill=Vegetation)) +
  ylab(label="pH") +
  xlab(label="Soil district") +
  theme_bw() +
  # scale_fill_manual(values= map.PdGn.Comb26.k9$branch.centroids.ord$colors,
  #                   name="Soil district")+ 
  facet_wrap(~ Depth_interval) 

ggplot(data=BASONET.PdGn[!is.na(BASONET.PdGn$EC) & 
                           !is.na(BASONET.PdGn$PdGn) &
                           !is.na(BASONET.PdGn$Vegetation),]) +  
  geom_boxplot(aes(x = Vegetation, y= EC))+
  ylab(label="EC") +
  xlab(label="Soil district") +
  theme_bw() #+
  # scale_fill_manual(values= map.PdGn.Comb26.k9$branch.centroids.ord$colors,
  #                   name="Soil district")+
  #facet_wrap(~ PdGn) 


### Bulk density
BASONET.SubsoilCompaction.PdGn$Basonet_Plot_Date <- paste0(BASONET.SubsoilCompaction.PdGn$myID,"_",BASONET.SubsoilCompaction.PdGn$Date)
BASONET.SubsoilCompaction.PdGn <- left_join(BASONET.SubsoilCompaction.PdGn, BASONET_spp)

ggplot(data=BASONET.SubsoilCompaction.PdGn[!is.na(BASONET.SubsoilCompaction.PdGn$Vegetation),]) +  
  geom_boxplot(aes(x =Vegetation ,
                   y= Bulk_Density, 
                   fill=as.factor(Condition))) +
  ylab(label="Bulk density") +
  xlab(label="Vegetation") +
  theme_bw()




# 3. PERMANOVA -----------------------------------------------------------

BASONET.PdGn.subset <- BASONET.PdGn[complete.cases(BASONET.PdGn[,c("PdGn","Vegetation","Depth_interval",
                                                                   "Silt","Clay","pH","TOC","TN_ppm","K_ppm",
                                                                   "Ca_ppm", "Mg_ppm","CECef")]) &
                                      BASONET.PdGn$Date == 2021, ]

Response <- BASONET.PdGn.subset[,c("Silt","Clay","pH","TOC",
                                   "TN_ppm","K_ppm","Ca_ppm", "Mg_ppm","CECef")]

#Response <- scale(Response,center=TRUE, scale=TRUE)

### PERMANOVA
library(vegan)
set.seed(1984)
permanova.BASONET <- adonis2(Response ~ PdGn*Vegetation+Depth_interval,
                             data=BASONET.PdGn.subset,
                             method = "mahalanobis",
                             permutations=9999)
permanova.BASONET
densityplot(permustats(permanova.BASONET))

plot(varpart(Y=Response, ~ Vegetation , ~ PdGn ,
             data = BASONET.PdGn.subset))

BASONET.dist <- vegdist(x = Response, method = "mahalanobis")

set.seed(1984)
BASONET_disper <- betadisper(d =BASONET.dist,
                             group=BASONET.PdGn.subset$Vegetation)
plot(BASONET_disper)


###############################################################################################