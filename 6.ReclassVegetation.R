############################################################################################################################################
###  Method for optimizing pedogenon classes as soil districts for the Basque Country
###  in the context of the European soil Monitoring Law
###  In this script: Differences between vegetation groups
###  Calculate basal area by species (m2/ha) and percentage of total BA
###  Calculate basal area by vegetation group and percentage of total BA
###  Identify dominant species by plot
###  Identify dominant vegetation group (define criteria for mixed)
###  Designate as managed / not managed

### To do: Calculate stand density index by dominant species

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

### Create unique ID including date for the Soil and Pedogenon data
BASONET.PdGn$Basonet_Plot_Date <- paste0(BASONET.PdGn$myID,"_",BASONET.PdGn$Date)


# 1.1 Load 2001 data -------------------------------------------------------

### Load stand information from BASONET plots
BASONET_dasom_2001 <- read_excel("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Stand_data_Basonet/BasonetDasomN.xlsx", 
                                 sheet = "pies2001")
### check with the species summary file
BASONET_spp_2001 <- read_excel("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Stand_data_Basonet/BasonetDasomN.xlsx", 
                               sheet = "parcelas2001")
### Create unique ID including date 
BASONET_spp_2001$Basonet_Plot_Date <- paste0("BASONET_", BASONET_spp_2001$th, "_", BASONET_spp_2001$parc, "_2001")

### And Create unique ID including date for 2001 stand data
BASONET_dasom_2001$Basonet_Plot_Date <- paste0("BASONET_", BASONET_dasom_2001$th, "_", BASONET_dasom_2001$parc, "_2001")

setdiff(sort(unique(BASONET_dasom_2001$Basonet_Plot_Date)),
        (sort(unique(BASONET_spp_2001$Basonet_Plot_Date))))

setdiff(sort(unique(BASONET_dasom_2001$Basonet_Plot_Date)),
        (sort(unique(BASONET.PdGn$Basonet_Plot_Date)))) ### I don´t have soil data from 21 plots

setdiff(sort(unique(BASONET.PdGn$Basonet_Plot_Date)),
        sort(unique(BASONET_dasom_2001$Basonet_Plot_Date))) 

setdiff(sort(unique(BASONET_spp_2001$Basonet_Plot_Date)), sort(unique(BASONET_dasom_2001$Basonet_Plot_Date)))
setdiff(sort(unique(BASONET_spp_2001$parc)), sort(unique(BASONET_dasom_2001$parc)))
length(setdiff(sort(unique(BASONET_spp_2001$Basonet_Plot_Date)), sort(unique(BASONET_dasom_2001$Basonet_Plot_Date))))
### 63 plots are many

### Subset columns
BASONET_spp_2001 <- BASONET_spp_2001[,c("th","parc","Basonet_Plot_Date", "Clase","Sp1","Sp2","Sp3")]
BASONET_spp_2001$Date <- 2001


# 1.2 Calculate Basal area 2001 by plot ---------------------------------------

#BASONET_dasom_2001$Diam_ave <- NA
summary(BASONET_dasom_2001$Diam1_mm)
summary(BASONET_dasom_2001$Diam2_mm)
BASONET_dasom_2001$Diam_ave_cm <- (BASONET_dasom_2001$Diam1_mm + BASONET_dasom_2001$Diam2_mm)/20
hist(BASONET_dasom_2001$Diam_ave_cm, breaks=30)
BASONET_dasom_2001$BArea_m2 <- pi*((BASONET_dasom_2001$Diam_ave_cm/200)^2)
Area_plot <- pi*(25^2)
Area_plot_ha <- Area_plot/10000

### What are the quantiles?
quantile(BASONET_dasom_2001$Diam_ave_cm, p=c(0.25,0.75,seq(from=0.1, to=0.9, by=0.1)))
#    10%    20%    30%    40%    50%    60%    70%    80%    90% 
# 12.400 15.000 17.600 21.450 25.100 29.050 34.305 42.950 50.685
summary(BASONET_dasom_2001$Diam_ave_cm)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7.50   16.20   25.10   28.73   38.09  140.00 

### Add species names
spp_names <- read_excel("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Stand_data_Basonet/spp_names.xlsx")
View(spp_names)

spp_names <- distinct(spp_names)
spp_names <- spp_names[spp_names$Spx != 10, ]
spp_names[duplicated(spp_names$Spx),]
spp_names[spp_names$Spx == 45,]$Name <- "Quercus ilex"
spp_names[spp_names$Spx == 92,]$Name <- "Robinia pseudoacacia"
spp_names[spp_names$Spx == 243,]$Name <-  "Quercus humilis" 
spp_names <- distinct(spp_names)

### Attach species
BASONET_dasom_2001 <- merge(BASONET_dasom_2001,spp_names, by.x="Species",by.y="Spx")

### Summarise total BA (m2/ha) and QMD (cm) by plot
BASONET_dasom_2001_summary_plot <- BASONET_dasom_2001 %>%
  group_by(., Basonet_Plot_Date) %>%
  summarise(BA_Total_m2ha = sum(BArea_m2, na.rm = TRUE)/Area_plot_ha,
            QMD_cm = sqrt(sum(Diam_ave_cm^2, na.rm=TRUE)/n()))

View(BASONET_dasom_2001_summary_plot)
hist(BASONET_dasom_2001_summary_plot$BA_Total_m2ha, breaks=30)
hist(BASONET_dasom_2001_summary_plot$QMD_cm, breaks=30)
plot(BASONET_dasom_2001_summary_plot$QMD_cm, BASONET_dasom_2001_summary_plot$BA_Total_m2ha)

### Auxiliary functions to apply by plot (subset data for plot outside function)

### Function to calculate Stand Density Index by the summation method
### Long, J.N., Daniel, T.W., 1990. Assessment of growing stock in uneven-aged stands. West. J. Appl. For. 5, 93–96.
### Woodall, C. W., Miles, P. D., & Vissage, J. S. (2005). Determining maximum stand density index in mixed species stands
### for strategic-scale stocking assessments. Forest Ecology and Management, 216(1-3), 367-377.

SDI_func <- function(df, Diam_ave_cm) {
  
  ### Assign diameter class
  df$D_class <- ifelse(df$Diam_ave_cm <= 10,"0-10", 
                       ifelse((df$Diam_ave_cm > 10 & df$Diam_ave_cm <= 20),"10-20",
                              ifelse((df$Diam_ave_cm > 20 & df$Diam_ave_cm <= 30),"20-30",
                                     ifelse((df$Diam_ave_cm > 30 & df$Diam_ave_cm <= 40),"30-40",
                                            ifelse((df$Diam_ave_cm > 40 & df$Diam_ave_cm <= 50), "40-50",
                                                   ifelse((df$Diam_ave_cm > 50 & df$Diam_ave_cm <= 60),"50-60",
                                                          ifelse((df$Diam_ave_cm > 60 & df$Diam_ave_cm <= 70), "60-70",
                                                                 ifelse((df$Diam_ave_cm > 70 & df$Diam_ave_cm <= 80),"70-80",
                                                                        ifelse((df$Diam_ave_cm > 80 & df$Diam_ave_cm <= 90), "80-90",
                                                                                ifelse((df$Diam_ave_cm > 90 & df$Diam_ave_cm <= 100),"90-100",
                                                                                       ifelse((df$Diam_ave_cm > 100 & df$Diam_ave_cm <=110),"100-110",
                                                                                              ifelse((df$Diam_ave_cm > 110 & df$Diam_ave_cm <= 120),"110-120",
                                                                                                     ifelse((df$Diam_ave_cm > 120 & df$Diam_ave_cm <= 130),"120-130",
                                                                                                            ifelse((df$Diam_ave_cm > 130 & df$Diam_ave_cm <= 140),"130-140",
                                                                                                                   ifelse((df$Diam_ave_cm > 140 & df$Diam_ave_cm <= 150),"140-150", NA )))))))))))))))
    
  ### Assign midpoint in cm
  df$DBH <- ifelse(df$Diam_ave_cm <= 10,5, 
                       ifelse((df$Diam_ave_cm > 10 & df$Diam_ave_cm <= 20),15,
                              ifelse((df$Diam_ave_cm > 20 & df$Diam_ave_cm <= 30),25,
                                     ifelse((df$Diam_ave_cm > 30 & df$Diam_ave_cm <= 40),35,
                                            ifelse((df$Diam_ave_cm > 40 & df$Diam_ave_cm <= 50), 45,
                                                   ifelse((df$Diam_ave_cm > 50 & df$Diam_ave_cm <= 60),55,
                                                          ifelse((df$Diam_ave_cm > 60 & df$Diam_ave_cm <= 70),65,
                                                                 ifelse((df$Diam_ave_cm > 70 & df$Diam_ave_cm <= 80),75,
                                                                        ifelse((df$Diam_ave_cm > 80 & df$Diam_ave_cm <= 90), 85,
                                                                               ifelse((df$Diam_ave_cm > 90 & df$Diam_ave_cm <= 100),95,
                                                                                      ifelse((df$Diam_ave_cm > 100 & df$Diam_ave_cm <=110),105,
                                                                                             ifelse((df$Diam_ave_cm > 110 & df$Diam_ave_cm <= 120),115,
                                                                                                    ifelse((df$Diam_ave_cm > 120 & df$Diam_ave_cm <= 130),125,
                                                                                                           ifelse((df$Diam_ave_cm > 130 & df$Diam_ave_cm <= 140),135,
                                                                                                                  ifelse((df$Diam_ave_cm > 140 & df$Diam_ave_cm <= 150),145, NA )))))))))))))))
  
  
  ### Count number of trees per diameter class
  df_Dclass <- df %>%
    group_by(., D_class, DBH) %>%
    summarise(tph = n()/Area_plot_ha) ### Calculate trees per hectare
  
  SDI <- sum(((df_Dclass$DBH/25)^1.6) * df_Dclass$tph, na.rm=TRUE)
  
  return(SDI)
  
}

Spp_BA_rank <- function(df, Species, Diam_ave_cm, BArea_m2, Area_plot_ha){
  
  
  #df <- BASONET_dasom_2001[BASONET_dasom_2001$Basonet_Plot_Date == unique_plots[4],]
  ### IF colname is not called Species, change the name
  colnames(df)[colnames(df)== Species] <- "Species"
  colnames(df)[colnames(df)== Diam_ave_cm] <- "Diam_ave_cm"
  
  ### In case BA is not calculated
  if(isFALSE(BArea_m2 %in% colnames(df))) {
    df$BArea_m2 <- pi*((df$Diam_ave_cm/200)^2)
  } else {
    colnames(df)[colnames(df)== BArea_m2] <- "BArea_m2"
  }
  
  BA_Total <- sum(df$BArea_m2,na.rm = TRUE)
  
  SDI_Total <- SDI_func(df = df, Diam_ave_cm = "Diam_ave_cm" )
  
  ### Summarise by species
  df_summary <- df %>%
    group_by(., Species) %>%
    summarise(BA_m2ha = sum(BArea_m2,na.rm = TRUE)/Area_plot_ha,
              BA_p = sum(BArea_m2,na.rm = TRUE)/BA_Total*100,
              QMD_cm = sqrt(sum(Diam_ave_cm^2, na.rm=TRUE)/n()) #,
              # Height_mean = mean(Height_m, na.rm = TRUE),
              # Height_median = median(Height_m, na.rm = TRUE),
              # Height_p05 = quantile(Height_m, probs=0.05, na.rm = TRUE),
              # Height_p95 = quantile(Height_m, probs=0.95, na.rm = TRUE),
              # Quality_mean = mean(Quality, na.rm = TRUE)
              )  
  
  ### Attach species
  df_summary <- merge(df_summary,spp_names, by.x="Species",by.y="Spx")
  df_summary <- arrange(df_summary,desc(BA_m2ha)) ### Sort by dominant species
  
  # Number of species in that plot
  n.Spp <- nrow(df_summary)
  
  ### Create empty dataframe for output
  df.out <- data.frame(matrix(ncol = 18, nrow = 1, data =NA))
  
  #provide column names
  colnames(df.out) <- c("BASp1_m2ha","BASp2_m2ha","BASp3_m2ha",
                        "BASp1_p","BASp2_p","BASp3_p",
                        "QMD_Sp1","QMD_Sp2","QMD_Sp3",
                        "SDI_Sp1","SDI_Sp2","SDI_Sp3",
                        "SDI_Sp1_p","SDI_Sp2_p","SDI_Sp3_p",
                        "Sp1", "Sp2", "Sp3")
  if(n.Spp == 1) {
    
    df.out$BASp1_m2ha <- df_summary[1,]$BA_m2ha
    df.out$BASp1_p <- df_summary[1,]$BA_p 
    df.out$QMD_Sp1 <- df_summary[1,]$QMD_cm 
    df.out$Sp1 <- df_summary[1,]$Name
    df.out$SDI_Sp1 <- SDI_func(df = df[df$Species == df_summary[1,]$Species,], Diam_ave_cm = "Diam_ave_cm")
    df.out$SDI_Sp1_p <- df.out$SDI_Sp1 / SDI_Total * 100
    
  } else if (n.Spp == 2) {
    
    df.out$BASp1_m2ha <- df_summary[1,]$BA_m2ha
    df.out$BASp1_p <- df_summary[1,]$BA_p 
    df.out$QMD_Sp1 <- df_summary[1,]$QMD_cm 
    df.out$Sp1 <- df_summary[1,]$Name
    df.out$SDI_Sp1 <- SDI_func(df = df[df$Species == df_summary[1,]$Species,], Diam_ave_cm = "Diam_ave_cm")
    df.out$SDI_Sp1_p <- df.out$SDI_Sp1 / SDI_Total * 100
    
    df.out$BASp2_m2ha <- df_summary[2,]$BA_m2ha
    df.out$BASp2_p <- df_summary[2,]$BA_p 
    df.out$QMD_Sp2 <- df_summary[2,]$QMD_cm
    df.out$Sp2 <- df_summary[2,]$Name
    df.out$SDI_Sp2 <- SDI_func(df = df[df$Species == df_summary[2,]$Species,], Diam_ave_cm = "Diam_ave_cm")
    df.out$SDI_Sp2_p <- df.out$SDI_Sp2 / SDI_Total * 100
    
    
  } else if (n.Spp >= 3) {
    
    df.out$BASp1_m2ha <- df_summary[1,]$BA_m2ha
    df.out$BASp1_p <- df_summary[1,]$BA_p 
    df.out$QMD_Sp1 <- df_summary[1,]$QMD_cm 
    df.out$Sp1 <- df_summary[1,]$Name
    df.out$SDI_Sp1 <- SDI_func(df = df[df$Species == df_summary[1,]$Species,], Diam_ave_cm = "Diam_ave_cm")
    df.out$SDI_Sp1_p <- df.out$SDI_Sp1 / SDI_Total * 100
    
    df.out$BASp2_m2ha <- df_summary[2,]$BA_m2ha
    df.out$BASp2_p <- df_summary[2,]$BA_p 
    df.out$QMD_Sp2 <- df_summary[2,]$QMD_cm
    df.out$Sp2 <- df_summary[2,]$Name
    df.out$SDI_Sp2 <- SDI_func(df = df[df$Species == df_summary[2,]$Species,], Diam_ave_cm = "Diam_ave_cm")
    df.out$SDI_Sp2_p <- df.out$SDI_Sp2 / SDI_Total * 100
    
    df.out$BASp3_m2ha <- df_summary[3,]$BA_m2ha
    df.out$BASp3_p <- df_summary[3,]$BA_p 
    df.out$QMD_Sp3 <- df_summary[3,]$QMD_cm
    df.out$Sp3 <- df_summary[3,]$Name
    df.out$SDI_Sp3 <- SDI_func(df = df[df$Species == df_summary[3,]$Species,], Diam_ave_cm = "Diam_ave_cm")
    df.out$SDI_Sp3_p <- df.out$SDI_Sp3 / SDI_Total * 100
    
  }
  
  df.out$SDI <- SDI_Total
  return(df.out)
  
  }
 
#BASONET_dasom_2001[BASONET_dasom_2001$Species == 45,]$Name <- "Quercus ilex"
unique_plots <- sort(unique(BASONET_dasom_2001_summary_plot$Basonet_Plot_Date))

for (i in 1:length(unique_plots)) {
  print(i)
  if(i == 1) {
    ### Select the rows for this plot
    plot.i <- BASONET_dasom_2001[BASONET_dasom_2001$Basonet_Plot_Date == unique_plots[i],]
    plot.i.out <- Spp_BA_rank(df = plot.i, 
                              Species = "Species",
                              Diam_ave_cm = "Diam_ave_cm",
                              BArea_m2 = "BArea_m2", 
                              Area_plot_ha = Area_plot_ha)
    ### Attach Plot ID
    plot.i.out$Basonet_Plot_Date <- unique_plots[i]
    } else if (i > 1) {
    ### Select the rows for this plot
    plot.j <- BASONET_dasom_2001[BASONET_dasom_2001$Basonet_Plot_Date == unique_plots[i],]
    plot.j.out <- Spp_BA_rank(df = plot.j, 
                              Species = "Species",
                              Diam_ave_cm = "Diam_ave_cm",
                              BArea_m2 = "BArea_m2", 
                              Area_plot_ha = Area_plot_ha)
    ### Attach Plot ID
    plot.j.out$Basonet_Plot_Date <- unique_plots[i]
    ### Bind to previous plot
    plot.i.out <- rbind(plot.i.out,plot.j.out)
    }
}

### Now JOIN with the Plot summary
BASONET_dasom_2001_summary_plot <- merge(BASONET_dasom_2001_summary_plot,plot.i.out, by = "Basonet_Plot_Date")


# 1.3 Assign vegetation group to 2001 data ------------------------------------

### Assign species to vegetation groups
#   8  15  16  18  21  24  25  26  28  33  34  35  37  41  42  43  44  45  48  51  54 
#  55  56  57  58  60  61  64  65  68
#  71  72  73  74  75  76  77  78  79  81 
#  92  95  97  99 215 255 256 278 373 578 657

invasive             <- c(7,11,92,207,292,307,392) ### We have 92 (Robinia pseudoacacia) in 5 plots
conifers_plantations <- c(28,18,30,33,34,35,235,335,435,319,229,230,17)
eucalyptus           <- c(60,61,62,63,64,264,364,464)
conifers_native      <- c(14,21,22,23,24,25,26,27,31,32,37,237,38,39,219,238)
broadleaved          <- c(13,41,243,42,43,44,244,45,46,47,48,50,51,53,54,57,59,
                        257,357,457,557,657,757,857,957,55,255,56,256,356,245,
                        58,65,66,67,68,70,75,77,275,277,377,71,72,73,273,373,
                        74,1,80,82,84,86,87,88,89,253,268,293,294,858,81,83,283,
                        93,94,2,3,4,5,6,8,9,12,15,16,40,49,52,76,78,90,95,96,
                        97,98,99,215,276,278,291,295,297,299,315,355,376,378,
                        395,399,415,476,478,499,576,578,676,69,469,495,
                        279,281,282,289,389,455,456,489,515,569,595,599,776,778,
                        258,678,79,901,902)

BASONET_dasom_2001$Vegetation <- ifelse(BASONET_dasom_2001$Species %in% invasive, "Invasive",
                                                ifelse(BASONET_dasom_2001$Species %in% conifers_plantations, "Conifer_plantation",
                                                       ifelse(BASONET_dasom_2001$Species %in% eucalyptus, "Eucalyptus",
                                                              ifelse(BASONET_dasom_2001$Species %in% conifers_native, "Conifer",
                                                                     ifelse(BASONET_dasom_2001$Species %in% broadleaved, "Broadleaved", 
                                                                            BASONET_dasom_2001$Species)))))



### Auxiliary function to calculate dominant Vegetation group

VegGroup_BA_rank <- function(df, Species, Area_plot_ha){
  
  BA_Total <- sum(df$BArea_m2,na.rm = TRUE)
  
  ### Summarise by Group
  df_summary <- df %>%
    group_by(., Vegetation) %>%
    summarise(BA_m2ha = sum(BArea_m2,na.rm = TRUE)/Area_plot_ha,
              BA_p = sum(BArea_m2,na.rm = TRUE)/BA_Total*100) %>%
    arrange(., desc(BA_m2ha))

  # Number of vegetation groups in that plot
  n.Groups <- nrow(df_summary)
  
  ### Create empty dataframe for output
  df.out <- data.frame(matrix(ncol = 10, nrow = 1, data =NA))
  
  #provide column names
  colnames(df.out) <- c("BAg1_m2ha","BAg2_m2ha","BAg3_m2ha",
                        "BAg1_p","BAg2_p","BAg3_p",
                        "Group1", "Group2", "Group3", "DomGroup")
  if(n.Groups == 1) {
    
    df.out$BAg1_m2ha <- df_summary[1,]$BA_m2ha
    df.out$BAg1_p <- df_summary[1,]$BA_p 
    df.out$Group1 <- df_summary[1,]$Vegetation
    df.out$DomGroup <- df_summary[1,]$Vegetation
    
  } else if (n.Groups == 2) {
    
    df.out$BAg1_m2ha <- df_summary[1,]$BA_m2ha
    df.out$BAg1_p <- df_summary[1,]$BA_p 
    df.out$Group1 <- df_summary[1,]$Vegetation
    
    df.out$BAg2_m2ha <- df_summary[2,]$BA_m2ha
    df.out$BAg2_p <- df_summary[2,]$BA_p 
    df.out$Group2 <- df_summary[2,]$Vegetation
    
    
  } else if (n.Groups >= 3) {
    
    df.out$BAg1_m2ha <- df_summary[1,]$BA_m2ha
    df.out$BAg1_p <- df_summary[1,]$BA_p 
    df.out$Group1 <- df_summary[1,]$Vegetation
    
    df.out$BAg2_m2ha <- df_summary[2,]$BA_m2ha
    df.out$BAg2_p <- df_summary[2,]$BA_p 
    df.out$Group2 <- df_summary[2,]$Vegetation
    
    df.out$BAg3_m2ha <- df_summary[3,]$BA_m2ha
    df.out$BAg3_p <- df_summary[3,]$BA_p 
    df.out$Group3 <- df_summary[3,]$Vegetation
    
  }
  
  ### If the dominant group has > 75% BA, then classify dominant vegetation as that group.
  df.out$DomGroup <- ifelse(df.out$BAg1_p >= 75, df.out$Group1, NA)
  
  ### Conditions for diagonal cases
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == df.out$Group2, df.out$Group1, df.out$DomGroup)
  
  ### Conditions for Group1 == Broadleaved
  ### "Broadleaved"        "Conifer_plantation" "Conifer"            "Eucalyptus"         "Invasive"   
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Broadleaved" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75 &
                              df.out$Group2 == "Conifer", 
                            "Mixed",
                            df.out$DomGroup)
  
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Broadleaved" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75 &
                              df.out$Group2 == "Conifer_plantation", 
                            "Mixed_plantation",
                            df.out$DomGroup)
  
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Broadleaved" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75 &
                              df.out$Group2 == "Eucalyptus", 
                            "Broadleaved_Eucalyptus",
                            df.out$DomGroup)
  
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Broadleaved" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75 &
                              df.out$Group2 == "Invasive", 
                            "Broadleaved",
                            df.out$DomGroup)
  
  ### Conditions for Group1 == Conifer
  ### "Broadleaved"        "Conifer_plantation" "Conifer"            "Eucalyptus"         "Invasive"   
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Conifer" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75 &
                              df.out$Group2 == "Broadleaved", 
                            "Mixed",
                            df.out$DomGroup)
  
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Conifer" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75 &
                              df.out$Group2 == "Conifer_plantation", 
                            "Conifer_plantation",
                            df.out$DomGroup)
  
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Conifer" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75 &
                              df.out$Group2 == "Eucalyptus", 
                            "Conifer_plantation",
                            df.out$DomGroup)
  
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Conifer" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75 &
                              df.out$Group2 == "Invasive", 
                            "Conifer",
                            df.out$DomGroup)
  
  ### Conditions for Group1 == Conifer_plantation
  ### "Broadleaved"        "Conifer_plantation" "Conifer"            "Eucalyptus"         "Invasive"   
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Conifer_plantation" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75 &
                              df.out$Group2 == "Broadleaved", 
                            "Mixed_plantation",
                            df.out$DomGroup)
  
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Conifer_plantation" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75 &
                              df.out$Group2 == "Conifer", 
                            "Conifer_plantation",
                            df.out$DomGroup)
  
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Conifer_plantation" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75 &
                              df.out$Group2 == "Eucalyptus", 
                            "Conifer_plantation",
                            df.out$DomGroup)
  
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Conifer_plantation" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75 &
                              df.out$Group2 == "Invasive", 
                            "Conifer_plantation",
                            df.out$DomGroup)
  

  ### Conditions for Group1 == Eucalyptus
  ### "Broadleaved"        "Conifer_plantation" "Conifer"            "Eucalyptus"         "Invasive"   
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Eucalyptus" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75 &
                              df.out$Group2 == "Broadleaved", 
                            "Broadleaved_Eucalyptus",
                            df.out$DomGroup)
  
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Eucalyptus" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75 &
                              df.out$Group2 == "Conifer_plantation", 
                            "Conifer_plantation",
                            df.out$DomGroup)
  
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Eucalyptus" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75 &
                              df.out$Group2 == "Conifer", 
                            "Conifer_plantation",
                            df.out$DomGroup)
  
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Eucalyptus" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75 &
                              df.out$Group2 == "Invasive", 
                            "Eucalyptus",
                            df.out$DomGroup)
  
  ### Conditions for Group1 == Invasive
  ### "Broadleaved"        "Conifer_plantation" "Conifer"            "Eucalyptus"         "Invasive"   
  df.out$DomGroup <- ifelse(!is.na(df.out$Group1) &
                              !is.na(df.out$Group2) &
                              df.out$Group1 == "Invasive" & 
                              df.out$BAg1_p >= 25 &
                              df.out$BAg1_p < 75, 
                            df.out$Group2,
                            df.out$DomGroup)
  
  return(df.out)
  
}

#BASONET_dasom_2001[BASONET_dasom_2001$Species == 45,]$Name <- "Quercus ilex"
unique_plots <- sort(unique(BASONET_dasom_2001_summary_plot$Basonet_Plot_Date))

for (i in 1:length(unique_plots)) {
  print(i)
  if(i == 1) {
    ### Select the rows for this plot
    plot.i <- BASONET_dasom_2001[BASONET_dasom_2001$Basonet_Plot_Date == unique_plots[i],]
    plot.i.out <- VegGroup_BA_rank(df = plot.i, 
                              Species = "Species",
                              Area_plot_ha = Area_plot_ha)
    ### Attach Plot ID
    plot.i.out$Basonet_Plot_Date <- unique_plots[i]
  } else if (i > 1) {
    ### Select the rows for this plot
    plot.j <- BASONET_dasom_2001[BASONET_dasom_2001$Basonet_Plot_Date == unique_plots[i],]
    plot.j.out <- VegGroup_BA_rank(df = plot.j, 
                              Species = "Species",
                              Area_plot_ha = Area_plot_ha)
    ### Attach Plot ID
    plot.j.out$Basonet_Plot_Date <- unique_plots[i]
    ### Bind to previous plot
    plot.i.out <- rbind(plot.i.out,plot.j.out)
  }
}

### Now JOIN with the Plot summary
BASONET_dasom_2001_summary_plot <- merge(BASONET_dasom_2001_summary_plot,plot.i.out, by = "Basonet_Plot_Date")


# 2.1 Load 2021 data ---------------------------------------

BASONET_dasom_2021 <- read_excel("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Stand_data_Basonet/BasonetDasomN.xlsx", 
                                 sheet = "pies202122")
View(BASONET_dasom_2021)

unique(BASONET_dasom_2021$Idparcela2) ### This is a mess
sort(unique(BASONET_dasom_2021$Idparcela))

BASONET_spp_2021 <- read_excel("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/Stand_data_Basonet/BasonetDasomN.xlsx", 
                               sheet = "parcelas2021222")
#View(BASONET_spp_2021)

### Compare columns BASONET_dasom_2021
sort(unique(BASONET_dasom_2021$Idparcela))
sort(unique(BASONET_spp_2021$Idparcela))

setdiff(sort(unique(BASONET_dasom_2021$Idparcela)),
        sort(unique(BASONET_spp_2021$Idparcela)))

length(setdiff(sort(unique(BASONET_dasom_2021$Idparcela)),
               sort(unique(BASONET_spp_2021$Idparcela)))) 
### Only 1 plot
setdiff(sort(unique(BASONET_spp_2021$Idparcela)),
        sort(unique(BASONET_dasom_2021$Idparcela))) ### 46 plots. Very strange naming...

BASONET_spp_2021$Basonet_Plot_Date <- paste0("BASONET_", BASONET_spp_2021$th, "_", BASONET_spp_2021$parc, "_2021")


### Left join for the 
BASONET_spp_2021[BASONET_spp_2021$Idparcela=="01_5074767",]

BASONET_dasom_2021[BASONET_dasom_2021$Idparcela=="48_5074767",]

BASONET_spp_2021j <- BASONET_spp_2021[,c("th","parc","Idparcela","Basonet_Plot_Date")]
BASONET_spp_2021j$Date <- 2021

### Add Plot ID information
BASONET_dasom_2021 <- left_join(BASONET_dasom_2021, BASONET_spp_2021j)
unique(BASONET_dasom_2021$Basonet_Plot_Date)


# 2.2 Calculate Basal area 2021 by plot ---------------------------------------

#BASONET_dasom_2021$Diam_ave <- NA
summary(BASONET_dasom_2021$Diam1_mm)
summary(BASONET_dasom_2021$Diam2_mm)
BASONET_dasom_2021$Diam_ave_cm <- (BASONET_dasom_2021$Diam1_mm + BASONET_dasom_2021$Diam2_mm)/20
hist(BASONET_dasom_2021$Diam_ave_cm, breaks=30)
BASONET_dasom_2021$BArea_m2 <- pi*((BASONET_dasom_2021$Diam_ave_cm/200)^2)
Area_plot <- pi*(25^2)
Area_plot_ha <- Area_plot/10000

### What are the quantiles?
quantile(BASONET_dasom_2021$Diam_ave_cm, p=c(0.25,0.75,seq(from=0.1, to=0.9, by=0.1)))
#    10%    20%    30%    40%    50%    60%    70%    80%    90% 
# 12.400 15.000 17.600 21.450 25.100 29.050 34.305 42.950 50.685
summary(BASONET_dasom_2021$Diam_ave_cm)

### Attach species
BASONET_dasom_2021 <- merge(BASONET_dasom_2021,spp_names, by.x="Species",by.y="Spx")

length(unique(BASONET_dasom_2021$Basonet_Plot_Date))

### Summarise total BA (m2/ha) and QMD (cm) by plot
BASONET_dasom_2021_summary_plot <- BASONET_dasom_2021 %>%
  group_by(., Basonet_Plot_Date) %>%
  summarise(BA_Total_m2ha = sum(BArea_m2, na.rm = TRUE)/Area_plot_ha,
            QMD_cm = sqrt(sum(Diam_ave_cm^2, na.rm=TRUE)/n()))

View(BASONET_dasom_2021_summary_plot)
hist(BASONET_dasom_2021_summary_plot$BA_Total_m2ha, breaks=30)
hist(BASONET_dasom_2021_summary_plot$QMD_cm, breaks=30)
plot(BASONET_dasom_2021_summary_plot$QMD_cm, BASONET_dasom_2021_summary_plot$BA_Total_m2ha)

### summarise

#BASONET_dasom_2021[BASONET_dasom_2021$Species == 45,]$Name <- "Quercus ilex"
unique_plots <- sort(unique(BASONET_dasom_2021_summary_plot$Basonet_Plot_Date))

rm(plot.i,plot.i.out,plot.j,plot.j.out)

for (i in 1:length(unique_plots)) {
  print(i)
  if(i == 1) {
    ### Select the rows for this plot
    plot.i <- BASONET_dasom_2021[BASONET_dasom_2021$Basonet_Plot_Date == unique_plots[i],]
    plot.i.out <- Spp_BA_rank(df = plot.i, 
                              Species = "Species",
                              Diam_ave_cm = "Diam_ave_cm",
                              BArea_m2 = "BArea_m2", 
                              Area_plot_ha = Area_plot_ha)
    ### Attach Plot ID
    plot.i.out$Basonet_Plot_Date <- unique_plots[i]
  } else if (i > 1) {
    ### Select the rows for this plot
    plot.j <- BASONET_dasom_2021[BASONET_dasom_2021$Basonet_Plot_Date == unique_plots[i],]
    plot.j.out <- Spp_BA_rank(df = plot.j, 
                              Species = "Species",
                              Diam_ave_cm = "Diam_ave_cm",
                              BArea_m2 = "BArea_m2", 
                              Area_plot_ha = Area_plot_ha)
    ### Attach Plot ID
    plot.j.out$Basonet_Plot_Date <- unique_plots[i]
    ### Bind to previous plot
    plot.i.out <- rbind(plot.i.out,plot.j.out)
  }
}

### Now JOIN with the Plot summary
BASONET_dasom_2021_summary_plot <- merge(BASONET_dasom_2021_summary_plot, plot.i.out, by = "Basonet_Plot_Date")


# 2.3 Assign vegetation group to 2021 data ---------------------------------------------

BASONET_dasom_2021$Vegetation <- ifelse(BASONET_dasom_2021$Species %in% invasive, "Invasive",
                                        ifelse(BASONET_dasom_2021$Species %in% conifers_plantations, "Conifer_plantation",
                                               ifelse(BASONET_dasom_2021$Species %in% eucalyptus, "Eucalyptus",
                                                      ifelse(BASONET_dasom_2021$Species %in% conifers_native, "Conifer",
                                                             ifelse(BASONET_dasom_2021$Species %in% broadleaved, "Broadleaved", 
                                                                    BASONET_dasom_2021$Species)))))

#BASONET_dasom_2021[BASONET_dasom_2021$Species == 45,]$Name <- "Quercus ilex"
unique_plots <- sort(unique(BASONET_dasom_2021_summary_plot$Basonet_Plot_Date))

rm(plot.i,plot.i.out,plot.j,plot.j.out)

for (i in 1:length(unique_plots)) {
  print(i)
  if(i == 1) {
    ### Select the rows for this plot
    plot.i <- BASONET_dasom_2021[BASONET_dasom_2021$Basonet_Plot_Date == unique_plots[i],]
    plot.i.out <- VegGroup_BA_rank(df = plot.i, 
                                   Species = "Species",
                                   Area_plot_ha = Area_plot_ha)
    ### Attach Plot ID
    plot.i.out$Basonet_Plot_Date <- unique_plots[i]
  } else if (i > 1) {
    ### Select the rows for this plot
    plot.j <- BASONET_dasom_2021[BASONET_dasom_2021$Basonet_Plot_Date == unique_plots[i],]
    plot.j.out <- VegGroup_BA_rank(df = plot.j, 
                                   Species = "Species",
                                   Area_plot_ha = Area_plot_ha)
    ### Attach Plot ID
    plot.j.out$Basonet_Plot_Date <- unique_plots[i]
    ### Bind to previous plot
    plot.i.out <- rbind(plot.i.out,plot.j.out)
  }
}

### Now JOIN with the Plot summary
BASONET_dasom_2021_summary_plot <- merge(BASONET_dasom_2021_summary_plot,plot.i.out, by = "Basonet_Plot_Date")



#3. ### Join 2001 and 2021 data, and soil data ---------------------------------------------

table(BASONET_dasom_2021_summary_plot$DomGroup)
table(BASONET_dasom_2001_summary_plot$DomGroup)

### check colnames
colnames(BASONET_dasom_2021_summary_plot)
colnames(BASONET_dasom_2001_summary_plot)

### Bring Year
BASONET_dasom_2021_summary_plot$Date <- 2021
BASONET_dasom_2001_summary_plot$Date <- 2001

### Join soil and stand data
BASONET_dasom <- rbind(BASONET_dasom_2021_summary_plot,BASONET_dasom_2001_summary_plot)
BASONET.PdGn <- left_join(BASONET.PdGn, BASONET_dasom)

### Drop plots with invasive vegetation from Broadleaved, mixed, conifer
exclude <- (BASONET.PdGn$DomGroup %in% c("Broadleaved","Mixed","Conifer") &
              BASONET.PdGn$Group1 == "Invasive") | 
  (!is.na(BASONET.PdGn$Group2) &
    BASONET.PdGn$DomGroup %in% c("Broadleaved","Mixed","Conifer") &
  BASONET.PdGn$Group2 == "Invasive")
BASONET.PdGn <- BASONET.PdGn[!exclude, ]
#BASONET.PdGn[787,c("DomGroup","Group2","Group1")]

BASONET.PdGn$DomGroup <- factor(BASONET.PdGn$DomGroup,
                                levels=c("Broadleaved","Mixed","Conifer", 
                                         "Mixed_plantation","Conifer_plantation",
                                         "Eucalyptus","Broadleaved_Eucalyptus","Invasive"))

### Split by depth
BASONET.PdGn.0_20 <- BASONET.PdGn[BASONET.PdGn$Depth_interval == "0-20 cm", ]
BASONET.PdGn.20_40 <- BASONET.PdGn[BASONET.PdGn$Depth_interval == "20-40 cm", ]
BASONET.PdGn.0_20$DomGroup <- factor(BASONET.PdGn.0_20$DomGroup)
BASONET.PdGn.20_40$DomGroup <- factor(BASONET.PdGn.20_40$DomGroup)



# ### 4. Prepare indicators by soil district -------------------------------------

ggplot(data=BASONET.PdGn[!is.na(BASONET.PdGn$TOC) & 
                           !is.na(BASONET.PdGn$Clay) &
                           !is.na(BASONET.PdGn$PdGn) &
                           !is.na(BASONET.PdGn$DomGroup)&
                           BASONET.PdGn$Date ==2021,]) +  
  geom_boxplot(aes(x = DomGroup, y = TOC, fill = DomGroup)) +
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
                           !is.na(BASONET.PdGn$DomGroup)&
                           BASONET.PdGn$Date ==2021,]) +  
  geom_boxplot(aes(x = as.factor(PdGn),
                   y = TOC_Clay, 
                   fill = DomGroup)) +
  ylab(label="TOC / Clay") +
  xlab(label="Vegetation") +
  theme_bw() +
  # scale_fill_manual(values= map.PdGn.Comb26.k9$branch.centroids.ord$colors,
  #                   name="Soil district") +
  geom_hline(yintercept =1/13,
             color = "red",
             lwd = 1,  
             linetype = "dashed")+
  ylim(0,0.4)


ggplot(data=BASONET.PdGn[!is.na(BASONET.PdGn$TOC) & 
                           !is.na(BASONET.PdGn$Clay) &
                           !is.na(BASONET.PdGn$PdGn) &
                           !is.na(BASONET.PdGn$DomGroup)&
                           BASONET.PdGn$Date ==2021,]) +  
  geom_boxplot(aes(x = DomGroup,
                   y = TOC_Clay, 
                   color = DomGroup),
               na.rm = FALSE,
               ## Note:
               position = position_dodge(preserve = 'single')) +
  ylab(label="TOC / Clay") +
  xlab(label="soil district") +
  theme_bw() +
  scale_color_viridis(discrete=TRUE)+
  # scale_fill_manual(values= map.PdGn.Comb26.k9$branch.centroids.ord$colors,
  #                   name="Soil district") +
  geom_hline(yintercept =1/13,
             color = "red",
             lwd = 1,  
             linetype = "dashed")+
  ylim(0,0.5)

ggplot(data=BASONET.PdGn[!is.na(BASONET.PdGn$TOC) & 
                           !is.na(BASONET.PdGn$PdGn) &
                           !is.na(BASONET.PdGn$DomGroup)&
                           BASONET.PdGn$Date ==2021,]) +  
  geom_point(aes(x = QMD_cm ,
                 y = TOC, 
                 color = as.factor(PdGn),
                 shape = Depth_interval)) +
  ylab(label="TOC") +
  #xlab(label="QMD (cm)") +
  theme_bw() +
  #scale_color_viridis(discrete=TRUE)+
  scale_color_manual(values= map.PdGn.Comb26.k9$branch.centroids.ord$colors,
                    name="Soil district") +
  facet_wrap(~ DomGroup)



### pH
ggplot(data=BASONET.PdGn.0_20[!is.na(BASONET.PdGn.0_20$pH) & 
                           !is.na(BASONET.PdGn.0_20$PdGn) &
                           !is.na(BASONET.PdGn.0_20$DomGroup),]) +  
  geom_boxplot(aes(y= pH, x=DomGroup, fill= DomGroup)) +
  ylab(label="pH") +
  xlab(label="Vegetation") +
  theme_bw() +
  # scale_fill_manual(values= map.PdGn.Comb26.k9$branch.centroids.ord$colors,
  #                   name="Soil district")+ 
  facet_wrap(~ as.factor(PdGn))







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

## Join with Pedogenon and Soil data

BASONET.PdGn.subset <- BASONET.PdGn[complete.cases(BASONET.PdGn[,c("PdGn","DomGroup","Group1","Depth_interval",
                                                                   "Silt","Clay","pH","TOC","TN_ppm","K_ppm",
                                                                   "Ca_ppm", "Mg_ppm","CECef")]) &
                                      BASONET.PdGn$Date == 2021, ]

Response <- BASONET.PdGn.subset[,c("Silt","Clay","pH","TOC",
                                   "TN_ppm","K_ppm","Ca_ppm", "Mg_ppm","CECef")]

#Response <- scale(Response,center=TRUE, scale=TRUE)

### PERMANOVA
library(vegan)
set.seed(1984)
permanova.BASONET <- adonis2(Response ~ PdGn*Group1+Depth_interval,
                             data=BASONET.PdGn.subset,
                             method = "mahalanobis",
                             permutations=999)
permanova.BASONET
densityplot(permustats(permanova.BASONET))

plot(varpart(Y=Response, ~ Group1 , ~ PdGn ,
             data = BASONET.PdGn.subset))

BASONET.dist <- vegdist(x = Response, method = "mahalanobis")

set.seed(1984)
BASONET_disper <- betadisper(d =BASONET.dist,
                             group=BASONET.PdGn.subset$Group1)
plot(BASONET_disper)

save.image("C:/Users/mercedes.roman/Desktop/SELVANS/WP1/R_output/6.SoilCondition/6.SoilConditionVegReclass24042024.RData")

###############################################################################################