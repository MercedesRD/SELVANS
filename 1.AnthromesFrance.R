#### Pedogenon maps for France

### Date 08/11/2023
### Author: Mercedes Roman Dobarco
### Objective: Plot the evolution of anthromes in France for identifying our reference point

# ### 6. Anthromes ------------------------------------------------------------

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

### mys study area
### I bring shapefile
france_regions <- read_sf("C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Administrative/regions-20180101-shp/regions-20180101.shp")

### Only metropolitan France
france_regions <- france_regions[france_regions$code_insee %in% c(11:93),]
plot(france_regions["nom"])
france <- st_union(france_regions)

# Import the precomputed Anthromes-12k-DGG dataset. These include fixed inputs like land area, region, and biome as well as the anthrome classifications for HYDE 3.2 baseline, upper, and lower scenarios.

# HYDE/Anthromes fixed inputs and baseline scenario
an12_dgg_inputs <- read_sf('C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Anthromes/an12_dgg_inputs/Anthromes-12k-DGG/an12_dgg_inputs.shp')
an12_dgg_baseline <- read_csv('C:/Users/mercedes.roman/Desktop/SELVANS/France/Covariates/Anthromes/an12_dgg_baseline.csv')

## Merge
an12_dgg <- an12_dgg_inputs %>%
  left_join(an12_dgg_baseline, by = 'id')

### remove unneeded files to save RAM
rm(an12_dgg_baseline,  an12_dgg_inputs, france_regions)

### subset for France
an12_dgg_Fr <- st_intersection(an12_dgg, france)

### I don´t want to rasterise just yet. Only to prepare the diagram of the evolution of anthromes in France
colnames(an12_dgg_Fr)
rm(an12_dgg); gc()

### what are the unique values present for anthromes?
names(an12_dgg_Fr)
unique_anthromes <- apply(X = as.data.frame(an12_dgg_Fr)[,6:80], MARGIN = 2, FUN=function(x) unique(x))
sort(unique(unlist(unique_anthromes)))
#  11 12    22 23 24 31 32 33 34 41 42 43 51 52 53 54
#  11 12 21 22 23 24 31 32 33 34 41 42 43 51 52 53 54
#  11 12    22 23 24 31 32 33 34 41 42 43 51 52 53 54 61 62
# 11 Urban
# 12 Mixed settlements
# 21 Rice villages
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
# 61 Wild forests
# 62 Wild treeless or barren

an12_dgg_Fr_df <- as.data.frame(an12_dgg_Fr)

## long format
my_cols <- colnames(an12_dgg_Fr_df[,6:80])
an12_dgg_Fr_df_long <- an12_dgg_Fr_df[,6:80] %>% 
  pivot_longer(., names_to = "year", cols= my_cols, values_to="Anthromes")

# Transform into factor
an12_dgg_Fr_df_long$Anthromes <- factor(an12_dgg_Fr_df_long$Anthromes, 
                                        levels=c("11","12", 
                                                 "22", "23", "24",
                                                 "31", "32", "33", "34",
                                                 "41", "42", "43",
                                                 "51", "52", "53", "54",
                                                 "61", "62"))

### Summarise, count number of cells by anthrome and year
anthromes_summary <- an12_dgg_Fr_df_long %>%
  group_by(., year, Anthromes) %>%
  summarise(Area = n())
anthromes_summary <- anthromes_summary %>% 
  mutate(time_step = ordered(year, levels = my_cols)) ### New variable, ordered levels of time step
area_total <- sum(anthromes_summary[anthromes_summary$year=="0AD",]$Area)
anthromes_summary$Area_perc <- round(anthromes_summary$Area/area_total*100, digits=1)


length(unique(anthromes_summary$Anthromes)) ### 18 colors I need
library(colorspace)

q18 <-c("11"="brown","12"="brown1",
        "22"="cornflowerblue","23"="darkorchid","24"="#F161AE", 
        "31"="#72F6D2","32"="#CBC937","33"="#F2F081","34"="#F0EFA9",
        "41"="#E09F2A","42"="#ECB776","43"="#F6DCB5",
        "51"="#19AA34","52"="#5FCE2C","53"="#B6E284","54"="#F6F8E3",
        "61"= "#b7edec", "62"="#f7f6f0")

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
# 61 Wild forests "#b7edec"
# 62 Wild treeless or barren "#f7f6f0"

anthromes_plot <- anthromes_summary %>% 
  ggplot(aes(time_step, Area_perc)) +
  geom_col(aes(fill =  Anthromes, color =  Anthromes), width = .8, size = .1)+
  geom_vline(xintercept = c('0AD', '800AD', "1860AD"), color = 'black') +
  geom_hline(yintercept = .5, linetype = 2, color = 'black') +
  scale_fill_manual(values = q18,
                    labels= c("Urban"," Mixed settlements", 
                              "Irrigated villages","Rainfed villages","Pastoral villages", 
                              "Residential irrigated croplands","Residential rainfed croplands",
                              "Populated croplands", "Remote croplands",
                              "Residential rangelands","Populated rangelands","Remote rangelands",
                              "Residential woodlands","Populated woodlands","Remote woodlands",
                              "Inhabited drylands",
                              "Wild forests", "Wild treeless" ))+
  scale_color_manual(values = q18,
                     labels = c("Urban"," Mixed settlements",
                                "Irrigated villages","Rainfed villages","Pastoral villages",
                                "Residential irrigated croplands","Residential rainfed croplands",
                                "Populated croplands", "Remote croplands",
                                "Residential rangelands","Populated rangelands","Remote rangelands",
                                "Residential woodlands","Populated woodlands","Remote woodlands",
                                "Inhabited drylands",
                                "Wild forests", "Wild treeless" )) +
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


########################################