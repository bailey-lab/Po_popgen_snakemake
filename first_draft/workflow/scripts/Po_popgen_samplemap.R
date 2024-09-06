## Generate maps representing sample origins and other metadata

setwd("~/Documents/P. Ovale Genomic Analysis")

library(tidyverse)
library(ggplot2)
library(rnaturalearth) 
#install.packages("ggpattern")
library(rnaturalearthdata)
library(sf)
library(ggspatial)
library(colorspace)
library(dplyr)
library(ggpattern)
#library(XML)
#library(tmaptools)
#library(magick)

#First, you need to download shape files. Ideally, you would be able to do this entirely within the R ecosystem. However, as you will see below, I ran into issues with the DRC regions.

admin10 <- ne_download(scale="large", type = "admin_1_states_provinces",
                       category = "cultural", returnclass = "sf")
rivers10 <- ne_download(scale = 10, type = 'rivers_lake_centerlines', 
                        category = 'physical', returnclass = "sf")
lakes10 <- ne_download(scale = "medium", type = 'lakes', 
                       category = 'physical', returnclass = "sf")
oceans10 <- ne_download(scale = "medium", type = "coastline",
                        category = 'physical', returnclass = "sf")
sov110 <- ne_download(scale="medium", type = "sovereignty",
                      category = "cultural", returnclass = "sf")
Africa <- ne_countries(scale="medium", type = "sovereignty", continent = "Africa", returnclass = "sf", )

#https://gadm.org/download_country.html then download the shapefile. You will get a zip file that you need to unzip. For our purposes, you want the _1.shp files which are at the state/region/province level.

DRC_regions <- st_read("shape_files/gadm41_COD_shp/gadm41_COD_1.shp")

Ethiopia_regions <- st_read("shape_files/gadm41_ETH_shp/gadm41_ETH_1.shp")

Tanzania_regions <- st_read("shape_files/gadm41_TZA_shp/gadm41_TZA_1.shp")

Cameroon_regions <- st_read("shape_files/gadm41_CMR_shp/gadm41_CMR_1.shp")

Senegal_regions <- st_read("shape_files/gadm41_SEN_shp/gadm41_SEN_1.shp")

IvoryCoast_regions <- st_read("shape_files/gadm41_CIV_shp/gadm41_CIV_1.shp")

Nigeria_regions <- st_read("shape_files/gadm41_NGA_shp/gadm41_NGA_1.shp")

all_regions <- rbind(DRC_regions, Ethiopia_regions, Tanzania_regions,Cameroon_regions,Senegal_regions,IvoryCoast_regions,Nigeria_regions)

focus_countries <- dplyr::filter(Africa, adm0_a3 == "COD" | adm0_a3 == "ETH" | adm0_a3 == "TZA" | adm0_a3 == "CMR" | adm0_a3 == "NGA" | adm0_a3 == "SEN" | adm0_a3 == "CIV")
po_countries <- dplyr::filter(Africa, adm0_a3 == "COD" | adm0_a3 == "ETH" | adm0_a3 == "TZA" | adm0_a3 == "CMR" | adm0_a3 == "SEN" | adm0_a3 == "CIV")
poc_countries <- dplyr::filter(Africa, adm0_a3 == "COD" | adm0_a3 == "ETH" | adm0_a3 == "TZA" | adm0_a3 == "CMR")
  
pow_countries <-  dplyr::filter(Africa, adm0_a3 == "COD" | adm0_a3 == "ETH" | adm0_a3 == "TZA" | adm0_a3 == "CMR" | adm0_a3 == "SEN" | adm0_a3 == "CIV")
pf_countries <- dplyr::filter(Africa, adm0_a3 == "COD" | adm0_a3 == "ETH" | adm0_a3 == "TZA" | adm0_a3 == "NGA" | adm0_a3 == "SEN" | adm0_a3 == "CIV")
po_regions <- dplyr::filter(all_regions, COUNTRY == "Democratic Republic of the Congo" & NAME_1 == "Kinshasa" |
                              COUNTRY == "Cameroon" & NAME_1 == "Nord-Ouest" |  
                              COUNTRY == "Ethiopia" & NAME_1 == "Amhara" |
                              COUNTRY == "Tanzania" & NAME_1 == "Pwani" |
                              COUNTRY == "Tanzania" & NAME_1 == "Songwe" )

poc_regions <- dplyr::filter(all_regions, COUNTRY == "Democratic Republic of the Congo" & NAME_1 == "Bas-Uele")
pow_regions <- dplyr::filter(all_regions, COUNTRY == "Democratic Republic of the Congo" & NAME_1 == "Sud-Kivu" |
                COUNTRY == "Senegal" & NAME_1 == "Tambacounda" |  
                GID_0 == "CIV" & NAME_1 == "Bas-Sassandra")

pf_regions <- dplyr::filter(all_regions, COUNTRY == "Democratic Republic of the Congo" & NAME_1 == "Kinshasa" |
                              COUNTRY == "Ethiopia" & NAME_1 == "Amhara" |
                              COUNTRY == "Tanzania" & NAME_1 == "Tanga" |
                              COUNTRY == "Tanzania" & NAME_1 == "Morogoro" |
                              COUNTRY == "Tanzania" & NAME_1 == "Kagera" |
                               COUNTRY == "Senegal" & NAME_1 == "Dakar" | 
                              GID_0 == "CIV" & NAME_1 == "Abidjan" |
                              COUNTRY == "Nigeria" & NAME_1 == "Lagos" |
                              ###Added for ease of visuals
                              COUNTRY == "Senegal" & NAME_1 == "Thiès" | 

                              GID_0 == "CIV" & NAME_1 == "Lagunes" |

                              COUNTRY == "Nigeria" & NAME_1 == "Ogun")

water <- st_read("shape_files/Africa_waterbody/Africa_waterbody.shp")


###Plot of all study locations (all Po)
ggplot() + geom_sf(data=Africa, fill="gray90") + 
  geom_sf(data = po_countries, fill="gray96", lwd=0.5) + 
  geom_sf(data = po_regions, fill="#E7298A") +
  geom_sf(data = poc_regions, fill="#E7298A") +
  geom_sf(data = pow_regions, fill="#E7298A") +
  geom_sf(data = water, fill="aliceblue") +
  coord_sf(ylim = c(-20,22), xlim= c(-18,47), expand = FALSE) +
  ggtitle("Study Sites")+
  theme_light() +
  theme(panel.background = element_rect(fill = 'aliceblue', colour = 'grey97')) +
  theme(axis.ticks = element_blank(),
        
        axis.text.x = element_blank(),
        
       axis.text.y = element_blank()) +
  
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  xlab("") +
  
  ylab("") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.title = element_text(hjust = 0.5))

library("ggpubr")
ggsave("maps/all_po_map.png", plot = last_plot(), device = "png", dpi = 600, width = 10, height = 5)


###Plot of specific species origins
ggplot() + geom_sf(data=Africa, fill="gray90") + 
  geom_sf(data = po_countries, fill="gray96", lwd=0.5) + 
  geom_sf(data = po_regions, fill="#E7298A") +
  geom_sf(data = poc_regions, fill="#7570B3") +
  geom_sf(data = pow_regions, fill="#66A61E") +
  geom_sf(data = water, fill="aliceblue") +
  coord_sf(ylim = c(-20,22), xlim= c(-18,47), expand = FALSE) +
  ggtitle(expression(paste(italic("P. ovale"), " WGS Sample Sites")))+
  theme_light() +
  theme(panel.background = element_rect(fill = 'aliceblue', colour = 'grey97')) +
  theme(axis.ticks = element_blank(),
        
        axis.text.x = element_blank(),
        
        axis.text.y = element_blank()) +
  
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  xlab("") +
  
  ylab("") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("maps/poc-pow_map.png", plot = last_plot(), device = "png", dpi = 600, width = 10, height = 5)


###Plot of P. falciparum origins
ggplot() + geom_sf(data=Africa, fill="gray90") + 
  geom_sf(data = pf_countries, fill="gray96", lwd=0.5) + 
  geom_sf(data = pf_regions, fill="red2") +
  geom_sf(data = water, fill="aliceblue") +
  coord_sf(ylim = c(-20,22),  xlim= c(-18,47),expand = FALSE) +
  ggtitle(expression(paste(italic("P. falciparum"), " Sample Sites")))+
  theme_light() +
  theme(panel.background = element_rect(fill = 'aliceblue', colour = 'grey97')) +
  theme(axis.ticks = element_blank(),
        
        axis.text.x = element_blank(),
        
        axis.text.y = element_blank()) +
  
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  xlab("") +
  
  ylab("") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("maps/pf_map.png", plot = last_plot(), device = "png", dpi = 600, width = 10, height = 5)


###specific regions for unique coloring
kinshasa <- dplyr::filter(all_regions, COUNTRY == "Democratic Republic of the Congo" & NAME_1 == "Kinshasa")
kinshasa_coord <- sf::st_point_on_surface(kinshasa)
dschang <- dplyr::filter(all_regions, COUNTRY == "Cameroon" & NAME_1 == "Nord-Ouest")
amhara  <- dplyr::filter(all_regions, COUNTRY == "Ethiopia" & NAME_1 == "Amhara")
bagamoyo <- dplyr::filter(all_regions, COUNTRY == "Tanzania" & NAME_1 == "Pwani")
songwe <- dplyr::filter(all_regions, COUNTRY == "Tanzania" & NAME_1 == "Songwe")
DRC_east_poc <- dplyr::filter(all_regions, COUNTRY == "Democratic Republic of the Congo" & NAME_1 == "Bas-Uele")

DRC_east_pow <- dplyr::filter(all_regions, COUNTRY == "Democratic Republic of the Congo" & NAME_1 == "Sud-Kivu")
senegal <- dplyr::filter(all_regions, COUNTRY == "Senegal" & NAME_1 == "Tambacounda")
ivory <- dplyr::filter(all_regions, GID_0 == "CIV" & NAME_1 == "Bas-Sassandra")

#values = c("DRC (West)" = "#D95F02", "DRC (East)"= "#7570B3","Tanzania"="#E7298A", "Ethiopia"="#66A61E", "Cameroon"="#E6AB02", "Senegal"="#A6761D","Ivory Coast"="#666666"


###Plot of poc samples, colored like PCA
ggplot() + geom_sf(data=Africa, fill="gray90") + 
  geom_sf(data = poc_countries, fill="gray96", lwd=0.5) + 
  geom_sf(data = kinshasa, fill="#D95F02") +
  geom_sf(data = dschang, fill="#E6AB02") +
  geom_sf(data = amhara, fill="#66A61E") +
  geom_sf(data = bagamoyo, fill="#E7298A") +
  geom_sf(data = DRC_east_poc, fill="#7570B3") +
  geom_sf(data = water, fill="aliceblue") +
  coord_sf(ylim = c(-20,22), xlim= c(-18,47), expand = FALSE) +
  ggtitle(expression(paste(italic("P. ovalecurtisi"), " PCA Sample Sites")))+
  theme_light() +
  theme(panel.background = element_rect(fill = 'aliceblue', colour = 'grey97')) +
  theme(axis.ticks = element_blank(),
        
        axis.text.x = element_blank(),
        
        axis.text.y = element_blank()) +
  
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  xlab("") +
  
  ylab("") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("maps/poc_map.png", plot = last_plot(), device = "png", dpi = 600, width = 7, height = 5)


###Plot of pow samples, colored like PCA
ggplot() + geom_sf(data=Africa, fill="gray90") + 
  geom_sf(data = po_countries, fill="gray96", lwd=0.5) + 
  geom_sf(data = kinshasa, fill="#D95F02") +
  geom_sf(data = dschang, fill="#E6AB02") +
  geom_sf(data = amhara, fill="#66A61E") +
  geom_sf(data = bagamoyo, fill="#E7298A") +
  geom_sf(data = songwe, fill="skyblue3") +
  geom_sf(data = DRC_east_pow, fill="#7570B3") +
  geom_sf(data = senegal, fill="#A6761D") +
  geom_sf(data = ivory, fill="#666666") +
  geom_sf(data = water, fill="aliceblue") +
  coord_sf(ylim = c(-20,22), xlim= c(-18,47), expand = FALSE) +
  ggtitle(expression(paste(italic("P. ovalewallikeri"), " PCA Sample Sites")))+
  theme_light() +
  theme(panel.background = element_rect(fill = 'aliceblue', colour = 'grey97')) +
  theme(axis.ticks = element_blank(),
        
        axis.text.x = element_blank(),
        
        axis.text.y = element_blank()) +
  
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  xlab("") +
  
  ylab("") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("maps/pow_map.png", plot = last_plot(), device = "png", dpi = 600, width = 7, height = 5)



###Plot of poc and pow samples, colored like PCA
### Using points over study sites rather than coloring in entire regions
ggplot() + geom_sf(data=Africa, fill="gray90") + 
  geom_sf(data = po_countries, fill="gray96", lwd=0.5) + 
  geom_sf(data = water, fill="aliceblue") +
  #Kinshasa
  geom_point(aes(x=15.3, y = -4.317), color="#D95F02", size = 4) +
  #Dschang
  geom_point(aes(x = 10.04, y = 5.27), 
             color="#E6AB02", size=4) +
  #Amhara
  geom_point(aes(x = 37.95, y = 11.66), 
             color="#66A61E", size=4) +
  #Bagamoyo
  geom_point(aes(x = 38.9, y = -6.433), 
             color="#E7298A", size=4) +
  #Songwe
  geom_point(aes(x = 32.533, y = -8.517), 
             color="skyblue3", size=4) +
  #Bas-Uele
  geom_point(aes(x = 24.733, y = 2.8), 
             color="#7570B3", size=4) +  
  #Sud-Kivu
  geom_point(aes(x = 28.867, y = -2.5), 
              color="#7570B3", size=4) +  
  #Tambacounda
  geom_point(aes(x = -13.667, y = 13.767), 
             color="#A6761D", size=4) + 
  #Soubre, Ivory Coast
  geom_point(aes(x = -6.6, y = 5.783), 
             color="#666666", size=4) + 
  coord_sf(ylim = c(-25,31), xlim= c(-18,48), expand = FALSE) +
  ggtitle(expression(paste(italic("P. ovale"), " PCA Sample Sites")))+
  theme_light() +
  theme(panel.background = element_rect(fill = 'aliceblue', colour = 'grey97')) +
  theme(axis.ticks = element_blank(),
        
        axis.text.x = element_blank(),
        
        axis.text.y = element_blank()) +
  
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  xlab("") +
  
  ylab("") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("maps/po_point_map.png", plot = last_plot(), device = "png", dpi = 600, width = 7, height = 5)







####Older version of Poc and Pow sample map, coloring in all regions rather than placing points over cities
ggplot() + geom_sf(data=Africa, fill="gray90") + 
  geom_sf(data = po_countries, fill="gray96", lwd=0.5) + 
  geom_sf(data = kinshasa, fill="#D95F02") +
  geom_sf(data = dschang, fill="#E6AB02") +
  geom_sf(data = amhara, fill="#66A61E") +
  geom_sf(data = bagamoyo, fill="#E7298A") +
  geom_sf(data = songwe, fill="skyblue3") +
  geom_sf(data = DRC_east_pow, fill="#7570B3") +
  geom_sf(data = DRC_east_poc, fill="#7570B3") +
  geom_sf(data = senegal, fill="#A6761D") +
  geom_sf(data = ivory, fill="#666666") +
  geom_sf(data = water, fill="aliceblue") +
  coord_sf(ylim = c(-25,31), xlim= c(-18,48), expand = FALSE) +
  ggtitle(expression(paste(italic("P. ovale"), " PCA Sample Sites")))+
  theme_light() +
  theme(panel.background = element_rect(fill = 'aliceblue', colour = 'grey97')) +
  theme(axis.ticks = element_blank(),
        
        axis.text.x = element_blank(),
        
        axis.text.y = element_blank()) +
  
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  xlab("") +
  
  ylab("") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("maps/po_map.png", plot = last_plot(), device = "png", dpi = 600, width = 7, height = 5)


### Zach's code beneath

sequencing_samples <- data.table::fread("Pm_Po_Sample_Locations_ALL.csv", header = TRUE, data.table = FALSE)

country_only <- sequencing_samples |> subset(Country == "Cameroon" | Country == "Nigeria")

by <- join_by("admin" == "Country")

country_only <- right_join(focus_countries, country_only, by, multiple = "all")

all_regions$NAME_1 <- all_regions$NAME_1 |> stringr::str_replace("Bas-Uele", "Bas-Uélé") |> stringr::str_replace("Gambela Peoples", "Gambela")

by <- join_by("NAME_1" == "Region")

regional_data <- sequencing_samples |> subset(Country != "Cameroon" & Country != "Nigeria")

regional_data <- right_join(all_regions, regional_data, by, multiple = "all")

water<-st_read("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/MSMT non-falciparum Projects/Africa_waterbody.shp")

focus_countries_pm <- focus_countries |> subset(adm0_a3 != "COD" & adm0_a3 != "ETH")

pm_countries <- country_only |> subset(Species == "Pm")

po_countries <- country_only |> subset(Species == "Po")

pm_regions <- regional_data |> subset(Species == "Pm")

po_regions <- regional_data |> subset(Species == "Po")

pmplot <- ggplot() + 
  geom_sf(data=Africa, fill="gray90")+
  geom_sf(data = focus_countries_pm, fill="gray96", lwd=0.5) +
  geom_sf(data = Tanzania_regions, fill="gray96") +
  geom_sf(data = pm_countries,aes(geometry=geometry, fill=Samples)) +
  geom_sf(data = pm_regions, aes(geometry = geometry, fill = Samples)) +
    scale_fill_distiller(type="seq", palette = "PuRd", direction = 1, name="Pm Samples", breaks = c(2, 4, 6, 8, 10, 12))+
  geom_sf(data = water, fill="aliceblue") +
  #annotation_scale(location = "bl", width_hint = 0.25) +
  #annotation_north_arrow(location = "bl", which_north = "true", 
  #                       pad_x = unit(0.85, "in"), pad_y = unit(0.2, "in"),
  #                       style = north_arrow_fancy_orienteering)+
  #annotate("text", x = 30.5, y = -0.7, label = "Uganda", color="grey30", size=4 , fontface="italic") +
  #annotate("text", x = 30., y = -3.2, label = "Burundi", color="grey30", size=3 , fontface="italic") +
  #annotate("text", x = 29.9, y = -1.9, label = "Rwanda", color="grey30", size=3 , fontface="italic") +
  #annotate("text", x = 38, y = -1.5, label = "Kenya", color="grey30", size=6 , fontface="italic") +
  #annotate("text", x = 38, y = -12, label = "Mozambique", color="grey30", size=5 , fontface="italic") +
  #annotate("text", x = 30.5, y = -10, label = "Zambia", color="grey30", size=6 , fontface="italic") +
  #annotate("text", x = 33.8, y = -11.5, label = "Malawi", color="grey30", size=3 , fontface="italic") +
  #annotate("text", x = 29.6, y = -7.5, label = "DRC", color="grey30", size=3 , fontface="italic") +
  coord_sf(ylim = c(-35,38), expand = FALSE) +
  #coord_sf(xlim = c(29.5, 40), ylim = c(-1,-12), expand = TRUE) +
  ggtitle(expression(paste(italic("P. malariae "), "WGS Samples")))+
  theme_light() +
  theme(panel.background = element_rect(fill = 'aliceblue', colour = 'grey97')) +
  theme(axis.ticks = element_blank(),
        
        axis.text.x = element_blank(),
        
        axis.text.y = element_blank()) +
  
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  xlab("") +
  
  ylab("") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.title = element_text(hjust = 0.5))

plot_ratio <- get_asp_ratio(Africa)

pmplotfile <- paste0("Pm_Samples.png")

ggsave(pmplotfile, plot = pmplot, dpi=600, width = plot_ratio * 5, height = 5)

magick_pmplot <- image_trim(image_read(pmplotfile)) #this crops whitespace which hasn't been a major issue for me before but is for some reason a big one here

image_write(magick_pmplot, pmplotfile)

ggplot() + geom_sf(data=Africa, fill="gray90") + geom_sf(data = focus_countries, fill="gray96", lwd=0.5)

poplot <- ggplot() + 
  geom_sf(data=Africa, fill="gray90")+
  geom_sf(data = focus_countries, fill="gray96", lwd=0.5) +
  geom_sf(data = all_regions, fill="gray96") +
  geom_sf(data = po_countries,aes(geometry=geometry, fill=Samples)) +
  geom_sf(data = po_regions, aes(geometry = geometry, fill = Samples)) +
  scale_fill_distiller(type="seq", palette = "Greens", direction = 1, name="Po Samples", breaks = c(2, 4, 6, 8))+
  geom_sf(data = water, fill="aliceblue") +
  #annotation_scale(location = "bl", width_hint = 0.25) +
  #annotation_north_arrow(location = "bl", which_north = "true", 
  #                       pad_x = unit(0.85, "in"), pad_y = unit(0.2, "in"),
  #                       style = north_arrow_fancy_orienteering)+
  #annotate("text", x = 30.5, y = -0.7, label = "Uganda", color="grey30", size=4 , fontface="italic") +
  #annotate("text", x = 30., y = -3.2, label = "Burundi", color="grey30", size=3 , fontface="italic") +
  #annotate("text", x = 29.9, y = -1.9, label = "Rwanda", color="grey30", size=3 , fontface="italic") +
  #annotate("text", x = 38, y = -1.5, label = "Kenya", color="grey30", size=6 , fontface="italic") +
  #annotate("text", x = 38, y = -12, label = "Mozambique", color="grey30", size=5 , fontface="italic") +
  #annotate("text", x = 30.5, y = -10, label = "Zambia", color="grey30", size=6 , fontface="italic") +
  #annotate("text", x = 33.8, y = -11.5, label = "Malawi", color="grey30", size=3 , fontface="italic") +
#annotate("text", x = 29.6, y = -7.5, label = "DRC", color="grey30", size=3 , fontface="italic") +
coord_sf(ylim = c(-35,38), expand = FALSE) +
  #coord_sf(xlim = c(29.5, 40), ylim = c(-1,-12), expand = TRUE) +
  ggtitle(expression(paste(italic("P. ovale "), "spp. WGS Samples")))+
  theme_light() +
  theme(panel.background = element_rect(fill = 'aliceblue', colour = 'grey97')) +
  theme(axis.ticks = element_blank(),
        
        axis.text.x = element_blank(),
        
        axis.text.y = element_blank()) +
  
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  xlab("") +
  
  ylab("") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.title = element_text(hjust = 0.5))

poplotfile <- paste0("Po_Samples.png")

ggsave(poplotfile, plot = poplot, dpi=600, width = plot_ratio * 5, height = 5)

magick_poplot <- image_trim(image_read(poplotfile)) #this crops whitespace which hasn't been a major issue for me before but is for some reason a big one here

image_write(magick_poplot, poplotfile)
