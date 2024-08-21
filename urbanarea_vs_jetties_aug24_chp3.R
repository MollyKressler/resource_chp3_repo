## Considering urban areas versus locations of jetties

## created by Molly Kressler, 20 August 2024

# response to comments from co-authors about the appropriateness of 'distance form the nearest jetty' as a metric of the human disturbance in Bimini


pacman::p_load(tidyverse,sf,ggplot2,cowplot,patchwork,flextable,sf,ggsn,dismo,gbm)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd')

j<-st_zm(st_read('boatlaunchesBimini.kml'),crs='WGS84')%>%
mutate(jetty=seq(1,15,1),.before='geometry')%>%
dplyr::select(-Description) # jetties visible from Google Earth in Bimini. NOT including small single residence jetties. 

h20 <- st_as_sf(st_read('winter2020habitat_polygonised_forSF_ECdata.shp'), crs = 'WGS84')



## Calculate the total urban area on land

h20 %>% mutate(area = st_area(.), .before = 'geometry') # 5,141,267 m^2 (approx. 5mill m^2)

## Calculate the distance of urban area that is adjacent to aquatic habitat

wet <- st_union(h20[1:3,])
urb <- st_as_sf(h20[5,])

i <- st_as_sf(st_intersection(wet, urb))

urb.plot <- ggplot()+
	geom_sf(data = i)+
	labs(subtitle = 'Where Bare/Urban intersects Aquatic')+
	theme_bw() ### urban/bare feature is not a better proxy for human disturbance because it rgisters bare pristine areas as urban, if we were to assign is as distance from urban area. 

ggsave(urb.plot, file = 'coauthors_feedback_response_UrbanVSJetty.png', device = 'png', unit = 'in', dpi = 850, height = 5, width = 5)

## Calculate the total urban area on land within 250m of each jetty location in the summer 2024 analysis of chapter 3 

jb <- j%>%
		st_buffer(500)%>%
		st_intersection(., urb)

intx.area<- st_area(jb)
sum(intx.area) # 1,322,185 m^2

intersecton_with_jetties<- ggplot()+
	geom_sf(data = urb, col = 'grey55', fill = 'grey85')+
	geom_sf(data = jb, fill = 'deeppink4', col = NA, alpha = 0.5)+
	geom_sf(data = j, pch = 19, size = 1, col = 'goldenrod2', fill = 'goldenrod2')+
	theme_void()+
	labs(justfication = 'center', caption = '500m radius around Jetties (pink) \n\ Overlap (500m radii) = 3.56mill m^2 \n\ Total Urban/Bare area = 5.14mill m^2')
ggsave(intersecton_with_jetties, file = 'coauthors_feedback_response_Urban_intersect_Jetty.png', device = 'png', unit = 'in', dpi = 850, height = 4, width = 3.5)

















