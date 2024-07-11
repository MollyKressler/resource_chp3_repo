## BRUVS data cleaner and tidier. 

## created by Molly M Kressler
## updated numerous times but major update in July 2024 

pacman::p_load(sf, tidyverse, dplyr, lubridate, ggplot2, cowplot, patchwork, readr, readxl)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd/')

############ July 2024 update: use multiple fish Families from 2014 and 2018 ############
### match habitat data to bruvs data by year 


######
## - Import data and convert to sf
######
d14 <- read_csv('BullockDedmanETAL2024AnimBehav_familyabundance_BRUVS_data.csv')%>%
	filter(!is.na(BRUV_ID)) # 2014 data from Bullock et al. Anim. Behav. 2024, available for download at https://data.mendeley.com/datasets/2vj2tw8jkw/1. 
s14 <- st_as_sf(d14, coords = c('GPS_W', 'GPS_N'), crs = 'WGS84')%>%
	mutate(source = 'Bullock', year = 2014, .before = geometry)



d18 <- read_excel('MF_calculatedFamilysabundanceperBRUV_JUNE2022_Driscoll_BRUVSdata.xlsx', sheet = 'Sheet1', range = 'A1:F103')
l18 <- read_excel('MF_calculatedFamilysabundanceperBRUV_JUNE2022_Driscoll_BRUVSdata.xlsx', sheet = 'BRUVS_final data', range = 'A1:C103') %>% rename(BRUV_ID = SiteName)
dd18 <- left_join(d18, l18, by = 'BRUV_ID')
s18 <- st_as_sf(dd18, coords = c('Longitude','Latitude'), crs = 'WGS84')%>%
	mutate(source = 'Driscoll', year = 2018, .before = geometry)


######
## - Match to habitat data 
######

h14 <- st_as_sf(st_read('winter2014habitat_hexagon_grid_NOland_ECdata.shp'), crs='WGS84')
h18 <- st_as_sf(st_read('summer2018habitat_hexagon_grid_NOland_ECdata.shp'), crs='WGS84')


	joined14<-st_join(s14, h14,join=st_nearest_feature)
	joined14

	joined18<-st_join(s18,h18,join=st_nearest_feature)
	joined18

	# save 
	write_sf(joined14,'bruvs2014_data_joinedWITHhabitat_winter2014hab_july2024.shp',driver='ESRI Shapefile')
	write_sf(joined18,'bruvs2018_data_joinedWITHhabitat_summer18hab_july2024.shp',driver='ESRI Shapefile')

	write_sf(joined14,'bruvs2014_data_joinedWITHhabitat_winter2014hab_july2024.csv',driver='CSV')
	write_sf(joined18,'bruvs2018_data_joinedWITHhabitat_summer18hab_july2024.csv',driver='CSV')

######
## - Join 2014 & 2018 data 
######

	a14 <- joined14 %>%
		dplyr::select(-Carangidae, -Hemiramphidae)%>%
		mutate(BRUV_ID = paste0(source,BRUV_ID))%>%
		rename(Gerridae=Gerreidae)

	a18 <- joined18 %>%
		dplyr::select(-Sphraenidae)%>%
		mutate(BRUV_ID = paste0(source,BRUV_ID))

	data <- bind_rows(a14,a18)
	stopifnot(nrow(data) == nrow(a14) + nrow(a18) & ncol(data) == ncol(a14) & ncol(data) == ncol(a18)) # check 


######
## - Save
######


	write_sf(data,'bruvs2014AND2018_data_joinedWITHhabitat_summer18habANDwinter2014hab_july2024.shp',driver='ESRI Shapefile')

	write_sf(data,'bruvs2014AND2018_data_joinedWITHhabitat_summer18habANDwinter2014hab_july2024.shp',driver='CSV')


























