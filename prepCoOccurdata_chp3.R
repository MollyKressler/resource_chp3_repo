## Assigning spatial data to co-occurence dataset (of large sharks and juveniles) for chapter 3 BSEM modeling

pacman::p_load(tidyverse,sf,lubridate,ggplot2,flextable)

setwd('/Users/mollykressler/Documents/data_phd')

#############################
## - Load datasets
#############################

## co-occurence data 

	c <- read.csv('lemonspredators_20192020blacktipsANDbulls/cooccurence_juvenilesANDlargesharks_within60minoflargedett.csv')%>%
		select(-X)%>%
		as_tibble()%>%
		mutate(rID=paste0('r',rID)) # habitat by buffers df has buffID/rID with 'r' before the number 

## habitat data by buffers

	b <- st_as_sf(st_read('buffers20192020_withEmilyCourmierhabitat_data4Winter2020.shp'),crs='WGS84')%>%
		mutate(rID=buffID)

#############################
## - Assign spatial geometry to co-occurence data
#############################

	data <- left_join(b,c)
	## NAs from 'b' are from where the receiver lies outside the known habitat extent. These are classified as prop_unknown=1, and are fine.
	summary(data)

	data %>% filter(is.na(lowdensg)) # 1 receiver without predator data (rID = 100 from construction of pressure metric)

	data <- data %>%
		filter(!is.na(pressure)) # removes the 1 row with NAs for pressure. See Note below for justification. 

#############################
## - Save data
#############################

	st_write(data,'coocurrencemetrics_withHabitat4winter2020_fromEC_nov23.csv',driver='CSV')

	st_write(data,'coocurrencemetrics_withHabitat4winter2020_fromEC_nov23.shp',driver='ESRI Shapefile')

#############################
## - Notes for further use of dataset 
#############################


## There will be fewer receivers in this set than in the juvenile detections set because they were detected at fewer receivers. That's fine. We are assigning pressure by nearest receiver so it's either going to be the same receiver or the closest. 