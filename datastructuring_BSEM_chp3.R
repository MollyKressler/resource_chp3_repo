## Structuring data for path analysis in the package NIMBLE 
## Data from BRUVS collected by Sarah Driscoll and Henriette Grimmel. Habitat data provided by the BBFSF and Matthew Smukall. Shark detection data provided by Evan Byrnes, Clemency White, the BBFSF, and Matthew Smukall, with support from Vital Heim.

## created by Molly Kressler, 12 April 2023


## Load workspace

pacman::p_load(tidyverse,sf,ggplot2,cowplot,patchwork,flextable,sf,ggsn)
setwd('/Users/mollykressler/Documents/data_phd')

## Load data

fish<-st_as_sf(st_read('resource_chp3/predictions_from_simplifiedBRTs_Gerreidae_and_SWSpecies_feb23.shp'),crs='WGS84')
	# shannon index for species, modelled by a BRT, simplified (low/medium/high densities of seagrass, sargassum, urban and rocky, deepwater, dist2shore): sppsmPD
	# Gerreidae counts, modelled by a BRT, simplified (low density seagrass and dist2shore ONLY): btGerrSimPD

sharksANDhabitat<-st_as_sf(st_read('cleaned_detections/detections20192020_NOghosts_withUPDATEDhabitat_sept22.shp'),crs='WGS84')# starting with the data from chapter 2 of MK thesis, Habitat or safety?. Sharks have PCLs from 60-100cm. The data has been filtered for ghosts detections at a threshold of 6hrs (see chp2 or the manuscript github for details). Habitat has been assigned by the effectve detection range buffer area which is centred around the receiver point location. Habtat data has been calculated from shapefiles provided by Matthew Smukall and the BBFSF. Buffer shape is missign from the shapefile, so you need to load that separately, in order to do the spatial maths on fish. 

buffs<-st_as_sf(st_read('buffers_2019_2020.shp'),crs='WGS84')%>%mutate(buff_ar=st_area(geometry))

	# check:
	ggplot()+geom_sf(data=sharksANDhabitat)+geom_sf(data=buffs,alpha=0.5,fill='violetred2',col=NA)+theme_bw()
	# SHOULD SEE: 35 black dots (receivers) and pale pink crcles around each black dot (the effective detection range of 211m)

### Modelling dataframe 
# three components: sharkiness, fishiness, habitat
	# point detections of sharks (receivers)
	# fishiness variable estmates modelled by boosted regression tree analysis and assigned by detection buffer area (average of hexagons in buffer)
	# habitat data assigned by detection point buffer area - only the variables that are in the BRT?

	## Deprecated approach
		# Step 1: cut fishiness for the buffer area 
		ggplot()+geom_sf(data=fishINbuffs)

			fishINbuffs<-st_intersection(buffs,fish) # cuts the hexagons but then there are multiple rows per buff. this old problem. Also: for Gerr, add the hexagon values because the predictions are counts, for SW-Spp average because those are not additive
			# struggling to do it as an sf - convert to a df then match to the sf - shark with habitat via buffer IDs

			a<-as_data_frame(fishINbuffs)%>%select(-geometry)
			
			b2<-a%>%
				group_by(buffID)%>%
				mutate(gerr_sum=mean(btGerSimPD),sppSW_avg=mean(sppsmPD))%>%
				filter(row_number()==1)%>%
				select(buffID,gerr_sum,sppSW_avg)
			# buffers 12 & 13 are completely outside of the known habitat (and therefore outside of the predictions for the BRT)
			# buffers &, partially outside/inside. These need reconstructed Gerreidae counts (sppSW is fine because it's an average)
				partials<-c('r10','r17')
				partial<-a%>%
					filter(buffID%in%partials)%>%
					group_by(buffID)%>%
					mutate(gerr_sum=(mean(btGerSimPD)),sppSW_avg=mean(sppsmPD))%>%
					filter(row_number()==1)%>%
					select(buffID,gerr_sum,sppSW_avg)
				partial$gerr_sum<-as.numeric(partial$gerr_sum)

			# join partials to mains, and add missing with 'NA's
				missing<-tibble(buffID=c('r12','r13'),gerr_sum=c(NA,NA),sppSW_avg=c(NA,NA))
				summary_fish_in_buffs<-bind_rows(b2%>%filter(!buffID%in%partials),partial,missing)

			# make summary_fish_in_buffs an sf again
				sf_fishinessBYbuffs<-left_join(buffs,summary_fish_in_buffs,by='buffID')


		# Step 2: Match fishiness to the sharksANDhabitat based on the buffer ID 
			sf_fishinessBYbuffs
			sharksANDhabitat

			data<-st_join(sharksANDhabitat,sf_fishinessBYbuffs)
			stopifnot(nrow(data)==nrow(sharksANDhabitat)) # check 
			summary(data$gerr_sum)
			summary(fish$btGerSimPD) # Gerries only predicted up to 197 per hexagon. maybe I should take the average. Take the average. Adjust above. 
			# average makes the buffer estimates more appropriate.

	#### DIFFERENT APPROACH, 2/5/2023 : predict the BRTs straight into a buffer+habitat df. 
	# this would be mathematically better, than taking lots of averages. 

		simple_gerr<-readRDS('resource_chp3/model_RDS/simplified_gerreidae_BRT_feb23_lowdensitySG_dist2shore.RDS')
		simple_spp<-readRDS('resource_chp3/model_RDS/simplified_speciesShannonIndex_BRT_feb23_baresandSeasonRemoved.RDS')

		simple.models<-c('simple_gerr','simple_spp')

		land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')

		buffswithhab<-st_join(buffs,sharksANDhabitat)

		sf4preds<-buffswithhab%>%group_by(buffID)%>%filter(row_number()==1)%>%dplyr::select('buffID','prop_ldsg','prop_medsg','prop_hdsg','prop_sarg','prop_urb_r','prop_deep','prop_marsh','prop_man','prop_spong','prop_inlvg','prop_unkwn','dist_cmg')%>%mutate(dist2shore=st_distance(st_centroid(geometry),st_union(land)))
			sf4preds$dist2shore<-as.numeric(sf4preds$dist2shore)
		df4preds<-as.data.frame(sf4preds)%>%select(-geometry)
		df4preds$dist2shore<-as.numeric(df4preds$dist2shore)
		### you left off here, you need to predict the brts into the buffer areas
		pacman::p_load(dismo,gbm)

		sf4<-left_join(df4preds,sf4preds,by='buffID',keep=FALSE)%>%
			dplyr::select(buffID,geometry)

		for(i in simple.models){
			mod<-get(i)
			predicts<-predict.gbm(mod,df4preds,n.trees=mod$gbm.call$best.trees,type='response')
			predicts<-as.data.frame(predicts)
			predicts.df<-bind_cols(df4preds,predicts)
			sf4b<-sf4%>%mutate(preds=predicts.df$predicts)
			new=paste0(i,'PD')
			sf4<-sf4b%>%dplyr::rename_with(.fn=~paste0(new),.cols='preds')
			}

			# back transform the gerreidae (poisson, log)
			sf4<-sf4%>%mutate(bt.gerr.simpPD=exp(.$simple_gerrPD),.before='geometry')%>%rename(sppsimpPD=simple_sppPD)

			# save the updated sf4preds and df4preds
				st_write(sf4,'sf_with_predictions_INTO_BUFFERS_fromBRTs_feb23.shp',append=FALSE) # append command set to FALSE replaces the old file with the new one. If set to TRUE it would add the new data as further layers.
				st_write(sf4,'sf_with_predictions_INTO_BUFFERS_fromBRTs_feb23.csv',delete_layer=TRUE,delete_dsn=TRUE)










fishiness_in_buffers<-st_as_sf(st_read('sf_with_predictions_INTO_BUFFERS_fromBRTs_feb23.shp'),crs='WGS84')%>%dplyr::select(-smpl_PD)

# spatially join sharks & habitat to fishiness in buffers 

data<-st_join(sharksANDhabitat,fishiness_in_buffers)
# yay!

	st_write(data,'resource_chp3/data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_may23.shp',driver='ESRI Shapefile')
	st_write(data,'resource_chp3/data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_may23.csv',driver='CSV')


