## Structuring data for path analysis in the package NIMBLE 
## Data from BRUVS collected by Sarah Driscoll and Henriette Grimmel. Habitat data provided by the BBFSF and Matthew Smukall. Shark detection data provided by Evan Byrnes, Clemency White, the BBFSF, and Matthew Smukall, with support from Vital Heim.

## created by Molly Kressler, 12 April 2023

########################################################
## Notes 
########################################################

## 24 May 2023: structuring second dataframe for path analysis. DF based on the hexagon predictions for shark use and fish occurence across the study site. Vars: pr(use) for sharks, BRT estimates of Gerreidae count (back transformed), BRT estmates of species diversity, habitat (raw). Calculated Var: cross metric for Gerries and Diversity (gerreidae * diversity) 'fishy'. 

## 11 July 2023: new parameters after collaboration meeting with Adrian Gleiss and Evan Byrnes, also RBS, on 7.7.23. New params: deltaTemp,depthReceiver,dist2urban, prPred. 
	## Comprehensive notes in Obsidian canvas for chapter 3, but to summarise: 
		## need data from EB for deltaTemp
		## depth, deltaTemp & prPred will be available only for the pointdata (dateframe 1 here)
		## waiting for data from Matt for prPred, but in the meantime I'm going to simulate a probability from a gamma distribution
		## dist2urbn can be calculated with habitat shapefile and sf
		## depth data is availabel in the receiver metadata 

## 25 July 2023: receievd predator species detection data from MS at BBFSF. cleaned in separate file. Here, integrated into df1 (point data) as a probability of detection of a predator. (1) calculate receiver specific detection probabilities; (2) join them to df1/pointdata. 

## 10 August 2023
	## Update data frames with new habitat data from Emily Courmier - processed in 'loading_cleaning_habitatdata_EmilyCourmier.R' & formatted into hexagons, buffers, and detections in 'joininghab2detections_ghostsremoved_EmilyChabitatdata.R'
	## add depth of NEAREST receiver to hexdata (fishiness)

## 10 October 2023: new metric for predation pressure
	## relPropPD doesn't make (the best) biological sense - a receiver with high Ns of large shark detections doesn't have a 'real' and high predation risk if no or few juvenile sharks also use that area. 
	## develop a cooccurence metric of juveniles with large sharks within an hour at each receiver
		## need filtered detections of juveniles
		## detections of large sharks
	## apply to pointdata, and join by nearest feature for hexdata

## 16 November 2023: new metric for predation pressure, 'pressure' = the co-occurence of juveniles and large sharks adjusted for the number of days of co-occurences relative to length of observation period
	## apply to pointdata, and join by nearest feature for hexdata
	## apply to hexdata, and join by nearest feature for hexdata


########################################################
########################################################

## Overview of dataframes
	## data frame 1: spatially organised by receivers. point detections of sharks, BRT estimates of Gerriedae count and Species diversity, and habitat. This is for the process model of sharkiness in path analysis. 
	## data frame 2: spatially organised by hexagons. Predictions of shark use (Pr(use)), Gerreidae back-transformed counts from BRT, Species diversity index from BRT, cross metric of 'fishy', and habitat. This is data for the process model of fishiness in path analysis. 

## Load workspace

pacman::p_load(tidyverse,sf,ggplot2,cowplot,patchwork,flextable,sf,ggsn,dismo,gbm)
setwd('/Users/mollykressler/Documents/data_phd')

## Load workspace, R Server 
pacman::p_load(MCMCvis,tidyverse,sf,nimble,devtools,flextable,webshot2,sfdep,sp,spdep,beepr, HDInterval)
setwd("~/resource/data_and_RDS_NOTforupload")



## Load data

	hab.new <-st_as_sf(st_read('winter2020habitat_hexagon_grid_NOland_ECdata.shp'),crs='WGS84')

	fish<-st_as_sf(st_read('resource_chp3/predictions_from_simplifiedBRTs_Gerreidae_and_SWSpecies_aug23.shp'),crs='WGS84')
		# shannon index for species, modelled by a BRT, simplified (low/medium/high densities of seagrass, sargassum, urban and rocky, deepwater, dist2shore): sppsmPD
		# Gerreidae counts, modelled by a BRT, simplified (low density seagrass and dist2shore ONLY): btGerrSimPD

	sharksANDhabitat<-st_as_sf(st_read('cleaned_detections/detections20192020_NOghosts_withEmilyCourmierhabitat_data4Winter2020_sept22.shp'),crs='WGS84')# starting with the data from chapter 2 of MK thesis, Habitat or safety?. Sharks have PCLs from 60-100cm. The data has been filtered for ghosts detections at a threshold of 6hrs (see chp2 or the manuscript github for details). Habitat has been assigned by the effectve detection range buffer area which is centred around the receiver point location. Habtat data has been calculated from shapefiles provided by Matthew Smukall and the BBFSF. Buffer shape is missign from the shapefile, so you need to load that separately, in order to do the spatial maths on fish. 

	buffs<-st_as_sf(st_read('buffers20192020_withEmilyCourmierhabitat_data4Winter2020.shp'),crs='WGS84')%>%
		mutate(buff_ar=st_area(geometry),buffID=as.factor(buffID))%>%
		distinct()%>%
		group_by(buffID)%>%slice(1) # for reasons I don't know why, r18 has two rows of the same value 
		summary(buffs)


		# check:
		ggplot()+geom_sf(data=sharksANDhabitat)+
			geom_sf(data=buffs,alpha=0.4,fill='violetred2',col=NA)+
			theme_bw()
		# SHOULD SEE: 35 black dots (receivers) and pale pink crcles around each black dot (the effective detection range of 211m)

	# land & centroid of mangroves

		land<-st_as_sf(st_union(st_read('bim_onlyland_noDots.kml')),crs='WGS84')
		cmg<-st_as_sf(st_read('habitat_model/cleaned_code/habitatorsafetylemonsharks_ms/data/centroid_mangroves_north.shp'),crs='WGS84')

	# hexagon df of Pr(use) with habitat BEFORE updating. This dataframe has information for low AND high tide, so you need to filter out for low only

		pruse<-st_as_sf(st_read('habitat_model/no ghosts of 6hr threshold/habitatmodel_preds_method5_withANDwithoutREFUGE_ghostsremoved_glmers_dec22_updatedhabitat.shp'),crs='WGS84')%>%
		filter(tidephs=='L')%>%
			dplyr::select(jcode,prp_brs,prp_lds,prp_mds,AIC_m5_,geometry)%>%
			rename(method5_PrUSE=AIC_m5_)

	# hexagon df of Gerreidae (BRT) & SW Species (BRT) - if these hexs dont align with 'pruse', re-predict the simplified BRTs into the Pr(use) hexagon df. 
		
		fishes<-st_as_sf(st_read('resource_chp3/predictions_from_simplifiedBRTs_Gerreidae_and_SWSpecies_aug23.shp'),crs='WGS84')%>%
			dplyr::select(jcode,prp_lds,prp_mds,prp_hds,dist_cmg,dist2shore,sppsimpPD,btGerSmPD,geometry)%>%
			mutate(fishy=(sppsimpPD*btGerSmPD),.before=geometry)
		fishes

	# depth of receiver metadata
		d <- st_as_sf(st_read('biminireceivers_withLocations_20190414to20201213.shp'),crs='WGS84')
		d$id <-as.factor(d$id)
		summary(d)
		d2 <- d%>%distinct()%>%st_as_sf()

### Modelling dataframe 1
# three components: sharkiness, fishiness, habitat

	# point detections of sharks (receivers)
	# fishiness variable estmates modelled by boosted regression tree analysis and assigned by detection buffer area (average of hexagons in buffer)
	# habitat data assigned by detection point buffer area - only the variables that are in the BRT?

	## Deprecated approach
		# Step 1: cut fishiness for the buffer area 
		# Step 2: Match fishiness to the sharksANDhabitat based on the buffer ID 

  #### DIFFERENT APPROACH, 2/5/2023 : predict the BRTs straight into a buffer+habitat df. 
	# this would be mathematically better, than taking lots of averages. 

		simple_gerr<-readRDS('resource_chp3/model_RDS/simplified_GerriesShannonIndex_BRT_aug23_highDenSgANDseasonremoved.RDS')
		simple_spp<-readRDS('resource_chp3/model_RDS/simplified_speciesShannonIndex_BRT_aug23_highDenSgANDseasonremoved.RDS')

		simple.models<-c('simple_gerr','simple_spp')

	### There's an issue with the buffers (r33 specifically) where it wasn't identifying the habitat in the buffer. Here's code to fix it in the document. 
		buffswithhab<-st_join(buffs,sharksANDhabitat)
		r33<-buffs%>%filter(buffID=='r33')
		ggplot()+geom_sf(data=r33withhab2)+geom_sf(data=land)

		sf_use_s2(FALSE) # turn speherical geometry off

		buffs2.join<-st_join(st_make_valid(r33),st_make_valid(hab.new))%>%mutate(hab_geo=st_geometry(st_intersection(st_make_valid(st_geometry(hab.new)),st_make_valid(st_geometry(r33)))))%>%mutate(hab_area=st_area(hab_geo))

		buffs2.spread<-buffs2.join%>%spread(key= 'habitat',value = 'hab_area',fill=0)%>%group_by(buffID)%>%summarise(across(c('buff_ar','highdensg','lowdensg','meddensg','sargassum','rocks.urb',),max))%>%mutate(prop_ldsg=lowdensg/buff_ar,prop_medsg=meddensg/buff_ar,prop_hdsg=highdensg/buff_ar,prop_sarg=sargassum/buff_ar,prop_urb_r=rocks.urb/buff_ar)%>%dplyr::select(buffID,buff_ar,prop_ldsg,prop_medsg,prop_hdsg,prop_sarg,prop_urb_r,geometry)%>%add_column(prop_unkwn=0,.before='geometry')%>%
			add_column(prop_marsh=0,prop_deep=0,prop_inlvg=0,prop_brs=0,prop_man=0,prop_spong=0)
			
		buffs.new<-buffs2.spread
			{
			buffs.new$prop_brs<-as.numeric(buffs.new$prop_brs)
			buffs.new$prop_ldsg<-as.numeric(buffs.new$prop_ldsg)
			buffs.new$prop_medsg<-as.numeric(buffs.new$prop_medsg)
			buffs.new$prop_hdsg<-as.numeric(buffs.new$prop_hdsg)
			buffs.new$prop_man<-as.numeric(buffs.new$prop_man)
			buffs.new$prop_sarg<-as.numeric(buffs.new$prop_sarg)
			buffs.new$prop_marsh<-as.numeric(buffs.new$prop_marsh)
			buffs.new$prop_urb_r<-as.numeric(buffs.new$prop_urb_r)
			buffs.new$prop_deep<-as.numeric(buffs.new$prop_deep)
			buffs.new$prop_inlvg<-as.numeric(buffs.new$prop_inlvg)
			buffs.new$prop_unkwn<-as.numeric(buffs.new$prop_unkwn)
			}
		
		cmg<-st_as_sf(st_read('habitat_model/centroid_mangroves_north.shp'))

		buffs.new2<-buffs.new%>%add_column(dist_cmg=st_distance(st_centroid(buffs.new$geometry),cmg),.before='geometry')

		r33withhab2<-buffs.new2%>%
			dplyr::select('buffID','prop_brs','prop_ldsg','prop_medsg','prop_hdsg','prop_sarg','prop_urb_r','prop_deep','prop_marsh','prop_man','prop_spong','prop_inlvg','prop_unkwn','dist_cmg')%>%
			mutate(dist2shore=st_distance(st_centroid(geometry),st_union(land)))%>%
			filter(row_number()==1)
		
			r33withhab2$dist2shore<-as.numeric(r33withhab2$dist2shore)
			r33withhab2$dist_cmg<-as.numeric(r33withhab2$dist_cmg)
		
			r33df<-as.data.frame(r33withhab2)%>%dplyr::select(-geometry)


	sf4preds<-buffs%>%group_by(buffID)%>%
		dplyr::select('buffID','prop_ldsg','prop_mdsg','prop_hdsg','prop_unkwn','dist_cmg')%>%
		mutate(dist2shore=as.numeric(st_distance(st_centroid(geometry),st_union(land))),.before='geometry')

	df4preds<-as.data.frame(sf4preds)%>%dplyr::select(-geometry)%>%rename(prp_lds=prop_ldsg,prp_mds=prop_mdsg,prp_hds=prop_hdsg)
	sapply(df4preds,class)

	sf4<-sf4preds%>%rename(prp_lds=prop_ldsg,prp_mds=prop_mdsg,prp_hds=prop_hdsg)

	for(i in simple.models){
		mod<-get(i)
		predicts<-predict.gbm(mod,df4preds,n.trees=mod$gbm.call$best.trees,type='response')
		predicts<-as.data.frame(predicts)
		predicts.df<-bind_cols(df4preds,predicts)
		sf4b<-left_join(sf4,predicts.df)
		new=paste0(i,'PD')
		sf4<-sf4b%>%dplyr::rename_with(.fn=~paste0(new),.cols='predicts')
		}

		# back transform the gerreidae (poisson, log)
		sf4<-sf4%>%mutate(bt.gerr.simPD=exp(simple_gerrPD),.before='geometry')%>%rename(sppsimPD=simple_sppPD)
		## sf4 has lost the geometry as a result from the left join. so join again to preserve the geometry
		sf44 <- left_join(sf4preds,sf4)
		summary(sf44)
		summary(sf4preds)
		summary(sf4) # checks 

		# save the updated sf4preds and df4preds
			st_write(sf44,'resource_chp3/sf_with_predictions_INTO_BUFFERS_fromBRTs_aug23.shp',append=FALSE) # append command set to FALSE replaces the old file with the new one. If set to TRUE it would add the new data as further layers.
			write_sf(sf44,'resource_chp3/sf_with_predictions_INTO_BUFFERS_fromBRTs_aug23.csv')



	fishiness_in_buffers<-st_as_sf(st_read('resource_chp3/sf_with_predictions_INTO_BUFFERS_fromBRTs_aug23.shp'),crs='WGS84')%>%
		dplyr::select(-smpl_PD)%>%
		rename('prop_ldsg'='prp_ldsg','prop_medsg'='prp_mdsg','prop_hdsg'='prp_hdsg','prp_nkw'='prp_nkw','dist2shore'='dst2shr','dist_cmg'='dst_cmg')


	# spatially join sharks & habitat to fishiness in buffers 
	names(sharksANDhabitat)
	not.needed.vars<-c("ghost" , "FID",'prop_brs','prop_ldsg','prop_medsg','prop_hdsg','prop_unkwn','dist_cmg')
	needed.vars <- c('time','GPS_W','GPS_N','PIT','pcl','sex','tidephs','dorn','month','ghost','buffID','geometry')
	sharksANDhabitat2<-sharksANDhabitat%>%
		dplyr::select(all_of(needed.vars))

	data5<-st_join(sharksANDhabitat2,fishiness_in_buffers%>%dplyr::select(-c(prp_lds,prp_mds,prp_hds)),join=st_nearest_feature)%>%
		rename(prop_unkn=prp_nkw) ## nearest fetaure because of the buffers cut 	off by new habitat shape
	nrow(sharksANDhabitat2) ## 33701
	summary(data5)

	st_write(data5,'resource_chp3/data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_aug23.shp',driver='ESRI Shapefile',append=FALSE)
	st_write(data5,'resource_chp3/data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_aug23.csv',driver='CSV',delete_layer=TRUE,delete_dsn=TRUE)


	#### STANDARDISE & MEAN-CENTRED
		data<-read.csv('resource_chp3/data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_aug23.csv')%>%
			rename(buffID=buffID.x)

		habitats<-c('prop_ldsg','prop_medsg','prop_hdsg','prop_unkn')
		fishmetrics<-c('bt_g_PD','sppsmPD')

		fishANDhab<-data%>%
			dplyr::select(buffID,dist2shore,dist_cmg,unlist(habitats),unlist(fishmetrics),-buffID.y)%>%
			group_by(buffID)%>%
			filter(row_number()==1)

		data2<- data %>%
			dplyr::select(PIT,buffID)%>%
			group_by(PIT,buffID)%>%
			tally()
		
		## modified data augmentation based on Chapter 7.8.4 in Kery & Royle 2016 (pages 358 - 362). 

		data3<-expand.grid(PIT=unique(data2$PIT),buffID=unique(data2$buffID))%>%
			left_join(.,data2)%>%
			mutate(n=replace_na(n,0))%>%
			left_join(.,fishANDhab,by='buffID')%>%
			mutate(gerries=bt_g_PD)%>%
			mutate(fishy=sppsmPD*gerries)%>%
			dplyr::arrange(.,n)%>% 
			mutate(z=case_when(n == 0 ~ 0,n > 0 ~ 1))%>%
			mutate(standard.shark=((n-mean(n))/sd(n)))%>% 
			mutate(standard.fish=((fishy-mean(fishy))/sd(fishy)))%>%
			mutate(standard.lds=((prop_ldsg-mean(prop_ldsg))/sd(prop_ldsg)))%>% 
			mutate(standard.mds=((prop_medsg-mean(prop_medsg))/sd(prop_medsg)))%>% 
			mutate(standard.hds=((prop_hdsg-mean(prop_hdsg))/sd(prop_hdsg)))%>% 
			mutate(standard.unkn=((prop_unkn-mean(prop_unkn))/sd(prop_unkn)))%>% 
			mutate(standard.dist2shore=((dist2shore-mean(dist2shore))/sd(dist2shore)))%>% 
			mutate(standard.distcmg=((dist_cmg-mean(dist_cmg))/sd(dist_cmg)))%>%
			mutate(buffIDnum=str_remove(.$buffID,"r"))%>%
			dplyr::select(PIT,buffIDnum,standard.shark,z,standard.fish,standard.dist2shore,standard.distcmg,standard.lds,standard.mds,standard.hds)

		data3$buffIDnum<-as.numeric(data3$buffIDnum)
		sapply(data3,class)

		nrow(data3%>%filter(z==1)) # 184 'n', 376 structural zeros: N = 560		


	#### Principal component analysis for medium and high density seagrass (hex)
		# examine correlation 
		cor(data3$standard.hds,data3$standard.mds) # - 0.414, weak/moderate correlation

		# calculate the principal components
		pca1<-prcomp(data3[,9:10])

		summary(pca1)	# PC1 describes 57.65% of the cumulative variation in the hgih and medium seagrass density data. Therefore proceeding with just the first PC is better than just picking one. 
		plot(pca1,type='lines') 
		pca1$sdev^2 # eigenvalues for each component
		pca1$rotation # look at loading coefficiations, eignenvectors. Loadings illustrate the association between each PC and the original vars. high loading = more contribution
		biplot(pca1,cex=.5) # hds and mds are relatively positivelly correlated, and are equidistant from the center (i.e. neither is much better than the other). 

		# seems like PCA1 is a reasonable choice. 

		data3$sgPCA1<-predict(pca1)[,1]
		summary(data3)	

		write.csv(data3,'resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_aug23.csv')

 # 11/7/23: adding: depthReceiver, deltaTemp and prPred
	# load dataframe 1 as made above. (df1:)
	pointdata<-read.csv('resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_aug23.csv')%>%dplyr::select(-X)
	pointdata$PIT<-as_factor(pointdata$PIT)
	pointdata$buffID<-as.numeric(pointdata$buffID)
	pointdata$buffIDnum<-as.numeric(pointdata$buffIDnum)
		
	stopifnot(nrow(pointdata)==560) # check 

	## depthReceiver - to pointdata
		d # loaded above, has the receivers deployed throughout the study with depths. Need to add the buffID in order to match it to the df1
		## dont need ot re-make it every time. can import finished fle: 'biminireceivers_withLocations_20190414to20201213.shp'
			e <- sharksANDhabitat%>%dplyr::select(FID,geometry)
			e2<-e%>%group_by(FID)%>%slice_head()%>%ungroup()%>%mutate(rID=seq(1,35,1),.before=FID)
			d2 <- st_join(e2,d,left=TRUE)%>%
			group_by(rID)%>%
			slice_head()%>%
			ungroup()
		st_write(d2,'biminireceivers_withLocations_20190414to20201213.shp',driver='ESRI Shapefile') 
		d2<-st_as_sf(st_read('biminireceivers_withLocations_20190414to20201213.shp'),crs='WGS84')

		pointdata

		pointdata$buffID <- as.character(pointdata$buffID)
		d2$rID <- as.character(d2$rID)
		pointdata2 <- left_join(pointdata,d2%>%rename(buffID=rID))%>%dplyr::select(-FID)
		head(pointdata2)
		p2 <-st_as_sf(pointdata2)

		## update file to include depth and receiver information
		st_write(p2,'resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_aug23.csv',driver='CSV',delete_layer=TRUE,delete_dsn=TRUE)
		st_write(p2,'resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_aug23.shp',driver='ESRI Shapefile')


	## dist2jetty - to pointdata
		land<-st_as_sf(st_union(st_read('bim_onlyland_noDots.kml')),crs='WGS84')
		hab.new # loaded above
		df1shp<-st_as_sf(st_read('resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_aug23.shp'),crs='WGS84') # load in df1 shp as updated by 'depthReceiver'
		stopifnot(nrow(df1shp)==560) # check 

		# import Jetty locations KML file

		j<-st_zm(st_read('boatlaunchesBimini.kml'),crs='WGS84')%>%
		mutate(jetty=seq(1,15,1),.before='geometry')%>%
		dplyr::select(-Description) # jetties visible from Google Earth in Bimini. NOT including small single residence jetties. 

		ronly<-df1shp%>%group_by(bffIDnm)%>%slice_head()%>%ungroup()%>%dplyr::select(bffIDnm,geometry)

		nearest2<-tibble(r=seq(1,35,1),jetty=st_nearest_feature(ronly,j))%>%mutate(r_geo = ronly$geometry)%>%
			mutate(jetty_geo=j[match(nearest2$jetty,j$jetty),'geometry'])%>%
			mutate(jetty_id=j[match(nearest2$jetty,j$jetty),'Name'])%>%
			mutate(dist2jetty=st_distance(nearest2$r_geo,nearest2$jetty_geo$geometry,by_element=TRUE))%>%
			rename(buffID=r)%>%
			mutate(buffID=as.character(buffID),dist2jetty=as.numeric(dist2jetty))

		nearest3<-nearest2%>%dplyr::select(buffID,jetty,dist2jetty)
		nearest3


		df1shp2<-left_join(df1shp,nearest3,by='buffID')


		## update file to include jetty and dist2jetty information
		st_write(df1shp2,'resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_aug23.csv',driver='CSV',delete_layer=TRUE,delete_dsn=TRUE)
		st_write(df1shp2,'resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_aug23.shp',driver='ESRI Shapefile')

	## DEPRECATED probPred - (1) calculate receiver specific detection probabilities; (2) join them to df1/pointdata. 
		## sum(predator dettections total in year) = dT
		## sum(predator detections at receiver 'x' in year) = dX
		### Pr(predators using the area over a year relative to other receivers) = dX / dT. 
		# spatially comparable, temporally flat 
		# for each species separately, and then both together

		workingpointdata <- st_as_sf(st_read('resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_aug23.shp'),crs='WGS84')%>%
				rename(buffIDnum=bffIDnm,standard.shark=stndrd_s,standard.fish=stndrd_f,standard.dist2shore=stndr_2,standard.distcmg=stndrd_d,standard.lds=stndrd_l,standard.mds=stndrd_m,standard.hds=stndrd_h)
		## calculate the relProbPD risk per receiver
		p <- st_as_sf(st_read('lemonspredators_20192020blacktipsANDbulls/predators_detectedwithinstudyarea_MKthesis2019to2020_aug2023.shp'))%>%
		  dplyr::select(-statn_x,-statn_y)
		p <- st_as_sf(st_read('predators_detectedwithinstudyarea_MKthesis2019to2020_aug2023.shp'))%>%
		  dplyr::select(-statn_x,-statn_y)
		
			wherePredsdontgo <- expand.grid(rID=c(10,22,24),Species=c('Carcharhinus limbatus','Negaprion brevirostris'),n=0)%>%
				left_join(.,d2)%>%
				st_as_sf()%>%
				dplyr::select(rID,Species,n,geometry)
			whereLimbatusdontgo <- expand.grid(rID=c(6,9,29,30,31),Species='Carcharhinus limbatus',n=0)%>%
				left_join(.,d2)%>%
				st_as_sf()%>%
				dplyr::select(rID,Species,n,geometry)

			p2<- p %>%
				dplyr::select(Species,rID)%>%
				group_by(Species,rID)%>%
				tally()%>%
				bind_rows(wherePredsdontgo,whereLimbatusdontgo)

			p3<-expand.grid(rID=seq(1,35,1))%>%
				left_join(.,p2,multiple='all')%>%
				mutate(n=replace_na(n,0),rID=as.character(rID))%>%
				st_as_sf()
			p3 # sf with one row per species per receiver, with the detections tallied for each row.

			p4 <- p3 %>%
				group_by(rID)%>%
				summarise(ndetts=sum(n))%>%
				st_join(.,p3b%>%dplyr::select(total),join=st_nearest_feature)%>%
				mutate(relPropPD = ndetts/total,relPropPD=replace_na(relPropPD,0)) # non-species specific relative proportion of risk - not going to try and format the species specific into the table atm. 

			p4<-st_join(p3,workingpointdata,join=st_nearest_feature)%>%
				dplyr::select(-PIT,-id,-locatin,-intervl)
			p4 # receiver, with relative proportion of predator detections, and geometry
			p5<-st_cast(p4,'POINT',group_or_split=FALSE) # for 3 receivers there was a duplication in the geometry, same point geometry but organised as multipoint. This just selects the first point feature form those multipoints. 
			summary(p5)
			write_sf(p5, 'lemonspredators_20192020blacktipsANDbulls/relativepredatorrisk_at_receivers_April2019December2020_lemonsANDblacktips.shp')

			relativePropPDdettsreceivers<-ggplot()+
			geom_sf(data=land,fill='grey70',col=NA)+
			geom_sf(data=p4,aes(size=relPropPD,col=relPropPD))+
				scale_colour_gradient(low='#8f5ab0',high='#e27c0e',limits=c(0.0,1.0),name='Relative Risk')+
				scale_size(guide='none')+
				theme_bw()
			
			ggsave(relativePropPDdettsreceivers,file='lemonspredators_20192020blacktipsANDbulls/descriptive_stats_and_figures/relativepredatorrisk_at_receivers_April2019December2020_lemonsANDblacktips.png',device='png',units='in',dpi=950,height=7,width=6)

			relp <- st_as_sf(st_read('relativepredatorrisk_at_receivers_April2019December2020_lemonsANDblacktips.shp'))
      relp
		
    
    ## join predator relProbPD data to pointdata df
			workingpointdata
			p4
			w2 <- st_join(workingpointdata,p4,left=TRUE,join=st_nearest_feature)%>%
				dplyr::select(-rID,-total,-ndetts)
			w2
		
		## standardise dist2jetty, relPropPD,depth
		summary(w2)
			w3 <- w2 %>% 				
				mutate(standard.relPropPD=((as.numeric(relPropPD)-mean(as.numeric(relPropPD)))/sd(as.numeric(relPropPD))),.before='geometry')%>%
				mutate(standard.dist2jetty=((as.numeric(dist2jetty)-mean(as.numeric(dist2jetty)))/sd(as.numeric(dist2jetty))),.before='geometry')%>%
				mutate(standard.depth=((as.numeric(depth)-mean(as.numeric(depth)))/sd(as.numeric(depth))),.before='geometry')%>%
				dplyr::select(-id,intervl,lctn_cd,locatin)



		## update file to include relProbPD
		st_write(w2,'resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_aug23.csv',driver='CSV',delete_layer=TRUE,delete_dsn=TRUE)
		st_write(w2,'resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_aug23.shp',driver='ESRI Shapefile',append=FALSE)

	## deltaTemp - by receiver 
		workingpointdata <- st_as_sf(st_read('standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_nov23.shp'),crs='WGS84')%>%
				rename(buffIDnum=bffIDnm,standard.shark=stndrd_s,standard.fish=stndrd_f,standard.dist2shore=stndr_2,standard.distcmg=stndrd_d,standard.lds=stndrd_l,standard.mds=stndrd_m,standard.hds=stndrd_h)
 # DEPRECATED November 2023 - 'pressure' metric of predator and juvenile co-occurence 
	# working pointdata df
	pp<-read.csv('resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_nov23.csv')%>%
		dplyr::select(buffIDnum,standard.distcmg)
	head(pp)
	pointdataTEST<-st_as_sf(st_read('resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_nov23.shp'),crs='WGS84')%>%
		rename(buffIDnum=bffIDnm,standard.shark=stndrd_s,standard.fish=stndrd_f,standard.dist2shore=stndr_2,standard.lds=stndrd_l,standard.mds=stndrd_m,standard.hds=stndrd_h, dist2jetty=dst2jtt,standard.depth=stndrd_dp,pressure=pressur,standard.distcmg=stndrd_ds)%>%
		dplyr::select(-pressure) # updating pressure
	## don't run every time, check: weirdly got rid of some vars, so adding them back 
		pointdata2 <- left_join(pointdata,pp,by='buffIDnum',relationship='many-to-many')%>%
			rename(std.2jetty=standard.dist2jetty,std.2cmg=standard.distcmg)
		pointdata2

	## pressure metric df
	c <- st_as_sf(st_read('coocurrencemetrics_withHabitat4winter2020_fromEC_nov23.shp'),crs='WGS84')%>%
		select(pressure, geometry)
	c

	## join by rID with st_nearest_feature

	j <- st_join(df22, c, join=st_nearest_feature)
	j

	## standardise and mean centre 
	jj <- j %>% 
		mutate(standard.press=((as.numeric(pressure)-mean(as.numeric(pressure)))/sd(as.numeric(pressure))),.before=geometry) %>%
		select(-standard.dist2jetty,-standard.press) ## too many file names which match when shortened by ESRI driver. so you need to do these at the load in. 



	######### 
	## save/update file to include pressure metric
	######### 
		st_write(jj,'resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_nov23.csv',driver='CSV',delete_layer=TRUE,delete_dsn=TRUE)
		write_sf(jj,'resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_nov23.shp',driver='ESRI Shapefile',append = FALSE)

	### 01/12/2023: There are duplicated rows in the final dataset - multipel rows of a reciever at PIT number. Post hoc removing it an will re-write this section of code at a later date. 
	### FIX HERE
	working<-st_as_sf(st_read('resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_nov23.shp'),crs='WGS84')%>%
		rename(buffIDnum=bffIDnm,
			standard.shark=stndrd_s,
			standard.fish=stndrd_f,
			standard.dist2shore=stndr_2,
			standard.lds=stndrd_l,
			standard.mds=stndrd_m,
			standard.hds=stndrd_h, 
			dist2jetty=dst2jtt,
			standard.depth=stndrd_dp,
			pressure=pressur,
			standard.distcmg=stndrd_ds)%>%
			mutate(standard.press=((as.numeric(pressure)-mean(as.numeric(pressure)))/sd(as.numeric(pressure))),.before=geometry) %>%
			mutate(standard.dist2jetty=((as.numeric(dist2jetty)-mean(as.numeric(dist2jetty)))/sd(as.numeric(dist2jetty))),.after=dist2jetty)

	noduplicates <- working %>%
		distinct(PIT,rID,.keep_all=T)%>%
		dplyr::select(-standard.press,-
standard.dist2jetty) # for some reason the driver wont write the file this big, so I remove these and recaluclate them in the data upload
	noduplicates

	st_write(noduplicates,'resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_dec23.csv',driver='CSV',delete_layer=TRUE,delete_dsn=TRUE)
	write_sf(noduplicates,'resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_dec23.shp',driver='ESRI Shapefile',append = FALSE)

  ## Error in dataframe - standard.distcmg and standard.dist2shore are the same values. 4/12/23023
	# 1. extract dist2shore and distcmg for each buffer from the august version of the data. call this the repair set. 
	august <- st_as_sf(st_read('resource_chp3/data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_aug23.shp'),crs='WGS84')
	august
	repairset <- august %>% 
		dplyr::select(buffID_x,dist_cmg,dist2shore) %>%
		rename(rID = buffID_x) %>%
		group_by(rID) %>%
		distinct()
	repairset # this will take a while to generate because of the nrow and distinct function 

	# 2. join modelling dataframe 1 as updated in december 2023 to repairset, and recaluclate the standardised versions of distcmg and dist2shore
	df2 <- st_as_sf(st_read('resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_dec23.shp'),crs='WGS84')%>%
		rename(buffIDnum=bffIDnm,
			standard.shark=stndrd_s,
			standard.fish=stndrd_f,
			standard.dist2shore=stndr_2,
			standard.lds=stndrd_l,
			standard.mds=stndrd_m,
			standard.hds=stndrd_h, 
			dist2jetty=dst2jtt,
			standard.depth=stndrd_dp,
			pressure=pressur,
			standard.distcmg=stndrd_ds)
	stopifnot(nrow(pointdata)==560) # check 
	df22 <- st_join(df2, repairset%>%dplyr::select(-rID))%>%
		mutate(standard.dist2shore = ((as.numeric(dist2shore)-mean(as.numeric(dist2shore)))/sd(as.numeric(dist2shore))))%>%
		mutate(standard.distcmg = ((as.numeric(dist_cmg)-mean(as.numeric(dist_cmg)))/sd(as.numeric(dist_cmg))),.before=sgPCA1)%>%
		dplyr::select(-dist_cmg,-dist2shore)	## dataframe gets too big with everything so remove none standardised vars 
	stopifnot(nrow(df22)==560) # check, dataframe is prone to duplicatign when 
	df22

	# 3. save df22 as CSV and SHP 

		st_write(df22,'resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_dec23.csv',driver='CSV',delete_layer=TRUE,delete_dsn=TRUE)
		
		write_sf(df22,'resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_dec23.shp',driver='ESRI Shapefile',append = FALSE)



	##################################################################
	## May 2024, adding relative predation pressure to the model. bind to data here
		relp <- as.data.frame(st_as_sf(st_read('lemonspredators_20192020blacktipsANDbulls/relativepredatorrisk_at_receivers_April2019December2020_lemonsANDblacktips.shp')))%>%
		    dplyr::select(-geometry)%>%
		    mutate(rID = paste0('r',rID))
		head(relp)
		
		pointdata <- left_join(pointdata, relp, by='rID')
    head(pointdata)    

        ## calculate standardised versions 
      
    pointdata <- pointdata %>%
        mutate(standard.relr = (relPropPD-mean(relPropPD))/sd(relPropPD))

    ## save this
      st_write(pointdata, 'resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_may24.shp', driver = 'ESRI Shapefile')
      st_write(pointdata, 'resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_may24.csv', driver = 'CSV')


###############################
### Modelling dataframe 2

	# join pruse and fishes
		## new habitat data area is smaller than the previous data, and therefore the fishes area is smaller than the pruse shapefile. This produces NAs in the jcode rows of pruse that have 0 matchign fishes jcodes. remove these, no need to st_intersection/difference. 

		ggplot()+geom_sf(data=pruse,aes(fill='violetred4',alpha=0.3))+
			geom_sf(data=fishes,aes(fill='orange2',alpha=0.3))+theme_void()

		combo.sf2 <- st_join(fishes,pruse%>%dplyr::select(method5_PrUSE),join=st_nearest_feature,left=TRUE)%>%
			dplyr::rename('mtd5_PrUSE'='method5_PrUSE')%>%
			drop_na() # these are the pruse grid cells that lie outside the new habitat extent

		summary(combo.sf2) # no NAs.
		stopifnot(nrow(fishes)==nrow(combo.sf2))

	## save as shp and csv
		st_write(combo.sf2,'data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_aug23.shp',driver='ESRI Shapefile',append=FALSE)
		st_write(combo.sf2,'resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_aug23.csv',driver='CSV',delete_layer=TRUE,delete_dsn=TRUE)

	#### STANDARDISE & MEAN-CENTRED
		data <- st_as_sf(st_read('resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_aug23.shp'))
		data2 <- data%>%
			mutate(standard.hexshark=((mtd5_PrUSE-mean(mtd5_PrUSE))/sd(mtd5_PrUSE)))%>% 
			mutate(standard.hexfish=((fishy-mean(fishy))/sd(fishy)))%>%
			mutate(standard.hexlds=((prp_lds-mean(prp_lds))/sd(prp_lds)))%>% 
			mutate(standard.hexmds=((prp_mds-mean(prp_mds))/sd(prp_mds)))%>% 
			mutate(standard.hexmds=((prp_hds-mean(prp_hds))/sd(prp_hds)))%>% 
			mutate(standard.hexdist2shore=((as.numeric(dist2shore)-mean(as.numeric(dist2shore)))/sd(as.numeric(dist2shore))))%>% 
			mutate(standard.hexdistcmg=((as.numeric(dist_cmg)-mean(as.numeric(dist_cmg)))/sd(as.numeric(dist_cmg))))%>%
			dplyr::select(jcode,standard.hexshark,standard.hexfish,standard.hexdist2shore,standard.hexdistcmg,standard.hexlds,standard.hexmds)
		stopifnot(nrow(combo.sf2)==nrow(data2))
		summary(data2)
		hist(data2$standard.hexdistcmg)

	####  reinstated 21/09/2023 - see notes in Phd notebook4.2, 15/8/2023 Principal component analysis for seagrass (hex)
		data2 <- read.csv('data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_aug23.csv')%>%mutate(jcode=as.numeric(jcode)) # for data import on 21/09/2023
		# examine correlation 
		cor(data2$standard.hexlds,data2$standard.hexmds) # - 0.454, moderate correlation

		# calculate the principal components - for BRTS, high got dropped, so medium and low are relevant.
		sapply(data2,class)
		data3<-as_tibble(data2)

		hexpca1<-prcomp(data3[,6:7])
         
		summary(hexpca1)	 	
		plot(hexpca1,type='lines') # pc2 eigen value is below the threshodl for meaningness. 
		hexpca1$sdev^2 # eigenvalues for each component
		hexpca1$rotation # look at loading coefficiations, eignenvectors. Loadings illustrate the association between each PC and the original vars. high loading = more contribution
		biplot(hexpca1,cex=0.5) # angle between them is greater than 90º indicating a negative association between the variables (makes sense)

		# PCA1 captures more of the variance, the variables have a negative association 

		data2$standard.hexsgPCA1<-predict(hexpca1)[,1]
		summary(data2)	


	## save

		st_write(data2,'resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_aug23.shp',driver='ESRI Shapefile',append=FALSE)
		st_write(data2,'resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_aug23.csv',driver='CSV',delete_layer=TRUE,delete_dsn=TRUE)


	## adding dist2urban 
		# load in df2s (fishy hexagons, standardised and mean-centred)
		df2<-read.csv('resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_aug23.csv')
		df2sf<-st_as_sf(st_read('resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_aug23.shp')) %>%
			rename(standard.hexshark=stndrd_hxs,standard.hexfish=stndrd_hxf,standard.hexdist2shore=stndr_2,standard.hexdistcmg=stndrd_hxd,standard.hexlds=stndrd_hxl,standard.hexmds=stndrd_hxm)
		
		# import Jetty locations KML file

		j<-st_zm(st_read('boatlaunchesBimini.kml'),crs='WGS84')%>%
		mutate(jetty=seq(1,15,1),.before='geometry')%>%
		dplyr::select(-Description) # jetties visible from Google Earth in Bimini. NOT including small single residence jetties. 

		nj1<-tibble(jcode=as.character(df2sf$jcode),nj=st_nearest_feature(df2sf,j))
		nj1sf<-left_join(df2sf,nj1,by='jcode')%>%
			rename(jetty=nj,jcode_geo=geometry)

			j3<-as.data.frame(nj1sf)%>%dplyr::select(jetty)
			j4<-left_join(j,j3,multiple='all',by='jetty')%>%rename(jID=jetty,jetty_geo=geometry) # jetties with geos in order of jcode for easy binding. an SF
			st_geometry(j4)

			nj2sf2<-bind_cols(nj1sf,j4)%>%
				dplyr::select(-jID)
			nj2sf2b<-nj2sf2%>%mutate(center=st_centroid(st_geometry(st_make_valid(nj2sf2))))
			summary(nj2sf2b)

			nj2sf3<-nj2sf2b%>%
				mutate(dist2jetty=as.numeric(st_distance(center,j4$jetty_geo, by_element=TRUE)),.before=jcode_geo)
				summary(nj2sf3)
				st_geometry(nj2sf3) # should be a multipolygon (aka the hexagon)


		## update file to include jetty and dist2jetty information
		st_write(nj2sf3,'resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_aug23.csv',driver='CSV',delete_layer=TRUE,delete_dsn=TRUE)	
		st_write(nj2sf3,'resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_aug23.shp',driver='ESRI Shapefile',append=FALSE)

	## deltaTemp - by nearest receiver 

	## 


	## DEPRECATED relProbPD and depthReceiver by nearest receiver, and standardise dist2jetty
		fd <- st_as_sf(st_read('resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_aug23.shp'),crs='WGS84') %>%
			rename(standard.hexshark=stndrd_hxs,standard.hexfish=stndrd_hxf,standard.hexdist2shore=stndrd_h2,standard.hexdistcmg=stndrd_hxd,standard.hexlds=stndrd_hxl,standard.hexmds=stndrd_hxm)
		hp <- st_as_sf(st_read('lemonspredators_20192020blacktipsANDbulls/relativepredatorrisk_at_receivers_April2019December2020_lemonsANDblacktips.shp'),crs='WGS84')
		d
		
		data2 <- st_join(fd,hp,left=TRUE,join=st_nearest_feature)%>%
				dplyr::select(-rlPrpPD,-rID,-total,-ndetts)%>%
				mutate(standard.relPropPD=((as.numeric(relPropPD)-mean(as.numeric(relPropPD)))/sd(as.numeric(relPropPD))),.before='geometry')%>%
				mutate(standard.dist2jetty=((as.numeric(dist2jetty)-mean(as.numeric(dist2jetty)))/sd(as.numeric(dist2jetty))),.before='geometry')%>%
				st_join(.,d,left=TRUE,join=st_nearest_feature)%>%
				mutate(standard.depth=((as.numeric(depth)-mean(as.numeric(depth)))/sd(as.numeric(depth))),.before='geometry')%>%
				dplyr::select(-locatin,-intervl,-lctn_cd,-depth,-rID,-id,-FID,-relPropPD,-dist2jetty)

		st_write(data2,'resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_aug23.shp',driver='ESRI Shapefile',append=FALSE)
		st_write(data2,'resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_aug23.csv',driver='CSV',delete_layer=TRUE,delete_dsn=TRUE)


	## 16 November 2023: new 'pressure' metric for co-occurence of juveniles and large sharks 
		## 17/11/2023 unclear where the files for the shp have gone - only the dbf remains. So re-spatially matchign using jcode and the empty hexagon grid from file 'winter2020habitat_hexagon_grid_NOland_ECdata.shp'
			hexdata<-read.csv('resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_aug23.csv')%>%mutate(jcode=as.character(jcode))
			head(hexdata)
			grid <- st_as_sf(st_read('winter2020habitat_hexagon_grid_NOland_ECdata.shp'),crs='WGS84')%>%
				select(jcode,geometry)
				grid

			## join the grid with the csv
				matched <- left_join(hexdata,grid, by='jcode')
			## convert to sf 
				hexsf <- st_as_sf(matched)
			## save
				st_write(hexsf,'resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_nov23.shp',driver='ESRI Shapefile',append=FALSE)

		## import data
		hexsf22 <- st_as_sf(st_read('resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_nov23.shp'),crs='WGS84')%>%
			rename(standard.hexshark = stndrd_hxs,
				standard.hexfish = stndrd_hxf,
				standard.hexdist2shore = stndrd_h2,
				standard.hexdistcmg = stndrd_hxd,
				standard.hexlowdensg = stndrd_hxl,
				standard.hexmeddensg = stndrd_hxm,
				standard.depth = stndrd_d,
				standard.dist2jetty = stndrd_d2,
				pressure = pressur,
				standard.sgPCA1 = st_PCA1)

        
		## join pressure data to hexsf
			## pressure metric df
			c <- st_as_sf(st_read('coocurrencemetrics_withHabitat4winter2020_fromEC_nov23.shp'),crs='WGS84')%>%
				select(pressure, geometry)
			c

			hexsf2 <- st_join(hexsf, c, join=st_nearest_feature)

		## standardise pressure 

			hexsf3 <- hexsf2 %>%
				mutate(standard.hexpress=((as.numeric(pressure)-mean(as.numeric(pressure)))/sd(as.numeric(pressure))),.before=geometry)
		########
		## save 
		########
		st_write(hexsf3,'resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_nov23.shp',driver='ESRI Shapefile',append=FALSE)
		st_write(hexsf3,'resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_nov23.csv',driver='CSV',delete_layer=TRUE,delete_dsn=TRUE)


	










##########################################
## Code graveyard
##########################################









### Df for spatial autocorrelation of hexagon grid cells 

pacman::p_load(sfdep,spdep)

data <- st_as_sf(st_read('resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_may23.shp'))%>%dplyr::select(jcode,geometry)
head(data)

# A CAR (conditional autoregressive) model requires 'queen neighborhood' information; sfdep package is supposed to help with this. 

	queen <- st_contiguity(data)
	 # neat but nimble only has tutorials for spdep (deprecated soon)
	indexneigh <- st_contiguity(data) # index of neighbours 
	weightsneigh <-st_weights(indexneigh)

	











