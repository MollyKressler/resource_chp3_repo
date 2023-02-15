## Models of teleost biodiversity distirbution across study area by habitat type 
#### initial modelling done with habtat data from SOSF pdf. 
#### models will need to be re-run with new year specific habitat data from Emily Courmier (when, TBD) - habitat will be matched to year appropriate remote sensing maps by source (HG vs SD, 2014 for HG and 2019 for SD)

## created 8 February 2023, by Molly Kressler
## updated 13 & 14 February 2023 - BRTs example for gbm.auto package. I need to be able to run the models and export the surfaces (spatially formatted predictions). I'll test the package functionality usign the dataset they provide. 


############ ############ ############ ############ ############ 
## NOTES ### ############ ############ ############ ############

# colors for figures: by season, wet = cadetblue3 and dry = tomato4


############ ############ ############ ############ ############ 
# DATASETS, load at the start ######## ############ ############ 

setwd('/Users/mollykressler/Documents/data_phd/resource_chp3')
pacman::p_load(tidyverse,sf,ggplot2,gridExtra,flextable,vegan,sf,ggsn,lme4,ggeffects,effects,beepr,modelsummary,performance)

# SHAPEFILE: bruvs data and habitat data (using sep22 updated habitat SOSF) joined together, with dist2shore calculated
	joined<-st_as_sf(st_read('bruvs_data_joinedWITHhabitat_feb23.shp'))%>%rename(spp_abundance=spp_bnd,spp_richness=spp_rch,SW_Species=SW_Spcs,SW_Families=SW_Fmls,prop_brs=prp_brs,prop_ldsg=prp_lds,prop_medsg=prp_mds,prop_hdsg=prp_hds,prop_sarg=prp_srg,prop_urb_r=prp_rb_,prop_deep=prop_dp,dist2shore=dst2shr,Sphyraenidae=Sphyrnd,Scaridae=Scarida,Haemulidae=Haemuld,Gerreidae=Gerreid,Belonidae=Belonid,Sparidae=Sparida)%>%dplyr::select(BRUV,jcode,spp_abundance,spp_richness,SW_Species,SW_Families,Sphyraenidae,Scaridae,Haemulidae,Gerreidae,Belonidae,Sparidae,prop_brs,prop_ldsg,prop_medsg,prop_hdsg,prop_sarg,prop_urb_r,prop_deep,dist2shore)

# DATE FRAME: bruvs data and habitat data (using sep22 updated habitat SOSF) joined together, with dist2shore calculated. Use in models. 

	joined_df<-read.csv('bruvs_DATAFRAME_joinedWITHhabitat_feb23.csv')%>%dplyr::select(-X)
	joined_df$BRUV<-as.factor(joined_df$BRUV)
	joined_df$Season<-as.factor(joined_df$Season)
	joined_df$SW_Species<-as.numeric(joined_df$SW_Species)
	joined_df$SW_Families<-as.numeric(joined_df$SW_Families)

# SHAPEFILE(S): land and water with grid

	setwd('/Users/mollykressler/Documents/data_phd')
	hab.grid<-st_as_sf(st_read('hexagon_grid_withUPDATEDhabitat_sept22.shp'),crs='WGS84')
	setwd('/Users/mollykressler/Documents/data_phd/resource_chp3')
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')

	hab.noland<-st_as_sf(st_difference(hab.grid,st_union(land))) # cut land out of hab.grid. 
	hab.nloand2<-st_cast(hab.noland,'MULTIPOLYGON')
	hab.atcentroids<-hab.nloand2%>%add_column(Longitude=st_coordinates(st_centroid(.))[,1],.before='prop_brs')%>%add_column(Latitude=st_coordinates(st_centroid(.))[,2],.before='Longitude')

############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ 

#### 9 - 10 February 2023 - (A) glm/glmms for response variables. 

### (A) GLMS for response vars: shannon index for families and species, and Gerreidae family counts. Explore some model selection, by dropping certain habitats. 
	## poisson with log link for Gerreidae
	## for SWs? gamma? with negative inverse link. 

	sapply(joined_df,class)
	summary(joined_df)

	# Shannon Index for Species
	Hspp<-glm(SW_Species~Season+dist2shore+prop_brs+prop_medsg+prop_sarg+prop_urb_r+prop_deep,data=joined_df,family=gaussian(link='identity'));beep(2)
	summary(Hspp)
	# Shannon Index for Families
	Hfam<-glm(SW_Families~Season+dist2shore+prop_brs+prop_medsg+prop_sarg+prop_urb_r+prop_deep,data=joined_df,family=gaussian(link='identity'));beep(2)
	summary(Hfam)

	# Gerreidae
	gerr<-glm(Gerreidae~Season+dist2shore+prop_brs+prop_medsg+prop_sarg+prop_urb_r+prop_deep,data=joined_df,family=poisson(link='log'));beep(2)
	summary(gerr)


	# save model RDS 

	saveRDS(gerr,'model_outputs/glm_gerreidae_feb23.RDS')
	saveRDS(Hspp,'model_outputs/glm_shannonindex_species_feb23.RDS')
	saveRDS(Hfam,'model_outputs/glm_shannonindex_families_feb23.RDS')


############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ 
			## Boosted Regression Tree Analysis ##
	## for Gerreidae family, SW Species, and SW Families ##
############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ 

## using this tutorial [https://rspatial.org/raster/sdm/9_sdm_brt.html#predicting-to-new-data], which is in associaton with the Elith et al. 2008 paper. 
	# not using gbm.auto. I tried and we didn't get along. 

		#sw1<-gbm.auto(grids=NULL,samples=samples,expvar=c( ,13:20),resvar=c('SW_Species'),tc=c(9),lr=c(0.001,0.001),bf=0.9,fam1='bernoulli',fam2='gaussian',savedir='/Users/mollykressler/Documents/data_phd/resource_chp3',savegbm=TRUE,varint=TRUE)

		#sw2<-gbm.loop(loops=10,samples=samples,expvar=c(2,13:20),resvar=c('SW_Species'),tc=c(9),lr=c(0.001,0.001),bf=0.9,fam1='bernoulli',fam2='gaussian',savedir='/Users/mollykressler/Documents/data_phd/resource_chp3',calcpreds=TRUE,savegbm=TRUE,multiplot=TRUE,)

## going to work through the tutorial top to bottom for SW Speces to see if I can get what I want: diversity index predictions into the hab.noland grid.

pacman::p_load(gbm,tidyverse,sf,ggplot2,dismo)

## using gbm.bfcheck to rule out bag fractions 
sw.fam.bf<-gbm.bfcheck(samples,resvar=c('SW_Families')) #  binary bag fraction must be at least 0.087. n = 242 "Gaussian bag fraction must be at least 0.284. n = 74"
sw.spp.bf<-gbm.bfcheck(samples,resvar=c('SW_Species')) # binary bag fraction must be at least 0.087. n = 242". "Gaussian bag fraction must be at least 0.091. n = 231"
gerr.bf<-gbm.bfcheck(samples,resvar=c('Gerreidae')) ## 0.009, extremely low. should be able to model without issues related to bf. 

## use gbm.step() to both cross validate the number of trees to use, and then run that model. 

	sw1<-gbm.step(joined_df,gbm.x=c('Season','prop_brs','prop_ldsg','prop_medsg','prop_hdsg','prop_sarg','prop_urb_r','prop_deep','dist2shore'),gbm.y=c('SW_Species'),tree.complexity=9,learning.rate=0.001,bag.fraction=0.75,family='gaussian',plot.main = TRUE) # 950 trees
	sw_species1<-sw1
	names(sw1)
	summary(sw1) # plots, and tables, the influence of each variable. 
	#saveRDS(sw_species1,'model_outputs/SW_Species_brt_tc9_lr001_gaussian.RDS')

## okay let's say we're happy with that model. And we're going to skip simplifying the model. Because we want all the explanatory vars to stay in =>
## Time to plot the functions and fitted values from the model 

gbm.plot(sw_species1, plot.layout=c(3,3),write.title=FALSE) 
	# the fitted functions from the BRT model
	# these can be misleading based on the distirbution of observations (can give the wrong indication of what different values of predictor do to the response if the data is sparse at certain predictor values)
	# the next plot function, gbm.plot.fits, can help investgate if this is the case.

gbm.plot.fits(sw_species1)
	# this plots the fitted values irt their predictor
	# in this data set (BRUVS, SW_Species) you can see that there are a lot of zeros, but otherwise many variables have a good spread of fitted values across the range for each predictor. 
	# above each graph is a 'wtm' value. This is the 'weighted mean' of the fitted values in relation to each non-factor predictor. 

## Interactions - investgate if pairwise interactiosn exits in the data as modelled by the BRT. And then you can plot them. 

find.int.swspp<-gbm.interactions(sw_species1) # run the function
find.int.swspp$interactions # print the interacton matrix
find.int.swspp$rank.list # prints just the list of non-zero interactions 

gbm.perspec(sw_species1,9,4) # plot pairwise interactions one at a time (or rather 2 covariates at a time with the response). Numbers correspond to the variables (easy way to identify them is the rank.list call from above). This also provides the maxmum fitted value by this interaction
gbm.perspec(sw_species1,9,6) 


## Predicting to new data. 









