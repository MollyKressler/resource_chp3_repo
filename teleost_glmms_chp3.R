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
pacman::p_load(tidyverse,sf,ggplot2,gridExtra,flextable,vegan,sf,ggsn,lme4,effects,forcats)

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
	## hab.noland is missing dist2shore
	hab.noland3<-hab.noland%>%mutate(dist2shore=as.numeric(st_distance(st_centroid(hab.noland),st_union(land))),.after='dist_cmg')
	#st_write(hab.noland3,'hexagon_grid_study_site_with_habitatFromSOSF_forchp3BRUVS.shp')

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

	saveRDS(gerr,'model_RDS/glm_gerreidae_feb23.RDS')
	saveRDS(Hspp,'model_RDS/glm_shannonindex_species_feb23.RDS')
	saveRDS(Hfam,'model_RDS/glm_shannonindex_families_feb23.RDS')


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
	#saveRDS(sw_species1,'model_RDS/SW_Species_brt_tc9_lr001_gaussian.RDS')
	sw_species1<-readRDS('model_RDS/SW_Species_brt_tc9_lr001_gaussian.RDS')
	summary(sw_species1)

	# extract the data on relative influence of each variable. and plot it. 
	hab.labels<-(c('dist2shore'='Dist. to Shore (m)', 'prop_brs'='Prop. of Bare Sand','prop_ldsg'='Prop. of \n\ Low Density \n\ Seagrass','prop_medsg'='Prop. of  \n\ Medium Density \n\ Seagrass','prop_hdsg'='Prop. of  \n\ High Density \n\ Seagrass','prop_sarg'='Prop. of Sargassum','prop_urb_r'='Prop. of \n\ Urban & Rocky','prop_deep'='Prop. of \n\ Deep Water'))
	
	infl.swsp<-sw_species1$contributions
	swspp.relinf<-infl.swsp%>%mutate(var=fct_reorder(var,rel.inf))%>%ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+geom_bar(stat='identity')+scale_fill_distiller(direction=1,palette='Oranges',limits=c(0,35),guide='none')+theme_bw()+coord_flip()+scale_x_discrete(labels=hab.labels)+ylab('Relative Influence')+xlab(NULL)
		#ggsave(swspp.relinf,file='relative_influence_vars_in_SWsppBRT_feb23.png',device='png',units='in',height=6,width=8,dpi=900)


## okay let's say we're happy with that model. And we're going to skip simplifying the model. Because we want all the explanatory vars to stay in =>
## Time to plot the functions and fitted values from the model 

fittedfx.swspp1<-gbm.plot(sw_species1, plot.layout=c(3,3),write.title=FALSE) 
	# the fitted functions from the BRT model
	# these can be misleading based on the distirbution of observations (can give the wrong indication of what different values of predictor do to the response if the data is sparse at certain predictor values)
	# the next plot function, gbm.plot.fits, can help investgate if this is the case.

	# you can just save these as the base R plot. Or you can export the values from the fittedfx... and plot them in ggplot2. 

gbm.plot.fits(sw_species1)
	# this plots the fitted values irt their predictor
	# in this data set (BRUVS, SW_Species) you can see that there are a lot of zeros, but otherwise many variables have a good spread of fitted values across the range for each predictor. 
	# above each graph is a 'wtm' value. This is the 'weighted mean' of the fitted values in relation to each non-factor predictor. 

## Interactions - investgate if pairwise interactiosn exits in the data as modelled by the BRT. And then you can plot them. 

find.int.swspp<-gbm.interactions(sw_species1) # run the function
find.int.swspp$interactions # print the interacton matrix
find.int.swspp$rank.list # prints just the list of non-zero interactions. Make these into flextables. 
swpp.interactions<-find.int.swspp$rank.list%>%dplyr::select(var1.names,var2.names,int.size)%>%rename('Variable 1'='var1.names','Variable 2' = 'var2.names','Size of Interaction'=int.size)%>%flextable()%>%theme_alafoli()%>%align(align = 'center', part = 'all')%>%font(fontname = 'Times', part = 'all')%>%fontsize(size = 11, part = 'all')%>%color(color='black',part='all')%>%autofit()
	save_as_image(swpp.interactions,'interactions_above05_SWsppBRT_feb23.png',webshot='webshot')
 
swspp_intx_distandmedsg<-gbm.perspec(sw_species1,9,4) # plot pairwise interactions one at a time (or rather 2 covariates at a time with the response). Numbers correspond to the variables (easy way to identify them is the rank.list call from above). This also provides the maxmum fitted value by this interaction


## Predicting to new data. 

# dataframe to predict into = data.frame version of hab.noland. Keep the grid codes (jcode), remove geometry and check that jcode is a factor. Then duplicate the dataset with half as W and half as D (Season)
# can use df4preds wth every model. 
hab.noland<-st_as_sf(st_read('hexagon_grid_study_site_with_habitatFromSOSF_forchp3BRUVS.shp'),crs='WGS84')%>%dplyr::select(-grid_ar,-nearest_mg,-dist_cmg)

df4preds<-as.data.frame(hab.noland)%>%dplyr::select(-geometry)%>%dplyr::select(-geometry)%>%slice(rep(1:n(),each=2))%>%mutate(Season=case_when(row_number()%%2==1 ~'D',.default='W'),.after=jcode)
sf4preds<-hab.noland%>%slice(rep(1:n(),each=2))%>%mutate(Season=case_when(row_number()%%2==1 ~'D',.default='W'),.after=jcode)

	df4preds$jcode<-as.factor(df4preds$jcode)
	df4preds$Season<-as.factor(df4preds$Season)
	sapply(df4preds,class)
	head(df4preds,6)
	df4preds<-df4preds%>%dplyr::select(-nearest_mg,-dist_cmg)
	#write.csv(df4preds,'df_for_prediction_chp3BRUVS_habitat_noland.csv')

# now the data frame is ready to predict with. not spatial, site names (jcode) with predictor variables. 

preds_SWspp<-predict.gbm(sw_species1,df4preds,n.trees=sw_species1$gbm.call$best.trees,type='response')
	# one item per jcode. Need to cbind.
p2spp<-as.data.frame(preds_SWspp)
preds2_SWspp<-bind_cols(df4preds,p2spp)
head(preds2_SWspp)

# match df4preds with predictions to the spatial geometry file
	# need to duplicate and add Seasons to spatial files. Then when plotting, I will facet by season. 

sf4preds<-hab.noland%>%slice(rep(1:n(),each=2))%>%mutate(Season=case_when(row_number()%%2==1 ~'D',.default='W'),.after=jcode) # use with all models 
# SW Species
spppreds.dry<-preds2_SWspp%>%filter(Season!='W')
spppreds.wet<-preds2_SWspp%>%filter(Season!='D')

sf4predsDRY<-sf4preds%>%filter(Season=='D')%>%mutate(preds_SWsp=spppreds.dry$preds_SWspp,.before='geometry')
sf4predsWET<-sf4preds%>%filter(Season=='W')%>%mutate(preds_SWsp=spppreds.wet$preds_SWspp,.before='geometry')

sf4preds<-bind_rows(sf4predsWET,sf4predsDRY)
	
# plot them
season.labels<-as_labeller(c('D'='Dry','W'='Wet'))

swspp.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=preds_SWsp),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0.5,1.7),guide=guide_colourbar(title=' Species \n\ Diversity \n\ Index'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
ggsave(swspp.plot,file='shannon_index_species_BRTpreds_feb23.png',device='png',units='in',height=8,width=10,dpi=1000)






