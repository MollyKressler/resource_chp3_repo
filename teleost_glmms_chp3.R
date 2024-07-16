## Models of teleost biodiversity distirbution across study area by habitat type 
#### initial modelling done with habtat data from SOSF pdf. 
#### models will need to be re-run with new year specific habitat data from Emily Courmier (when, TBD) - habitat will be matched to year appropriate remote sensing maps by source (HG vs SD, 2014 for HG and 2019 for SD)

## created 8 February 2023, by Molly Kressler
## updated 13 & 14 February 2023 - BRTs example for gbm.auto package. I need to be able to run the models and export the surfaces (spatially formatted predictions). I'll test the package functionality usign the dataset they provide. 
## updated 10 August 2023, new habitat data from Emily COurmier. Need to re-run full BRTs to run simplified BRTs in 'comparing_simmplifiedBRTS_chp3_neater.R'


############ ############ ############ ############ ############ 
## NOTES ### ############ ############ ############ ############

# colors for figures: by season, wet = cadetblue3 and dry = tomato4

pacman::p_load(tidyverse,sf,ggplot2,gridExtra,flextable,vegan,sf,ggsn,lme4,effects,forcats,gbm3, dismo)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd/')


############ ############ ############ ############ ############ 
# BRTS - updated August 2023 (new hab data) & June 2024 (Driscoll BRUVs only) & July 2024 (publically available 2014 and Driscoll 2018 with updated habitat specific to sampling year), copy and pasted from before
############ ############ ############ ############ ############ 

	n <-read.csv('bruvs2014AND2018_data_joinedWITHhabitat_summer18habANDwinter2014hab_july2024.csv')
	n


## RUN MODELS (or READ IN model objects)
## use gbm.step() to both cross validate the number of trees to use, and then run that model. 
	# Poisson because count data


	n1 <-gbm.step(n,gbm.x=c('prp_lds','prp_mds','prp_hds','dist2shore'),gbm.y='maxN',tree.complexity=2,learning.rate=0.001,bag.fraction=0.95,family='poisson',plot.main = TRUE)  ## resid dev 16.4, 3700 trees
	
	summary(n1)

	saveRDS(n1,'resource_chp3/model_RDS/maxN_gbm_poisson_july2024.RDS')


	## DEPRECATED - Species richness is also count data so also needs a poisson
		rich2<-gbm.step(joined_df,gbm.x=c('prp_lds','prp_mds','prp_hds','dist2shore'),gbm.y=c('spp_rch'),tree.complexity=5,learning.rate=0.0001,bag.fraction=0.95,family='poisson',plot.main = TRUE)  ## resid dev 1.8, 1750 trees, NO Season.
		rich1<-gbm.step(joined_df,gbm.x=c('prp_lds','prp_mds','prp_hds','dist2shore'),gbm.y=c('spp_rch'),tree.complexity=5,learning.rate=0.0001,bag.fraction=0.5,family='poisson',plot.main = TRUE)  ## resid dev 1.8, 4050 trees, NO Season.
		
		#rich3<-gbm.step(joined_df,gbm.x=c('prp_lds','prp_mds','prp_hds','dist2shore'),gbm.y=c('spp_rch'),tree.complexity=2,learning.rate=0.0001,bag.fraction=0.5,family='poisson',plot.main = TRUE)  ## resid dev 1.8, 4050 trees, NO Season.

		saveRDS(rich3,'resource_chp3/model_RDS/Species_Rchness_Poisson_brt_tc3_lr0001_poisson_NoSeason_jun24.RDS')



	# bag fraction = the fraction of the training set observations randomly selected to propose the next tree in the expansion
## output of gbm.step explained 
	# $fitted = the fitted values form the final tree. on the response scale. 
	# $fitted.vars = the variance of the fitted values, on the resposne scale
	# $residuals = the residuals of the fitted values, on the response scale. 
	# $contributons  = relative importance of the variables. 
	# $self.statistics = evaluation statistics, calculated on fitted values. Shows 'fit' of the model on the training data. NOT to be reported as model performance. 
	# $cv.statistics = most appropriate for evaluation. 
	# $weights = the weights used in fitting the models. defualt is 1. 
	# cv.values = the mean of CV (cross validating) estmates of predictive deviance. Calculated at each stagewise step. This is used in the auto-generated plots of tress versus deviance. 
	# $cv.loss.ses = standard errors in CV estimates of predictive deviance at each step in the stagewise process
	# $cv.roc.matrix = matrix of the values for area under the curve estimated on the excluded data, instead of deviance in the cv.loss.matrix.


	# read in model objects 
		n1 <- readRDS('resource_chp3/model_RDS/maxN_gbm_poisson_july2024.RDS')

		# deprecated
			richfull <- readRDS('resource_chp3/model_RDS/spp_richness_brt_tc5_lr001_gaussian_jun24.RDS')
			gerrfull <- readRDS('resource_chp3/model_RDS/gerrPoisson_brt_tc3_lr0001_poisson_NoSeason_jun24.RDS')
			sppfull<-readRDS('resource_chp3/model_RDS/SW_Species_brt_tc3_lr001_gaussian_aug23.RDS')
			famfull<-readRDS('resource_chp3/model_RDS/SW_Families_brt_tc9_lr001_gaussian_aug23.RDS')
			gerrfull<-readRDS('resource_chp3/model_RDS/gerrPoisson_brt_tc3_lr0001_poisson_aug23.RDS')
##########
## EVALUATE MODELS
# extract the data on relative influence of each variable. and plot it. 
	hab.labels<-(c('dist2shore'='Dist. to Shore (m)', 'prp_lds'='Prop. of \n\ Low Density \n\ Seagrass','prp_mds'='Prop. of  \n\ Medium Density \n\ Seagrass','prp_hds'='Prop. of  \n\ High Density \n\ Seagrass'))
	# maxN 
		infl.n<-n1$contributions
		n.relinf<-infl.n%>%
			mutate(var=fct_reorder(var,rel.inf))%>%
			ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+
				geom_bar(stat='identity')+
				scale_fill_distiller(direction=1,palette='Blues',limits=c(15,35),guide='none')+
				theme_bw()+
				coord_flip()+
				scale_x_discrete(labels=hab.labels)+
				ylab('Relative Influence')+
				xlab(NULL)+
				theme(text = element_text(size = 14))

		ggsave(n.relinf,file='resource_chp3/BRTS_outputs/BRT_maxN_july2024/relative_influence_vars_maxN_july2024.png',device='png',units='in',height=4,width=5.5,dpi=900)
		
		res_n1 <- as_tibble(resid(n1))%>%rename(resids = value)
		F1 <- as_tibble(predict(n1))%>%rename(fitted = value)
		diag_n1 <- bind_cols(F1, res_n1)

		resVfit_n1 <- ggplot(data = diag_n1,aes(x= fitted, y = resids))+
			geom_point(col = '#3A6C74', fill = '#3A6C74')+ 
			labs(subtitle = 'Residuals vs. Fitted', y = 'Residuals', x = 'Fitted')+
			theme_bw()
		res_n1_hist <- ggplot(data = res_n1, aes(x = resids))+ 
			geom_histogram(binwidth = nrow(res_n1)/100,col = '#3A6C74', fill = '#3A6C74')+ 
			theme_bw()+
			labs(subtitle = 'Residuals', y = 'Frequency', x = 'Residuals')
		res_n1_qq <- ggplot(data=F1, aes(sample = fitted))+
			stat_qq(size=1,pch=21,col = '#3A6C74', fill = '#3A6C74')+
			labs(subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
			stat_qq_line(linetype=2, col='red')+
			theme_bw()

		diagnostics_n1 <- resVfit_n1+res_n1_hist+res_n1_qq
		
		ggsave(diagnostics_n1, file = 'resource_chp3/BRTS_outputs/BRT_maxN_july2024/diagnostics_BRT_maxN_july2024.png', device = 'png', unit = 'in', height = 4, width = 8, dpi = 850)	
	# deprecated - Gerr 
		infl.gerr<-gerrfull$contributions
		gerr.relinf<-infl.gerr%>%
			mutate(var=fct_reorder(var,rel.inf))%>%
			ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+
				geom_bar(stat='identity')+
				scale_fill_distiller(direction=1,palette='Oranges',limits=c(0,80),guide='none')+
				theme_bw()+
				coord_flip()+
				scale_x_discrete(labels=hab.labels)+
				ylab('Relative Influence')+
				xlab(NULL)	
		ggsave(gerr.relinf,file='resource_chp3/figures+tables/BRT_gerreidae_jun24/relative_influence_vars_in_GerreidaeBRT_jun24.png',device='png',units='in',height=6,width=8,dpi=900)
	# deprecated - Richness
		infl.rich<-richfull$contributions
		rich.relinf<-infl.rich%>%
			mutate(var=fct_reorder(var,rel.inf))%>%
			ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+
				geom_bar(stat='identity')+
				scale_fill_distiller(direction=1,palette='Blue',limits=c(0,50),guide='none')+
				theme_bw()+
				coord_flip()+
				scale_x_discrete(labels=hab.labels)+
				ylab('Relative Influence')+
				xlab(NULL)	
		ggsave(gerr.relinf,file='resource_chp3/figures+tables/BRT_species_richness_jun24/relative_influence_vars_in_sppRichnessBRT_jun24.png',device='png',units='in',height=6,width=8,dpi=900)
	# deprecated - SW Species 
		infl.swsp<-sw_species2$contributions
		swspp.relinf<-infl.swsp%>%
			mutate(var=fct_reorder(var,rel.inf))%>%
			ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+
				geom_bar(stat='identity')+scale_fill_distiller(direction=1,palette='Oranges',limits=c(0,65),guide='none')+
				theme_bw()+
				coord_flip()+
				scale_x_discrete(labels=hab.labels)+
				ylab('Relative Influence')+
				xlab(NULL)
			ggsave(swspp.relinf,file='resource_chp3/figures+tables/BRT_sw_species_feb23/relative_influence_vars_in_SWsppBRT_aug23.png',device='png',units='in',height=6,width=8,dpi=900)

	# deprecated - SW Families 
		infl.swf<-sw_families2$contributions
		swf.relinf<-infl.swf%>%
			mutate(var=fct_reorder(var,rel.inf))%>%
				ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+
				geom_bar(stat='identity')+
				scale_fill_distiller(direction=1,palette='Oranges',limits=c(0,70),guide='none')+
				theme_bw()+
				coord_flip()+
				scale_x_discrete(labels=hab.labels)+ylab('Relative Influence')+
				xlab(NULL)
			ggsave(swf.relinf,file='resource_chp3/figures+tables/BRT_sw_families_feb23/relative_influence_vars_in_SWfamiliesBRT_aug23.png',device='png',units='in',height=6,width=8,dpi=900)

	# deprecated - Gerreidae
		infl.gerr<-gerr2$contributions
		gerr.relinf<-infl.gerr%>%
			mutate(var=fct_reorder(var,rel.inf))%>%
			ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+
				geom_bar(stat='identity')+
				scale_fill_distiller(direction=1,palette='Oranges',limits=c(0,70),guide='none')+
				theme_bw()+
				coord_flip()+
				scale_x_discrete(labels=hab.labels)+
				ylab('Relative Influence')+
				xlab(NULL)	
		ggsave(gerr.relinf,file='resource_chp3/figures+tables/BRT_gerreidae_feb23/relative_influence_vars_in_GerreidaeBRT_aug23.png',device='png',units='in',height=6,width=8,dpi=900)


##########
## Before you predict you need to do gbm.simplify - don't go further. 

	df4preds<-read_csv('winter2020habitat_hexagon_grid_NOland_ECdata.csv')

	sf4preds<-st_as_sf(st_read('winter2020habitat_hexagon_grid_NOland_ECdata.shp'),crs='WGS84')

full.models<- n1 

for(i in full.models){
	mod<-get(i)
	predicts<-predict.gbm(mod,df4preds,n.trees=mod$gbm.call$best.trees,type='response')
	predicts<-as.data.frame(predicts)
	predicts.df<-bind_cols(df4preds,predicts)

	dry<-predicts.df%>%filter(Season!='W')
	wet<-predicts.df%>%filter(Season=='W')
	sf.dry<-sf4preds%>%filter(Season=='D')%>%mutate(preds=dry$predicts)
	sf.wet<-sf4preds%>%filter(Season=='W')%>%mutate(preds=wet$predicts)
	sf4preds3<-bind_rows(sf.dry,sf.wet)
	new=paste0(i,'PD')
	sf4preds<-sf4preds3%>%dplyr::rename_with(.fn=~paste0(new),.cols='preds')
	}

	# back transform the gerreidae (poisson, log)
	sf4preds<-sf4preds%>%mutate(bt.gerr.fullPD=exp(.$gerrfullPD),.before='geometry')

	# save the updated sf4preds and df4preds
		st_write(sf4preds,'sf_for_predictions_fromBRTs_aug23.shp',append=FALSE) # append command set to FALSE replaces the old file with the new one. If set to TRUE it would add the new data as further layers.
		st_write(sf4preds,'sf_for_predictions_fromBRTs_aug23.csv',delete_layer=TRUE,delete_dsn=TRUE)









############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ 
############ ############CODE GRAVEYARD############ ############ 
############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ 
## deprecated, pre-August 2023 (see  below for August 2023 updates)
# DATASETS, load at the start ######## ############ ############ 

############ ############ ############ ############ ############ 
pacman::p_load(tidyverse,sf,ggplot2,gridExtra,flextable,vegan,sf,ggsn,lme4,effects,forcats,gbm3)

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

	hab.noland<-st_as_sf(st_read('hexagon_grid_study_site_with_habitatFromSOSF_forchp3BRUVS.shp'),crs='WGS84')
		# how I made it: <-st_as_sf(st_difference(hab.grid,st_union(land))) # cut land out of hab.grid. 
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

## RUN MODELS (or READ IN model objects)
## use gbm.step() to both cross validate the number of trees to use, and then run that model. 

	# SW species
	sw1<-gbm.step(joined_df,gbm.x=c('Season','prop_brs','prop_ldsg','prop_medsg','prop_hdsg','prop_sarg','prop_urb_r','prop_deep','dist2shore'),gbm.y=c('SW_Species'),tree.complexity=9,learning.rate=0.001,bag.fraction=0.75,family='gaussian',plot.main = TRUE) # 950 trees
	sw_species2<-gbm.step(joined_df,gbm.x=c('Season','prop_brs','prop_ldsg','prop_medsg','prop_hdsg','prop_sarg','prop_urb_r','prop_deep','dist2shore'),gbm.y=c('SW_Species'),tree.complexity=3,learning.rate=0.001,bag.fraction=0.75,family='gaussian',plot.main = TRUE) # 3000 trees and residual dev (mean) at 0.209 at tc=3
		sw_species1<-sw1
	names(sw_species2)
	summary(sw_species2$n.trees)
	summary(sw1) # plots, and tables, the influence of each variable. 
	#saveRDS(sw_species1,'model_RDS/SW_Species_brt_tc9_lr001_gaussian.RDS')
	#saveRDS(sw_species2,'model_RDS/SW_Species_brt_tc3_lr001_gaussian.RDS')
	
	# SW Families 
	sw_families1<-gbm.step(joined_df,gbm.x=c('Season','prop_brs','prop_ldsg','prop_medsg','prop_hdsg','prop_sarg','prop_urb_r','prop_deep','dist2shore'),gbm.y=c('SW_Species'),tree.complexity=9,learning.rate=0.001,bag.fraction=0.75,family='gaussian',plot.main = TRUE) # 950 trees
	sw_families2<-gbm.step(joined_df,gbm.x=c('Season','prop_brs','prop_ldsg','prop_medsg','prop_hdsg','prop_sarg','prop_urb_r','prop_deep','dist2shore'),gbm.y=c('SW_Families'),tree.complexity=3,learning.rate=0.001,bag.fraction=0.75,family='gaussian',plot.main = TRUE) # 1100 trees
	names(sw_families1)
	summary(sw_families1)
	summary(sw_families2)
	plot(sw_families2$residuals)
	# test out simplifying 
		#sw_families.simp<-gbm.simplify(sw_families1,n.drops=7) # suggests removing 1 var. 
		#sw_families1.reduced<-gbm.step(joined_df,gbm.x=sw_families.simp$pred.list[[1]],gbm.y=c('Gerreidae'),tree.complexity=3,learning.rate=0.0001,bag.fraction=0.5,family='poisson',plot.main = TRUE) # 950 trees
		# I dont think the model improves much really. 
	#saveRDS(sw_families1,'model_RDS/SW_Families_brt_tc9_lr001_gaussian.RDS')
	#saveRDS(sw_families2,'model_RDS/SW_Families_brt_tc3_lr001_gaussian.RDS')

	# Gerreidae (counts)  - predicting on the total wrong scale...
	#gerr1<-gbm.step(joined_df,gbm.x=c('Season','prop_brs','prop_ldsg','prop_medsg','prop_hdsg','prop_sarg','prop_urb_r','prop_deep','dist2shore'),gbm.y=c('Gerreidae'),tree.complexity=2,learning.rate=0.0001,bag.fraction=0.65,family='poisson',plot.main = TRUE) # 950 trees
	gerr2<-gbm.step(joined_df,gbm.x=c('Season','prop_brs','prop_ldsg','prop_medsg','prop_hdsg','prop_sarg','prop_urb_r','prop_deep','dist2shore'),gbm.y=c('Gerreidae'),tree.complexity=3,learning.rate=0.0001,bag.fraction=0.65,family='poisson',plot.main = TRUE) # 950 trees
	names(gerr2)
	summary(gerr2)
	gerr2$n.trees
		# really broad CIs in defining the optimal fit.

	#saveRDS(gerr2,'model_RDS/Gerreidae_brt_tc3_lr0001_poisson.RDS') 

	gerr2<-readRDS('model_RDS/Gerreidae_brt_tc3_lr0001_poisson.RDS')
	gerr2$gbm.call$


	# Gerreidae - 27 February 2023 - as gamma distributed. Gamma distributions don't allow for zeros, so add a constant to all data points for Gerreidae. You also have to use the github version of gbm to access this update (to run a gamma distribution)
	
		# gbm.step is from dismo. 
		# need to use gbm function from gbm3/gbm to use gamma distribution

	library('devtools')
	install_github("gbm-developers/gbm3")
	pacman::p_load(gbm3)
	available_distributions()

	joined_df<-joined_df%>%mutate(Gerr_Gamma=Gerreidae+0.01)

	gerr3<-gbm(joined_df,gbm.x=c('Season','prop_brs','prop_ldsg','prop_medsg','prop_hdsg','prop_sarg','prop_urb_r','prop_deep','dist2shore'),gbm.y=c('Gerr_Gamma'),tree.complexity=3,learning.rate=0.01,bag.fraction=0.65,family='Gamma',plot.main = TRUE) 
	
	gerr3<-gbm(Gerr_Gamma~Season+prop_brs+prop_ldsg+prop_medsg+prop_hdsg+prop_sarg+prop_urb_r+prop_deep+dist2shore, data=joined_df,n.trees=1000,interaction.depth=3,train.fraction=0.3,shrinkage=0.01,bag.fraction=0.65,distribution='Gamma',keep.data = TRUE) 

	gerr3$variables
	summary.gbm(gerr3)


	#saveRDS(sw_families1,'model_RDS/SW_Families_brt_tc9_lr001_gaussian.RDS')

	# bag fraction = the fraction of the training set observations randomly selected to propose the next tree in the expansion
## output of gbm.step explained 
	# $fitted = the fitted values form the final tree. on the response scale. 
	# $fitted.vars = the variance of the fitted values, on the resposne scale
	# $residuals = the residuals of the fitted values, on the response scale. 
	# $contributons  = relative importance of the variables. 
	# $self.statistics = evaluation statistics, calculated on fitted values. Shows 'fit' of the model on the training data. NOT to be reported as model performance. 
	# $cv.statistics = most appropriate for evaluation. 
	# $weights = the weights used in fitting the models. defualt is 1. 
	# cv.values = the mean of CV (cross validating) estmates of predictive deviance. Calculated at each stagewise step. This is used in the auto-generated plots of tress versus deviance. 
	# $cv.loss.ses = standard errors in CV estimates of predictive deviance at each step in the stagewise process
	# $cv.roc.matrix = matrix of the values for area under the curve estimated on the excluded data, instead of deviance in the cv.loss.matrix.


	# read in model objects 
	sw_species2<-readRDS('model_RDS/SW_Species_brt_tc3_lr001_gaussian.RDS')
	sw_families2<-readRDS('model_RDS/SW_Families_brt_tc3_lr001_gaussian.RDS')
	gerr2<-readRDS('model_RDS/Gerreidae_brt_tc3_lr0001_poisson.RDS')

	sw_families2$gbm.call$gbm.y

## EVALUATE MODELS
# extract the data on relative influence of each variable. and plot it. 
	hab.labels<-(c('dist2shore'='Dist. to Shore (m)', 'prop_brs'='Prop. of Bare Sand','prop_ldsg'='Prop. of \n\ Low Density \n\ Seagrass','prop_medsg'='Prop. of  \n\ Medium Density \n\ Seagrass','prop_hdsg'='Prop. of  \n\ High Density \n\ Seagrass','prop_sarg'='Prop. of Sargassum','prop_urb_r'='Prop. of \n\ Urban & Rocky','prop_deep'='Prop. of \n\ Deep Water'))
	# SW Species 
	infl.swsp<-sw_species2$contributions
	swspp.relinf<-infl.swsp%>%mutate(var=fct_reorder(var,rel.inf))%>%ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+geom_bar(stat='identity')+scale_fill_distiller(direction=1,palette='Oranges',limits=c(0,35),guide='none')+theme_bw()+coord_flip()+scale_x_discrete(labels=hab.labels)+ylab('Relative Influence')+xlab(NULL)
		ggsave(swspp.relinf,file='relative_influence_vars_in_SWsppBRT_feb23.png',device='png',units='in',height=6,width=8,dpi=900)

	# SW Families 
	infl.swf<-sw_families2$contributions
	swf.relinf<-infl.swf%>%mutate(var=fct_reorder(var,rel.inf))%>%ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+geom_bar(stat='identity')+scale_fill_distiller(direction=1,palette='Oranges',limits=c(0,35),guide='none')+theme_bw()+coord_flip()+scale_x_discrete(labels=hab.labels)+ylab('Relative Influence')+xlab(NULL)
		ggsave(swf.relinf,file='figures+tables/BRT_sw_families_feb23/relative_influence_vars_in_SWfamiliesBRT_feb23.png',device='png',units='in',height=6,width=8,dpi=900)

	# Gerreidae
	infl.gerr<-gerr2$contributions
	gerr.relinf<-infl.gerr%>%mutate(var=fct_reorder(var,rel.inf))%>%ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+geom_bar(stat='identity')+scale_fill_distiller(direction=1,palette='Oranges',limits=c(0,35),guide='none')+theme_bw()+coord_flip()+scale_x_discrete(labels=hab.labels)+ylab('Relative Influence')+xlab(NULL)
		ggsave(gerr.relinf,file='relative_influence_vars_in_GerreidaeBRT_feb23.png',device='png',units='in',height=6,width=8,dpi=900)


		## okay let's say we're happy with that model. And we're going to skip simplifying the model. Because we want all the explanatory vars to stay in =>

## Time to plot the functions and fitted values from the model 

	#SW Species 
	fittedfx.swspp1<-gbm.plot(sw_species2, plot.layout=c(3,3),write.title=FALSE) 
		# the fitted functions from the BRT model
		# these can be misleading based on the distirbution of observations (can give the wrong indication of what different values of predictor do to the response if the data is sparse at certain predictor values)
		# the next plot function, gbm.plot.fits, can help investgate if this is the case.

		# use base R save functionto save this. 
		

	# SW Families
	fittedfx.swf1<-gbm.plot(sw_families2, plot.layout=c(3,3),write.title=FALSE) 

	#SW Species 
	gbm.plot.fits(sw_species2)
		# this plots the fitted values irt their predictor
		# in this data set (BRUVS, SW_Species) you can see that there are a lot of zeros, but otherwise many variables have a good spread of fitted values across the range for each predictor. 
		# above each graph is a 'wtm' value. This is the 'weighted mean' of the fitted values in relation to each non-factor predictor. 
		# extract the fitted values and plot them
			a1<-plot(sw_species2,'dist2shore',return.grid=TRUE)
			a2<-plot(sw_species2,'prop_brs',return.grid=TRUE)
			a3<-plot(sw_species2,'prop_ldsg',return.grid=TRUE)
			a4<-plot(sw_species2,'prop_medsg',return.grid=TRUE)
			a5<-plot(sw_species2,'prop_hdsg',return.grid=TRUE)
			a6<-plot(sw_species2,'prop_sarg',return.grid=TRUE)
			a7<-plot(sw_species2,'prop_urb_r',return.grid=TRUE)
			a8<-plot(sw_species2,'prop_deep',return.grid=TRUE)
			a9<-plot(sw_species2,'Season',return.grid=TRUE)
			hab.colors<-c('prop_brs'='navy','prop_ldsg'='violetred4','prop_medsg'='cadetblue4','prop_hdsg'='cadetblue3','prop_sarg'='goldenrod3','prop_urb_r'='orange3','prop_deep'='hotpink2')	
			hab.labels.nodist<-(c('prop_brs'='Bare Sand','prop_ldsg'='Low Density Seagrass','prop_medsg'='Medium Density Seagrass','prop_hdsg'=' High Density Seagrass','prop_sarg'='Sargassum','prop_urb_r'='Urban & Rocky','prop_deep'=' Deep Water'))
			swspp_fittedvalues<-ggplot()+scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits=c(0.9,1.3))+geom_line(data=a2,aes(x=prop_brs,y=y,col='prop_brs',),lwd=1)+geom_line(data=a3,aes(x=prop_ldsg,y=y,col='prop_ldsg'),lwd=1)+geom_line(data=a4,aes(x=prop_medsg,y=y,col='prop_medsg'),lwd=1)+geom_line(data=a5,aes(x=prop_hdsg,y=y,col='prop_hdsg'),lwd=1)+geom_line(data=a6,aes(x=prop_sarg,y=y,col='prop_sarg'),lwd=1)+geom_line(data=a7,aes(x=prop_urb_r,y=y,col='prop_urb_r'),lwd=1)	+geom_line(data=a8,aes(x=prop_deep,y=y,col='prop_deep'),lwd=1)+scale_color_manual(values=hab.colors,labels=hab.labels.nodist)+labs(x='Proportion',y='Fitted value',color='Habitat Type')+theme_bw()
			ggsave(swspp_fittedvalues,file='swspp_fittedvalues_proportionHabTypes_BRT_feb23.png',device='png',unit='in',height=6,width=9,dpi=800)

	# SW Families
		gbm.plot.fits(sw_families2)
		# extract the fitted values and plot them
			a1<-plot(sw_families2,'dist2shore',return.grid=TRUE)
			a2<-plot(sw_families2,'prop_brs',return.grid=TRUE)
			a3<-plot(sw_families2,'prop_ldsg',return.grid=TRUE)
			a4<-plot(sw_families2,'prop_medsg',return.grid=TRUE)
			a5<-plot(sw_families2,'prop_hdsg',return.grid=TRUE)
			a6<-plot(sw_families2,'prop_sarg',return.grid=TRUE)
			a7<-plot(sw_families2,'prop_urb_r',return.grid=TRUE)
			a8<-plot(sw_families2,'prop_deep',return.grid=TRUE)
			a9<-plot(sw_families2,'Season',return.grid=TRUE)
			hab.colors<-c('prop_brs'='navy','prop_ldsg'='violetred4','prop_medsg'='cadetblue4','prop_hdsg'='cadetblue3','prop_sarg'='goldenrod3','prop_urb_r'='orange3','prop_deep'='hotpink2')	
			hab.labels.nodist<-(c('prop_brs'='Bare Sand','prop_ldsg'='Low Density Seagrass','prop_medsg'='Medium Density Seagrass','prop_hdsg'=' High Density Seagrass','prop_sarg'='Sargassum','prop_urb_r'='Urban & Rocky','prop_deep'=' Deep Water'))
			swfam_fittedvalues<-ggplot()+scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits=c(0.9,1.3))+geom_line(data=a2,aes(x=prop_brs,y=y,col='prop_brs',),lwd=1)+geom_line(data=a3,aes(x=prop_ldsg,y=y,col='prop_ldsg'),lwd=1)+geom_line(data=a4,aes(x=prop_medsg,y=y,col='prop_medsg'),lwd=1)+geom_line(data=a5,aes(x=prop_hdsg,y=y,col='prop_hdsg'),lwd=1)+geom_line(data=a6,aes(x=prop_sarg,y=y,col='prop_sarg'),lwd=1)+geom_line(data=a7,aes(x=prop_urb_r,y=y,col='prop_urb_r'),lwd=1)	+geom_line(data=a8,aes(x=prop_deep,y=y,col='prop_deep'),lwd=1)+scale_color_manual(values=hab.colors,labels=hab.labels.nodist)+labs(x='Proportion',y='Fitted value',color='Habitat Type')+theme_bw()
			ggsave(swfam_fittedvalues,file='figures+tables/BRT_sw_families_feb23/swFamilies_fittedvalues_proportionHabTypes_BRT_feb23.png',device='png',unit='in',height=6,width=9,dpi=800)


	# Gerreidae
		gbm.plot.fits(gerr2)

		# extract the fitted values and plot them
			a1<-plot(gerr2,'dist2shore',return.grid=TRUE)
			a2<-plot(gerr2,'prop_brs',return.grid=TRUE)
			a3<-plot(gerr2,'prop_ldsg',return.grid=TRUE)
			a4<-plot(gerr2,'prop_medsg',return.grid=TRUE)
			a5<-plot(gerr2,'prop_hdsg',return.grid=TRUE)
			a6<-plot(gerr2,'prop_sarg',return.grid=TRUE)
			a7<-plot(gerr2,'prop_urb_r',return.grid=TRUE)
			a8<-plot(gerr2,'prop_deep',return.grid=TRUE)
			a9<-plot(gerr2,'Season',return.grid=TRUE)
			hab.colors<-c('prop_brs'='navy','prop_ldsg'='violetred4','prop_medsg'='cadetblue4','prop_hdsg'='cadetblue3','prop_sarg'='goldenrod3','prop_urb_r'='orange3','prop_deep'='hotpink2')	
			hab.labels.nodist<-(c('prop_brs'='Bare Sand','prop_ldsg'='Low Density Seagrass','prop_medsg'='Medium Density Seagrass','prop_hdsg'=' High Density Seagrass','prop_sarg'='Sargassum','prop_urb_r'='Urban & Rocky','prop_deep'=' Deep Water'))
			gerr2_fittedvalues<-ggplot()+scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits=c(0.99,1.2))+geom_line(data=a2,aes(x=prop_brs,y=y,col='prop_brs',),lwd=1)+geom_line(data=a3,aes(x=prop_ldsg,y=y,col='prop_ldsg'),lwd=1)+geom_line(data=a4,aes(x=prop_medsg,y=y,col='prop_medsg'),lwd=1)+geom_line(data=a5,aes(x=prop_hdsg,y=y,col='prop_hdsg'),lwd=1)+geom_line(data=a6,aes(x=prop_sarg,y=y,col='prop_sarg'),lwd=1)+geom_line(data=a7,aes(x=prop_urb_r,y=y,col='prop_urb_r'),lwd=1)	+geom_line(data=a8,aes(x=prop_deep,y=y,col='prop_deep'),lwd=1)+scale_color_manual(values=hab.colors,labels=hab.labels.nodist)+labs(x='Proportion',y='Fitted value',color='Habitat Type')+theme_bw()
			ggsave(gerr2_fittedvalues,file='gerr2_fittedvalues_proportionHabTypes_BRT_feb23.png',device='png',unit='in',height=6,width=9,dpi=800)


## Interactions - investgate if pairwise interactiosn exits in the data as modelled by the BRT. And then you can plot them. 
	# SW Species 
		find.int.swspp<-gbm.interactions(sw_species2) # run the function, calculates
		find.int.swspp$interactions # print the interacton matrix
		find.int.swspp$rank.list # prints just the list of non-zero interactions. Make these into flextables. 
		# make a table from it 
		swpp.interactions<-find.int.swspp$rank.list%>%dplyr::select(var1.names,var2.names,int.size)%>%rename('Variable 1'='var1.names','Variable 2' = 'var2.names','Size of Interaction'=int.size)%>%flextable()%>%theme_alafoli()%>%align(align = 'center', part = 'all')%>%font(fontname = 'Times', part = 'all')%>%fontsize(size = 11, part = 'all')%>%color(color='black',part='all')%>%autofit()
			save_as_image(swpp.interactions,'interactions_above05_SWsppBRT_feb23.png',webshot='webshot')
	 
		swspp_intx_distandmedsg<-gbm.perspec(sw_species2,9,4) # plot pairwise interactions one at a time (or rather 2 covariates at a time with the response). Numbers correspond to the variables (easy way to identify them is the rank.list call from above). This also provides the maxmum fitted value by this interaction

	# SW Families 
		find.int.swf<-gbm.interactions(sw_families2) # run the function, calculates
		find.int.swf$interactions # print the interacton matrix
		find.int.swf$rank.list # prints just the list of non-zero interactions. Make these into flextables. 
		# make a table from it 
		swf.interactions<-find.int.swf$rank.list%>%dplyr::select(var1.names,var2.names,int.size)%>%rename('Variable 1'='var1.names','Variable 2' = 'var2.names','Size of Interaction'=int.size)%>%flextable()%>%theme_alafoli()%>%align(align = 'center', part = 'all')%>%font(fontname = 'Times', part = 'all')%>%fontsize(size = 11, part = 'all')%>%color(color='black',part='all')%>%autofit()
			save_as_image(swf.interactions,'figures+tables/BRT_sw_families_feb23/interactions_above05_SWfamiliesBRT_feb23.png',webshot='webshot')
	 
		swspp_intx_distandmedsg<-gbm.perspec(sw_families2,9,4)

	# Gerreidae 
		find.int.gerr<-gbm.interactions(gerr2) # run the function, calculates
		find.int.gerr$interactions # print the interacton matrix
		find.int.gerr$rank.list # prints just the list of non-zero interactions. Make these into flextables. 
		# make a table from it 
		gerr.interactions<-find.int.gerr$rank.list%>%dplyr::select(var1.names,var2.names,int.size)%>%rename('Variable 1'='var1.names','Variable 2' = 'var2.names','Size of Interaction'=int.size)%>%flextable()%>%theme_alafoli()%>%align(align = 'center', part = 'all')%>%font(fontname = 'Times', part = 'all')%>%fontsize(size = 11, part = 'all')%>%color(color='black',part='all')%>%autofit()
			save_as_image(swf.interactions,'figures+tables/BRT_gerreidae_feb23/interactions_above05_GerreidaeBRT_feb23.png',webshot='webshot')
	 

## PREDICT FROM MODELS
#### Predicting to new data. 

# dataframe to predict into = data.frame version of hab.noland. Keep the grid codes (jcode), remove geometry and check that jcode is a factor. Then duplicate the dataset with half as W and half as D (Season)
# can use df4preds and sf4preds wth every model. 

hab.noland<-st_as_sf(st_read('hexagon_grid_study_site_with_habitatFromSOSF_forchp3BRUVS.shp'),crs='WGS84')%>%dplyr::select(-grid_ar,-nearest_mg,-dist_cmg)

# How I made the dataframe for predicting: 
	df4preds<-as.data.frame(hab.noland)%>%dplyr::select(-geometry)%>%dplyr::select(-geometry)%>%slice(rep(1:n(),each=2))%>%mutate(Season=case_when(row_number()%%2==1 ~'D',.default='W'),.after=jcode)%>%dplyr::select(-nearest_mg,-dist_cmg)
# Import this at start: 
	df4preds<-read.csv('df_for_prediction_chp3BRUVS_habitat_noland.csv')%>%dplyr::select(-X)
		df4preds$jcode<-as.factor(df4preds$jcode)
		df4preds$Season<-as.factor(df4preds$Season)
# How I made the spatial file for matching to df4preds with predictions: 
	sf4preds<-hab.noland%>%slice(rep(1:n(),each=2))%>%mutate(Season=case_when(row_number()%%2==1 ~'D',.default='W'),.after=jcode)

# Import updated each time: 
	sf4preds<-st_as_sf(st_read('sf_for_predictions_fromBRTs_feb23.shp'),crs='WGS84')%>%dplyr::select(-preds_SWsp,-preds_SWf)
		sf4preds$jcode<-as.factor(sf4preds$jcode)
		sf4preds$Season<-as.factor(sf4preds$Season)


# Dataframes and Simple Features are ready. Use predict.gbm() to create a list f predictions for each row in df4preds. 

	# SW Species
	preds_SWspp<-predict.gbm(sw_species2,df4preds,n.trees=sw_species2$gbm.call$best.trees,type='response')	# one item per jcode. 
	#make it a df, and bind to df4preds
	p2spp<-as.data.frame(preds_SWspp)
	preds2_SWspp<-bind_cols(df4preds,p2spp)

	# SW Families
	preds_swf<-predict.gbm(sw_families2,df4preds,n.trees=sw_families2$gbm.call$best.trees,type='response')	# one item per jcode. 
	#make it a df, and bind to df4preds
	p2f<-as.data.frame(preds_swf)
	preds2_SWf<-bind_cols(df4preds,p2f)

	# Gerreidae
	preds_gerr<-predict.gbm(gerr2,df4preds,n.trees=gerr2$gbm.call$best.trees,type='response')	# one item per jcode. 
	#make it a df, and bind to df4preds
	# back transform
	p2g<-as.data.frame(preds_gerr)
	preds2_gerr<-bind_cols(df4preds,p2g)

# match df4preds with predictions to the spatial geometry file. Remember Seasons - need to match. can't just cbind bcause sf stored data differently 

# SW Species
	# filter predictions by season (W/D)
		spppreds.dry<-preds2_SWspp%>%filter(Season!='W')
		spppreds.wet<-preds2_SWspp%>%filter(Season!='D')
	# add predictions to sf4preds, filtered by season.
		sf4predsDRY<-sf4preds%>%filter(Season=='D')%>%mutate(preds_SWsp=spppreds.dry$preds_SWspp,.before='geometry')
		sf4predsWET<-sf4preds%>%filter(Season=='W')%>%mutate(preds_SWsp=spppreds.wet$preds_SWspp,.before='geometry')
	# (re)combine separate seasons
		sf4preds<-bind_rows(sf4predsWET,sf4predsDRY)
		sf4preds$Season<-as.factor(sf4preds$Season)

# SW Family
	# filter predictions by season (W/D)
		fpreds.dry<-preds2_SWf%>%filter(Season!='W')
		fpreds.wet<-preds2_SWf%>%filter(Season!='D')
	# add predictions to sf4preds, filtered by season.
		sf4predsDRY<-sf4preds%>%filter(Season=='D')%>%mutate(preds_SWf=fpreds.dry$preds_swf,.before='geometry')
		sf4predsWET<-sf4preds%>%filter(Season=='W')%>%mutate(preds_SWf=fpreds.wet$preds_swf,.before='geometry')
	# (re)combine separate seasons
		sf4preds<-bind_rows(sf4predsWET,sf4predsDRY)
		sf4preds$Season<-as.factor(sf4preds$Season)

# Gerreidae
	# filter predictions by season (W/D)
		gerr.preds.dry<-preds2_gerr%>%filter(Season!='W')
		gerr.preds.wet<-preds2_gerr%>%filter(Season!='D')
	# add predictions to sf4preds, filtered by season.
		sf4predsDRY<-sf4preds%>%filter(Season=='D')%>%mutate(preds_Gerr=gerr.preds.dry$preds_gerr,.before='geometry')
		sf4predsWET<-sf4preds%>%filter(Season=='W')%>%mutate(preds_Gerr=gerr.preds.wet$preds_gerr,.before='geometry')
	# (re)combine separate seasons
		sf4preds<-bind_rows(sf4predsWET,sf4predsDRY)
		sf4preds$Season<-as.factor(sf4preds$Season)

	# back transform
		sf4preds<-sf4preds%>%mutate(bt.gerrPD=exp(.$preds_Gerr),.before='geometry')


	# save the updated sf4preds and df4preds
		st_write(sf4preds,'sf_for_predictions_fromBRTs_feb23.shp',append=FALSE) # append command set to FALSE replaces the old file with the new one. If set to TRUE it would add the new data as further layers.
		st_write(sf4preds,'sf_for_predictions_fromBRTs_feb23.csv',delete_layer=TRUE,delete_dsn=TRUE)

## PLOT
## some helpful formatting code for all models/plots
	season.labels<-as_labeller(c('D'='Dry','W'='Wet'))

	{sf4preds<-st_as_sf(st_read('sf_for_predictions_fromBRTs_feb23.shp'),crs='WGS84')%>%rename(btGerPD='bt_grPD',btGerSimPD='bt_g_PD')
		sf4preds$jcode<-as.factor(sf4preds$jcode)
		sf4preds$Season<-as.factor(sf4preds$Season)}
		# predictions for full models are stored in sf4preds, as preds_SWf, _SWsp, _gerr.


# SW Species
	swspp.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=preds_SWsp),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0.5,1.7),guide=guide_colourbar(title=' Species \n\ Diversity \n\ Index'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
	ggsave(swspp.plot,file='figures+tables/shannon_index_species_BRTpreds_feb23.png',device='png',units='in',height=8,width=10,dpi=1000)

# SW Families 
	swf.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=preds_SWf),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0.5,1.7),guide=guide_colourbar(title=' Family \n\ Diversity \n\ Index'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
	ggsave(swf.plot,file='figures+tables/BRT_sw_families_feb23/shannon_index_families_BRTpreds_feb23.png',device='png',units='in',height=8,width=10,dpi=1000)

# Gerreidae
	gerr.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=btGerPD),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0,50),guide=guide_colourbar(title='Predicted Abundance \n\ of Gerreidae'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
	ggsave(gerr.plot,file='figures+tables/BRT_gerreidae_feb23/abundance_Gerreidae_BRTpreds_feb23.png',device='png',units='in',height=8,width=10,dpi=1000)













