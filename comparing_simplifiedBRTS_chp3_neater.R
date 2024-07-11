#### Neat Version of BRTS for bruvs teleost data for chapter 3 
## boosted regression trees, simplified

## created by Molly Kressler, 21 February 2023

## updated 10 August 2023, new habitat data from Emily Courmier. re-running the simplified BRTs process 
## updated 17 June 2024, only S Driscoll BRUVS


## Load Workspace and Files 
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd')
pacman::p_load(tidyverse,sf,ggplot2,gridExtra,flextable,sf,ggsn,dismo,gbm,patchwork)

# shapefiles, useful for plotting
	land<-st_as_sf(st_union(st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')))
	
	# data post August 2023 - use this
	hab.grid<-st_as_sf(st_read('winter2020habitat_hexagon_grid_ECdata.shp'),crs='WGS84')%>%mutate(jcode=as_factor(jcode))
	hab.noland<-st_as_sf(st_read('winter2020habitat_hexagon_grid_noLAND_ECdata.shp'),crs='WGS84')%>%
		mutate(jcode=as_factor(jcode))
	summary(hab.noland)
		# check
			ggplot()+geom_sf(data=hab.noland,alpha=0.5,col='violetred2')+
				theme_bw()

		
# Data frames and shape files with data 
	## need a df and a sf for predicting into: jcodes and habitat data, and Season - W and D
		df4preds<-read.csv('winter2020habitat_hexagon_grid_NOland_ECdata.csv')
			
			summary(df4preds)

		sf4preds<-st_as_sf(st_read('winter2020habitat_hexagon_grid_NOland_ECdata.shp'),crs='WGS84')
			summary(sf4preds)
			# predictions for full models are stored in sf4preds, as preds_SWf, _SWsp, _gerr.

		joined_df<-read.csv('bruvs_data_SarahDriscollonly_joinedWITHhabitat_winter2020EChabitat_june24.csv')%>%
			rename(Gerreidae = Gerreid, dist2shore = dst2shr)
			joined_df$BRUV<-as.factor(joined_df$BRUV)	
		summary(joined_df)		


# Model RDS - BRTS, with all 9 predictor variables, 'full'
	richfull <- readRDS('resource_chp3/model_RDS/Species_Rchness_Poisson_brt_tc5_lr0001_poisson_NoSeason_jun24.RDS')
	richfull_tc3 <- readRDS('resource_chp3/model_RDS/Species_Rchness_Poisson_brt_tc3_lr0001_poisson_NoSeason_jun24.RDS')
	gerrfull <- readRDS('resource_chp3/model_RDS/gerrPoisson_brt_tc3_lr0001_poisson_NoSeason_jun24.RDS')


	# deprecated june 2024
	sppfull<-readRDS('resource_chp3/model_RDS/SW_Species_brt_tc3_lr001_gaussian_aug23.RDS')
	famfull<-readRDS('resource_chp3/model_RDS/SW_Families_brt_tc9_lr001_gaussian_aug23.RDS')
	gerrfull<-readRDS('resource_chp3/model_RDS/gerrPoisson_brt_tc3_lr0001_poisson_aug23.RDS')

########################

## Instructions

# Run through the following code for each model. In the final step, saving the end outputs (an updated sf4preds and a plot of the predicted values), the file name needs to be changed to reference the appropriate model. 

########################


# Simplify - 

fullmodel<-gerrfull

simple<-gbm.simplify(fullmodel,n.drops=2) # 4 predictors, must have at least 2, so n.drops = 2

simple.model<-gbm.step(joined_df,gbm.x=simple$pred.list[[2]],gbm.y=fullmodel$gbm.call$gbm.y,tree.complexity=fullmodel$gbm.call$tree.complexity,learning.rate=as.numeric(fullmodel$gbm.call$learning.rate),bag.fraction=fullmodel$gbm.call$bag.fraction,family=fullmodel$gbm.call$family,plot.main = TRUE) # 1.8, 3950 trees - selectign this one
simple1.model<-gbm.step(joined_df,gbm.x=simple$pred.list[[1]],gbm.y=fullmodel$gbm.call$gbm.y,tree.complexity=fullmodel$gbm.call$tree.complexity,learning.rate=as.numeric(fullmodel$gbm.call$learning.rate),bag.fraction=fullmodel$gbm.call$bag.fraction,family=fullmodel$gbm.call$family,plot.main = TRUE) # 1.7, 5900 trees 

simple.model$var.names
summary(simple1.model)

saveRDS(simple1.model,'resource_chp3/model_RDS/simplified_Gerries_tc3_BRT_jun24_lowdensg_removed.RDS')

	# sppsimple<-simple.model # high density and Season removed
	# gerrsimple<-simple1.model # lowdensityseagrass removed


# Get n.trees and mean deviance of original & simplified models 
	richsimple<-readRDS('resource_chp3/model_RDS/simplified_sppRichness_tc5_BRT_jun24_lowandhigh_densg_removed.RDS')
	gerrsimple<-readRDS('resource_chp3/model_RDS/simplified_Gerries_tc3_BRT_jun24_lowdensg_removed.RDS')

	info<-as.data.frame(matrix(ncol=3,nrow=4))
	colnames(info)=c('Model','Deviance','n.trees')
	info$Model<-c('gerrfull', 'gerrsimple', 'richfull', 'richsimple')
	for(i in info$Model){
		mod<-get(i)
		dev<-as.numeric(mod$cv.statistics$deviance.mean[1])
		info$Deviance[info$Model==i]<-dev
		trees<-as.numeric(mod$n.trees)
		info$n.trees[info$Model==i]<-trees
			}

	info2<-info%>%
		mutate('Explanatory variables'=c('low, medium and high density seagrass & dist. to shore','medium and high density seagrass & dist. to shore','low, medium and high density seagrass & dist. to shore','medium density seagrass & dist. to shore'))%>%
			mutate_if(is.numeric, round, digits = 3)%>%
			flextable()%>%
			theme_zebra()%>%
			autofit()

	save_as_image(info2,'resource_chp3/BRTS_outputs/ntrees_and_deviance_of_full_and_simpleBRTS_chp3_jun24.png',webshot='webshot')
	save_as_docx(info2, path='resource_chp3/BRTS_outputs/ntrees_and_deviance_of_full_and_simpleBRTS_chp3_jun24.docx')


	###########
	## without family 

		sppsimple<-readRDS('resource_chp3/model_RDS/simplified_speciesShannonIndex_BRT_aug23_highDenSgANDseasonremoved.RDS')
		gerrsimple<-readRDS('resource_chp3/model_RDS/simplified_GerriesShannonIndex_BRT_aug23_highDenSgANDseasonremoved.RDS')
	    sppfull<-readRDS('resource_chp3/model_RDS/SW_Species_brt_tc3_lr001_gaussian_aug23.RDS')
	    gerrfull<-readRDS('resource_chp3/model_RDS/gerrPoisson_brt_tc3_lr0001_poisson_aug23.RDS')


		info<-as.data.frame(matrix(ncol=4,nrow=4))
		colnames(info)=c('model','modelnames','deviance','n.trees')
		info$modelnames<-c('Species Full','Gerridae Full','Species Simplified','Gerridae Simplified')
		info$model<-c('sppfull','gerrfull','sppsimple','gerrsimple')
		for(i in info$model){
			mod<-get(i)
			dev<-as.numeric(mod$cv.statistics$deviance.mean[1])
			info$deviance[info$model==i]<-dev
			trees<-as.numeric(mod$n.trees)
			info$n.trees[info$model==i]<-trees
				}
		
		info2<-info%>%
			dplyr::select(-model)%>%
			rename('Model' = modelnames)%>%
			mutate('Parameters Included'=c('low, medium, & high density seagrasses, season (W/D), distance to shore','low, medium, & high density seagrasses, season (W/D), distance to shore','low & medium density seagrasses, distance to shore','low & medium density seagrasses,  distance to shore'))%>%
			mutate_if(is.numeric, round, digits = 3)%>%
			flextable()%>%
			set_header_labels(model='Model', deviance = 'Deviance')%>%
			theme_zebra()%>%
			align(align = 'center', part = 'all')%>%
			font(fontname = 'Arial', part = 'all')%>%
			fontsize(size = 11, part = 'all')%>%
			color(color='black',part='all')%>%
			autofit()

		save_as_image(info2,'resource_chp3/BRTS_outputs/ntrees_and_deviance_of_full_and_simpleBRTS_chp3.png',webshot='webshot',res=650)



# Predict 
	gerr.predicts<-predict.gbm(gerrsimple,sf4preds,n.trees=gerrsimple$gbm.call$best.trees,type='response')
	gerr.predicts<-as.data.frame(gerr.predicts)

	rich.predicts<-predict.gbm(richsimple,sf4preds,n.trees=richsimple$gbm.call$best.trees,type='response')
	rich.predicts<-as.data.frame(rich.predicts)
	
	predicts.sf<-bind_cols(sf4preds,rich.predicts, gerr.predicts)

	# save the updated sf4preds and df4preds
		st_write(predicts.sf,'sf_for_predictions_fromBRTs_jun24.shp',driver='ESRI Shapefile',append=FALSE) # append command set to FALSE replaces the old file with the new one. If set to TRUE it would add the new data as further layers.
		st_write(predicts.sf,'sf_for_predictions_fromBRTs_jun24.csv', driver = 'CSV')


# Plot 
	predsf<-st_as_sf(st_read('sf_for_predictions_fromBRTs_jun24.shp'),crs='WGS84')
	predsf$jcode<-as.factor(sf4preds$jcode)

	summary(predsf) # rich 6.7, 8.5. gerr 0.3, 1.4


	spp.simple.plot<-ggplot()+geom_sf(data=predsf,aes(fill=rch_prd),col=NA)+
		scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(6,9),guide=guide_colourbar(title=' Species \n\ Richness'))+
		geom_sf(data=land,fill='gray98')+
		theme_bw()

	gerr.simple.plot<-ggplot()+
		geom_sf(data=predsf,aes(fill=grr_prd),col=NA)+
		scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0,3),guide=guide_colourbar(title='Gerreidae\n\ Abundance\n\ (counts)'))+
		geom_sf(data=land,fill='gray98')+
		theme_bw()+
		theme(axis.text.y = element_blank())

	hz <- spp.simple.plot + gerr.simple.plot

	ggsave(gerr.simple.plot,file='resource_chp3/BRTS_outputs/simplified_BRT_Gerreidae_preds_jun24.png',device='png',units='in',height=8,width=5,dpi=400)

	ggsave(spp.simple.plot,file='resource_chp3/BRTS_outputs/simplified_BRT_SpeciesRichness_preds_jun24.png',device='png',units='in',height=8,width=5,dpi=400)
	
	ggsave(hz,file='resource_chp3/BRTS_outputs/simplified_BRTs_richness_and_gerries_preds_jun24_forpub.png',device='png',units='in',height=8,width=10,dpi=400)


	## deprecated june 2024
		spp.simple.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=sppsimpPD),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0.5,1.7),guide=guide_colourbar(title=' Species \n\ Diversity \n\ Index'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
		fam.simple.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=famsimpPD),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0.1,1.4),guide=guide_colourbar(title=' Family \n\ Diversity \n\ Index'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
		gerr.simple.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=btGerSmPD),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(3,200),guide=guide_colourbar(title='Gerreidae\n\ Abundance\n\ (counts)'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
				sppsimple$contributions

		ggsave(gerr.simple.plot,file='resource_chp3/BRTS_outputs/simplified_BRT_Gerreidaebacktransformed_preds_aug23.png',device='png',units='in',height=8,width=10,dpi=1000)
		ggsave(spp.simple.plot,file='resource_chp3/BRTS_outputs/simplified_BRT_SpeciesSWbacktransformed_preds_aug23.png',device='png',units='in',height=8,width=10,dpi=1000)
		ggsave(fam.simple.plot,file='resource_chp3/BRTS_outputs/simplified_BRT_FamilySWbacktransformed_preds_aug23.png',device='png',units='in',height=8,width=10,dpi=1000)

# compare the predicted full vs predicted simple versus real distirbutions with boxplots (predicted) and points from raw data
# real = joined_df
# predicted = sf4preds

# dataframe with one col for each models predictions
	pred<-as_tibble(predsf)%>%dplyr::select(jcode,prp_mds, prp_hds, dst2shr, rch_prd, grr_prd) 
	summary(pred)
	pred.long<- pred%>% pivot_longer(pred,cols=rch_prd:grr_prd,names_to='model',values_to='preds')%>%
		mutate(model = as.factor(model))
	summary(pred.long)

# dataframe with original data 
	raw<-joined_df%>%dplyr::select(jcode,prp_mds, prp_hds, dist2shore, spp_rch, Gerreidae)
	raw.long<- raw %>% pivot_longer(cols=spp_rch:Gerreidae,names_to='model',values_to='preds') %>%
		mutate(jcode = as.factor(jcode))

# join 
	data <- left_join(pred.long, raw.long, by = 'jcode', relationship = 'many-to-many')


# plot 
	season.labels<-as_labeller(c('D'='Dry','W'='Wet'))
	indices.labels<-as_labeller(c('sppsimpPD'='SW Families (simplified)',
		'sppflPD'='SW Species (full)'))#,
		#'famsimpPD'='SW Families (simplified)',
		#'sppsimpPD'='SW Species (simplified)'))
	gerr.labels<-as_labeller(c('bt_g_PD'='Gerreidae (full)','
		btGerSmPD'='Gerreidae (simplified)'))

	famspp_preds<-pred.long%>%dplyr::filter(model!=c('bt_g_PD','btGerSmPD'))			
	famspp_raw<-raw.long2%>%dplyr::filter(model!=c('bt_g_PD'))%>%dplyr::filter(model!=c('btGerSmPD'))
	summary(famspp_preds)

	indices.plot<-ggplot(pred.long%>%dplyr::filter(model!=c('bt_g_PD','btGerSmPD')),aes(x=model,y=preds))+
		geom_boxplot()+
		geom_jitter(data=famspp_raw,aes(x=model,y=preds),pch=19,alpha=0.7,col='#3a6c74')+
		ylab(lab='Diversity Index Value')+
		xlab(lab='Model')+
		scale_x_discrete(labels=indices.labels)+
		theme_bw()+
		theme(axis.text.x = element_text(angle = 30,vjust=0.5,size=10))

	gfullpredlong <- pred.long %>% filter(model=='bt_g_PD')
	gsimppredlong <- pred.long %>% filter(model=='btGerSimPD')

	gerr.plot<-	ggplot(pred.long%>%dplyr::filter(model==c('bt_g_PD','btGerSmPD')),aes(x=model,y=preds))+
		geom_boxplot()+
		geom_jitter(data=raw.long2%>%dplyr::filter(model==c('bt_g_PD','btGerSmPD')),aes(x=model,y=preds),pch=19,alpha=0.7,col='#3a6c74')+
		ylab(lab='Values (count)')+xlab(lab='Model')+
		scale_x_discrete(labels=c('Gerreidae (full)','Gerreidae (simplified)'))+
		theme_bw()+ 
		theme(axis.text.x = element_text(angle = 30,vjust=0.5,size=10))

	comparison <- indices.plot + gerr.plot + plot_annotation( title = 'Predicted Values versus Observations, for Full and Simplified Models')

	ggsave(comparison,file='resource_chp3/BRTS_outputs/Predicted Values versus Observations, for Full and Simplified Models_AUGUST2023.png',device='png',units='in',height=6,width=8.5,dpi=850)

	# Now - not moving forward with SW Families. So remake plots with Spp Index and re-combien with gridExtra with Gerreidae

	spponly.plot<-ggplot(famspp_preds%>%filter(model!='famsimpPD')%>%filter(model!='fmfllPD'),aes(x=model,y=preds))+geom_boxplot()+
		geom_jitter(data=famspp_raw%>%filter(model!='famsimpPD')%>%filter(model!='fmfllPD'),aes(x=model,y=preds),pch=19,alpha=0.7,col='#3a6c74')+
		ylab(lab='Diversity Index Value')+
		xlab(lab='Model')+
		scale_x_discrete(labels=indices.labels)+
		theme_bw()+ 
		theme(axis.text.x = element_text(angle = 30,vjust=0.5,size=10))

	comparison2 <- spponly.plot + gerr.plot + plot_annotation( title = 'Predicted Values versus Observations, for Full and Simplified Models')

		ggsave(comparison2,file='resource_chp3/BRTS_outputs/Predicted Values versus Observations, for Full and Simplified Models, no Family Diversity Index_AUGUST2023.png',device='png',units='in',height=6,width=12,dpi=1000)
	
	# only simplified models

	spponly.simpleonly.plot<-ggplot(famspp_preds%>%filter(model=='sppsimpPD')%>%filter(model!='fmfllPD'),aes(x=model,y=preds))+
		geom_boxplot()+geom_jitter(data=famspp_raw%>%dplyr::filter(model=='sppsimpPD'),aes(x=model,y=preds),pch=19,alpha=0.7,col='orange3')+
		ylab(lab='Diversity Index Value')+
		xlab(lab=NULL)+
		scale_x_discrete(labels=indices.labels)+
		theme_bw()

	gerr.simpleonly.plot<-ggplot(pred.long%>%dplyr::filter(model=='btGerSmPD'),aes(x=model,y=preds))+
		geom_boxplot()+
		geom_jitter(data=raw.long2%>%dplyr::filter(model=='btGerSmPD'),aes(x=model,y=preds),pch=19,alpha=0.7,col='orange3')+
		ylab(lab='Values (count)')+
		xlab(lab=NULL)+
		scale_x_discrete(labels=gerr.labels)+
		theme_bw()

		comparison3<-grid.arrange(spponly.simpleonly.plot,gerr.simpleonly.plot,ncol=2,nrow=1,top=textGrob('Predicted Values versus Observations for Simplified Models',gp=gpar(fontsize=12,font=2)))

		ggsave(comparison3,file='BRTS_outputs/Predicted Values versus Observations, Simplified models of Species Diveristy and Gerreidae counts.png',device='png',units='in',height=6,width=7,dpi=1000)



#### Save predictions in a format for modelling with further in chapter 3 - PREY LANDSCAPE MODEL outputs. 

	data <- st_as_sf(st_read('resource_chp3/sf_for_predictions_fromBRTs_jun24.shp'), crs='WGS84')

# need them in the right file. 
st_write(data,'resource_chp3/predictions_from_simplifiedBRTs_Gerreidae_and_sppRichness_june24.shp')
# as a dataframe
st_write(data,'resource_chp3/predictions_from_simplifiedBRTs_Gerreidae_and_spprRichness_june24.csv')


########################

## Evaluate simplified models - Gerreidae and SW Species

########################

### import data 

	richsimple<-readRDS('resource_chp3/model_RDS/simplified_sppRichness_tc5_BRT_jun24_lowandhigh_densg_removed.RDS')
	gerrsimple<-readRDS('resource_chp3/model_RDS/simplified_Gerries_tc3_BRT_jun24_lowdensg_removed.RDS')

  	richfull <- readRDS('resource_chp3/model_RDS/Species_Rchness_Poisson_brt_tc5_lr0001_poisson_NoSeason_jun24.RDS')
	gerrfull <- readRDS('resource_chp3/model_RDS/gerrPoisson_brt_tc3_lr0001_poisson_NoSeason_jun24.RDS')

	spp.hab.labels<-(c('dist2shore'='Dist. to Shore (m)', 'prp_lds'='Prop. of \n\ Low Density \n\ Seagrass','prp_mds'='Prop. of  \n\ Medium Density \n\ Seagrass','prp_hds'='Prop. of  \n\ High Density \n\ Seagrass'))
	gerr.hab.labels<-(c('prp_lds'='Prop. of \n\ Low Density \n\ Seagrass','dist2shore'='Dist. to Shore (m)','prp_mds'='Prop. of  \n\ Medium Density \n\ Seagrass','prp_hds'='Prop. of  \n\ High Density \n\ Seagrass'))

# Influence of vars 

	infl.gerr<-gerrsimple$contributions
	infl.rich<-richsimple$contributions

	rich.relinf<-infl.rich%>%mutate(var=fct_reorder(var,rel.inf))%>%
		ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+
		geom_bar(stat='identity')+
		scale_fill_gradient(low='#3cbcfc',high='#3a6c74', space='Lab',guide='none', aesthetics='fill')+
		theme_bw()+
		coord_flip()+
		scale_x_discrete(labels=spp.hab.labels)+
		ylab('Relative Influence')+
		ggtitle('Species Diversity Simplified BRT ')+
		xlab(NULL)

		ggsave(rich.relinf,file='resource_chp3/BRTS_outputs/relative_influence_vars_in_simpifiedsppRichnessBRT_jun24.png',device='png',units='in',height=4,width=6,dpi=450)
	
	gerr.relinf<-infl.gerr%>%mutate(var=fct_reorder(var,rel.inf))%>%ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+
		geom_bar(stat='identity')+
		scale_fill_gradient(low='#3cbcfc',high='#3a6c74', space='Lab',guide='none', aesthetics='fill')+
		theme_bw()+
		coord_flip()+
		ggtitle('Gerreidae Simplified BRT ')+
		scale_x_discrete(labels=gerr.hab.labels)+
		ylab('Relative Influence')+xlab(NULL)
	gerr.relinf
		ggsave(gerr.relinf,file='resource_chp3/BRTS_outputs/relative_influence_vars_in_simplifiedGerreidaeBRT_june24.png',device='png',units='in',height=4,width=6,dpi=450)
 
# Fitted values for dist2shore

	preds <- read.csv('resource_chp3/predictions_from_simplifiedBRTs_Gerreidae_and_spprRichness_june24.csv')

	dist2shorefitted_values_simplifiedmodels<-ggplot()+
			geom_smooth(data=preds,aes(x=dist2shore,y=sppsimpPD),fill='3a6c74',col='grey26',alpha=0.5)+
			geom_smooth(data=preds,aes(x=dist2shore,y=btGerSmPD),fill='3cbcfc',alpha=0.5,col='grey26')+
			xlab('Distance to Shore (m)')+
			ylab('Fitted Values')+
			scale_x_continuous(limits=c(0,1800))+
			scale_y_continuous()+
			theme_bw()
	dist2shorefitted_values_simplifiedmodels
		ggsave(dist2shorefitted_values_simplifiedmodels,file='resource_chp3/BRTS_outputs/smoothed_fittedvalues_simplifiedSWSpeciesINorange_and_SimplifiedDist2ShoreINpurple_fittedvalues_BRT_feb23.png',device='png',unit='in',height=6,width=9,dpi=800)


	# compare ldsg and dist2shore in gerrsimple
	gbm.perspec(sppsimple,1,2)

	# compare 2 vars effects in sppsimple
	sppsimple.gbmint<-gbm.interactions(sppsimple)
	sppsimple.gbmint$rank.list
	gbm.perspec(sppsimple,3,2)

# Diagnostic plots - QQs 

preds <- read.csv('resource_chp3/predictions_from_simplifiedBRTs_Gerreidae_and_spprRichness_june24.csv')

gerr.simple.qq <- ggplot(data=preds, aes(sample = gerr.predicts))+
	stat_qq(size=1,pch=21)+
	labs(title='Simplified Gerridae BRT', subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
	stat_qq_line(linetype=2, col='red')+
	theme_bw()

	gerr.simple.qq

	ggsave(gerr.simple.qq,file = 'resource_chp3/BRTS_outputs/quantilequantile_plot_gerrsimple_jun24.png', dpi=450, unit= 'in', height = 4, width = 4)

richness.qq <- ggplot(data=preds, aes(sample = rich.predicts))+
	stat_qq(size=1,pch=21)+
	labs(title='Simplified Species Richness BRT', subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
	stat_qq_line(linetype=2, col='red')+
	theme_bw()
	richness.qq
	ggsave(richness.qq,file = 'resource_chp3/BRTS_outputs/quantilequantile_plot_sppRichness_jun24.png', dpi=450, unit= 'in', height = 4, width = 4)


































