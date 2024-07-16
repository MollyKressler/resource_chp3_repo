#### Neat Version of BRTS for bruvs teleost data for chapter 3 
## boosted regression trees, simplified

## created by Molly Kressler, 21 February 2023

## updated 10 August 2023, new habitat data from Emily Courmier. re-running the simplified BRTs process 
## updated 17 June 2024, only S Driscoll BRUVS


## Load Workspace and Files 
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd')
pacman::p_load(tidyverse,sf,ggplot2,gridExtra,flextable,sf,ggsn,dismo,gbm,patchwork)


#####
## - Load data 
#####

	# shapefiles, useful for plotting
		land<-st_as_sf(st_union(st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')))
		
		# data post August 2023 - use this
		hab.grid<-st_as_sf(st_read('winter2020habitat_hexagon_grid_ECdata.shp'),crs='WGS84')%>%mutate(jcode=as_factor(jcode))
		hab.noland<-st_as_sf(st_read('winter2020habitat_hexagon_grid_noLAND_ECdata.shp'),crs='WGS84')%>%
			mutate(jcode=as_factor(jcode))

			
	# Data frames and shape files with data 
		## need a df and a sf for predicting into: jcodes and habitat data, and Season - W and D
		sf4preds<-st_as_sf(st_read('winter2020habitat_hexagon_grid_NOland_ECdata.shp'),crs='WGS84')
		b <-read.csv('bruvs2014AND2018_data_joinedWITHhabitat_summer18habANDwinter2014hab_july2024.csv')


	# Model RDS - BRTS, with all 9 predictor variables, 'full'
		n1 <- readRDS('resource_chp3/model_RDS/maxN_gbm_poisson_july2024.RDS')

		# deprecated june 2024 & july 2024
			sppfull<-readRDS('resource_chp3/model_RDS/SW_Species_brt_tc3_lr001_gaussian_aug23.RDS')
			famfull<-readRDS('resource_chp3/model_RDS/SW_Families_brt_tc9_lr001_gaussian_aug23.RDS')
			gerrfull<-readRDS('resource_chp3/model_RDS/gerrPoisson_brt_tc3_lr0001_poisson_aug23.RDS')		
			richfull <- readRDS('resource_chp3/model_RDS/Species_Rchness_Poisson_brt_tc5_lr0001_poisson_NoSeason_jun24.RDS')
			richfull_tc3 <- readRDS('resource_chp3/model_RDS/Species_Rchness_Poisson_brt_tc3_lr0001_poisson_NoSeason_jun24.RDS')
			gerrfull <- readRDS('resource_chp3/model_RDS/gerrPoisson_brt_tc3_lr0001_poisson_NoSeason_jun24.RDS')


########################

## Instructions

# Run through the following code for each model. In the final step, saving the end outputs (an updated sf4preds and a plot of the predicted values), the file name needs to be changed to reference the appropriate model. 

########################

	# Simplify - 

	fullmodel<-n1

	simple<-gbm.simplify(fullmodel,n.drops=2) # 4 predictors, must have at least 2, so n.drops = 2

	simple.model<-gbm.step(b,gbm.x=simple$pred.list[[2]],gbm.y=fullmodel$gbm.call$gbm.y,tree.complexity=fullmodel$gbm.call$tree.complexity,learning.rate=as.numeric(fullmodel$gbm.call$learning.rate),bag.fraction=fullmodel$gbm.call$bag.fraction,family=fullmodel$gbm.call$family,plot.main = TRUE) # 1.8, 3950 trees - selectign this one
	simple1.model<-gbm.step(b,gbm.x=simple$pred.list[[1]],gbm.y=fullmodel$gbm.call$gbm.y,tree.complexity=fullmodel$gbm.call$tree.complexity,learning.rate=as.numeric(fullmodel$gbm.call$learning.rate),bag.fraction=fullmodel$gbm.call$bag.fraction,family=fullmodel$gbm.call$family,plot.main = TRUE) # 1.7, 5900 trees 

	simple.model$var.names # marginally better bc more trees, and marginal change in the standard residual deviance 
	summary(simple.model)

	saveRDS(simple.model,'resource_chp3/model_RDS/simplified_maxN_gbm_poisson_HDS_dist2shore_july2024.RDS')

	# Get n.trees and mean deviance of original & simplified models 
		mod <- list(n1, simple.model)
		
		info<-as.data.frame(matrix(ncol=3,nrow=2))
		colnames(info)=c('Model','Deviance','n.trees')
		info$Model<-c('n1', 'simple.model')
		for(i in info$Model){
			mod<-get(i)
			dev<-as.numeric(mod$cv.statistics$deviance.mean[1])
			info$Deviance[info$Model==i]<-dev
			trees<-as.numeric(mod$n.trees)
			info$n.trees[info$Model==i]<-trees
				}

		info2<-info%>%
			mutate('Explanatory variables'=c('low, medium and high density seagrass & dist. to shore','high density seagrass & dist. to shore'))%>%
				mutate_if(is.numeric, round, digits = 3)%>%
				flextable()%>%
				theme_zebra()%>%
				autofit()

		save_as_image(info2,'resource_chp3/BRTS_outputs/BRT_maxN_july2024/ntrees_and_deviance_of_full_and_simpleBRTS_chp3_july2024.png',webshot='webshot')
		save_as_docx(info2, path='resource_chp3/BRTS_outputs/BRT_maxN_july2024/ntrees_and_deviance_of_full_and_simpleBRTS_chp3_july2024.docx')


		###########
		## deprecated - without family 

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


######
## - Evaluate simplified model 
######
		simple.model <- readRDS('resource_chp3/model_RDS/simplified_maxN_gbm_poisson_HDS_dist2shore_july2024.RDS')
		hab.labels<-(c('dist2shore'='Dist. to Shore (m)', 'prp_lds'='Prop. of \n\ Low Density \n\ Seagrass','prp_mds'='Prop. of  \n\ Medium Density \n\ Seagrass','prp_hds'='Prop. of  \n\ High Density \n\ Seagrass'))

		infl.simple<-simple.model$contributions
		simple.relinf<-infl.simple%>%
			mutate(var=fct_reorder(var,rel.inf))%>%
			ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+
				geom_bar(stat='identity')+
				scale_fill_distiller(direction=1,palette='Blues',limits=c(35,65),guide='none')+
				theme_bw()+
				coord_flip()+
				scale_x_discrete(labels=hab.labels)+
				ylab('Relative Influence')+
				xlab(NULL)	+
				theme(text = element_text(size = 14))
		ggsave(simple.relinf,file='resource_chp3/BRTS_outputs/BRT_maxN_july2024/relative_influence_vars_maxN_BRT_SIMPLIFIED_july2024.png',device='png',units='in',height=4,width=5.5,dpi=900)
		
		res_simple <- as_tibble(resid(simple.model))%>%rename(resids = value)
		F1 <- as_tibble(predict(simple.model))%>%rename(fitted = value)
		diag_simple <- bind_cols(F1, res_simple)

		resVfit_simple <- ggplot(data = diag_simple,aes(x= fitted, y = resids))+
			geom_point(col = '#3A6C74', fill = '#3A6C74')+ 
			labs(subtitle = 'Residuals vs. Fitted', y = 'Residuals', x = 'Fitted')+
			theme_bw()
		res_simple_hist <- ggplot(data = res_simple, aes(x = resids))+ 
			geom_histogram(binwidth = nrow(res_n1)/100,col = '#3A6C74', fill = '#3A6C74')+ 
			theme_bw()+
			labs(subtitle = 'Residuals', y = 'Frequency', x = 'Residuals')
		res_simple_qq <- ggplot(data=F1, aes(sample = fitted))+
			stat_qq(size=1,pch=21,col = '#3A6C74', fill = '#3A6C74')+
			labs(subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
			stat_qq_line(linetype=2, col='red')+
			theme_bw()

		diagnostics_simple <- resVfit_simple+res_simple_hist+res_simple_qq
		
		ggsave(diagnostics_simple, file = 'resource_chp3/BRTS_outputs/BRT_maxN_july2024/diagnostics_BRT_SIMPLIFIED_maxN_july2024.png', device = 'png', unit = 'in', height = 4, width = 8, dpi = 850)	

######
## - Predict into 2020 habitat data - hexagon grid
######

	maxN.preds<-predict.gbm(simple.model,sf4preds,n.trees=simple.model$gbm.call$best.trees,type='response')
	maxN.preds<-as.data.frame(maxN.preds) 

	predicts.sf<-bind_cols(sf4preds, maxN.preds)

	#####
	## - save
	#####
		st_write(predicts.sf,'predictions_from_simplifiedBRTs_Gerreidae_and_sppRichness_june24.shp',driver='ESRI Shapefile',append=FALSE) # append command set to FALSE replaces the old file with the new one. If set to TRUE it would add the new data as further layers.
		st_write(predicts.sf,'predictions_from_simplifiedBRTs_Gerreidae_and_sppRichness_june24.csv', driver = 'CSV')

	## plot

	preds <- st_as_sf(st_read('predictions_from_simplifiedBRTs_Gerreidae_and_sppRichness_june24.shp'), crs = 'WGS84')
		
	maxN.plotLOG<-ggplot()+geom_sf(data=preds,aes(fill=log(maxN_preds)),col=NA)+
		scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0,5),guide=guide_colourbar(title=' MaxN of 4\n\ Families'))+
		geom_sf(data=land,fill='gray98')+
		theme_bw()
	ggsave(maxN.plotLOG,file='resource_chp3/BRTS_outputs/simplified_BRT_LOG_maxN_preds_july24.png',device='png',units='in',height=8,width=5,dpi=850)		

	maxN.plot<-ggplot()+geom_sf(data=preds,aes(fill=maxN_preds),col=NA)+
		scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0,100),guide=guide_colourbar(title=' MaxN of 4\n\ Families'))+
		geom_sf(data=land,fill='gray98')+
		theme_bw()
	ggsave(maxN.plot,file='resource_chp3/BRTS_outputs/simplified_BRT_maxN_preds_july24.png',device='png',units='in',height=8,width=5,dpi=850)

	## deprecated june 2024
		spp.simple.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=sppsimpPD),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0.5,1.7),guide=guide_colourbar(title=' Species \n\ Diversity \n\ Index'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
		fam.simple.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=famsimpPD),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0.1,1.4),guide=guide_colourbar(title=' Family \n\ Diversity \n\ Index'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
		gerr.simple.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=btGerSmPD),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(3,200),guide=guide_colourbar(title='Gerreidae\n\ Abundance\n\ (counts)'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
				sppsimple$contributions

		ggsave(gerr.simple.plot,file='resource_chp3/BRTS_outputs/simplified_BRT_Gerreidaebacktransformed_preds_aug23.png',device='png',units='in',height=8,width=10,dpi=1000)
		ggsave(spp.simple.plot,file='resource_chp3/BRTS_outputs/simplified_BRT_SpeciesSWbacktransformed_preds_aug23.png',device='png',units='in',height=8,width=10,dpi=1000)
		ggsave(fam.simple.plot,file='resource_chp3/BRTS_outputs/simplified_BRT_FamilySWbacktransformed_preds_aug23.png',device='png',units='in',height=8,width=10,dpi=1000)


######
## - Predict into 2020 habitat data - buffers
######

	buffs.sf <- st_as_sf(st_read('resource_chp3/sf_with_predictions_INTO_BUFFERS_fromBRTs_aug23.shp'),crs='WGS84')%>%
		dplyr::select(-bt_g_PD, -smpl_PD, -sppsmPD)%>%
		rename(dist2shore = dst2shr)

	simple.model <- readRDS('resource_chp3/model_RDS/simplified_maxN_gbm_poisson_HDS_dist2shore_july2024.RDS')

	maxN.preds<-predict.gbm(simple.model,buffs.sf,n.trees=simple.model$gbm.call$best.trees,type='response')
	maxN.preds<-as.data.frame(maxN.preds)

	predicts.sf<-bind_cols(buffs.sf, maxN.preds)

	#####
	## - save
	##### 
		st_write(predicts.sf,'sf_with_predictions_INTO_BUFFERS_fromBRTs_july2024.shp',driver='ESRI Shapefile',append=FALSE) 

		st_write(predicts.sf,'sf_with_predictions_INTO_BUFFERS_fromBRTs_july2024.csv', driver = 'CSV')


######
## - Compare raw data to predictions 
######

	pred <- st_as_sf(st_read('predictions_from_simplifiedBRTs_Gerreidae_and_sppRichness_june24.shp'), crs = 'WGS84')%>% 
		dplyr::select(jcode, maxN_preds)%>%
		mutate(jcode = as.numeric(jcode))
	b <-read_csv('bruvs2014AND2018_data_joinedWITHhabitat_summer18habANDwinter2014hab_july2024.csv')%>%
		dplyr::select(jcode, maxN)

	# join 
		data <- left_join(pred, b, by = 'jcode', relationship = 'many-to-many')%>%
			filter


	# plot
	boxplot<- ggplot()+
		geom_boxplot(data = pred, aes(y = log(maxN_preds+1)), fill = '#3A6C74', alpha = 0.2)+
		geom_jitter(data = b, aes(x = 0, y = log(1+maxN)), pch = 19, alpha = 0.75,,col = '#3A6C74')+
		theme_bw()+
		theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
		labs(y = 'MaxN (log)', x = NULL)+
		ggtitle('Simplified')

	ggsave(boxplot,file='resource_chp3/BRTS_outputs/BRT_maxN_july2024/predsVSraw_maxN_simplifiedBRT.png',device='png',units='in',height=4,width=4.5,dpi=850)


















########################

## Code graveyard - deprecated figures

########################

### Evaluate simplified models - Gerreidae and SW Species
##import data 

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
































