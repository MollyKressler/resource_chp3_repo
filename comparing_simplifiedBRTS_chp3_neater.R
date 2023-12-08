#### Neat Version of BRTS for bruvs teleost data for chapter 3 
## boosted regression trees, simplified

## created by Molly Kressler, 21 February 2023

## updated 10 August 2023, new habitat data from Emily Courmier. re-running the simplified BRTs process 


## Load Workspace and Files 
setwd('/Users/mollykressler/Documents/data_phd')
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

	# deprecated: data from MS, pre-dating AUgust 2023
		setwd('/Users/mollykressler/Documents/data_phd')
		hab.grid<-st_as_sf(st_read('hexagon_grid_withUPDATEDhabitat_sept22.shp'),crs='WGS84')
	
# Data frames and shape files with data 
	## need a df and a sf for predicting into: jcodes and habitat data, and Season - W and D
		{df4preds<-read.csv('winter2020habitat_hexagon_grid_NOland_ECdata.csv')%>%
			slice(rep(1:n(),each=2))%>%
			mutate(Season=case_when(row_number()%%2==1 ~'D',.default='W'),.after=jcode)
			df4preds$jcode<-as.factor(df4preds$jcode)
			df4preds$Season<-as.factor(df4preds$Season)}
			summary(df4preds)

		{sf4preds<-st_as_sf(st_read('winter2020habitat_hexagon_grid_NOland_ECdata.shp'),crs='WGS84')%>%
			slice(rep(1:n(),each=2))%>%
			mutate(Season=case_when(row_number()%%2==1 ~'D',.default='W'),.after=jcode)
			sf4preds$jcode<-as.factor(sf4preds$jcode)
			sf4preds$Season<-as.factor(sf4preds$Season)}
			summary(sf4preds)
			# predictions for full models are stored in sf4preds, as preds_SWf, _SWsp, _gerr.

		joined_df<-read.csv('bruvs_data_joinedWITHhabitat_winter2020EChabitat_aug23.csv')
			joined_df$BRUV<-as.factor(joined_df$BRUV)
			joined_df$Season<-as.factor(joined_df$Season)
			joined_df$SW_Species<-as.numeric(joined_df$SW_Species)
			joined_df$SW_Families<-as.numeric(joined_df$SW_Families)	
		summary(joined_df)		


# Model RDS - BRTS, with all 9 predictor variables, 'full'

	sppfull<-readRDS('resource_chp3/model_RDS/SW_Species_brt_tc3_lr001_gaussian_aug23.RDS')
	famfull<-readRDS('resource_chp3/model_RDS/SW_Families_brt_tc9_lr001_gaussian_aug23.RDS')
	gerrfull<-readRDS('resource_chp3/model_RDS/gerrPoisson_brt_tc3_lr0001_poisson_aug23.RDS')

########################

## Instructions

# Run through the following code for each model. In the final step, saving the end outputs (an updated sf4preds and a plot of the predicted values), the file name needs to be changed to reference the appropriate model. 

########################


# Simplify - 

fullmodel<-gerrfull

simple<-gbm.simplify(fullmodel,n.drops=3) # 5 predictors, must have at least 2, so n.drops = 3

simple.model<-gbm.step(joined_df,gbm.x=simple$pred.list[[2]],gbm.y=fullmodel$gbm.call$gbm.y,tree.complexity=fullmodel$gbm.call$tree.complexity,learning.rate=as.numeric(fullmodel$gbm.call$learning.rate),bag.fraction=fullmodel$gbm.call$bag.fraction,family=fullmodel$gbm.call$family,plot.main = TRUE) 
simple.model$var.names
saveRDS(simple.model,'resource_chp3/model_RDS/simplified_GerriesShannonIndex_BRT_aug23_highDenSgANDseasonremoved.RDS')

	# sppsimple<-simple.model # high density and Season removed
	# famsimple<-simple.model # medium density and Season removed
	# gerrsimple<-simple.model # high density and Season removed


# Get n.trees and mean deviance of original & simplified models 

	sppsimple<-readRDS('resource_chp3/model_RDS/simplified_speciesShannonIndex_BRT_aug23_highDenSgANDseasonremoved.RDS')
	famsimple<-readRDS('resource_chp3/model_RDS/simplified_familyShannonIndex_BRT_aug23_mediumDenSgANDseasonremoved.RDS')
	gerrsimple<-readRDS('resource_chp3/model_RDS/simplified_GerriesShannonIndex_BRT_aug23_highDenSgANDseasonremoved.RDS')

	info<-as.data.frame(matrix(ncol=3,nrow=6))
	colnames(info)=c('model','deviance','n.trees')
	info$model<-c('sppfull','famfull','gerrfull','sppsimple','famsimple','gerrsimple')
	for(i in info$model){
		mod<-get(i)
		dev<-as.numeric(mod$cv.statistics$deviance.mean[1])
		info$deviance[info$model==i]<-dev
		trees<-as.numeric(mod$n.trees)
		info$n.trees[info$model==i]<-trees
			}
	
	info2<-info%>%mutate(note=c('low, medium, & high density seagrasses, season (W/D), distance to shore','low, medium, & high density seagrasses, season (W/D), distance to shore','low, medium, & high density seagrasses, season (W/D), distance to shore','low & medium density seagrasses, distance to shore','low & high density seagrasses,  distance to shore','low & medium density seagrasses,  distance to shore'))%>%flextable()%>%theme_alafoli()%>%align(align = 'center', part = 'all')%>%font(fontname = 'Times', part = 'all')%>%fontsize(size = 11, part = 'all')%>%color(color='black',part='all')%>%autofit()

	save_as_image(info2,'resource_chp3/BRTS_outputs/ntrees_and_deviance_of_full_and_simpleBRTS_chp3.png',webshot='webshot')


	###########
	## without family 

	sppsimple<-readRDS('resource_chp3/model_RDS/simplified_speciesShannonIndex_BRT_aug23_highDenSgANDseasonremoved.RDS')
	gerrsimple<-readRDS('resource_chp3/model_RDS/simplified_GerriesShannonIndex_BRT_aug23_highDenSgANDseasonremoved.RDS')
    sppfull<-readRDS('resource_chp3/model_RDS/SW_Species_brt_tc3_lr001_gaussian_aug23.RDS')
    gerrfull<-readRDS('resource_chp3/model_RDS/gerrPoisson_brt_tc3_lr0001_poisson_aug23.RDS')


	info<-as.data.frame(matrix(ncol=3,nrow=4))
	colnames(info)=c('model','deviance','n.trees')
	info$model<-c('sppfull','gerrfull','sppsimple','gerrsimple')
	for(i in info$model){
		mod<-get(i)
		dev<-as.numeric(mod$cv.statistics$deviance.mean[1])
		info$deviance[info$model==i]<-dev
		trees<-as.numeric(mod$n.trees)
		info$n.trees[info$model==i]<-trees
			}
	
	info2<-info%>%
		mutate('Parameters Included'=c('low, medium, & high density seagrasses, season (W/D), distance to shore','low, medium, & high density seagrasses, season (W/D), distance to shore','low & medium density seagrasses, distance to shore','low & medium density seagrasses,  distance to shore'))%>%
		flextable()%>%
		set_header_labels(model='Model', deviance = 'Deviance')%>%
		theme_zebra()%>%
		align(align = 'center', part = 'all')%>%
		font(fontname = 'Arial', part = 'all')%>%
		fontsize(size = 11, part = 'all')%>%
		color(color='black',part='all')%>%
		autofit()

	save_as_image(info2,'resource_chp3/BRTS_outputs/ntrees_and_deviance_of_full_and_simpleBRTS_chp3.png',webshot='webshot')

#(skip evaluating, because for the purpose of comparing we just need the end plots)

# Predict 

simple.models<-c('sppsimple','famsimple','gerrsimple')

for(i in simple.models){
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
	sf4preds<-sf4preds%>%mutate(btGerSmPD=exp(.$gerrsimplePD),.before='geometry')%>%
		rename(sppsimpPD=sppsimplePD,famsimpPD=famsimplePD,gerrsimpPD=gerrsimplePD)

	# save the updated sf4preds and df4preds
		st_write(sf4preds,'sf_for_predictions_fromBRTs_aug23.shp',driver='ESRI Shapefile',append=FALSE) # append command set to FALSE replaces the old file with the new one. If set to TRUE it would add the new data as further layers.
		st_write(sf4preds,'sf_for_predictions_fromBRTs_aug23.csv',delete_layer=TRUE,delete_dsn=TRUE)


# Plot 
	{sf4preds<-st_as_sf(st_read('sf_for_predictions_fromBRTs_aug23.shp'),crs='WGS84')
	sf4preds$jcode<-as.factor(sf4preds$jcode)
	sf4preds$Season<-as.factor(sf4preds$Season)}

	season.labels<-as_labeller(c('D'='Dry','W'='Wet'))

	summary(sf4preds) # spp 0.5-1.6. fam 0.1-1.5. gerr 14-30

	spp.simple.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=sppsimpPD),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0.5,1.7),guide=guide_colourbar(title=' Species \n\ Diversity \n\ Index'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
	fam.simple.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=famsimpPD),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0.1,1.4),guide=guide_colourbar(title=' Family \n\ Diversity \n\ Index'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
	gerr.simple.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=btGerSmPD),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(3,200),guide=guide_colourbar(title='Gerreidae\n\ Abundance\n\ (counts)'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
			sppsimple$contributions

	ggsave(gerr.simple.plot,file='resource_chp3/BRTS_outputs/simplified_BRT_Gerreidaebacktransformed_preds_aug23.png',device='png',units='in',height=8,width=10,dpi=1000)
	ggsave(spp.simple.plot,file='resource_chp3/BRTS_outputs/simplified_BRT_SpeciesSWbacktransformed_preds_aug23.png',device='png',units='in',height=8,width=10,dpi=1000)
	ggsave(fam.simple.plot,file='resource_chp3/BRTS_outputs/simplified_BRT_FamilySWbacktransformed_preds_aug23.png',device='png',units='in',height=8,width=10,dpi=1000)



###YOU LEFT OFF HERE



	# compare the predicted full vs predicted simple versus real distirbutions with boxplots (predicted) and points from raw data
		# real = joined_df
		# predicted = sf4preds

		# dataframe with one col for each models predictions
			{
			pred<-as.data.frame(sf4preds)%>%dplyr::select(Season,fmfllPD,sppflPD,bt_g_PD,famsimpPD,sppsimpPD,btGerSmPD) # has wet and dry season
			summary(pred)
			pred.long<-as.data.frame(pivot_longer(pred,cols=fmfllPD:btGerSmPD,names_to='model',values_to='preds'))
			pred.long$model<-as.factor(pred.long$model)
			summary(pred.long)
			}

		# dataframe with original data 
			{
			raw<-joined_df%>%dplyr::select(Season,SW_Species,SW_Families,Gerreidae)%>%rename(sppflPD=SW_Species,fmfllPD=SW_Families,bt_g_PD=Gerreidae)
			raw.long<-as.data.frame(pivot_longer(raw,cols=sppflPD:bt_g_PD,names_to='model',values_to='preds'))

			for.simples<-raw%>%rename(sppsimpPD=sppflPD,famsimpPD=fmfllPD,btGerSmPD=bt_g_PD)
			for.simples.long<-as.data.frame(pivot_longer(for.simples,cols=sppsimpPD:btGerSmPD,names_to='model',values_to='preds'))
			
			raw.long2<-bind_rows(raw.long,for.simples.long)
			raw.long2$model<-as.factor(raw.long2$model)
			summary(raw.long2)
			}


		# plot 
			season.labels<-as_labeller(c('D'='Dry','W'='Wet'))
			indices.labels<-as_labeller(c('fmfllPD'='SW Families (full)',
				'sppflPD'='SW Species (full)',
				'famsimpPD'='SW Families (simplified)',
				'sppsimpPD'='SW Species (simplified)'))
			gerr.labels<-as_labeller(c('bt_g_PD'='Gerreidae (full)','
				btGerSmPD'='Gerreidae (simplified)'))

			famspp_preds<-pred.long%>%dplyr::filter(model!=c('bt_g_PD','btGerSmPD'))			
			famspp_raw<-raw.long2%>%dplyr::filter(model!=c('bt_g_PD'))%>%dplyr::filter(model!=c('btGerSmPD'))
			summary(famspp_raw)

			indices.plot<-ggplot(pred.long%>%dplyr::filter(model!=c('bt_g_PD','btGerSmPD')),aes(x=model,y=preds))+
				geom_boxplot()+
				geom_jitter(data=famspp_raw,aes(x=model,y=preds),pch=19,alpha=0.7,col='orange3')+
				ylab(lab='Diversity Index Value')+
				xlab(lab='Model')+
				scale_x_discrete(labels=indices.labels)+
				theme_bw()+
				theme(axis.text.x = element_text(angle = 30,vjust=0.5,size=10))

			gerr.plot<-	ggplot(pred.long%>%dplyr::filter(model==c('bt_g_PD','btGerSimPD')),aes(x=model,y=preds))+geom_boxplot()+geom_jitter(data=raw.long2%>%dplyr::filter(model==c('bt_g_PD','btGerSimPD')),aes(x=model,y=preds),pch=19,alpha=0.7,col='orange3')+ylab(lab='Values (count)')+xlab(lab='Model')+scale_x_discrete(labels=gerr.labels)+theme_bw()+ theme(axis.text.x = element_text(angle = 30,vjust=0.5,size=10))

			comparison <- indices.plot + gerr.plot + plot_annotation( title = 'Predicted Values versus Observations, for Full and Simplified Models')

			ggsave(comparison,file='resource_chp3/BRTS_outputs/Predicted Values versus Observations, for Full and Simplified Models_AUGUST2023.png',device='png',units='in',height=6,width=12,dpi=1000)

			# Now - not moving forward with SW Families. So remake plots with Spp Index and re-combien with gridExtra with Gerreidae

			spponly.plot<-ggplot(famspp_preds%>%filter(model!='famsimpPD')%>%filter(model!='fmfllPD'),aes(x=model,y=preds))+geom_boxplot()+
				geom_jitter(data=famspp_raw%>%filter(model!='famsimpPD')%>%filter(model!='fmfllPD'),aes(x=model,y=preds),pch=19,alpha=0.7,col='orange3')+
				ylab(lab='Diversity Index Value')+
				xlab(lab='Model')+scale_x_discrete(labels=indices.labels)+
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

head(df4preds,2)
simpledf4preds<-df4preds%>%filter(Season!='W')%>%dplyr::select(-Season)
nrow(df4preds)
simplesf4preds<-sf4preds%>%filter(Season!='W')%>%dplyr::select(-Season)%>%dplyr::select(jcode,sppsimpPD,btGerSmPD)

simplified.predictions<-left_join(simpledf4preds,simplesf4preds,by='jcode')
head(simplified.predictions)


# as a simple feature
st_write(simplified.predictions,'resource_chp3/predictions_from_simplifiedBRTs_Gerreidae_and_SWSpecies_aug23.shp')
# as a dataframe
st_write(simplified.predictions,'resource_chp3/predictions_from_simplifiedBRTs_Gerreidae_and_SWSpecies_aug23.csv')


########################

## Evaluate simplified models - Gerreidae and SW Species

########################

### YOU LEFT OFF HERE 
	spp.hab.labels<-(c('dist2shore'='Dist. to Shore (m)', 'prp_lds'='Prop. of \n\ Low Density \n\ Seagrass','prp_mds'='Prop. of  \n\ Medium Density \n\ Seagrass','prp_hds'='Prop. of  \n\ High Density \n\ Seagrass')
	gerr.hab.labels<-(c('prp_lds'='Prop. of \n\ Low Density \n\ Seagrass','dist2shore'='Dist. to Shore (m)','prp_mds'='Prop. of  \n\ Medium Density \n\ Seagrass'))

# Influence of vars 

	infl.gerr<-gerrsimple$contributions
	infl.spp<-sppsimple$contributions

	spp.relinf<-infl.spp%>%mutate(var=fct_reorder(var,rel.inf))%>%
		ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+
		geom_bar(stat='identity')+
		scale_fill_distiller(direction=1,palette='Oranges',limits=c(0,70),guide='none')+
		theme_bw()+
		coord_flip()+
		scale_x_discrete(labels=hab.labels)+
		ylab('Relative Influence')+
		xlab(NULL)

		ggsave(spp.relinf,file='resource_chp3/BRTS_outputs/relative_influence_vars_in_simpifiedSWSpeciesBRT_aug23.png',device='png',units='in',height=6,width=8,dpi=900)
	
	gerr.relinf<-infl.gerr%>%mutate(var=fct_reorder(var,rel.inf))%>%ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+geom_bar(stat='identity')+scale_fill_distiller(direction=1,palette='Oranges',guide='none',limits=c(0,60))+theme_bw()+coord_flip()+scale_x_discrete(labels=gerr.hab.labels)+ylab('Relative Influence')+xlab(NULL)
		ggsave(gerr.relinf,file='resource_chp3/BRTS_outputs/relative_influence_vars_in_simplifiedGerreidaeBRT_aug23.png',device='png',units='in',height=6,width=8,dpi=900)

# Fitted values for dist2shore

		a1<-plot(sppsimple,'dist2shore',return.grid=TRUE)

		b1<-plot(gerrsimple,'dist2shore',return.grid=TRUE)

		dist2shorefitted_values_simplifiedmodels
		ggplot()+
			geom_smooth(data=a1,aes(x=dist2shore,y=y),fill='orange1',col='grey26',alpha=0.5)+
			geom_smooth(data=b1,aes(x=dist2shore,y=y),fill='darkorchid4',alpha=0.5,col='grey26')+
			xlab('Distance to Shore (m)')+
			ylab('Fitted Values')+
			scale_x_continuous(limits=c(0,1800))+
			scale_y_continuous()+
			theme_bw()

		ggsave(dist2shorefitted_values_simplifiedmodels,file='resource_chp3/BRTS_outputs/smoothed_fittedvalues_simplifiedSWSpeciesINorange_and_SimplifiedDist2ShoreINpurple_fittedvalues_BRT_feb23.png',device='png',unit='in',height=6,width=9,dpi=800)


	# compare ldsg and dist2shore in gerrsimple
	gbm.perspec(sppsimple,1,2)

	# compare 2 vars effects in sppsimple
	sppsimple.gbmint<-gbm.interactions(sppsimple)
	sppsimple.gbmint$rank.list
	gbm.perspec(sppsimple,3,2)





































