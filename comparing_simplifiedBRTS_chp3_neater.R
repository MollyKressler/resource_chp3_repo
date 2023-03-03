

## created by Molly Kressler, 21 February 2023


## Load Workspace and Files 

pacman::p_load(tidyverse,sf,ggplot2,gridExtra,flextable,sf,ggsn,dismo,gbm)

# shapefiles, useful for plotting

	setwd('/Users/mollykressler/Documents/data_phd')
	hab.grid<-st_as_sf(st_read('hexagon_grid_withUPDATEDhabitat_sept22.shp'),crs='WGS84')
	setwd('/Users/mollykressler/Documents/data_phd/resource_chp3')
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')
	hab.noland<-st_as_sf(st_read('hexagon_grid_study_site_with_habitatFromSOSF_forchp3BRUVS.shp'),crs='WGS84')
# Data frames and shape files with data 

	setwd('/Users/mollykressler/Documents/data_phd/resource_chp3')

	{df4preds<-read.csv('df_for_prediction_chp3BRUVS_habitat_noland.csv')%>%dplyr::select(-X)
		df4preds$jcode<-as.factor(df4preds$jcode)
		df4preds$Season<-as.factor(df4preds$Season)}

	{sf4preds<-st_as_sf(st_read('sf_for_predictions_fromBRTs_feb23.shp'),crs='WGS84')%>%dplyr::select(-sppsmPD,-fmsmpPD,-grrsmPD,-bt_g_PD)
		sf4preds$jcode<-as.factor(sf4preds$jcode)
		sf4preds$Season<-as.factor(sf4preds$Season)}
		# predictions for full models are stored in sf4preds, as preds_SWf, _SWsp, _gerr.

	{joined_df<-read.csv('bruvs_DATAFRAME_joinedWITHhabitat_feb23.csv')%>%dplyr::select(-X)
		joined_df$BRUV<-as.factor(joined_df$BRUV)
		joined_df$Season<-as.factor(joined_df$Season)
		joined_df$SW_Species<-as.numeric(joined_df$SW_Species)
		joined_df$SW_Families<-as.numeric(joined_df$SW_Families)	}
		
# Model RDS - BRTS, with all 9 predictor variables, 'full'

	sppfull<-readRDS('model_RDS/SW_Species_brt_tc3_lr001_gaussian.RDS')
	famfull<-readRDS('model_RDS/SW_Families_brt_tc3_lr001_gaussian.RDS')
	gerrfull<-readRDS('model_RDS/Gerreidae_brt_tc3_lr0001_poisson.RDS')

# Model RDS - BRTS, simplified, 'simple'. Not moving forward with Family Shannon Index. Had to re-run both on 27 Feb, bc computer crashed before I could save the RDS. 
	# sppsimple<-simple.model # baresand removed
	# gerrsimple<-simple.model #only ldsg and dist2shore left in.



########################

## Instructions

# Run through the following code for each model. In the final step, saving the end outputs (an updated sf4preds and a plot of the predicted values), the file name needs to be changed to reference the appropriate model. 

########################


# Simplify - 

fullmodel<-famfull

simple<-gbm.simplify(fullmodel,n.drops=7)

simple.model<-gbm.step(joined_df,gbm.x=simple$pred.list[[2]],gbm.y=fullmodel$gbm.call$gbm.y,tree.complexity=fullmodel$gbm.call$tree.complexity,learning.rate=as.numeric(fullmodel$gbm.call$learning.rate),bag.fraction=fullmodel$gbm.call$bag.fraction,family=fullmodel$gbm.call$family,plot.main = TRUE) 

#saveRDS(simple.model,'model_RDS/simplified_familyShannonIndex_BRT_feb23_baresandSeasonremoved.RDS')

	# sppsimple<-simple.model # baresand and Season removed
	# famsimple<-simple.model # baresand and Season removed
	# gerrsimple<-simple.model #only ldsg and dist2shore left in.


# Get n.trees and mean deviance of original & simplified models 

	sppsimple<-readRDS('model_RDS/simplified_speciesShannonIndex_BRT_feb23_baresandSeasonremoved.RDS')
	famsimple<-readRDS('model_RDS/simplified_familyShannonIndex_BRT_feb23_baresandSeasonremoved.RDS')
	gerrsimple<-readRDS('model_RDS/simplified_gerreidae_BRT_feb23_lowdensitySG_dist2shore.RDS')

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
	info2<-info%>%mutate(note=c(NA,NA,NA,'baresand & Season removed','baresand & Season removed','only low den. sg and dist2shore in'))%>%flextable()%>%theme_alafoli()%>%align(align = 'center', part = 'all')%>%font(fontname = 'Times', part = 'all')%>%fontsize(size = 11, part = 'all')%>%color(color='black',part='all')%>%autofit()

	save_as_image(info2,'figures+tables/ntrees_and_deviance_of_full_and_simpleBRTS_chp3.png',webshot='webshot')

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
	sf4preds<-sf4preds%>%mutate(bt.gerr.simpPD=exp(.$gerrsimplePD),.before='geometry')%>%rename(sppsimpPD=sppsimplePD,famsimpPD=famsimplePD,gerrsimpPD=gerrsimplePD)

	# save the updated sf4preds and df4preds
		st_write(sf4preds,'sf_for_predictions_fromBRTs_feb23.shp',append=FALSE) # append command set to FALSE replaces the old file with the new one. If set to TRUE it would add the new data as further layers.
		st_write(sf4preds,'sf_for_predictions_fromBRTs_feb23.csv',delete_layer=TRUE,delete_dsn=TRUE)



# Plot 

	{sf4preds<-st_as_sf(st_read('sf_for_predictions_fromBRTs_feb23.shp'),crs='WGS84')%>%rename(btGerPD='bt_grPD',btGerSimPD='bt_g_PD')
	sf4preds$jcode<-as.factor(sf4preds$jcode)
	sf4preds$Season<-as.factor(sf4preds$Season)}
	# predictions for full models are stored in sf4preds, as preds_SWf, _SWsp, _gerr.

	season.labels<-as_labeller(c('D'='Dry','W'='Wet'))

	summary(sf4preds) # spp 0.5-1.6. fam 0.1-0.5. gerr 3-200

	spp.simple.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=sppsmPD),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0.5,1.7),guide=guide_colourbar(title=' Species \n\ Diversity \n\ Index'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
	fam.simple.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=fmsmpPD),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0.1,0.5),guide=guide_colourbar(title=' Family \n\ Diversity \n\ Index'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
	gerr.simple.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=btGerSimPD),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(3,200),guide=guide_colourbar(title='Gerreidae\n\ Abundance\n\ (counts)'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
			sppsimple$contributions

	ggsave(gerr.simple.plot,file='figures+tables/simplified_onltLDSGanddist2shore_BRT_Gerreidaebacktransformed_preds_feb23.png',device='png',units='in',height=8,width=10,dpi=1000)



	# compare the predicted full vs predicted simple versus real distirbutions with boxplots (predicted) and points from raw data
		# real = joined_df
		# predicted = sf4preds

		# dataframe with one col for each models predictions
			{
			pred<-as.data.frame(sf4preds)%>%dplyr::select(Season,prds_SWf,prds_SWs,btGerPD,fmsmpPD,sppsmPD,btGerSimPD) # has wet and dry season
			summary(pred)
			pred.long<-as.data.frame(pivot_longer(pred,cols=prds_SWf:btGerSimPD,names_to='model',values_to='preds'))
			pred.long$model<-as.factor(pred.long$model)
			summary(pred.long)
			}

		# dataframe with original data 
			{
				raw<-joined_df%>%dplyr::select(Season,SW_Species,SW_Families,Gerreidae)%>%rename(prds_SWs=SW_Species,prds_SWf=SW_Families,btGerPD=Gerreidae)
			raw.long<-as.data.frame(pivot_longer(raw,cols=prds_SWs:btGerPD,names_to='model',values_to='preds'))

			for.simples<-raw%>%rename(sppsmPD=prds_SWs,fmsmpPD=prds_SWf,btGerSimPD=btGerPD)
			for.simples.long<-as.data.frame(pivot_longer(for.simples,cols=sppsmPD:btGerSimPD,names_to='model',values_to='preds'))
			
			raw.long2<-bind_rows(raw.long,for.simples.long)
			raw.long2$model<-as.factor(raw.long2$model)
			summary(raw.long2)
			}


		# plot 
			season.labels<-as_labeller(c('D'='Dry','W'='Wet'))
			indices.labels<-as_labeller(c('prds_SWf'='SW Families (full)','prds_SWs'='SW Species (full)','fmsmpPD'='SW Families (simplified)','sppsmPD'='SW Species (simplified)'))
			gerr.labels<-as_labeller(c('btGerPD'='Gerreidae (full)','btGerSimPD'='Gerreidae (simplified)'))

		
			indices.plot<-ggplot(pred.long%>%dplyr::filter(model!='btGerPD')%>%filter(model!='btGerSimPD'),aes(x=model,y=preds))+geom_boxplot()+geom_jitter(data=raw.long2%>%dplyr::filter(model!='btGerPD')%>%filter(model!='btGerSimPD'),aes(x=model,y=preds),pch=19,alpha=0.7,col='orange3')+ylab(lab='Diversity Index Value')+xlab(lab='Model')+scale_x_discrete(labels=indices.labels)+theme_bw()+ theme(axis.text.x = element_text(angle = 30,vjust=0.5,size=10))

			gerr.plot<-	ggplot(pred.long%>%dplyr::filter(model==c('btGerPD','btGerSimPD')),aes(x=model,y=preds))+geom_boxplot()+geom_jitter(data=raw.long2%>%dplyr::filter(model==c('btGerPD','btGerSimPD')),aes(x=model,y=preds),pch=19,alpha=0.7,col='orange3')+ylab(lab='Values (count)')+xlab(lab='Model')+scale_x_discrete(labels=gerr.labels)+theme_bw()+ theme(axis.text.x = element_text(angle = 30,vjust=0.5,size=10))

			comparison<-grid.arrange(indices.plot,gerr.plot,ncol=2,nrow=1,top=textGrob('Predicted Values versus Observations, for Full and Simplified Models',gp=gpar(fontsize=12,font=2)))

			ggsave(comparison,file='figures+tables/Predicted Values versus Observations, for Full and Simplified Models.png',device='png',units='in',height=6,width=12,dpi=1000)

			# Now - not moving forward with SW Families. So remake plots with Spp Index and re-combien with gridExtra with Gerreidae

			spponly.plot<-ggplot(pred.long%>%dplyr::filter(model!='btGerPD')%>%filter(model!='btGerSimPD')%>%filter(model!='fmsmpPD')%>%filter(model!='prds_SWf'),aes(x=model,y=preds))+geom_boxplot()+geom_jitter(data=raw.long2%>%dplyr::filter(model!='btGerPD')%>%filter(model!='btGerSimPD')%>%filter(model!='fmsmpPD')%>%filter(model!='prds_SWf'),aes(x=model,y=preds),pch=19,alpha=0.7,col='orange3')+ylab(lab='Diversity Index Value')+xlab(lab='Model')+scale_x_discrete(labels=indices.labels)+theme_bw()+ theme(axis.text.x = element_text(angle = 30,vjust=0.5,size=10))


				comparison2<-grid.arrange(spponly.plot,gerr.plot,ncol=2,nrow=1,top=textGrob('Predicted Values versus Observations, for Full and Simplified Models',gp=gpar(fontsize=12,font=2)))

				ggsave(comparison2,file='figures+tables/Predicted Values versus Observations, for Full and Simplified Models, no Family Diversity Index.png',device='png',units='in',height=6,width=12,dpi=1000)
			
			# only simplified models

			spponly.simpleonly.plot<-ggplot(pred.long%>%dplyr::filter(model=='sppsmPD'),aes(x=model,y=preds))+geom_boxplot()+geom_jitter(data=raw.long2%>%dplyr::filter(model=='sppsmPD'),aes(x=model,y=preds),pch=19,alpha=0.7,col='orange3')+ylab(lab='Diversity Index Value')+xlab(lab=NULL)+scale_x_discrete(labels=indices.labels)+theme_bw()

			gerr.simpleonly.plot<-ggplot(pred.long%>%dplyr::filter(model=='btGerSimPD'),aes(x=model,y=preds))+geom_boxplot()+geom_jitter(data=raw.long2%>%dplyr::filter(model=='btGerSimPD'),aes(x=model,y=preds),pch=19,alpha=0.7,col='orange3')+ylab(lab='Values (count)')+xlab(lab=NULL)+scale_x_discrete(labels=gerr.labels)+theme_bw()

				comparison3<-grid.arrange(spponly.simpleonly.plot,gerr.simpleonly.plot,ncol=2,nrow=1,top=textGrob('Predicted Values versus Observations for Simplified Models',gp=gpar(fontsize=12,font=2)))

				ggsave(comparison3,file='figures+tables/Predicted Values versus Observations, Simplified models of Species Diveristy and Gerreidae counts.png',device='png',units='in',height=6,width=7,dpi=1000)



#### Save predictions in a format for modelling with further in chapter 3 - PREY LANDSCAPE MODEL outputs. 

head(df4preds2,2)
simpledf4preds<-df4preds%>%filter(Season!='W')%>%dplyr::select(-Season)
nrow(df4preds2)
simplesf4preds<-sf4preds%>%filter(Season!='W')%>%dplyr::select(-Season)%>%dplyr::select(jcode,sppsmPD,btGerSimPD)

simplified.predictions<-left_join(simpledf4preds,simplesf4preds,by='jcode')
head(simplified.predictions)


# as a simple feature
st_write(simplified.predictions,'predictions_from_simplifiedBRTs_Gerreidae_and_SWSpecies_feb23.shp')
# as a dataframe
st_write(simplified.predictions,'predictions_from_simplifiedBRTs_Gerreidae_and_SWSpecies_feb23.csv')


########################

## Evaluate simplified models - Gerreidae and SW Species

########################


	spp.hab.labels<-(c('dist2shore'='Dist. to Shore (m)', 'prop_ldsg'='Prop. of \n\ Low Density \n\ Seagrass','prop_medsg'='Prop. of  \n\ Medium Density \n\ Seagrass','prop_hdsg'='Prop. of  \n\ High Density \n\ Seagrass','prop_sarg'='Prop. of Sargassum','prop_urb_r'='Prop. of \n\ Urban & Rocky','prop_deep'='Prop. of \n\ Deep Water'))
	gerr.hab.labels<-(c('prop_ldsg'='Prop. of \n\ Low Density \n\ Seagrass','dist2shore'='Dist. to Shore (m)'))

# Influence of vars 

	infl.gerr<-gerrsimple$contributions
	infl.spp<-sppsimple$contributions

	spp.relinf<-infl.spp%>%mutate(var=fct_reorder(var,rel.inf))%>%ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+geom_bar(stat='identity')+scale_fill_distiller(direction=1,palette='Oranges',limits=c(0,45),guide='none')+theme_bw()+coord_flip()+scale_x_discrete(labels=hab.labels)+ylab('Relative Influence')+xlab(NULL)

		ggsave(spp.relinf,file='relative_influence_vars_in_simpifiedSWSpeciesBRT_feb23.png',device='png',units='in',height=6,width=8,dpi=900)
	
	gerr.relinf<-infl.gerr%>%mutate(var=fct_reorder(var,rel.inf))%>%ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+geom_bar(stat='identity')+scale_fill_distiller(direction=1,palette='Oranges',guide='none',limits=c(0,60))+theme_bw()+coord_flip()+scale_x_discrete(labels=gerr.hab.labels)+ylab('Relative Influence')+xlab(NULL)
		ggsave(gerr.relinf,file='relative_influence_vars_in_simplifiedGerreidaeBRT_feb23.png',device='png',units='in',height=6,width=8,dpi=900)

# Fitted values for dist2shore

		a1<-plot(sppsimple,'dist2shore',return.grid=TRUE)

		b1<-plot(gerrsimple,'dist2shore',return.grid=TRUE)
		dist2shorefitted_values_simplifiedmodels<-ggplot()+geom_smooth(data=a1,aes(x=dist2shore,y=y),fill='orange1',alpha=0.5,col='grey26')+geom_smooth(data=b1,aes(x=dist2shore,y=y),fill='darkorchid4',alpha=0.5,col='grey26')+xlab('Distance to Shore (m)')+ylab('Fitted Values')+scale_x_continuous(limits=c(0,1800))+scale_y_continuous()+theme_bw()

		ggsave(dist2shorefitted_values_simplifiedmodels,file='smoothed_fittedvalues_simplifiedSWSpeciesINorange_and_SimplifiedDist2ShoreINpurple_fittedvalues_BRT_feb23.png',device='png',unit='in',height=6,width=9,dpi=800)


	# compare ldsg and dist2shore in gerrsimple
	gbm.perspec(gerrsimple,1,2)

	# compare 2 vars effects in sppsimple
	sppsimple.gbmint<-gbm.interactions(sppsimple)
	sppsimple.gbmint$rank.list
	gbm.perspec(sppsimple,7,2)





































