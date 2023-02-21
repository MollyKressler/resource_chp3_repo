

## created by Molly Kressler, 21 February 2023


## Load Workspace and Files 

pacman::p_load(tidyverse,sf,ggplot2,labeller,gridExtra,flextable,sf,ggsn,dismo,gbm)

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

	{sf4preds<-st_as_sf(st_read('sf_for_predictions_fromBRTs_feb23.shp'),crs='WGS84')
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


########################

## Instructions

# Run through the following code for each model. In the final step, saving the end outputs (an updated sf4preds and a plot of the predicted values), the file name needs to be changed to reference the appropriate model. 

########################


# Simplify - 

fullmodel<-gerrfull

simple<-gbm.simplify(fullmodel,n.drops=7)

simple.model<-gbm.step(joined_df,gbm.x=simple$pred.list[[7]],gbm.y=fullmodel$gbm.call$gbm.y,tree.complexity=fullmodel$gbm.call$tree.complexity,learning.rate=as.numeric(fullmodel$gbm.call$learning.rate),bag.fraction=fullmodel$gbm.call$bag.fraction,family=fullmodel$gbm.call$family,plot.main = TRUE) 


	# sppsimple<-simple.model # baresand removed
	# famsimple<-simple.model #baresand removed
	# gerrsimple<-simple.model #only ldsg and dist2shore left in
# Get n.trees and mean deviance of original & simplified models 

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
	info<-info%>%mutate(note=c(NA,NA,NA,'baresand removed','baresand removed','only low den. sg and dist2shore in'))%>%flextable()%>%theme_alafoli()%>%align(align = 'center', part = 'all')%>%font(fontname = 'Times', part = 'all')%>%fontsize(size = 11, part = 'all')%>%color(color='black',part='all')%>%autofit()

	save_as_image(info,'ntrees_and_deviance_of_full_and_simpleBRTS_chp3.png',webshot='webshot')

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
	sf4preds<-sf4preds%>%mutate(bt.gerrPD=exp(.$preds_Gerr),bt.gerr.simpPD=exp(.$gerrsimplePD),.before='geometry')

	# save the updated sf4preds and df4preds
		st_write(sf4preds,'sf_for_predictions_fromBRTs_feb23.shp',append=FALSE) # append command set to FALSE replaces the old file with the new one. If set to TRUE it would add the new data as further layers.
		st_write(sf4preds,'sf_for_predictions_fromBRTs_feb23.csv',delete_layer=TRUE,delete_dsn=TRUE)



# Plot 

	{sf4preds<-st_as_sf(st_read('sf_for_predictions_fromBRTs_feb23.shp'),crs='WGS84')%>%rename(btGerPD='bt_grPD',btGerSimPD='bt_g_PD')
	sf4preds$jcode<-as.factor(sf4preds$jcode)
	sf4preds$Season<-as.factor(sf4preds$Season)}
	# predictions for full models are stored in sf4preds, as preds_SWf, _SWsp, _gerr.

	season.labels<-as_labeller(c('D'='Dry','W'='Wet'))

	summary(sf4preds) # spp 0.5-1.6. fam 0.1-0.5. gerr 2.5-3.6

	spp.simple.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=sppsmPD),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0.5,1.7),guide=guide_colourbar(title=' Species \n\ Diversity \n\ Index'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
	fam.simple.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=fmsmpPD),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0.1,0.5),guide=guide_colourbar(title=' Family \n\ Diversity \n\ Index'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
	gerr.simple.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=btGerSimPD),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0,50),guide=guide_colourbar(title='Gerreidae\n\ Abundance\n\ (counts)'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
	

	ggsave(fam.simple.plot,file='figures+tables/simplified_removedBareSand_BRT_SWfamilies_preds_feb23.png',device='png',units='in',height=8,width=10,dpi=1000)



	# compare the predicted full vs predicted simple versus real distirbutions with boxplots (predicted) and points from raw data
		# real = joined_df
		# predicted = sf4preds

		# dataframe with one col for each models predictions

			pred<-as.data.frame(sf4preds)%>%dplyr::select(Season,prds_SWf,prds_SWs,btGerPD,fmsmpPD,sppsmPD,btGerSimPD) # has wet and dry season
			summary(pred)
			pred.long<-as.data.frame(pivot_longer(pred,cols=prds_SWf:btGerSimPD,names_to='model',values_to='preds'))
			pred.long$model<-as.factor(pred.long$model)
			summary(pred.long)

		# dataframe with original data 
			raw<-joined_df%>%dplyr::select(Season,SW_Species,SW_Families,Gerreidae)%>%rename(prds_SWs=SW_Species,prds_SWf=SW_Families,btGerPD=Gerreidae)
			raw.long<-as.data.frame(pivot_longer(raw,cols=prds_SWs:btGerPD,names_to='model',values_to='preds'))

			for.simples<-raw%>%rename(sppsmPD=prds_SWs,fmsmpPD=prds_SWf,btGerSimPD=btGerPD)
			for.simples.long<-as.data.frame(pivot_longer(for.simples,cols=sppsmPD:btGerSimPD,names_to='model',values_to='preds'))
			
			raw.long2<-bind_rows(raw.long,for.simples.long)
			raw.long2$model<-as.factor(raw.long2$model)
			summary(raw.long2)


		# plot 
			season.labels<-as_labeller(c('D'='Dry','W'='Wet'))
			indices.labels<-as_labeller(c('prds_SWf'='SW Families (full)','prds_SWs'='SW Species (full)','fmsmpPD'='SW Families (simplified)','sppsmPD'='SW Species (simplified)'))
			gerr.labels<-as_labeller(c('btGerPD'='Gerreidae (full)','btGerSimPD'='Gerreidae (simplified)'))

		
			indices.plot<-ggplot(pred.long%>%dplyr::filter(model!='btGerPD')%>%filter(model!='btGerSimPD'),aes(x=model,y=preds))+geom_boxplot()+geom_jitter(data=raw.long2%>%dplyr::filter(model!='btGerPD')%>%filter(model!='btGerSimPD'),aes(x=model,y=preds),pch=19,alpha=0.7,col='orange3')+ylab(lab='Diversity Index Value')+xlab(lab='Model')+scale_x_discrete(labels=indices.labels)+theme_bw()

			gerr.plot<-ggplot(pred.long%>%dplyr::filter(model==c('btGerPD','btGerSimPD')),aes(x=model,y=preds))+geom_boxplot()+geom_jitter(data=raw.long2%>%dplyr::filter(model==c('btGerPD','btGerSimPD')),aes(x=model,y=preds),pch=19,alpha=0.7,col='orange3')+ylab(lab='Values (count)')+xlab(lab='Model')+scale_x_discrete(labels=gerr.labels)+theme_bw()

			grid.arrange(indices.plot,gerr.plot,ncol=2,nrow=1,top=textGrob('Predicted Values versus Observations, for Full and Simplified Models',gp=gpar(fontsize=12,font=2)))













