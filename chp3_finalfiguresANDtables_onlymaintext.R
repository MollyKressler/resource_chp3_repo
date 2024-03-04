### Final Figures code for Chapter 3: Risk or Resource

## compiling code for figures in final manuscript, not including supporting information tables and figures

#####
pacman::p_load(sf, tidyverse, lubridate, ggplot2, patchwork, cowplot, flextable,MuMIn,ggplot2,lme4,stats,ggeffects,gtsummary, MCMCvis)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd')
#####

###################
## Metadata on shark detections
###################

	meta <- read.csv('habitat_model/no ghosts of 6hr threshold/biometric_and_metadata_for_habsel_glmers_dec22.csv')%>%dplyr::select(-X)

	meta_flex<-flextable(meta%>%dplyr::select(-pit, -period)%>%
      mutate_if(is.numeric, round, digits = 3))%>%
		set_header_labels(sex = 'Sex',average_duration_hours = 'Average Duration\nbetween Detection (hours)',min_duration_hours = 'Minimum Duration\nbetween Detection (hours)',max_duration_days = 'Maximum Duration\nbetween Detection (days)',days = 'Days at\nLiberty',pcl = 'Pre-Caudal\nLength (cm)')%>%
		theme_zebra()%>%
		align(align = 'center', part = 'all')%>%
		font(fontname = 'Arial', part = 'all')%>%
		fontsize(size = 10, part = 'all')%>%
		autofit()
	meta_flex

	save_as_image(meta_flex, 'metadata_tags_inchapter3.png',webshot=webshot, res=850)


 ## metadata table on sizes of individuals 
	m <- read_csv('lemonspredators_20192020blacktipsANDbulls/all_predators.csv', col_select=c(Acoustic.ID:Longitude))%>%
		mutate(PredatorID = seq(1,49,1))
	d <- read_csv('lemonspredators_20192020blacktipsANDbulls/predators_detectedwithinstudyarea_MKthesis2019to2020_aug2023.csv')%>%
		mutate(Acoustic.ID = as.character(elasmo))
		
	a<- left_join(d, m%>%dplyr::select(PredatorID, Acoustic.ID), relationship = 'many-to-one')%>%
		filter(!is.na(PredatorID))%>%
		group_by(PredatorID)%>%
		tally()

	m2 <- left_join(m,a)%>%	
		filter(!is.na(n))%>%
		mutate(PredatorID = seq(1,21,1)) # resetting this for table readability only
	m2	# metadata for only sharks detected within period and area of juveniles detection set 

	meta.flextable <- m2 %>%
		mutate('Date (ymd)' = str_split(Date,' ', simplify=TRUE)[,1])%>%
		dplyr::select(PredatorID, Species, Sex, 'Date (ymd)', PCL, Project.Location,n)%>%
		flextable()%>%
		set_header_labels(PredatorID = 'Predator ID', Project.Location = 'Original Project \n\ Location', n = 'Number \n\ of Detections')%>%
		theme_zebra()%>%
		align(align = 'center', part = 'all')%>%
		font(fontname = 'Arial', part = 'all')%>%
		fontsize(size = 10, part = 'all')%>%
		autofit()
	meta.flextable
	save_as_image(meta.flextable, 'lemonspredators_20192020blacktipsANDbulls/onlypredatorsinstudyareaandtime_predators_PCL_tagDate_tagLocation.png', webshot='webshot')

###################
## Habitat data map - from winter 2020, from Emily Courmier 
###################

	## import polygonised data from Winter 2020
	data <- st_as_sf(st_read('winter2020habitat_polygonised_forSF_ECdata.shp'), crs='WGS84')
	labels <- as_labeller(c(lowdensg='Low Density Seagrass',mediumdensg='Medium Density Seagrass',highdensg='High Density Seagrass',vegetated='Vegetated (terrestrial)',bare.urban='Bare &/or Urban (terrestrial)'))
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	
	cmg<-st_as_sf(st_read('habitat_model/centroid_mangroves_north.shp'),crs='WSG84')
	r <- st_as_sf(st_read('receivers_in_thesis_data_dec2023.shp'),crs='WGS84')
	jettys<- st_as_sf(st_read('boatlaunchesBimini.kml'),crs='WGS84')
	nrow(jettys)

	winter2020map<-
	ggplot()+geom_sf(data=data,aes(fill=feature,col=feature),lwd=0)+
		scale_fill_manual(values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'),labels = labels)+
		scale_colour_manual(values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'),labels = labels)+
		geom_sf(data=r%>%filter(sharks=='no'),fill='grey82', col='white',pch=22,size=3)+
		geom_sf(data=r%>%filter(sharks=='yes'),fill='white', col='grey50',pch=21,size=3)+
		geom_sf(data=land,alpha=0,col='grey30')+
		geom_sf(data=cmg,pch=23,size=4,col='#FFFFFF',fill='#032D59')+
		geom_sf(data = jettys, pch = 24, size = 2, col='goldenrod1',fill='goldenrod1')+
		theme_bw()+ 
		theme(axis.text.x = element_text(angle = 90), legend.position = c(.75,0.17))
	winter2020map
		# save
		ggsave(winter2020map,file='habbimini2020_emilycourmier_seagrasses_urban_vegetated.png',device='png',unit='in',height=8,width=8,dpi=650)


###################
## Descriptive plots for predator metric
###################
	
	dd <- read_csv('lemonspredators_20192020blacktipsANDbulls/cooccurence_juvenilesANDlargesharks_within60minoflargedett.csv',col_select = c(rID:pressure))
	## individually
	a <- ggplot(data=dd%>%mutate(rID=as.character(rID)),aes(x=reorder(rID,-n.cooc),y=n.cooc))+
			geom_bar(stat='identity',fill='#3a6c74')+
			ylab('Co-occurences (total, count)')+
			xlab('Receiver')+
			theme_bw()	

	b <- ggplot(data=dd%>%mutate(rID=as.character(rID)),aes(x=reorder(rID,-cooc.days),y=cooc.days))+
			geom_bar(stat='identity',fill='#3a6c74')+
			ylab('Days of Co-occurences (total, count)')+
			xlab('Receiver')+
			theme_bw()

	pred.plots.forEquation <- a / b 
	pred.plots.forEquation		

	ggsave(pred.plots.forEquation,file='resource_chp3/cooccurences_anddaysof_predationpressureMetric_forsupplmentary.png',device=png,units='in',height=6,width=6,dpi=850)

	## or, in one plot, 2 y-axes
	ddd <- dd %>% 
		dplyr::select('rID','n.cooc','cooc.days')%>% 
		mutate(cooc.days = cooc.days*10)%>%
		pivot_longer(-rID) %>% 
		mutate(rID = as.character(rID))

	coeff <- 10
	colors <- c('#b6e5fc','#3a6c74')

	cooccur.plots <- ggplot(ddd, aes(x=reorder(rID,-value),y = value ,fill = name))+
		geom_col(position='dodge2')+
		scale_y_continuous('Co-occurrences (count)', sec.axis = sec_axis(~.*.10, name = 'Days of Co-occurrences (count * 10)'))+
		scale_fill_manual(values=colors,labels = c('Days of Co-oc.', 'Co-oc. Count'),name='')+
		xlab('Receiver')+
		theme_bw()+
		theme(legend.position = 'bottom')
	
	ggsave(cooccur.plots,file='resource_chp3/cooccurences_anddaysof_predationpressureMetric_forsupplmentary.png',device=png,units='in',height=4,width=6,dpi=850)

	c<-ggplot(data=dd%>%mutate(rID=as.character(rID)),aes(x=reorder(rID,-pressure),y=pressure))+
			geom_bar(stat='identity',fill='#3a6c74')+
			ylab('Pressure')+
			xlab('Receiver')+
			theme_bw()

	ggsave(c,file='resource_chp3/histogram_PRESSURE_byreceiver.png',device=png,units='in',height=3,width=6,dpi=850)


###################
## Descriptive plots for fishes
###################

 #032d59 dark blue

##### BEFORE BRTS - at each site for SuppInfo
	joined_df<-read.csv('resource_chp3/bruvs_DATAFRAME_joinedWITHhabitat_feb23.csv')
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	

	swspp.histogram <- ggplot(data = joined_df, aes(x=SW_Species))+
		geom_histogram(binwidth = 0.1,alpha=0.75, fill = '#032d59')+
			xlab('Species Diversity Index')+
			theme_bw()
	gerries.histogram <- ggplot(data = joined_df, aes(x=Gerreidae))+
		geom_histogram(binwidth = 10,alpha=0.75, fill = '#032d59')+
			xlab('Gerredae (count)')+
			theme_bw()
	bruvs.hist <- swspp.histogram / gerries.histogram

	ggsave(bruvs.hist,file='resource_chp3/histogram_SpeciesShannonIndex_and_Gerries_BRUVS_forSuppInfo.png',device=png,units='in',height=6,width=6,dpi=850)


##### AFTER BRTS - at each hexagon grid cell for results
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	
	after<-st_as_sf(st_read('resource_chp3/predictions_from_simplifiedBRTs_Gerreidae_and_SWSpecies_aug23.shp'),crs='WGS84')%>%
			dplyr::select(jcode,prp_lds,prp_mds,prp_hds,dist_cmg,dist2shore,sppsimpPD,btGerSmPD,geometry)%>%
			mutate(fishy=(sppsimpPD*btGerSmPD),.before=geometry)

	gerries <-	ggplot()+
		geom_sf(data=after,aes(fill=btGerSmPD),col=NA)+
		scale_fill_gradient(low='#b6e5fc',high='#3a6c74', space='Lab', aesthetics='fill',limits = c(0,40),guide=guide_colourbar(title='Predicted \n\ Gerridae'))+
		theme_bw()+
		theme(axis.text.x = element_text(angle = 90), legend.position="bottom",legend.direction = "horizontal")+
		geom_sf(data=land,col='grey30',fill='grey82', alpha = 0.7)

	spp <-	ggplot()+
		geom_sf(data=after,aes(fill=sppsimpPD),col=NA)+
		scale_fill_gradient(low='#b6e5fc',high='#3a6c74', space='Lab', aesthetics='fill',limits = c(0.5,1.5),guide=guide_colourbar(title='Predicted\n\ Species \n\ Diversity'))+
		theme_bw()+
		theme(axis.text.x = element_text(angle = 90), legend.position="bottom",legend.direction = "horizontal",axis.text.y=element_blank(),axis.ticks.y = element_blank())+
		geom_sf(data=land,col='grey30',fill='grey82', alpha = 0.7)

	fishy <-	ggplot()+
		geom_sf(data=after,aes(fill=fishy),col=NA)+
		scale_fill_gradient(low='#b6e5fc',high='#3a6c74', space='Lab', aesthetics='fill',limits = c(0,40),guide=guide_colourbar(title='Calculated \n\ Fish Metric'))+
		theme_bw()+
		theme(axis.text.x = element_text(angle = 90), legend.position="bottom",legend.direction = "horizontal",axis.text.y=element_blank(),axis.ticks.y = element_blank())+
		geom_sf(data=land,col='grey30',fill='grey82', alpha = 0.7)


	after.plots <- gerries | spp | fishy
	after.plots

	ggsave(after.plots,file='resource_chp3/postBRTpredictions_SpeciesShannonIndex_and_Gerries_and_FishyMetric_BRUVS_.png',device=png,units='in',height=8,width=10,dpi=850)
 ## QQ plots 
	gerr.simple.qq <- ggplot(data=preds, aes(sample = btGerSmPD))+
	stat_qq(size=1,pch=21)+
	labs(title='Simplified Gerridae BRT', subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
	stat_qq_line(linetype=2, col='red')+
	theme_bw()
	gerr.simple.qq
	ggsave(gerr.simple.qq,file = 'resource_chp3/BRTS_outputs/quantilequantile_plot_gerrsimple.png', dpi=850, unit= 'in', height = 4, width = 4)

	spp.simple.qq <- ggplot(data=preds, aes(sample = sppsimpPD))+
		stat_qq(size=1,pch=21)+
		labs(title='Simplified Species Diversity Index BRT', subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
		stat_qq_line(linetype=2, col='red')+
		theme_bw()
		spp.simple.qq
		ggsave(spp.simple.qq,file = 'resource_chp3/BRTS_outputs/quantilequantile_plot_sppsimple.png', dpi=850, unit= 'in', height = 4, width = 4)


###################
## Descriptive plots for Seagrasses PCA (raw data is in habitat map)
###################
	## seagrass PCA is standardised and mean-centred in hexdata
	withpca<-st_as_sf(st_read('resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_nov23.shp'),crs='WGS84')
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	

	sgPCAplot <- ggplot()+
		geom_sf(data=withpca,aes(fill=st_PCA1),col=NA)+
		scale_fill_gradient(low='#b6e5fc',high='#3a6c74', space='Lab', aesthetics='fill',limits = c(-3,3),guide=guide_colourbar(title='Seagrass PCA (standardised)'))+
		theme_bw()+
		theme(axis.text.x = element_text(angle = 90), legend.position="bottom",legend.direction = "horizontal")+
		geom_sf(data=land,col='grey30',fill='grey82', alpha = 0.7)
	
	ggsave(sgPCAplot,file='resource_chp3/seagrassesPCA_distributionplot_forchpater3writeup_.png',device=png,units='in',height=8,width=4.2,dpi=850)

###################
## Patterns of fish, predation pressure and seagrass based on abiotic covariates from 'best' dredge models
###################
	
	# read in model RDS
	sgm <- readRDS('resource_chp3/hypotesting_dredge_results/seagrasses_glm_hypotesting_distancemetrics_sept23.RDS')
	fm <- readRDS('resource_chp3/hypotesting_dredge_results/fishesmetric_glm_hypotesting_distancemetrics_sept23.RDS')
	pm <- readRDS('resource_chp3/hypotesting_dredge_results/largesharks_pressure_glm_hypotesting_distancemetrics_dec23.RDS')
	# get data for predicting into hexagons 
	hexdata<-read.csv('resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_nov23.csv')%>%
		mutate(jcode=as.numeric(jcode))%>%
		rename(standard.hexsgPCA1 = st_PCA1)
	pointdata<-read.csv('resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_dec23.csv')%>%
			mutate(standard.press=((as.numeric(pressure)-mean(as.numeric(pressure)))/sd(as.numeric(pressure)))) %>%
			mutate(standard.dist2jetty=((as.numeric(dist2jetty)-mean(as.numeric(dist2jetty)))/sd(as.numeric(dist2jetty))),.after=dist2jetty)


	navy: #3a6c74
	
	##################
	## seagrasses: dist2jetty, dist2shore, distcmg, dst2shore*distcmg
		
		sgm.hex.predicts.dist2jetty<-as.data.frame(ggpredict(sgm,terms='standard.dist2jetty[all]',type='fixed'))
		sgm.plot1<-	ggplot()+
			geom_point(data=hexdata, aes(x=standard.dist2jetty, y=standard.hexsgPCA1),pch=21,cex=0.5,col='#3a6c74')+
			ylab('Seagrass PCA (standardised)')+xlab('Dist. to Jetty (standardised)')+
			geom_line(data=sgm.hex.predicts.dist2jetty,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
			geom_ribbon(data=sgm.hex.predicts.dist2jetty,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
			theme_bw()

		sgm.hex.predicts.dist2shore<-as.data.frame(ggpredict(sgm,terms='standard.hexdist2shore[all]',type='fixed'))
		sgm.plot2<- ggplot()+
			geom_point(data=hexdata, aes(x=standard.hexdist2shore, y=standard.hexsgPCA1,xmin=min(hexdata$standard.hexdist2shore),xmax=max(standard.hexdist2shore)),pch=21,cex=0.5,col='#3a6c74')+
			ylab('Seagrass PCA (standardised)')+xlab('Dist. to Shore (standardised)')+
			geom_line(data=sgm.hex.predicts.dist2shore,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
			geom_ribbon(data=sgm.hex.predicts.dist2shore,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
			theme_bw()

		sgm.hex.predicts.dist2cmg<-as.data.frame(ggpredict(sgm,terms='standard.hexdistcmg[all]',type='fixed'))
		sgm.plot3<- ggplot()+
			geom_point(data=hexdata, aes(x=standard.hexdistcmg, y=standard.hexsgPCA1,xmin=min(hexdata$standard.hexdistcmg),xmax=max(standard.hexdistcmg)),pch=21,cex=0.5,col='#3a6c74')+
			ylab('Seagrass PCA (standardised)')+xlab('Dist. to Refuge (standardised)')+
			geom_line(data=sgm.hex.predicts.dist2cmg,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
			geom_ribbon(data=sgm.hex.predicts.dist2cmg,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
			theme_bw()

		## you have to calculate and plot the interaction terms differently. ggpredict puts the first term on the x-axis and then shows the second term at representative values - here three: the mean, and +- a SD unit 
		sgm.hex.predicts.intx.dist2.shore.cmg<-as.data.frame(ggpredict(sgm,terms=c("standard.hexdist2shore", "standard.hexdistcmg"),type='fixed'))
		
		intx.col.sgm <- c('#3a6c74','#708d8e','#3cbcfc')
		sgm.plot4<-
		ggplot()+
			geom_line(data=sgm.hex.predicts.intx.dist2.shore.cmg,aes(x=x,y=predicted, group=group, col=group),size=.6,linetype=1)+
			scale_color_manual(values=intx.col.sgm,name='Levels')+
			geom_ribbon(data=sgm.hex.predicts.intx.dist2.shore.cmg,aes(group=group,x=x,ymin=conf.low,ymax=conf.high,fill=group),alpha=0.3)+		
			scale_fill_manual(values=intx.col.sgm,name='Levels')+
			ylab('Seagrass PCA (standardised)')+
			xlab('Dist. to Shore (standardised')+
			theme_bw()+
		theme(legend.position='bottom')


		# put them together in a box
		sgm.prediction.plots.formatted <- (sgm.plot1 + sgm.plot2)/(sgm.plot3 + sgm.plot4)

		ggsave(sgm.prediction.plots.formatted,file='hypotesting_dredge_results/seescapescolortheme_seagrassesGLM_fromDredge_fourplots_predictions.png',device='png',units='in',dpi=350,height=6,width=7)


	##################
	## fishes: sgPCA1, dist2shore, distcmg, dst2shore*distcmg
		
		fm.hex.predicts.sgPCA<-as.data.frame(ggpredict(fm,terms='standard.hexsgPCA1[all]',type='fixed'))
		fm.plot1<- ggplot()+
			geom_point(data=hexdata, aes(x=standard.hexsgPCA1, y=standard.hexfish),pch=21,cex=0.5,col='#3a6c74')+
			ylab('Fishes (standardised)')+xlab('Seagrasses PCA (standardised)')+
			geom_line(data=fm.hex.predicts.sgPCA,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
			geom_ribbon(data=fm.hex.predicts.sgPCA,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
			theme_bw()

		fm.hex.predicts.distcmg<-as.data.frame(ggpredict(fm,terms='standard.hexdistcmg[all]',type='fixed'))
		fm.plot2<- ggplot()+
			geom_point(data=hexdata, aes(x=standard.hexdistcmg, y=standard.hexfish,xmin=min(hexdata$standard.hexdistcmg),xmax=max(standard.hexdistcmg)),pch=21,cex=0.5,col='#3a6c74')+
			ylab('Fishes (standardised)')+xlab('Dist. to Refuge (standardised)')+
			geom_line(data=fm.hex.predicts.distcmg,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
			geom_ribbon(data=fm.hex.predicts.distcmg,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
			theme_bw()

		fm.hex.predicts.dist2shore<-as.data.frame(ggpredict(fm,terms='standard.hexdist2shore[all]',type='fixed'))
		fm.plot3<- ggplot()+
			geom_point(data=hexdata, aes(x=standard.hexdist2shore, y=standard.hexfish,xmin=min(hexdata$standard.hexdist2shore),xmax=max(standard.hexdist2shore)),pch=21,cex=0.5,col='#3a6c74')+
			ylab('Fishes (standardised)')+xlab('Dist. to Shore (standardised)')+
			geom_line(data=fm.hex.predicts.dist2shore,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
			geom_ribbon(data=fm.hex.predicts.dist2shore,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
			theme_bw()

		## you have to calculate and plot the interaction terms differently. ggpredict puts the first term on the x-axis and then shows the second term at representative values - here three: the mean, and +- a SD unit 
		fm.hex.predicts.intx.dist2.shore.cmg<-as.data.frame(ggpredict(fm,terms=c("standard.hexdist2shore", "standard.hexdistcmg"),type='fixed'))
		
		intx.col.fm <- c('#3a6c74','#708d8e','#3cbcfc')
		fm.plot4<-
		ggplot()+
			geom_line(data=fm.hex.predicts.intx.dist2.shore.cmg,aes(x=x,y=predicted, group=group, col=group),size=.6,linetype=1)+
			scale_color_manual(values=intx.col.fm,name='Levels')+
			geom_ribbon(data=fm.hex.predicts.intx.dist2.shore.cmg,aes(group=group,x=x,ymin=conf.low,ymax=conf.high,fill=group),alpha=0.3)+		
			scale_fill_manual(values=intx.col.fm,name='Levels')+
			ylab('Fishes (standardised)')+
			xlab('Dist. to Shore (standardised)')+
			theme_bw()+
		theme(legend.position='bottom')


		# put them together in a box
		fm.prediction.plots.formatted <- (fm.plot1 + fm.plot2)/(fm.plot3 + fm.plot4)

		ggsave(fm.prediction.plots.formatted,file='hypotesting_dredge_results/seescapescolortheme_fishmetricGLM_fromDredge_fourplots_predictions.png',device='png',units='in',dpi=350,height=6,width=7)


	##################
	## large sharks: depth, dist2jetty, distcmg, depth*dist2jetty

	pm..predicts.dist2jetty<-as.data.frame(ggpredict(pm,terms='standard.dist2jetty[all]',type='fixed'))
	pm.plot1<- ggplot()+
		geom_point(data=pointdata, aes(x=standard.dist2jetty, y=standard.press),pch=21,col='#3a6c74')+
		ylab('Predation Pressure (standardised)')+xlab('Dist. to Jetty (standardised)')+
		geom_line(data=pm..predicts.dist2jetty,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
		geom_ribbon(data=pm..predicts.dist2jetty,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
		theme_bw()
		
	pm..predicts.distcmg<-as.data.frame(ggpredict(pm,terms='standard.dist2shore[all]',type='fixed'))
	pm.plot2<- ggplot()+
		geom_point(data=pointdata, aes(x=standard.dist2shore, y=standard.press,xmin=min(pointdata$standard.dist2shore),xmax=max(standard.dist2shore,ymin=min(pointdata$standard.press),ymax=max(pointdata$standard.press))),pch=21,col='#3a6c74')+
		ylab('Predation Pressure (standardised)')+xlab('Dist. to Shore (standardised)')+
		geom_line(data=pm..predicts.distcmg,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
		geom_ribbon(data=pm..predicts.distcmg,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
		theme_bw()


	pm..predicts.depth<-as.data.frame(ggpredict(pm,terms='standard.depth[all]',type='fixed'))
	pm.plot3<- ggplot()+
		geom_point(data=pointdata, aes(x=standard.depth, y=standard.press,xmin=min(pointdata$standard.depth),xmax=max(standard.depth)),pch=21,col='#3a6c74')+
		ylab('Predation Pressure (standardised)')+xlab('Depth (standardised)')+
		geom_line(data=pm..predicts.depth,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
		geom_ribbon(data=pm..predicts.depth,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
		theme_bw()

	## you have to calculate and plot the interaction terms differently. ggpredict puts the first term on the x-axis and then shows the second term at representative values - here three: the mean, and +- a SD unit 
	
	pm..predicts.intx.dist2.shore.cmg<-as.data.frame(ggpredict(pm,terms=c("standard.depth", "standard.dist2jetty"),type='fixed'))
	
	intx.col.pm <- c('#3a6c74','#708d8e','#3cbcfc')
	pm.plot4<-
	ggplot()+
		geom_line(data=pm..predicts.intx.dist2.shore.cmg,aes(x=x,y=predicted, group=group, col=group),size=.6,linetype=1)+
		scale_color_manual(values=intx.col.pm,name='Levels')+
		geom_ribbon(data=pm..predicts.intx.dist2.shore.cmg,aes(group=group,x=x,ymin=conf.low,ymax=conf.high,fill=group),alpha=0.3)+		
		scale_fill_manual(values=intx.col.pm,name='Levels')+
		ylab('Predation Pressure (standardised)')+
		xlab('Depth (standardised)')+
		theme_bw()+
		theme(legend.position='bottom')


	# put them together in a box
	pm.prediction.plots.formatted <- (pm.plot1 + pm.plot2)/(pm.plot3 + pm.plot4)

	ggsave(pm.prediction.plots.formatted,file='hypotesting_dredge_results/seescapescolortheme_largesharks_metricGLM_fromDredge_fourplots_predictions.png',device='png',units='in',dpi=350,height=6,width=7)



	## alternative formatting - 3 columns, four rows

	all.in.one.3col.4rows <- (sgm.plot1 / sgm.plot2 / sgm.plot3 / sgm.plot4) |	(fm.plot1 / fm.plot2 / fm.plot3 / fm.plot4) | 	(pm.plot1 / pm.plot2 / pm.plot3 / pm.plot4)

	ggsave(all.in.one.3col.4rows,file='hypotesting_dredge_results/seescapescolortheme_diffformatting_allin1_fishSGlargesharks_fromDredge_predictions.png',device='png',units='in',dpi=500,height=10,width=8.5)

  ## QQ plots 

	sgm.preds<-as.data.frame(predict(sgm))%>%rename(preds='predict(sgm)')
	fm.preds<-as.data.frame(predict(fm))%>%rename(preds='predict(fm)')
	pm.preds<-as.data.frame(predict(pm))%>%rename(preds='predict(pm)')

	sgm.qq <- ggplot(data=sgm.preds, aes(sample = preds))+
		stat_qq(size=1,pch=21)+
		labs(title='Seagrass PCA, Best Dredged Model', subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
		stat_qq_line(linetype=2, col='red')+
		theme_bw()
		sgm.qq
		ggsave(sgm.qq,file = 'resource_chp3/hypotesting_dredge_results/quantilequantile_plot_sgm.png', dpi=850, unit= 'in', height = 4, width = 4)

	fm.qq <- ggplot(data=fm.preds, aes(sample = preds))+
		stat_qq(size=1,pch=21)+
		labs(title='Prey Metric, Best Dredged Model', subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
		stat_qq_line(linetype=2, col='red')+
		theme_bw()
		fm.qq
		ggsave(fm.qq,file = 'resource_chp3/hypotesting_dredge_results/quantilequantile_plot_fm.png', dpi=850, unit= 'in', height = 4, width = 4)

	pm.qq <- ggplot(data=pm.preds, aes(sample = preds))+
		stat_qq(size=1,pch=21)+
		labs(title='Predation Pressure, Best Dredged Model', subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
		stat_qq_line(linetype=2, col='red')+
		theme_bw()
		pm.qq
		ggsave(pm.qq,file = 'resource_chp3/hypotesting_dredge_results/quantilequantile_plot_pm.png', dpi=850, unit= 'in', height = 4, width = 4)



	##############################################
	## Table: pathway description, estimates, CI and Support ##
	##############################################
		
		## For local macbook
		samplesList3b <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model3b_niter10000_burn2000_chains3_4dec2023.RDS')

		## For R Server
		pacman::p_load(tidybayes,bayesplot,MCMCvis,ggdist,nlist,forcats,patchwork)
		pacman::p_load(MCMCvis)
		
		samplesList3b <- readRDS('~/resource/data_and_RDS_NOTforupload/mcmcsamples_model3b_niter10000_burn2000_chains3_4dec2023.RDS')
		##
		

		samplesList3b %>% filter(Rowname=='path')
		
		

		pathwayresults_table <- MCMCsummary(samplesList3b,round=5,pg0=TRUE,params='path', probs=c(0.05,0.95))%>%
				tibble::rownames_to_column()%>%
				rename_with(str_to_title)%>%
			rename('pg0'='P>0')%>%
		  	mutate(pg00 = case_when(Mean > 0 ~ as.numeric(pg0), Mean < 0 ~ 1-as.numeric(pg0), .default = as.numeric(pg0)))%>%
			arrange(-pg00)%>%
			rename(Pathway = Rowname, 'Prop. of posterior with \n\ same sign as estimate' = 'pg00', Estimate = 'Mean','lower'='5%',upper='95%')%>%
			mutate('Path' = parse_number(Pathway), .before = 'Pathway')%>%  
			mutate(Estimate = round(Estimate, 3))%>%
			mutate(lower = case_when(Path != 1 ~ round(lower,3), Path == 1 ~ round(lower,5)),upper = case_when(Path != 1 ~ round(upper,3), Path == 1 ~ round(upper,5)))%>%
			mutate(CI = paste0('[',lower,',',upper,']'),.after='Estimate')%>%
			mutate('Standardised Estimate \n\ (95% credible interval)' = paste0(Estimate,' ',CI),.after='Pathway')%>%
			mutate('Path' = parse_number(Pathway), .before = 'Pathway')%>%  
			dplyr::select(-N.eff,-Rhat,-lower,-upper,-Estimate,-CI,-Sd, -pg0)%>%	
			flextable()%>%
				compose(i=1,j=2, as_paragraph('Juvenile sharks ~ Dist. to Refuge + Dist. to Shore + \n\ Seagrasses + Teleost fish'))%>%
				compose(i=2,j=2, as_paragraph('Juvenile sharks ~ Depth + Dist. to Shore +  \n\ Dist. to Jetty + Predator Pressure'))%>%
				compose(i=3,j=2, as_paragraph('Juvenile sharks ~ Dist. to Jetty + Dist. to Shore + \n\ Dist. to Refuge + Seagrasses'))%>%
				compose(i=4,j=2, as_paragraph('Juvenile sharks ~ Dist. to Refuge + Dist. to Jetty'))%>%
				theme_zebra()%>%
				align(j=3:4, align = 'center', part = 'all')%>%
				font(fontname = 'Arial', part = 'all')%>%
				color(color='black',part='all')%>%
				fontsize(size = 10, part = 'all')%>%
				autofit()
		pathwayresults_table

		save_as_image(pathwayresults_table,path='resource_chp3/nimblemodel_outputs/pathwayresultssummary_model3b_niter20000_burn2000_chains3_4dec2023.png')	



































##