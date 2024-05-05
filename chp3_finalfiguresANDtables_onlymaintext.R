### Final Figures code for Chapter 3: Risk or Resource

## compiling code for figures in final manuscript, not including supporting information tables and figures

#####
pacman::p_load(sf, tidyverse, lubridate, ggplot2, patchwork, cowplot, flextable,MuMIn,ggplot2,lme4,stats,ggeffects,gtsummary, MCMCvis)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd')
#####


###################
## Common dataframes and formatting
###################
    hexdata <- read.csv('resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_may24.csv')
    hexsf <- st_as_sf(st_read('resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_may24.shp'),crs='WGS84')%>%
      dplyr::select(-stndrd_, -stndrd_r,-sqzrisk)%>%
      rename(standard.hexshark = stndrd_hxs,
        standard.hexfish = stndrd_hxf,
        standard.hexdist2shore = stndrd_h2,
        standard.hexdistcmg = stndrd_hxd,
        standard.hexlowdensg = stndrd_hxl,
        standard.hexmeddensg = stndrd_hxm,
        standard.dist2jetty = stndrd_d2,
        standard.depth = stndrd_d,
        standard.sgPCA1 = st_PCA1,
        logit.sqzrisk = lgt_sqz,
        zlogit.sqzrisk = zlgt_sq,
        relPropPD = rlPrpPD
        )

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
## Descriptive plots for predator metric: relative risk at receivers
###################

	relp <- st_as_sf(st_read('lemonspredators_20192020blacktipsANDbulls/relativepredatorrisk_at_receivers_April2019December2020_lemonsANDblacktips.shp'))
    relp
   	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	

	
	relativePropPDdettsreceivers<-ggplot()+
	geom_sf(data=land,col='grey75',lwd=0.5)+
	geom_sf(data=relp,aes(size=ndetts,col=relPropPD))+
		scale_colour_gradient(low='navy',high='goldenrod2',limits=c(0.0,0.5),name='Relative Risk')+
		scale_size(name = 'No. of Detections')+
		theme_bw()
	relativePropPDdettsreceivers
	
	ggsave(relativePropPDdettsreceivers,file='lemonspredators_20192020blacktipsANDbulls/descriptive_stats_and_figures/relativepredatorrisk_at_receivers_April2019December2020_lemonsANDblacktips.png',device='png',units='in',dpi=950,height=7,width=6)


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



##############################################
## Table: pathway description, estimates, CI and Support ##
##############################################
	
	## For local macbook
	samplesList4 <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model4_niter12000_burn4000_chains3_4may2024.RDS')


	pathwayresults_table <- MCMCsummary(samplesList4,round=5,pg0=TRUE,params='path', probs=c(0.05,0.95))%>%
			tibble::rownames_to_column()%>%
			rename_with(str_to_title)%>%
		rename('pg0'='P>0')%>%
	  	mutate(pg00 = case_when(Mean > 0 ~ as.numeric(pg0), Mean < 0 ~ 1-as.numeric(pg0), .default = as.numeric(pg0)))%>%
		arrange(-pg00)%>%
		rename(Pathway = Rowname, 'Prop. of posterior with \n\ same sign as estimate' = 'pg00', Estimate = 'Mean','lower'='5%',upper='95%')%>%
		mutate('Path' = parse_number(Pathway), .before = 'Pathway')%>%  
		mutate(Estimate = round(Estimate, 4))%>%
		mutate(lower = case_when(Path != 1 ~ round(lower,3), Path == 1 ~ round(lower,5)),upper = case_when(Path != 1 ~ round(upper,4), Path == 1 ~ round(upper,4)))%>%
		mutate(CI = paste0('[',lower,',',upper,']'),.after='Estimate')%>%
		mutate('Standardised Estimate \n\ (95% credible interval)' = paste0(Estimate,' ',CI),.after='Pathway')%>%
		mutate('Path' = parse_number(Pathway), .before = 'Pathway')%>%  
		dplyr::select(-N.eff,-Rhat,-lower,-upper,-Estimate,-CI,-Sd, -pg0)%>%	
		flextable()%>%
			compose(i=1,j=2, as_paragraph('Juvenile sharks ~ Dist. to Refuge + Dist. to Shore + \n\ Seagrasses + Teleost fish'))%>%
			compose(i=2,j=2, as_paragraph('Juvenile sharks ~ Depth + Dist. to Shore +  \n\ Dist. to Jetty + Depth*Dist. to Jetty\n\ + Relative Risk'))%>%
			compose(i=3,j=2, as_paragraph('Juvenile sharks ~ Dist. to Jetty + Dist. to Shore + \n\ Dist. to Refuge + Seagrasses'))%>%
			compose(i=4,j=2, as_paragraph('Juvenile sharks ~ Dist. to Refuge + Dist. to Jetty'))%>%
			theme_zebra()%>%
			align(j=3:4, align = 'center', part = 'all')%>%
			font(fontname = 'Arial', part = 'all')%>%
			color(color='black',part='all')%>%
			fontsize(size = 10, part = 'all')%>%
			autofit()
	pathwayresults_table

	save_as_image(pathwayresults_table,path='resource_chp3/nimblemodel_outputs/pathwayresultssummary_model4_niter12000_burn4000_chains3_4may2024.png')	





































##