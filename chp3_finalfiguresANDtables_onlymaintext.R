### Final Figures code for Chapter 3: Risk or Resource

## compiling code for figures in final manuscript, not including supporting information tables and figures

#####
pacman::p_load(sf, tidyverse, lubridate, ggplot2, patchwork, cowplot, flextable,MuMIn,ggplot2,lme4,stats,ggeffects,gtsummary, MCMCvis, viridis, tidybayes)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd')
#####


###################
## Common dataframes and formatting
###################
    hexdata <- read.csv('hexdata_juvLemonsPrUse_withAllCov_MKthesis.csv')%>%
	    mutate(st_shark = (m5_PUSE - mean(m5_PUSE))/sd(m5_PUSE))%>%
	    mutate(tide = case_when(tidephs == 'L' ~ 1, tidephs == 'H' ~ 0))%>%
      rename(standard.hexshark = st_shark,
        standard.hexfish = st_maxN,
        standard.hexdist2shore = st_d2sh,
        standard.hexdistcmg = st_dstc,
        standard.hexhighdensg = st_hds,
        standard.hexmeddensg = st_mds,
        standard.dist2jetty = st_d2jetty,
        standard.depth = st_depth,
        standard.SG = st_SG,
        zlogit.sqzrisk = st_risk
        )
    hexsf <- st_as_sf(st_read('hexdata_juvLemonsPrUse_withAllCov_MKthesis.shp'),crs='WGS84')%>%
	    mutate(st_shark = (m5_PUSE - mean(m5_PUSE))/sd(m5_PUSE))%>%
	    mutate(tide = case_when(tidephs == 'L' ~ 1, tidephs == 'H' ~ 0))%>%
      rename(standard.hexshark = st_shark,
        standard.hexfish = st_maxN,
        standard.hexdist2shore = st_d2sh,
        standard.hexdistcmg = st_dstc,
        standard.hexhighdensg = st_hds,
        standard.hexmeddensg = st_mds,
        standard.dist2jetty = st_d2jetty,
        standard.depth = st_depth,
        standard.SG = st_SG,
        zlogit.sqzrisk = st_risk
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

	save_as_docx(meta_flex, path = 'resource_chp3/metadata_tags_inchapter3.docx')

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
	labels <- as_labeller(c(lowdensg='Low Density Seagrass',mediumdensg='Medium Density Seagrass',highdensg='High Density Seagrass',vegetated='Vegetated (terrestrial)',bare.urban='Bare &/or Urban (terrestrial)'))
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	
	cmg<-st_as_sf(st_read('habitat_model/centroid_mangroves_north.shp'),crs='WSG84')
	r <- st_as_sf(st_read('receivers_in_thesis_data_dec2023.shp'),crs='WGS84')
	jettys<- st_as_sf(st_read('boatlaunchesBimini.kml'),crs='WGS84')
	nrow(jettys)

	h14 <- st_as_sf(st_read('winter2014habitat_polygonised_forSF_ECdata.shp'), crs='WGS84')
	h18 <- st_as_sf(st_read('summer2018habitat_polygonised_forSF_ECdata.shp'), crs='WGS84')
	h20 <- st_as_sf(st_read('winter2020habitat_polygonised_forSF_ECdata.shp'), crs = 'WGS84')


	## coverage of each type of habitat feature by area 
		total.area <- st_area(h20)/sum(st_area(h20))
		# values left to right: lds, mds, hds, vegetated, bare/urban
		# 0.271, 0.184, 0.127, 0.315, 0.100


	# hab maps 
	plot1 <- ggplot()+
		geom_sf(data=h14,aes(fill=feature,col=feature),lwd=0)+
		scale_fill_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den. Seagrass', 'Low Den. Seagrass', 'Medium Den. Seagrass', 'Vegetated'))+
		scale_colour_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den. Seagrass', 'Low Den. Seagrass', 'Medium Den. Seagrass', 'Vegetated'))+
		geom_sf(data = land, col = 'grey32',fill = NA, lwd = 0.25)+
		theme_bw()+
		theme(legend.position = 'bottom', legend.direction = 'horizontal')+
		ggtitle('2014')+
      scale_y_continuous(limits = c(25.665,25.78), breaks = c(25.70, 25.75))+
      scale_x_continuous(limits = c(-79.31,-79.23), breaks = c(-79.30,-79.25))
	plot2 <- ggplot()+
		geom_sf(data=h18,aes(fill=feature,col=feature),lwd=0)+
		scale_fill_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den. Seagrass', 'Low Den. Seagrass', 'Medium Den. Seagrass', 'Vegetated'))+
		scale_colour_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den. Seagrass', 'Low Den. Seagrass', 'Medium Den. Seagrass', 'Vegetated'))+
		geom_sf(data = land, col = 'grey32',fill = NA, lwd = 0.25)+
		theme_bw()+	
		theme(axis.text.y = element_blank(),legend.position = 'bottom', legend.direction = 'horizontal')+
		ggtitle('2018')+
      scale_y_continuous(limits = c(25.665,25.78), breaks = c(25.70, 25.75))+
      scale_x_continuous(limits = c(-79.31,-79.23), breaks = c(-79.30,-79.25))
	plot3 <- ggplot()+
		geom_sf(data=h20,aes(fill=feature,col=feature),lwd=0)+
		scale_fill_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den. Seagrass', 'Low Den. Seagrass', 'Medium Den. Seagrass', 'Vegetated'))+
		scale_colour_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den. Seagrass', 'Low Den. Seagrass', 'Medium Den. Seagrass', 'Vegetated'))+
		geom_sf(data = land, col = 'grey32', fill = NA, lwd = 0.25)+
		theme_bw()+
		theme(axis.text.y = element_blank(),legend.position = 'bottom', legend.direction = 'horizontal')+
		ggtitle('2020')+
      scale_y_continuous(limits = c(25.665,25.78), breaks = c(25.70, 25.75))+
      scale_x_continuous(limits = c(-79.31,-79.23), breaks = c(-79.30,-79.25))

	plot <- (plot1 | plot2 | plot3) + plot_layout(guides = "collect") & theme(legend.position = 'bottom', text = element_text(size = 11))
	ggsave(plot,file='chp3_ms_v6_figures_tables/habitat_bimini_wint14_summ18_wint20_ECdata.png',device='png',unit='in',height=5.5,width=8.5,dpi=950)

	# hab map, same hab data as above for 2020, but with additional information on jettys, receivers etc.
	plot3b <- ggplot()+
		geom_sf(data=h20,aes(fill=feature,col=feature),lwd=0)+
		geom_sf(data=r%>%filter(sharks=='no'),fill='grey82', col='white',pch=22,size=2)+
		geom_sf(data=r%>%filter(sharks=='yes'),fill='white', col='grey50',pch=21,size=2)+
		geom_sf(data = land, col = 'grey32', fill = NA, lwd = 0.25)+
		geom_sf(data=cmg,pch=23,size=3,col='#FFFFFF',fill='#032D59')+
		geom_sf(data = jettys, pch = 24, size = 1, col='goldenrod1',fill='goldenrod1')+
		scale_fill_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den. Seagrass', 'Low Den. Seagrass', 'Medium Den. Seagrass', 'Vegetated'))+
		scale_colour_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den. Seagrass', 'Low Den. Seagrass', 'Medium Den. Seagrass', 'Vegetated'))+
		theme_bw()+
		theme(legend.position = 'bottom', legend.direction = 'horizontal')+
		ggtitle('2020')+
      scale_y_continuous(limits = c(25.645,25.79), breaks = c(25.65,25.70, 25.75))+
      scale_x_continuous( breaks = c(-79.30,-79.25))

		plot <- (plot1 | plot2 | plot3b) + plot_layout(guides = "collect") & theme(legend.position = 'bottom', legend.justification = 'centre',text = element_text(size = 11))
		# save
		ggsave(plot,file='chp3_ms_v6_figures_tables/habitat_bimini_wint14_summ18_wint20_withCoVdata4BSEM_ECdata.png',device='png',unit='in',height=5, width = 8.5, dpi=800)

	# just 2020 alone with covariates data points 
	plot3c <- ggplot()+
		geom_sf(data=h20,aes(fill=feature,col=feature),lwd=0)+
		geom_sf(data=r%>%filter(sharks=='no'),fill='grey82', col='white',pch=22,size=2)+
		geom_sf(data=r%>%filter(sharks=='yes'),fill='white', col='grey50',pch=21,size=2)+
		geom_sf(data=cmg,pch=23,size=3,col='#FFFFFF',fill='#032D59')+
		geom_sf(data = jettys, pch = 24, size = 1, col='goldenrod1',fill='goldenrod1')+
		scale_fill_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den.\n\ Seagrass', 'Low Den.\n\ Seagrass', 'Medium Den.\n\ Seagrass', 'Vegetated'))+
		scale_colour_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den.\n\ Seagrass', 'Low Den.\n\ Seagrass', 'Medium Den.\n\ Seagrass', 'Vegetated'))+
		geom_sf(data = land, col = 'grey32', fill = NA, lwd = 0.25)+
		theme_bw()+
		theme(legend.position = 'bottom', legend.direction = 'horizontal')+
		ggtitle('2020')+
		guides(fill=guide_legend(nrow=2,byrow=TRUE))+
    scale_y_continuous(limits = c(25.645,25.79), breaks = c(25.65,25.70, 25.75))+
    scale_x_continuous( breaks = c(-79.30,-79.25))

		ggsave(plot3c,file='chp3_ms_v6_figures_tables/solohabbimini2020_withCoVforBSEM_ECdata.png',device='png',unit='in',width=5, height = 7,dpi=850)

###################
## BRUVS locations map
###################
		
		land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	
		data <- st_as_sf(st_read('bruvs2014AND2018_data_joinedWITHhabitat_summer18habANDwinter2014hab_july2024.shp'), crs= 'WGS84')
		### SD and HG data 
		bruvslocations <- ggplot()+
			geom_sf(data = land, col = 'grey75', lwd = 0.5)+
				geom_sf(data = data, aes(col = source, pch = source))+
				scale_color_manual(values = c('#020f75', '#7c1447'))+
				theme_bw()

		withsize <- ggplot()+
			geom_sf(data = land, col = 'grey75', lwd = 0.5)+
				geom_sf(data = data, aes(col = source, pch = source, size = maxN))+
				scale_color_manual(values = c('#020f75', '#7c1447'))+
				theme_bw()
		ggsave(withsize,file='bruvsdata_eval/bruvs_locations_hg_sd_bysource_withsize.png',device='png',unit='in',height=5, width = 5,dpi=850)


###################
## MaxN predictions
###################
	##### BEFORE BRTS - at each site for SuppInfo
	joined_df<-read.csv('bruvs2014AND2018_data_joinedWITHhabitat_summer18habANDwinter2014hab_july2024.csv')
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	

	maxN.histogram <- ggplot(data = joined_df, aes(x=log(maxN+.1)))+
		geom_histogram(binwidth = 0.1,alpha=0.75, fill = '#2b0052')+
			xlab('MaxN (log+0.1)')+
			theme_bw()

	ggsave(maxN.histogram,file='resource_chp3/histogram_maxN_BRUVS_forSuppInfo.png',device=png,units='in',height=3,width=6,dpi=850)
	ggsave(maxN.histogram,file='chp3_ms_v6_figures_tables/histogram_maxN_BRUVS_forSuppInfo.png',device='png',units='in',dpi=650, height = 3, width = 6)		
 	
 	#### BRT predictions

	data <- hexsf %>% dplyr::select(jcode, tidephs,standard.hexfish,mxN_prd, geometry)
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	
	labels = c('High', 'Low')
  names(labels) = c('H','L')

	maxN.plot<-ggplot()+		
		geom_sf(data = land, col = 'grey32', fill = 'grey75', lwd = 0.25)+
		geom_sf(data=data,aes(fill=mxN_prd),col=NA)+
		scale_fill_gradientn(colors = c('#fafcfc', '#2b0052'), limits = c(0,70), name = 'MaxN', oob = scales::squish)+
		theme_bw()+
		ggtitle('MaxN')+
	 	facet_wrap(~factor(tidephs, levels = c('H','L'), labels = labels), ncol=2)+
	  theme(legend.position = 'right',strip.background = element_rect(color = NA, fill = NA), strip.text = element_text(size = 10, hjust = 0), axis.text = element_text(size =10), plot.title = element_text(size = 10))+
	  guides(size = 'none')+
    scale_y_continuous(limits = c(25.665,25.78), breaks = c(25.70, 25.75))+
    scale_x_continuous(limits = c(-79.31,-79.23), breaks = c(-79.30,-79.25))

	ggsave(maxN.plot,file='chp3_ms_v6_figures_tables/simplified_BRT_LOG_maxN_preds_july24.png',device='png',units='in',dpi=650, height = 3.,width=5)		
 
	# calculate CI (method 1)
		result <- t.test(preds$maxN_preds)
		result$conf.int

		# (method 2) R program to find the confidence interval
		 
		# Calculate the mean of the sample data
		mean_value <- mean(preds$maxN_preds)
		 
		# Compute the size
		n <- length(preds$maxN_preds)
		 
		# Find the standard deviation
		standard_deviation <- sd(preds$maxN_preds)
		 
		# Find the standard error
		standard_error <- standard_deviation / sqrt(n)
		alpha = 0.05
		degrees_of_freedom = n - 1
		t_score = qt(p=alpha/2, df=degrees_of_freedom,lower.tail=F)
		margin_error <- t_score * standard_error
		 
		# Calculating lower bound and upper bound
		lower_bound <- mean_value - margin_error
		upper_bound <- mean_value + margin_error
		 
		# Print the confidence interval
		print(c(lower_bound,upper_bound)) # this is the same as call a t.test of the data column. 


###################
## Descriptive plots for predator metric: relative risk at receivers
###################
  N = 35
	relp <- st_as_sf(st_read('predators_pointdatasums_studyareaonly_withTidewithHab_MKthesis20192020.shp'))%>%
    mutate(geometry = st_centroid(geometry))%>%
		dplyr::select(buffID, tidephs, n, st_risk,geometry)
    relp
  land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	

	labels = c('High', 'Low')
  names(labels) = c('H','L')

	relativePropPDdettsreceivers<-ggplot()+
		geom_sf(data = land, col = 'grey32', fill = 'grey75', lwd = 0.25)+
	  geom_sf(data = relp, col = 'grey32', pch = 4, size = 0.75)+
	  geom_sf(data = relp%>%filter(n!=0), col = '#b27409', pch = 1, aes(size = n))+
	  scale_size_continuous(range = c(0,12))+
	 	facet_wrap(~factor(tidephs, levels = c('H','L'), labels = labels), ncol=2)+
	  theme_bw()+
		ggtitle('Predation Risk')+
	  theme(legend.position = 'bottom',strip.background = element_rect(color = NA, fill = NA), strip.text = element_text(size = 10, hjust = 0), axis.text = element_text(size =9), plot.title = element_text(size = 10))+
	  guides(size = 'none')+
    scale_y_continuous(limits = c(25.645,25.79), breaks = c(25.65,25.70, 25.75))+ 
    scale_x_continuous( breaks = c(-79.30,-79.25))
 
	ggsave(relativePropPDdettsreceivers,file='chp3_ms_v6_figures_tables/relativepredatorrisk_at_receivers.png',device='png',units='in',dpi=650, height = 3.25,width=4)


###################
## Descriptive plots for Seagrasses PCA (raw data is in habitat map)
###################
	
	## seagrass 
	datasg <- hexsf %>% dplyr::select(jcode, tidephs,prp_SG,prp_mds,prp_hds, geometry)
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	

	sgplot <- ggplot()+
		geom_sf(data = land, col = 'grey32', fill = 'grey75', lwd = 0.25)+
		geom_sf(data=datasg,aes(fill=prp_SG),col=NA)+
		scale_fill_gradientn(colors = c('#fafcfc', '#3A6C74'), limits = c(0,1), name = 'Proportion', oob = scales::squish)+
		theme_bw()+
		ggtitle('Seagrasses')+
	  theme(legend.position = 'bottom',axis.text = element_text(size =10), plot.title = element_text(size = 10))+
    scale_y_continuous(limits = c(25.665,25.78), breaks = c(25.70,25.75))+
    scale_x_continuous(limits = c(-79.31,-79.23), breaks = c(-79.30,-79.25))

	ggsave(sgplot,file='chp3_ms_v6_figures_tables/prop_of_seagrass_forchpater3writeup_.png',device=png,units='in',height=4, width=3,dpi=950)


###################
## Three mid level predictors: seagrasses PCA, fish, relative risk 
###################
	{
		#pipe for Figure 4 - remove some axes
	fish <- ggplot()+		
		geom_sf(data = land, col = 'grey32', fill = 'grey75', lwd = 0.25)+
		geom_sf(data=data,aes(fill=mxN_prd),col=NA)+
		scale_fill_gradientn(colors = c('#fafcfc', '#2b0052'), limits = c(0,70), name = 'MaxN', oob = scales::squish)+
		geom_sf(data=land,fill='gray75')+
		theme_bw()+
		#ggtitle('Prey Availability (maxN)')+
	 	facet_wrap(~factor(tidephs, levels = c('H','L'), labels = labels), ncol=2)+
	  theme(legend.position = 'bottom',strip.background = element_rect(color = NA, fill = NA), strip.text = element_text(size = 10, hjust = 0), axis.text = element_blank(), plot.title = element_text(size = 10))+
    scale_y_continuous(limits = c(25.645,25.79), breaks = c(25.65,25.70, 25.75))+
    scale_x_continuous( breaks = c(-79.30,-79.25), limits = c(-79.33,-79.21))
  sg <-  ggplot()+
		geom_sf(data = land, col = 'grey32', fill = 'grey75', lwd = 0.25)+
		geom_sf(data=datasg,aes(fill=prp_SG),col=NA)+
		scale_fill_gradientn(colors = c('#fafcfc', '#3A6C74'), limits = c(0,1), name = 'Proportion', oob = scales::squish)+
		theme_bw()+
		#ggtitle('Seagrasses')+
	  theme(legend.position = 'bottom',axis.text = element_text(size =10), plot.title = element_text(size = 10))+
    scale_y_continuous(limits = c(25.645,25.79), breaks = c(25.65,25.70, 25.75))+
    scale_x_continuous( breaks = c(-79.30,-79.25), limits = c(-79.33,-79.21))
  pd <- ggplot()+
		geom_sf(data = land, col = 'grey32', fill = 'grey75', lwd = 0.25)+
	  geom_sf(data = relp, col = 'grey32', pch = 4, size = 0.75)+
	  geom_sf(data = relp%>%filter(n!=0), col = '#b27409', pch = 1, aes(size = n))+
	  scale_size_continuous(range = c(0,12))+
	 	facet_wrap(~factor(tidephs, levels = c('H','L'), labels = labels), ncol=2)+
	  theme_bw()+
		#ggtitle('Predation Risk')+
	  theme(legend.position = 'bottom',strip.background = element_rect(color = NA, fill = NA), strip.text = element_text(size = 10, hjust = 0), axis.text.x = element_text(size =9), axis.text.y = element_blank(), plot.title = element_text(size = 10))+
	  guides(size = 'none')+
    scale_y_continuous(limits = c(25.645,25.79), breaks = c(25.65,25.70, 25.75))+ 
    scale_x_continuous( breaks = c(-79.30,-79.25), limits = c(-79.33,-79.21))

	theme_legend = theme(legend.box.margin = margin(-10,-5,0,-5), legend.direction = 'vertical', legend.box = 'vertical',legend.box.just = 'left', legend.spacing = unit(0, 'pt'), legend.margin = margin(0,0,0,0))

	sgleg <- get_legend(sg + theme_legend)
	fishleg <- get_legend(fish + theme_legend) 

	leg12 <- plot_grid(sgleg,fishleg,nrow=1, align = 'hv')
	sg2 <- sg + guides(fill = 'none')
	fish2 <- fish + guides(fill = 'none')
	layout <- "ABB
						 FCC"
	plot <-  sg2 + fish2 + pd + leg12 + plot_layout(design = layout,ncol = 3, nrow = 2, width = c(1,1,1), height = c(1,1)) & theme(plot.margin = unit(c(0,0,0,0),'cm'))

	ggsave(plot,file='chp3_ms_v6_figures_tables/sg_maxN_and_relRisk_chp3.png',device='png',units='in',dpi=1000)		
 	}  



##############################################
## Table: pathway description, estimates, CI and Support ##
##############################################

  samplesList6a <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model6_niter5000_burn1000_chains3_jan2025.RDS')

	pathwayresults_table <- MCMCsummary(samplesList6a,round=5,pg0=TRUE,params='path', probs=c(0.05,0.95))%>%
			tibble::rownames_to_column()%>%
			rename_with(str_to_title)%>%
		rename('pg0'='P>0')%>%
	  	mutate(pg00 = case_when(Mean > 0 ~ as.numeric(pg0), Mean < 0 ~ 1-as.numeric(pg0), .default = as.numeric(pg0)))%>%
		arrange(-pg00)%>%
		rename(Pathway = Rowname, 'Prop. of posterior with \n\ same sign as estimate' = 'pg00', Estimate = 'Mean','lower'='5%',upper='95%')%>%
		mutate('Path' = parse_number(Pathway), .before = 'Pathway')%>%  
		mutate(Estimate = round(Estimate, 4))%>%
		mutate(lower = case_when(Path == 3 ~ round(lower,5), Path !=3 ~ round(lower,5)),upper = case_when(Path == 3 ~ round(upper,4), Path !=3 ~ round(upper,5)))%>%
		mutate(CI = paste0('[',lower,',',upper,']'),.after='Estimate')%>%
		mutate('Standardised Estimate \n\ (95% credible interval)' = paste0(Estimate,' ',CI),.after='Pathway')%>%
		mutate('Path' = parse_number(Pathway), .before = 'Pathway')%>%  
		dplyr::select(-N.eff,-Rhat,-lower,-upper,-Estimate,-CI,-Sd, -pg0)%>%	
		flextable()%>%
			compose(i=1,j=2, as_paragraph('Juvenile sharks ~ Dist. to Refuge + Dist. to Shore + \n\ Depth + Tide State + Relative Risk'))%>%
			compose(i=2,j=2, as_paragraph('Juvenile sharks ~ Medium Density Seagrass +  \n\ High Density Seagrass +  Dist. to Jetty + Dist. to Shore\n\ + Tide State + Teleost fish'))%>%
			compose(i=3,j=2, as_paragraph('Juvenile sharks ~ Dist. to Jetty + Dist. to Shore + \n\ Dist. to Refuge + Seagrasses'))%>%
			align(j=3:4, align = 'center', part = 'all')%>%
			align(j=1, align = 'center', part = 'all')%>%
			font(fontname = 'Arial', part = 'all')%>%
			color(color='black',part='all')%>%
			fontsize(size = 10, part = 'all')%>%
			autofit()%>%
			theme_zebra()
	pathwayresults_table

	save_as_image(pathwayresults_table,path='resource_chp3/nimblemodel_outputs/pathwayresultssummary_model5a_niter5000_burn1000_chains3_Jan2025.png')	
	save_as_docx(pathwayresults_table,path='resource_chp3/nimblemodel_outputs/pathwayresultssummary_model5a_niter5000_burn1000_chains3_Jan2025.docx')	
	# for manuscript v6
	save_as_image(pathwayresults_table,path='chp3_ms_v6_figures_tables/pathwayresultssummary_model5a_niter5000_burn1000_chains3_Jan2025.png'
	save_as_docx(pathwayresults_table,path='chp3_ms_v6_figures_tables/pathwayresultssummary_model5a_niter5000_burn1000_chains3_Jan2025.docx')	

##############################################
## Caterpillars ##
##############################################
 # grab draws with gather_draws and create label for paths based on iterations and sequence of paths minN to maxN. 
    
    d6 <- gather_draws(samplesList6a,path[])%>%
      group_by(.chain)%>%
      mutate(pathID = paste0('path',rep(1:3, each=8000)))%>% 
      mutate(pathIDnum = rep(1:3, each=8000))%>%
      ungroup() # each path estimate for each chain (of 3) in order, starting with path[1] first estimate in chain 1 
    
    # use tidybayes to plot 
    
    caterpillars <- ggplot(d6, aes(y=reorder(pathID,pathIDnum,decreasing=T),x=.value, col = pathID))+
      stat_pointinterval(.width=c(.50,.95),point_size=2)+
      scale_color_manual(values = c('#3A6C74','#2b0052','#b27409'))+
      ylab('Path ID')+
      xlab('Estimate (mean) & CI (.5,.95)')+
      geom_vline(xintercept=0,linetype=3)+
      guides(col = 'none')+
      theme_bw()
    caterpillars  
    ggsave(caterpillars,file='resource_chp3/nimblemodel_outputs/caterpillarsPlot_mcmcsamples_model6_niter5000_burn1000_chains3_jan2025.png',device='png',dpi=650,width=4,height=3.5,units='in')
    ##for manuscript v6
    ggsave(caterpillars,file='chp3_ms_v6_figures_tables/caterpillarsPlot_mcmcsamples_model6_jan2025.png',device='png',dpi=650,width=4,height=3.5,units='in')
 

##############################################
## Plots from predictions of effects from path 3 into 2020 hexdata ##
##############################################









##