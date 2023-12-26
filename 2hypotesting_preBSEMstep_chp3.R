## Hypothesis testing for chapter 3 (resource and risk BSEM) 

## inspired by methods in Bennett et al. Ecology & Evolution

## created 18 September 2023 by Molly M Kressler 

# Three hypotheses tested with the 'dredge' function in the package 'MuMIn' (Multi-Model inference)
https://rdrr.io/cran/MuMIn/man/dredge.html 

## Updated 6 October 2023: post-meeting with Rich. Fundamentals approved, but some tweaks based on how bad the seagrass PCA models are (residuals and linear assumptions not met), and some ideas about fish predictors. 


#############################
## - Load Workspace and Data
#############################

pacman::p_load(tidyverse,MuMIn,ggplot2,flextable,cowplot,patchwork,lme4,stats,ggeffects,gtsummary)
  options(na.action = "na.fail")

setwd('/Users/mollykressler/Documents/data_phd/resource_chp3')

pointdata<-read.csv('standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_dec23.csv')%>%mutate(standard.press=((as.numeric(pressure)-mean(as.numeric(pressure)))/sd(as.numeric(pressure)))) %>%
			mutate(standard.dist2jetty=((as.numeric(dist2jetty)-mean(as.numeric(dist2jetty)))/sd(as.numeric(dist2jetty))),.after=dist2jetty)

hexdata<-read.csv('data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_aug23.csv')%>%mutate(jcode=as.numeric(jcode))%>%
	mutate(standard.squaredSG = standard.hexsgPCA1^2,log.stand.dist2shore = log(standard.hexdist2shore + 1))
names(hexdata)

	# the seagrass PCA data was clusterng in plots. explore:

	mdsldscomparedtoPCA.densities<-ggplot(data=hexdata)+
		geom_density(aes(x=standard.hexsgPCA1),fill='cadetblue2',alpha=0.3)+
		geom_density(aes(x=standard.hexmds),fill='violetred2',alpha=0.3)+
		geom_density(aes(x=standard.hexlds),fill='goldenrod3',alpha=0.3)+
		theme_bw()

		ggsave(mdsldscomparedtoPCA.densities,file='hypotesting_dredge_results/mdsldscomparedtoPCA_hexdata_geomDensities_oct23.png',device='png',units='in',dpi=600,height=3,width=3.5)
		# beta-binomial transformation on mds and lds and then pca? to maybe meet normality criteria? 

		



sub_sample<-hexdata%>%slice_sample(n=50)
write.csv(sub_sample,'subsample_hexdata_formodeltesting.csv')

#############################
## - Define Hypotheses
#############################

	##Define some hypotheses to test:
	## On the basis that we've established a preference for distance form the central mangrove forest and medium density seagrass...
	#	1.  Seagrasses and distance metrics: with increasing distance from jetties, medium density seagrass will increase. And with increasing depth and distance from central refuge, medium density seagrass will decrease.
		# global model for MuMIn:    mds ~ depth x distcmg x dist2jetty x dist2shore 

	#	2. Habitat metrics on sharks, using receiver point data for sharks: juvenile sharks will select for shallow habitats close to the central mangroves, but distant from man-made jetties. SHorelines will have null results because the relationship is not linear.  
		## Fish will selet for shallow habitats close to shorelines. Jetties will have no effet and the distance from the central mangroves will also be irrelevant.
		# global model:   sharks ~ depth * dist2jetty * dist2shore * distcmg 
		# global model:   fish ~ depth * dist2jetty * dist2shore * distcmg 

	#	3. Large Shark detections: predict there will be more large sharks with inceasing depth and distance from the shoreline and the central mangrove forest. Predict that distance to jetty will be inversely related to large shark presence. 




#############################
## - Hypothesis 1: seagrasses PCA
#############################

	sg <- glm(standard.hexsgPCA1 ~ standard.depth * standard.hexdist2shore * standard.hexdistcmg * standard.dist2jetty, data=hexdata) #define the global model
	
	sg.dredge <- dredge(sg, beta='sd', evaluate=TRUE, trace=FALSE, extra='R^2',m.lim=c(0,4)) # dredge from the global model
	
	sg.modelsranked.tabled <- sg.dredge %>%
	  as_tibble %>%
	  mutate(weight=round(weight,5),
	         model = 1:n()) %>%
	  mutate(Null = ifelse(df == 2, 1, NA))%>%
	  pivot_longer(cols=starts_with(c('st','Null','exp')),values_to='estimate')%>%	  
	  filter(!is.na(estimate)) %>%
	  group_by(pick(2, 3,4,5,6,7,8)) %>%
	  summarise(model = paste(name, collapse = ' + ')) %>%
	  ungroup() %>%
	  select(model, everything())%>% 
	  arrange(-weight) %>%
	  flextable()%>%
	  theme_zebra()%>%
	  set_header_labels(model = 'Model',delta = 'dAICc')%>%
	  align(align = 'left', part = 'all')%>%
	  color(color='black',part='all')%>%
	  fontsize(size = 10, part = 'all')%>%
	  autofit()

	  save_as_image(sg.modelsranked.tabled,'dredged_results_seagrassesPCAresponse_depth_dist2shore_distcmg_dist2jetty.png',webshot='webshot')


#############################
## - Hypothesis 2 : Juveniles and fish, distances only as predictors.
#############################

  ## sharks, point
	sharkP.glm <- glm(standard.shark ~ standard.depth * standard.dist2shore * standard.distcmg * standard.dist2jetty, data=pointdata) #define the global model

	sharkP.dredge <- dredge(sharkP.glm, beta='sd', evaluate=TRUE, trace=FALSE, extra='R^2',m.lim=c(0,4)) # dredge from the global model

	sharkP.modelsranked.tabled <- sharkP.dredge%>%
		as_tibble %>%
	  mutate(weight=round(weight,9), model = 1:n()) %>%
	  mutate(Null = ifelse(df == 2, 1, NA))%>%
	  pivot_longer(cols=starts_with(c('st','Null')),values_to='estimate')%>%	  
	  filter(!is.na(estimate)) %>%
	  group_by(pick(2,3,4,5,6,7,8)) %>%
	  summarise(model = paste(name, collapse = ' + ')) %>%
	  ungroup() %>%
	  select(model, everything())%>% 
	  arrange(-weight) %>%
	  flextable()%>%	  
	  theme_zebra()%>%
	  set_header_labels(model = 'Model',delta = 'dAICc')%>%
	  align(align = 'left', part = 'all')%>%
	  color(color='black',part='all')%>%
	  fontsize(size = 10, part = 'all')%>%
	  autofit()

	  save_as_image(sharkP.modelsranked.tabled,'dredged_results_SharksPOINTresponse_depth_dist2shore_distcmg_dist2jetty.png',webshot='webshot')

	summary(get.models(shark.dredge,1)[[1]]) #"best" model
	best.shark <- get.models(shark.dredge,1)[[1]]

	## try with a random effect of site 
	sharkP.glmer <- lmer(standard.shark ~ standard.depth * standard.dist2shore * standard.distcmg * standard.dist2jetty + (1|buffIDnum), data=pointdata) #define the global model

	sharkP.glmer.dredge <- dredge(sharkP.glmer, beta='sd', evaluate=TRUE, trace=FALSE, extra='R^2',m.lim=c(0,5)) # dredge from the global model

	sharkP.glmer.modelsranked.tabled <- sharkP.glmer.dredge%>%
		as_tibble %>%
	  mutate(weight=round(weight,9), model = 1:n()) %>%
	  mutate(Null = ifelse(df == 3, 1, NA))%>%
	  pivot_longer(cols=starts_with(c('st','Null')),values_to='estimate')%>%	  
	  filter(!is.na(estimate)) %>%
	  group_by(pick(2,3,4,5,6,7,8)) %>%
	  summarise(model = paste(name, collapse = ' + ')) %>%
	  ungroup() %>%
	  select(model, everything())%>% 
	  arrange(-weight) %>%
	  flextable()%>%	  
	  theme_zebra()%>%
	  set_header_labels(model = 'Model',delta = 'dAICc')%>%
	  align(align = 'left', part = 'all')%>%
	  color(color='black',part='all')%>%
	  fontsize(size = 10, part = 'all')%>%
	  autofit()

	  summary(get.models(sharkP.glmer.dredge,1)[[1]]) #"best" model
	best.shark <- get.models(sharkP.glmer.dredge,1)[[1]]
 
	  save_as_image(sharkP.glmer.modelsranked.tabled,'dredged_GLMER_1BUFFID_results_SharksPOINTresponse_depth_dist2shore_distcmg_dist2jetty.png',webshot='webshot')


  ## fish, hex
	fish.glm <- glm(standard.hexfish ~ standard.depth * standard.hexdist2shore * standard.hexdistcmg * standard.dist2jetty + standard.hexsgPCA1, data=hexdata) #define the global model
	
	fish.glm2 <- glm(standard.hexfish ~ standard.depth * standard.hexdist2shore * standard.hexdistcmg * standard.dist2jetty + standard.hexsgPCA1+standard.squaredSG, data=hexdata) #define the global model, 6/10/2023: quadratic for the seagrass.
	fish.glm3 <- glm(standard.hexfish ~ standard.depth * log(standard.hexdist2shore+1) * standard.hexdistcmg * standard.dist2jetty + standard.hexsgPCA1+standard.squaredSG, data=hexdata) #define the global model, 6/10/2023: quadratic for the seagrass, and log normal for dist2shore.


	fish.dredge <- dredge(fish.glm, beta='sd', evaluate=TRUE, trace=FALSE, extra='R^2',m.lim=c(0,4)) # dredge from the global model
	fish.dredge2 <- dredge(fish.glm2, beta='sd', evaluate=TRUE, trace=FALSE, extra='R^2',m.lim=c(0,4),subset=(standard.hexsgPCA1|!standard.squaredSG)) # dredge from the global model,quadratic for seagrass
	fish.dredge3 <- dredge(fish.glm3, beta='sd', evaluate=TRUE, trace=FALSE, extra='R^2',m.lim=c(0,4),subset=(standard.hexsgPCA1|!standard.squaredSG)) # dredge from the global model,quadratic for seagrass

	fish.modelsranked.tabled3 <- fish.dredge3%>%
		as_tibble %>%
	  mutate(weight=round(weight,9), model = 1:n()) %>%
	  mutate(Null = ifelse(df == 2, 1, NA))%>%
	  pivot_longer(cols=starts_with(c('st','Null','log(sta')),values_to='estimate')%>%	  
	 filter(!is.na(estimate)) %>%
	  group_by(pick(2,3,4,5,6,7,8)) %>%
	  summarise(model = paste(name, collapse = ' + ')) %>%
	  ungroup() %>%
	  select(model, everything())%>% 
	  dplyr::arrange(-weight) %>%
	  flextable()%>%	  
	  theme_zebra()%>%
	  set_header_labels(model = 'Model',delta = 'dAICc')%>%
	  align(align = 'left', part = 'all')%>%
	  color(color='black',part='all')%>%
	  fontsize(size = 10, part = 'all')%>%
	  autofit()

	  save_as_image(fish.modelsranked.tabled3,'hypotesting_dredge_results/dredged_results_fishHEXresponse_depth_quadraticSGpca_LOGNORMALdist2shorePLUS1_distcmg_dist2jetty.png',webshot='webshot')
	get.models(fish.dredge,1)[[1]]


#############################
## - Hypothesis 3:  Large sharks
#############################

  ## predators, point
	pred.glmer <- glm(standard.relPropPD ~ standard.depth * standard.dist2shore * standard.distcmg * standard.dist2jetty , data=pointdata) #define the global model

	pred.dredged <- dredge(pred.glmer, beta='sd', evaluate=TRUE, trace=FALSE, extra='R^2',m.lim=c(0,4)) # dredge from the global model

	pred.glmer.modelsranked.tabled <- pred.dredged%>%
		as_tibble %>%
	  mutate(weight=round(weight,9), model = 1:n()) %>%
	  mutate(Null = ifelse(df == 2, 1, NA))%>%
	  pivot_longer(cols=starts_with(c('st','Null')),values_to='estimate')%>%	  
	  filter(!is.na(estimate)) %>%
	  group_by(pick(2,3,4,5,6,7,8)) %>%
	  summarise(model = paste(name, collapse = ' + ')) %>%
	  ungroup() %>%
	  select(model, everything())%>% 
	  arrange(-weight) %>%
	  flextable()%>%	  
	  theme_zebra()%>%
	  set_header_labels(model = 'Model',delta = 'dAICc')%>%
	  align(align = 'left', part = 'all')%>%
	  color(color='black',part='all')%>%
	  fontsize(size = 10, part = 'all')%>%
	  autofit()
	  
	  save_as_image(pred.glmer.modelsranked.tabled,'dredged_results_RelPropPDresponse_depth_dist2shore_distcmg_dist2jetty.png',webshot='webshot')
	get.models(fish.dredge,1)[[1]]



#############################
## - GLM/GLMMs of 'best' models: the effect of distance metrics on predictor variables
#############################

# H1: seagrasses and distance metrics

	sgm <- glm(standard.hexsgPCA1 ~ standard.dist2jetty + standard.hexdist2shore + standard.hexdistcmg + standard.hexdist2shore*standard.hexdistcmg, data=hexdata, family=gaussian)

# H2: fish and distance metrics, updated 6/10/23

	fm <- glm(standard.hexfish ~ standard.hexdist2shore + standard.hexdistcmg + standard.hexsgPCA1 + standard.hexdist2shore*standard.hexdistcmg, data=hexdata, family=gaussian)

	fm2 <- glm(standard.hexfish ~ standard.depth + log.stand.dist2shore + standard.hexdistcmg + standard.dist2jetty + standard.hexsgPCA1 + standard.squaredSG, data=hexdata, family=gaussian)
	# depth is data defficient because it only has four groups.
	fm2 <- glm(standard.hexfish ~ standard.depth + log.stand.dist2shore + standard.hexdistcmg + standard.dist2jetty + standard.hexsgPCA1 + standard.squaredSG, data=hexdata, family=gaussian)

# H3: large sharks and distance metrics

	pm <- glm(standard.relPropPD ~ standard.depth + standard.dist2jetty + standard.distcmg + standard.depth*standard.dist2jetty, data=pointdata)
	## not a glmer because of rank defficiency, wont converge. but in the bsem it is a glmer. 

# save the models as RDS 
	saveRDS(sgm,'seagrasses_glm_hypotesting_distancemetrics_sept23.RDS')
	saveRDS(fm,'fishesmetric_glm_hypotesting_distancemetrics_sept23.RDS')
	saveRDS(fm2,'hypotesting_dredge_results/fishesmetric_glm_lognormalSHOREplus1_quadraticSG_hypotesting_distancemetrics_sept23.RDS')
	saveRDS(pm,'largesharks_glm_hypotesting_distancemetrics_sept23.RDS')

#############################
## - Summary table of GLM/GLMMs of 'best' models: the effect of distance metrics on predictor variables
#############################

# read in model RDS
	sgm <- readRDS('hypotesting_dredge_results/seagrasses_glm_hypotesting_distancemetrics_sept23.RDS')
	fm <- readRDS('hypotesting_dredge_results/fishesmetric_glm_hypotesting_distancemetrics_sept23.RDS')
	fm2 <- readRDS('hypotesting_dredge_results/fishesmetric_glm_lognormalSHOREplus1_quadraticSG_hypotesting_distancemetrics_sept23.RDS')
	pm <- readRDS('hypotesting_dredge_results/largesharks_glm_hypotesting_distancemetrics_sept23.RDS')

# make table
	models <- c(sgm, fm2, pm)
	
	sgf <- flextable::as_flextable(sgm)
	pf <- flextable::as_flextable(pm)
	ff <- flextable::as_flextable(fm2)

	s1 <- tbl_regression(sgm,exp=FALSE,conf_level=0.95,label=list(standard.hexdist2shore='Dist. to Shore',standard.hexdistcmg='Dist. to Central Mangroves',standard.dist2jetty='Dist. to Jetty','standard.hexdist2shore*standard.hexdistcmg'='Dist. to Shore * Dist. to Central Mangroves'))%>%
		bold_p(t=0.05)
	f1 <- tbl_regression(fm2, exp=FALSE,conf_level=0.95,label=list(standard.depth='Depth',standard.dist2jetty='Dist. to Jetty','log.stand.dist2shore'   ='Dist. to Shore, (log+1)',standard.hexdistcmg='Dist. to Central Mangroves',standard.hexsgPCA1='Seagrass PCA',standard.squaredSG = 'Squared SG'))%>%
		bold_p(t=0.05)	
	p1 <- tbl_regression(pm, exp=FALSE,conf_level=0.95, label=list(standard.depth='Depth',standard.dist2jetty='Dist. to Jetty',standard.distcmg='Dist. to Central Mangroves','standard.depth*standard.dist2jetty'='Depth * Dist. to Jetty'))%>%
		bold_p(t=0.05)

	# side by side
	tbl_merge(tbls = list(s1,f1), tab_spanner = c('Hypothesis 1', 'Hypothesis 2'))
	# stacked 
	stacked <- tbl_stack(list(s1,f1,p1),group_header=c('1','2','3'))%>%
	bold_levels()
	
	show_header_names(stacked)
 
	responses <- as_tibble(x=c('Seagrasses','Fish','Large Sharks'))%>%rename(Respones=value)

	stacked.summary<-stacked%>%
		modify_header(groupname_col = '**Hypothesis**',label='**Predictor**')%>%
		as_flex_table()%>%	
		autofit()

	save_as_image(stacked.summary,'hypotesting_dredge_results/stackedsummarytables_hypotheses_chp3_seagrasses_fishDredge3_largesharks_glms_oct2023.png' ,webshot='webshot')

	## check resiudals of models

		# sgm
			par(mfrow=c(2,2))
			qq.sgm<-plot(sgm)
		# fm2
			par(mfrow=c(2,2))
			qq.fm2<-plot(fm2)
		# pm
			par(mfrow=c(2,2))
			qq.pm<-plot(pm)


#############################
## - Prediction plots GLM/GLMMs of 'best' models: the effect of distance metrics on predictor variables
#############################

# read in model RDS
	sgm <- readRDS('hypotesting_dredge_results/seagrasses_glm_hypotesting_distancemetrics_sept23.RDS')
	fm <- readRDS('hypotesting_dredge_results/fishesmetric_glm_hypotesting_distancemetrics_sept23.RDS')
	fm2 <- readRDS('hypotesting_dredge_results/fishesmetric_glm_lognormalSHOREplus1_quadraticSG_hypotesting_distancemetrics_sept23.RDS')
	pm <- readRDS('hypotesting_dredge_results/largesharks_pressure_glm_hypotesting_distancemetrics_dec23.RDS')

# find min and max of preditors and store in a df 
	 	hexlims <- as.data.frame(matrix(ncol=3,nrow=5))
	 		colnames(hexlims)=c('var','min','max')
		 	hexlims$var=c('standard.dist2jetty','standard.hexdist2shore','standard.hexdistcmg','standard.depth','standard.hexsgPCA1')
		 	hexlims$min=c(min(hexdata$standard.dist2jetty),min(hexdata$standard.hexdist2shore),min(hexdata$standard.hexdistcmg),min(hexdata$standard.depth),min(hexdata$standard.hexsgPCA1))
		 	hexlims$max=c(max(hexdata$standard.dist2jetty),max(hexdata$standard.hexdist2shore),max(hexdata$standard.hexdistcmg),max(hexdata$standard.depth),max(hexdata$standard.hexsgPCA1))
	 	hexlims

	 	pointlims <- as.data.frame(matrix(ncol=3,nrow=5))
	 		colnames(pointlims)=c('var','min','max')
		 	pointlims$var=c('standard.dist2shore','standard.distcmg','standard.depth','standard.dist2jetty','sgPCA1')
		 	pointlims$min=c(min(pointdata$standard.dist2shore),min(pointdata$standard.distcmg),min(pointdata$standard.depth),min(pointdata$standard.dist2jetty),min(pointdata$sgPCA1))
		 	pointlims$max=c(max(pointdata$standard.dist2shore),max(pointdata$standard.distcmg),max(pointdata$standard.depth),max(pointdata$standard.dist2jetty),max(pointdata$sgPCA1))
	 	pointlims

## used ggpredict to calculate estimate values for responses when a focal predictor is variable and non-focal are held constant.
	## colours for graphs 
	 	# seagrasses
	 	 '#653515' # points
	 	 '#6d2e35' # lines
	 	 '#3e2736' # ribbons

	 	# large sharks 
	 	 '#d13474' # points
	 	 '#a22558' # lines
	 	 '#e280a8' # ribbons
	 	
	 	# fishes 
	 	 '#836304' # points
	 	 '#fad361' # ribbons
	 	 '#e3ab07' # lines
	
	##################
	## seagrasses: dist2jetty, dist2shore, distcmg, dst2shore*distcmg
		
		sgm.hex.predicts.dist2jetty<-as.data.frame(ggpredict(sgm,terms='standard.dist2jetty[all]',type='fixed'))
		sgm.plot1<- ggplot()+
			geom_point(data=hexdata, aes(x=standard.dist2jetty, y=standard.hexsgPCA1),pch=21,cex=0.5,col='#3a1d20')+
			ylab('Seagrass PCA (standardised)')+xlab('Distance from Jetty (standardised)')+
			geom_line(data=sgm.hex.predicts.dist2jetty,aes(x=x,y=predicted),col='#653515',size=.6,linetype=1)+
			geom_ribbon(data=sgm.hex.predicts.dist2jetty,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#6d2e35')+
			theme_bw()

		sgm.hex.predicts.dist2shore<-as.data.frame(ggpredict(sgm,terms='standard.hexdist2shore[all]',type='fixed'))
		sgm.plot2<- ggplot()+
			geom_point(data=hexdata, aes(x=standard.hexdist2shore, y=standard.hexsgPCA1),pch=21,cex=0.5,col='#3a1d20')+
			geom_line(data=sgm.hex.predicts.dist2shore,aes(x=x,y=predicted),col='#653515',size=.6,linetype=1)+
			geom_ribbon(data=sgm.hex.predicts.dist2shore,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#6d2e35')+
			ylab('Seagrass PCA (standardised)')+xlab('Distance from Shore (standardised)')+
			theme_bw()

		sgm.hex.predicts.dist2cmg<-as.data.frame(ggpredict(sgm,terms='standard.hexdistcmg[all]',type='fixed'))
		sgm.plot3<- ggplot()+
			geom_point(data=hexdata, aes(x=standard.hexdistcmg, y=standard.hexsgPCA1,xmin=min(hexdata$standard.hexdistcmg),xmax=max(standard.hexdistcmg)),pch=21,cex=0.5,col='#3a1d20')+
			ylab('Seagrass PCA (standardised)')+xlab('Distance from CM (standardised)')+
			geom_line(data=sgm.hex.predicts.dist2cmg,aes(x=x,y=predicted),col='#653515',size=.6,linetype=1)+
			geom_ribbon(data=sgm.hex.predicts.dist2cmg,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#6d2e35')+
			theme_bw()

		## you have to calculate and plot the interaction terms differently. ggpredict puts the first term on the x-axis and then shows the second term at representative values - here three: the mean, and +- a SD unit 
		sgm.hex.predicts.intx.dist2.shore.cmg<-as.data.frame(ggpredict(sgm,terms=c("standard.hexdist2shore", "standard.hexdistcmg"),type='fixed'))
		
		intx.col.sgm <- c('#653515','#6d2e35','#2a4858')
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

		ggsave(sgm.prediction.plots.formatted,file='hypotesting_dredge_results/seagrassesGLM_fromDredge_fourplots_predictions.png',device='png',units='in',dpi=350,height=6,width=7)


	##################
	## fishes
		
		# fishes, f1: glm(standard.hexfish ~ standard.depth * standard.hexdist2shore * standard.hexdistcmg * standard.dist2jetty + standard.hexsgPCA1, data=hexdata) 
			fm.hex.predicts.sgPCA<-as.data.frame(ggpredict(fm,terms='standard.hexsgPCA1[all]',type='fixed'))
			fm.plot1<- ggplot()+
				geom_point(data=hexdata, aes(x=standard.hexsgPCA1, y=standard.hexfish),pch=21,cex=0.5,col='#836304')+
				ylab('Fishes (standardised)')+xlab('Seagrasses PCA (standardised)')+
				geom_line(data=fm.hex.predicts.dist2jetty,aes(x=x,y=predicted),col='#e3ab07',size=.6,linetype=1)+
				geom_ribbon(data=fm.hex.predicts.sgPCA,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#fad361')+
				theme_bw()

			fm.hex.predicts.distcmg<-as.data.frame(ggpredict(fm,terms='standard.hexdistcmg[all]',type='fixed'))
			fm.plot2<- ggplot()+
				geom_point(data=hexdata, aes(x=standard.hexdistcmg, y=standard.hexfish,xmin=min(hexdata$standard.hexdistcmg),xmax=max(standard.hexdistcmg)),pch=21,cex=0.5,col='#836304')+
				ylab('Fishes (standardised)')+xlab('Distance from CM (standardised)')+
				geom_line(data=fm.hex.predicts.distcmg,aes(x=x,y=predicted),col='#e3ab07',size=.6,linetype=1)+
				geom_ribbon(data=fm.hex.predicts.distcmg,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#fad361')+
				theme_bw()

			fm.hex.predicts.dist2shore<-as.data.frame(ggpredict(fm,terms='standard.hexdist2shore[all]',type='fixed'))
			fm.plot3<- ggplot()+
				geom_point(data=hexdata, aes(x=standard.hexdist2shore, y=standard.hexfish,xmin=min(hexdata$standard.hexdist2shore),xmax=max(standard.hexdist2shore)),pch=21,cex=0.5,col='#836304')+
				ylab('Fishes (standardised)')+xlab('Distance from Shore (standardised)')+
				geom_line(data=fm.hex.predicts.dist2shore,aes(x=x,y=predicted),col='#e3ab07',size=.6,linetype=1)+
				geom_ribbon(data=fm.hex.predicts.dist2shore,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#fad361')+
				theme_bw()

			## you have to calculate and plot the interaction terms differently. ggpredict puts the first term on the x-axis and then shows the second term at representative values - here three: the mean, and +- a SD unit 
			fm.hex.predicts.intx.dist2.shore.cmg<-as.data.frame(ggpredict(fm,terms=c("standard.hexdist2shore", "standard.hexdistcmg"),type='fixed'))
			
			intx.col.fm <- c('#9e8d59','#e3ab07','#aa4827')
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

			ggsave(fm.prediction.plots.formatted,file='hypotesting_dredge_results/fishmetricGLM_fromDredge_fourplots_predictions.png',device='png',units='in',dpi=350,height=6,width=7)

		# fishes, f2: lognormal dist2shore + 1, quadratic for seagrsses PCA 
			fm2.hex.predicts.sgPCA<-as.data.frame(ggpredict(fm2,terms=c('standard.hexsgPCA1[all]','standard.squaredSG'),type='fixed'))
			fm2.hex.predictsAVERAGE.sgPCA<-as.data.frame(ggeffect(fm2,terms=c('standard.hexsgPCA1[all]','standard.squaredSG'),type='fixed'))
			fm2.plot1<- ggplot()+
				geom_point(data=hexdata, aes(x=standard.hexsgPCA1, y=standard.hexfish),pch=21,cex=0.5,col='#836304')+
				ylab('Fishes (standardised)')+xlab('Seagrasses PCA (quadratic,standardised)')+
				geom_line(data=fm2.hex.predicts.sgPCA,aes(x=x,y=predicted,group=group),col='#e3ab07',size=.6,linetype=1)+
				geom_ribbon(data=fm2.hex.predicts.sgPCA,aes(x=x,ymin=conf.low,ymax=conf.high,group=group),alpha=0.3, fill='#fad361')+
				theme_bw()


			fm2.hex.predicts.distcmg<-as.data.frame(ggpredict(fm2,terms='standard.hexdistcmg[all]',type='fixed'))
			fm2.plot2<- ggplot()+
				geom_line(data=fm2.hex.predicts.distcmg,aes(x=x,y=predicted),col='#e3ab07',size=.6,linetype=1)+
				geom_ribbon(data=fm2.hex.predicts.distcmg,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#fad361')+
				geom_point(data=hexdata, aes(x=standard.hexdistcmg, y=standard.hexfish),pch=21,cex=0.5,col='#836304')+
				ylab('Fishes (standardised)')+xlab('Distance from CM (standardised)')+
				theme_bw()

			fm2.hex.predicts.dist2shore<-as.data.frame(ggpredict(fm2,terms='log.stand.dist2shore[all]',type='fixed'))
			fm2.plot3<- ggplot()+
				geom_line(data=fm2.hex.predicts.dist2shore,aes(x=x,y=predicted),col='#e3ab07',size=.6,linetype=1)+
				geom_ribbon(data=fm2.hex.predicts.dist2shore,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#fad361')+
				geom_point(data=hexdata, aes(x=log.stand.dist2shore, y=standard.hexfish),pch=21,cex=0.5,col='#836304')+
				ylab('Fishes (standardised)')+xlab('Distance from Shore (log, standardised)')+
				theme_bw()

			fm2.hex.predicts.jetty<-as.data.frame(ggpredict(fm2,terms='standard.dist2jetty[all]',type='fixed'))
			fm2.plot4<-ggplot()+
				geom_line(data=fm2.hex.predicts.jetty,aes(x=x,y=predicted),col='#e3ab07',size=.6,linetype=1)+
				geom_ribbon(data=fm2.hex.predicts.jetty,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#fad361')+
				geom_point(data=hexdata, aes(x=standard.dist2jetty, y=standard.hexfish),pch=21,cex=0.5,col='#836304')+
				ylab('Fishes (standardised)')+xlab('Dist. to Jetty (standardised)')+
				theme_bw()

			fm2.hex.predicts.depth<-as.data.frame(ggpredict(fm2,terms='standard.depth[all]',type='fixed'))
			fm2.plot5<-ggplot()+
				geom_point(data=hexdata, aes(x=standard.dist2jetty, y=standard.hexfish),pch=21,cex=0.5,col='#836304')+
				geom_line(data=fm2.hex.predicts.depth,aes(x=x,y=predicted),col='#e3ab07',size=.3,linetype=1)+
				ylab('Fishes (standardised)')+xlab('Depth (standardised)')+
				theme_bw()


			# put them together in a box
			fm2.prediction.plots.formatted <- (fm2.plot1 + fm2.plot2)/(fm2.plot3 + fm2.plot4)

			ggsave(fm2.prediction.plots.formatted,file='hypotesting_dredge_results/fishmetricGLM_fromDredge_logDist2shore_quadraticSG_fourplots_predictions_oct23.png',device='png',units='in',dpi=350,height=6,width=7)



	##################
	## large sharks: depth, dist2jetty, distcmg, depth*dist2jetty

	pm..predicts.dist2jetty<-as.data.frame(ggpredict(pm,terms='standard.dist2jetty[all]',type='fixed'))
	pm.plot1<- ggplot()+
		geom_point(data=pointdata, aes(x=standard.dist2jetty, y=standard.relPropPD),pch=21,col='#d13474')+
		ylab('Large Sharks (standardised)')+xlab('Distance from Jetty (standardised)')+
		geom_line(data=pm..predicts.dist2jetty,aes(x=x,y=predicted),col='#a22558',size=.6,linetype=1)+
		geom_ribbon(data=pm..predicts.dist2jetty,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#e280a8')+
		theme_bw()

	pm..predicts.distcmg<-as.data.frame(ggpredict(pm,terms='standard.distcmg[all]',type='fixed'))
	pm.plot2<- ggplot()+
		geom_point(data=pointdata, aes(x=standard.distcmg, y=standard.relPropPD,xmin=min(pointdata$standard.distcmg),xmax=max(standard.distcmg,ymin=min(pointdata$standard.relPropPD),ymax=max(pointdata$standard.relPropPD))),pch=21,col='#d13474')+
		ylab('Large Sharks (standardised)')+xlab('Distance from CM (standardised)')+
		geom_line(data=pm..predicts.distcmg,aes(x=x,y=predicted),col='#a22558',size=.6,linetype=1)+
		geom_ribbon(data=pm..predicts.distcmg,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#e280a8')+
		theme_bw()

	pm..predicts.depth<-as.data.frame(ggpredict(pm,terms='standard.depth[all]',type='fixed'))
	pm.plot3<- ggplot()+
		geom_point(data=pointdata, aes(x=standard.depth, y=standard.relPropPD,xmin=min(pointdata$standard.depth),xmax=max(standard.depth)),pch=21,col='#d13474')+
		ylab('Large Sharks (standardised)')+xlab('Depth (standardised)')+
		geom_line(data=pm..predicts.depth,aes(x=x,y=predicted),col='#a22558',size=.6,linetype=1)+
		geom_ribbon(data=pm..predicts.depth,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#e280a8')+
		theme_bw()

	## you have to calculate and plot the interaction terms differently. ggpredict puts the first term on the x-axis and then shows the second term at representative values - here three: the mean, and +- a SD unit 
	
	pm..predicts.intx.dist2.shore.cmg<-as.data.frame(ggpredict(pm,terms=c("standard.depth", "standard.dist2jetty"),type='fixed'))
	
	intx.col.pm <- c('#e280a8','#a22558','#310b1a')
	pm.plot4<-
	ggplot()+
		geom_line(data=pm..predicts.intx.dist2.shore.cmg,aes(x=x,y=predicted, group=group, col=group),size=.6,linetype=1)+
		scale_color_manual(values=intx.col.pm,name='Levels')+
		geom_ribbon(data=pm..predicts.intx.dist2.shore.cmg,aes(group=group,x=x,ymin=conf.low,ymax=conf.high,fill=group),alpha=0.3)+		
		scale_fill_manual(values=intx.col.pm,name='Levels')+
		ylab('Large Sharks (standardised)')+
		xlab('Depth (standardised)')+
		theme_bw()+
		theme(legend.position='bottom')


	# put them together in a box
	pm.prediction.plots.formatted <- (pm.plot1 + pm.plot2)/(pm.plot3 + pm.plot4)

	ggsave(pm.prediction.plots.formatted,file='hypotesting_dredge_results/largesharks_metricGLM_fromDredge_fourplots_predictions.png',device='png',units='in',dpi=350,height=6,width=7)



	## alternative formatting - 3 columns, four rows

	all.in.one.3col.4rows <- (sgm.plot1 / sgm.plot2 / sgm.plot3 / sgm.plot4) |	(fm.plot1 / fm.plot2 / fm.plot3 / fm.plot4) | 	(pm.plot1 / pm.plot2 / pm.plot3 / pm.plot4) 

	ggsave(all.in.one.3col.4rows,file='hypotesting_dredge_results/diffformatting_allin1_fishSGlargesharks_fromDredge_predictions.png',device='png',units='in',dpi=500,height=10,width=8.5)



grid <- st_as_sf(st_read('winter2020habitat_hexagon_grid_NOland_ECdata.shp'),crs='WGS84')
land <- st_as_sf(st_read('bim_onlyland.kml'),crs='WGS84')
grid.plot<-ggplot()+geom_sf(data=st_difference(grid,st_union(land)),fill='white',col='grey32',pch=.8)+theme_void()
grid.plot<-ggplot()+geom_sf(data=grid,fill='white',col='grey32',pch=.8)+theme_void()
grid.plot
ggsave(grid.plot, file='griddedhabitatdata_sept23_EChabdata_noland.png',device='png',units='in',dpi=600,height=6,width=4)



###########################
###########################
## same plots, all navy 

	navy: #3a6c74
	
	##################
	## seagrasses: dist2jetty, dist2shore, distcmg, dst2shore*distcmg
		
		sgm.hex.predicts.dist2jetty<-as.data.frame(ggpredict(sgm,terms='standard.dist2jetty[all]',type='fixed'))
		sgm.plot1<- ggplot()+
			geom_point(data=hexdata, aes(x=standard.dist2jetty, y=standard.hexsgPCA1),pch=21,cex=0.5,col='#3a6c74')+
			ylab('Seagrass PCA (standardised)')+xlab('Distance from Jetty (standardised)')+
			geom_line(data=sgm.hex.predicts.dist2jetty,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
			geom_ribbon(data=sgm.hex.predicts.dist2jetty,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
			theme_bw()

		sgm.hex.predicts.dist2shore<-as.data.frame(ggpredict(sgm,terms='standard.hexdist2shore[all]',type='fixed'))
		sgm.plot2<- ggplot()+
			geom_point(data=hexdata, aes(x=standard.hexdist2shore, y=standard.hexsgPCA1,xmin=min(hexdata$standard.hexdist2shore),xmax=max(standard.hexdist2shore)),pch=21,cex=0.5,col='#3a6c74')+
			ylab('Seagrass PCA (standardised)')+xlab('Distance from Shore (standardised)')+
			geom_line(data=sgm.hex.predicts.dist2shore,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
			geom_ribbon(data=sgm.hex.predicts.dist2shore,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
			theme_bw()

		sgm.hex.predicts.dist2cmg<-as.data.frame(ggpredict(sgm,terms='standard.hexdistcmg[all]',type='fixed'))
		sgm.plot3<- ggplot()+
			geom_point(data=hexdata, aes(x=standard.hexdistcmg, y=standard.hexsgPCA1,xmin=min(hexdata$standard.hexdistcmg),xmax=max(standard.hexdistcmg)),pch=21,cex=0.5,col='#3a6c74')+
			ylab('Seagrass PCA (standardised)')+xlab('Distance from CM (standardised)')+
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


	
		
		# fishes, f1: glm(standard.hexfish ~ standard.depth * standard.hexdist2shore * standard.hexdistcmg * standard.dist2jetty + standard.hexsgPCA1, data=hexdata) 
			fm.hex.predicts.sgPCA<-as.data.frame(ggpredict(fm,terms='standard.hexsgPCA1[all]',type='fixed'))
			fm.plot1<- ggplot()+
				geom_point(data=hexdata, aes(x=standard.hexsgPCA1, y=standard.hexfish),pch=21,cex=0.5,col='#3a6c74')+
				ylab('Fishes (standardised)')+xlab('Seagrasses PCA (standardised)')+
				geom_line(data=fm.hex.predicts.dist2jetty,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
				geom_ribbon(data=fm.hex.predicts.sgPCA,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
				theme_bw()

			fm.hex.predicts.distcmg<-as.data.frame(ggpredict(fm,terms='standard.hexdistcmg[all]',type='fixed'))
			fm.plot2<- ggplot()+
				geom_point(data=hexdata, aes(x=standard.hexdistcmg, y=standard.hexfish,xmin=min(hexdata$standard.hexdistcmg),xmax=max(standard.hexdistcmg)),pch=21,cex=0.5,col='#3a6c74')+
				ylab('Fishes (standardised)')+xlab('Distance from CM (standardised)')+
				geom_line(data=fm.hex.predicts.distcmg,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
				geom_ribbon(data=fm.hex.predicts.distcmg,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
				theme_bw()

			fm.hex.predicts.dist2shore<-as.data.frame(ggpredict(fm,terms='standard.hexdist2shore[all]',type='fixed'))
			fm.plot3<- ggplot()+
				geom_point(data=hexdata, aes(x=standard.hexdist2shore, y=standard.hexfish,xmin=min(hexdata$standard.hexdist2shore),xmax=max(standard.hexdist2shore)),pch=21,cex=0.5,col='#3a6c74')+
				ylab('Fishes (standardised)')+xlab('Distance from Shore (standardised)')+
				geom_line(data=fm.hex.predicts.dist2shore,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
				geom_ribbon(data=fm.hex.predicts.dist2shore,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
				theme_bw()

			## you have to calculate and plot the interaction terms differently. ggpredict puts the first term on the x-axis and then shows the second term at representative values - here three: the mean, and +- a SD unit 
			fm.hex.predicts.intx.dist2.shore.cmg<-as.data.frame(ggpredict(fm,terms=c("standard.hexdist2shore", "standard.hexdistcmg"),type='fixed'))
			
		intx.col.sgm <- c('#3a6c74','#708d8e','#3cbcfc')
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

			ggsave(fm.prediction.plots.formatted,file='hypotesting_dredge_results/fishmetricGLM_fromDredge_fourplots_predictions.png',device='png',units='in',dpi=350,height=6,width=7)

	##################
	## large sharks: depth, dist2jetty, distcmg, depth*dist2jetty
	pm..predicts.dist2jetty<-as.data.frame(ggpredict(pm,terms='standard.dist2jetty[all]',type='fixed'))
	pm.plot1<- ggplot()+
		geom_point(data=pointdata, aes(x=standard.dist2jetty, y=standard.press),pch=21,col='#3a6c74')+
		ylab('Large Sharks (standardised)')+xlab('Distance from Jetty (standardised)')+
		geom_line(data=pm..predicts.dist2jetty,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
		geom_ribbon(data=pm..predicts.dist2jetty,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
		theme_bw()

	pm..predicts.distcmg<-as.data.frame(ggpredict(pm,terms='standard.distcmg[all]',type='fixed'))
	pm.plot2<- ggplot()+
		geom_point(data=pointdata, aes(x=standard.distcmg, y=standard.press,xmin=min(pointdata$standard.distcmg),xmax=max(standard.distcmg,ymin=min(pointdata$standard.press),ymax=max(pointdata$standard.press))),pch=21,col='#3a6c74')+
		ylab('Large Sharks (standardised)')+xlab('Distance from CM (standardised)')+
		geom_line(data=pm..predicts.distcmg,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
		geom_ribbon(data=pm..predicts.distcmg,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
		theme_bw()

	pm..predicts.depth<-as.data.frame(ggpredict(pm,terms='standard.depth[all]',type='fixed'))
	pm.plot3<- ggplot()+
		geom_point(data=pointdata, aes(x=standard.depth, y=standard.press,xmin=min(pointdata$standard.depth),xmax=max(standard.depth)),pch=21,col='#3a6c74')+
		ylab('Large Sharks (standardised)')+xlab('Depth (standardised)')+
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
		ylab('Large Sharks (standardised)')+
		xlab('Depth (standardised)')+
		theme_bw()+
		theme(legend.position='bottom')



	# put them together in a box
	pm.prediction.plots.formatted <- (pm.plot1 + pm.plot2)/(pm.plot3 + pm.plot4)

	ggsave(pm.prediction.plots.formatted,file='hypotesting_dredge_results/seescapescolortheme_largesharks_metricGLM_fromDredge_fourplots_predictions.png',device='png',units='in',dpi=350,height=6,width=7)



	## alternative formatting - 3 columns, four rows

	all.in.one.3col.4rows <- (sgm.plot1 / sgm.plot2 / sgm.plot3 / sgm.plot4 ) |	(fm2.plot1 / fm2.plot2 / fm2.plot3 / fm2.plot4 / fm2.plot5) | 	(pm.plot1 / pm.plot2 / pm.plot3 / pm.plot4) 

	ggsave(all.in.one.3col.4rows,file='hypotesting_dredge_results/seescapescolortheme_diffformatting_allin1_fishSGlargesharks_fromDredge_predictions.png',device='png',units='in',dpi=500,height=10,width=8.5)




































