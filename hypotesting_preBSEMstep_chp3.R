## Hypothesis testing for chapter 3 (resource and risk BSEM) 

## inspired by methods in Bennett et al. Ecology & Evolution

## created 18 September 2023 by Molly M Kressler 

# Three hypotheses tested with the 'dredge' function in the package 'MuMIn' (Multi-Model inference)
https://rdrr.io/cran/MuMIn/man/dredge.html 


#############################
## - Load Workspace and Data
#############################

pacman::p_load(tidyverse,MuMIn,ggplot2,flextable,cowplot,patchwork,lme4,stats,ggeffects,gtsummary)
  options(na.action = "na.fail")

setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd/resource_chp3')

pointdata<-read.csv('standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_may24.csv')

hexdata<-read.csv('data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_may24.csv')%>%
		mutate(jcode=as.numeric(jcode))%>%
		rename(standard.hexsgPCA1 = standard.sgPCA1)

names(hexdata)

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
	
	sg.dredged.models <- get.models(sg.dredge,subset=TRUE)
	summary(get.models(sg.dredge,1)[[1]]) #"best" model
	best.sg <- get.models(sg.dredge,1)[[1]]

	sg.modelsranked.tabled <- sg.dredge %>%
	  as_tibble %>%
	  round(3)%>%
	  mutate(weight=round(weight,3),
	         model = 1:n()) %>%
	  mutate(Null = ifelse(df == 2, 1, NA))%>%
	  pivot_longer(cols=starts_with(c('st','Null')),values_to='estimate')%>%	  
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

  ## sharks, point - updated 17/11/2023
	sharkP.glm.abiotics <- glm(standard.shark ~ standard.depth * standard.dist2shore * standard.distcmg * standard.dist2jetty, data=pointdata) #define the global model
	sharkP.glm.fishplusabiotics <- glm(standard.shark ~ standard.fish + standard.depth * standard.dist2shore * standard.distcmg * standard.dist2jetty, data=pointdata) #define the global model
	sharkP.glm.pressplusabiotics <- glm(standard.shark ~ standard.press + standard.depth * standard.dist2shore * standard.distcmg * standard.dist2jetty, data=pointdata) #define the global model

	sharkP.dredge.abiotics <- dredge(sharkP.glm.abiotics, beta='sd', evaluate=TRUE, trace=FALSE, extra='R^2',m.lim=c(0,4)) # dredge from the global model
	sharkP.dredge.fishplusabiotics <- dredge(sharkP.glm.fishplusabiotics, beta='sd', evaluate=TRUE, trace=FALSE, extra='R^2',m.lim=c(0,4)) # dredge from the global model
	sharkP.dredge.pressplusabiotics <- dredge(sharkP.glm.pressplusabiotics, beta='sd', evaluate=TRUE, trace=FALSE, extra='R^2',m.lim=c(0,4)) # dredge from the global model

	abiotics.sharkP.modelsranked.tabled <- sharkP.dredge.abiotics%>%
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

	fishplusbiotics.sharkP.modelsranked.tabled <- sharkP.dredge.fishplusabiotics%>%
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

	pressplusbiotics.sharkP.modelsranked.tabled <- sharkP.dredge.pressplusabiotics%>%
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

	save_as_image(abiotics.sharkP.modelsranked.tabled,'hypotesting_dredge_results/dredged_results_abiotics.sharkP.modelsranked.png',webshot='webshot')
	save_as_image(fishplusbiotics.sharkP.modelsranked.tabled,'hypotesting_dredge_results/dredged_results_fishplusbiotics.sharkP.modelsranked.png',webshot='webshot')
	save_as_image(pressplusbiotics.sharkP.modelsranked.tabled,'hypotesting_dredge_results/dredged_results_pressplusbiotics.sharkP.modelsranked.png',webshot='webshot')


# make table
	a1 <- glm(standard.shark ~ standard.depth + standard.dist2shore + standard.distcmg * standard.dist2jetty, data=pointdata)
	a2 <- get.models(sharkP.dredge.fishplusabiotics,1)[[1]]
	a3 <- get.models(sharkP.dredge.pressplusabiotics,1)[[1]]
	models <- c(sgm, fm, pm)
	
	a1f <- flextable::as_flextable(a1)
	a2f <- flextable::as_flextable(a2)
	a3f <- flextable::as_flextable(a3)

	a1t <- tbl_regression(a1,exp=FALSE,conf_level=0.95,label=list(standard.dist2shore='Dist. to Shore',
		standard.dist2jetty='Dist. to Jetty',
		standard.distcmg='Dist. to Mangrove',
		standard.depth='Depth',
		'standard.distcmg*standard.dist2jetty'='Dist. to Jetty * Dist. to Central Mangroves'))%>%
		bold_p(t=0.05)

	a2t <- tbl_regression(a2, exp=FALSE,conf_level=0.95,label=list(standard.dist2shore='Dist. to Shore',standard.distcmg='Dist. to Central Mangroves',standard.hexsgPCA1='Seagrass PCA','standard.dist2shore*standard.distcmg'='Dist. to Shore * Dist. to Central Mangroves'))%>%
		bold_p(t=0.05)	
	a3t <- tbl_regression(a3, exp=FALSE,conf_level=0.95, label=list(standard.depth='Depth',standard.dist2jetty='Dist. to Jetty',standard.distcmg='Dist. to Central Mangroves','standard.depth*standard.dist2jetty'='Depth * Dist. to Jetty'))%>%
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


	save_as_image(stacked.summary,'hypotesting_dredge_results/stackedsummarytables_hypotheses_chp3_seagrasses_fish_largesharks_glms_sept2023.png' ,webshot='webshot')




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

	fish.dredge <- dredge(fish.glm, beta='sd', evaluate=TRUE, trace=FALSE, extra='R^2',m.lim=c(0,4)) # dredge from the global model

	fish.modelsranked.tabled <- fish.dredge%>%
		as_tibble %>%
		round(3)%>%
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

	  save_as_image(fish.modelsranked.tabled,'dredged_results_fishHEXresponse_depth__MDS_LDS_dist2shore_distcmg_dist2jetty.png',webshot='webshot')
	get.models(fish.dredge,1)[[1]]


#############################
## - Hypothesis 3:  Large sharks
#############################

  ## predators, point
	pred.glm <- glm(zlogit.sqzrisk ~ standard.depth * standard.dist2shore * standard.distcmg * standard.dist2jetty, data=pointdata) #define the global model
	
	pred.dredged <- dredge(pred.glm, beta='sd', evaluate=TRUE, trace=FALSE, extra='R^2',m.lim=c(0,5)) # dredge from the global model

	press.glm.modelsranked.tabled <- pred.dredged%>%
		as_tibble %>%
		round(3) %>%
	  mutate(weight=round(weight,9), model = 1:n()) %>%
	  mutate(Null = ifelse(df == 2, 1, NA))%>%
	  pivot_longer(cols=starts_with(c('st','Null')),values_to='estimate')%>%	  
	  filter(!is.na(estimate)) %>%
	  group_by(pick(2,3,4,5,6,7,8)) %>%
	  summarise(model = paste(name, collapse = ' + ')) %>%
	  ungroup() %>%
	  dplyr::select(model, everything())%>% 
	  arrange(-weight) %>%
	  flextable()%>%	  
	  theme_zebra()%>%
	  set_header_labels(model = 'Model',delta = 'dAICc')%>%
	  align(align = 'left', part = 'all')%>%
	  color(color='black',part='all')%>%
	  fontsize(size = 10, part = 'all')%>%
	  autofit()
	press.glm.modelsranked.tabled
	  
	  save_as_image(press.glm.modelsranked.tabled,'resource_chp3/hypotesting_dredge_results/dredged_results_RelPropPDresponse_depth_dist2shore_distcmg_dist2jetty.png',webshot='webshot')
	
		pm1<-get.models(pred.dredged,1)[[1]]
		pm2<-get.models(pred.dredged,2)[[1]]
		pm3 <-glm(standard.press ~ standard.depth + standard.dist2jetty + standard.distcmg + standard.depth*standard.dist2jetty + standard.depth*standard.distcmg, data=pointdata) 
		summary(pm3)

		with(summary(pm3), 1 - deviance/null.deviance) R2 = 0.139

#############################
## - GLM/GLMMs of 'best' models: the effect of distance metrics on predictor variables
#############################

# H1: seagrasses and distance metrics

	sgm <- glm(standard.sgPCA1 ~ standard.dist2jetty + standard.hexdist2shore + standard.hexdistcmg + standard.hexdist2shore*standard.hexdistcmg, data=hexdata, family=gaussian)

# H2: fish and distance metrics

	fm <- glm(standard.hexfish ~ standard.hexdist2shore + standard.hexdistcmg + standard.sgPCA1 + standard.hexdist2shore*standard.hexdistcmg, data=hexdata, family=gaussian)

# H3: large sharks and distance metrics

	pm <- 	get.models(pred.dredged,1)[[1]]
	pmm<-glm(standard.press ~ standard.depth + standard.dist2jetty + standard.distcmg + standard.depth*standard.dist2jetty, data=pointdata)

# save the models as RDS 
	saveRDS(sgm,'seagrasses_glm_hypotesting_distancemetrics_sept23.RDS')
	saveRDS(fm,'fishesmetric_glm_hypotesting_distancemetrics_sept23.RDS')
	saveRDS(pm,'hypotesting_dredge_results/largesharks_relativerisk_glm_hypotesting_distancemetrics_may24.RDS')

#############################
## - Summary table of GLM/GLMMs of 'best' models: the effect of distance metrics on predictor variables
#############################

# read in model RDS
	sgm <- readRDS('hypotesting_dredge_results/seagrasses_glm_hypotesting_distancemetrics_sept23.RDS')
	fm <- readRDS('hypotesting_dredge_results/fishesmetric_glm_hypotesting_distancemetrics_sept23.RDS')
	pm <- readRDS('hypotesting_dredge_results/largesharks_relativerisk_glm_hypotesting_distancemetrics_may24.RDS')

# make table
	models <- c(sgm, fm, pm)
	
	sgf <- flextable::as_flextable(sgm)
	pf <- flextable::as_flextable(pm)
	ff <- flextable::as_flextable(fm)

	s1 <- tbl_regression(sgm,exp=FALSE,conf_level=0.95,label=list(standard.hexdist2shore='Dist. to Shore',standard.hexdistcmg='Dist. to Central Mangroves',standard.dist2jetty='Dist. to Jetty','standard.hexdist2shore*standard.hexdistcmg'='Dist. to Shore * Dist. to Central Mangroves'))%>%
		bold_p(t=0.05)
	f1 <- tbl_regression(fm, exp=FALSE,conf_level=0.95,label=list(standard.hexdist2shore='Dist. to Shore',standard.hexdistcmg='Dist. to Central Mangroves',standard.hexsgPCA1='Seagrass PCA','standard.hexdist2shore*standard.hexdistcmg'='Dist. to Shore * Dist. to Central Mangroves'))%>%
		bold_p(t=0.05)	
	p1 <- tbl_regression(pm, exp=FALSE,conf_level=0.95, label=list(standard.depth='Depth',standard.dist2jetty='Dist. to Jetty',standard.dist2shore='Dist. to Shore',standard.distcmg='Dist. to Central Mangroves', 'standard.depth:standard.dist2jetty' = ' Depth * Dist. to Jetty'))%>%
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

	## AIC of three. mdoels relative to each other
	AIC(sgm, fm, pm)

	save_as_image(stacked.summary,'hypotesting_dredge_results/stackedsummarytables_hypotheses_chp3_seagrasses_fish_largesharks_glms_may2024.png' ,webshot='webshot')

#############################
## - Residuals checking: histograms of residuals, and standardsised residuals versus fitted, looking for heteroscedasticity
## 19 Feb 2024, following feedback from RBS
	
	pacman::p_load(tidyverse, sf, ggplot2, cowplot, patchwork,lme4)
	# color guide
	  navy: '#3a6c74'; 
	  #for multiple colors: 
	  c('#3a6c74','#708d8e','#3cbcfc')
	
	# R Server, read in model RDS
	sgm <- readRDS('~/resource/data_and_RDS_NOTforupload/seagrasses_glm_hypotesting_distancemetrics_sept23.RDS')
	fm <- readRDS('~/resource/data_and_RDS_NOTforupload/fishesmetric_glm_hypotesting_distancemetrics_sept23.RDS')
	pm <- readRDS('~/resource/data_and_RDS_NOTforupload/largesharks_relativerisk_glm_hypotesting_distancemetrics_may24.RDS')

	hexdata<-read.csv('~/resource/data_and_RDS_NOTforupload/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_may24.csv')%>%
	  mutate(jcode=as.numeric(jcode))%>%
	  rename(standard.hexsgPCA1 = st_PCA1)

	pointdata <- read.csv('~/resource/data_and_RDS_NOTforupload/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_dec23.csv')%>%
	  mutate(standard.press=((as.numeric(pressure)-mean(as.numeric(pressure)))/sd(as.numeric(pressure)))) %>%
	  mutate(standard.dist2jetty=((as.numeric(dist2jetty)-mean(as.numeric(dist2jetty)))/sd(as.numeric(dist2jetty))),.after=dist2jetty)
	
  # histogram of residuals
	sgm.histresiduals <- ggplot(data = hexdata, aes(x = sgm$residuals)) +
	  geom_histogram(fill = '#708d8e', color = '#3a6c74') +
	  labs(title = 'Seagrass PCA', x = 'Residuals', y = 'Frequency')+
	  theme_bw()
	fm.histresiduals <- ggplot(data = hexdata, aes(x = fm$residuals)) +
	  geom_histogram(fill = '#708d8e', color = '#3a6c74') +
	  labs(title = 'Fish Metric', x = 'Residuals', y = 'Frequency')+
	  theme_bw()
	pm.histresiduals <- ggplot(data = pointdata, aes(x = pm$residuals)) +
	  geom_histogram(fill = '#708d8e', color = '#3a6c74') +
	  labs(title = 'Pressure', x = 'Residuals', y = 'Frequency')+
	  theme_bw()
	
	stacked.bedtdredge.residualsHistograms <- sgm.histresiduals / fm.histresiduals / pm.histresiduals
	
	ggsave(stacked.bedtdredge.residualsHistograms, file = 'hypotesting_dredge_results/dredgebestmodel_HISTresiduals_fish_sg_relativerisk_may24.png', device = 'png', units = 'in',height = 6.5, width = 3, dpi=950)
	
	## residuals versus fitted, plots
    sgm.resids <- resid(sgm)
    fm.resids <- resid(fm)
    pm.resids <- resid(pm)
    
    sgm.fitted.resids <- ggplot(data = sgm, aes(x=.fitted, y=.resid))+
      geom_point(color = '#3a6c74') +
      geom_hline(yintercept = 0,)+
      labs(title='', x='Fitted', y='Residuals')+
      theme_bw()
    fm.fitted.resids <- ggplot(data = fm, aes(x=.fitted, y=.resid))+
      geom_point(color = '#3a6c74') +
      geom_hline(yintercept = 0,)+
      labs(title='', x='Fitted', y='Residuals')+
      theme_bw()
    pm.fitted.resids <- ggplot(data = pm, aes(x=.fitted, y=.resid))+
      geom_point(color = '#3a6c74') +
      geom_hline(yintercept = 0,)+
      labs(title='', x='Fitted', y='Residuals')+
      theme_bw()

    
    stacked.bestdredge.residualsVSfitted <- sgm.fitted.resids / fm.fitted.resids / pm.fitted.resids
    
    ggsave(stacked.bestdredge.residualsVSfitted, file = 'hypotesting_dredge_results/dredgebestmodel_residualsVSfitted_fish_sg_relativerisk_may24.png', device = 'png', units = 'in',height = 6.5, width = 3, dpi=950)
    
    
    
	
#############################
## - Prediction plots GLM/GLMMs of 'best' models: the effect of distance metrics on predictor variables
#############################

# read in model RDS
	sgm <- readRDS('hypotesting_dredge_results/seagrasses_glm_hypotesting_distancemetrics_sept23.RDS')
	fm <- readRDS('hypotesting_dredge_results/fishesmetric_glm_hypotesting_distancemetrics_sept23.RDS')
	pm <- readRDS('hypotesting_dredge_results/largesharks_pressure_glm_hypotesting_distancemetrics_dec23.RDS')

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
			geom_point(data=hexdata, aes(x=standard.hexdist2shore, y=standard.hexsgPCA1,xmin=min(hexdata$standard.hexdist2shore),xmax=max(standard.hexdist2shore)),pch=21,cex=0.5,col='#3a1d20')+
			ylab('Seagrass PCA (standardised)')+xlab('Distance from Shore (standardised)')+
			geom_line(data=sgm.hex.predicts.dist2shore,aes(x=x,y=predicted),col='#653515',size=.6,linetype=1)+
			geom_ribbon(data=sgm.hex.predicts.dist2shore,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#6d2e35')+
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
	## fishes: sgPCA1, dist2shore, distcmg, dst2shore*distcmg
		
		fm.hex.predicts.sgPCA<-as.data.frame(ggpredict(fm,terms='standard.hexsgPCA1[all]',type='fixed'))
		fm.plot1<- ggplot()+
			geom_point(data=hexdata, aes(x=standard.hexsgPCA1, y=standard.hexfish),pch=21,cex=0.5,col='#836304')+
			ylab('Fishes (standardised)')+xlab('Seagrasses PCA (standardised)')+
			geom_line(data=fm.hex.predicts.sgPCA,aes(x=x,y=predicted),col='#e3ab07',size=.6,linetype=1)+
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


	##################
	## large sharks: depth, dist2jetty, distcmg, depth*dist2jetty

	pm..predicts.dist2jetty<-as.data.frame(ggpredict(pm,terms='standard.dist2jetty[all]',type='fixed'))
	pm.plot1<- ggplot()+
		geom_point(data=pointdata, aes(x=standard.dist2jetty, y=standard.press),pch=21,col='#d13474')+
		ylab('Large Sharks (standardised)')+xlab('Distance from Jetty (standardised)')+
		geom_line(data=pm..predicts.dist2jetty,aes(x=x,y=predicted),col='#a22558',size=.6,linetype=1)+
		geom_ribbon(data=pm..predicts.dist2jetty,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#e280a8')+
		theme_bw()

	pm..predicts.dist2shore<-as.data.frame(ggpredict(pm,terms='standard.dist2shore[all]',type='fixed'))
	pm.plot2<- ggplot()+
		geom_point(data=pointdata, aes(x=standard.dist2shore, y=standard.press,xmin=min(pointdata$standard.dist2shore),xmax=max(standard.dist2shore,ymin=min(pointdata$standard.press),ymax=max(pointdata$standard.press))),pch=21,col='#d13474')+
		ylab('Large Sharks (standardised)')+xlab('Distance from Shore (standardised)')+
		geom_line(data=pm..predicts.dist2shore,aes(x=x,y=predicted),col='#a22558',size=.6,linetype=1)+
		geom_ribbon(data=pm..predicts.dist2shore,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#e280a8')+
		theme_bw()

	pm..predicts.depth<-as.data.frame(ggpredict(pm,terms='standard.depth[all]',type='fixed'))
	pm.plot3<- ggplot()+
		geom_point(data=pointdata, aes(x=standard.depth, y=standard.press,xmin=min(pointdata$standard.depth),xmax=max(standard.depth)),pch=21,col='#d13474')+
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
		ylab('Predation Pressure (standardised)')+xlab('Distance from Jetty (standardised)')+
		geom_line(data=pm..predicts.dist2jetty,aes(x=x,y=predicted),col='#3a6c74',size=.6,linetype=1)+
		geom_ribbon(data=pm..predicts.dist2jetty,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.3, fill='#3a6c74')+
		theme_bw()
		
	pm..predicts.distcmg<-as.data.frame(ggpredict(pm,terms='standard.dist2shore[all]',type='fixed'))
	pm.plot2<- ggplot()+
		geom_point(data=pointdata, aes(x=standard.dist2shore, y=standard.press,xmin=min(pointdata$standard.dist2shore),xmax=max(standard.dist2shore,ymin=min(pointdata$standard.press),ymax=max(pointdata$standard.press))),pch=21,col='#3a6c74')+
		ylab('Predation Pressure (standardised)')+xlab('Distance to Shore (standardised)')+
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



































