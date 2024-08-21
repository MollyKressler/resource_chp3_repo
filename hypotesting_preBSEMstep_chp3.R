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

setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd')

pointdata<-read.csv('standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_july24.csv') 

hexdata<-read.csv('data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_andBRTs_july24.csv')%>%
		mutate(jcode=as.numeric(jcode))

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

	sg <- glm(standard.sgPCA1 ~ standard.depth * st.dist2shore * standard.distcmg * st.dist2jetty, data=hexdata) #define the global model

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
## - Hypothesis 2 : Fish and abiotics 
#############################

	 ## brief foray into juvenile sharks 
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

	######
  ## fish, hex

	#####
	fish.glm <- glm(standard.hexfish ~ standard.depth * st.dist2shore * st.dist2jetty * standard.distcmg + standard.sgPCA1, data=hexdata) #define the global model

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

	  save_as_image(fish.modelsranked.tabled,'resource_chp3/hypotesting_dredge_results/dredged_results_fishHEXresponse_July2024.png',webshot='webshot')
	get.models(fish.dredge,1)[[1]]


	## deprecated - gerries only, SD only 
		gerr.glm <- glm(standard.hexgerr ~ standard.depth * standard.hexdist2shore * standard.hexdistcmg * standard.hexdist2jetty + standard.sgPCA, data=hexdata) #define the global model

		gerr.dredge <- dredge(gerr.glm, beta='sd', evaluate=TRUE, trace=FALSE, extra='R^2',m.lim=c(0,4)) # dredge from the global model

		gerr.modelsranked.tabled <- gerr.dredge%>%
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

		  save_as_image(gerr.modelsranked.tabled,'dredged_results_Gerreidae_SDonly_HEXresponse_June2024.png',webshot='webshot')
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

	sgm <- get.models(sg.dredge,1)[[1]]

# H2: fish and distance metrics
		fm <- get.models(fish.dredge,1)[[1]]
		fm
	# deprecated
		fm 	<- glm(standard.hexfish ~ standard.hexdist2jetty + standard.hexdistcmg + standard.sgPCA + standard.hexdist2jetty*standard.hexdistcmg, data=hexdata, family=gaussian)

		gm <- glm(standard.hexgerr ~ standard.hexdist2jetty + standard.hexdistcmg + standard.sgPCA + standard.hexdist2jetty*standard.hexdistcmg, data=hexdata, family=gaussian) # SD gerries only 

# H3: large sharks and distance metrics

	pm <- 	get.models(pred.dredged,1)[[1]]
	pmm<-glm(standard.press ~ standard.depth + standard.dist2jetty + standard.distcmg + standard.depth*standard.dist2jetty, data=pointdata)

# save the models as RDS 
	saveRDS(sgm,'resource_chp3/hypotesting_dredge_results/seagrasses_glm_hypotesting_distancemetrics_july24.RDS')
	saveRDS(fm,'resource_chp3/hypotesting_dredge_results/fishesmetric_glm_hypotesting_distancemetrics_MaxN_July24.RDS')
	#saveRDS(gm,'resource_chp3/hypotesting_dredge_results/gerries_sdonly_glm_hypotesting_distancemetrics_june24.RDS')
	saveRDS(pm,'resource_chp3/hypotesting_dredge_results/largesharks_relativerisk_glm_hypotesting_distancemetrics_may24.RDS')

#############################
## - Summary table of GLM/GLMMs of 'best' models: the effect of distance metrics on predictor variables
#############################

# read in model RDS
	sgm <- readRDS('resource_chp3/hypotesting_dredge_results/seagrasses_glm_hypotesting_distancemetrics_july24.RDS')
	fm <- readRDS('resource_chp3/hypotesting_dredge_results/fishesmetric_glm_hypotesting_distancemetrics_MaxN_July24.RDS')
	pm <- readRDS('resource_chp3/hypotesting_dredge_results/largesharks_relativerisk_glm_hypotesting_distancemetrics_may24.RDS')

# make table
	models <- c(sgm, fm, pm)
	
	sgf <- flextable::as_flextable(sgm)
	pf <- flextable::as_flextable(pm)
	ff <- flextable::as_flextable(fm)
	gf <- flextable::as_flextable(gm)

	s1 <- tbl_regression(sgm,exp=FALSE,conf_level=0.95,label=list(standard.depth='Depth',standard.distcmg='Dist. to Central Mangroves',st.dist2jetty='Dist. to Jetty','st.dist2jetty:standard.distcmg'='Dist. to Shore * Dist. to Central Mangroves'))%>%
		bold_p(t=0.05)
	f1 <- tbl_regression(fm, exp=FALSE,conf_level=0.95,label=list(st.dist2shore='Dist. to Shore',st.dist2jetty='Dist. to Nearest Jetty',standard.sgPCA1='Seagrass PCA','st.dist2jetty:st.dist2shore '='Dist. to Nearest Jetty * Dist. to Shore'))%>%
		bold_p(t=0.05)	
	p1 <- tbl_regression(pm, exp=FALSE,conf_level=0.95, label=list(standard.depth='Depth',standard.dist2jetty='Dist. to Jetty',standard.dist2shore='Dist. to Shore',standard.distcmg='Dist. to Central Mangroves', 'standard.depth:standard.dist2jetty' = ' Depth * Dist. to Jetty'))%>%
		bold_p(t=0.05)

	# side by side
	tbl_merge(tbls = list(s1,f1,p1), tab_spanner = c('Hypothesis 1', 'Hypothesis 2', 'Hypothesis 3'))
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

	save_as_image(stacked.summary,'resource_chp3/hypotesting_dredge_results/stackedsummarytables_hypotheses_chp3_seagrasses_fish_largesharks_glms_July2024.png' ,webshot='webshot')
	save_as_docx(stacked.summary,path = 'resource_chp3/hypotesting_dredge_results/stackedsummarytables_hypotheses_chp3_seagrasses_fish_largesharks_glms_July2024.docx')

#############################
## - Residuals checking: histograms of residuals, and standardsised residuals versus fitted, looking for heteroscedasticity
#############################
    
    # calculate fitted and pull residuals 
		R1_sgm <- as_tibble(resid(sgm))%>%rename(resids = value)
		F1_sgm <- as_tibble(predict(sgm))%>%rename(fitted = value)
		diag_sgm <- bind_cols(F1_sgm, R1_sgm)

		R1_fm <- as_tibble(resid(fm))%>%rename(resids = value)
		F1_fm <- as_tibble(predict(fm))%>%rename(fitted = value)
		diag_fm <- bind_cols(F1_fm, R1_fm)

		R1_pm <- as_tibble(resid(pm))%>%rename(resids = value)
		F1_pm <- as_tibble(predict(pm))%>%rename(fitted = value)
		diag_pm <- bind_cols(F1_pm, R1_pm)

		# plots 
		resVfit_sgm <- ggplot(data = diag_sgm,aes(x= fitted, y = resids))+
			geom_point(col = '#3A6C74', fill = '#3A6C74')+ 
			labs(subtitle = 'Residuals vs. Fitted', y = 'Residuals', x = 'Fitted')+
			theme_bw()+
			ggtitle('Seagrasses')
		res_sgm_hist <- ggplot(data = R1_sgm, aes(x = resids))+ 
			geom_histogram(binwidth = nrow(R1_sgm)/10000,col = '#3A6C74', fill = '#3A6C74')+ 
			theme_bw()+
			labs(subtitle = 'Residuals', y = 'Frequency', x = 'Residuals')
		res_sgm_qq <- ggplot(data=F1_sgm, aes(sample = fitted))+
			stat_qq(size=2,pch=21,col = '#3A6C74')+
			labs(subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
			stat_qq_line(linetype=2, col='red')+
			theme_bw()

		resVfit_fm <- ggplot(data = diag_fm,aes(x= fitted, y = resids))+
			geom_point(col = '#3A6C74', fill = '#3A6C74')+ 
			labs(subtitle = 'Residuals vs. Fitted', y = 'Residuals', x = 'Fitted')+
			theme_bw()+
			ggtitle('Prey MaxN')
		res_fm_hist <- ggplot(data = R1_fm, aes(x = resids))+ 
			geom_histogram(binwidth = nrow(R1_fm)/10000,col = '#3A6C74', fill = '#3A6C74')+ 
			theme_bw()+
			labs(subtitle = 'Residuals', y = 'Frequency', x = 'Residuals')
		res_fm_qq <- ggplot(data=F1_fm, aes(sample = fitted))+
			stat_qq(size=2,pch=21,col = '#3A6C74')+
			labs(subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
			stat_qq_line(linetype=2, col='red')+
			theme_bw()

		resVfit_pm <- ggplot(data = diag_pm,aes(x= fitted, y = resids))+
			geom_point(col = '#3A6C74', fill = '#3A6C74')+ 
			labs(subtitle = 'Residuals vs. Fitted', y = 'Residuals', x = 'Fitted')+
			theme_bw()+
			ggtitle('Predation Risk')
		res_pm_hist <- ggplot(data = R1_pm, aes(x = resids))+ 
			geom_histogram(binwidth = nrow(R1_pm)/1000,col = '#3A6C74', fill = '#3A6C74')+ 
			theme_bw()+
			labs(subtitle = 'Residuals', y = 'Frequency', x = 'Residuals')
		res_pm_qq <- ggplot(data=F1_pm, aes(sample = fitted))+
			stat_qq(size=2,pch=21,col = '#3A6C74')+
			labs(subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
			stat_qq_line(linetype=2, col='red')+
			theme_bw()

		diagnostics_sgm <- resVfit_sgm+res_sgm_hist+res_sgm_qq
		diagnostics_fm <- resVfit_fm+res_fm_hist+res_fm_qq
		diagnostics_pm <- resVfit_pm+res_pm_hist+res_pm_qq

		ggsave(diagnostics_sgm, file = 'resource_chp3/hypotesting_dredge_results/diagnostic_plots_sagrasses_hypotestingGLM_july2024.png', device = 'png', unit = 'in', width = 8, height = 4, dpi = 850)
		ggsave(diagnostics_fm, file = 'resource_chp3/hypotesting_dredge_results/diagnostic_plots_fishMaxN_hypotestingGLM_july2024.png', device = 'png', unit = 'in', width = 8, height = 4, dpi = 850)
		ggsave(diagnostics_pm, file = 'resource_chp3/hypotesting_dredge_results/diagnostic_plots_relrisk_hypotestingGLM_july2024.png', device = 'png', unit = 'in', width = 8, height = 4, dpi = 850)

    
	
#############################
## - (as of July'24, not in the manuscript) Prediction plots GLM/GLMMs of 'best' models: the effect of distance metrics on predictor variables
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



































