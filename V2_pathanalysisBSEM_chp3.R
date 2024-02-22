## Path analysis for Risk & Resource BSEM (chapter 3)

## Code and approach here heavily inspired by modelling code from paper Bennett, S., Harris, M. P., Wanless, S., Green, J. A., Newell, M. A., Searle, K. R., & Daunt, F. (2022). Earlier and more frequent occupation of breeding sites during the non‐breeding season increases breeding success in a colonial seabird. Ecology and Evolution, 12(9). https://doi.org/10.1002/ece3.9213.

## created 10 May 2023 by Molly M Kressler 


###################################
########## RUN AT OPEN ############
###################################

## Load Workspace 

pacman::p_load(MCMCvis,tidyverse,sf,nimble,devtools,flextable,webshot2,sfdep,sp,spdep,beepr)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd')

#####################################################
###### DF 1, for process model for sharkiness  ######
### Define Sharkiness
## counts of detections per individual per Site (buffer/receiver)

pointdata<-st_as_sf(st_read('resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_dec23.shp'),crs='WGS84')%>%
		rename(buffIDnum=bffIDnm,
			standard.shark=stndrd_s,
			standard.fish=stndrd_f,
			standard.dist2shore=stndr_2,
			standard.lds=stndrd_l,
			standard.mds=stndrd_m,
			standard.hds=stndrd_h, 
			dist2jetty=dst2jtt,
			standard.depth=stndrd_dp,
			pressure=pressur,
			standard.distcmg=stndrd_ds)%>%
			mutate(standard.press=((as.numeric(pressure)-mean(as.numeric(pressure)))/sd(as.numeric(pressure))),.before=geometry) %>%
			mutate(standard.dist2jetty=((as.numeric(dist2jetty)-mean(as.numeric(dist2jetty)))/sd(as.numeric(dist2jetty))),.after=dist2jetty)
stopifnot(nrow(pointdata)==560) # check 
sapply(pointdata,class)
summary(pointdata)

	pits <- unique(pointdata$PIT)
	tags <- list(paste0('Tag',seq(1:16)))
	both <- bind_cols(pits,tags) %>%rename(PIT=...1,tagID=...2)
	pointdata <- left_join(pointdata,both) %>%relocate(tagID, .after=PIT)
	head(pointdata)
	pointdata%>%filter(PIT=='molly5') # check that PIT is gettin assigned the same Tag ID the whole way through

	# Table of random sample of data frame for graphical methods figure 
	pointflex <- pointdata%>% 
		as_tibble%>%
		dplyr::select(-z,-buffIDnum,-PIT,-geometry)%>%
		rename('Tag' = 'tagID')%>%
		slice(1:5)%>% 
		flextable()%>%
		theme_zebra()%>%
		set_header_labels(rowname = 'Coefficient',SD='Sd')%>%
		align(align = 'center', part = 'all')%>%
		font(fontname = 'Arial', part = 'all')%>%
		color(color='black',part='all')%>%
		fontsize(size = 10, part = 'all')%>%
		autofit()
		pointflex
		save_as_image(pointflex, 'resource_chp3/forgraphicalmethodsfigure_chapter3_datapreparation_POINTDATA.png',webshot='webshot')

####################################################
###### DF 2, for process model for fishiness  ######

hexdata<-read.csv('resource_chp3/data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_nov23.csv')%>%mutate(jcode=as.numeric(jcode))
summary(hexdata)
stopifnot(nrow(hexdata)==2663) # check 

	# Table of random sample of data frame for graphical methods figure 
	hexflex <- hexdata%>%
	dplyr::rename('Hex.No' = 'jcode')%>%
		slice(1:5)%>% 
		flextable()%>%
		theme_zebra()%>%
		set_header_labels(rowname = 'Coefficient',SD='Sd')%>%
		align(align = 'center', part = 'all')%>%
		font(fontname = 'Arial', part = 'all')%>%
		color(color='black',part='all')%>%
		fontsize(size = 10, part = 'all')%>%
		autofit()
		hexflex
		save_as_image(hexflex, 'resource_chp3/forgraphicalmethodsfigure_chapter3_datapreparation_HEXDATA.png',webshot='webshot')

##################################################################

#####################	
## - MODEL3B:  Paths informed by hypothesis testng and dredge models - Active model, 17/11/2023

	modelCode3b<-nimbleCode({

	  ##########################
	  ######### priors #########
	  
	  #### prior for sharkiness #### 
	  for(i in 1:7){
		   a[i] ~ dnorm(0,.01) 
			}
		  # prior for intercept - sharks 
		  b ~ dnorm(0,.001)
		  # prior for residual variance - sharks 
		  tau.shark ~ dgamma(0.01,0.01) # 
		  sigma.shark <- sqrt(1/tau.shark)
	  
	  #### prior for predators #### 
	  for(i in 1:4){
		   e[i] ~ dnorm(0,.001) 
			}
		  # prior for intercept - sharks 
		  f ~ dnorm(0,.001)
		  # prior for residual variance - sharks 
		  tau.pred ~ dgamma(0.01,0.01) # 
		  sigma.pred <- sqrt(1/tau.pred)
	  
	  # group-level effect of buffID on sharkiness and predators.
	  for(i in 1:B){
	    epsi_shark[i]~ dnorm(0,tau.epsi_shark)
		}
		  tau.epsi_shark ~ dgamma(0.01,0.01)
	  	 sigma.epsi_shark <- sqrt(1/tau.epsi_shark)
	  
	  #### prior for fishiness (cross metric, hexagons) ####
	  for(i in 1:4){
		    c[i] ~ dnorm(0,0.001) 
			}
		  # prior for intercept - fishiness
		  d ~ dnorm(0,0.001)
		  
		  # prior for residual variance (precision) - fishiness
		  tau.fish ~ dgamma(0.01,0.01) 
		  sigma.fish <- sqrt(1/tau.fish) # gives sd 
	  
	  #### prior for seagrasses PCA @ hexagon level ####
	  for(i in 1:4){
			j[i] ~ dnorm(0,0.01) 
					}
		  # prior for intercept - hex 
			k ~ dnorm(0,0.001) # low
					
		  # prior for residual variance (precision) - hex sg PCA 
			tau.hexsg ~ dgamma(0.01,0.01)
		 	sigma.hexsg <- sqrt(1/tau.hexsg)
	

	  #########################################################
	  ######### Likelihoods - data and process models #########	

		 ## informed by hypothesis exploration 

	  ### data model for seagrasses - hexagons
	    for(i in 1:hex.N){
	    standard.hexsg[i] ~ dnorm(hexsg.mu[i],tau.hexsg)

	    hexsg.mu[i] <- k + j[1]*standard.hexdist2shore[i] + j[2]*standard.hexdistcmg[i] + j[3]*standard.hexdist2jetty[i] + j[4]*standard.hexdist2shore[i]*standard.hexdistcmg[i] 
		}
	  
	  ### data model for fishiness - hexagons 
	   for(i in 1:hex.N){
	    standard.hexfish[i] ~ dnorm(hexfish.mu[i],tau.fish) 

	    hexfish.mu[i] <- d + c[1]*standard.hexdist2shore[i] + c[2]*standard.hexdistcmg[i] + c[3]*standard.hexsg[i]+ c[4]*standard.hexdist2shore[i]*standard.hexdistcmg[i]
		}
	
	  ### data model for large shark detectons - point data 
	   for(i in 1:point.N){
	   	standard.pointpress[i] ~ dnorm(pred.mu[i], tau.pred)

	   	pred.mu[i] <- f + e[1]*standard.pointdist2shore[i] + e[2]*standard.pointdepth[i]+ e[3]*standard.pointdist2jetty[i] + e[4]*standard.pointdist2jetty[i]*standard.pointdepth[i]
	   }

	  ### data model for sharkiness - pointdata
	   for(i in 1:point.N){
	    standard.shark[i] ~ dnorm(shark.mu[i],tau.shark)   

	    shark.mu[i] <- b + epsi_shark[buffID[i]] + a[1]*standard.fish.pred[i] + a[2]*standard.pointdist2shore[i] + a[3]*standard.pointdistcmg[i] + a[4]*sgPCA1[i] + a[5]*standard.pointdepth[i]+ a[6]*standard.pointdist2jetty[i]+ a[7]*standard.pointpress[i]
	  	}  

	  ######### Derived Parameters #########
	  # for estmating total pathways 
		## coefficients for distance metrics are from process models of those predictors
		## paths 1-3 assess the support of the effect of the primary parameter, e.g. pressure, on juvenile spatial baheviour as that primary parameter is determined/influenced by the abiotic habitat features, e.g. depth  

		path[1] <-  j[1] * j[2] * j[3] * j[4] * a[4]  # seagrass dredge informed path: sg dist2jetty dist2shore distcmg 

		path[2] <-  j[4] * c[3] * a[1] # fish dredge informed path: fish sg distcmg dist2shore

		path[3] <-  e[1] * e[2] * e[3] * e[4] * a[7] # predator pressure informed path: predators depth dist2shore dist2jetty

		path[4] <- a[3] * a[6] # refugia and anthropocene - could directly impact movement, so ecologically, a path (disturbance)

		## need to calculate the values from the abiotics along the paths to the initial parameter, e.g. abtiocs to fishes in path 1. 

		value[1] <-  j[1] * j[2] * j[3] * j[4] # path 2, shore * refuge * jetty
		value[2] <- e[1] * e[2] * e[3] * e[4] # path 3, shore * jetty * depth 

		### predictions based on paths 2 & 3

		for(v in 1:30){
		p2.mu[v] <- j[4]*shore.p2.dum[v]*cmg.dum[v] * c[3]*sg.dum[v] * a[1]*fish.dum[v]
			} # problem(?) here is that the values of 'v' for different things dont relate. so is this the nesting thing? where I need to have it run through all 30 of shore.p2.dum, holding c3 and a1 to means; then repeat for c3 and a1? Do do I actually need three dataframes for path 2, one for each predictor/coeff where it runs from min to max and the others hold constant? but then how to bring them together. 22/2/2024

		for(z in 1:30){
		p3.mu[z] <- e[1]*shore.p3.dum[z] * e[2]*standard.pointdepth[z]* e[3]*standard.pointdist2jetty[z] * e[4]*standard.pointdist2jetty[z]*standard.pointdepth[z]* a[7]*standard.pointpress[z]
			}

	}) # end of model code 


	##########################
	## Compile the model code
	##########################

		myConstants3b<-list(point.N=560,hex.N=2663,B=max(pointdata$buffIDnum),buffID=pointdata$buffID, v = nrow(path2.dum), z = nrow(path3.dum))

		myData3b<-list(
		  # tell nimble the covariates 
		  	# point data
		  standard.shark = pointdata$standard.shark,
		  standard.pointdist2shore = pointdata$standard.dist2shore,
		  standard.pointdistcmg = pointdata$standard.distcmg,
		  standard.fish.pred = pointdata$standard.fish,
		  sgPCA1 = pointdata$sgPCA1,
		  standard.pointdepth= pointdata$standard.depth,
		  standard.pointdist2jetty= pointdata$standard.dist2jetty,
		  standard.pointpress = pointdata$standard.press,
		  	# hex data 
		  standard.hexfish = hexdata$standard.hexfish,
		  standard.hexdist2shore = hexdata$standard.hexdist2shore,
		  standard.hexdistcmg = hexdata$standard.hexdistcmg,
		  standard.hexdist2jetty = hexdata$standard.dist2jetty,
		  standard.hexsg = hexdata$st_PCA1,
		  standard.hexdist2jetty = hexdata$standard.dist2jetty,
		  sg.dum = path2.dum$sg.dum,
		  cmg.dum = path2.dum$cmg.dum,
		  shore.p2.dum = path2.dum$shore.dum,
		  shore.p3.dum = path3.dum$shore.dum,
		  depth.dum = path3.dum$depth.dum,
		  jetty.dum = path3.dum$jetty.dum,
		  fish.dum = path2.dum$fish.dum,
		  press.dum = path3.dum$press.dum
			)

		init.values3b<-list(a=rnorm(7,0,1),
			b=rnorm(1),
			c=rnorm(4,0,1),
			d=rnorm(1),
			e=rnorm(4,0,1),
			f=rnorm(1),
			j=rnorm(4,0,1),
			k=rnorm(1),
			path=rnorm(4,0,.05),
			value=rnorm(2,0,.05),
			epsi_shark=rgamma(35,0.01,0.01),
			tau.shark=rgamma(1,0.01,0.01),
			tau.pred=rgamma(1,0.01,0.01),
			tau.fish=rgamma(1,0.01,0.01),
			tau.epsi_shark=rgamma(1,0.01,0.01),
			tau.hexsg=rgamma(1,0.01,0.01)
				)
			
		model3b<-nimbleModel(code=modelCode3b, name="model3b",data=myData3b,constants = myConstants3b,inits=init.values3b) #define the model

			model3b$calculate() # if = NA, indicates missing or invalid initial values, and you have to fix the model until it is numeric.
			model3b$initializeInfo()

		Cm3b<-compileNimble(model3b);beep(2) # compile the model


	##########################
	## Compile & Run the MCMC
	##########################

			conf3b <- configureMCMC(model3b,monitors=c('a','b','c','d','j','k','e','f','tau.epsi_shark','tau.fish','tau.shark','tau.pred','tau.hexsg','path','value'),onlySlice=FALSE) 
			MCMC_model3b <- buildMCMC(conf3b,na.rm=TRUE)
			ccMCMC3b <-compileNimble(MCMC_model3b, project = model3b)
			samples3b <- runMCMC(ccMCMC3b,niter=10000, nburnin=2000, nchains=3,samplesAsCodaMCMC = TRUE);beep(2)

		saveRDS(samples3b,'resource_chp3/nimblemodel_outputs/mcmcsamples_model3b_niter10000_burn2000_chains3_4dec2023.RDS')
		
		saveRDS(samples3b,'resource_chp3/nimblemodel_outputs/mcmcsamples_model3b_niter10000_burn2000_chains3_22feb2024.RDS')
		
		mcmc_summary_Cmodel3b<-MCMCsummary(samples3b,round=4,pg0=TRUE,prob=c(0.05,0.95))%>%
				tibble::rownames_to_column()%>%
				rename_with(str_to_title)%>%
				rename(Parameter = Rowname)%>%
				rename('% of posterior with \n\ same sign as estimate' = 'P>0', Estimate = 'Mean','lower'='5%',upper='95%')%>%
				mutate(CI = paste0('[',lower,',',upper,']'),.after='Estimate')%>%
				dplyr::select(-lower,-upper,-Sd)%>%	
				flextable()%>%
				theme_zebra()%>%
				set_header_labels(rowname = 'Coefficient',SD='Sd')%>%
				align(align = 'center', part = 'all')%>%
				font(fontname = 'Arial', part = 'all')%>%
				color(color='black',part='all')%>%
				fontsize(size = 10, part = 'all')%>%
				autofit()
		mcmc_summary_Cmodel3b
		
		# caterpillar plots 

		MCMCplot(samples3b,ci=c(50,95),params=c('path')) # point = median, thick line = 50% CI, thin line = 95% CI 

		# trace and density plots

		MCMCtrace(samplesList3b,pdf=TRUE,ind=TRUE, Rhat=TRUE, n.eff=TRUE) # ind = TRUE, separate density lines per chain. # pdf = FALSE, don't export to a pdf automatically. 


	###########################################################
	## Summary Table & Caterpillar plots with MCMCvis & tidybayes to show small values ##
	###########################################################

		## Feb 2024 - need to amend the pg0

		pacman::p_load(tidybayes,bayesplot,ggdist,nlist,forcats,patchwork)

		# import RDS

		samplesList3b <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model3b_niter10000_burn2000_chains3_4dec2023.RDS')
			summary(samplesList3b)

		mcmc_summary_Cmodel3b_samplesListfromRDS<-MCMCsummary(samplesList3b,round=3,probs=c(0.05,0.95),pg0=TRUE)%>%
				tibble::rownames_to_column()%>%
				rename_with(str_to_title)%>%
			rename('pg0'='P>0')%>%
		  mutate(pg00 = case_when(Mean > 0 ~ as.numeric(pg0), Mean < 0 ~ 1-as.numeric(pg0), .default = as.numeric(pg0)))%>%
				rename(Parameter = Rowname, 'Prop. of posterior with \n\ same sign as estimate' = 'pg00', Estimate = 'Mean','lower'='5%',upper='95%')%>%
				mutate(CI = paste0('[',lower,',',upper,']'),.after='Estimate')%>%
				dplyr::select(-lower,-upper,-Sd, -pg0)%>%	
				filter(Parameter!=c('value[1]','value[2]'))%>%
				mutate_if(is.numeric, round, digits = 3)%>%
				flextable()%>%
				theme_zebra()%>%
				set_header_labels(rowname = 'Coefficient',SD='Sd')%>%
				align(align = 'center', part = 'all')%>%
				font(fontname = 'Arial', part = 'all')%>%
				color(color='black',part='all')%>%
				fontsize(size = 10, part = 'all')%>%
				autofit()
		mcmc_summary_Cmodel3b_samplesListfromRDS

		save_as_image(mcmc_summary_Cmodel3b_samplesListfromRDS,path='resource_chp3/nimblemodel_outputs/mcmcsamples_model3b_niter20000_burn2000_chains3_4dec2023.png',res=850)	

		# grab draws with gather_draws and create label for paths based on iterations and sequence of paths minN to maxN. 

		d3b <- gather_draws(samplesList3b,path[])%>%
				group_by(.chain)%>%
				mutate(pathID = paste0('path',rep(1:4, each=8000)))%>% 
				mutate(pathIDnum = rep(1:4, each=8000))%>%
				ungroup() # each path estimate for each chain (of 3) in order, starting with path[1] first estimate in chain 1 
			head(d3b)
			summary(d3b)
			unique(d3b$pathID) # check

			write.csv(d3b,'resource_chp3/nimblemodel_outputs/samples3b_spreadlong_model3b_niter10000_burn2000_chains3_4dec2023.csv')

			draws_3b_a <- gather_draws(samplesList3b,a[])%>%
				group_by(.chain)%>%
				mutate(a.index = rep(1:7, each=8000))%>%
				mutate(aa = paste0(.variable,a.index))%>%
				ungroup()
			draws_3b_a




		## pivot_wider to spread the pathID column into multiple columns, and fill with .value. 

			w3b <- pivot_wider(d3b%>%dplyr::select(-pathIDnum), names_from=pathID, values_from=.value)
			nrow(w3b)
			head(w3b)
			names(w3b)
			summary(w3b)

			write.csv(w3b,'resource_chp3/nimblemodel_outputs/samples3a_pivotwider_model3a_niter20000_burn2000_chains3_4dec2023.csv')

		# use tidybayes to plot 

		caterpillars <- ggplot(d3b, aes(y=reorder(pathID,pathIDnum,decreasing=T),x=.value))+
			stat_pointinterval(.width=c(.50,.95),point_size=2)+
			ylab('Path ID')+
			xlab('Estimate (mean) & CI (.5,.95)')+
			geom_vline(xintercept=0,linetype=3)+
			theme_bw()
		caterpillars	
		ggsave(caterpillars,file='resource_chp3/nimblemodel_outputs/caterpillarsPlot_mcmcsamples_model3b_niter10000_burn2000_chains3_4dec2023.png',device='png',dpi=400,width=5,height=5,units='in')


	###########################################################
	## Table: pathway description, estimates, CI and Support ##
	###########################################################
		## For local macbook
		samplesList3b <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model3b_niter10000_burn2000_chains3_4dec2023.RDS')

		## For R Server
		pacman::p_load(tidybayes,bayesplot,MCMCvis,ggdist,nlist,forcats,patchwork)
		pacman::p_load(MCMCvis)
		
		samplesList3b <- readRDS('~/resource/data_and_RDS_NOTforupload/mcmcsamples_model3b_niter10000_burn2000_chains3_4dec2023.RDS')
		##
		samplesList3b %>% filter(Rowname=='path')
		
		pathwayresults_table

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
				compose(i=2,j=2, as_paragraph('Juvenile sharks ~ Dist. to Jetty + Dist. to Shore + \n\ Dist. to Refuge + Seagrasses'))%>%
				compose(i=3,j=2, as_paragraph('Juvenile sharks ~ Depth + Dist. to Shore +  \n\ Dist. to Jetty + Predator Pressure'))%>%
				compose(i=4,j=2, as_paragraph('Juvenile sharks ~ Dist. to Refuge + Dist. to Jetty'))%>%
				theme_zebra()%>%
				align(j=3:4, align = 'center', part = 'all')%>%
				font(fontname = 'Arial', part = 'all')%>%
				color(color='black',part='all')%>%
				fontsize(size = 10, part = 'all')%>%
				autofit()
		pathwayresults_table

		save_as_image(pathwayresults_table,path='resource_chp3/nimblemodel_outputs/pathwayresultssummary_model3b_niter20000_burn2000_chains3_4dec2023.png')	


	###########################################################
	## effects plots attempts ##
	###########################################################
	## For local macbook
		samplesList3b <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model3b_niter10000_burn2000_chains3_4dec2023.RDS')
		pacman::p_load(tidybayes)

		## make dummy data frames 
			sg.dum <- seq(min(pointdata$sgPCA1), max(pointdata$sgPCA1), (max(pointdata$sgPCA1)-min(pointdata$sgPCA1))/30)
			cmg.dum <- seq(min(pointdata$standard.distcmg), max(pointdata$standard.distcmg), (max(pointdata$standard.distcmg)-min(pointdata$standard.distcmg))/30)
			shore.dum <- seq(min(pointdata$standard.dist2shore), max(pointdata$standard.dist2shore), (max(pointdata$standard.dist2shore)-min(pointdata$standard.dist2shore))/30)
			depth.dum <- seq(min(pointdata$standard.depth), max(pointdata$standard.depth), (max(pointdata$standard.depth)-min(pointdata$standard.depth))/30)
			jetty.dum <- seq(min(pointdata$standard.dist2jetty), max(pointdata$standard.dist2jetty), (max(pointdata$standard.dist2jetty)-min(pointdata$standard.dist2jetty))/30)
			fish.dum <- seq(min(pointdata$standard.fish), max(pointdata$standard.fish), (max(pointdata$standard.fish)-min(pointdata$standard.fish))/30)
			press.dum <- seq(min(pointdata$standard.pres), max(pointdata$standard.pres), (max(pointdata$standard.pres)-min(pointdata$standard.pres))/30)

			path2.dum <- as.data.frame(cbind(sg.dum,cmg.dum,shore.dum,fish.dum))
			path3.dum <- as.data.frame(cbind(shore.dum,depth.dum,jetty.dum,press.dum))


		## path[2] <-  j[4] * c[3] * a[1] # fish dredge informed path: fish sg distcmg dist2shore
		## path[3] <-  e[1] * e[2] * e[3] * e[4] * a[7] # predator pressure informed path: predators depth dist2shore dist2jetty
			model3b$getVarNames() # prints variable names
			model3b$getNodeNames() # prints nodes 
			model3b$expandNodeNames('a') # prints all nodes of the variable 'a'
			# you can swap model3b and Cm3b in these node/var calls

	



	###########################################################
	## Residuals of paths using tidybayes ##
	###########################################################
	## For local macbook
		samplesList3b <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model3b_niter10000_burn2000_chains3_4dec2023.RDS')

		## For R Server
		pacman::p_load(tidybayes,bayesplot,MCMCvis,ggdist,nlist,forcats,patchwork)
		pacman::p_load(MCMCvis)
		
		samplesList3b <- readRDS('~/resource/data_and_RDS_NOTforupload/mcmcsamples_model3b_niter10000_burn2000_chains3_4dec2023.RDS')
		##


	## add_residual_draws
		pointdata %>%
			add_residual_draws(samplesList3b)%>%
			ggplot(aes(x=.row,y = .residual))+
			stat_pointinterval()


