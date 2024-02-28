 ## Path analysis for Risk & Resource BSEM (chapter 3)

## Code and approach here heavily inspired by modelling code from paper Bennett, S., Harris, M. P., Wanless, S., Green, J. A., Newell, M. A., Searle, K. R., & Daunt, F. (2022). Earlier and more frequent occupation of breeding sites during the non‐breeding season increases breeding success in a colonial seabird. Ecology and Evolution, 12(9). https://doi.org/10.1002/ece3.9213.

## created 10 May 2023 :: Molly M Kressler 

## inference from paths 2 & 3, 27 February 2024 :: Molly M Kressler & Rich B Sherley 


###################################
########## RUN AT OPEN ############
###################################

## Load Workspace, local macbook

pacman::p_load(MCMCvis,tidyverse,sf,nimble,devtools,flextable,webshot2,sfdep,sp,spdep,beepr, HDInterval)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd')

## Load workspace, R Server 
pacman::p_load(MCMCvis,tidyverse,sf,nimble,devtools,flextable,webshot2,sfdep,sp,spdep,beepr, HDInterval)
setwd("~/resource/data_and_RDS_NOTforupload")

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
		  standard.hexdist2jetty = hexdata$standard.dist2jetty
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

		# import RDS, R Server
		samplesList3b <- readRDS('mcmcsamples_model3b_niter10000_burn2000_chains3_4dec2023.RDS')
		
		# import RDS, locacl macbook 
		samplesList3b <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model3b_niter10000_burn2000_chains3_4dec2023.RDS')
			summary(samplesList3b)

			head(samplesList3b$chain1)
			summary(samplesList3b)
			str(samplesList3b)

			## How to make an object with all draws of e.g. j4. 
			j4.ch<-c(samplesList3b$chain1[,22],samplesList3b$chain2[,22],samplesList3b$chain3[,22])
			head(c3.ch)
			hdi(j4.ch)[2]
			pacman::p_load(HDInterval)

			c3.ch<-c(samplesList3b$chain1[,11],samplesList3b$chain2[,11],samplesList3b$chain3[,11])

			a1.ch<-c(samplesList3b$chain1[,1],samplesList3b$chain2[,1],samplesList3b$chain3[,1])


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

			j4_draws <- gather_draws(samplesList3b,j[])%>%
				group_by(.chain)%>%
				ungroup()
			j4<- as.data.frame(j4_draws[,5])
			head(j4)
			samplesList3b$chain

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


##################################################
## Inference from Paths ##
##################################################

		## For local macbook
		samplesList3b <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model3b_niter10000_burn2000_chains3_4dec2023.RDS')
		## For R Server
		samplesList3b <- readRDS('mcmcsamples_model3b_niter10000_burn2000_chains3_4dec2023.RDS')
		hexdata <- read.csv('data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_nov23.csv')
	  head(hexdata)	
	
	  summary(samplesList3b)
	head(samplesList3b$chain1)
	str(samplesList3b)

## Make test sample data frames - based on the hexagon df

	hexsamp <- hexdata %>%
		sample_n(5)%>%
		as.data.frame() 

## From model3b samplesList, make objects with all draws of each coeefficient from path 2 and path 3, e.g. j4. 
	# path 2: j4, c3, a1
	j4.ch<-c(samplesList3b$chain1[,22],samplesList3b$chain2[,22],samplesList3b$chain3[,22])
	c3.ch<-c(samplesList3b$chain1[,11],samplesList3b$chain2[,11],samplesList3b$chain3[,11])
	a1.ch<-c(samplesList3b$chain1[,1],samplesList3b$chain2[,1],samplesList3b$chain3[,1])

	# path 3: e1, e2, e3, e4, a7
	e1.ch <-c(samplesList3b$chain1[,14], samplesList3b$chain2[,14],samplesList3b$chain3[,14] )
	e2.ch <-c(samplesList3b$chain1[,15], samplesList3b$chain2[,15],samplesList3b$chain3[,15] )
	e3.ch <-c(samplesList3b$chain1[,16], samplesList3b$chain2[,16],samplesList3b$chain3[,16] )
	e4.ch <-c(samplesList3b$chain1[,17], samplesList3b$chain2[,17],samplesList3b$chain3[,17] )
	a7.ch <-c(samplesList3b$chain1[1:3,7], samplesList3b$chain2[,7],samplesList3b$chain3[,7] )

## Run loops for each path, to calculate path estimate given a hexagon cell

	## Variable-Coefficient key
		# j4 = dist2shore * distcmg
		# c3 = seagrasses
		# a1 = fish
		# e1 = dist2shore
		# e2 = depth
		# e3 = dist2jetty
		# e4 = dist2jetty*depth
		# a7 = predator pressure

	## Define size of objects/matrices
	J = 3 * 8000 # chains x iterations
	n.hex = nrow(hexdata)	

	## set up prediction df

		preds.path2 = as.data.frame(matrix(NA,ncol = J, nrow = n.hex))
		head(preds.path2)

	## write progress bar function
		pb <- txtProgressBar(min = 1, max = n.hex, style = 3)

	## path 2 loop
		for(i in 1:n.hex){
		  for(j in 1:J){
		    preds.path2[i,j] <- (j4.ch[j]*hexdata$standard.hexdist2shore[i]*hexdata$standard.hexdistcmg[i]) + (c3.ch[j]*hexdata$st_PCA1[i]) + (a1.ch[j]*hexdata$standard.hexfish[i])
		  } 
		  	setTxtProgressBar(pb, i)
		}
	
	## path 3 loop
	for(i in 1:n.hex){
		for(j in 1:J{
			preds.path3[i,j] <-  (e1.ch[j]*hexdata$standard.pointdist2shore[i]) + (e2.ch[j]*hexdata$standard.pointdepth[i]) + (e3.ch[j]*hexdata$standard.pointdist2jetty[i]) + (e4.ch[j]*hexdata$standard.pointdist2jetty[i]*hexdata$standard.pointdepth[i]) + (a7.ch[j]*hexdata$standard.pointpress) 
		}
	}

## Calculate marginal means, and HDI (highest density intervals)

	mean(as.numeric(preds.path2[1,1:24000]))

	cols <- c('jcode','mean', 'lower', 'upper')
	p2pred <- as.data.frame(matrix(ncol=4, nrow = n.hex))
	p3pred <- as.data.frame(matrix(ncol=4, nrow = n.hex))
	colnames(p2pred) = cols
	colnames(p3pred) = cols
	p2pred$jcode <- as.character(hexdata$jcode)
	p3pred$jcode <- as.character(hexdata$jcode)
	head(p2pred)

	pb <- txtProgressBar(min = 1, max = n.hex, style = 3)
	
	for(i in 1:n.hex){
		p2pred[i,2] <- mean(as.numeric(preds.path2[i,]))
		p2pred[i,3] <- hdi(preds.path2[i,])[2]
		p2pred[i,4] <- hdi(preds.path2[i,])[1]
		setTxtProgressBar(pb, i)
		#p3pred[i,2] <- mean(as.numeric(preds.path3[i,]))
		#p3pred[i,3] <- hdi(preds.path3[i,])[2]
		#p3pred[i,4] <- hdi(preds.path3[i,])[1]
	};beep(3)

	head(p2pred)
	head(p3pred)


## Save path estimates and path means + HDI dfs

	write.csv(preds.path2, 'path2_estimates_at_hexagons_model3bdec2023_calcFeb2024.csv')
	write.csv(preds.path3, 'path3_estimates_at_hexagons_model3bdec2023_calcFeb2024.csv')
	write.csv(p2preds, 'path2_means_andHDI_at_heaxgons_model3bdec2023_calcFeb2024.csv')
	write.csv(p3preds, 'path3_means_andHDI_at_heaxgons_model3bdec2023_calcFeb2024.csv')

##################################################
## Path predictions, spatial ##
##################################################

	hexsf <- st_as_sf(st_read('resource_chp3/standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_dec23.shp'), crs = 'WGS84')
	land <- st_as_sf(st_read('bim_onlyland_noDots.kml'), crs = 'WGS84')

	## join spatial geomtry by jcode to preds
		p2.sf <- st_as_sf(left_join(p2pred, hexsf, by='jcode'))%>%
		  dplyr::select(jcode, mean, lower, upper, geometry)
		p3.sf <- st_as_sf(left_join(p3pred, hexsf, by='jcode'))%>%
		  dplyr::select(jcode, mean, lower, upper, geometry)

	## path 2 plots - mean, upper, lower 

	p2.mean <- ggplot()+
		geom_sf(data = p2.sf, aes(fill = mean, col = mean), lwd=0)+
	  scale_fill_gradient(low='#b6e5fc',high='#3a6c74', space='Lab', aesthetics='fill',limits = c(-1,1),guide=guide_colourbar(title='Mean'))+
	  scale_colour_gradient(low='#b6e5fc',high='#3a6c74', space='Lab', aesthetics='col',limits = c(-1,1),guide=guide_colourbar(title='Mean'))+
	  theme_bw()+
		geom_sf(data = land, aes(alpha = 0, col = 'grey30'))

	p2.lower <- ggplot()+
		geom_sf(data = p2.sf, aes(fill = lower, col = NA))+
		scale_fill_gradient(low='#b6e5fc',high='#3a6c74', space='Lab', aesthetics='fill',limits = c(,),guide=guide_colourbar(title='HDI Lower'))+
		theme_bw()+
		geom_sf(data = land, aes(alpha = 0, col = 'grey30'))
aq
	p2.upper <- ggplot()+
		geom_sf(data = p2.sf, aes(fill = upper, col = NA))+
		scale_fill_gradient(low='#b6e5fc',high='#3a6c74', space='Lab', aesthetics='fill',limits = c(,),guide=guide_colourbar(title='HDI Upper'))+
		theme_bw()+
		geom_sf(data = land, aes(alpha = 0, col = 'grey30'))


	## path 3 plots - mean, upper, lower 

	p3.mean <- ggplot()+
		geom_sf(data = p3.sf, aes(fill = mean, col = NA))+
		scale_fill_gradient(low='#b6e5fc',high='#3a6c74', space='Lab', aesthetics='fill',limits = c(,),guide=guide_colourbar(title='Mean'))+
		theme_bw()+
		geom_sf(data = land, aes(alpha = 0, col = 'grey30'))

	p3.lower <- ggplot()+
		geom_sf(data = p3.sf, aes(fill = lower, col = NA))+
		scale_fill_gradient(low='#b6e5fc',high='#3a6c74', space='Lab', aesthetics='fill',limits = c(,),guide=guide_colourbar(title='HDI Lower'))+
		theme_bw()+
		geom_sf(data = land, aes(alpha = 0, col = 'grey30'))

	p3.upper <- ggplot()+
		geom_sf(data = p3.sf, aes(fill = upper, col = NA))+
		scale_fill_gradient(low='#b6e5fc',high='#3a6c74', space='Lab', aesthetics='fill',limits = c(,),guide=guide_colourbar(title='HDI Upper'))+
		theme_bw()+
		geom_sf(data = land, aes(alpha = 0, col = 'grey30'))


	## various arrangements for pub
		## just the means, side by side
		means.. <- p2.mean + p3.mean + plot.layout(guides = 'collect')

		## means with their HDIs, horizontal
		p2.mean.hdis.wide <- p2.mean + p2.lower + p2.upper + plot.layout(guides = 'collect')
		p3.mean.hdis.wide <- p3.mean + p3.lower + p3.upper + plot.layout(guides = 'collect')
		
		## means with their HDIs, stacked
		p2.mean.hdis.long <- p2.mean / p2.lower / p2.upper + plot.layout(guides = 'collect')
		p3.mean.hdis.long <- p3.mean / p3.lower / p3.upper + plot.layout(guides = 'collect')

		## means with their HDIs, means = 2 cols and big, hdis = 1 col small under means
		p2.mean.hdis.square <- p2.mean / (p2.lower + p2.upper) + plot.layout(ncol =2, nrow =2, guides = 'collect')
		p3.mean.hdis.square <- p3.mean / (p3.lower + p3.upper) + plot.layout(ncol =2, nrow =2, guides = 'collect')




