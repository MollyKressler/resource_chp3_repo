## Path analysis for fish abundance, shark presence, and habitat composition, on each other.

## Modelling code from paper Bennett, S., Harris, M. P., Wanless, S., Green, J. A., Newell, M. A., Searle, K. R., & Daunt, F. (2022). Earlier and more frequent occupation of breeding sites during the non‐breeding season increases breeding success in a colonial seabird. Ecology and Evolution, 12(9). https://doi.org/10.1002/ece3.9213.

## created 10 May 2023 by Molly M Kressler 


#########################################################################
######################### RUN AT OPEN ###################################

## Load Workspace 

pacman::p_load(tidyverse,sf,nimble,MCMCvis,flextable,webshot2)

setwd('/Users/mollykressler/Documents/data_phd/resource_chp3')

#####################################################
###### DF 1, for process model for sharkiness  ######
### Define Sharkiness
## counts of detections per individual per Site (buffer/receiver)

pointdata<-read.csv('standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_may23.csv')%>%dplyr::select(-X)#%>%rename('prop_brs'='prp_','prop_ldsg'='prp_lds','prop_medsg'='prp_mds','prop_hdsg'='prp_hds','prop_sarg'='prp_srg','prop_marsh'='prp_mrs','prop_urb_r'='prp_rb_','prop_deep'='prop_dp','prop_spong'='prp_spn','prop_unkwn'='prp_nkw')
pointdata$PIT<-as_factor(pointdata$PIT)
pointdata$buffID<-as.numeric(pointdata$buffID)
pointdata$buffIDnum<-as.numeric(pointdata$buffIDnum)
	
stopifnot(nrow(pointdata)==560) # check 
names(pointdata)

summary(pointdata)

####################################################
###### DF 2, for process model for fishiness  ######

hexdata<-read.csv('standardisedmeancentred_data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_fromPRuse_andBRTs_may23.csv')%>%dplyr::select(-X)%>%
	mutate(jcode=as.numeric(jcode))
	summary(hexdata)
## nimble needs only numeric data
sapply(hexdata,class)


######################### END ###################################
#################################################################



########################################################
################# BUGS models with NIMBLE ##############

#### overview of models 
	## model1a - uses only DF1, and doesn't work as ntended becuase site explains all the variance in fish (bc of how the data is structured).
	## model2a - uses DF1 for modelling sharks and coefficients for shark predictors; and DF2 for fishes and coefficients for fish predictors. 

	## for habitat variables, we used the types identified in previous modelling as influencing habitat selection. For sharks, this is distance to the central mangrove point and medium density seagrass; for Gerreidae, this is distance to shore and low density seagrass; and for the diversity index for species, we included only the habitats with relative influence in the BRTs > 15%, which are distance to the shore and low density seagrass.

	## in all models, variables are mean centred and standardised 

##########################################
########### MODEL 1A, RECEIVER LEVEL MODEL - shark point detectons (counts per receiver for each individual, sample size N=35), Gerreidae abundance predictions (BRT simplified model), habitat 

	## random effect for SITE but not for shark individual
	# all vars as normal, because they are + and - and centred on 0. 


	##########################

	modelCode1a<-nimbleCode({
	  
	  ##########################
	  ######### priors #########
	  
	  #### prior for sharkiness #### 
	  for(i in 1:5){
	    a[i] ~ dnorm(0,.001) 
		}
	  # prior for intercept - sharks 
	  b ~ dnorm(0,.001)
	  # prior for residual variance - sharks 
	  tau.shark ~ dgamma(0.01,0.01) # 
	  sigma.shark <- sqrt(1/tau.shark)
	  
	  # group-level effect of buffID on sharkiness. prior for variance (epsi_site) calculated from inverse of precision (epsi_site)
	  for(i in 1:B){
	    epsi_shark[i]~ dnorm(0,tau.epsi_shark)
		}
	  tau.epsi_shark ~ dgamma(0.01,0.01)
	  sigma.epsi_shark <- sqrt(1/tau.epsi_shark)
	  
	  
	  #### prior for fishiness (Gerreidae) ####
	  for(i in 1:5){
	    c[i] ~ dnorm(0,0.001) 
		}
	  # prior for intercept - fishiness
	  d ~ dnorm(0,0.001)
	  
	  # prior for residual variance (precision) - fishiness
	  tau.fish ~ dgamma(0.01,0.01) 
	  sigma.fish <- sqrt(1/tau.fish) # gives sd 
	  
	  # group-level effect for site ID on fishiness (epsi_fish)
	  # prior for precision of random effect (buffID), then calculated to variance 
	  for(i in 1:B){ # B is no. of buffers
	    epsi_fish[i] ~ dnorm(0,tau.epsi_fish)
		}
	  tau.epsi_fish ~ dgamma(0.01,0.01)
	  sigma.epsi_fish <- sqrt(1/tau.epsi_fish)
	  
	  #########################################################
	  ######### Likelihoods - data and process models #########	
	  
	  ### data model for sharkiness 
	  for(i in 1:N){
	    standard.shark[i] ~ dnorm(shark.mu[i],tau.shark)   
	  # process model for sharkiness ~ fish + dist2shore + dist_cmg + prop_ldsg + prop_medsg ? + epsi_site?
	    shark.mu[i] <- b + epsi_shark[buffID[i]] + a[1]*standard.fish.pred[i] + a[2]*standard.dist2shore[i] + a[3]*standard.distcmg[i] + a[4]*standard.lds[i] + a[5]*standard.mds[i] 
	  ### data model for fishiness 
	    standard.fish[i] ~ dnorm(fish.mu[i],tau.fish) 
	  # process model for fishiness ~ shark + dist2shore + dist_cmg + prop_ldsg + prop_medsg  + epsi_fish
	    fish.mu[i] <- d + epsi_fish[buffID[i]]+ c[1]*standard.shark.pred[i] + c[2]*standard.dist2shore[i] + c[3]*standard.distcmg[i] + c[4]*standard.lds[i] + c[5]*standard.mds[i] 
		}

	  ######### derived parameters #########
	  # for estmating total pathways 
	  # write one for each route through the pathway diagram
	  
	  path[1]<- a[2]*a[3]*a[4]*a[5] 
	  # shark ~ dist2shore + dist_cmg + prop_ldsg + prop_medsg 
	  path[2]<- c[2]*c[3]*c[4]*c[5] 
	  # fish ~ dist2shore + dist_cmg + prop_ldsg + prop_medsg 
	  path[3]<- a[1] 
	  # shark ~ fish 
	  path[4]<- c[1] 
	  # fish ~ shark 
	  path[5]<- a[1]*c[2]*c[3]*c[4]*c[5] 
	  # shark ~ fish + dist2shore + dist_cmg + prop_ldsg + prop_medsg # with coeff for habitat as they went into fish and came through
	  path[6]<- c[1]*a[2]*a[3]*a[4]*a[5]
	  # fish ~ shark + dist2shore + dist_cmg + prop_ldsg + prop_medsg 
	  
	  
	}) # end of model code 

	##########################
	## Compile the model code
	##########################

		myCovs<-c('standard.fish','standard.shark','standard.dist2shore','standard.distcmg','standard.lds','standard.mds','buffID')
		myConstants<-list(N=560,B=max(pointdata$buffIDnum),buffID=pointdata$buffIDnum) 
		myData<-list(
		  # tell nimble the covariates 
		  standard.fish = pointdata$standard.fish,
		  standard.shark = pointdata$standard.shark,
		  standard.dist2shore = pointdata$standard.dist2shore,
		  standard.distcmg = pointdata$standard.distcmg,
		  standard.lds = pointdata$standard.lds,
		  standard.mds = pointdata$standard.mds,
		  standard.shark.pred = pointdata$standard.shark,
		  standard.fish.pred = pointdata$standard.fish
		  # tell nimble the constants
		  #N=560,
		  #buffID=as.numeric(factor(pointdata$buffID))
		)

		init.values<-list(a=rnorm(5,0,1),b=rnorm(1),c=rnorm(5,0,1),d=rnorm(1),tau.shark=rgamma(1,0.01,0.01),tau.fish=rgamma(1,0.01,0.01),tau.epsi_shark=rgamma(1,0.01,0.01),tau.epsi_fish=rgamma(1,0.01,0.01),standard.fish=rnorm(myConstants$N,0,1),standard.shark=rnorm(myConstants$N,0,1))
			
		model1a<-nimbleModel(code=modelCode1a, name="model1a",data=myData,constants = myConstants,inits=init.values) #define the model
		Cm1a<-compileNimble(model1a) # compile the model

		#saveRDS(Cm1a,'nimblemodelCompiled_model1a_BSEM_resourceChp3.RDS')
		#	Cm1a<-readRDS('nimblemodelCompiled_model1a_BSEM_resourceChp3.RDS')

	##########################
	## Compile the MCMC, Run and Evaluate samples
	##########################

		C1a.MCMC.outputB<-nimbleMCMC(code=modelCode1a,data=myData,constants= myConstants, inits=init.values,niter=10000,nburnin=1000,nchains=4) #,summary=TRUE,samples=TRUE

		samples_chain1of4_C1aMCMC<-as.data.frame(C1a.MCMC.output$chain1)%>%mutate(chain='1')%>%mutate(iter=seq(1,4000,1))
		samples_chain2of4_C1aMCMC<-as.data.frame(C1a.MCMC.output$chain2)%>%mutate(chain='2')%>%mutate(iter=seq(1,4000,1))
		samples_chain3of4_C1aMCMC<-as.data.frame(C1a.MCMC.output$chain3)%>%mutate(chain='3')%>%mutate(iter=seq(1,4000,1))
		samples_chain4of4_C1aMCMC<-as.data.frame(C1a.MCMC.output$chain4)%>%mutate(chain='4')%>%mutate(iter=seq(1,4000,1))

		samples_fourchains_C1aMCMC<-bind_rows(samples_chain1of4_C1aMCMC,samples_chain2of4_C1aMCMC,samples_chain3of4_C1aMCMC,samples_chain4of4_C1aMCMC)%>%
			rename(a1='a[1]',a2='a[2]',a3='a[3]',a4='a[4]',a5='a[5]',c1='c[1]',c2='c[2]',c3='c[3]',c4='c[4]',c5='c[5]')

		head(samples_fourchains_C1aMCMC)


		## evaluate chain sampling

		paramnames<-colnames(samples_fourchains_C1aMCMC%>%group_by(chain))
		paramnames<-paramnames[1:24]
		plot_list = list()
		for(i in paramnames){
			var<-sym(i)
			p = ggplot(data=samples_fourchains_C1aMCMC%>%group_by(chain),aes(x=iter,col=chain,y=!!var))+
			geom_line()+
			ylab(NULL)+
			ggtitle(paste(var))+
			theme_bw()
			#print(p) # include if you want it to print each one.
			plot_list[[i]] = p
			filename = paste('fourchainsplot_model1a_nimbleBSEM_', var, '.png',sep="")
			png(filename, units='in', res=1100, width=7,height=3)
			print(plot_list[[i]])
			dev.off()
		}

	##########################
	## Posterior Inference 
	##########################

		# numeric summaries

		mcmc_summary_Cmodel1a<-MCMCsummary(C1a.MCMC.outputB,round=3)%>%
				tibble::rownames_to_column()%>%
				rename_with(str_to_title)%>%
				flextable()%>%
				theme_alafoli()%>%
				set_header_labels(rowname = 'Coefficient',SD='Sd')%>%
				align(align = 'center', part = 'all')%>%
				font(fontname = 'Times', part = 'all')%>%
				color(color='black',part='all')%>%
				fontsize(size = 10, part = 'all')%>%
				autofit()
			save_as_image(mcmc_summary_Cmodel1a,'mcmc_summary_4chains_5000iter_1000burnin_model1a_nimbleBSEM.png',webshot='webshot')

		# caterpillar plots 


		MCMCplot(C1a.MCMC.output,ci=c(50,95)) # point = median, thick line = 50% CI, thin line = 95% CI 
		MCMCsummary(C1a.MCMC.output,round=2)%>%
				tibble::rownames_to_column()%>%
				rename_with(str_to_title)%>%
				ggplot()+
				geom_pointrange(aes(x=Rowname,ymin='2.5%',ymax='97.5%'))

		ggsave(,file='mcmc_medianANDeffectTails_4chains_5000iter_1000burnin_model1a_nimbleBSEM.png',device=png,dpi=800, units='in',height=5,width=3)

		# trace and density plots

		coeffNintercept<-colnames(samples_fourchains_C1aMCMC,do.NULL=TRUE,prefix='row')


		MCMCtrace(C1a.MCMC.output,pdf=TRUE,ind=TRUE, Rhat=TRUE, n.eff=TRUE) # ind = TRUE, separate density lines per chain. # pdf = FALSE, don't export to a pdf automatically. 


		# DAG 
		#Cm1a$plotGraph()





##########################################
########### MODEL 2A, HEXAGON LEVEL MODEL

##########################

modelCode2a<-nimbleCode({
  
  ##########################
  ######### priors #########
  
  #### prior for sharkiness #### 
  for(i in 1:5){
    a[i] ~ dnorm(0,.001) 
	}
  # prior for intercept - sharks 
  b ~ dnorm(0,.001)
  # prior for residual variance - sharks 
  tau.shark ~ dgamma(0.01,0.01) # 
  sigma.shark <- sqrt(1/tau.shark)
  
  # group-level effect of buffID on sharkiness. prior for variance (epsi_site) calculated from inverse of precision (epsi_site)
  for(i in 1:B){
    epsi_shark[i]~ dnorm(0,tau.epsi_shark)
	}
  tau.epsi_shark ~ dunif(0,100)
  sigma.epsi_shark <- sqrt(1/tau.epsi_shark)
  
  
  #### prior for fishiness (cross metric, hexagons) ####
  for(i in 1:5){
	    c[i] ~ dnorm(0,0.001) 
		}
	  # prior for intercept - fishiness
	  d ~ dnorm(0,0.001)
	  
	  # prior for residual variance (precision) - fishiness
	  tau.fish ~ dgamma(0.01,0.01) 
	  sigma.fish <- sqrt(1/tau.fish) # gives sd 
  

  #### prior for seagrass PCA @ hexagon level ####
  for(i in 1:2){
	  	e[i] ~ dnorm(0,.001)
	  }

	  # prior for intercept - hexagon sg PCA
	  f ~ dnorm(0,.001)

	  # prior for residual variance (precision) - hexagon sg PCA 
	  tau.hexsgPCA ~ dgamma(.01,.01)
	  sigma.hexsgPCA <- sqrt(1/tau.hexsgPCA)

	  # prior for autocorrelation structure - hexagon sg PCA

  
  #### prior for seagrass PCA @ pointdata/receiver level ####
	for(i in 1:2){
		g[i] ~ dnorm(0,.001)
	}
	  # prior for intercept - point sg PCA
		h ~ dnorm(0,.001)
	  # prior for residual variance (precision) - point sg PCA 
		tau.sgPCA ~ dgamma(0.01,0.01)
	  sigma.sgPCA <- sqrt(1/tau.sgPCA)
  
  # group-level effect of buffID on sg PCA. prior for variance (epsi_site) calculated from inverse of precision (epsi_site)
  for(i in 1:B){
    epsi_sgPCA[i]~ dnorm(0,tau.epsi_sgPCA)
	}
  tau.epsi_sgPCA ~ dunif(0,100)
  sigma.epsi_sgPCA <- sqrt(1/tau.epsi_sgPCA)




  #########################################################
  ######### Likelihoods - data and process models #########	
  
  ### data model for sgPCA for pointdata/receiver level (sharks)
  for(i in 1:B){
    standard.sgPCA1[i] ~ dnorm(sgPCA.mu[i],tau.sgPCA)   
  # process model for point-sgPCA ~ dist2shore + dist_cmg + epsi_site
    sgPCA.mu[i] <- h + epsi_sgPCA[buffID[i]] + g[1]*standard.dist2shore[i] + g[2]*standard.distcmg[i] 
  }

  ### data model for sgPCA for hexdata/hexagon level (fishes) 
  for(i in 1:hex.N){
    standard.hexsgPCA1[i] ~ dnorm(hexsgPCA.mu[i],tau.hexsgPCA) 
  # process model for fishiness ~ shark + dist2shore + dist_cmg + prop_ldsg + prop_medsg  + epsi_fish
    hexsgPCA.mu[i] <- f + e[1]*standard.hexdist2shore[i] + e[2]*standard.hexdistcmg[i] 
	}


  ### data model for sharkiness 
  for(i in 1:point.N){
    standard.shark[i] ~ dnorm(shark.mu[i],tau.shark)   
  # process model for sharkiness ~ fish + dist2shore + dist_cmg + prop_ldsg + prop_medsg ? + epsi_site?
    shark.mu[i] <- b + epsi_shark[buffID[i]] + a[1]*standard.fish.pred[i] + a[2]*standard.dist2shore[i] + a[3]*standard.distcmg[i] + a[4]*standard.sgPCA1[i] # could be problematic that the value here is meant to come from data but is a estimated parameter above.
  }
  ### data model for fishiness 
  for(i in 1:hex.N){
    standard.hexfish[i] ~ dnorm(hexfish.mu[i],tau.fish) 
  # process model for fishiness ~ shark + dist2shore + dist_cmg + prop_ldsg + prop_medsg  + epsi_fish
    hexfish.mu[i] <- d + c[1]*standard.hexshark[i] + c[2]*standard.hexdist2shore[i] + c[3]*standard.hexdistcmg[i] + c[4]*standard.hexsgPCA1[i]
	}

  ######### derived parameters #########
  # for estmating total pathways 
  # write one for each route through the pathway diagram
  
  #path[1]<- g[1]*g[2] # point SG PCA ~ dist2shore + dist_cmg 
  #path[2]<- e[1]*e[2] # hex SG PCA ~ dist2shore + dist_cmg
  
  ## 26 May 2023: technically, paths 1 & 2 here should not be included. these paths should be paths that end at shark and fish. so it should actually be just the hab pathways that arrive at shark/fish, and the hab + shark/fish pathways that end at fish/shark respectively. Keeping for now, and not stopping the current MCMC run (model2a 200000iter, 25000burnin, 3 chains) because not sure they effect the runs of the other pathways but in next model run, remove. 

  path[1]<- g[1]*g[2]*a[4] # shark ~ dist2shore + dist_cmg + sgPCA
  path[2]<- e[1]*e[2]*c[4] # fish ~ dist2shore + dist_cmg + sgPCA

  path[3]<- c[1] # fish ~ shark 
  path[4]<- a[1] # shark ~ fish

  path[5]<- g[1]*g[2]*a[4]*a[1] # shark ~ fish + dist2shore + dist_cmg + prop_ldsg + prop_medsg # with coeff for habitat as they went into fish and came through
  path[6]<- e[1]*e[2]*c[4]*c[1] # fish ~ shark + dist2shore + dist_cmg + hexsgPCA: e2+3 are the effect of habtat on shark VIA their effect on seagrass, and c4 is the effect of seagrass, and c1 is the effect of sharks, on fish.    

  
  
}) # end of model code 


##########################
## Compile the model code
##########################

	myConstants2<-list(point.N=560,hex.N=nrow(hexdata),B=max(pointdata$buffIDnum),buffID=pointdata$buffIDnum) 
	myData2<-list(
	  # tell nimble the covariates 
	  	# point data
	  standard.shark = pointdata$standard.shark,
	  standard.dist2shore = pointdata$standard.dist2shore,
	  standard.distcmg = pointdata$standard.distcmg,
	  #standard.lds = pointdata$standard.lds,
	  #standard.mds = pointdata$standard.mds,
	  standard.fish.pred = pointdata$standard.fish,
	  standard.sgPCA1 = pointdata$standard.sgPCA1,
	  	# hex data 
	  standard.hexfish = hexdata$standard.hexfish,
	  standard.hexshark = hexdata$standard.hexshark,
	  standard.hexdist2shore = hexdata$standard.hexdist2shore,
	  standard.hexdistcmg = hexdata$standard.hexdistcmg,
	  #standard.hexlds = hexdata$standard.hexlds,
	  #standard.hexmds = hexdata$standard.hexmds,
	  standard.hexsgPCA1 = hexdata$standard.hexsgPCA1
	)

	init.values2<-list(a=rnorm(5,0,1),b=rnorm(1),c=rnorm(5,0,1),d=rnorm(1),e=rnorm(2,0,1),f=rnorm(1),g=rnorm(2,0,1),h=rnorm(1),tau.shark=rgamma(1,0.01,0.01),tau.fish=rgamma(1,0.01,0.01),tau.epsi_shark=runif(1,0,100),tau.epsi_sgPCA=runif(1,0,100))
		
	model2a<-nimbleModel(code=modelCode2a, name="model2a",data=myData2,constants = myConstants2,inits=init.values2) #define the model

	Cm2a<-compileNimble(model2a) # compile the model

  # saveRDS(Cm2a,'nimblemodelCompiled_model2a_BSEM_resourceChp3.RDS')
	#	Cm2a<-readRDS('nimblemodelCompiled_model2a_BSEM_resourceChp3.RDS')

##########################
## Compile & Run the MCMC
##########################

	C2a.MCMC.output<-nimbleMCMC(code=modelCode2a,data=myData2,constants= myConstants2, inits=init.values2,niter=50000,nburnin=2000,nchains=3,monitors=c('a','b','c','d','e','f','g','h','path','tau.epsi_shark','tau.fish','tau.shark','tau.epsi_sgPCA')) 

##########################
## Posterior Inference 
##########################

	# numeric summaries

	mcmc_summary_Cmodel2a<-MCMCsummary(C2a.MCMC.output,round=3,pg0=TRUE)%>%
			tibble::rownames_to_column()%>%
			rename_with(str_to_title)%>%
			flextable()%>%
			theme_alafoli()%>%
			set_header_labels(rowname = 'Coefficient',SD='Sd')%>%
			align(align = 'center', part = 'all')%>%
			font(fontname = 'Times', part = 'all')%>%
			color(color='black',part='all')%>%
			fontsize(size = 10, part = 'all')%>%
			autofit()

		save_as_image(mcmc_summary_Cmodel2a,'mcmc_summary_4chains_5000iter_1000burnin_model2a_nimbleBSEM.png',webshot='webshot')

	# caterpillar plots 


	MCMCplot(C2a.MCMC.output,ci=c(50,95),params=c('a','c','e','g','path')) # point = median, thick line = 50% CI, thin line = 95% CI 



	# trace and density plots

	coeffNintercept<-colnames(samples_fourchains_C1aMCMC,do.NULL=TRUE,prefix='row')


	MCMCtrace(C2a.MCMC.outputB,pdf=TRUE,ind=TRUE, Rhat=TRUE, n.eff=TRUE) # ind = TRUE, separate density lines per chain. # pdf = FALSE, don't export to a pdf automatically. 






##############################################################
######################### code graveyard #####################
##############################################################


	########### This is before I remembered I need to mean centre and standardse the variables

		##MODEL 1B, RECEIVER LEVEL MODEL: shark point detectons (counts per receiver for each individual, sample size N=35), Gerreidae abundance predictions (BRT simplified model), habitat 
			## for habitat variables, we used the types identified in previous modelling as influencing habitat selection. For sharks, this is distance to the central mangrove point and medium density seagrass; for Gerreidae, this is distance to shore and low density seagrass; and for the diversity index for species, we included only the habitats with relative influence in the BRTs > 15%, which are distance to the shore and low density seagrass. 
			## random effect for SITE but not for shark individual

			## 15/5: fishy variable with gamma distribution; 'n' for sharks with a neg-binomial distribution

		# arguments to pass later to nimble in nimbleModel()

			myCovs<-c('n','buffID','fishy','dist2shore','dist_cmg','prop_ldsg','prop_medsg')
			myConstants<-list(n.obsv=560,numGroups=35)
			myData<-list(y=pointdata,x=myCovs)
			myInits<-list(alpha = ,beta = ,theta = )

			nimbleModel(constants=myConstants,data=myData, inits=list()) ## need to define myCovs

		model1<-nimblecode({})

		##########################
		######### priors #########

		#### prior for sharkiness #### 
		for(i in 1:3){
			a[i] ~ dnorm(0,.001) # trying negative binomial for sharks 'n', dnegbin(prob = p, size = r), a mean and a dispersion
		}
		# prior for intercept - sharks 
			b ~ dnorm(0,.001)
		# prior for residual variance - sharks 
			tau.shark ~ exp(1) # 
			sigma.shark <- 1/tau.shark # dispersion 


		#### prior for fishiness (Gerreidae) ####
		for(i in 1:3){
			c[i] ~ dgamma(0.01,0.01) 
		}
		# prior for intercept - fishiness
			d ~ dgamma(0.01,0.01)

		# prior for residual variance (precision) - fishiness
		tau.fish ~ dinvgamma(.01,.01) # inverse gamma? 
		sigma.fish ~ 1/tau.fish^2 #(shape/rate^2)

		# group-level effect for site ID on fishiness (epsi_fish)
		# prior for precision of random effect (buffID), then calculated to variance 
		for(i in 1:n.obsv){
			epsi_fish[i] ~ dnorm(0,tau.epsi_fish)
		}
		tau.epsi_fish ~ dgamma(.001,.001)
		sigma.epsi_fish ~ 1/tau.epsi_fish


		#### prior for habitats: distance to shore, distance to refuge, low density seagrass, medium density seagrass ####

		# priors - habitat: dist2shore
			e ~ dbeta(0,.1) # slope
			ee ~ dbeta(0,.1) # intercept

		# priors - habitat: dist_cmg
			f ~ dbeta(0,.1) # slope
			ff ~ dbeta(0,.1) # intercept

		# priors - habitat: prop_ldsg
			g ~ dbeta(0,.1) # slope
			gg ~ dbeta(0,.1) # intercept

		# priors - habitat: prop_medsg
			h ~ dbeta(0,.1) # slope
			hh ~ dbeta(0,.1) # intercept


		# group-level effect of buffID on habitat values. prior for variance (epsi_site) calculated from inverse of precision (epsi_site)
		for(i in 1:n.obsv){
			epsi_site[i]~dnorm(0,tau.epsi_site)
		}
		tau.epsi_site ~ dgamma(.001,.001)
		sigma.epsi_site ~ 1/tau.epsi_site

		# prior for residual variance - site of hab sample
		tau.site ~ dgamma(.001,.001)
		sigma.site ~ 1/tau.site


		###########################################
		######### data and process models #########	

		### data model for sharkiness 
			for(i in 1:n.obsv){
				shark[y] ~ dnbinom(shark.mu[i],sigma.shark) # should this be tau.shark or sigma.shark? tau is the dispersion or the precision? blogs say the dispersion
				# dnbinom, needs a mean (mu) and a dispersion (sigma?tau?) 
				# lost in the sauce
			}

		# process model for sharkiness ~ fish + dist2shore + dist_cmg + prop_ldsg + prop_medsg ? + epsi_site?
			for(i in 1:n.obsv){
				shark.mu[y] ~ b + a[1]*gerries[y] + a[2]dist2shore[y] + a[3]dist_cmg[y] + a[4]prop_ldsg[y] + a[5]prop_medsg[y] + epsi_site[buffID[y]] # ?? also epsi_fish? no, only in the one that specifies fish (fish on the LHS). 
			}

		### data model for fishiness 
			for(i in 1:n.obsv){
				fish[y] ~ dgamma(fish.shape[i],tau.fish) # tau.fish is actually the rate (?..i think) currently, its got a gamma distribution, and is used to calculate the precision 
			}
		# process model for fishiness ~ shark + dist2shore + dist_cmg + prop_ldsg + prop_medsg  + epsi_fish
			for(i in 1:n.obsv){
				fish.shape[y] ~ d + c[1]*n[y] + c[2]dist2shore[y] + c[3]dist_cmg[y] + c[4]prop_ldsg[y] + c[5]prop_medsg[y] + epsi_fish[buffID[y]]
			}


		######### derived parameters #########
		# for estmating total pathways 
		# write one for each route through the pathway diagram

			habs: dist2shore + dist_cmg + prop_ldsg + prop_medsg 
		 
		### you left off here on Friday (12/5), about to re-work this section so it is correct and updated with the multiple hab types and their coefficients. 

		path[1]<-  # shark ~ dist2shore + dist_cmg + prop_ldsg + prop_medsg 
		path[2]<-  # fish ~ dist2shore + dist_cmg + prop_ldsg + prop_medsg + epsi_site
		path[3]<-  # shark ~ fish 
		path[4]<-  # fish ~ shark + epsi_site
		path[5]<- # shark ~ fish + dist2shore + dist_cmg + prop_ldsg + prop_medsg 
		path[6]<- # fish ~ shark + dist2shore + dist_cmg + prop_ldsg + prop_medsg + epsi_site

		# ^^ some should have ntercepts? when? why?


		######################################################
		#### Code for model with individual random effect ####

		## random effect for individual on sharks 
		# prior for precision of random effect (individual), then calculated to variance 
			for(i in 1:PIT){
				epsi_shark[i] ~ dnorm(0,tau.epsi_shark)
			}
			tau.epsi_shark ~ dgamma(.001,.001)
			sigma.epsi_shark ~ 1/tau.epsi_shark

			# prior for residual variance - sharkiness 
			tau.shark ~ dgamma(.001,.001)
			sigma.shark <- 1/tau.shark


