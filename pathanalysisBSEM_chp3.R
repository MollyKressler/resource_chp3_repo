## Path analysis for fish abundance, shark presence, and habitat composition, on each other.

## Modelling code from paper Bennett, S., Harris, M. P., Wanless, S., Green, J. A., Newell, M. A., Searle, K. R., & Daunt, F. (2022). Earlier and more frequent occupation of breeding sites during the non‐breeding season increases breeding success in a colonial seabird. Ecology and Evolution, 12(9). https://doi.org/10.1002/ece3.9213.

## created 10 May 2023 by Molly M Kressler 


#########################################################################
######################### RUN AT OPEN ###################################

## Load Workspace 

pacman::p_load(tidyverse,sf,nimble)

#setwd('/Users/mollykressler/Documents/data_phd/resource_chp3')

data<-read.csv('data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_may23.csv')#%>%rename('prop_brs'='prp_','prop_ldsg'='prp_lds','prop_medsg'='prp_mds','prop_hdsg'='prp_hds','prop_sarg'='prp_srg','prop_marsh'='prp_mrs','prop_urb_r'='prp_rb_','prop_deep'='prop_dp','prop_spong'='prp_spn','prop_unkwn'='prp_nkw')
data$PIT<-as.character(data$PIT)
data$buffID<-as.character(data$buffID)
	
stopifnot(nrow(data)==33701) # check 

################################
###### Define Sharkiness  ######
## counts of detections per individual per Site (buffer/receiver)

habitats<-c('prop_brs','prop_ldsg','prop_medsg','prop_hdsg','prop_sarg','prop_marsh','prop_urb_r','prop_deep','prop_unkwn')
fishmetrics<-c('bt_g_PD','sppsmPD')

fishANDhab<-data%>%
	dplyr::select(buffID,dist2shore,dst_cmg,unlist(habitats),unlist(fishmetrics))%>%
	group_by(buffID)%>%
	filter(row_number()==1)

data2<-data%>%
	dplyr::select(PIT,buffID)%>%
	group_by(PIT,buffID)%>%
	tally()
data3<-expand.grid(PIT=unique(data2$PIT),buffID=unique(data2$buffID))%>%
	left_join(.,data2)%>%
	mutate(n=replace_na(n,0))%>%
	left_join(.,fishANDhab,by='buffID')%>%
	mutate(gerries=bt_g_PD)%>%
	mutate(fishy=sppsmPD*gerries)%>% 
	mutate(standard.shark=((n-mean(n))/sd(n)))%>% 
	mutate(standard.fish=((fishy-mean(fishy))/sd(fishy)))%>%
	mutate(standard.lds=((prop_ldsg-mean(prop_ldsg))/sd(prop_ldsg)))%>% 
	mutate(standard.mds=((prop_medsg-mean(prop_medsg))/sd(prop_medsg)))%>% 
	mutate(standard.unkn=((prop_unkwn-mean(prop_unkwn))/sd(prop_unkwn)))%>% 
	mutate(standard.dist2shore=((dist2shore-mean(dist2shore))/sd(dist2shore)))%>% 
	mutate(standard.distcmg=((dst_cmg-mean(dst_cmg))/sd(dst_cmg)))%>%
	dplyr::select(PIT,buffID,standard.shark,standard.fish,standard.dist2shore,standard.distcmg,standard.lds,standard.mds)

data3$PIT<-as.factor(data3$PIT)
data3$buffID<-as.factor(data3$buffID)

stopifnot(!is.na(data3$standard.lds) & !is.na(data3$prop_brs) & nrow(data3)==560)

summary(data3)
head(data3,5)


######################### END ###################################
#################################################################


#####################################################################
######################### CHECKS  ###################################


###
## Compare Fishy (cross metric) and Diversity and Gerreidae

	ggplot(data=data3)+geom_histogram(aes(x=gerries),fill='violetred2',alpha=0.5, binwidth=5)+geom_histogram(aes(x=fishy),fill='cadetblue3',alpha=0.5, binwidth=5)+theme_bw()+xlab('Value')+ylab('Count')+theme(legend.position='topright')

	ggplot(data=data3)+geom_histogram(aes(x=standard.fish),fill='goldenrod3',alpha=0.5, binwidth=.5)+theme_bw()

###
## DISTRIBUTION CHECKS 
	## for standardised fish cross metric variable 'fishy' (count of gerries x shannon index for species) 
	summary(data3$standard.fish) # min=-0.83, max = 3.395

		ggplot(data=data3)+geom_histogram(aes(x=(standard.fish)),fill='cadetblue4',alpha=0.8,binwidth=0.5)+theme_bw()
		#normal I guess


	## for standardised shark variable 'n' (sum(detections) per individual @ each receiver)

		ggplot(data=data3)+geom_histogram(aes(x=standard.shark),fill='cadetblue4',alpha=0.5,binwidth=.5)+theme_bw()

######################### END ###################################
#################################################################


########################################################
################# BUGS models with NIMBLE ##############


##########################################
########### MODEL 1A, RECEIVER LEVEL MODEL: MEAN CENTRED AND STANDARDISED DATA - shark point detectons (counts per receiver for each individual, sample size N=35), Gerreidae abundance predictions (BRT simplified model), habitat 
	## for habitat variables, we used the types identified in previous modelling as influencing habitat selection. For sharks, this is distance to the central mangrove point and medium density seagrass; for Gerreidae, this is distance to shore and low density seagrass; and for the diversity index for species, we included only the habitats with relative influence in the BRTs > 15%, which are distance to the shore and low density seagrass. 
	## random effect for SITE but not for shark individual

# all vars as normal, because they are + and - and centred on 0. 

# arguments to pass later to nimble in nimbleModel()

	myCovs<-c('n','buffID','fishy','dist2shore','dist_cmg','prop_ldsg','prop_medsg')
	myConstants<-list(N=560,numGroups=35)
	myData<-list(y=data3,x=myCovs)
	myInits<-list(alpha = ,beta = ,theta = )

	nimbleModel(constants=myConstants,data=myData, inits=list()) ## need to define myCovs

##########################

modelCode1a<-nimbleCode({

##########################
######### priors #########

#### prior for sharkiness #### 
for(i in 1:5){
	a[i] ~ dnorm(0,.001) }
# prior for intercept - sharks 
	b ~ dnorm(0,.001)
# prior for residual variance - sharks 
	tau.shark ~ dgamma(0.001,0.001) # 
	sigma.shark <- 1/tau.shark 

# group-level effect of buffID on sharkiness. prior for variance (epsi_site) calculated from inverse of precision (epsi_site)
for(i in 1:N){
	epsi_shark[i]~ dnorm(0,tau.epsi_shark)}
tau.epsi_shark ~ dgamma(.001,.001)
sigma.epsi_shark ~ 1/tau.epsi_shark


#### prior for fishiness (Gerreidae) ####
for(i in 1:5){
	c[i] ~ dnorm(0,0.001) }
# prior for intercept - fishiness
	d ~ dnorm(0,0.001)

# prior for residual variance (precision) - fishiness
tau.fish ~ dgamma(0.001,0.001) 
sigma.fish ~ 1/tau.fish 

# group-level effect for site ID on fishiness (epsi_fish)
# prior for precision of random effect (buffID), then calculated to variance 
for(i in 1:N){ #n.obsv? or should it related to the grouping
	epsi_fish[i] ~ dnorm(0,tau.epsi_fish)}
tau.epsi_fish ~ dgamma(0.001,0.001)
sigma.epsi_fish ~ 1/tau.epsi_fish


#### prior for habitats: distance to shore, distance to refuge, low density seagrass, medium density seagrass ####
 # not even sure I need to specify these since they aren't on the LHS
# priors - habitat: dist2shore
	e ~ dnorm(0.001,0.001) # slope
	ee ~ dnorm(0.001,0.001) # intercept

# priors - habitat: dist_cmg
	f ~ dnorm(0.001,0.001) # slope
	ff ~ dnorm(0.001,0.001) # intercept

# priors - habitat: prop_ldsg
	g ~ dnorm(0.001,0.001) # slope
	gg ~ dnorm(0.001,0.001) # intercept

# priors - habitat: prop_medsg
	h ~ dnorm(0.001,0.001) # slope
	hh ~ dnorm(0.001,0.001) # intercept



#########################################################
######### Likelihoods - data and process models #########	

### data model for sharkiness 
	for(i in 1:N){
		shark[i] ~ dnorm(shark.mu[i],tau.shark) }

# process model for sharkiness ~ fish + dist2shore + dist_cmg + prop_ldsg + prop_medsg ? + epsi_site?
	for(i in 1:N){
		shark.mu[i] ~ b + a[1]*standard.fish[i] + a[2]*standard.dist2shore[i] + a[3]*standard.distcmg[i] + a[4]*standard.lds[i] + a[5]*standard.mds[i] + epsi_shark[buffID[i]]}

### data model for fishiness 
	for(i in 1:N){
		fish[i] ~ dnorm(fish.mu[i],tau.fish) }
# process model for fishiness ~ shark + dist2shore + dist_cmg + prop_ldsg + prop_medsg  + epsi_fish
	for(i in 1:N){
		fish.mu[i] ~ d + c[1]*standard.shark[i] + c[2]*standard.dist2shore[i] + c[3]*standard.distcmg[i] + c[4]*standard.lds[i] + c[5]*standard.mds[i] + epsi_fish[buffID[i]]}


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

		# alternate derived paths 
			path[1]<- e*f*g*h*a[2]*a[3]*a[4]*a[5] 
				# shark ~ dist2shore + dist_cmg + prop_ldsg + prop_medsg 
			path[2]<- e*f*g*h*c[2]*c[3]*c[4]*c[5] 
				# fish ~ dist2shore + dist_cmg + prop_ldsg + prop_medsg 
			path[3]<- a[1] 
				# shark ~ fish 
			path[4]<- c[1] 
				# fish ~ shark 
			path[5]<- e*f*g*h*a[1]*c[2]*c[3]*c[4]*c[5] 
				# shark ~ fish + dist2shore + dist_cmg + prop_ldsg + prop_medsg # with coeff for habitat as they went into fish and came through
			path[6]<- e*f*g*h*c[1]*a[2]*a[3]*a[4]*a[5]
				# fish ~ shark + dist2shore + dist_cmg + prop_ldsg + prop_medsg 


##########################
## Compile the model code
##########################


	myCovs<-c('standard.fish','standard.shark','standard.dist2shore','standard.distcmg','standard.lds','standard.mds','buffID')
	myConstants<-list(N=560,buffID=35) # if above, for group-level effects, the 'N' should be related to the numbe of groups, then the 'N' should be changed to numGroups, or numGroups changed to buffID. 
	myData<-list(
		# tell nimble the covariates 
		standard.fish = data3$standard.fish,
		standard.shark = data3$standard.shark,
		standard.dist2shore = data3$standard.dist2shore,
		standard.distcmg = data3$standard.distcmg,
		standard.lds = data3$standard.lds,
		standard.mds = data3$standard.mds,
		# tell nimble the constants
		N=560,
		buffID=as.numeric(factor(data3$buffID))
		)
	initial.values <- function() list(theta = runif(1,0,1))
	init.values<-list(alpha=1,beta=1,theta= rep(0.1,myConstants$N))

model1a<-nimbleModel(code=modelCode1a, name="model1a",data=myData)
	
	model1a$getNodeNames

model1aMCMC.output<-nimbleMCMC(code=modelCode1a,data=myData,inits=initial.values,niter=5000,nburnin=1000,nchains=1)

##########################
## Compile the MCMC
##########################
##########################
## Run the MCMC
##########################
##########################
## Evaluate the samples from the mcmc 
##########################





######################### END ###################################
#################################################################













##########################################
##########################################
##########################################
########### This is before I remembered I need to mean centre and standardse the variables

	##MODEL 1B, RECEIVER LEVEL MODEL: shark point detectons (counts per receiver for each individual, sample size N=35), Gerreidae abundance predictions (BRT simplified model), habitat 
		## for habitat variables, we used the types identified in previous modelling as influencing habitat selection. For sharks, this is distance to the central mangrove point and medium density seagrass; for Gerreidae, this is distance to shore and low density seagrass; and for the diversity index for species, we included only the habitats with relative influence in the BRTs > 15%, which are distance to the shore and low density seagrass. 
		## random effect for SITE but not for shark individual

		## 15/5: fishy variable with gamma distribution; 'n' for sharks with a neg-binomial distribution

	# arguments to pass later to nimble in nimbleModel()

		myCovs<-c('n','buffID','fishy','dist2shore','dist_cmg','prop_ldsg','prop_medsg')
		myConstants<-list(n.obsv=560,numGroups=35)
		myData<-list(y=data3,x=myCovs)
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


