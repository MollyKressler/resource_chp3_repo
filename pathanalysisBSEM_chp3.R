## Path analysis for fish abundance, shark presence, and habitat composition, on each other.

## Modelling code from paper Bennett, S., Harris, M. P., Wanless, S., Green, J. A., Newell, M. A., Searle, K. R., & Daunt, F. (2022). Earlier and more frequent occupation of breeding sites during the non‐breeding season increases breeding success in a colonial seabird. Ecology and Evolution, 12(9). https://doi.org/10.1002/ece3.9213.

## created 10 May 2023 by Molly M Kressler 


#########################################################################
######################### RUN AT OPEN ###################################

## Load Workspace 

pacman::p_load(tidyverse,sf,nimble,MCMCvis,flextable,webshot2)

setwd('/Users/mollykressler/Documents/data_phd/resource_chp3')

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
	mutate(buffIDnum=str_remove(data3$buffID,"r"))%>%
	dplyr::select(PIT,buffIDnum,standard.shark,standard.fish,standard.dist2shore,standard.distcmg,standard.lds,standard.mds)

data3$buffIDnum<-as.numeric(data3$buffIDnum)

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
  tau.shark ~ dgamma(0.001,0.001) # 
  sigma.shark <- sqrt(1/tau.shark)
  
  # group-level effect of buffID on sharkiness. prior for variance (epsi_site) calculated from inverse of precision (epsi_site)
  for(i in 1:B){
    epsi_shark[i]~ dnorm(0,tau.epsi_shark)
}
  tau.epsi_shark ~ dgamma(.001,.001)
  sigma.epsi_shark <- sqrt(1/tau.epsi_shark)
  
  
  #### prior for fishiness (Gerreidae) ####
  for(i in 1:5){
    c[i] ~ dnorm(0,0.001) 
}
  # prior for intercept - fishiness
  d ~ dnorm(0,0.001)
  
  # prior for residual variance (precision) - fishiness
  tau.fish ~ dgamma(0.001,0.001) 
  sigma.fish <- sqrt(1/tau.fish) # gives sd 
  
  # group-level effect for site ID on fishiness (epsi_fish)
  # prior for precision of random effect (buffID), then calculated to variance 
  for(i in 1:B){ # B is no. of buffers
    epsi_fish[i] ~ dnorm(0,tau.epsi_fish)
}
  tau.epsi_fish ~ dgamma(0.001,0.001)
  sigma.epsi_fish <- sqrt(1/tau.epsi_fish)
  
  #########################################################
  ######### Likelihoods - data and process models #########	
  
  ### data model for sharkiness 
  for(i in 1:N){
    standard.shark[i] ~ dnorm(shark.mu[i],tau.shark)   
  # process model for sharkiness ~ fish + dist2shore + dist_cmg + prop_ldsg + prop_medsg ? + epsi_site?
    shark.mu[i] <- b + epsi_shark[buffID[i]] #+ a[1]*standard.fish.pred[i] + a[2]*standard.dist2shore[i] + a[3]*standard.distcmg[i] + a[4]*standard.lds[i] + a[5]*standard.mds[i] 
  ### data model for fishiness 
    standard.fish[i] ~ dnorm(fish.mu[i],tau.fish) 
  # process model for fishiness ~ shark + dist2shore + dist_cmg + prop_ldsg + prop_medsg  + epsi_fish
    fish.mu[i] <- d + epsi_fish[buffID[i]]#+ c[1]*standard.shark.pred[i] + c[2]*standard.dist2shore[i] + c[3]*standard.distcmg[i] + c[4]*standard.lds[i] + c[5]*standard.mds[i] 
}

  ######### derived parameters #########
  # for estmating total pathways 
  # write one for each route through the pathway diagram
  
  path[1]<- b*a[2]*a[3]*a[4]*a[5] 
  # shark ~ dist2shore + dist_cmg + prop_ldsg + prop_medsg 
  path[2]<- d*c[2]*c[3]*c[4]*c[5] 
  # fish ~ dist2shore + dist_cmg + prop_ldsg + prop_medsg 
  path[3]<- b*a[1] 
  # shark ~ fish 
  path[4]<- d*c[1] 
  # fish ~ shark 
  path[5]<- b*a[1]*c[2]*c[3]*c[4]*c[5] 
  # shark ~ fish + dist2shore + dist_cmg + prop_ldsg + prop_medsg # with coeff for habitat as they went into fish and came through
  path[6]<- d*c[1]*a[2]*a[3]*a[4]*a[5]
  # fish ~ shark + dist2shore + dist_cmg + prop_ldsg + prop_medsg 
  
  
}) # end of model code 

##########################
## Compile the model code
##########################

myCovs<-c('standard.fish','standard.shark','standard.dist2shore','standard.distcmg','standard.lds','standard.mds','buffID')
myConstants<-list(N=560,B=max(data3$buffIDnum),buffID=data3$buffIDnum) 
myData<-list(
  # tell nimble the covariates 
  standard.fish = data3$standard.fish,
  standard.shark = data3$standard.shark,
  standard.dist2shore = data3$standard.dist2shore,
  standard.distcmg = data3$standard.distcmg,
  standard.lds = data3$standard.lds,
  standard.mds = data3$standard.mds,
  standard.shark.pred = data3$standard.shark,
  standard.fish.pred = data3$standard.fish
  # tell nimble the constants
  #N=560,
  #buffID=as.numeric(factor(data3$buffID))
)

init.values<-list(a=rnorm(5,0,1),b=rnorm(1),c=rnorm(5,0,1),d=rnorm(1),tau.shark=rgamma(1,0.001,0.001),tau.fish=rgamma(1,0.001,0.001))
	
model1a<-nimbleModel(code=modelCode1a, name="model1a",data=myData,constants = myConstants,inits=init.values) #define the model
Cm1a<-compileNimble(model1a) # compile the model

#saveRDS(Cm1a,'nimblemodelCompiled_model1a_BSEM_resourceChp3.RDS')
#	Cm1a<-readRDS('nimblemodelCompiled_model1a_BSEM_resourceChp3.RDS')

##########################
## Compile the MCMC, Run and Evaluate samples
##########################

C1a.MCMC.output<-nimbleMCMC(code=modelCode1a,data=myData,constants= myConstants, inits=init.values,niter=5000,nburnin=1000,nchains=4) #,summary=TRUE,samples=TRUE

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

mcmc_summary_Cmodel1a<-MCMCsummary(C1a.MCMC.output,round=2)%>%
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














# ? calculate the path estimates - take the estimates for each coeff for each chain, then multiply across the path

a1<-c(C1a.MCMC.output$chain1[,'a[1]'],
	C1a.MCMC.output$chain2[,'a[1]'],
	C1a.MCMC.output$chain3[,'a[1]'],
	C1a.MCMC.output$chain4[,'a[1]'])
a2<-c(C1a.MCMC.output$chain1[,'a[2]'],
	C1a.MCMC.output$chain2[,'a[2]'],
	C1a.MCMC.output$chain3[,'a[2]'],
	C1a.MCMC.output$chain4[,'a[2]'])
a3<-c(C1a.MCMC.output$chain1[,'a[3]'],
	C1a.MCMC.output$chain2[,'a[3]'],
	C1a.MCMC.output$chain3[,'a[3]'],
	C1a.MCMC.output$chain4[,'a[3]'])
a4<-c(C1a.MCMC.output$chain1[,'a[4]'],
	C1a.MCMC.output$chain2[,'a[4]'],
	C1a.MCMC.output$chain3[,'a[4]'],
	C1a.MCMC.output$chain4[,'a[4]'])
b<-c(C1a.MCMC.output$chain1[,'b'],
	C1a.MCMC.output$chain2[,'b'],
	C1a.MCMC.output$chain3[,'b'],
	C1a.MCMC.output$chain4[,'b'])
path1<-a1*a2*a3*a4*b # gives a massive number..not right. missing something
mean(path1)
quantile(path1)
path1%>%as_tibble()%>%ggplot()+geom_histogram(aes(x=value))+theme_bw()

### 18/5 you left off here, you want to make a flextable of chains summary table rounded to 2 decimals places. you need to refresh on what the nEff is. additionally, just go through the whole Gimenez tutorial. this should help and be the figures/tables that rich is expecting. 













######################### END ###################################
#################################################################


  ### data model for sharkiness 
  for(i in 1:N){
    standard.shark[i] ~ dnorm(shark.mu[i],tau.shark) }
  
  # process model for sharkiness ~ fish + dist2shore + dist_cmg + prop_ldsg + prop_medsg ? + epsi_site?
  for(i in 1:N){
    shark.mu[i] <- b + a[1]*standard.fish[i] + a[2]*standard.dist2shore[i] + a[3]*standard.distcmg[i] + a[4]*standard.lds[i] + a[5]*standard.mds[i] + epsi_shark[buffID[i]]}
  
  ### data model for fishiness 
  for(i in 1:N){
    standard.fish[i] ~ dnorm(fish.mu[i],tau.fish) }
  # process model for fishiness ~ shark + dist2shore + dist_cmg + prop_ldsg + prop_medsg  + epsi_fish
  for(i in 1:N){
    fish.mu[i] <- d + c[1]*standard.shark[i] + c[2]*standard.dist2shore[i] + c[3]*standard.distcmg[i] + c[4]*standard.lds[i] + c[5]*standard.mds[i] + epsi_fish[buffID[i]]}
  










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


