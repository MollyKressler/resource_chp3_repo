## Path analysis for directionality of fish abundance, shark presence, and habitat composition, on each other.

## Modelling code from paper Bennett, S., Harris, M. P., Wanless, S., Green, J. A., Newell, M. A., Searle, K. R., & Daunt, F. (2022). Earlier and more frequent occupation of breeding sites during the non‐breeding season increases breeding success in a colonial seabird. Ecology and Evolution, 12(9). https://doi.org/10.1002/ece3.9213.

## created 10 May 2023 by Molly M Kressler 

## Load Workspace 

pacman::p_load(tidyverse,sf,nimble)

setwd('/Users/mollykressler/Documents/data_phd/resource_chp3')

not.needed.vars<-c("difftim","sgmnt_s","sgmnt_n","SgmntID","ghost" , "FID","station")
data<-read.csv('data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_may23.csv')%>%dplyr::select(-not.needed.vars)
data$PIT<-as.factor(data$PIT)
	
stopifnot(nrow(data)==33701) # check 
#############################################
###### Write model in BUGS with NIMBLE ######

# parts of nimble model: data, params, distirbutiions, logical operators/process models, derived params/generated quantities (optional)

### variables 
# sharkiness = detections (points), presence.
# fishiness = either Gerreidae abundance or Shannon Index of Species, as predicted by simplified BRTs.
# habitat = 12 proportions of habitat types.
# individual = PIT, factor in data

########### MODEL 1: shark point detectons, Gerreidae abundance predictions (BRT simplified model), habitat (12 variables, proportions of habitat)

model1<-nimblecode({})

##########################
######### priors #########

#### prior for sharkiness ####
for(i in 1:3){
	a[i] ~ dnorm(0,.0001)
}
# prior for intercept - sharks 
	b ~ dnorm(0,.0001)

## random effect for individual on sharks 
# prior for precision of random effect (individual), then calculated to variance 
for(i in 1:individuals){
	epsi_shark[i] ~ dnorm(0,tau.epsi_shark)
}
tau.epsi_shark ~ dgamma(.001,.001)
sigma.epsi_shark ~ 1/tau.epsi_shark

# prior for residual variance - sharkiness 
tau.shark ~ dgamma(.001,.001)
sigma.shark <- 1/tau.shark


#### prior for fishiness ####
for(i in 1:3){
	c[i] ~ dnorm(0,.0001)
}
# prior for intercept - fishiness
	d ~ dnorm(0,.0001)

# random effect for site ID on fishiness 
# prior for precision of random effect (buffID), then calculated to variance 
for(i in 1:sites){
	epsi_fish[i] ~ dnorm(0,tau.epsi_fish)
}
tau.epsi_fish ~ dgamma(.001,.001)
sigma.epsi_fish ~ 1/tau.epsi_fish

# prior for residual variance - fishiness
tau.fish ~ dgamma(.001,.001)
sigma.fish ~ 1/tau.fish


#### prior for habitat ####

# prior for intercept - habitat
e ~ dnorm(0,0.0001)

# prior for slope - habitat 
f ~ dnorm(0,0.0001)

# random effect of buffID on habitat values. prior for variance (sigma) calculated from inverse of precision (tau)
for(i in 1:sites){
	epsi.site[i]~dnorm(0,tau.epsi_site)
}
tau.epsi_site ~ dgamma(.001,.001)
sigma.epsi_site ~ 1/tau.epsi_site

# prior for residual variance - site of hab sample
tau.site ~ dgamma(.001,.001)
sigma.site ~ 1/tau.site


###########################################
######### data and process models #########	

### data model for sharkiness ~ hab + fish + indiv

# process model for sharkiness ~ hab + fish + indiv


### data model for fishiness ~ shark+ hab + indiv 

# process model for fishiness ~ shark + hab + indiv


######### derived parameters #########
# for estmating total pathways 
# write one for each route through the pathway diagram

path[1]<- # shark ~ hab + indiv
path[2]<- # fish ~ hab 
path[3]<- # shark ~ fish + indiv
path[4]<- # fish ~ shark
path[5]<- # shark ~ fish + indiv + hab
path[6]<- # fish ~ shark + hab 






