 ## Path analysis for Risk & Resource BSEM (chapter 3)

## Code and approach here heavily inspired by modelling code from paper Bennett, S., Harris, M. P., Wanless, S., Green, J. A., Newell, M. A., Searle, K. R., & Daunt, F. (2022). Earlier and more frequent occupation of breeding sites during the non‐breeding season increases breeding success in a colonial seabird. Ecology and Evolution, 12(9). https://doi.org/10.1002/ece3.9213.

## created 10 May 2023 :: Molly M Kressler 

## inference from paths 2 & 3, 27 February 2024 :: Molly M Kressler & Rich B Sherley 
## model4 and model4b, 3 May 2024 :: Molly M Kressler & Rich B Sherley

###################################
########## RUN AT OPEN ############
###################################

## Load Workspace, local macbook

pacman::p_load(MCMCvis,tidyverse,sf,nimble,devtools,flextable,arm,webshot2,sfdep,sp,spdep,beepr, HDInterval, patchwork, cowplot, tidybayes)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd')

## Load workspace, R Server 
pacman::p_load(MCMCvis,tidyverse,sf,nimble,devtools,flextable,webshot2,sfdep,sp,spdep,beepr, HDInterval)
setwd("~/resource/data_and_RDS_NOTforupload")


###################################
##########     END     ############
###################################

#####################################################
###### DF 1, for process model for sharkiness  ######
### Define Sharkiness
## counts of detections per individual per Site (buffer/receiver)

pointdata<-read.csv('pointdata_juvlemons_withAllCoV_MKthesis20192020.csv')%>%
    mutate(tide = case_when(tidephs == 'L' ~ 1, tidephs == 'H' ~ 0))%>%
    mutate(buffIDnum = parse_number(buffID))
stopifnot(nrow(pointdata)==560*2) # check 
sapply(pointdata,class)
summary(pointdata)
  
	# Table of random sample of data frame for graphical methods figure 
  pits <- unique(pointdata$PIT)
  tags <- list(paste0('Tag',seq(1:16)))
  both <- bind_cols(pits,tags) %>%rename(PIT=...1,tagID=...2)
  pointdata <- left_join(pointdata,both) %>%relocate(tagID, .after=PIT)
  head(pointdata)
  pointdata%>%filter(PIT=='molly5') # check that PIT is gettin assigned the same Tag ID the whole way through

	pointflex <- pointdata%>% 
		as_tibble%>%
    dplyr::select(tagID, starts_with('st_'))%>%
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
		save_as_image(pointflex, 'resource_chp3/forgraphicalmethodsfigure_chapter3_datapreparation_POINTDATA_oct2024.png',webshot='webshot')

####################################################
###### DF 2, for process model for fishiness  ######

hexdata<-read.csv('hexdata_juvLemonsPrUse_withAllCov_MKthesis.csv')%>%
    mutate(st_shark = (m5_PUSE - mean(m5_PUSE))/sd(m5_PUSE))%>%
    mutate(tide = case_when(tidephs == 'L' ~ 1, tidephs == 'H' ~ 0))
stopifnot(nrow(hexdata)==5296) # check 
summary(hexdata)

	# Table of random sample of data frame for graphical methods figure 
	hexflex <- hexdata%>%
	dplyr::rename('Hex.No' = 'jcode')%>%
    dplyr::select(Hex.No, starts_with('st_'))%>%
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
		save_as_image(hexflex, 'resource_chp3/forgraphicalmethodsfigure_chapter3_datapreparation_HEXDATA_oct2024.png',webshot='webshot')



##################################################################
##################################################################

##################### 
## - MODEL6: tidal version of model5a. Covariates, predicted or observed, were done so with respect to tide (high/low).

  modelCode6<-nimbleCode({
      

    ##########################
    ######### priors #########
    
    #### prior for sharkiness #### 
    for(i in 1:8){
       a[i] ~ dnorm(0,.01) 
      }
    for(i in 1:2){
       g[i] ~ dnorm(0,.01)# for pred and fish only models, for interpretting coefficients 
      }
    for(i in 1:1){
      h[i] ~ dnorm(0, 0.01) # for tide on sharks - the difference between high and low tide 
    }
      # prior for residual variance - sharks 
      tau.shark ~ dgamma(0.01,0.01) 
      sigma.shark <- sqrt(1/tau.shark)
      tau.shark.hex ~ dgamma(0.01, 0.01) # for glm of fish on sharks, at hexagon level
      sigma.shark.hex <- sqrt(1/tau.shark.hex)
    
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
      # prior for intercept - sharks 
      b ~ dnorm(0,.001)
      tau.epsi_shark ~ dgamma(0.01,0.01)
      sigma.epsi_shark <- sqrt(1/tau.epsi_shark)
    
    #### prior for fishiness (cross metric, hexagons) ####
    for(i in 1:6){
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

    ## glmm for the effect of tide (only), informed by local knowledge of the impact of tides 
      for(i in 1: point.N){
        standard.shark4[i] ~ dnorm(tideonly.mu[i], tau.shark)

        tideonly.mu[i] <- epsi_shark[buffID[i]] + h[1]*pointtide[i] 
      }

     ## informed by hypothesis exploration 

    ### data model for seagrasses - hexagons
      for(i in 1:hex.N){
      standard.hexsg[i] ~ dnorm(hexsg.mu[i],tau.hexsg)

      hexsg.mu[i] <- k + j[1]*standard.hexdist2shore[i] + j[2]*standard.hexdistcmg[i] + j[3]*standard.hexdist2jetty[i] + j[4]*standard.hexdist2shore[i]*standard.hexdistcmg[i] 
    }
    
    ### data model for fishiness - hexagons 
     for(i in 1:hex.N){
      standard.hexfish[i] ~ dnorm(fish.mu[i],tau.fish) 

      fish.mu[i] <- d + c[1]*standard.hexdist2jetty[i] + c[2]*standard.hexdist2shore[i] + c[3]*standard.hexmds[i]+ c[4]*standard.hexhds[i]+ c[5]*standard.hexdist2jetty[i]*standard.hexdist2shore[i] + c[6]*hextide[i]

      # glm for the effect of fish (only)
      standard.hexshark[i] ~ dnorm(fishonly.mu[i], tau.shark.hex)
      fishonly.mu[i] <- d + g[1]*standard.hexfish[i]
    }
  
    ### data model for large shark detectons - point data 
     for(i in 1:point.N){
      zlogit.sqzrisk[i] ~ dnorm(pred.mu[i], tau.pred)

      pred.mu[i] <- f + e[1]*standard.pointdist2shore[i] + e[2]*standard.pointdepth[i] + e[3]*standard.pointdistcmg[i] + e[4]*standard.pointdepth[i]*standard.pointdistcmg[i]

      ## glm for the effect of predation (only) 
      standard.shark3[i] ~ dnorm(predonly.mu[i], tau.shark)
      predonly.mu[i] <-  epsi_shark[buffID[i]] + g[2]*zlogit.sqzrisk.pred[i]

     }

    ### data model for sharkiness - pointdata
     for(i in 1:point.N){
      standard.shark[i] ~ dnorm(shark.mu[i],tau.shark)   
      shark.mu[i] <- epsi_shark[buffID[i]] + a[1]*standard.fish.pred[i] + a[2]*standard.pointdist2shore[i] + a[3]*standard.pointdistcmg[i] + a[4]*standard.sg[i] + a[5]*standard.pointdepth[i]+ a[6]*standard.pointdist2jetty[i]+ a[7]*zlogit.sqzrisk[i] + a[8]*pointtide[i]

      }  

    ######### Derived Parameters #########
    # for estmating total pathways 
    ## coefficients for distance metrics are from process models of those predictors. leave out coefficients for interactions if the fixed effects coefficients are already included.

    path[1] <-  a[4] + h[1] # seagrass and tide

    path[2] <-  a[1] + c[1] + c[2] + c[3] + c[4] + a[8] # fish, the things that effect fish, and tide

    path[3] <-  a[7] + e[1] + e[2] + e[3] + a[8] # predator risk, the things that effect risk, and tide

    path[4] <-  a[8] + a[2] + a[3] + a[4] + a[5] + a[6] # tide - should it have abiotics in it? tide in the context of the habitat? Because similiarly you have 'plain' coefficiednt estimates g1 an g2 for fish and pred. 

    path[5] <-  a[7] + e[1] + e[2] + e[3] # predator risk, the things that effect risk, no tide

    path[6] <-  a[1] + c[1] + c[2] + c[3] + c[4] # fish, the things that effect fish, no tide

    ## need to calculate the values from the abiotics along the paths to the initial parameter, e.g. abtiocs to fishes in path 1. 

    value[1] <-  j[1] * j[2] * j[3]  # path 2, shore * refuge * jetty
    value[2] <- e[1] * e[2] * e[3] * e[4] # path 3, shore * jetty * depth 
    value[3] <- a[2] + a[3] + a[4] + a[5] + a[6]  # abiotics in path 4 for tide

  }) # end of model code 


    ##########################
    ## Compile the model code
    ##########################
    
    myConstants6<-list(point.N=1120,hex.N=5296,B=max(pointdata$buffIDnum),buffID=pointdata$buffIDnum)
    
    myData6<-list(
      # tell nimble the covariates 
      # point data
      standard.shark = pointdata$st_shark,
      standard.pointdist2shore = pointdata$st_d2shore,
      standard.pointdistcmg = pointdata$st_distcmg,
      standard.fish.pred = pointdata$st_maxN,
      standard.sg = pointdata$st_SG,
      standard.pointdepth= pointdata$st_depth,
      standard.pointdist2jetty= pointdata$st_d2jetty,
      zlogit.sqzrisk = pointdata$st_risk,
      zlogit.sqzrisk.pred = pointdata$st_risk,
      standard.shark3 = pointdata$st_shark,
      standard.shark4 = pointdata$st_shark,
      pointtide = pointdata$tide,
      # hex data 
      standard.hexfish = hexdata$st_maxN,
      standard.hexdist2shore = hexdata$st_d2sh,
      standard.hexdistcmg = hexdata$st_dstc,
      standard.hexdist2jetty = hexdata$st_d2jetty,
      standard.hexmds = hexdata$st_mds,
      standard.hexhds = hexdata$st_hds,
      standard.hexshark = hexdata$st_shark,
      hextide = hexdata$tide,
      standard.hexsg = hexdata$st_SG
      )
    
    init.values6<-list(a=rnorm(8,0,1),
                        b=rnorm(1),
                        c=rnorm(6,0,1),
                        d=rnorm(1),
                        e=rnorm(4,0,1),
                        f=rnorm(1),
                        j=rnorm(4,0,1),
                        k=rnorm(1),
                        g=rnorm(2,0,1),
                        h=rnorm(1,0,.1),
                        path=rnorm(6,0,1),
                        value=rnorm(3,0,.05),
                        epsi_shark=rgamma(35,1,0.1),
                        tau.shark=rgamma(1,0.01,0.01),
                        tau.fish=rgamma(1,0.01,0.01),
                        tau.pred=rgamma(1,0.01,0.01),
                        tau.epsi_shark=rgamma(1,0.01,0.01),
                        tau.hexsg=rgamma(1,0.01,0.01),
                        tau.shark.hex = rgamma(1, 0.01, 0.01)
    )
    
    ## model4 define and compile
    model6<-nimbleModel(code=modelCode6, name="model6",data=myData6,constants = myConstants6,inits=init.values6) #define the model
    
    model6$calculate() # if = NA, indicates missing or invalid initial values, and you have to fix the model until it is numeric.
    model6$initializeInfo()
    
    Cm6<-compileNimble(model6) # compile the model

    
    ##########################
    ## Compile & Run the MCMC
    ##########################
    
    ## model5
    conf6<- configureMCMC(model6,monitors=c('a','b','c','d','j','k','e','f','g', 'h','tau.epsi_shark','tau.fish','tau.shark','tau.pred','tau.hexsg','path','value'),onlySlice=FALSE) 
    MCMC_model6 <- buildMCMC(conf6,na.rm=TRUE)
    ccMCMC6 <-compileNimble(MCMC_model6, project = model6)
    samples6 <- runMCMC(ccMCMC6,niter=5000, nburnin=1000, nchains=3,samplesAsCodaMCMC = TRUE) 

    summary(samples6)
    
    saveRDS(samples6,'resource_chp3/nimblemodel_outputs/mcmcsamples_model6_niter5000_burn1000_chains3_oct2024.RDS')

    mcmc_summary_Cmodel6<-MCMCsummary(samplesList6a,round=4,pg0=TRUE,prob=c(0.05,0.95))%>%
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
    mcmc_summary_Cmodel6
    
    # caterpillar plots 
    
    MCMCplot(samples6,ci=c(50,95),params=c('path')) # point = median, thick line = 50% CI, thin line = 95% CI 
    
    # trace and density plots
    
    MCMCtrace(samplesList6,pdf=TRUE,ind=TRUE, Rhat=TRUE, n.eff=TRUE) # ind = TRUE, separate density lines per chain. # pdf = FALSE, don't export to a pdf automatically. 
    
    
    ###########################################################
    ## Summary Table & Caterpillar plots with MCMCvis & tidybayes to show small values ##
    ###########################################################
    
    # import RDS, local macbook 
    samplesList6a <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model6_niter5000_burn1000_chains3_oct2024.RDS')

    mcmc_summary_Cmodel6_samplesListfromRDS<-MCMCsummary(samplesList6a,round=4,probs=c(0.05,0.95),pg0=TRUE)%>%
      tibble::rownames_to_column()%>%
      rename_with(str_to_title)%>%
      rename('pg0'='P>0')%>%
      mutate(pg00 = case_when(Mean >= 0 ~ as.numeric(pg0), Mean < 0 ~ 1-as.numeric(pg0), .default = as.numeric(pg0)))%>%
      rename(Parameter = Rowname, 'Prop. of posterior with \n\ same sign as estimate' = 'pg00', Estimate = 'Mean','lower'='5%',upper='95%')%>%
      mutate(CI = paste0('[',lower,',',upper,']'),.after='Estimate')%>%
      dplyr::select(-lower,-upper,-Sd, -pg0)%>% 
      filter(Parameter!='value[1]', Parameter!='value[2]', Parameter!='g[2]', Parameter!='g[1]')%>%
      flextable()%>%
      theme_zebra()%>%
      set_header_labels(rowname = 'Coefficient',SD='Sd')%>%
      align(align = 'center', part = 'all')%>%
      font(fontname = 'Arial', part = 'all')%>%
      color(color='black',part='all')%>%
      fontsize(size = 10, part = 'all')%>%
      autofit()
    mcmc_summary_Cmodel6_samplesListfromRDS
    
    save_as_image(mcmc_summary_Cmodel6_samplesListfromRDS,path='resource_chp3/nimblemodel_outputs/mcmcsamples__model6_niter5000_burn1000_chains3_oct2024.png',res=850)  
    save_as_docx(mcmc_summary_Cmodel6_samplesListfromRDS,path='resource_chp3/nimblemodel_outputs/mcmcsamples__model6_niter5000_burn1000_chains3_oct2024.docx')  
    
    # grab draws with gather_draws and create label for paths based on iterations and sequence of paths minN to maxN. 
    
    d6 <- gather_draws(samplesList6a,path[])%>%
      group_by(.chain)%>%
      mutate(pathID = paste0('path',rep(1:6, each=4000)))%>% 
      mutate(pathIDnum = rep(1:6, each=4000))%>%
      ungroup() # each path estimate for each chain (of 3) in order, starting with path[1] first estimate in chain 1 
    
    # use tidybayes to plot 
    
    caterpillars <- ggplot(d6, aes(y=reorder(pathID,pathIDnum,decreasing=T),x=.value))+
      stat_pointinterval(.width=c(.50,.8),point_size=2)+
      ylab('Path ID')+
      xlab('Parameter Estimate & CI')+
      geom_vline(xintercept=0,linetype=3)+
      guides(col = 'none')+
      theme_bw()
    caterpillars  
    ggsave(caterpillars,file='resource_chp3/nimblemodel_outputs/caterpillarsPlot__model6_niter5000_burn1000_chains3_oct2024.png',device='png',dpi=400,width=5,height=4,units='in')
    

    





    
##################### 
## - MODEL5A:  Paths informed by hypothesis testng and dredge models - Relative predation risk not pressure. Updated July 2024 with Driscoll BRUVS from 2018 and Bullock BRUVS from 2014 with maxN calculated for the four families. 
    
    modelCode5a<-nimbleCode({
      

    ##########################
    ######### priors #########
    
    #### prior for sharkiness #### 
    for(i in 1:7){
       a[i] ~ dnorm(0,.01) 
      }
    for(i in 1:2){
       g[i] ~ dnorm(0,.01)# for pred and fish only models, for interpretting coefficients 
      }
      # prior for residual variance - sharks 
      tau.shark ~ dgamma(0.01,0.01) 
      sigma.shark <- sqrt(1/tau.shark)
      tau.shark.hex ~ dgamma(0.01, 0.01) # for glm of fish on sharks, at hexagon level
      sigma.shark.hex <- sqrt(1/tau.shark.hex)
    
    #### prior for predators #### 
    for(i in 1:5){
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
      # prior for intercept - sharks 
      b ~ dnorm(0,.001)
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

      hexsg.mu[i] <- k + j[1]*standard.hexdepth[i] + j[2]*standard.hexdistcmg[i] + j[3]*standard.hexdist2jetty[i] + j[4]*standard.hexdist2jetty[i]*standard.hexdistcmg[i] 
    }
    
    ### data model for fishiness - hexagons 
     for(i in 1:hex.N){
      standard.hexfish[i] ~ dnorm(fish.mu[i],tau.fish) 

      fish.mu[i] <- d + c[1]*standard.hexdist2jetty[i] + c[2]*standard.hexdist2shore[i] + c[3]*standard.hexsg[i]+ c[4]*standard.hexdist2jetty[i]*standard.hexdist2shore[i]

      # glm for the effect of fish (only)
      standard.shark2[i] ~ dnorm(fishonly.mu[i], tau.shark.hex)
      fishonly.mu[i] <- d + g[1]*standard.hexfish[i]
    }
  
    ### data model for large shark detectons - point data 
     for(i in 1:point.N){
      zlogit.sqzrisk[i] ~ dnorm(pred.mu[i], tau.pred)

      pred.mu[i] <- f + e[1]*standard.pointdist2shore[i] + e[2]*standard.pointdepth[i]+ e[3]*standard.pointdist2jetty[i] + e[4]*standard.pointdistcmg[i] + e[5]*standard.pointdepth[i]*standard.pointdist2jetty[i]

      ## glm for the effect of predation (only) 
      standard.shark3[i] ~ dnorm(predonly.mu[i], tau.shark)
      predonly.mu[i] <-  epsi_shark[buffID[i]] + g[2]*zlogit.sqzrisk.pred[i]

     }

    ### data model for sharkiness - pointdata
     for(i in 1:point.N){
      standard.shark[i] ~ dnorm(shark.mu[i],tau.shark)   
      shark.mu[i] <- epsi_shark[buffID[i]] + a[1]*standard.fish.pred[i] + a[2]*standard.pointdist2shore[i] + a[3]*standard.pointdistcmg[i] + a[4]*sgPCA1[i] + a[5]*standard.pointdepth[i]+ a[6]*standard.pointdist2jetty[i]+ a[7]*zlogit.sqzrisk[i]

      }  

    ######### Derived Parameters #########
    # for estmating total pathways 
    ## coefficients for distance metrics are from process models of those predictors
    ## paths 1-3 assess the support of the effect of the primary parameter, e.g. pressure, on juvenile spatial baheviour as that primary parameter is determined/influenced by the abiotic habitat features, e.g. depth  

    path[1] <-  j[1] * j[2] * j[3] * a[4]  # seagrass dredge informed path: sg dist2jetty dist2shore distcmg 

    path[2] <-  j[1] * j[3] * c[2] * c[3] * a[1] # fish dredge informed path: maxN sg (depth + jetty) dist2shore dist2jetty (covered with sg)

    path[3] <-  e[1] * e[2] * e[3] * e[4] * a[7] # predator pressure informed path: predators depth dist2shore dist2jetty distcmg

    ## need to calculate the values from the abiotics along the paths to the initial parameter, e.g. abtiocs to fishes in path 1. 

    value[1] <-  j[1] * j[2] * j[3]  # path 2, shore * refuge * jetty
    value[2] <- e[1] * e[2] * e[3] * e[4] # path 3, shore * jetty * depth 

  }) # end of model code 


    ##########################
    ## Compile the model code
    ##########################
    
    myConstants5a<-list(point.N=560,hex.N=2663,B=max(pointdata$buffID),buffID=pointdata$buffID)
    
    myData5a<-list(
      # tell nimble the covariates 
      # point data
      standard.shark = pointdata$standard.shark,
      standard.pointdist2shore = pointdata$standard.dist2shore,
      standard.pointdistcmg = pointdata$standard.distcmg,
      standard.fish.pred = pointdata$standard.fish,
      sgPCA1 = pointdata$sgPCA1,
      standard.pointdepth= pointdata$standard.depth,
      standard.pointdist2jetty= pointdata$standard.dist2jetty,
      zlogit.sqzrisk = pointdata$zlogit.sqzrisk,
      zlogit.sqzrisk.pred = pointdata$zlogit.sqzrisk,
      standard.shark3 = pointdata$standard.shark,
      # hex data 
      standard.hexfish = hexdata$standard.hexfish,
      standard.hexdist2shore = hexdata$standard.hexdist2shore,
      standard.hexdistcmg = hexdata$standard.hexdistcmg,
      standard.hexdist2jetty = hexdata$standard.dist2jetty,
      standard.hexdepth = hexdata$standard.depth,
      standard.hexsg = hexdata$standard.sgPCA1,
      standard.shark2 = hexdata$standard.hexshark
      )
    
    init.values5a<-list(a=rnorm(7,0,1),
                        b=rnorm(1),
                        c=rnorm(4,0,1),
                        d=rnorm(1),
                        e=rnorm(5,0,1),
                        f=rnorm(1),
                        j=rnorm(4,0,1),
                        k=rnorm(1),
                        g=rnorm(2,0,1),
                        path=rnorm(3,0,.05),
                        value=rnorm(2,0,.05),
                        epsi_shark=rgamma(35,0.01,0.01),
                        tau.shark=rgamma(1,0.01,0.01),
                        tau.fish=rgamma(1,0.01,0.01),
                        tau.pred=rgamma(1,0.01,0.01),
                        tau.epsi_shark=rgamma(1,0.01,0.01),
                        tau.hexsg=rgamma(1,0.01,0.01),
                        tau.shark.hex = rgamma(1, 0.01, 0.01)
    )
    
    ## model4 define and compile
    model5a<-nimbleModel(code=modelCode5a, name="model5a",data=myData5a,constants = myConstants5a,inits=init.values5a) #define the model
    
    model5a$calculate() # if = NA, indicates missing or invalid initial values, and you have to fix the model until it is numeric.
    model5a$initializeInfo()
    
    Cm5a<-compileNimble(model5a) # compile the model

    
    ##########################
    ## Compile & Run the MCMC
    ##########################
    
    ## model5
    conf5a<- configureMCMC(model5a,monitors=c('a','b','c','d','j','k','e','f','g','tau.epsi_shark','tau.fish','tau.shark','tau.pred','tau.hexsg','path','value'),onlySlice=FALSE) 
    MCMC_model5a <- buildMCMC(conf5a,na.rm=TRUE)
    ccMCMC5a <-compileNimble(MCMC_model5a, project = model5a)
    samples5a <- runMCMC(ccMCMC5a,niter=5000, nburnin=1000, nchains=3,samplesAsCodaMCMC = TRUE) 

    summary(samples5a)
    
    saveRDS(samples5a,'resource_chp3/nimblemodel_outputs/mcmcsamples_model5a_niter5000_burn1000_chains3_July2024.RDS')

    mcmc_summary_Cmodel5a<-MCMCsummary(samples5a,round=4,pg0=TRUE,prob=c(0.05,0.95))%>%
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
    mcmc_summary_Cmodel5a
    
    # caterpillar plots 
    
    MCMCplot(samples5a,ci=c(50,95),params=c('path')) # point = median, thick line = 50% CI, thin line = 95% CI 
    
    # trace and density plots
    
    MCMCtrace(samplesList5a,pdf=TRUE,ind=TRUE, Rhat=TRUE, n.eff=TRUE) # ind = TRUE, separate density lines per chain. # pdf = FALSE, don't export to a pdf automatically. 
    
    
    ###########################################################
    ## Summary Table & Caterpillar plots with MCMCvis & tidybayes to show small values ##
    ###########################################################
    
    # import RDS, local macbook 
    samplesList5a <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model5a_niter5000_burn1000_chains3_July2024.RDS')

    mcmc_summary_Cmodel5a_samplesListfromRDS<-MCMCsummary(samplesList5a,round=4,probs=c(0.05,0.95),pg0=TRUE)%>%
      tibble::rownames_to_column()%>%
      rename_with(str_to_title)%>%
      rename('pg0'='P>0')%>%
      mutate(pg00 = case_when(Mean >= 0 ~ as.numeric(pg0), Mean < 0 ~ 1-as.numeric(pg0), .default = as.numeric(pg0)))%>%
      rename(Parameter = Rowname, 'Prop. of posterior with \n\ same sign as estimate' = 'pg00', Estimate = 'Mean','lower'='5%',upper='95%')%>%
      mutate(CI = paste0('[',lower,',',upper,']'),.after='Estimate')%>%
      dplyr::select(-lower,-upper,-Sd, -pg0)%>% 
      filter(Parameter!='value[1]', Parameter!='value[2]', Parameter!='g[2]', Parameter!='g[1]')%>%
      flextable()%>%
      theme_zebra()%>%
      set_header_labels(rowname = 'Coefficient',SD='Sd')%>%
      align(align = 'center', part = 'all')%>%
      font(fontname = 'Arial', part = 'all')%>%
      color(color='black',part='all')%>%
      fontsize(size = 10, part = 'all')%>%
      autofit()
    mcmc_summary_Cmodel5a_samplesListfromRDS
    
    save_as_image(mcmc_summary_Cmodel5a_samplesListfromRDS,path='resource_chp3/nimblemodel_outputs/mcmcsamples_model5a_niter20000_burn12000_chains3_july2024.png',res=850)  
    save_as_docx(mcmc_summary_Cmodel5a_samplesListfromRDS,path='resource_chp3/nimblemodel_outputs/mcmcsamples_model5a_niter20000_burn12000_chains3_july2024.docx')  
    
    # grab draws with gather_draws and create label for paths based on iterations and sequence of paths minN to maxN. 
    
    d5a <- gather_draws(samplesList5a,path[])%>%
      group_by(.chain)%>%
      mutate(pathID = paste0('path',rep(1:3, each=4000)))%>% 
      mutate(pathIDnum = rep(1:3, each=4000))%>%
      ungroup() # each path estimate for each chain (of 3) in order, starting with path[1] first estimate in chain 1 
    
    #write.csv(d3b,'resource_chp3/nimblemodel_outputs/samples3b_spreadlong_model3b_niter10000_burn2000_chains3_4dec2023.csv')
    
    j4_draws <- gather_draws(samplesList5a,j[])%>%
      group_by(.chain)%>%
      ungroup()
    j4<- as.data.frame(j4_draws[,5])
    head(j4)
    samplesList4$chain
    
    ## pivot_wider to spread the pathID column into multiple columns, and fill with .value. 
    
    w3b <- pivot_wider(d3b%>%dplyr::select(-pathIDnum), names_from=pathID, values_from=.value)
    nrow(w3b)
    head(w3b)
    names(w3b)
    summary(w3b)
    
    #write.csv(w3b,'resource_chp3/nimblemodel_outputs/samples3a_pivotwider_model3a_niter20000_burn2000_chains3_4dec2023.csv')
    
    # use tidybayes to plot 
    
    caterpillars <- ggplot(d5a, aes(y=reorder(pathID,pathIDnum,decreasing=T),x=.value, col = pathID))+
      stat_pointinterval(.width=c(.50,.95),point_size=2)+
      scale_color_manual(values = c('#3A6C74','#FFBD59','#020F75'))+
      ylab('Path ID')+
      xlab('Estimate (mean) & CI (.5,.95)')+
      geom_vline(xintercept=0,linetype=3)+
      guides(col = 'none')+
      theme_bw()
    caterpillars  
    ggsave(caterpillars,file='resource_chp3/nimblemodel_outputs/caterpillarsPlot_mcmcsamples__model5a_niter20000_burn12000_chains3_july2024.png',device='png',dpi=400,width=5,height=4,units='in')
    


#####################	
## - MODEL4:  Paths informed by hypothesis testng and dredge models - Relative predation risk not pressure. Updated June 2024 with only Driscoll BRUVs
    
    modelCode4<-nimbleCode({
      

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
    for(i in 1:5){
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

      hexfish.mu[i] <- d + c[1]*standard.hexdist2jetty[i] + c[2]*standard.hexdistcmg[i] + c[3]*standard.hexsg[i]+ c[4]*standard.hexdist2jetty[i]*standard.hexdistcmg[i]
    }
  
    ### data model for large shark detectons - point data 
     for(i in 1:point.N){
      zlogit.sqzrisk[i] ~ dnorm(pred.mu[i], tau.pred)

      pred.mu[i] <- f + e[1]*standard.pointdist2shore[i] + e[2]*standard.pointdepth[i]+ e[3]*standard.pointdist2jetty[i] + e[4]*standard.pointdistcmg[i] + e[5]*standard.pointdepth[i]*standard.pointdist2jetty[i]
     }

    ### data model for sharkiness - pointdata
     for(i in 1:point.N){
      standard.shark[i] ~ dnorm(shark.mu[i],tau.shark)   

      shark.mu[i] <- b + epsi_shark[buffID[i]] + a[1]*standard.fish.pred[i] + a[2]*standard.pointdist2shore[i] + a[3]*standard.pointdistcmg[i] + a[4]*sgPCA1[i] + a[5]*standard.pointdepth[i]+ a[6]*standard.pointdist2jetty[i]+ a[7]*zlogit.sqzrisk[i]
      }  

    ######### Derived Parameters #########
    # for estmating total pathways 
    ## coefficients for distance metrics are from process models of those predictors
    ## paths 1-3 assess the support of the effect of the primary parameter, e.g. pressure, on juvenile spatial baheviour as that primary parameter is determined/influenced by the abiotic habitat features, e.g. depth  

    path[1] <-  j[1] * j[2] * j[3] * j[4] * a[4]  # seagrass dredge informed path: sg dist2jetty dist2shore distcmg 

    path[2] <-  j[2] * j[3] * c[3] * a[1] # fish dredge informed path: fish sg distcmg dist2jetty

    path[3] <-  e[1] * e[2] * e[3] * e[4] * e[5] * a[7] # predator pressure informed path: predators depth dist2shore dist2jetty distcmg

    path[4] <- a[3] * a[6] # refugia and anthropocene - could directly impact movement, so ecologically, a path (disturbance)

    ## need to calculate the values from the abiotics along the paths to the initial parameter, e.g. abtiocs to fishes in path 1. 

    value[1] <-  j[1] * j[2] * j[3] * j[4] # path 2, shore * refuge * jetty
    value[2] <- e[1] * e[2] * e[3] * e[4] * e[5] # path 3, shore * jetty * depth 

  }) # end of model code 


    ##########################
    ## Compile the model code
    ##########################
    
    myConstants4<-list(point.N=560,hex.N=2663,B=max(pointdata$buffIDnum),buffID=pointdata$buffID)
    
    myData4<-list(
      # tell nimble the covariates 
      # point data
      standard.shark = pointdata$standard.shark,
      standard.pointdist2shore = pointdata$standard.dist2shore,
      standard.pointdistcmg = pointdata$standard.distcmg,
      standard.fish.pred = pointdata$standard.fish,
      sgPCA1 = pointdata$sgPCA1,
      standard.pointdepth= pointdata$standard.depth,
      standard.pointdist2jetty= pointdata$standard.dist2jetty,
      zlogit.sqzrisk = pointdata$zlogit.sqzrisk,
      standard.larges = pointdata$standard.larges,
      # hex data 
      standard.hexfish = hexdata$standard.hexfish,
      standard.hexdist2shore = hexdata$standard.hexdist2shore,
      standard.hexdistcmg = hexdata$standard.hexdistcmg,
      standard.hexdist2jetty = hexdata$standard.hexdist2jetty,
      standard.hexsg = hexdata$standard.sgPCA,
      standard.hexdist2jetty = hexdata$standard.hexdist2jetty
    )
    
    init.values4<-list(a=rnorm(7,0,1),
                        b=rnorm(1),
                        c=rnorm(4,0,1),
                        d=rnorm(1),
                        e=rnorm(5,0,1),
                        f=rnorm(1),
                        j=rnorm(4,0,1),
                        k=rnorm(1),
                        path=rnorm(4,0,.05),
                        value=rnorm(2,0,.05),
                        epsi_shark=rgamma(35,0.01,0.01),
                        tau.shark=rgamma(1,0.01,0.01),
                       tau.fish=rgamma(1,0.01,0.01),
                       tau.pred=rgamma(1,0.01,0.01),
                       tau.epsi_shark=rgamma(1,0.01,0.01),
                        tau.hexsg=rgamma(1,0.01,0.01)
    )
    
    ## model4 define and compile
    model4<-nimbleModel(code=modelCode4, name="model4",data=myData4,constants = myConstants4,inits=init.values4) #define the model
    
    model4$calculate() # if = NA, indicates missing or invalid initial values, and you have to fix the model until it is numeric.
    model4$initializeInfo()
    
    Cm4<-compileNimble(model4);beep(2) # compile the model

    
    ##########################
    ## Compile & Run the MCMC
    ##########################
    
    ## model4
    conf4 <- configureMCMC(model4,monitors=c('a','b','c','d','j','k','e','f','tau.epsi_shark','tau.fish','tau.shark','tau.pred','tau.hexsg','path','value'),onlySlice=FALSE) 
    MCMC_model4 <- buildMCMC(conf4,na.rm=TRUE)
    ccMCMC4 <-compileNimble(MCMC_model4, project = model4)
    samples4 <- runMCMC(ccMCMC4,niter=20000, nburnin=12000, nchains=3,samplesAsCodaMCMC = TRUE);beep(2) # predation pressure in model4 = z-scored logit transformed lemon squeezed relative risk (proporiton 0 to 1)

    summary(samples4)
    
    saveRDS(samples4,'resource_chp3/nimblemodel_outputs/mcmcsamples_model4_niter20000_burn12000_chains3_June2024.RDS')

    mcmc_summary_Cmodel4<-MCMCsummary(samples4,round=4,pg0=TRUE,prob=c(0.05,0.95))%>%
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
    mcmc_summary_Cmodel4
    
    # caterpillar plots 
    
    MCMCplot(samples4,ci=c(50,95),params=c('path')) # point = median, thick line = 50% CI, thin line = 95% CI 
    
    # trace and density plots
    
    MCMCtrace(samplesList4,pdf=TRUE,ind=TRUE, Rhat=TRUE, n.eff=TRUE) # ind = TRUE, separate density lines per chain. # pdf = FALSE, don't export to a pdf automatically. 
    
    
    ###########################################################
    ## Summary Table & Caterpillar plots with MCMCvis & tidybayes to show small values ##
    ###########################################################
    
    # import RDS, R Server
    samplesList4 <- readRDS('mcmcsamples_model4_niter20000_burn12000_chains3_4may2024.RDS')
    
    # import RDS, local macbook 
    samplesList4 <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model4_niter20000_burn12000_chains3_June2024.RDS')
    summary(samplesList4)

    
    head(samplesList4$chain1)
    summary(samplesList4)
    str(samplesList4)
    
    mcmc_summary_Cmodel4_samplesListfromRDS<-MCMCsummary(samplesList4,round=4,probs=c(0.05,0.95),pg0=TRUE)%>%
      tibble::rownames_to_column()%>%
      rename_with(str_to_title)%>%
      rename('pg0'='P>0')%>%
      mutate(pg00 = case_when(Mean >= 0 ~ as.numeric(pg0), Mean < 0 ~ 1-as.numeric(pg0), .default = as.numeric(pg0)))%>%
      rename(Parameter = Rowname, 'Prop. of posterior with \n\ same sign as estimate' = 'pg00', Estimate = 'Mean','lower'='5%',upper='95%')%>%
      mutate(CI = paste0('[',lower,',',upper,']'),.after='Estimate')%>%
      dplyr::select(-lower,-upper,-Sd, -pg0)%>%	
      filter(Parameter!='value[1]', Parameter!='value[2]')%>%
      flextable()%>%
      theme_zebra()%>%
      set_header_labels(rowname = 'Coefficient',SD='Sd')%>%
      align(align = 'center', part = 'all')%>%
      font(fontname = 'Arial', part = 'all')%>%
      color(color='black',part='all')%>%
      fontsize(size = 10, part = 'all')%>%
      autofit()
    mcmc_summary_Cmodel4_samplesListfromRDS
    
    save_as_image(mcmc_summary_Cmodel4_samplesListfromRDS,path='resource_chp3/nimblemodel_outputs/mcmcsamples_model4_niter20000_burn12000_chains3_june2024.png',res=850)	
    
    # grab draws with gather_draws and create label for paths based on iterations and sequence of paths minN to maxN. 
    
    d3b <- gather_draws(samplesList4,path[])%>%
      group_by(.chain)%>%
      mutate(pathID = paste0('path',rep(1:4, each=8000)))%>% 
      mutate(pathIDnum = rep(1:4, each=8000))%>%
      ungroup() # each path estimate for each chain (of 3) in order, starting with path[1] first estimate in chain 1 
    head(d3b)
    summary(d3b)
    unique(d3b$pathID) # check
    
    write.csv(d3b,'resource_chp3/nimblemodel_outputs/samples3b_spreadlong_model3b_niter10000_burn2000_chains3_4dec2023.csv')
    
    j4_draws <- gather_draws(samplesList4,j[])%>%
      group_by(.chain)%>%
      ungroup()
    j4<- as.data.frame(j4_draws[,5])
    head(j4)
    samplesList4$chain
    
    ## pivot_wider to spread the pathID column into multiple columns, and fill with .value. 
    
    w3b <- pivot_wider(d3b%>%dplyr::select(-pathIDnum), names_from=pathID, values_from=.value)
    nrow(w3b)
    head(w3b)
    names(w3b)
    summary(w3b)
    
    write.csv(w3b,'resource_chp3/nimblemodel_outputs/samples3a_pivotwider_model3a_niter20000_burn2000_chains3_4dec2023.csv')
    
    # use tidybayes to plot 
    
    caterpillars <- ggplot(d3b, aes(y=reorder(pathID,pathIDnum,decreasing=T),x=.value))+
      stat_pointinterval(.width=c(.50,.95),point_size=1.8)+
      ylab('Path ID')+
      xlab('Estimate (mean) & CI (.5,.95)')+
      geom_vline(xintercept=0,linetype=3)+
      theme_bw()
    caterpillars	
    ggsave(caterpillars,file='resource_chp3/nimblemodel_outputs/caterpillarsPlot_mcmcsamples_model4_niter12000_burn4000_chains3_june2024.png',device='png',dpi=400,width=5,height=5,units='in')
    

    
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
				compose(i=2,j=2, as_paragraph('Juvenile sharks ~ Depth + Dist. to Shore +  \n\ Dist. to Jetty + Predator Pressure'))%>%
				compose(i=3,j=2, as_paragraph('Juvenile sharks ~ Dist. to Jetty + Dist. to Shore + \n\ Dist. to Refuge + Seagrasses'))%>%
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
    samplesList4 <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model4_niter20000_burn12000_chains3_4may2024.RDS')
    
    samplesList5a <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model5a_niter5000_burn1000_chains3_July2024.RDS')

    pointdata<-read.csv('standardised_meancentred_data_for_bayes_structural_EQ_modelling_optionC_sharkiness_fishiness_habitat_JULY24.shp')
    head(pointdata)

    hexdata <- read.csv('data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_andBRTs_july24.csv')
    hexsf <- st_as_sf(st_read('data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_andBRTs_july24.shp'),crs='WGS84')%>%
      rename(standard.hexshark = stndrd_hxs,
        standard.hexfish = stndrd_hxf,
        standard.hexdist2shore = stndrd_h2,
        standard.hexdistcmg = stndrd_hxd,
        standard.hexlowdensg = stndrd_hxl,
        standard.hexmeddensg = stndrd_hxm,
        standard.dist2jetty = stndrd_d2,
        standard.depth = stndrd_d,
        standard.sgPCA1 = st_PCA1,
        zlogit.sqzrisk = zlgt_sq,
        relPropPD = rlPrpPD
        )

	
  ## Make test sample data frames - based on the hexagon df

  	hexsamp <- hexdata %>%
  		sample_n(5)%>%
  		as.data.frame() 


  ## From model4 samplesList, make objects with all draws of each coeefficient from path 2 and path 3, e.g. j4. 
  	# path 2: j4, c3, a1
  	j4.ch<-c(samplesList4$chain1[,22],samplesList4$chain2[,22],samplesList4$chain3[,22])
  	c3.ch<-c(samplesList4$chain1[,11],samplesList4$chain2[,11],samplesList4$chain3[,11])
  	a1.ch<-c(samplesList4$chain1[,1],samplesList4$chain2[,1],samplesList4$chain3[,1])

  	# path 3: e1, e2, e3, e4, e5, a7
  	e1.ch <-c(samplesList5a$chain1[,14], samplesList5a$chain2[,14],samplesList5a$chain3[,14] )
  	e2.ch <-c(samplesList5a$chain1[,15], samplesList5a$chain2[,15],samplesList5a$chain3[,15] )
  	e3.ch <-c(samplesList5a$chain1[,16], samplesList5a$chain2[,16],samplesList5a$chain3[,16] )
    e4.ch <-c(samplesList5a$chain1[,17], samplesList5a$chain2[,17],samplesList5a$chain3[,17] )
    # not in july2024 model5a. e5.ch <-c(samplesList5a$chain1[,18], samplesList5a$chain2[,18],samplesList5a$chain3[,18] )
  	a7.ch <-c(samplesList5a$chain1[,7], samplesList5a$chain2[,7],samplesList5a$chain3[,7] )

  ## Run loops for each path, to calculate path estimate given a hexagon cell

    	## Variable-Coefficient key
    		# j4 = dist2shore * distcmg
    		# c3 = seagrasses
    		# a1 = fish
    		# e1 = dist2shore
    		# e2 = depth
    		# e3 = dist2jetty
        # e4 = distcmg
    		# e5 = dist2jetty*depth
    		# a7 = predator pressure

    	## Define size of objects/matrices
    	J = 3 * 4000 # chains x iterations, 8000 in June, 4000 in July
    	n.hex = nrow(hexdata)	

    	## set up prediction df

    		preds.path2 = as.data.frame(matrix(NA,ncol = J, nrow = n.hex))
    		head(preds.path2)

    		preds.path3 = as.data.frame(matrix(NA,ncol = J, nrow = n.hex))
    		head(preds.path3)


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
    		setTxtProgressBar(pb, i)
    		for(j in 1:J){
    			preds.path3[i,j] <-  (e1.ch[j]*hexdata$standard.hexdist2shore[i]) + (e2.ch[j]*hexdata$standard.depth[i]) + (e3.ch[j]*hexdata$standard.dist2jetty[i]) + (e4.ch[j]*hexdata$standard.hexdistcmg[i]) + (a7.ch[j]*hexdata$zlogit.sqzrisk[i]) 
    		}
    	};beep(3)
          ### not in july2024:  (e5.ch[j]*hexdata$standard.depth[i]*hexdata$standard.hexdistcmg[i])

    	## save calculations
        saveRDS(preds.path3,'resource_chp3/path_inference/path3_estimates_at_hexagons_model5aJuly2024_calcJuly2024.RData')
        saveRDS(preds.path2,'resource_chp3/path_inference/path2_estimates_at_hexagons_model4may2024_calcMay2024.RData')

  ## Calculate marginal means, and HDI (highest density intervals)

  	## deprecated and inaccurate way to calculate HDIs - HDInterval::hdi needs values in columns. 

  		mean(as.numeric(preds.path3[1,]))

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
  			p3pred[i,2] <- mean(as.numeric(preds.path3[i,]))
  			p3pred[i,3] <- hdi(preds.path3[i,])[2]
  			p3pred[i,4] <- hdi(preds.path3[i,])[1]
  		};beep(3)

  		head(p2pred)
  		nrow(p2pred)
  		head(p3pred)

  	## Updated approach: 5 March 2024
  	p2 <- readRDS('resource_chp3/path_inference/path2_estimates_at_hexagons_model4may2024_calcMay2024.RData')
  	p3 <- readRDS('resource_chp3/path_inference/path3_estimates_at_hexagons_model5aJuly2024_calcJuly2024.RData')
  	
    hexdata <- read.csv('data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_andBRTs_july24.csv')%>%mutate(jcode=as.character(jcode))
      head(hexdata)

  	## re-format data
  	dat.p2 <- as.matrix(p2)
  	dat.p2<-dat.p2[,-1]
  	dat.p2[1:10,1:10]
  	dim(dat.p2)[2]
  			# rowMeans(dat)[1]
  			# mean(dat[1,])
  			# hdi(dat[1,])[1]

  	out.p2 <- data.frame(
  		  Index = rep(0,dim(dat.p2)[1]),
  		  Mean = rep(0,dim(dat.p2)[1]),
  		  Low = rep(0,dim(dat.p2)[1]),
  		  Upp = rep(0,dim(dat.p2)[1]))

  	dat.p3 <- as.matrix(p3)
  	out.p3 <- data.frame(
  		  Index = rep(0,dim(dat.p3)[1]),
  		  Mean = rep(0,dim(dat.p3)[1]),
  		  Low = rep(0,dim(dat.p3)[1]),
  		  Upp = rep(0,dim(dat.p3)[1]))


  	## calculate means and HDIs for each hexagon
  	for(i in 1:dim(dat.p2)[1]){
  		  out.p2$Index[i] = i
  		  out.p2$Mean[i] = mean(dat.p2[i,])
  		  out.p2$Low[i] = HDInterval::hdi(dat.p2[i,])[1]
  		  out.p2$Upp[i] = HDInterval::hdi(dat.p2[i,])[2]
  		  } # takes a beat
    	head(out.p2)
  	
  	for(i in 1:dim(dat.p3)[1]){
  		  out.p3$Index[i] = i
  		  out.p3$Mean[i] = mean(dat.p3[i,])
  		  out.p3$Low[i] = HDInterval::hdi(dat.p3[i,])[1]
  		  out.p3$Upp[i] = HDInterval::hdi(dat.p3[i,])[2]
  		  } # takes a beat
    	head(out.p3)

    	## attach jcodes back 
    	out.p2 <- cbind(hexdata%>%dplyr::select('jcode'), out.p2)%>%
    		dplyr::select(-Index)
    	head(out.p2)
    	
    	out.p3 <- cbind(hexdata%>%dplyr::select('jcode'), out.p3)%>%
    		dplyr::select(-Index)
    	head(out.p3)

  ## Save path estimates and path means + HDI dfs

  	saveRDS(out.p2, 'resource_chp3/path_inference/path2_means_andHDI_at_heaxgons_model4_may2024.RData')
  	saveRDS(out.p3, 'resource_chp3/path_inference/path3_means_andHDI_at_hexagons_model5a_july2024.RData')

##################################################
## Path inference diagnostics ##
##################################################

	## local R, macbook 
	p2pred <- readRDS('resource_chp3/path_inference/path2_means_andHDI_at_heaxgons_model4_may2024.RData')%>%
		mutate(jcode = as.character(jcode))
	p3pred <- readRDS('resource_chp3/path_inference/path3_means_andHDI_at_hexagons_model5a_july2024.RData')%>%
		mutate(jcode = as.character(jcode))

## histograms of HDIs
	p2up <- ggplot(data=p2pred, aes(x=Upp))+geom_histogram(binwidth=.2, fill='#143B43', col='#143B43',lwd=0.05) + 
    theme_bw()  + 
    labs(subtitle = 'A') + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())+xlim(-4,4)

	p2low <- ggplot(data=p2pred, aes(x=Low))+geom_histogram(binwidth=.2, fill='#143B43', col='#143B43',lwd=0.05) + 
    theme_bw() + 
    labs(subtitle = 'B') + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())+xlim(-4,4)
	
  p3up <- ggplot(data=p3pred, aes(x=Upp))+geom_histogram(binwidth=.2, fill='#020F75', col='#020F75',lwd=0.05) + 
    theme_bw() + 
    labs(subtitle = 'Upper') + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),plot.subtitle = element_text(face = "italic"))+
    xlim(-4,4)
	p3low <- ggplot(data=p3pred, aes(x=Low))+
    geom_histogram(binwidth=.2, fill='#020F75', col='#020F75',lwd=0.05) + 
    theme_bw() + 
    labs(subtitle = 'Lower') + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),plot.subtitle = element_text(face = "italic"))+
    xlim(-4,4)

  hist.HDIs.path2.path3 <- p2up | p2low | p3up | p3low 
  hist.HDIs.path3 <-  p3up | p3low 

  ggsave(hist.HDIs.path2.path3, file = 'resource_chp3/path_inference/histograms_HDIs_path2and3.png', dpi = 850, units = 'in', height = 3, width = 8.5)
  ggsave(hist.HDIs.path3, file = 'resource_chp3/path_inference/histograms_HDIs_path3_model5_july24.png', dpi = 850, units = 'in', height = 3, width = 8.5)
 

##################################################
## Path predictions, spatial ##
##################################################
	pacman::p_load(sf, ggplot2, patchwork,tidyverse)
  
    hexsf <- st_as_sf(st_read('data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_andBRTs_july24.shp'),crs='WGS84')%>%
      rename(standard.hexshark = stndrd_hxs,
        standard.hexfish = stndrd_hxf,
        standard.hexdist2shore = stndrd_h2,
        standard.hexdistcmg = stndrd_hxd,
        standard.hexlowdensg = stndrd_hxl,
        standard.hexmeddensg = stndrd_hxm,
        standard.dist2jetty = stndrd_d2,
        standard.depth = stndrd_d,
        standard.sgPCA1 = st_PCA1,
        zlogit.sqzrisk = zlgt_sq,
        relPropPD = rlPrpPD
        )

	land <- st_as_sf(st_read('bim_onlyland_noDots.kml'), crs = 'WGS84')
	p2pred <- readRDS('resource_chp3/path_inference/path2_means_andHDI_at_heaxgons_model4_may2024.RData')%>%
		mutate(jcode = as.character(jcode))
	p3pred <- readRDS('resource_chp3/path_inference/path3_means_andHDI_at_hexagons_model5a_july2024.RData')

	head(p3pred)

	## cut out the splattering of stuff to the south for plotting
		 cropbox <- st_as_sf(st_bbox(c(xmin = 79.24, xmax=79.31, ymin = 25.68, ymax = 25.78), crs='WGS84'))
		 hexsf2 <- st_crop(hexsf, c(xmin = -79.23, xmax=-79.31, ymin = 25.68, ymax = 25.78))
		 southerninset <- st_crop(hexsf, c(xmin = -79.30, xmax=-79.31, ymin = 25.66, ymax = 25.68))
	
	## join spatial geomtry by jcode to preds
		p2.sf <- st_as_sf(left_join(hexsf, p2pred, by='jcode'))%>%
		  dplyr::select(jcode, Mean, Low, Upp, geometry)
		
    p2.sf.cropped <- st_as_sf(left_join(hexsf2, p2pred, by='jcode'))%>%
		  dplyr::select(jcode, Mean, Low, Upp, geometry)
		
		p2.southerninset <- st_as_sf(left_join(southerninset, p2pred, by='jcode'))%>%
		  dplyr::select(jcode, Mean, Low, Upp, geometry)
		
		p3.sf <- st_as_sf(left_join(hexsf, p3pred, by='jcode'))%>%
		  dplyr::select(jcode, Mean, Low, Upp, geometry)

		p3.sf.cropped <- st_as_sf(left_join(hexsf2, p3pred, by='jcode'))%>%
		  dplyr::select(jcode, Mean, Low, Upp, geometry)

		p3.southerninset <- st_as_sf(left_join(southerninset,p3pred, by='jcode'))%>%
		  dplyr::select(jcode, Mean, Low, Upp, geometry)

	## path 2 plots - mean, upper, lower 

  colours = c('#ffffff', '#666282','#3a6c74','#143B43')
  colours.diff = c('#d2c9cf','#a097a6','#858093','#696b81','#4c566e','#2d435a','#063146') 
	
	 #scale_fill_steps(name = 'Mean', low = colours[1], high = colours[3],right = TRUE,space = 'Lab', guide = 'coloursteps', aesthetics = 'fill', breaks = c(-3,0,.2,.6),show.limits = TRUE, labels=scales::label_number(accuracy=0.1))+


	p2.mean <- ggplot()+
		theme_bw()+
		geom_sf(data = p2.sf, aes(fill = Mean), lwd=0, col = 'white')+
		geom_sf(data = land, col = 'grey75', lwd=0.5)+
		scale_fill_gradientn(colors = c('#C8D9DA', '#073B46'), limits = c(-1,1), oob = scales::squish, name = 'Mean')+
	  	theme(axis.text.x = element_text(angle=45, hjust = 1))
	 p2.mean

	p2.lower <- ggplot()+
		geom_sf(data = p2.sf, aes(fill = Low), lwd=0)+
		geom_sf(data = land, col = 'grey75', lwd=0.5)+
	  theme_bw()+
		scale_fill_gradientn(colors = c('#C8D9DA', '#073B46'), limits = c(-1,1), oob = scales::squish, name = 'Lower HDI')+
	  	theme(axis.text.x = element_text(angle=45, hjust = 1),legend.direction = "horizontal", legend.position = "bottom")
	p2.lower

	p2.upper <- ggplot()+
		geom_sf(data = p2.sf, aes(fill = Upp), lwd=0)+
		geom_sf(data = land, col = 'grey75', lwd=0.5)+
		scale_fill_gradientn(colors = c('#C8D9DA', '#073B46'), limits = c(-1,1), oob = scales::squish, name = 'Upper HDI')+
	  theme_bw()+
	  	theme(axis.text.x = element_text(angle=45, hjust = 1),  legend.direction = "horizontal", legend.position = "bottom")
  p2.upper
  
    ## what if...plot the lower/uppder HDIs as the difference between mean and the HDI
      p2.upper.diff <- ggplot()+
        geom_sf(data = p2.sf, aes(fill = abs(Upp-Mean)), lwd=0)+
		geom_sf(data = land, col = 'grey75', lwd=0.5)+
		scale_fill_gradientn(colors = c('#C8D9DA', '#073B46'), limits = c(0,0.75), oob = scales::squish, name = 'Upper HDI\n\ Absolute\n\ difference')+
        theme_bw()+
        theme(axis.text.x = element_text(angle=45, hjust = 1))
      p2.upper.diff
      
      p2.lower.diff <- ggplot()+
        geom_sf(data = p2.sf, aes(fill = abs(Mean-Low)), lwd=0)+
		geom_sf(data = land, col = 'grey75', lwd=0.5)+
		scale_fill_gradientn(colors = c('#C8D9DA', '#073B46'), limits = c(0,0.75), oob = scales::squish, name = 'Lower HDI\n\ Absolute\n\ difference')+
        theme_bw()+
        theme(axis.text.x = element_text(angle=45, hjust = 1))
      p2.lower.diff
      

      
	## path 3 plots - mean, upper, lower 

  p3.mean <- ggplot()+
    geom_sf(data = p3.sf, aes(fill = abs(Mean)), lwd=0)+
    geom_sf(data = land, col = 'grey75', lwd=0.5)+
    theme_bw()+
    scale_fill_gradientn(colors = c('#fafcfc', '#073B46'), limits = c(0,5), oob = scales::squish, name = 'Absolute \n\ Effect')+
      theme(axis.text.x = element_text(angle=45, hjust = 1))
   ggsave(p3.mean, file = 'resource_chp3/path_inference/path3_absolute_magnitude_of_effect_estimates_spatial_model5_july24.png', device = 'png', unit = 'in', dpi = 900, width = 5)
 

  p3.mean <- ggplot()+
    geom_sf(data = p3.sf, aes(fill = (Mean)), lwd=0)+
    geom_sf(data = land, col = 'grey75', lwd=0.5)+
    theme_bw()+
    scale_fill_gradientn(colors = c('#fafcfc', '#073B46'), limits = c(0,5), oob = scales::squish, name = 'Mean')+
      theme(axis.text.x = element_text(angle=45, hjust = 1))
   ggsave(p3.mean, file = 'resource_chp3/path_inference/path3_mean_estimates_spatial_model5_july24.png', device = 'png', unit = 'in', dpi = 900, width = 5)
 

	p3.lower <- ggplot()+
		geom_sf(data = p3.sf, aes(fill = abs(Low)), lwd=0)+
		geom_sf(data = land, col = 'grey75', lwd=0.5)+
	  theme_bw()+
    scale_fill_gradientn(colors = c('#fafcfc', '#073B46'), limits = c(0,5), oob = scales::squish, name = 'Absolute \n\ Effect')+
	  theme(axis.text.x = element_text(angle=45, hjust = 1),  legend.direction = "horizontal", legend.position = "bottom")+
    labs(subtitle = 'Lower')
	 p3.lower

	p3.upper <- ggplot()+
		geom_sf(data = p3.sf, aes(fill = abs(Upp)), lwd=0)+
		geom_sf(data = land, col = 'grey75', lwd=0.5)+
	  theme_bw()+
    scale_fill_gradientn(colors = c('#fafcfc', '#073B46'), limits = c(0,5), oob = scales::squish, name = 'Absolute \n\ Effect')+
	  theme(axis.text.x = element_text(angle=45, hjust = 1), legend.direction = "horizontal", legend.position = "bottom")+
    labs(subtitle = 'Upper')
    p3.upper

	  ## what if...plot the lower/uppder HDIs as the difference between mean and the HDI
        p3.upper.diff <- ggplot()+
          geom_sf(data = p3.sf, aes(fill = abs(Upp-Mean)), lwd=0)+
   		geom_sf(data = land, col = 'grey75', lwd=0.5)+
  		scale_fill_gradientn(colors = c('#C8D9DA', '#073B46'), limits = c(0,1), oob = scales::squish, name = 'Lower')+
          theme_bw()+
          theme(axis.text.x = element_text(angle=45, hjust = 1))
        p3.upper.diff
        
        p3.lower.diff <- ggplot()+
          geom_sf(data = p3.sf, aes(fill = abs(Mean-Low)), lwd=0)+
   		geom_sf(data = land, col = 'grey75', lwd=0.5)+
  		scale_fill_gradientn(colors = c('#C8D9DA', '#073B46'), limits = c(0,1), oob = scales::squish, name = 'Lower')+
          theme_bw()+
          theme(axis.text.x = element_text(angle=45, hjust = 1),)
        p3.lower.diff


	## various arrangements for pub
		## just the means, side by side
		p3mean.noyaxis <- p3.mean + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())
    
    means.. <- p2.mean + p3mean.noyaxis + plot_layout(guides = 'collect')
    means..

    p3upper.noyaxis <- p3.upper + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())
    p2upper.noyaxis <- p2.upper + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())

    p2.hdis.wide.axescollected <- p2.lower + p2upper.noyaxis 
    p3.hdis.wide.axescollected <- (p3.lower | p3upper.noyaxis) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')


		## means with their HDIs, horizontal
    p3lower.noyaxis <- p3.lower + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())
    p2lower.noyaxis <- p2.lower + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())

    p2.mean.hdis.wide <- p2.mean + p3lower.noyaxis + p2upper.noyaxis 
		p3.mean.hdis.wide <- p3.mean + p2lower.noyaxis + p3upper.noyaxis 
		
		## means with their HDIs, stacked
		p2.mean.hdis.long <- p2.mean / p2.lower / p2.upper 
		p3.mean.hdis.long <- p3.mean / p3.lower / p3.upper 

		## means with their HDIs, means = 2 cols and big, hdis = 1 col small under means
		p2.mean.diffshdis.square <- p2.mean / (p2.lower.diff | p2.upper.diff) + plot_layout(widths = c(1), heights = c(2,1))
		p2.mean.hdis.square <- p2.mean / (p2.lower | p2.upper) + plot_layout(widths = c(1), heights = c(2,1))

		p3.mean.hdis.square <- p3.mean / (p3.lower | p3.upper) + plot_layout(widths = c(1), heights = c(2,1))
		p3.mean.diffshdis.square <- p3.mean / (p3.lower.diff | p3.upper.diff) + plot_layout(widths = c(1), heights = c(2,1))

		## mean with difference
		 diffmean.p2 <- p2.lower.diff | p2.upper.diff
		 diffmean.p3 <- p3.lower.diff | p3.upper.diff

		

	  ## save plots 
	  ##### path 2 plots
		ggsave(p2.mean, file = 'resource_chp3/path_inference/path2_mean_estimates_spatial_may24.png', device = 'png', unit = 'in', dpi = 900, height = 6, width = 4.5)
		ggsave(p2.lower, file = 'resource_chp3/path_inference/path2_lowerHDI_estimates_spatial_may24.png', device = 'png', unit = 'in', dpi = 900, height = 6, width = 4.5)
		ggsave(p2.upper, file = 'resource_chp3/path_inference/path2_upperHDI_estimates_spatial_may24.png', device = 'png', unit = 'in', dpi = 900, height = 6, width = 4.5)
		ggsave(p2.mean.hdis.wide, file = 'resource_chp3/path_inference/path2_mean_andHDI_WIDEpanel_estimates_spatial_may24.png', device = 'png', unit = 'in', dpi = 900, height = 6, width = 13)
		ggsave(p2.mean.hdis.long, file = 'resource_chp3/path_inference/path2_mean_andHDI_LONGpanel_estimates_spatial_may24.png', device = 'png', unit = 'in', dpi = 900, height = 12, width = 4.5)
		ggsave(p2.mean.diffshdis.square, file = 'resource_chp3/path_inference/path2_mean_and_diffHDI_SQUAREpanel_estimates_spatial_may24.png', device = 'png', unit = 'in', dpi = 900, height = 8, width = 6)
		ggsave(p2.mean.hdis.square, file = 'resource_chp3/path_inference/path2_mean_and_UppandLowHDI_SQUAREpanel_estimates_spatial_may24.png', device = 'png', unit = 'in', dpi = 900, height = 8, width = 6)
		ggsave(p2.hdis.wide.axescollected, file = 'resource_chp3/path_inference/path2_HDI_WIDEpanel_estimates_spatial_may24.png', device = 'png', unit = 'in', dpi = 900, height = 5, width = 9)

	  ####### path 3 plots 
 
 		ggsave(p3.mean, file = 'resource_chp3/path_inference/path3_mean_estimates_spatial_july24.png', device = 'png', unit = 'in', dpi = 900, height = 7)
		ggsave(p3.lower, file = 'resource_chp3/path_inference/path3_lowerHDI_estimates_spatial_may24.png', device = 'png', unit = 'in', dpi = 900, height = 6, width = 4.5)
		ggsave(p3.upper, file = 'resource_chp3/path_inference/path3_upperHDI_estimates_spatial_may24.png', device = 'png', unit = 'in', dpi = 900, height = 6, width = 4.5)
		ggsave(p3.mean.hdis.wide, file = 'resource_chp3/path_inference/path3_mean_andHDI_WIDEpanel_estimates_spatial_may24.png', device = 'png', unit = 'in', dpi = 900, height = 6, width = 13)
		ggsave(p3.mean.hdis.long, file = 'resource_chp3/path_inference/path3_mean_andHDI_LONGpanel_estimates_spatial_may24.png', device = 'png', unit = 'in', dpi = 900, height = 12, width = 4.5)
		ggsave(p3.mean.hdis.square, file = 'resource_chp3/path_inference/path3_mean_andHDI_SQUAREpanel_estimates_spatial_may24.png', device = 'png', unit = 'in', dpi = 900, height = 8, width = 6)
		ggsave(p3.mean.diffshdis.square, file = 'resource_chp3/path_inference/path3_mean_and_diffsHDI_SQUAREpanel_estimates_spatial_may24.png', device = 'png', unit = 'in', dpi = 900, height = 8, width = 6)
		ggsave(diffmean.p3, file = 'resource_chp3/path_inference/path3_HDIdifference_estimates_spatial_may24.png', device = 'png', unit = 'in', dpi = 900, height = 4, width = 8)
 		ggsave(p3.hdis.wide.axescollected, file = 'resource_chp3/path_inference/path3_HDI_WIDEpanel_Absolute_effect_estimates_spatial_july24.png', device = 'png', unit = 'in', dpi = 900, width = 8)

		
 		ggsave(means.., file = 'resource_chp3/path_inference/path2_and_path3_mean_estimates_spatial_forpubs_collectedaxes_may24.png', device = 'png', unit = 'in', dpi = 900, height = 5, width = 7)


	


