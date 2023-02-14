getwd()
setwd('/Users/mollykressler/Documents/data_phd')

################### BRUVS data cleaner masterfile ################


# compare SD and HGRB data for duplicate entries - use GPS and season for potental copies
	
	# import libraries 
	library('tidyverse')
	library('dplyr')
	library('purrr')
	
	#import files
	HGbig<-read.csv('HGRBdatabig&env.csv', header=T, sep=',')
		HGbig[is.na(HGbig)]<-0
	HGsmall<-read.csv('HGRBdatasmall.csv', header=T, sep=',')
	SDorig<-read.csv('driscolldata_sum&env.csv', header=T, sep=',')
		head(HGbig)
	
	#factor-ise appropriate variables using 'lapply'
	head(HGbig)
		HGbigfactors<-c(1,4:7,11:13,17:18) #define object of columns in HGbig that you want to make a factor
		HGbig[HGbigfactors]<-lapply(HGbig[HGbigfactors],factor) #use lapply to change the class of the columns in the col-object into factors
		sapply(HGbig,class) #check that it worked, it did
	head(HGsmall)
		HGsmall$BRUV_ID<-as.factor(HGsmall$BRUV_ID)
		sapply(HGsmall,class)
	head(SDorig)
		SDorigfactors<-c(1,11:12)
		SDorig[SDorigfactors]<-lapply(SDorig[SDorigfactors],factor)
		sapply(SDorig, class)

	#rename SDorig headers to match HGRB data - BRUV_ID, 
		#site name, BRUV_ID
		SDorig <- SDorig %>% rename(BRUV_ID=SiteName)
			head(SDorig) # all good 
	
	# FOUND ZERO DUPLICATES. trying to identify duplicates between SD and HG datasets. used this code to identify zero duplicates: 
		# I think I can do this with dplyr, using a filtering join -  'anti-join': return rows of 'x' that do not have a match in 'y'. Use to see what is unique in x. So here it iwll tell me which od SD is not already in HG latlongs.
		#rename some columns to make SD and HG match up
	#		SD2<-data.frame(SDorig$Latitude, SDorig$Longitude)
	#			SD2 <- SD2 %>% rename(Latitude=SDorig.Latitude, Longitude=SDorig.Longitude)
	#		HG2<-data.frame(HGsmall$Latitude, HGsmall$Longitude)
	#			HG2 <- HG2 %>% rename(Latitude=HGsmall.Latitude, Longitude=HGsmall.Longitude)
				#just do it with the Lats and longs
	#		antiSDHG<-anti_join(SD2, HG2,by=NULL,copy=FALSE)
	#		summary(antiSDHG)
	#		summary(SD2)
			#they are exactly the same, which would suggest the SD and HG data are separate. cross-checked in excel, no duplicate LAT-LONG combos found (2 duplicate Lats found in SD, but thats fine.)




	#kepp tidal phase state information in HG, and when I bind HG and SD, put NAs in for SD. I will post hoc reinsert tidal phase information. 


#want to join HG and SD - combine cases to combine only summative stats by binding rows - first need to add SD/HG to BRUV_ID to differentiate between HG1 and SD1. 
	##need to cut them down to only what they both have
	## make a data file from HGbig that has environmentals from locations - can use this as an abiotic layer separate from the maxN layers 

	summary(SDorig) 
	head(HGbig) 

	# add SD/HG to identify source
	SDorig$source<- 'SD'
	HGbig$source<- 'HG'

	# create dataframes with  the summative columns: bruv_id, lat, long, spp_abundance (hg MaxN_Species needs rename), spp_rich, shannon (SD shannon needs rename to SW_species), source & tidal state and phase in HG
	SDorig<- SDorig%>%rename(SW_species=Shannon)
	SDorig<- SDorig%>%rename(SW_Species=SW_species)
	SDorig<-SDorig%>%rename(spp_richness=spp_rich)

	HGbig<-HGbig%>%rename(spp_abundance=MaxN_Species)

	SDsumm<-SDorig%>%select(BRUV_ID, Latitude, Longitude, spp_abundance, spp_richness,SW_Species, Season,source)
	HGsumm<-HGbig%>%select(BRUV_ID, Latitude, Longitude, spp_abundance, spp_richness,SW_Species, Season,source,ST_Time,Time_Midpoint,Tidal_State,Tidal_Phase,Tidal_State_Phase,Av_Temperature, Av_Depth,Bottom_type,Av_Visibility,SG_Coverage)
	head(HGsumm)
	SDHG<-bind_rows(HGsumm,SDsumm,row.names=NULL)
	head(SDHG)
	summary(SDHG) #It worked and there are NAs where SD data is missing environmental data

# export cleaned and consolidated bruvs data from SD and HG-RB sources. 
	# data from SDHGsumm
	write.csv(SDHG, '/Users/mollykressler/Documents/data_phd/bruvs_SDHG_summary&environs.csv')

plot(SDHG$SG_Coverage,SDHG$spp_richness,pch=19)
m1<-lm(SDHG$spp_richness~SDHG$SG_Coverage)
	abline(m1, col='blue')








##############################################################
# separate abitics from HG, new file, by location
library('lubridate')
head(HGbig)
	#New variable: average flow
	average flow= (end flow - start flow)/(end time - start time) = flow/minute
#####so after finishing this section I realised ST_flow_D is NOT velocity but degree direction (caridnal), so this was unneeded, but a very good exercise so it's here for notes. I've commented-out the buts that aren't useful for keeping in HGbig. 
	# good website: https://data.library.virginia.edu/working-with-dates-and-time-in-r-using-the-lubridate-package/ 
		timeST<-strptime(HGbig$ST_Time, format="%H:%M")
		timeE<-strptime(HGbig$E_Time, format="%H:%M")
		end<-HGbig$E_Time
		int<-interval(timeST,timeE)
		durationHG<-as.duration(int) ##object with diraton between start and end times
	str(durationHG)
	## add duration values to table
		HGbig$duration<-durationHG
			head(HGbig)
	## difference in flow between start and stop 
		#HGbig$flow_D_diff<-abs((HGbig$ST_Flow_D)-HGbig$E_Flow_D) #absolute value of the difference in flow degree between the start and end
	## divide location flow change (flow_diff) by duration
		#HGbig$avg_flow_D_change<-(HGbig$flow_diff)/(as.numeric(HGbig$duration/dminutes(1)))

##include: Lat/Long, season, Time_midpoint, Tidal-state, tidal_phase, tidal_state_phase, Av_Flow_Velocity, Av-Depth, Sg-coverage, av-saliinity, av-visibility, av-do2, time-midpoint (can assign night/day based on this)

	HGabiotic<-HGbig %>% select(Latitude, Longitude, BRUV_ID, Season, Time_Midpoint, Tidal_State, Tidal_Phase, Tidal_State_Phase, Degree_Diff, Bottom_type, Av_Temperature, Av_Flow_Velocity, Av_Depth, SG_Coverage, Av_Salinity,Av_Visibility, Av_DO2)

	summary(HGabiotic)
	##rename columns, then export 
	HGabiotic <- HGabiotic %>% rename(time_midpoint=Time_Midpoint, tide_state=Tidal_State,tide_phase=Tidal_Phase, tide_state_phase=Tidal_State_Phase, flowD_diff=Degree_Diff, bottom_type=Bottom_type, avg_temp=Av_Temperature, avg_flowV=Av_Flow_Velocity, avg_depth=Av_Depth, sg_cover=SG_Coverage, avg_salt=Av_Salinity, avg_vis=Av_Visibility, avg_do2=Av_DO2)
		head(HGabiotic)
	write.csv(HGabiotic, '/Users/mollykressler/Documents/data_phd/abioticBYlocation_HG-RB.csv')
























