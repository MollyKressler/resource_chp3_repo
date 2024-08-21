## Updated approach for structuring data for path analysis in the package NIMBLE based on co-author feedback in August 2024 
# corresponds to the analysis done in chapter 3 of my thesis.

## created by Molly Kressler, 20 August 2024

## Major changes include: deparating out low and high tide

## Structure: 
	# Add habitat to large receiver array (n=54) 
	# Collate predator data set, with habitat and tides (receivers) - done 21/08
	# Collate juveniles pointdata data frame - one row per receiver per individual with counts, with tide - done 21/08
	# Collate BRUVs with temporally matched habitat data and tides (hexagons) - 2014-L, 2014-H, 2018-L, 2018-H. Then import BRT results (having used the first this dataset) and cbind MaxN predictions. 
	
	## Collate receiver specfic values of predation pressure and juvenile shark detections at low and high tide - done 21/08
	
	# put together hexagon modelling dataframes at low and high tide, standardise and mean centre variables, z.logit transform predation. 
	# put together receiver/pointdata modelling dataframes at low and high tide, standardise and mean centre variables, z.logit transform predation. 

	
	#### YOU LEFT OF ON 21/08, having processed predator and uvenile data into full dataframes. what's missing from those are BRT prediction of prey MaxN. You have not done any of the  BRUVs data handling because you are waiting to get data on tide from Sarah Driscoll. You can find this information in her thesis if you want to extract it. but give her 1-2 weeks and then crcle back to this. 

pacman::p_load(tidyverse,fuzzyjoin,lubridate,sf,ggplot2,cowplot,patchwork,flextable,sf,ggsn,dismo,gbm)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd')

#####
## - Data sets 
#####

## passive acoustic telemetry

p.df <- read_csv('lemonspredators_20192020blacktipsANDbulls/predators_detectedwithinstudyarea_MKthesis2019to2020_aug2024.csv')%>%
	mutate(time = as_datetime(time), tz = 'EST')
p.sf <- st_as_sf(st_read('lemonspredators_20192020blacktipsANDbulls/predators_detectedwithinstudyarea_MKthesis2019to2020_aug2024.shp'))%>%
		dplyr::select(-time)%>%
		mutate(time = as_datetime(p.df$time),.before = 'statn_x')## this doesn't have full date times. 

j <- st_as_sf(st_read('cleaned_detections/detections20192020_NOghosts_withEmilyCourmierhabitat_data4Winter2020_sept22.shp'), crs = 'WGS84')

r <- st_as_sf(st_read('receivers_in_thesis_data_dec2023.shp'),crs='WGS84') # receiver array used in thesis (n = 54)
r.new4 <- st_as_sf(st_read('receivers_in_thesis_data_dec2023_withEC2020Habdata.shp'), crs = 'WGS84') # made below, in Add habitat data to receivers based on detection radius (211m) - fewer receivers than in 'r' becaus eit removes the ones way outside thestudy area. 
## BRUVs data 

f14 <- 

f18 <- 

## habitat and features

h20 <- st_as_sf(st_read('winter2020habitat_polygonised_forSF_ECdata.shp'), crs = 'WGS84')

jt<-st_zm(st_read('boatlaunchesBimini.kml'),crs='WGS84')%>%
mutate(jetty=seq(1,15,1),.before='geometry')%>%
dplyr::select(-Description) # jetties visible from Google Earth in Bimini. NOT including small single residence jetties. 

labels <- as_labeller(c(lowdensg='Low Density Seagrass',mediumdensg='Medium Density Seagrass',highdensg='High Density Seagrass',vegetated='Vegetated (terrestrial)',bare.urban='Bare &/or Urban (terrestrial)'))
land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	
cmg<-st_as_sf(st_read('habitat_model/centroid_mangroves_north.shp'),crs='WSG84')

## data organisation strutcures- buffers or hexagons

buffs<-st_as_sf(st_read('buffers20192020_withEmilyCourmierhabitat_data4Winter2020.shp'),crs='WGS84')%>%
	mutate(buff_ar=st_area(geometry),buffID=as.factor(buffID))%>%
	distinct()%>%
	group_by(buffID)%>%slice(1) # for reasons I don't know why, r18 has two rows of the same value 
	summary(buffs)

hex <-st_as_sf(st_read('winter2020habitat_hexagon_grid_NOland_ECdata.shp'),crs='WGS84')

## tides

tides2019<-read.csv('tidesbimini/tides2019.csv') # dataset downloaded from NOAA
tides2020<-read.csv('tidesbimini/tides2020.csv') # dataset downloaded from NOAA
tall<-bind_rows(tides2019,tides2020)
tall$Date.Time<-as_datetime(tall$Date.Time,tz="EST")
#create tidephase df - this describes the interval during which the tide phase is peak low or high tide. 18,000 = 2 hours either side of peak time. 
tidephasesdata<-tall%>%dplyr::select(Date.Time,Prediction,Type) %>%
 dplyr::rename(Tidepeak_time=Date.Time,Tidepeak_type=Type, Tide_height=Prediction) %>%
  mutate(Phase_start=Tidepeak_time-dseconds(9000))%>% 
  mutate(Phase_end=Tidepeak_time+dseconds(9000))

#############################################

##########
## - Add habitat data to receivers based on detection radius (211m)
##########

	r2<-r %>%
		st_buffer(211)%>%
		mutate(buff_ar=st_area(.))%>%
		mutate(r.id = seq(1,54,1), .before = 'sharks')
	
	box<-c(xmin=-79.30786,ymin=25.66316,xmax=-79.23485,ymax=25.7796) # values from bounding box of 'h20'
	r2.crop<-st_crop(r2,box)
	
	# Identify which buffers/receivers got cropped out FULLY, visually:

	within <- st_covered_by(r2,st_union(hex))
	within %>% print(n=54) # entirely within, r = 1:10, 14, 22, 24, 26, 27,29, 30, 32, 46
	within <- r2[c(1:10, 14, 22, 24, 26, 27,29, 30, 32, 46),]

	outwith <- st_disjoint(r2, st_union(hex))
	outwith %>% print(n = 54) # r = 11:13, 15:19, 21, 33:45, 47:54
	outside <- r2[c(11:13, 15:19, 21, 33:45, 47:54),]
	outside.r<-outside%>%
				mutate(prop_ldsg=0,prop_mdsg=0,prop_hdsg=0,prop_vegg=0,prop_brurb=0,.before='geometry')%>%
				add_column(prop_unkwn=1,.before='geometry')

 
	partial <- r2[c(20,23,25,28,31),]#partially in, partially out: 20,23,25,28,31
	## use partial and within in the caluclations of proportions of habiitat. after you will subtract the proportions of known from 1 to calculate the proprotion of unknown. Not used in analysis but good to know. 

 	r3 <- r2[c(1:10, 14, 20, 22, 23, 24, 25, 26, 27, 28,29, 30, 31, 32, 46),]

	sf_use_s2(FALSE) # turn speherical geometry off

	r2.join<-st_join(st_make_valid(r3),st_make_valid(h20))%>%
		mutate(hab_geo=st_geometry(st_intersection(st_make_valid(st_geometry(h20)),st_make_valid(st_geometry(r3)))))%>%
		mutate(hab_area=st_area(hab_geo))

	r2.spread<-r2.join%>%
		spread(key= 'feature',value = 'hab_area',fill=0)%>%
		group_by(r.id)%>%
		summarise(across(c('buff_ar','lowdensg','mediumdensg','highdensg','vegetated','bare.urban'),max))%>%
		mutate(prop_ldsg=as.numeric(lowdensg/buff_ar),
			prop_mdsg=as.numeric(mediumdensg/buff_ar),
			prop_hdsg=as.numeric(highdensg/buff_ar),
			prop_vegg=as.numeric(vegetated/buff_ar),
			prop_brurb=as.numeric(bare.urban/buff_ar),
			prop_unkwn=0)%>%
		dplyr::select('r.id','buff_ar','lowdensg','mediumdensg','highdensg','vegetated','bare.urban','prop_ldsg','prop_mdsg','prop_hdsg','prop_vegg','prop_brurb','prop_unkwn',geometry)
	

	r.new2<-bind_rows(r2.spread,outside.r) ## add the r back that were cropped (because outside the extent of known habitat)

	## Calculate distance to the centre point of the mangroves in the north 
	sf_use_s2(TRUE) # turn speherical geometry off
	r.new2<-r.new2%>%
		 mutate(dist_cmg=as.numeric(st_distance(st_centroid(r.new2$geometry),cmg)),.before='geometry')%>%
		 rename(meddensg=mediumdensg)

	## Caluclate the distance to the shore 
	r.new3<-r.new2%>%
		 mutate(dist2shore=as.numeric(st_distance(st_centroid(r.new2$geometry),st_union(land))),.before='geometry')

	## Calculate the distance to the nearest jetty 
	ronly <- r %>%
		mutate(r.id = seq(1,54,1), .before = 'sharks')%>%
		dplyr::select(-sharks)
	nearest2<- tibble(r=seq(1,54,1),jetty=st_nearest_feature(ronly,jt))%>%mutate(r_geo = ronly$geometry)%>%
			mutate(jetty_geo=jt[match(nearest2$jetty,jt$jetty),'geometry'])%>%
			mutate(jetty_id=jt[match(nearest2$jetty,jt$jetty),'Name'])%>%
			st_as_sf()%>%		
			mutate(dist2jetty=st_distance(nearest2$r_geo,nearest2$jetty_geo$geometry,by_element=TRUE))%>%
			rename(r.id=r)%>%
			mutate(r.id=as.character(r.id),dist2jetty=as.numeric(dist2jetty))

		nearest3<-nearest2%>%dplyr::select(r.id,jetty,dist2jetty)
	r.new4 <- r.new3 %>%
		cbind(., nearest3)%>%
		dplyr::select(-r_geo, -sharks, -r.id.1, -lowdensg, -meddensg, -highdensg, -vegetated, -bare.urban)%>%
		mutate(buff_ar = as.numeric(buff_ar))%>%
		rename(arbitr.ID = r.id)

	## Clean up and save 	 

		st_write(r.new4, 'receivers_in_thesis_data_dec2023_withEC2020Habdata.shp', driver = 'ESRI Shapefile')


##########
## - Constructing Predatory Sharks data set 
##########

## Cut 'p' to receivers wihtin study area ('r.set', which includes some with unknown habitat, but removes a handful of way far away receivers). At the same time, join habitat data.
	r.set <- st_as_sf(st_read('receivers_in_thesis_data_dec2023_withEC2020Habdata.shp'), crs = 'WGS84')%>%
		dplyr::select(-arbitr_ID)

	p.sf <- p.sf %>% dplyr::select(time, statn_x, elasmo, locatin, Species, Sex, PCL, depth, geometry)%>%
		rename(station = statn_x)
	pp <- st_intersection(p.sf,r.set) #takes a minute
	ppp <- st_as_sf(pp)
	p <- as_tibble(ppp)

	write_csv(p, 'tempfile_predators_instudysite_withhab_aug2024.csv')

## Assign Tide by time
	# join the detections (data) and the tidephase df, and classify whether the deteciton is inside a low or high tide tidephase interval

	p_withtide<-interval_left_join(p,tidephasesdata,by=c("time"="Phase_start","time"="Phase_end"))%>%st_as_sf()

	# detections can fall outside the window of peak times, classify these as Tidepeak_type='M'
		p_withtide<-p_withtide%>%replace_na(list(Tidepeak_type='M'))
		p_withtide$Tidepeak_type<-as.factor(p_withtide$Tidepeak_type)
		summary(p_withtide) # majority occur within L or H tide. 

	# get rid of what is now unnescessary
		p_withtide<-p_withtide%>%
			dplyr::select(-Tidepeak_time,-Phase_end,-Phase_start,-Tide_height)

## Descriptive statistics on predation pressure at tide times
	p_withtide %>% group_by(Tidepeak_type)%>% tally()
	# 68068 at low tide, all species 
	# 103724 at high tide, all species
	spp_tide <- p_withtide %>%
		filter(Tidepeak_type != 'M')%>%
		group_by(Species,Tidepeak_type)%>%
		tally()
	spp_tide_plot <- ggplot(data=spp_tide,aes(x=Species,y=n,fill=Tidepeak_type)
		)+
		geom_bar(position='stack',stat='identity')+
		scale_fill_manual(values=c('deepskyblue4','deepskyblue1'),name='Tide',labels=c('High','Low'),guide=guide_legend(nrow=1,direction='horizontal'))+
		labs(y = 'Detections (count)')+
		theme_bw()+
		theme(legend.position = 'bottom')
	ggsave(spp_tide_plot,file='lemonspredators_20192020blacktipsANDbulls/descriptive_stats_and_figures/detectionsovertime_atTidepeaks_lemonpredators_20192020.png',device='png',dpi=900,unit='in',height=4,width=4)

######
## - save 
######
	st_write(p_withtide, 'predators_studyareaonly_withTidewithHab_MKthesis20192020.csv', driver = 'CSV')
	st_write(p_withtide, 'predators_studyareaonly_withTidewithHab_MKthesis20192020.shp', driver = 'ESRI Shapefile')



##########
## - Juveniles pointdata - intx: PIT, receiver, grouped by tide 
##########

# simplify some dataframes 
jj <- j %>% dplyr::select(time, PIT, buffID, tidephs)

j2 <- j %>% distinct(buffID, .keep_all = TRUE)
j2 <- j2 %>% dplyr::select(buffID, geometry)

#### update r.set so it have the same buffID values as the juvenile data. Only need to run once, don't run at open, just import the dataframe.
r.set <- st_as_sf(st_read('receivers_in_thesis_data_dec2023_withEC2020Habdata.shp'), crs = 'WGS84')%>%
		dplyr::select(-arbitr_ID)

rr <- st_join(r.set, j2)
r.na <- rr %>% 
		filter(is.na(buffID))%>%
		mutate(buffID = paste0('r', seq(36,54,1)))
r.set <- rr %>% 
		filter(!is.na(buffID))%>%
		bind_rows(., r.na)

st_write(r.set, 'receivers_in_thesis_data_dec2023_withEC2020Habdata_withcorrectBuffID.shp')
r.set <- st_as_sf(st_read('receivers_in_thesis_data_dec2023_withEC2020Habdata_withcorrectBuffID.shp'), crs = 'WGS84')

# remove M tide

j.noM <- jj %>% filter(tidephs != 'M') 

## expand grid by PIT and buffID

tally <- j.noM %>%
		dplyr::select(PIT,buffID,tidephs)%>%
		group_by(PIT,buffID,tidephs)%>%
		tally()

expand <- expand.grid(PIT=unique(j.noM$PIT),buffID=unique(j2$buffID),tidephs=unique(j.noM$tidephs))%>% ## j2$buffIDs, because 2 are used only at M-tides. 
		left_join(.,tally)%>%
		mutate(n=replace_na(n,0))%>%
		dplyr::select(-geometry)%>%
		left_join(.,r.set,by='buffID')%>%
		st_as_sf()
		expand
	stopifnot(nrow(expand) ==560*2) 

######
## - save 
######
	st_write(expand, 'juvenilelemon_pointdatasums_detts_n35receivers_withTIDE_withEC20hab.csv', driver = 'CSV')
	st_write(expand, 'juvenilelemon_pointdatasums_detts_n35receivers_withTIDE_withEC20hab.shp', driver = 'ESRI Shapefile')


##########
## - Pointdata: collate preds and juvs at receivers 
##########

jv <- st_as_sf(st_read('juvenilelemon_pointdatasums_detts_n35receivers_withTIDE_withEC20hab.shp'), crs = 'WGS84')

pr <- st_as_sf(st_read('predators_studyareaonly_withTidewithHab_MKthesis20192020.shp'), crs = 'WGS84')%>% rename(tidephs = Tdpk_ty)

r.set <- st_as_sf(st_read('receivers_in_thesis_data_dec2023_withEC2020Habdata_withcorrectBuffID.shp'), crs = 'WGS84')
r.df <- as_tibble(r.set) %>% dplyr::select(-geometry)

# add consistent buffID to preds from r.set and exclude M-tide
	r <- r.set %>% dplyr::select(buffID)
	pr2 <- st_join(pr, r)%>%
		filter(tidephs != 'M') 
	pr2
	unique(r.set$buffID)

# expand preds to full array and calculate relative predation pressure
		## sum(predator dettections total in year) = dT
		## sum(predator detections at receiver 'x' in year) = dX
		### Pr(predators using the area over a year relative to other receivers) = dX / dT. 
		# spatially comparable, temporally flat 

	tally <- pr2 %>%
		dplyr::select(buffID,tidephs)%>%
		group_by(buffID,tidephs)%>%
		tally()

	pr.expand <- expand.grid(buffID=unique(r.set$buffID),tidephs=unique(pr2$tidephs))%>% 
		left_join(.,tally)%>%
		mutate(n=replace_na(n,0))%>%
		dplyr::select(-geometry)%>%
		left_join(., r.set, by = 'buffID')%>%
		st_as_sf()

		## left off here - next step is calculating rel pr pr.  
	N = nrow(pr.expand)
	pacman::p_load(arm)
	pr3 <- pr.expand %>%
		mutate(
			total = sum(n),
			relPr = as.numeric(n/total), #proportion of detectons at receiver at given tide
			sqzrisk = ((relPr*(N-1))+0.5)/N,
			logit.sqzrisk = logit(sqzrisk),
			zlogit.sqzrisk = (logit.sqzrisk-mean(logit.sqzrisk))/sd(logit.sqzrisk), # standardise and mean centre all variables
			st.lds = (prop_ldsg - mean(prop_ldsg))/sd(prop_ldsg),
			st.mds = (prop_mdsg - mean(prop_mdsg))/sd(prop_mdsg),
			st.hds = (prop_hdsg - mean(prop_hdsg))/sd(prop_hdsg),
			st.distcmg = (dist_cmg - mean(dist_cmg))/sd(dist_cmg),
			st.d2shore = (dist2shore - mean(dist2shore))/sd(dist2shore),
			st.d2jetty = (dist2jetty - mean(dist2jetty))/sd(dist2jetty),
			.before = geometry)%>%
		rename(st.risk = zlogit.sqzrisk)
	pr3

	#####
	## - save
	#####
		st_write(pr3, 'predators_pointdatasums_studyareaonly_withTidewithHab_MKthesis20192020.shp', driver = 'ESRI Shapefile')
		st_write(pr3, 'predators_pointdatasums_studyareaonly_withTidewithHab_MKthesis20192020.csv', driver = 'CSV')

## combine juvs-pointdatasums with preds-pointdatasums - take zlogit.sqzrisk from predators and join based on buffID and tidephs
	names(jv)
	names(pr3)

	pr4 <- as_tibble(pr3) %>% dplyr::select(buffID, tidephs, zlogit.sqzrisk)

	pointdata <- jv%>%
		rename(n.juv = n)%>%
		left_join(.,pr4, by = c('buffID', 'tidephs'))%>%
		mutate( # standardise and mean centre all variables
			st.shark = (n.juv - mean(n.juv))/sd(n.juv),
			st.lds = (prop_ldsg - mean(prop_ldsg))/sd(prop_ldsg),
			st.mds = (prop_mdsg - mean(prop_mdsg))/sd(prop_mdsg),
			st.hds = (prop_hdsg - mean(prop_hdsg))/sd(prop_hdsg),
			st.distcmg = (dist_cmg - mean(dist_cmg))/sd(dist_cmg),
			st.d2shore = (dist2shore - mean(dist2shore))/sd(dist2shore),
			st.d2jetty = (dist2jetty - mean(dist2jetty))/sd(dist2jetty),
			.before = geometry)%>%
		rename(st.risk = zlogit.sqzrisk)
		summary(pointdata) # should be zero NAs

	#####
	## - save
	#####
		st_write(pointdata, 'pointdata_juvlemons_withTidewithHab_MKthesis20192020.shp', driver = 'ESRI Shapefile')
		st_write(pointdata, 'pointdata_juvlemons_withTidewithHab_MKthesis20192020.csv', driver = 'CSV')
		test<-st_as_sf(st_read('pointdata_juvlemons_withTidewithHab_MKthesis20192020.shp'))
		test



