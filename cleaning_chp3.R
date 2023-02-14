### For cleaning and processing of data for chapter 3 - resource selection model, based on habittat seleciton of prey and hypoerprior of shark habitat selection. 

## created 3 Februry 2023, by Molly Kressler
## updated: 7 Feb. 2023 
## updated: 8 Feb. 2023

#### HABITAT JOINING will need updating when I get the updated year-specific remote sensed habitat maps from Emily Courmier. Acquisition TBD. 


pacman::p_load(tidyverse,sf,ggplot2,gridExtra,flextable,vegan,sf,ggsn)
setwd('/Users/mollykressler/Documents/data_phd/resource_chp3')
# colors for figures: by season, wet = cadetblue3 and dry = tomato4

############ ############ ############ ############ ############ 
############ ############ # DATASETS # ############ ############ 
############ ########## LOAD @ THE START  ######### ############ 
############ ############ ############ ############ ############ 

# SHAPEFILE: bruvs data and habitat data (using sep22 updated habitat SOSF) joined together, with dist2shore calculated
	joined<-st_as_sf(st_read('bruvs_data_joinedWITHhabitat_feb23.shp'))%>%rename(spp_abundance=spp_bnd,spp_richness=spp_rch,SW_Species=SW_Spcs,SW_Families=SW_Fmls,prop_brs=prp_brs,prop_ldsg=prp_lds,prop_medsg=prp_mds,prop_hdsg=prp_hds,prop_sarg=prp_srg,prop_urb_r=prp_rb_,prop_deep=prop_dp,dist2shore=dst2shr,Sphyraenidae=Sphyrnd,Scaridae=Scarida,Haemulidae=Haemuld,Gerreidae=Gerreid,Belonidae=Belonid,Sparidae=Sparida)%>%select(BRUV,jcode,spp_abundance,spp_richness,SW_Species,SW_Families,Sphyraenidae,Scaridae,Haemulidae,Gerreidae,Belonidae,Sparidae,prop_brs,prop_ldsg,prop_medsg,prop_hdsg,prop_sarg,prop_urb_r,prop_deep,dist2shore)

# DATE FRAME: bruvs data and habitat data (using sep22 updated habitat SOSF) joined together, with dist2shore calculated

	joined_df<-read.csv('bruvs_DATAFRAME_joinedWITHhabitat_feb23.csv')

# SHAPEFILE(S): land and water with grid

	setwd('/Users/mollykressler/Documents/data_phd')
	hab.grid<-st_as_sf(st_read('hexagon_grid_withUPDATEDhabitat_sept22.shp'),crs='WGS84')
	setwd('/Users/mollykressler/Documents/data_phd/resource_chp3')
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')

	hab.noland<-st_difference(hab.grid,st_union(land)) # cut land out of hab.grid. 


############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ 


############ ############ ############ ############ ############ 
############ 3 February 2023 - collating and combining BRUVS data from Sarah Driscoll and Henriette Grimmel 

## Combining BRUVS data 
	## BRUVS data from Sarah Driscoll

		sd<-read.csv('raw_and_outdated/driscoll_bruvs_data_raw.csv')%>% rename(BRUV_ID=SiteName,SW_Species=Shannon,spp_richness=spp_rich)%>%add_column(source='SD')
		names(sd)

		sd2<-sd%>%select(BRUV_ID, Latitude, Longitude, spp_abundance, spp_richness,SW_Species, Season,source) #common variables to sd and hg

		spp.sd<-names(sd[13:90])

		## BRUVS data from Henriette Grimmel and Rob Bullock 

		hg<-read.csv('raw_and_outdated/grimmel_bruvs_data_raw.csv')%>%add_column(source='HG')%>%rename(spp_abundance=MaxN_Species)
		hg<-hg[,c(1:128,257)]
		names(hg)
		hg2<-hg%>%select(BRUV_ID, Latitude, Longitude, spp_abundance, spp_richness,SW_Species, Season,source) # common variables

	## BRUVS data combined from SD and HG, with families. 

		hg.fam<-hg%>%select(BRUV_ID,Latitude, Longitude, spp_abundance, spp_richness,SW_Species, Season,source,Sphyraenidae,Scaridae,Haemulidae,Gerreidae,Belonidae,Sparidae)%>%add_column(BRUV=paste0(.$source,.$BRUV_ID),.before='BRUV_ID')
		head(hg.fam,1)

		sd.fam<-sd%>%dplyr::select(BRUV_ID, Latitude, Longitude, spp_abundance, spp_richness,SW_Species, Season,source,Gerres.sp.,Sphyraena.sp.,Scaridae.sp.,Belonidae.sp.,Sparidae.sp.)%>%mutate(Haemulidae=rowSums(across(starts_with('Haemulon'))))%>%rename(Gerreidae=Gerres.sp.,Sphyraenidae=Sphyraena.sp.,Belonidae=Belonidae.sp.,Sparidae=Sparidae.sp.,Scaridae=Scaridae.sp.)%>%add_column(BRUV=paste0(.$source,.$BRUV_ID),.before='BRUV_ID')
		head(sd.fam,1)

		data<-bind_rows(hg.fam,sd.fam)
		head(data,2)
		names(data)

	write.csv(data,'bruvs_data_cleaned_sd_and_hg.csv')

## summary plots - updated 8 February 2023 - Shannon index (H) plots 

	data<-read.csv('bruvs_data_cleaned_sd_and_hg.csv')%>%select(-X)

	richness<-ggplot(data=data,aes(x=spp_richness))+geom_histogram(binwidth=1,col='violetred4',fill='violetred3')+xlab('Richness')+ylab(NULL)+theme_bw()
	abundance<-ggplot(data=data,aes(x=spp_abundance))+geom_histogram(binwidth=20,col='violetred4',fill='violetred3')+xlab('Abundance')+ylab(NULL)+theme_bw()
	Hspp<-ggplot(data=data,aes(x=SW_Species))+geom_histogram(binwidth=0.2,col='violetred4',fill='violetred3')+xlab('H, Species')+ylab(NULL)+theme_bw()
	Hfam<-ggplot(data=data,aes(x=SW_Families))+geom_histogram(binwidth=0.1,col='violetred4',fill='violetred3')+xlab('H, Families')+ylab(NULL)+theme_bw()
	sphrnaenidae<-ggplot(data=data,aes(x=Sphyraenidae))+geom_histogram(binwidth=1,col='cadetblue4',fill='cadetblue3')+ylab(NULL)+theme_bw()
	gerreidae<-ggplot(data=data,aes(x=Gerreidae))+geom_histogram(binwidth=5,col='cadetblue4',fill='cadetblue3')+ylab(NULL)+theme_bw()
	belonidae<-ggplot(data=data,aes(x=Belonidae))+geom_histogram(binwidth=5,col='cadetblue4',fill='cadetblue3')+ylab(NULL)+theme_bw()
	sparidae<-ggplot(data=data,aes(x=Sparidae))+geom_histogram(binwidth=1,col='cadetblue4',fill='cadetblue3')+ylab(NULL)+theme_bw()
	scaridae<-ggplot(data=data,aes(x=Scaridae))+geom_histogram(binwidth=5,col='cadetblue4',fill='cadetblue3')+theme_bw()+ylab(NULL)
	haemulidae<-ggplot(data=data,aes(x=Haemulidae))+geom_histogram(binwidth=5,col='cadetblue4',fill='cadetblue3')+ylab(NULL)+theme_bw()

	histograms.bruvs<-grid.arrange(richness,abundance,Hspp,Hfam,sphrnaenidae,gerreidae,belonidae,sparidae,scaridae,haemulidae,ncol=5)
	ggsave(histograms.bruvs,file='figures+tables/histograms_bruvs_families_richness_abundance.png',device='png',units='in',height=8,width=16,dpi=800)


## table that summarises what was seen on wach BRUV and other 'biometric' equivalent stats: family counts, richness and abundance, season, shannon.

	table.info<-data%>%select(BRUV,spp_abundance,spp_richness,SW_Species,Sphyraenidae,Scaridae,Haemulidae,Gerreidae,Belonidae,Sparidae)%>%filter(spp_abundance!='NA')
	table<-flextable(table.info)%>%set_header_labels(BRUV='ID',spp_abundance='Abundance',spp_richness='Richness',SW_Species='Shannon Index')%>%theme_alafoli()%>%align(align = 'center', part = 'all')%>%font(fontname = 'Times', part = 'all')%>%fontsize(size = 10, part = 'all')%>%color(color='black',part='all')%>%autofit()
	#save_as_image(table,'figures+tables/bruvs_table_richness_abundance_familiescounts.png',webshot='webshot2')



############ ############ ############ ############ ############ 
############ 7 February 2023 - (A) calculate Shannon index for families. SPATIAL: (B) import BRUVS data into sf. (C) join habitat to BRUVS locations. (D) calculate distance from nearest shore. (E) some figures of BRUV location and a table with BRUV data and habitat data 

## import cleaned combined BRUVS data

	data<-read.csv('bruvs_data_cleaned_sd_and_hg.csv')%>%select(-X)%>%filter(BRUV_ID!='NA')
	names(data)


## (A) calculate Shannon index for families per bruvs using 'vegan' package - the proportional abundance of families in a given location. Use the diversity function. 

	data(BCI)
	head(BCI)

	# subset the data to only have the families per row.
	fams<-data[10:15]
	head(fams,1)
	nrow(data)

	data2<-data%>%mutate('SW_Families'=diversity(fams,'shannon'),.after='SW_Species')
	head(data2,2)
	summary(data2)

	## update the saved data file
	#write.csv(data2,'bruvs_data_cleaned_sd_and_hg.csv')


## (B) import BRUVS data into sf. 

	data<-read.csv('bruvs_data_cleaned_sd_and_hg.csv')%>%select(-X)
		head(data,2)
	br<-st_as_sf(data,coords=c('Longitude','Latitude'),crs='WGS84')
	br
	ggplot()+geom_sf(data=br)+theme_bw() # all good. 


## (C) join habitat data to BRUVS location [temporary - updated year specific data coming from Emily soon]
### some BRUVS occur close to the edge of the land mass which is cut out of the habitat grid. So need to use an sf function which joins the point to the nearest grid cell -- try st_overlaps, st_touches, st_crosses. 

	setwd('/Users/mollykressler/Documents/data_phd')
	hab.grid<-st_as_sf(st_read('hexagon_grid_withUPDATEDhabitat_sept22.shp'),crs='WGS84')
	setwd('/Users/mollykressler/Documents/data_phd/resource_chp3')
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')

	hab.noland<-st_difference(hab.grid,st_union(land)) # cut land out of hab.grid. 

	#ggplot()+geom_sf(data=hab.noland,alpha=0.3,col='goldenrod2')+geom_sf(data=br,fill='cadetblue3',col='cadetblue3')+geom_sf(data=land,fill='grey82',col='grey70',alpha=0.3)+theme_bw() # all good. 

	joined<-st_join(br,hab.noland,join=st_nearest_feature)
	joined

## (D) calculate distance to nearest shore 

	joined<-joined%>%mutate(dist2shore=as.numeric(st_distance(st_centroid(joined),st_union(land))),.after='dist_cmg')
	# save 
	st_write(joined,'bruvs_data_joinedWITHhabitat_feb23.shp',driver='ESRI Shapefile')
		#joined<-st_as_sf(st_read('bruvs_data_joinedWITHhabitat_feb23.shp'))


## (E) some figures of BRUV location and a table with BRUV data and habitat data 

	ggplot()+geom_sf(data=land,fill='grey72',col='grey72')+geom_sf(data=joined,col='cadetblue4',fill='cadetblue4',alpha=0.7,size=1.5,pch=21)+theme(legend.position=c(.5,.5),legend.text=element_text(size=20),axis.text=element_text(size=20))+north(data=joined,location='bottomright',scale=0.11,symbol=12)+theme_bw() 
	# wet and dry seasons are different colours
	wetanddry<-ggplot()+geom_sf(data=land,fill='grey72',col='grey72')+geom_sf(data=joined,aes(fill=Season),col='grey45',alpha=0.5,size=1.5,pch=21)+theme(legend.position=c(.5,.5),legend.text=element_text(size=20),axis.text=element_text(size=20))+scale_fill_manual(values=c('tomato4','cadetblue3'),labels=c('Dry','Wet'))+north(data=joined,location='bottomright',scale=0.11,symbol=12)+theme_bw() 

	ggsave(wetanddry,file='figures+tables/bruvs_locations_wetanddry_seasons_onmap.png',device='png',units='in',height=4,width=4,dpi=1000)

	table.info<-joined%>%select(BRUV,Season,spp_abundance,spp_richness,SW_Species,SW_Families,Sphyraenidae,Scaridae,Haemulidae,Gerreidae,Belonidae,Sparidae,prop_brs,prop_ldsg,prop_medsg,prop_hdsg,prop_sarg,prop_urb_r,prop_deep,dist2shore)
	tab2<-as.data.frame(table.info)%>%select(-geometry)
	names(tab2)
	table<-flextable(tab2)%>%set_header_labels(BRUV='ID',spp_abundance='Abundance',spp_richness='Richness',SW_Species='H, species',SW_Families='H, families')%>%theme_alafoli()%>%align(align = 'center', part = 'all')%>%font(fontname = 'Times', part = 'all')%>%fontsize(size = 10, part = 'all')%>%color(color='black',part='all')%>%colformat_double(digits=3)%>%set_header_labels(prop_brs='Bare Sand',prop_ldsg='Low Density Seagrass',prop_medsg='Medium Density Seagrass',prop_hdsg='High Density Seagrass',prop_sarg='Sargassum',prop_urb_r='Urban & Rocky',prop_deep='Deep Water',dist2shore='Distance to Nearest Shore (m)')%>%autofit()
	save_as_image(table,'figures+tables/bruvs_table_data_joinedwithhabitat.png',zoom=2,webshot='webshot2')



############ 










































