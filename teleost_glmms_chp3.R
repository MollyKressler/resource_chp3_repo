## Models of teleost biodiversity distirbution across study area by habitat type 
#### initial modelling done with habtat data from SOSF pdf. 
#### models will need to be re-run with new year specific habitat data from Emily Courmier (when, TBD) - habitat will be matched to year appropriate remote sensing maps by source (HG vs SD, 2014 for HG and 2019 for SD)

## created 8 February 2023, by Molly Kressler
## updated 13 & 14 February 2023 - BRTs example for gbm.auto package. I need to be able to run the models and export the surfaces (spatially formatted predictions). I'll test the package functionality usign the dataset they provide. 


############ ############ ############ ############ ############ 
## NOTES ### ############ ############ ############ ############

# colors for figures: by season, wet = cadetblue3 and dry = tomato4


############ ############ ############ ############ ############ 
# DATASETS, load at the start ######## ############ ############ 

setwd('/Users/mollykressler/Documents/data_phd/resource_chp3')
pacman::p_load(tidyverse,sf,ggplot2,gridExtra,flextable,vegan,sf,ggsn,lme4,ggeffects,effects,beepr,modelsummary,performance)

# SHAPEFILE: bruvs data and habitat data (using sep22 updated habitat SOSF) joined together, with dist2shore calculated
	joined<-st_as_sf(st_read('bruvs_data_joinedWITHhabitat_feb23.shp'))%>%rename(spp_abundance=spp_bnd,spp_richness=spp_rch,SW_Species=SW_Spcs,SW_Families=SW_Fmls,prop_brs=prp_brs,prop_ldsg=prp_lds,prop_medsg=prp_mds,prop_hdsg=prp_hds,prop_sarg=prp_srg,prop_urb_r=prp_rb_,prop_deep=prop_dp,dist2shore=dst2shr,Sphyraenidae=Sphyrnd,Scaridae=Scarida,Haemulidae=Haemuld,Gerreidae=Gerreid,Belonidae=Belonid,Sparidae=Sparida)%>%select(BRUV,jcode,spp_abundance,spp_richness,SW_Species,SW_Families,Sphyraenidae,Scaridae,Haemulidae,Gerreidae,Belonidae,Sparidae,prop_brs,prop_ldsg,prop_medsg,prop_hdsg,prop_sarg,prop_urb_r,prop_deep,dist2shore)

# DATE FRAME: bruvs data and habitat data (using sep22 updated habitat SOSF) joined together, with dist2shore calculated. Use in models. 

	joined_df<-read.csv('bruvs_DATAFRAME_joinedWITHhabitat_feb23.csv')%>%select(-X)
	joined_df$BRUV<-as.factor(joined_df$BRUV)
	joined_df$Season<-as.factor(joined_df$Season)
	joined_df$SW_Species<-as.numeric(joined_df$SW_Species)
	joined_df$SW_Families<-as.numeric(joined_df$SW_Families)

# SHAPEFILE(S): land and water with grid

	setwd('/Users/mollykressler/Documents/data_phd')
	hab.grid<-st_as_sf(st_read('hexagon_grid_withUPDATEDhabitat_sept22.shp'),crs='WGS84')
	setwd('/Users/mollykressler/Documents/data_phd/resource_chp3')
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')

	hab.noland<-st_difference(hab.grid,st_union(land)) # cut land out of hab.grid. 

############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ 

#### 9 - 10 February 2023 - (A) glm/glmms for response variables. 

### (A) GLMS for response vars: shannon index for families and species, and Gerreidae family counts. Explore some model selection, by dropping certain habitats. 
	## poisson with log link for Gerreidae
	## for SWs? gamma? with negative inverse link. 

sapply(joined_df,class)
summary(joined_df)

# Shannon Index for Species
Hspp<-glm(SW_Species~Season+dist2shore+prop_brs+prop_medsg+prop_sarg+prop_urb_r+prop_deep,data=joined_df,family=gaussian(link='identity'));beep(2)
summary(Hspp)
# Shannon Index for Families
Hfam<-glm(SW_Families~Season+dist2shore+prop_brs+prop_medsg+prop_sarg+prop_urb_r+prop_deep,data=joined_df,family=gaussian(link='identity'));beep(2)
summary(Hfam)

# Gerreidae
gerr<-glm(Gerreidae~Season+dist2shore+prop_brs+prop_medsg+prop_sarg+prop_urb_r+prop_deep,data=joined_df,family=poisson(link='log'));beep(2)
summary(gerr)


# save model RDS 

saveRDS(gerr,'model_outputs/glm_gerreidae_feb23.RDS')
saveRDS(Hspp,'model_outputs/glm_shannonindex_species_feb23.RDS')
saveRDS(Hfam,'model_outputs/glm_shannonindex_families_feb23.RDS')



############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ 

#### 13 & 14 February 2023 - BRTs example for gbm.auto package. I need to be able to run the models and export the surfaces (spatially formatted predictions). I'll test the package functionality usign the dataset they provide. 




















