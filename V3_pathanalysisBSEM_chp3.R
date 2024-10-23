## Path analysis for Risk & Resource BSEM (chapter 3)

## Code and approach here heavily inspired by modelling code from paper Bennett, S., Harris, M. P., Wanless, S., Green, J. A., Newell, M. A., Searle, K. R., & Daunt, F. (2022). Earlier and more frequent occupation of breeding sites during the non‐breeding season increases breeding success in a colonial seabird. Ecology and Evolution, 12(9). https://doi.org/10.1002/ece3.9213.

## created 21 August 2024 :: Molly M Kressler 

# model runs separately for low and high tide. 

###################################
########## RUN AT OPEN ############
###################################

## Load Workspace, local macbook

pacman::p_load(MCMCvis,tidyverse,sf,nimble,devtools,flextable,arm,webshot2,sfdep,sp,spdep,beepr, HDInterval, patchwork, cowplot, tidybayes)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd')

###################################
##########     END     ############
###################################

#############################
### - Pointdata dataframes

## for modelling juvenile shark behaviour 
	pointdata <- read_csv('pointdata_juvlemons_withTidewithHab_MKthesis20192020.csv')%>%
		rename(st.risk = zlogit.sqzrisk)


## for modelling predatory large shark behaviour 

preds_pointdata <- read_csv('predators_pointdatasums_studyareaonly_withTidewithHab_MKthesis20192020.csv')
	# uses the N=54 receiver array, has rows for low and high tide for each receiver. Includes habitat information and predation pressure standardised and mean centred (zlogit.sqzrisk)

#############################
### - Hexagon grid dataframes

## for modelling prey fish distributions & seagrass distributions




























