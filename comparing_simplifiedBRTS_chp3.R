
## gbm and dismo packages allow you to produce simplified models. gbm.simplified function assesses BRTs from the permutations of your total variables, and identifies what level of reduction produces the most parsimonious model.
## While I am unsure of the reasoning of using simplified models in a complex ecological context where I am not sure of the exact interconnectedness of the predctor variables, I am insterested in usign this simplified model option to compare the performance of the full models against simplified models. I'm thinkign this might serve as a check for whether the models are overrun with noise when full. I'm hoping this can help me evaluate the performance of the full models, and the influence of certain habitat predictors on the outcome.

## created by Molly Kressler, 20 February 2023.


########################

###### Start up workspace

setwd('/Users/mollykressler/Documents/data_phd/resource_chp3')
pacman::p_load(tidyverse,sf,ggplot2,gridExtra,flextable,sf,ggsn,dismo,gbm)

###### 'Full' models, from file 'teleost_glmms_chp3.R' ######

	sppfull<-readRDS('model_RDS/SW_Species_brt_tc3_lr001_gaussian.RDS')
	famfull<-readRDS('model_RDS/SW_Families_brt_tc3_lr001_gaussian.RDS')
	gerrfull<-readRDS('model_RDS/Gerreidae_brt_tc3_lr0001_poisson.RDS')


###### DATE FRAME for modelling: bruvs data and habitat data (using sep22 updated habitat SOSF) joined together, with dist2shore calculated. Use in models. 

	joined_df<-read.csv('bruvs_DATAFRAME_joinedWITHhabitat_feb23.csv')%>%dplyr::select(-X)
	joined_df$BRUV<-as.factor(joined_df$BRUV)
	joined_df$Season<-as.factor(joined_df$Season)
	joined_df$SW_Species<-as.numeric(joined_df$SW_Species)
	joined_df$SW_Families<-as.numeric(joined_df$SW_Families)

###### Useful SHAPEFILE(S): land and water with grid

	setwd('/Users/mollykressler/Documents/data_phd')
	hab.grid<-st_as_sf(st_read('hexagon_grid_withUPDATEDhabitat_sept22.shp'),crs='WGS84')
	setwd('/Users/mollykressler/Documents/data_phd/resource_chp3')
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')

	hab.noland<-st_as_sf(st_read('hexagon_grid_study_site_with_habitatFromSOSF_forchp3BRUVS.shp'),crs='WGS84')
		# how I made it: <-st_as_sf(st_difference(hab.grid,st_union(land))) # cut land out of hab.grid. 

###### Dateframe and Simple Feature for predicting #####

	df4preds<-read.csv('df_for_prediction_chp3BRUVS_habitat_noland.csv')%>%dplyr::select(-X)
		df4preds$jcode<-as.factor(df4preds$jcode)
		df4preds$Season<-as.factor(df4preds$Season)

	sf4preds<-st_as_sf(st_read('sf_for_predictions_fromBRTs_feb23.shp'),crs='WGS84')
		sf4preds$jcode<-as.factor(sf4preds$jcode)
		sf4preds$Season<-as.factor(sf4preds$Season)
		# predictions for full models are stored in sf4preds, as preds_SWf, _SWsp, _gerr.

########################

###### Simplify full models using gbm.simplify. n.drops = maximum vars the function can drop. ######
# as a starting point use the same settigns for lr, tc, and bf

	# SW Spp
		sp_simp<-gbm.simplify(sppfull,n.drops=7) # deviance in full model = 0.2681
		# suggests dropping two vars: 
		sp_simp$pred.list[[2] # this notation prints the specifric list of variables it suggests 
		# it suggested dropping bare sand and season. 
		sp_simpby2<-gbm.step(joined_df,gbm.x=sp_simp$pred.list[[2]],gbm.y=c('SW_Species'),tree.complexity=3,learning.rate=0.001,bag.fraction=0.75,family='gaussian',plot.main = TRUE) 
			# 2500 trees, mean deviance = 0.336 - 100 more trees than the full model.
		# what if I only drop bare sande (SW Familiy simplified drops bare sand)
		sp_simpby1<-gbm.step(joined_df,gbm.x=sp_simp$pred.list[[1]],gbm.y=c('SW_Species'),tree.complexity=3,learning.rate=0.001,bag.fraction=0.75,family='gaussian',plot.main = TRUE) 
			# new deviance = 0.213, n.trees=220
		saveRDS(sp_simpby1,'model_RDS/simplifiedby1_BRT_SWspecies_feb23_tc3_lr001_bf75.RDS')
			
	# SW F
		fam_simp<-gbm.simplify(famfull,n.drops=7) # deviance in full model = 0.2675, 2300 trees
		# suggests dropping two vars: 
		fam_simp$pred.list[[2]]
		# it suggested dropping bare sand. 
		fam_simpby2<-gbm.step(joined_df,gbm.x=fam_simp$pred.list[[2]],gbm.y=c('SW_Families'),tree.complexity=3,learning.rate=0.001,bag.fraction=0.75,family='gaussian',plot.main = TRUE) 
			# new mean deviance = 0.13, n.tress=950, full model has 2250. 
		saveRDS(fam_simpby2,'model_RDS/simplifiedby1_BRT_SWFamilies_feb23_tc3_lr001_bf75.RDS')

########### you left off here. just ran fam_simpby2. RDS not yet saved because the ntrees is so low and you might want to play with the LR. Will need to re-do all the below plotting and evaluating of simplifed by 2 families model, bc it had the wrong RDS. 


	# Gerr



	gerr.simp<-gbm.simplify(gerr2,n.drops=5) # n.drops defines how many the simplify function CAN drop at a MAX.
	# estimated that I should drop 2 variables. 
	# re-run gbm.step but instead of listing the vars, I refer the gbm.x argument to the listed object form the simp model. 
	gerr.reducedby5<-gbm.step(joined_df,gbm.x=gerr.simp$pred.list[[5]],gbm.y=c('Gerreidae'),tree.complexity=3,learning.rate=0.0001,bag.fraction=0.5,family='poisson',plot.main = TRUE) # 2800 trees, mean resid dev = 10.591
	summary(gerr.reduced)

###### Evaluate models, plot fitted values, and explore the interactions 

	hab.labels<-(c('dist2shore'='Dist. to Shore (m)', 'prop_brs'='Prop. of Bare Sand','prop_ldsg'='Prop. of \n\ Low Density \n\ Seagrass','prop_medsg'='Prop. of  \n\ Medium Density \n\ Seagrass','prop_hdsg'='Prop. of  \n\ High Density \n\ Seagrass','prop_sarg'='Prop. of Sargassum','prop_urb_r'='Prop. of \n\ Urban & Rocky','prop_deep'='Prop. of \n\ Deep Water'))


	# SW Sp
		infl.simp.spp<-sp_simpby2$contributions

		swspp.relinf<-infl.simp.spp%>%mutate(var=fct_reorder(var,rel.inf))%>%ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+geom_bar(stat='identity')+scale_fill_distiller(direction=1,palette='Oranges',limits=c(0,35),guide='none')+theme_bw()+coord_flip()+scale_x_discrete(labels=hab.labels)+ylab('Relative Influence')+xlab(NULL)
		ggsave(swspp.relinf,file='figures+tables/simplifiedBRTs/simplifiedby2_SWspp_BRT_feb23_relativeinfluence_ofvars.png',device='png',units='in',width=6,height=6,dpi=900)

		fittedfx.swspsimp<-gbm.plot(sp_simpby2, plot.layout=c(3,3),write.title=FALSE) 

		find.int.swspp<-gbm.interactions(sp_simpby2) # run the function, calculates
		find.int.swspp$interactions # print the interacton matrix
		find.int.swspp$rank.list # prints just the list of non-zero interactions. Make these into flextables. 
		# make a table from it 
		swpp.interactions<-find.int.swspp$rank.list%>%dplyr::select(var1.names,var2.names,int.size)%>%rename('Variable 1'='var1.names','Variable 2' = 'var2.names','Size of Interaction'=int.size)%>%flextable()%>%theme_alafoli()%>%align(align = 'center', part = 'all')%>%font(fontname = 'Times', part = 'all')%>%fontsize(size = 11, part = 'all')%>%color(color='black',part='all')%>%autofit()
			save_as_image(swpp.interactions,'interactions_above05_SWsppBRT_simplifiedby2_feb23.png',webshot='webshot')
		# only dropping baresand 
		infl.simp.spp<-sp_simpby1$contributions

		swspp.relinf<-infl.simp.spp%>%mutate(var=fct_reorder(var,rel.inf))%>%ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+geom_bar(stat='identity')+scale_fill_distiller(direction=1,palette='Oranges',limits=c(0,35),guide='none')+theme_bw()+coord_flip()+scale_x_discrete(labels=hab.labels)+ylab('Relative Influence')+xlab(NULL)
		ggsave(swspp.relinf,file='figures+tables/simplifiedBRTs/simplifiedby1_SWspp_BRT_feb23_relativeinfluence_ofvars.png',device='png',units='in',width=6,height=6,dpi=900)

		fittedfx.swspsimp<-gbm.plot(sp_simpby1, plot.layout=c(3,3),write.title=FALSE) 

		find.int.swspp<-gbm.interactions(sp_simpby1) # run the function, calculates
		find.int.swspp$interactions # print the interacton matrix
		find.int.swspp$rank.list # prints just the list of non-zero interactions. Make these into flextables. 
		# make a table from it 
		swpp.interactions<-find.int.swspp$rank.list%>%dplyr::select(var1.names,var2.names,int.size)%>%rename('Variable 1'='var1.names','Variable 2' = 'var2.names','Size of Interaction'=int.size)%>%flextable()%>%theme_alafoli()%>%align(align = 'center', part = 'all')%>%font(fontname = 'Times', part = 'all')%>%fontsize(size = 11, part = 'all')%>%color(color='black',part='all')%>%autofit()
			save_as_image(swpp.interactions,'interactions_above05_SWsppBRT_simplifiedby1_feb23.png',webshot='webshot')

	# SW Families
		infl.simp.fam<-fam_simpby1$contributions

		fam.relinf<-infl.simp.fam%>%mutate(var=fct_reorder(var,rel.inf))%>%ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+geom_bar(stat='identity')+scale_fill_distiller(direction=1,palette='Oranges',limits=c(0,35),guide='none')+theme_bw()+coord_flip()+scale_x_discrete(labels=hab.labels)+ylab('Relative Influence')+xlab(NULL)
		ggsave(fam.relinf,file='figures+tables/simplifiedBRTs/simplifiedby1_SWfam_BRT_feb23_relativeinfluence_ofvars.png',device='png',units='in',width=6,height=6,dpi=900)

		fittedfx.famsimp<-gbm.plot(fam_simpby1, plot.layout=c(3,3),write.title=FALSE) 

		find.int.fam<-gbm.interactions(fam_simpby1) # run the function, calculates
		find.int.fam$interactions # print the interacton matrix
		find.int.fam$rank.list # prints just the list of non-zero interactions. Make these into flextables. 
		# make a table from it 
		fam.interactions<-find.int.fam$rank.list%>%dplyr::select(var1.names,var2.names,int.size)%>%rename('Variable 1'='var1.names','Variable 2' = 'var2.names','Size of Interaction'=int.size)%>%flextable()%>%theme_alafoli()%>%align(align = 'center', part = 'all')%>%font(fontname = 'Times', part = 'all')%>%fontsize(size = 11, part = 'all')%>%color(color='black',part='all')%>%autofit()
			save_as_image(fam.interactions,'interactions_above05_SWfamBRT_simplifiedby1_feb23.png',webshot='webshot')



###### PREDICT FROM MODELS to new data, spatial
	season.labels<-as_labeller(c('D'='Dry','W'='Wet'))

	# SW Sp
		# simpby2 --this model removed Season so need to remove Season column and reduce back to only one row for each jcode. 
			sf4preds.noSeasons<-sf4preds%>%filter(Season!='D')%>%dplyr::select(-Season)
			df4preds.noSeasons<-df4preds%>%filter(Season!='D')%>%dplyr::select(-Season)

			simpSws_pr<-predict.gbm(sp_simpby2,df4preds.noSeasons,n.trees=sp_simpby2$gbm.call$best.trees,type='response')	# one item per jcode. 
		#make it a df, and bind to df4preds
			p2spp<-as.data.frame(simpSws_pr)
			preds2_SWspp<-bind_cols(df4preds.noSeasons,p2spp)
			head(preds2_SWspp)
		# this model removed Season so need to remove Season column and reduce back to only one row for each jcode. 
		# add predictions to sf4preds
			sf4preds.noSeasons<-sf4preds.noSeasons%>%mutate(simpSws_pr=preds2_SWspp$simpSws_pr,.before='geometry')
		# plot
			simple.SWsp.plot<-ggplot()+geom_sf(data=sf4preds.noSeasons,aes(fill=simpSws_pr),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0.5,1.7),guide=guide_colourbar(title=' Species \n\ Diversity \n\ Index'))+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
			ggsave(simple.SWsp.plot,file='figures+tables/simplifiedby2_BRTSWspecies_preds_feb23.png',device='png',units='in',height=4,width=5,dpi=1000)


		# simpby1 -- this one has season, but not baresand. 

				preds_simpby1sp<-predict.gbm(sp_simpby1,df4preds,n.trees=sp_simpby1$gbm.call$best.trees,type='response')	# one item per jcode. 
				#make it a df, and bind to df4preds
				p2spp<-as.data.frame(preds_simpby1sp)
				preds2_SWspp<-bind_cols(df4preds,p2spp)
			# filter predictions by season (W/D)
				spppreds.dry<-preds2_SWspp%>%filter(Season!='W')
				spppreds.wet<-preds2_SWspp%>%filter(Season!='D')
			# add predictions to sf4preds, filtered by season.
				sf4predsDRY<-sf4preds%>%filter(Season=='D')%>%mutate(pr_simp1sp=spppreds.dry$preds_simpby1sp,.before='geometry')
				sf4predsWET<-sf4preds%>%filter(Season=='W')%>%mutate(pr_simp1sp=spppreds.wet$preds_simpby1sp,.before='geometry')
			# (re)combine separate seasons
				sf4preds<-bind_rows(sf4predsWET,sf4predsDRY)
				sf4preds$Season<-as.factor(sf4preds$Season)
			# plot
			swspp.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=pr_simp1sp),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0.5,1.7),guide=guide_colourbar(title=' Species \n\ Diversity \n\ Index'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
			ggsave(swspp.plot,file='figures+tables/simplifiedby1_BRTSWspecies_preds_feb23.png',device='png',units='in',height=8,width=10,dpi=1000)

	# Sw Families

		preds_simpby1fam<-predict.gbm(fam_simpby1,df4preds,n.trees=fam_simpby1$gbm.call$best.trees,type='response')	# one item per jcode. 
		#make it a df, and bind to df4preds
		p2spp<-as.data.frame(preds_simpby1fam)
		preds2_SWfam<-bind_cols(df4preds,p2spp)

		# filter predictions by season (W/D)
			spppreds.dry<-preds2_SWfam%>%filter(Season!='W')
			spppreds.wet<-preds2_SWfam%>%filter(Season!='D')
		# add predictions to sf4preds, filtered by season.
			sf4predsDRY<-sf4preds%>%filter(Season=='D')%>%mutate(pr_simp1F=spppreds.dry$preds_simpby1fam,.before='geometry')
			sf4predsWET<-sf4preds%>%filter(Season=='W')%>%mutate(pr_simp1F=spppreds.wet$preds_simpby1fam,.before='geometry')
		# (re)combine separate seasons
			sf4preds<-bind_rows(sf4predsWET,sf4predsDRY)
			sf4preds$Season<-as.factor(sf4preds$Season)
		# plot
		swFAM.plot<-ggplot()+geom_sf(data=sf4preds,aes(fill=pr_simp1F),col=NA)+scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0,0.5),guide=guide_colourbar(title=' Families \n\ Diversity \n\ Index'))+facet_grid(cols=vars(Season),labeller=season.labels)+theme_minimal()+geom_sf(data=land,fill='gray98')+theme_bw()
		ggsave(swspp.plot,file='figures+tables/simplifiedby1_BRTSWfamilies_preds_feb23.png',device='png',units='in',height=8,width=10,dpi=1000)

















