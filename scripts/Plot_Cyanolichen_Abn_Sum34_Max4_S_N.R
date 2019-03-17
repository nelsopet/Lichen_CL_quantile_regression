############################################################################################################################################
#Abundance of Cyanolichen vs Nitrogen 
################################################################################################################

	## Input files are generated in Calculate Abundances script. I have to switch each of these sources depending on which measure of Abundance wanted.
	## We decided to use the maximum Abundance because it was the most straightforward to explain with little to no loss in model fit.
## Instead of changing SOURCE, I should build a separate pipeline for each measure of Abundance (summing 3s and 4s, summing log transformed values and only max value)
	SOURCE<-LichDb_sFncGrpSens_All_Abun_N

XVAR<-c("N"                                                               
 ,"poly(N, 2)"                                                      
 ,"N+maxaug_c+mindec_c+precip_cm+CMD"         
 ,"poly(N, 2)+maxaug_c+mindec_c+precip_cm+CMD"
 ,"maxaug_c+mindec_c+precip_cm+CMD"
 ,"N+maxaug_c"                                                      
 ,"N+mindec_c"                                                      
 ,"N+precip_cm"                                                     
 ,"N+continen"
,"N+CMD"
,"maxaug_c"
,"mindec_c"
,"precip_cm"
,"continen"
,"CMD") 

XVAR_names<-c("N","N_poly","N_clim","N_poly_clim","clim","N_maxaug_c","N_mindec_c","N_precip_cm","N_continen","N_CMD", "Max C Aug Temp","Min C Dec Temp","Precip_cm","Continentality","CMD");

Abn_Cya_sum34_max4_N_Mods<-paste("Abn_Cya_sum34_max4 ~ ",paste(XVAR),sep="")

Abn_Cya_sum34_max4_N_Formula<-lapply(Abn_Cya_sum34_max4_N_Mods,as.formula)

Abn_Cya_sum34_max4_N_Mods_run<-lapply(Abn_Cya_sum34_max4_N_Formula,function (x) rq(x, tau=0.9, data=SOURCE))


	## Make null models
Abn_Cya_sum34_max4_N_intercept<-rq(formula = Abn_Cya_sum34_max4 ~ 1, tau = 0.9, data = SOURCE, model = T)


	## Grab each model fit but the intercept is hard coded 

Abn_Cya_sum34_max4_N_Mods_fit<-mapply(function (y) 1-Abn_Cya_sum34_max4_N_Mods_run[[y]]$rho/Abn_Cya_sum34_max4_N_intercept$rho, 1:length(XVAR))


	## makes a list of mod names

Abn_Cya_sum34_max4_N_Mod_names<-lapply(XVAR_names, function (x) paste("Abn_Cya_sum34_max4_",paste(x),sep=""))

names(Abn_Cya_sum34_max4_N_Mods_run)<-Abn_Cya_sum34_max4_N_Mod_names


	## Save model stats as a spreadsheet of manuscript

Abn_Cya_sum34_max4_N_Mods_AIC<-lapply(Abn_Cya_sum34_max4_N_Mods_run,AIC)
Abn_Cya_sum34_max4_N_Lichen_Abundance_vs_N_Stats<-cbind(XVAR_names,Abn_Cya_sum34_max4_N_Mods_AIC,Abn_Cya_sum34_max4_N_Mods_fit,paste(Abn_Cya_sum34_max4_N_Formula))
write.csv(Abn_Cya_sum34_max4_N_Lichen_Abundance_vs_N_Stats,"Abundance_Abn_Cya_sum34_max4_N_all_models.csv")

Reduc_Mods_forTable<-c("N","N_poly","N_poly_clim")
Abn_Cya_sum34_max4_N_Lichen_Abundance_vs_N_Stats_reduc<-cbind(XVAR_names,Abn_Cya_sum34_max4_N_Mods_fit,Abn_Cya_sum34_max4_N_Mods_AIC,paste(Abn_Cya_sum34_max4_N_Formula))
Abn_Cya_sum34_max4_N_Lichen_Abundance_vs_N_Stats_reduc<-subset(Abn_Cya_sum34_max4_N_Lichen_Abundance_vs_N_Stats_reduc, XVAR_names %in% Reduc_Mods_forTable)
write.csv(Abn_Cya_sum34_max4_N_Lichen_Abundance_vs_N_Stats_reduc,"Abundance_Abn_Cya_sum34_max4_N_all_models_reduc.csv")


	## Make predicted data files
	## Use a model file ... Abn_Cya_sum34_max4_Mods_run as input

MyFunc2<-function (x) {as.data.frame(predict(x, SOURCE, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000))}

	##This produces a matrix with a list of 10 objects with 3 observations, each with 5322 entries and 10 columns. 
Abn_Cya_sum34_max4_N_Mod_pred<-mapply(MyFunc2,Abn_Cya_sum34_max4_N_Mods_run)%>%as.data.frame()

	## This works but is hard coded for each predictor file. Have to change working directory each time

#pdf(paste("Abn_Cya_sum34_max4_N_Lichen_Abundance_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),".pdf",sep=""))
#
#MyPlot<-function (x)
#    
#{
#    x<-Abn_Cya_sum34_max4_N_Mod_pred[i]
#    x<-unlist(x, recursive = F, use.names =T)
#    fit<-unlist(x[1], recursive = F, use.names =T)
#    low<-unlist(x[2], recursive = F, use.names =T)
#    high<-unlist(x[3], recursive = F, use.names =T)
#    plot(Abn_Cya_sum34_max4~N, data=SOURCE, xlab="CMAQ N kg/ha/yr", ylab="Abundance Abn_Cya_sum34_max4 Lichens", main="Abundance of  Abn_Cya_sum34_max4 lichen species vs N deposition")
#    points(SOURCE$N, fit, col='blue', pch= 19)
#    points(SOURCE$N, low, col='aquamarine4', pch= 4, cex=0.25) 
#    points(SOURCE$N, high, col='aquamarine4', pch= 4, cex=0.25)	
#    text(max(SOURCE$N)*0.4,max(SOURCE$Abn_Cya_sum34_max4)*0.8,labels=paste("R1=",print(round(Abn_Cya_sum34_max4_N_Mods_fit[i],2))))
#    text(max(SOURCE$N)*0.4,max(SOURCE$Abn_Cya_sum34_max4)*0.75,labels=paste(as.character(Abn_Cya_sum34_max4_N_Formula[i])), cex=0.7)
#    legend('topright', legend=c("Raw Data","Fitted Values", "Prediction Interval"), col=c("black","blue", "aquamarine4"), pch=16)
# 
#}	
#
#for (i in 1:length(XVAR)) MyPlot()
#
#dev.off()
    #tiff(paste("Abn_Cya_sum34_max4_N_Cyanolichen_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),".tiff", sep=""), width=1200,  height=1200, units="px", pointsize = 24)
	  #jpeg(paste("Abn_Cya_sum34_max4_N_Cyanolichen_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),".jpeg", sep=""))
	  #jpeg(paste("Abn_Cya_sum34_max4_N_Cyanolichen_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.jpeg", sep=""))
    tiff(paste("Abn_Cya_sum34_max4_N_Cyanolichen_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.tiff", sep=""), width=1200,  height=1200, units="px", pointsize = 24)
    
    x<-Abn_Cya_sum34_max4_N_Mod_pred[2]
    x<-unlist(x, recursive = F, use.names =T)
    y<-Abn_Cya_sum34_max4_N_Mod_pred[4]
	  clim<-unlist(y[1], recursive = F, use.names =T)
	  clim_pred<-unlist(clim[1], recursive = F, use.names =T)
    fit<-unlist(x[1], recursive = F, use.names =T)
    low<-unlist(x[2], recursive = F, use.names =T)
    high<-unlist(x[3], recursive = F, use.names =T)
    SOURCE_N_fPlot<-cbind(SOURCE, fit, low, high, clim_pred)
	  SOURCE_N_fPlot<-subset(SOURCE_N_fPlot, N < 12.5)    
	  palette(c("darkorange","darkgreen"))
	  plot(Abn_Cya_sum34_max4~N, data=SOURCE, xlab="CMAQ N kg/ha/yr", ylab="", main="90th Quantile Regression of \n Cyanolichen Abundance vs N deposition", col=Area, cex.main = 1.6, cex.lab=1.7, cex.axis = 1.7)
    title(ylab ="Cyanolichen Abundance", cex.lab=1.7, line=2.5)
    points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$clim_pred, col=alpha('grey',0.25), pch= 19)
    points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$fit, col='blue', pch= 19)
    points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$low, col='black', pch= 4, cex=0.25) 
    points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$high, col='black', pch= 4, cex=0.25)	
    points(0.08,13.49563,col='red', pch=20,cex=2)
		#points(0.4,12.82085,col='red', pch=20,cex=2)
		points(0.7 ,12.14607,col='red', pch=20,cex=2)
		points(1.3 ,10.7965 ,col='red', pch=20,cex=2)
		points(3.5 ,6.747815,col='red', pch=20,cex=2)
		points(6.6 ,2.699126,col='red', pch=20,cex=2)
    #text(max(SOURCE$N)*0.55,max(SOURCE$Abn_Cya_sum34_max4)*0.6,labels=paste("R1=",print(round(Abn_Cya_sum34_max4_N_Mods_fit[2],2))))
    ## Old way with formula pasted but without coefficients and model fit stats, which were moved to a table
    ## text(max(SOURCE$N)*0.75,max(SOURCE$Abn_Cya_sum34_max4)*0.75,labels=paste(as.character(Abn_Cya_sum34_max4_N_Formula[2])), cex=1)
    ## text(max(SOURCE$N)*0.75,max(SOURCE$Abn_Cya_sum34_max4)*0.7,labels=paste("AIC=",print(round(as.numeric(Abn_Cya_sum34_max4_N_Mods_AIC[2],2)))))
	  ## text(max(SOURCE$N)*0.75,max(SOURCE$Abn_Cya_sum34_max4)*0.85,labels=paste("90th quantile regression model"))
  
    ## New way pasting formula directly
	     ##text(max(SOURCE$N)*0.55,max(SOURCE$Abn_Cya_sum34_max4)*0.75,labels=paste("Cyanolichen Abundance = 5.87254-(246.66380*N)+(107.45163*N2)", cex=1))
    
		
		#text(max(SOURCE$N)*0.5,max(SOURCE$Abn_Cya_sum34_max4)*0.6,labels=bquote(atop(.("Cyanolichen Abundance = 5.87-"),.("246.66*N+107.45*")*N^2)), cex=1.5)
		text(max(SOURCE$N)*0.5,max(SOURCE$Abn_Cya_sum34_max4)*0.9,labels=bquote(atop(.("Cyanolichen Abundance = 5.87-"),.("246.66*N+107.45*")*N^2)), cex=1.5)
		#legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1.47)
		#legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1)
		
    dev.off()


###### Make final version of figure with grey cloud of points in the background of the climate fitted values, the raw data as open black circles
###### and the fitted polynomial line with the predicted fitted 

				#five_pct_Abn_Cya_sum34_max4  <- max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*0.05) %>% as.data.frame()
				#ten_pct_Abn_Cya_sum34_max4   <- max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*0.1) %>% as.data.frame()
				#twenty_pct_Abn_Cya_sum34_max4<- max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*0.2) %>% as.data.frame()
				#fifty_pct_Abn_Cya_sum34_max4 <- max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*0.5) %>% as.data.frame()
				#eighty_pct_Abn_Cya_sum34_max4 <- max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*0.8) %>% as.data.frame()
				#
				#colnames(five_pct_Abn_Cya_sum34_max4)<-"N"
				#colnames(ten_pct_Abn_Cya_sum34_max4)<-"N"   
				#colnames(twenty_pct_Abn_Cya_sum34_max4)<-"N"
				#colnames(fifty_pct_Abn_Cya_sum34_max4)<-"N"
				#colnames(eighty_pct_Abn_Cya_sum34_max4)<-"N"

	
				#> five_pct_Abn_Cya_sum34_max4  
				#1 12.82085
				#> #
				#> ten_pct_Abn_Cya_sum34_max4   
				#1 12.14607
				#> #
				#> twenty_pct_Abn_Cya_sum34_max4
				#1 10.7965
				#> #
				#> fifty_pct_Abn_Cya_sum34_max4 
				#1 6.747815
				#> #
				#> eighty_pct_Abn_Cya_sum34_max4
				#1 2.699126
		
				## x - 5
				## y - Abn_Cya_sum34_max4_N_Mods_run
				## ModRun - Abn_Cya_sum34_max4_N_poly
				## ModPred Abn_Cya_sum34_max4_N_Mod_pred
			
	#GetYforX<-function (x,y,ModRun,ModPred)
	#{		
	#test_n<-as.data.frame(x)
	#colnames(test_n)<-"N"
	#(predict(y$ModRun, test_n))/max(ModPred$ModRun)
	#}
	#GetYforX(5,Abn_Cya_sum34_max4_N_Mods_run,Abn_Cya_sum34_max4_N_poly,Abn_Cya_sum34_max4_N_Mod_pred)
			
			test_n<-as.data.frame(7.5)
			colnames(test_n)<-"N"
			1-(predict(Abn_Cya_sum34_max4_N_Mods_run$Abn_Cya_sum34_max4_N_poly, test_n))/max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)
		
		## Heuristically solved for X given a Y using the polynomial 
			## Heuristically solved for X given a Y ...	34.58681
				## 0 % max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit) 
				##					  Heuristically solved for X given a Y ...13.49563 Abn_Cya_sum34_max4 0.08 N kg/ha/yr
				## 5 % loss of N. Cya Heuristically solved for X given a Y ...12.82085 Abn_Cya_sum34_max4 0.4  N kg/ha/yr
				## 10% loss of N. Cya Heuristically solved for X given a Y ...12.14607 Abn_Cya_sum34_max4 0.7  N kg/ha/yr
				## 20% loss of N. Cya Heuristically solved for X given a Y ...10.7965  Abn_Cya_sum34_max4 1.3  N kg/ha/yr
				## 50% loss of N. Cya Heuristically solved for X given a Y ...6.747815 Abn_Cya_sum34_max4 3.5  N kg/ha/yr
				## 80% loss of N. Cya Heuristically solved for X given a Y ...2.699126 Abn_Cya_sum34_max4 6.6  N kg/ha/yr

				new_n_1kg	<-1
				new_n_1.5kg	<-1.5
				new_n_2kg	<-2
				new_n_2.5kg	<-2.5
				new_n_3kg	<-3
				new_n_5kg   <-5
				new_n_1kg	<-as.data.frame(new_n_1kg)
				new_n_1.5kg	<-as.data.frame(new_n_1.5kg)
				new_n_2kg	<-as.data.frame(new_n_2kg)
				new_n_2.5kg	<-as.data.frame(new_n_2.5kg)
				new_n_3kg	<-as.data.frame(new_n_3kg)
				new_n_5kg   <-as.data.frame(new_n_5kg)
				
				colnames(new_n_1kg)<- "N"
				colnames(new_n_1.5kg)<- "N"
				colnames(new_n_2kg)<- "N"
				colnames(new_n_2.5kg)<- "N"
				colnames(new_n_3kg)<- "N"
				colnames(new_n_5kg)<- "N"

				c(
				 1-(predict(Abn_Cya_sum34_max4_N_Mods_run$Abn_Cya_sum34_max4_N_poly, new_n_1kg))/max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)
				,1-(predict(Abn_Cya_sum34_max4_N_Mods_run$Abn_Cya_sum34_max4_N_poly,new_n_1.5kg))/max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)
				,1-(predict(Abn_Cya_sum34_max4_N_Mods_run$Abn_Cya_sum34_max4_N_poly,new_n_2kg))/max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)
				,1-(predict(Abn_Cya_sum34_max4_N_Mods_run$Abn_Cya_sum34_max4_N_poly,new_n_2.5kg))/max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)
				,1-(predict(Abn_Cya_sum34_max4_N_Mods_run$Abn_Cya_sum34_max4_N_poly,new_n_3kg))/max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)
				,1-(predict(Abn_Cya_sum34_max4_N_Mods_run$Abn_Cya_sum34_max4_N_poly,new_n_5kg))/max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)
				)
				
				c(
				max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit) ## 0% change
				,max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*0.05) ## 5% change
				,max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*0.10) ## 10% change
				,max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*0.20) ## 20% change
				,max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*0.50) ## 50% change
				,max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*0.80) ## 80% change
				)
				
	################################################################################################################
#Abundance of cyanolichens vs sulphur deposition
################################################################################################################

	## Input files are generated in Calculate Abundances script. I have to switch each of these sources depending on which measure of Abundance wanted.
	## We decided to use the maximum Abundance because it was the most straightforward to explain with little to no loss in model fit.
## Instead of changing SOURCE, I should build a separate pipeline for each measure of Abundance (summing 3s and 4s, summing log transformed values and only max value)
	SOURCE<-LichDb_sFncGrpSens_All_Abun_S
	##SOURCE<-subset(LichDb_sFncGrpSens_All_Abun_S, S<20)
XVAR<-c("S"                                                               
 ,"poly(S, 2)"                                                      
 ,"S+maxaug_c+mindec_c+precip_cm+CMD"         
 ,"poly(S, 2)+maxaug_c+mindec_c+precip_cm+CMD"
 ,"maxaug_c+mindec_c+precip_cm+CMD"
 ,"S+maxaug_c"                                                      
 ,"S+mindec_c"                                                      
 ,"S+precip_cm"                                                     
 ,"S+continen"
,"S+CMD"
,"maxaug_c"
,"mindec_c"
,"precip_cm"
,"continen"
,"CMD") 

XVAR_names<-c("S","S_poly","S_clim","S_poly_clim","clim","S_maxaug_c","S_mindec_c","S_precip_cm","S_continen","S_CMD", "Max C Aug Temp","Min C Dec Temp","Precip_cm","Continentality","CMD");

Abn_Cya_S_sum34_max4_Mods<-paste("Abn_Cya_sum34_max4 ~ ",paste(XVAR),sep="")

Abn_Cya_S_sum34_max4_Formula<-lapply(Abn_Cya_S_sum34_max4_Mods,as.formula)

Abn_Cya_S_sum34_max4_Mods_run<-lapply(Abn_Cya_S_sum34_max4_Formula,function (x) rq(x, tau=0.9, data=SOURCE))
	## Make null models
Abn_Cya_S_sum34_max4_S_intercept<-rq(formula = Abn_Cya_sum34_max4 ~ 1, tau = 0.9, data = SOURCE, model = T)
	## Grab each model fit but the intercept is hard coded 

Abn_Cya_S_sum34_max4_Mods_fit<-mapply(function (y) 1-Abn_Cya_S_sum34_max4_Mods_run[[y]]$rho/Abn_Cya_S_sum34_max4_S_intercept$rho, 1:length(XVAR))
	## makes a list of mod names

Abn_Cya_S_sum34_max4_Mod_names<-lapply(XVAR_names, function (x) paste("Abn_Cya_sum34_max4_",paste(x),sep=""))

names(Abn_Cya_S_sum34_max4_Mods_run)<-Abn_Cya_S_sum34_max4_Mod_names

	## Save model stats as a spreadsheet of manuscript

Abn_Cya_S_sum34_max4_Mods_AIC<-lapply(Abn_Cya_S_sum34_max4_Mods_run,AIC)

Abn_Cya_sum34_max4_Lichen_Abundance_vs_S_Stats<-cbind(XVAR_names,Abn_Cya_S_sum34_max4_Mods_AIC,Abn_Cya_S_sum34_max4_Mods_fit,paste(Abn_Cya_S_sum34_max4_Formula))
write.csv(Abn_Cya_sum34_max4_Lichen_Abundance_vs_S_Stats,"Abundance_Abn_Cya_sum34_max4_S_all_models.csv")

#Reduc_Mods_forTable_S<-c("S","S_poly","S_poly_clim")
#Abn_Cya_sum34_max4_LicheS_Abundance_vs_S_Stats_reduc<-cbind(XVAR_names,Abn_Cya_S_sum34_max4_Mods_fit,Abn_Cya_S_sum34_max4_Mods_AIC,paste(Abn_Cya_S_sum34_max4_Formula))
#Abn_Cya_sum34_max4_LicheS_Abundance_vs_S_Stats_reduc<-subset(Abn_Cya_sum34_max4_LicheS_Abundance_vs_S_Stats_reduc, XVAR_names %in% Reduc_Mods_forTable_S)
#write.csv(Abn_Cya_sum34_max4_Lichen_Abundance_vs_S_Stats_reduc,"Abundance_Abn_Cya_sum34_max4_S_all_models_reduc.csv")

	## Make predicted data files
	## Use a model file ... Abn_Cya_S_sum34_max4_Mods_run as input
	## Make a SOURCE where there aren't wild estimates at higher S values based on a visual estimate of the 
MyFunc2<-function (x) {as.data.frame(predict(x, SOURCE, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000))}

	##This produces a matrix with a list of 10 objects with 3 observations, each with 5322 entries and 10 columns. 
Abn_Cya_S_sum34_max4_Mod_pred<-mapply(MyFunc2,Abn_Cya_S_sum34_max4_Mods_run)%>%as.data.frame()
  
  #tiff(paste("Abn_Cya_sum34_max4_S_Cyanolichen_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),".tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
	jpeg(paste("Abn_Cya_sum34_max4_S_Cyanolichen_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),".jpeg",sep=""))
	#jpeg(paste("Abn_Cya_sum34_max4_S_Cyanolichen_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.jpeg",sep=""))
  #tiff(paste("Abn_Cya_sum34_max4_S_Cyanolichen_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
  
	## This works but is hard coded for each predictor file. Have to change working directory each time
	## Make final figure with grey dots in the background
	## Need to add red dots for percentage change
	x<-Abn_Cya_S_sum34_max4_Mod_pred[2]
	y<-Abn_Cya_S_sum34_max4_Mod_pred[4]
	clim<-unlist(y[1], recursive = F, use.names =T)
	x<-unlist(x, recursive = F, use.names =T)
	fit<-unlist(x[1], recursive = F, use.names =T)
	low<-unlist(x[2], recursive = F, use.names =T)
	high<-unlist(x[3], recursive = F, use.names =T)
	clim_pred<-unlist(clim[1], recursive = F, use.names =T)
	SOURCE_S_fPlot<-cbind(SOURCE, fit, low, high, clim_pred)
	SOURCE_S_fPlot<-subset(SOURCE_S_fPlot, S < 20)
	palette(c("darkorange","darkgreen"))
	plot(Abn_Cya_sum34_max4~S, data=SOURCE, xlab="CMAQ S kg/ha/yr", ylab="", main="90th Quantile Regression of \n Cyanolichen Abundance vs S deposition", col=Area, cex.main = 1.6, cex.lab=1.7, cex.axis = 1.7)
	title(ylab ="Cyanolichen Abundance", cex.lab=1.7, line=2.5)
	points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$clim_pred, col=alpha('grey',0.25), pch= 19)
	points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$fit, col='blue', pch= 19)
	points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$low, col='black', pch= 4, cex=0.25) 
	points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$high, col='black', pch= 4, cex=0.25)	
		  points(0.21, 9.733974, col='red', pch=20,cex=2)
		  #points(0.7, 9.247276 ,col='red', pch=20,cex=2)
		  points(1.2 , 8.760577 ,col='red', pch=20,cex=2)
		  points(2.3 , 7.78718  ,col='red', pch=20,cex=2)
		  points(5.9 , 4.866987 ,col='red', pch=20,cex=2)
		  points(10.9, 1.946795 ,col='red', pch=20,cex=2)
    	#text(max(SOURCE$S)*0.75,max(SOURCE$Abn_Cya_sum34_max4)*0.65,labels=paste("R1=",print(round(Abn_Cya_S_sum34_max4_Mods_fit[2],2))))
    	#text(max(SOURCE$S)*0.75,max(SOURCE$Abn_Cya_sum34_max4)*0.7,labels=paste("AIC=",print(round(as.numeric(Abn_Cya_S_sum34_max4_Mods_AIC[2],2)))))
    	#text(max(SOURCE$S)*0.75,max(SOURCE$Abn_Cya_sum34_max4)*0.75,labels=paste(as.character(Abn_Cya_S_sum34_max4_Formula[2])), cex=1)
	
	text(max(SOURCE$S)*0.5,max(SOURCE$Abn_Cya_sum34_max4)*0.6,labels=bquote(atop(.("Cyanolichen Abundance = .08 -"), .("204.74*S + 94.51*")*S^2)), cex=1.5)
	#text(max(SOURCE$S)*0.5,max(SOURCE$Abn_Cya_sum34_max4)*0.9,labels=bquote(atop(.("Cyanolichen Abundance = .08 -"), .("204.74*S + 94.51*")*S^2)), cex=1.5)
	#legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1.47)
	legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1)
	dev.off()
##pdf(paste("Abn_Cya_sum34_max4_Lichen_Abundance_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),".pdf",sep=""))

#MyPlot<-function (x)
#    
#{
#    x<-Abn_Cya_S_sum34_max4_Mod_pred[i]
#    x<-unlist(x, recursive = F, use.names =T)
#    fit<-unlist(x[1], recursive = F, use.names =T)
#    low<-unlist(x[2], recursive = F, use.names =T)
#    high<-unlist(x[3], recursive = F, use.names =T)
#    SOURCE_S_fPlot<-cbind(SOURCE, fit, low, high)
#    SOURCE_S_fPlot<-subset(SOURCE_S_fPlot, S < 20)
#    plot(Abn_Cya_sum34_max4~S, data=SOURCE, xlab="CMAQ S kg/ha/yr", ylab="Abundance Abn_Cya_sum34_max4 Lichens", main="Abundance of  Abn_Cya_sum34_max4 lichen species vs S deposition")
#    points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$fit, col='blue', pch= 19)
#    points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$low, col='aquamarine4', pch= 4, cex=0.25) 
#    points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$high, col='aquamarine4', pch= 4, cex=0.25)	
#    text(max(SOURCE$S)*0.4,max(SOURCE$Abn_Cya_sum34_max4)*0.8,labels=paste("R1=",print(round(Abn_Cya_S_sum34_max4_Mods_fit[i],2))))
#    text(max(SOURCE$S)*0.4,max(SOURCE$Abn_Cya_sum34_max4)*0.75,labels=paste(as.character(Abn_Cya_S_sum34_max4_Formula[i])), cex=0.7)
#    legend('topright', legend=c("Raw Data","Fitted Values", "Prediction Interval"), col=c("black","blue", "aquamarine4"), pch=16)
# 
#}	
#
#for (i in 1:length(XVAR)) MyPlot()
#
#dev.off()

##### Make final version of figure with grey cloud of points in the background of the climate fitted values, the raw data as open black circles
###### and the fitted polynomial line with the predicted fitted 
	
	## The used the original model prediction but the range	of predictions threw off using this way since, at the high end
	## predictions using S went up (so max Cyan Abundance was at max S dep). Had to fix this by using an object from the plotting section
	## for cyanolichens stopping predictions at 20 kg S ha yr.
	
	#five_pct_Abn_Cya_sum34_max4_N  <- max(Abn_Cya_S_sum34_max4_Mod_pred$Abn_Cya_sum34_max4_S_poly$fit)-(max(Abn_Cya_S_sum34_max4_Mod_pred$Abn_Cya_sum34_max4_S_poly$fit)*0.05) %>% as.data.frame()
	#ten_pct_Abn_Cya_sum34_max4_N   <- max(Abn_Cya_S_sum34_max4_Mod_pred$Abn_Cya_sum34_max4_S_poly$fit)-(max(Abn_Cya_S_sum34_max4_Mod_pred$Abn_Cya_sum34_max4_S_poly$fit)*0.1) %>% as.data.frame()
	#twenty_pct_Abn_Cya_sum34_max4_N<- max(Abn_Cya_S_sum34_max4_Mod_pred$Abn_Cya_sum34_max4_S_poly$fit)-(max(Abn_Cya_S_sum34_max4_Mod_pred$Abn_Cya_sum34_max4_S_poly$fit)*0.2) %>% as.data.frame()
	#fifty_pct_Abn_Cya_sum34_max4_N <- max(Abn_Cya_S_sum34_max4_Mod_pred$Abn_Cya_sum34_max4_S_poly$fit)-(max(Abn_Cya_S_sum34_max4_Mod_pred$Abn_Cya_sum34_max4_S_poly$fit)*0.5) %>% as.data.frame()
	#eighty_pct_Abn_Cya_sum34_max4_N <- max(Abn_Cya_S_sum34_max4_Mod_pred$Abn_Cya_sum34_max4_S_poly$fit)-(max(Abn_Cya_S_sum34_max4_Mod_pred$Abn_Cya_sum34_max4_S_poly$fit)*0.8) %>% as.data.frame()
	
	five_pct_Abn_Cya_sum34_max4_S  <- max(SOURCE_S_fPlot$fit)-(max(SOURCE_S_fPlot$fit)*0.05) %>% as.data.frame()
	ten_pct_Abn_Cya_sum34_max4_S   <- max(SOURCE_S_fPlot$fit)-(max(SOURCE_S_fPlot$fit)*0.1) %>% as.data.frame()
	twenty_pct_Abn_Cya_sum34_max4_S<- max(SOURCE_S_fPlot$fit)-(max(SOURCE_S_fPlot$fit)*0.2) %>% as.data.frame()
	fifty_pct_Abn_Cya_sum34_max4_S <- max(SOURCE_S_fPlot$fit)-(max(SOURCE_S_fPlot$fit)*0.5) %>% as.data.frame()
	eighty_pct_Abn_Cya_sum34_max4_S<- max(SOURCE_S_fPlot$fit)-(max(SOURCE_S_fPlot$fit)*0.8) %>% as.data.frame()
	
	colnames(five_pct_Abn_Cya_sum34_max4_S)<-"S"
	colnames(ten_pct_Abn_Cya_sum34_max4_S)<-"S"   
	colnames(twenty_pct_Abn_Cya_sum34_max4_S)<-"S"
	colnames(fifty_pct_Abn_Cya_sum34_max4_S)<-"S"
	colnames(eighty_pct_Abn_Cya_sum34_max4_S)<-"S"

	## Find the value of Abn_Cya_sum34_max4 at zero deposition
		## Would it be zero deposition
			c(
			five_pct_Abn_Cya_sum34_max4_S  
			,ten_pct_Abn_Cya_sum34_max4_S   
			,twenty_pct_Abn_Cya_sum34_max4_S
			,fifty_pct_Abn_Cya_sum34_max4_S 
			,eighty_pct_Abn_Cya_sum34_max4_S
			)
			
			## x - 5
			## y - Abn_Cya_S_sum34_max4_Mods_run
			## ModRun - Abn_Cya_sum34_max4_S_poly
			## ModPred Abn_Cya_S_sum34_max4_Mod_pred
			
			    
	#Started to automate the function of finding a given X for a value of Y but ran into issues referencing elements using 
	# the characters "$". Do I need to unlist list first?
			
	#GetYforX<-function (x,y,ModRun,ModPred)
	#{		
	#test_s<-as.data.frame(x)
	#colnames(test_s)<-"S"
	#(predict(ModRun, test_s))/max(ModPred)
	#}
	#
	#GetYforX(3,Abn_Cya_S_sum34_max4_Mods_run,Abn_Cya_sum34_max4_S_poly,Abn_Cya_S_sum34_max4_Mod_pred)
	
	
		## Heuristically solved for X given a Y using the polynomial 
		## NOTE that the denominator needs to be changed for reponses wherew max Y value is at max X value due to poor fit 
		## at higher X values. Need to change denomnitor to reduced SOURCE used for plotting.
			test_s<-as.data.frame(0.5)
			colnames(test_s)<-"S"
			1-(predict(Abn_Cya_S_sum34_max4_Mods_run$Abn_Cya_sum34_max4_S_poly, test_s))/max(SOURCE_S_fPlot$fit)

				## 0 % (max Abn_Cya_sum34_max4 ) 
									# Heuristically solved for X S kg/ha/yr given a Y  9.733974 ... 0.21 S kg/ha/yr
				## 5 % loss of N. Cya Heuristically solved for X S kg/ha/yr given a Y  9.247276 ... 0.7  S kg/ha/yr
				## 10% loss of N. Cya Heuristically solved for X S kg/ha/yr given a Y  8.760577 ... 1.2  S kg/ha/yr
				## 20% loss of N. Cya Heuristically solved for X S kg/ha/yr given a Y  7.78718  ... 2.3  S kg/ha/yr
				## 50% loss of N. Cya Heuristically solved for X S kg/ha/yr given a Y  4.866987 ... 5.9  S kg/ha/yr
				## 80% loss of N. Cya Heuristically solved for X S kg/ha/yr given a Y  1.946795 ... 10.9 S kg/ha/yr
  
				new_S_1kg	<-1
				new_S_1.5kg	<-1.5
				new_S_2kg	<-2
				new_S_2.5kg	<-2.5
				new_S_3kg	<-3
				new_S_5kg   <-5
				new_S_1kg	<-as.data.frame(new_S_1kg)
				new_S_1.5kg	<-as.data.frame(new_S_1.5kg)
				new_S_2kg	<-as.data.frame(new_S_2kg)
				new_S_2.5kg	<-as.data.frame(new_S_2.5kg)
				new_S_3kg	<-as.data.frame(new_S_3kg)
				new_S_5kg   <-as.data.frame(new_S_5kg)
				
				colnames(new_S_1kg)<- "S"
				colnames(new_S_1.5kg)<- "S"
				colnames(new_S_2kg)<- "S"
				colnames(new_S_2.5kg)<- "S"
				colnames(new_S_3kg)<- "S"
				colnames(new_S_5kg)<- "S"

				c(
				 1-(predict(Abn_Cya_S_sum34_max4_Mods_run$Abn_Cya_sum34_max4_S_poly, new_S_1kg))/max(SOURCE_S_fPlot$fit)
				,1-(predict(Abn_Cya_S_sum34_max4_Mods_run$Abn_Cya_sum34_max4_S_poly,new_S_1.5kg))/max(SOURCE_S_fPlot$fit)
				,1-(predict(Abn_Cya_S_sum34_max4_Mods_run$Abn_Cya_sum34_max4_S_poly,new_S_2kg))/max(SOURCE_S_fPlot$fit)
				,1-(predict(Abn_Cya_S_sum34_max4_Mods_run$Abn_Cya_sum34_max4_S_poly,new_S_2.5kg))/max(SOURCE_S_fPlot$fit)
				,1-(predict(Abn_Cya_S_sum34_max4_Mods_run$Abn_Cya_sum34_max4_S_poly,new_S_3kg))/max(SOURCE_S_fPlot$fit)
				,1-(predict(Abn_Cya_S_sum34_max4_Mods_run$Abn_Cya_sum34_max4_S_poly,new_S_5kg))/max(SOURCE_S_fPlot$fit)
				)
				
				c(
				max(SOURCE_S_fPlot$fit) ## 0% change
				,max(SOURCE_S_fPlot$fit)-(max(SOURCE_S_fPlot$fit)*0.05) ## 5% change
				,max(SOURCE_S_fPlot$fit)-(max(SOURCE_S_fPlot$fit)*0.10) ## 10% change
				,max(SOURCE_S_fPlot$fit)-(max(SOURCE_S_fPlot$fit)*0.20) ## 20% change
				,max(SOURCE_S_fPlot$fit)-(max(SOURCE_S_fPlot$fit)*0.50) ## 50% change
				,max(SOURCE_S_fPlot$fit)-(max(SOURCE_S_fPlot$fit)*0.80) ## 80% change
				)