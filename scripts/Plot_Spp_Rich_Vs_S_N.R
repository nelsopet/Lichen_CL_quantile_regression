## Load required packages
require(plyr)
require(quantreg)
require(tidyverse)


############################################################################################################################################
#Lichen species richness vs Nitrogen 
################################################################################################################

	## Input files are generated in Calculate Counts script. I have to switch each of these sources depending on which measure of Count wanted.
	## We decided to use the maximum Count because it was the most straightforward to explain with little to no loss in model fit.
## Instead of changing SOURCE, I should build a separate pipeline for each measure of Count (summing 3s and 4s, summing log transformed values and only max value)
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

spp_rich_N_Mods<-paste("spp_rich ~ ",paste(XVAR),sep="")

spp_rich_N_Formula<-lapply(spp_rich_N_Mods,as.formula)

spp_rich_N_Mods_run<-lapply(spp_rich_N_Formula,function (x) rq(x, tau=0.9, data=SOURCE))


	## Make null models
spp_rich_N_intercept<-rq(formula = spp_rich ~ 1, tau = 0.9, data = SOURCE, model = T)


	## Grab each model fit but the intercept is hard coded 

spp_rich_N_Mods_fit<-mapply(function (y) 1-spp_rich_N_Mods_run[[y]]$rho/spp_rich_N_intercept$rho, 1:length(XVAR))


	## makes a list of mod names

spp_rich_N_Mod_names<-lapply(XVAR_names, function (x) paste("spp_rich_",paste(x),sep=""))

names(spp_rich_N_Mods_run)<-spp_rich_N_Mod_names


	## Save model stats as a spreadsheet of manuscript

spp_rich_N_Mods_AIC<-lapply(spp_rich_N_Mods_run,AIC)

spp_rich_N_Lichen_Count_vs_N_Stats<-cbind(XVAR_names,spp_rich_N_Mods_AIC,spp_rich_N_Mods_fit,paste(spp_rich_N_Formula))
write.csv(spp_rich_N_Lichen_Count_vs_N_Stats,"Count_spp_rich_N_all_models.csv")


	## Make predicted data files
	## Use a model file ... spp_rich_Mods_run as input

MyFunc2<-function (x) {as.data.frame(predict(x, SOURCE, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000))}

	##This produces a matrix with a list of 10 objects with 3 observations, each with 5322 entries and 10 columns. 
spp_rich_N_Mod_pred<-mapply(MyFunc2,spp_rich_N_Mods_run)%>%as.data.frame()

	## This works but is hard coded for each predictor file. Have to change working directory each time

pdf(paste("spp_rich_N_Lichen_Count_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),".pdf",sep=""))

MyPlot<-function (x)
    
{
    x<-spp_rich_N_Mod_pred[i]
    x<-unlist(x, recursive = F, use.names =T)
    fit<-unlist(x[1], recursive = F, use.names =T)
    low<-unlist(x[2], recursive = F, use.names =T)
    high<-unlist(x[3], recursive = F, use.names =T)
    plot(spp_rich~N, data=SOURCE, xlab="CMAQ N kg/ha/yr", ylab="Lichen species richness", main="Lichen species richness vs N deposition")
    points(SOURCE$N, fit, col='blue', pch= 19)
    points(SOURCE$N, low, col='aquamarine4', pch= 4, cex=0.25) 
    points(SOURCE$N, high, col='aquamarine4', pch= 4, cex=0.25)	
		#points(0.08,5.63, col='red', pch=20)
		#points(0.44,5.35, col='red', pch=20)
		#points(0.80,5.07, col='red', pch=20)
		#points(1.58,4.51, col='red', pch=20)
		#points(4.27,2.82, col='red', pch=20)
		#points(7.95,1.13, col='red', pch=20)
    text(max(SOURCE$N)*0.4,max(SOURCE$spp_rich)*0.8,labels=paste("R1=",print(round(spp_rich_N_Mods_fit[i],2))))
    text(max(SOURCE$N)*0.4,max(SOURCE$spp_rich)*0.75,labels=paste(as.character(spp_rich_N_Formula[i])), cex=0.7)
    legend('topright', legend=c("Raw Data","Fitted Values", "Prediction Interval"), col=c("black","blue", "aquamarine4"), pch=16)
 
}	

for (i in 1:length(XVAR)) MyPlot()

dev.off()
  #tiff(paste("spp_rich_N_Lichen_Count_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),".tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
	jpeg(paste("spp_rich_N_Lichen_Count_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),".jpeg",sep=""))
	#jpeg(paste("spp_rich_N_Lichen_Count_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.jpeg",sep=""))
  #tiff(paste("spp_rich_N_Lichen_Count_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
  
    x<-spp_rich_N_Mod_pred[2]
    x<-unlist(x, recursive = F, use.names =T)
	y<-spp_rich_N_Mod_pred[4]
	clim<-unlist(y[1], recursive = F, use.names =T)
	clim_pred<-unlist(clim[1], recursive = F, use.names =T)

    fit<-unlist(x[1], recursive = F, use.names =T)
    low<-unlist(x[2], recursive = F, use.names =T)
    high<-unlist(x[3], recursive = F, use.names =T)
    SOURCE_N_fPlot<-cbind(SOURCE, fit, low, high, clim_pred)
	  SOURCE_N_fPlot<-subset(SOURCE_N_fPlot, N < 12.5)    
    palette(c("darkorange","darkgreen"))
	  plot(spp_rich~N, data=SOURCE, xlab="CMAQ N kg/ha/yr", ylab="", main="90th Quantile Regresssion of Lichen \n Species Richness vs N deposition", col=Area, cex.main = 1.6, cex.lab=1.7, cex.axis = 1.7)
    title(ylab ="Lichen Species Richnesss", cex.lab=1.7, line=2.5)
    points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$clim_pred, col=alpha('grey',0.25), pch= 19)
    points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$fit, col='blue', pch= 19)
    points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$low, col='black', pch= 4, cex=0.25) 
    points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$high, col='black', pch= 4, cex=0.25)	
	points(0.08,32.80, col='red',     pch=20,cex=2)
	#points(0.86,31.16,col='red',     pch=20,cex=2)
	points(1.70,29.52, col='red',     pch=20,cex=2)
	points(3.50,26.24, col='red',     pch=20,cex=2)
	#points(6.65,5.203122, col='red', pch=20,cex=2)
	#points(12.8,2.081249, col='red', pch=20,cex=2)
     #text(max(SOURCE$N)*0.8,max(SOURCE$spp_rich)*0.70,labels=paste("R1=",print(round(spp_rich_N_Mods_fit[2],2))))
     #text(max(SOURCE$N)*0.8,max(SOURCE$spp_rich)*0.75,labels=paste(as.character(spp_rich_N_Formula[2])), cex=1)
     #text(max(SOURCE$N)*0.8,max(SOURCE$spp_rich)*0.65,labels=paste("AIC=",print(round(as.numeric(spp_rich_N_Mods_AIC[2],2)))))
	
	text(max(SOURCE$N)*0.6,max(SOURCE$spp_rich)*0.65,labels=bquote(atop(.("Lichen Species richness = 24.68-"),.("313.26*N + 81.44*")*N^2)), cex=1.5)
	#text(max(SOURCE$N)*0.6,max(SOURCE$spp_rich)*0.9,labels=bquote(atop(.("Lichen Species richness = 24.68-"),.("313.26*N + 81.44*")*N^2)), cex=1.5)
	#legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1.47)
	legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1)
	dev.off()

###### Make final version of figure with grey cloud of points in the background of the climate fitted values, the raw data as open black circles
###### and the fitted polynomial line with the predicted fitted 
		
		## x - 5
		## y - spp_rich_N_Mods_run
		## ModRun - spp_rich_N_poly
		## ModPred spp_rich_N_Mod_pred
			
	#GetYforX<-function (x,y,ModRun,ModPred)
	#{		
	#test_n<-as.data.frame(x)
	#colnames(test_n)<-"N"
	#(predict(y$ModRun, test_n))/max(ModPred$ModRun)
	#}
	#GetYforX(5,spp_rich_N_Mods_run,spp_rich_N_poly,spp_rich_N_Mod_pred)
			
			test_n<-as.data.frame(7.5)
			colnames(test_n)<-"N"
			1-(predict(spp_rich_N_Mods_run$spp_rich_N_poly, test_n))/max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)
		
		## Heuristically solved for X given a Y using the polynomial 
			## Heuristically solved for X given a Y ...	34.58681
				## 0 % (max spp_rich ) Heuristically solved for X given a Y   10.41...    0.08  N kg/ha/yr
				## 5 % loss of N. Cya Heuristically solved for X given a Y 9.89...	   0.6 N kg/ha/yr
				## 10% loss of N. Cya Heuristically solved for X given a Y 9.37...	   1.25 N kg/ha/yr
				## 20% loss of N. Cya Heuristically solved for X given a Y 8.324995... 2.45 N kg/ha/yr
				## 50% loss of N. Cya Heuristically solved for X given a Y 5.203122... 6.65 N kg/ha/yr
				## 80% loss of N. Cya Heuristically solved for X given a Y 2.081249... 12.8 N kg/ha/yr

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
				 1-(predict(spp_rich_N_Mods_run$spp_rich_N_poly, new_n_1kg))/max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)
				,1-(predict(spp_rich_N_Mods_run$spp_rich_N_poly,new_n_1.5kg))/max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)
				,1-(predict(spp_rich_N_Mods_run$spp_rich_N_poly,new_n_2kg))/max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)
				,1-(predict(spp_rich_N_Mods_run$spp_rich_N_poly,new_n_2.5kg))/max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)
				,1-(predict(spp_rich_N_Mods_run$spp_rich_N_poly,new_n_3kg))/max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)
				,1-(predict(spp_rich_N_Mods_run$spp_rich_N_poly,new_n_5kg))/max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)
				)
				
				c(
				max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit) ## 0% change
				,max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)-(max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)*0.05) ## 5% change
				,max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)-(max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)*0.10) ## 10% change
				,max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)-(max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)*0.20) ## 20% change
				,max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)-(max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)*0.50) ## 50% change
				,max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)-(max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)*0.80) ## 80% change
				)
								
################################################################################################################
#Lichen species richnesss vs sulphur deposition
################################################################################################################

	## Input files are generated in Calculate Counts script. I have to switch each of these sources depending on which measure of Count wanted.
	## We decided to use the maximum Count because it was the most straightforward to explain with little to no loss in model fit.
## Instead of changing SOURCE, I should build a separate pipeline for each measure of Count (summing 3s and 4s, summing log transformed values and only max value)
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

spp_rich_Mods<-paste("spp_rich ~ ",paste(XVAR),sep="")

spp_rich_Formula<-lapply(spp_rich_Mods,as.formula)

spp_rich_Mods_run<-lapply(spp_rich_Formula,function (x) rq(x, tau=0.9, data=SOURCE))
	## Make null models
spp_rich_S_intercept<-rq(formula = spp_rich ~ 1, tau = 0.9, data = SOURCE, model = T)
	## Grab each model fit but the intercept is hard coded 

spp_rich_Mods_fit<-mapply(function (y) 1-spp_rich_Mods_run[[y]]$rho/spp_rich_S_intercept$rho, 1:length(XVAR))
	## makes a list of mod names

spp_rich_Mod_names<-lapply(XVAR_names, function (x) paste("spp_rich_",paste(x),sep=""))

names(spp_rich_Mods_run)<-spp_rich_Mod_names

	## Save model stats as a spreadsheet of manuscript

spp_rich_Mods_AIC<-lapply(spp_rich_Mods_run,AIC)

spp_rich_Lichen_Count_vs_S_Stats<-cbind(XVAR_names,spp_rich_Mods_AIC,spp_rich_Mods_fit,paste(spp_rich_Formula))
write.csv(spp_rich_Lichen_Count_vs_S_Stats,"Count_spp_rich_S_all_models.csv")

	## Make predicted data files
	## Use a model file ... spp_rich_Mods_run as input
	## Make a SOURCE where there aren't wild estimates at higher S values based on a visual estimate of the 
MyFunc2<-function (x) {as.data.frame(predict(x, SOURCE, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000))}

	##This produces a matrix with a list of 10 objects with 3 observations, each with 5322 entries and 10 columns. 
spp_rich_Mod_pred<-mapply(MyFunc2,spp_rich_Mods_run)%>%as.data.frame()
	## 
	## This works but is hard coded for each predictor file. Have to change working directory each time
	## Make final figure with grey dots in the background
	## Need to add red dots for percentage change


#jpeg(paste("spp_rich_Lichen_Count_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),".jpeg",sep=""))
#tiff(paste("spp_rich_Lichen_Count_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),".tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
#jpeg(paste("spp_rich_Lichen_Count_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.jpeg",sep=""))
#jpeg(paste("spp_rich_Lichen_Count_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.jpeg",sep=""))
#jpeg(paste("spp_rich_Lichen_Count_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_noyaxis_nolegend.jpeg",sep=""))
#tiff(paste("spp_rich_Lichen_Count_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
tiff(paste("spp_rich_Lichen_Count_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_noyaxis_nolegend.tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)

	
x<-spp_rich_Mod_pred[2]
y<-spp_rich_Mod_pred[4]
clim<-unlist(y[1], recursive = F, use.names =T)
x<-unlist(x, recursive = F, use.names =T)
fit_1<-unlist(x[1], recursive = F, use.names =T)
low_1<-unlist(x[2], recursive = F, use.names =T)
high_1<-unlist(x[3], recursive = F, use.names =T)
clim_pred<-unlist(clim[1], recursive = F, use.names =T)
SOURCE_fPlot<-cbind(SOURCE, fit_1, low_1, high_1, clim_pred)
SOURCE_fPlot<-subset(SOURCE_fPlot, S < 20)
palette(c("darkorange","darkgreen"))
plot(spp_rich~S, data=SOURCE, ylab="" ,xlab="CMAQ S kg/ha/yr", main="90th Quantile Regression of Lichen \n Species Richness vs S deposition", col=Area, cex.main = 1.6, cex.lab=1.7, cex.axis = 1.7)
#title(ylab ="Lichen Species Richnesss", cex.lab=1.7, line=2.5)
points(SOURCE_fPlot$S, SOURCE_fPlot$clim_pred, col=alpha('grey',0.25), pch= 19)
points(SOURCE_fPlot$S, SOURCE_fPlot$fit_1, col='blue', pch= 19)
points(SOURCE_fPlot$S, SOURCE_fPlot$low_1, col='black', pch= 4, cex=0.25) 
points(SOURCE_fPlot$S, SOURCE_fPlot$high_1, col='black', pch= 4, cex=0.25)	
	points(0.21,29.361145, col='red', pch=20,cex=2)
	#points(1.51,27.893088,col='red', pch=20,cex=2)
	points(2.90,26.42503, col='red', pch=20,cex=2)
	points(6.00,23.488916,col='red', pch=20,cex=2)
	#points(,4.54, col='red',        pch=20,cex=2)
	#points(,1.82, col='red',        pch=20,cex=2)
    #text(max(SOURCE$S)*0.8,max(SOURCE$spp_rich)*0.70,labels=paste("R1=",print(round(spp_rich_Mods_fit[2],2))))
    #text(max(SOURCE$S)*0.8,max(SOURCE$spp_rich)*0.75,labels=paste(as.character(spp_rich_Formula[2])), cex=1)
    #text(max(SOURCE$S)*0.8,max(SOURCE$spp_rich)*0.65,labels=paste("AIC=",print(round(as.numeric(spp_rich_Mods_AIC[2],2)))))

text(max(SOURCE$S)*0.55,max(SOURCE$spp_rich)*0.65,labels=bquote(atop(.("Lichen Species Richness = 25.84 -"), .("285.70*S+ 83.95*")*S^2*.("- 18.51*")*S^3)), cex=1.5)
#text(max(SOURCE$S)*0.55,max(SOURCE$spp_rich)*0.9,labels=bquote(atop(.("Lichen Species Richness = 25.84 -285.70*S"), .("+ 83.95*")*S^2*.("- 18.51*")*S^3)), cex=1.5)
#legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1)
#legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1.47)
dev.off()

#pdf(paste("spp_rich_Lichen_Count_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),".pdf",sep=""))
#
#MyPlot<-function (x)
#    
#{
#    x<-spp_rich_Mod_pred[i]
#    x<-unlist(x, recursive = F, use.names =T)
#    fit<-unlist(x[1], recursive = F, use.names =T)
#    low<-unlist(x[2], recursive = F, use.names =T)
#    high<-unlist(x[3], recursive = F, use.names =T)
#    SOURCE_fPlot<-cbind(SOURCE, fit, low, high)
#    SOURCE_fPlot<-subset(SOURCE_fPlot, S < 20)
#    plot(spp_rich~S, data=SOURCE, xlab="CMAQ S kg/ha/yr", ylab="Count spp_rich Lichens", main="Count of  spp_rich lichen species vs S deposition")
#    points(SOURCE_fPlot$S, SOURCE_fPlot$fit, col='blue', pch= 19)
#    points(SOURCE_fPlot$S, SOURCE_fPlot$low, col='aquamarine4', pch= 4, cex=0.25) 
#    points(SOURCE_fPlot$S, SOURCE_fPlot$high, col='aquamarine4', pch= 4, cex=0.25)	
#    text(max(SOURCE$S)*0.4,max(SOURCE$spp_rich)*0.8,labels=paste("R1=",print(round(spp_rich_Mods_fit[i],2))))
#    text(max(SOURCE$S)*0.4,max(SOURCE$spp_rich)*0.75,labels=paste(as.character(spp_rich_Formula[i])), cex=0.7)
#    legend('topright', legend=c("Raw Data","Fitted Values", "Prediction Interval"), col=c("black","blue", "aquamarine4"), pch=16)
# 
#}	
#
#for (i in 1:length(XVAR)) MyPlot()
#
#dev.off()

##### Make final version of figure with grey cloud of points in the background of the climate fitted values, the raw data as open black circles
###### and the fitted polynomial line with the predicted fitted 
				
			## x - 5
			## y - spp_rich_Mods_run
			## ModRun - spp_rich_S_poly
			## ModPred spp_rich_Mod_pred
						    
	#Started to automate the function of finding a given X for a value of Y but ran into issues referencing elements using 
	# the characters "$". Do I need to unlist list first?
			
	#GetYforX<-function (x,y,ModRun,ModPred)
	#{		
	#test_s<-as.data.frame(x)
	#colnames(test_s)<-"S"
	#(predict(ModRun, test_s))/max(ModPred)
	#}
	#
	#GetYforX(3,spp_rich_Mods_run,spp_rich_S_poly,spp_rich_Mod_pred)
	
		## Heuristically solved for X given a Y using the polynomial 
			test_s<-as.data.frame(0.5)
			colnames(test_s)<-"S"
			1-(predict(spp_rich_Mods_run$spp_rich_S_poly, test_s))/max(SOURCE_fPlot$fit)

				## 0 % (max spp_rich ) Heuristically solved for  X  S kg/ha/yr given a Y 9.08... 0.21 S kg/ha/yr
				## 5 % loss of N. Cya Heuristically solved for X S kg/ha/yr given a Y 8.62... 0.82 S kg/ha/yr
				## 10% loss of N. Cya Heuristically solved for X S kg/ha/yr given a Y 8.17... 1.44 S kg/ha/yr
				## 20% loss of N. Cya Heuristically solved for X S kg/ha/yr given a Y 7.26... 2.80 S kg/ha/yr
				## 50% loss of N. Cya Heuristically solved for X S kg/ha/yr given a Y 4.54... 5.70 S kg/ha/yr
				## 80% loss of N. Cya Heuristically solved for X S kg/ha/yr given a Y 1.82 ...14.60 S kg/ha/yr

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
				 1-(predict(spp_rich_Mods_run$spp_rich_S_poly, new_S_1kg))/max(SOURCE_fPlot$fit)
				,1-(predict(spp_rich_Mods_run$spp_rich_S_poly,new_S_1.5kg))/max(SOURCE_fPlot$fit)
				,1-(predict(spp_rich_Mods_run$spp_rich_S_poly,new_S_2kg))/max(SOURCE_fPlot$fit)
				,1-(predict(spp_rich_Mods_run$spp_rich_S_poly,new_S_2.5kg))/max(SOURCE_fPlot$fit)
				,1-(predict(spp_rich_Mods_run$spp_rich_S_poly,new_S_3kg))/max(SOURCE_fPlot$fit)
				,1-(predict(spp_rich_Mods_run$spp_rich_S_poly,new_S_5kg))/max(SOURCE_fPlot$fit)
				)
				
				c(
				max(SOURCE_fPlot$fit) ## 0% change
				,max(SOURCE_fPlot$fit)-(max(SOURCE_fPlot$fit)*0.05) ## 5% change
				,max(SOURCE_fPlot$fit)-(max(SOURCE_fPlot$fit)*0.10) ## 10% change
				,max(SOURCE_fPlot$fit)-(max(SOURCE_fPlot$fit)*0.20) ## 20% change
				,max(SOURCE_fPlot$fit)-(max(SOURCE_fPlot$fit)*0.50) ## 50% change
				,max(SOURCE_fPlot$fit)-(max(SOURCE_fPlot$fit)*0.80) ## 80% change
				)
				
