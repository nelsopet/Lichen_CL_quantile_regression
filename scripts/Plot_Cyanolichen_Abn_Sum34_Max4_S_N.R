## Install required packages
#Function taken from https://nhsrcommunity.com/blog/a-simple-function-to-install-and-load-packages-in-r/
install_or_load_pack <- function(pack){
  
  create.pkg <- pack[!(pack %in% installed.packages()[, "Package"])]
  
  if (length(create.pkg))
    
    install.packages(create.pkg, dependencies = TRUE)
  
  sapply(pack, require, character.only = TRUE)
  
  #I know I should be using purr here, but this is before the Tidyverse is loaded. I know you Tidyverse trend setters will have me here.
  
}

#Make the list of packages needed here
packs<-c("plyr","quantreg","tidyverse","scales","vegan","fmsb")

install_or_load_pack(packs)

############################################################################################################################################
#Abundance of Cyanolichen vs Nitrogen 
################################################################################################################

  ## Input files are generated in LichenCL_QuantReg_Script.R script. I have to switch each of these sources depending on which measure of Abundance wanted.	SOURCE<-LichDb_sFncGrpSens_All_Abun_N
SOURCE<-LichDb_sFncGrpSens_All_Abun_N

  ## Calculate correlations amongst predictors and response
cor_vars<-c("Abn_Cya_sum34_max4","N","maxaug_c","mindec_c","precip_cm","continen","CMD")
select(SOURCE, cor_vars) %>% cor()

  #Make a list of models with different combinations of predictors	
XVAR<-c("N"                                                               
 ,"poly(N, 2, raw=T)"                                                      
 ,"N+maxaug_c+mindec_c+precip_cm+CMD"         
 ,"poly(N, 2, raw=T)+maxaug_c+mindec_c+precip_cm+CMD"
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

  ## Variance Inflation Factor for polynomial model, which ended up being the best model 
VIF(lm(Abn_Cya_sum34_max4~poly(N, 2, raw=T)+maxaug_c+mindec_c+precip_cm+CMD, data=SOURCE))

  #Give each set of predictors a short name
XVAR_names<-c("N","N_poly","N_clim","N_poly_clim","clim","N_maxaug_c","N_mindec_c","N_precip_cm","N_continen","N_CMD", "Max C Aug Temp","Min C Dec Temp","Precip_cm","Continentality","CMD");

  #Combine each set of predictors with the response, turn it into a formula and fit 90th quantile regression model
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
  
  ##Calculate AIC for each model and combine R1 statistic
Abn_Cya_sum34_max4_N_Mods_AIC<-lapply(Abn_Cya_sum34_max4_N_Mods_run,AIC)
Abn_Cya_sum34_max4_N_Lichen_Abundance_vs_N_Stats<-cbind(XVAR_names,Abn_Cya_sum34_max4_N_Mods_AIC,Abn_Cya_sum34_max4_N_Mods_fit,paste(Abn_Cya_sum34_max4_N_Formula))
  
  ## Save model stats as an output file
write.csv(Abn_Cya_sum34_max4_N_Lichen_Abundance_vs_N_Stats,"output/Abundance_Abn_Cya_sum34_max4_N_all_models.csv")

	## Make predicted data files
	## Use a model file ... Abn_Cya_sum34_max4_Mods_run as input
MyFunc2<-function (x) {as.data.frame(predict(x, SOURCE, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000))}

	##This produces a matrix with a list of 10 objects with 3 observations, each with 5322 entries and 10 columns. 
Abn_Cya_sum34_max4_N_Mod_pred<-mapply(MyFunc2,Abn_Cya_sum34_max4_N_Mods_run)%>%as.data.frame()

  #Creates a pdf of all models for a given lichen group. Only used for comparing models. Not for final model plotting.
pdf(paste("output/Abn_Cya_sum34_max4_N_Lichen_Abundance_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),".pdf",sep=""))
MyPlot<-function (x)
{
    x<-Abn_Cya_sum34_max4_N_Mod_pred[i]
    x<-unlist(x, recursive = F, use.names =T)
    fit<-unlist(x[1], recursive = F, use.names =T)
    low<-unlist(x[2], recursive = F, use.names =T)
    high<-unlist(x[3], recursive = F, use.names =T)
    plot(Abn_Cya_sum34_max4~N, data=SOURCE, xlab="CMAQ N kg/ha/yr", ylab="Abundance Abn_Cya_sum34_max4 Lichens", main="Abundance of  Abn_Cya_sum34_max4 lichen species vs N deposition")
    points(SOURCE$N, fit, col='blue', pch= 19)
    points(SOURCE$N, low, col='aquamarine4', pch= 4, cex=0.25) 
    points(SOURCE$N, high, col='aquamarine4', pch= 4, cex=0.25)	
    text(max(SOURCE$N)*0.4,max(SOURCE$Abn_Cya_sum34_max4)*0.8,labels=paste("R1=",print(round(Abn_Cya_sum34_max4_N_Mods_fit[i],2))))
    text(max(SOURCE$N)*0.4,max(SOURCE$Abn_Cya_sum34_max4)*0.75,labels=paste(as.character(Abn_Cya_sum34_max4_N_Formula[i])), cex=0.7)
    legend('topright', legend=c("Raw Data","Fitted Values", "Prediction Interval"), col=c("black","blue", "aquamarine4"), pch=16)
}	
for (i in 1:length(XVAR)) MyPlot()
dev.off()

  #Calculate confidence intervals for percent decline in lichen metric for a bunch of deposition increments
new_n_Abn_Cya_sum34_max4<-c(1,1.5,2,2.5,3,5,7.5,10,12.5,15,17.5,20) %>% as.data.frame()
colnames(new_n_Abn_Cya_sum34_max4)<-"N"
100*round(1-(predict(Abn_Cya_sum34_max4_N_Mods_run$Abn_Cya_sum34_max4_N_poly, new_n_Abn_Cya_sum34_max4))/max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit),2)
new_n_CI_Abn_Cya_sum34_max4<-predict(Abn_Cya_sum34_max4_N_Mods_run$Abn_Cya_sum34_max4_N_poly, new_n_Abn_Cya_sum34_max4, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000)
stuff<-as.data.frame(new_n_CI_Abn_Cya_sum34_max4)
cbind(
  round((max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(stuff$higher))/max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*100,0)
  ,round((max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(stuff$fit))/max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*100,0)
  ,round((max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(stuff$lower))/max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*100,0)
)

  ##Calculate incremental percent decline of each lichen index along the 90th quantile fitted line
c(
max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit) ## 0% change
,max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*0.05) ## 5% change
,max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*0.10) ## 10% change
,max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*0.20) ## 20% change
,max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*0.50) ## 50% change
,max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)-(max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)*0.80) ## 80% change
)

  #Heuristically solved for X given a Y using the polynomial to find deposition associated with each incremental percent decline
test_n<-as.data.frame(0.2)
colnames(test_n)<-"N"
  
  #returns percent decline
1-(predict(Abn_Cya_sum34_max4_N_Mods_run$Abn_Cya_sum34_max4_N_poly, test_n))/max(Abn_Cya_sum34_max4_N_Mod_pred$Abn_Cya_sum34_max4_N_poly$fit)
  
  #returns absolute decline
predict(Abn_Cya_sum34_max4_N_Mods_run$Abn_Cya_sum34_max4_N_poly, test_n)

  #Take deposition associated with incremental percent decline and calculate confidence intervals
n_pct_Abn_Cya_sum34_max4_N<-c(0.08,0.4, 0.7 ,1.3 ,3.5 ,6.6)%>% as.data.frame()
colnames(n_pct_Abn_Cya_sum34_max4_N)<-"N"
n_pct_CI_Abn_Cya_sum34_max4_N<-predict(Abn_Cya_sum34_max4_N_Mods_run$Abn_Cya_sum34_max4_N_poly, n_pct_Abn_Cya_sum34_max4_N, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000)
  
  # print upper and lower bounds of confidence intervals for incremental percent decline
n_pct_CI_Abn_Cya_sum34_max4_N
  # 0 % 13.495629 13.434098 15.877696
  # 5 % 12.770348 12.739062 14.966161
  # 10% 12.108583 12.104228 14.136051
  # 20% 10.837854 10.770273 12.564827
  # 50%  6.780830  6.515438  7.635893
  # 80%  2.670614  2.322223  2.956652

#Use same heuristic solving approach as above to find deposition associated with upper and lower confidence intervals of incremental precent decline
  ## 0 % max	of N. Cya Heuristically solved for X given a Y ...13.49563 Abn_Cya_sum34_max4 undef-0.08-0.1 N kg/ha/yr
  ## 5 % loss of N. Cya Heuristically solved for X given a Y ...12.82085 Abn_Cya_sum34_max4 undef-0.4- 0.414  N kg/ha/yr
  ## 10% loss of N. Cya Heuristically solved for X given a Y ...12.14607 Abn_Cya_sum34_max4 undef-0.7- 0.7   N kg/ha/yr
  ## 20% loss of N. Cya Heuristically solved for X given a Y ...10.7965  Abn_Cya_sum34_max4 0.49 -1.3- 1.33  N kg/ha/yr
  ## 50% loss of N. Cya Heuristically solved for X given a Y ...6.747815 Abn_Cya_sum34_max4 2.99 -3.5- 3.67  N kg/ha/yr
  ## 80% loss of N. Cya Heuristically solved for X given a Y ...2.699126 Abn_Cya_sum34_max4 6.33 -6.6- 6.95  N kg/ha/yr

### Scatterplots with fitted lines, confidence intervals and percent decline point estimates with and without legends and
### some with and without one of the axis labels. Uncomment the one image driver and one text label and legend (if desired).
### Label and legend positions may need to be adjusted.

        # Get coeffecients to plot in figures 
        round(coefficients(Abn_Cya_sum34_max4_N_Mods_run$Abn_Cya_sum34_max4_N_poly),2)

				#tiff(paste("output/Abn_Cya_sum34_max4_N_Cyanolichen_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),".tiff", sep=""), width=1200,  height=1200, units="px", pointsize = 24)
				#jpeg(paste("output/Abn_Cya_sum34_max4_N_Cyanolichen_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),".jpeg", sep=""))
				#jpeg(paste("output/Abn_Cya_sum34_max4_N_Cyanolichen_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.jpeg", sep=""))
				 tiff(paste("output/Abn_Cya_sum34_max4_N_Cyanolichen_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.tiff", sep=""), width=1200,  height=1200, units="px", pointsize = 24)
				
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
				plot(Abn_Cya_sum34_max4~N, data=SOURCE, xlab=expression(Nitrogen ~ (kg ~ N ~ ha^{-1} ~ y^{-1})), ylab="",  col=Area, cex.main = 1.6, cex.lab=1.7, cex.axis = 1.7)#main="90th Quantile Regression of \n Cyanolichen Abundance vs N deposition",
				title(ylab ="Cyanolichen Abundance", cex.lab=1.7, line=2.5)
				points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$clim_pred, col=alpha('grey',0.25), pch= 19)
				points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$fit, col='blue', pch= 19)
				points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$low, col='black', pch= 4, cex=0.25) 
				points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$high, col='black', pch= 4, cex=0.25)	
				#Red dot coordinates are determined by heuristically solving for incremental percent decline (right value) to determine associated deposition (left value)
				points(0.08,13.49563,col='red', pch=20,cex=2) #0% decline
				#points(0.4, 12.82085,col='red', pch=20,cex=2) #5% decline
				points(0.7 ,12.14607,col='red', pch=20,cex=2) #10% decline
				points(1.3 ,10.7965 ,col='red', pch=20,cex=2) #20% decline
				points(3.5 ,6.747815,col='red', pch=20,cex=2) #50% decline
				points(6.6 ,2.699126,col='red', pch=20,cex=2)  #80% decline
				#text(max(SOURCE$N)*0.75,max(SOURCE$Abn_Cya_sum34_max4)*0.85,labels=paste("90th quantile regression model"))
				#text(max(SOURCE$N)*0.5,max(SOURCE$Abn_Cya_sum34_max4)*0.6,labels=bquote(atop(.("Cyanolichen Abundance = 13.68-"),.("2.31*N+0.10*")*N^2)), cex=1.5)
				text(max(SOURCE$N)*0.5,max(SOURCE$Abn_Cya_sum34_max4)*0.9,labels=bquote(atop(.("Cyanolichen Abundance = 13.68-"),.("2.31*N+0.10*")*N^2)), cex=1.5)
				#legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1.47)
				#legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1)
				
				dev.off()
				
################################################################################################################
#Abundance of cyanolichens vs sulphur deposition
################################################################################################################

  ## Input files are generated in LichenCL_QuantReg_Script.R script. I have to switch each of these sources depending on which measure of Abundance wanted.	SOURCE<-LichDb_sFncGrpSens_All_Abun_S
SOURCE<-LichDb_sFncGrpSens_All_Abun_S

  ## Calculate correlations amongst predictors and response
cor_vars<-c("Abn_Cya_sum34_max4","S","maxaug_c","mindec_c","precip_cm","continen","CMD")
select(SOURCE,cor_vars) %>% cor()				
	
  #Make a list of models with different combinations of predictors			
XVAR<-c("S"                                                               
 ,"poly(S, 2, raw=T)"                                                      
 ,"S+maxaug_c+mindec_c+precip_cm+CMD"         
 ,"poly(S, 2, raw=T)+maxaug_c+mindec_c+precip_cm+CMD"
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

  ## Variance Inflation Factor for polynomial model, which ended up being the best model 
VIF(lm(Abn_Cya_sum34_max4~poly(S, 2, raw=T)+maxaug_c+mindec_c+precip_cm+CMD, data=SOURCE))

  #Give each set of predictors a short name
XVAR_names<-c("S","S_poly","S_clim","S_poly_clim","clim","S_maxaug_c","S_mindec_c","S_precip_cm","S_continen","S_CMD", "Max C Aug Temp","Min C Dec Temp","Precip_cm","Continentality","CMD");

  #Combine each set of predictors with the response, turn it into a formula and fit 90th quantile regression model
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

  ##Calculate AIC for each model and combine R1 statistic
Abn_Cya_S_sum34_max4_Mods_AIC<-lapply(Abn_Cya_S_sum34_max4_Mods_run,AIC)
Abn_Cya_sum34_max4_Lichen_Abundance_vs_S_Stats<-cbind(XVAR_names,Abn_Cya_S_sum34_max4_Mods_AIC,Abn_Cya_S_sum34_max4_Mods_fit,paste(Abn_Cya_S_sum34_max4_Formula))
    
  ## Save model stats as an output file
write.csv(Abn_Cya_sum34_max4_Lichen_Abundance_vs_S_Stats,"output/Abundance_Abn_Cya_sum34_max4_S_all_models.csv")

	## Make predicted data files
	## Use a model file ... Abn_Cya_S_sum34_max4_Mods_run as input
MyFunc2<-function (x) {as.data.frame(predict(x, SOURCE, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000))}

	##This produces a matrix with a list of 10 objects with 3 observations, each with 5322 entries and 10 columns. 
Abn_Cya_S_sum34_max4_Mod_pred<-mapply(MyFunc2,Abn_Cya_S_sum34_max4_Mods_run)%>%as.data.frame()

  #Creates a pdf of all models for a given lichen group. Only used for comparing models. Not for final model plotting.
pdf(paste("output/Abn_Cya_sum34_max4_Lichen_Abundance_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),".pdf",sep=""))

MyPlot<-function (x)
{
    x<-Abn_Cya_S_sum34_max4_Mod_pred[i]
    x<-unlist(x, recursive = F, use.names =T)
    fit<-unlist(x[1], recursive = F, use.names =T)
    low<-unlist(x[2], recursive = F, use.names =T)
    high<-unlist(x[3], recursive = F, use.names =T)
    SOURCE_S_fPlot<-cbind(SOURCE, fit, low, high)
    SOURCE_S_fPlot<-subset(SOURCE_S_fPlot, S < 20)
    plot(Abn_Cya_sum34_max4~S, data=SOURCE, xlab="CMAQ S kg/ha/yr", ylab="Abundance Abn_Cya_sum34_max4 Lichens", main="Abundance of  Abn_Cya_sum34_max4 lichen species vs S deposition")
    points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$fit, col='blue', pch= 19)
    points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$low, col='aquamarine4', pch= 4, cex=0.25) 
    points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$high, col='aquamarine4', pch= 4, cex=0.25)	
    text(max(SOURCE$S)*0.4,max(SOURCE$Abn_Cya_sum34_max4)*0.8,labels=paste("R1=",print(round(Abn_Cya_S_sum34_max4_Mods_fit[i],2))))
    text(max(SOURCE$S)*0.4,max(SOURCE$Abn_Cya_sum34_max4)*0.75,labels=paste(as.character(Abn_Cya_S_sum34_max4_Formula[i])), cex=0.7)
    legend('topright', legend=c("Raw Data","Fitted Values", "Prediction Interval"), col=c("black","blue", "aquamarine4"), pch=16)
}	
for (i in 1:length(XVAR)) MyPlot()
dev.off()

  #Make object with fitted values for polynomial model for use  in calculating response to deposition increments
x<-Abn_Cya_S_sum34_max4_Mod_pred[2]
x<-unlist(x, recursive = F, use.names =T)
fit<-unlist(x[1], recursive = F, use.names =T)
low<-unlist(x[2], recursive = F, use.names =T)
high<-unlist(x[3], recursive = F, use.names =T)
SOURCE_S_fPlot<-cbind(SOURCE, fit, low, high)
SOURCE_S_fPlot<-subset(SOURCE_S_fPlot, S < 20)

  #Calculate confidence intervals for percent decline in lichen metric for a bunch of deposition increments
new_s_Abn_Cya_sum34_max4_S<-c(1,1.5,2,2.5,3,5,7.5,10,12.5,15,17.5,20) %>% as.data.frame()
colnames(new_s_Abn_Cya_sum34_max4_S)<-"S"
100*round(1-(predict(Abn_Cya_S_sum34_max4_Mods_run$Abn_Cya_sum34_max4_S_poly,new_s_Abn_Cya_sum34_max4_S))/max(SOURCE_S_fPlot$fit),2)
new_s_CI_Abn_Cya_sum34_max4_S<-predict(Abn_Cya_S_sum34_max4_Mods_run$Abn_Cya_sum34_max4_S_poly, new_s_Abn_Cya_sum34_max4_S, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000)
stuff<-as.data.frame(new_s_CI_Abn_Cya_sum34_max4_S)
cbind(
   round((max(SOURCE_S_fPlot$fit)-(stuff$higher))/max(SOURCE_S_fPlot$fit)*100,0)
  ,round((max(SOURCE_S_fPlot$fit)-(stuff$fit))/max(SOURCE_S_fPlot$fit)*100,0)
  ,round((max(SOURCE_S_fPlot$fit)-(stuff$lower))/max(SOURCE_S_fPlot$fit)*100,0)
)

  ##Calculate incremental percent decline of each lichen index along the 90th quantile fitted line			
c(
max(SOURCE_S_fPlot$fit) ## 0% change
,max(SOURCE_S_fPlot$fit)-(max(SOURCE_S_fPlot$fit)*0.05) ## 5% change
,max(SOURCE_S_fPlot$fit)-(max(SOURCE_S_fPlot$fit)*0.10) ## 10% change
,max(SOURCE_S_fPlot$fit)-(max(SOURCE_S_fPlot$fit)*0.20) ## 20% change
,max(SOURCE_S_fPlot$fit)-(max(SOURCE_S_fPlot$fit)*0.50) ## 50% change
,max(SOURCE_S_fPlot$fit)-(max(SOURCE_S_fPlot$fit)*0.80) ## 80% change
)

  #Heuristically solved for X given a Y using the polynomial to find deposition associated with each incremental percent decline
test_s<-as.data.frame(0.05)
colnames(test_s)<-"S"
  
  #returns percent decline
1-(predict(Abn_Cya_S_sum34_max4_Mods_run$Abn_Cya_sum34_max4_S_poly, test_s))/max(SOURCE_S_fPlot$fit)
  
  #returns absolute decline
predict(Abn_Cya_S_sum34_max4_Mods_run$Abn_Cya_sum34_max4_S_poly, test_s)

  #Take deposition associated with incremental percent decline and calculate confidence intervals
s_pct_Abn_Cya_sum34_max4_S<-c(0.21,0.7, 1.2 ,2.3 ,5.9 ,10.9)%>% as.data.frame()
colnames(s_pct_Abn_Cya_sum34_max4_S)<-"S"
s_pct_CI_Abn_Cya_sum34_max4_S<-predict(Abn_Cya_S_sum34_max4_Mods_run$Abn_Cya_sum34_max4_S_poly, s_pct_Abn_Cya_sum34_max4_S, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000)

  # print upper and lower bounds of confidence intervals for incremental percent decline
s_pct_CI_Abn_Cya_sum34_max4_S
  # 0 % 9.733974 9.237599 10.329497
  # 5 % 9.242846 8.788535  9.797999
  # 10% 8.754785 8.341838  9.275963
  # 20% 7.727596 7.301401  8.177313
  # 50% 4.813352 4.463998  5.116947
  # 80% 1.902970 1.584362  2.187304

#Use same heuristic solving approach as above to find deposition associated with upper and lower confidence intervals of incremental precent decline
  ## 0 % loss of N. Cya Heuristically solved for X S kg/ha/yr given a Y  9.733974 ... undef-0.21-0.71 S kg/ha/yr
  ## 5 % loss of N. Cya Heuristically solved for X S kg/ha/yr given a Y  9.247276 ... 0.15 -0.7 -1.17  S kg/ha/yr
  ## 10% loss of N. Cya Heuristically solved for X S kg/ha/yr given a Y  8.760577 ... 0.67 -1.2 -1.63  S kg/ha/yr
  ## 20% loss of N. Cya Heuristically solved for X S kg/ha/yr given a Y  7.78718  ... 1.81 -2.3 -2.77  S kg/ha/yr
  ## 50% loss of N. Cya Heuristically solved for X S kg/ha/yr given a Y  4.866987 ... 5.49 -5.9 -6.4  S kg/ha/yr
  ## 80% loss of N. Cya Heuristically solved for X S kg/ha/yr given a Y  1.946795 ... 10.3 -10.9-11.65 S kg/ha/yr

### Scatterplots with fitted lines, confidence intervals and percent decline point estimates with and without legends and
### some with and without one of the axis labels. Uncomment the one image driver and one text label and legend (if desired).
### Label and legend positions may need to be adjusted.

      # Get coeffecients to plot in figures 
			round(coefficients(Abn_Cya_S_sum34_max4_Mods_run$Abn_Cya_sum34_max4_S_poly),2)
      #9.95                -1.03                 0.03 
				
      #tiff(paste("output/Abn_Cya_sum34_max4_S_Cyanolichen_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),".tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
      #jpeg(paste("output/Abn_Cya_sum34_max4_S_Cyanolichen_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),".jpeg",sep=""))
      #jpeg(paste("output/Abn_Cya_sum34_max4_S_Cyanolichen_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.jpeg",sep=""))
      tiff(paste("output/Abn_Cya_sum34_max4_S_Cyanolichen_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
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
      #Red dot coordinates are determined by heuristically solving for incremental percent decline (right value) to determine associated deposition (left value)
      points(0.21, 9.733974, col='red', pch=20,cex=2) #0% decline
      #points(0.7, 9.247276 ,col='red', pch=20,cex=2) #5% decline
      points(1.2 , 8.760577 ,col='red', pch=20,cex=2) #10% decline
      points(2.3 , 7.78718  ,col='red', pch=20,cex=2) #20% decline
      points(5.9 , 4.866987 ,col='red', pch=20,cex=2) #50% decline
      points(10.9, 1.946795 ,col='red', pch=20,cex=2) #80% decline
      
      #text(max(SOURCE$S)*0.5,max(SOURCE$Abn_Cya_sum34_max4)*0.6,labels=bquote(atop(.("Cyanolichen Abundance = 9.95 -"), .("1.03*S + 0.03*")*S^2)), cex=1.5)
      text(max(SOURCE$S)*0.5,max(SOURCE$Abn_Cya_sum34_max4)*0.9,labels=bquote(atop(.("Cyanolichen Abundance = 9.95 -"), .("1.03*S + 0.03*")*S^2)), cex=1.5)
      #legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1.47)
      #legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1)
      dev.off()			