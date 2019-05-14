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
#Lichen species richness vs Nitrogen 
################################################################################################################

  ## Input files are generated in LichenCL_QuantReg_Script.R script. I have to switch each of these sources depending on which measure of Count wanted.
SOURCE<-LichDb_sFncGrpSens_All_Abun_N

  ## Calculate correlations amongst predictors and response
cor_vars<-c("spp_rich","N","maxaug_c","mindec_c","precip_cm","continen","CMD")
select(SOURCE,cor_vars) %>% cor()

  #Make a list of models with different combinations of predictors
XVAR<-c("N"                                                               
 ,"poly(N, 2, raw = T)"                                                      
 ,"N+maxaug_c+mindec_c+precip_cm+CMD"         
 ,"poly(N, 2, raw = T)+maxaug_c+mindec_c+precip_cm+CMD"
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
VIF(lm(spp_rich~poly(N, 2, raw=T)+maxaug_c+mindec_c+precip_cm+CMD, data=SOURCE))

  #Give each set of predictors a short name
XVAR_names<-c("N","N_poly","N_clim","N_poly_clim","clim","N_maxaug_c","N_mindec_c","N_precip_cm","N_continen","N_CMD", "Max C Aug Temp","Min C Dec Temp","Precip_cm","Continentality","CMD");

  #Combine each set of predictors with the response, turn it into a formula and fit 90th quantile regression model
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

  ##Calculate AIC for each model and combine R1 statistic
spp_rich_N_Mods_AIC<-lapply(spp_rich_N_Mods_run,AIC)
spp_rich_N_Lichen_Count_vs_N_Stats<-cbind(XVAR_names,spp_rich_N_Mods_AIC,spp_rich_N_Mods_fit,paste(spp_rich_N_Formula))
  
  ## Save model stats as an output file
write.csv(spp_rich_N_Lichen_Count_vs_N_Stats,"Count_spp_rich_N_all_models.csv")

	## Make predicted data files
	## Use a model file ... spp_rich_Mods_run as input
MyFunc2<-function (x) {as.data.frame(predict(x, SOURCE, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000))}

	##This produces a matrix with a list of 10 objects with 3 observations, each with 5322 entries and 10 columns. 
spp_rich_N_Mod_pred<-mapply(MyFunc2,spp_rich_N_Mods_run)%>%as.data.frame()

  #Creates a pdf of all models for a given lichen group. Only used for comparing models. Not for final model plotting.
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
    text(max(SOURCE$N)*0.4,max(SOURCE$spp_rich)*0.8,labels=paste("R1=",print(round(spp_rich_N_Mods_fit[i],2))))
    text(max(SOURCE$N)*0.4,max(SOURCE$spp_rich)*0.75,labels=paste(as.character(spp_rich_N_Formula[i])), cex=0.7)
    legend('topright', legend=c("Raw Data","Fitted Values", "Prediction Interval"), col=c("black","blue", "aquamarine4"), pch=16)
}	
for (i in 1:length(XVAR)) MyPlot()
dev.off()

  ##Calculate confidence intervals for percent decline in lichen metric for a bunch of deposition increments
new_n_spp_rich_N<-c(1,1.5,2,2.5,3,5,7.5,10,12.5,15,17.5,20) %>% as.data.frame()
colnames(new_n_spp_rich_N)<-"N"
100*round(1-(predict(spp_rich_N_Mods_run$spp_rich_N_poly, new_n_spp_rich_N))/max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit),2)
new_n_CI_spp_rich_N<-predict(spp_rich_N_Mods_run$spp_rich_N_poly, new_n_spp_rich_N, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000)
stuff<-as.data.frame(new_n_CI_spp_rich_N)
cbind(
  round((max(SOURCE_N_fPlot$fit)-(stuff$higher))/max(SOURCE_N_fPlot$fit)*100,0)
  ,round((max(SOURCE_N_fPlot$fit)-(stuff$fit))/max(SOURCE_N_fPlot$fit)*100,0)
  ,round((max(SOURCE_N_fPlot$fit)-(stuff$lower))/max(SOURCE_N_fPlot$fit)*100,0)
)
  
  ##Calculate incremental percent decline of each lichen index along the 90th quantile fitted line
c(
max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit) ## 0% change
,max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)-(max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)*0.05) ## 5% change
,max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)-(max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)*0.10) ## 10% change
,max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)-(max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)*0.20) ## 20% change
,max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)-(max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)*0.50) ## 50% change
,max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)-(max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)*0.80) ## 80% change
)

  ## Heuristically solved for X given a Y using the polynomial to find deposition associated with each incremental percent decline
test_n<-as.data.frame(6.65)
colnames(test_n)<-"N"
  
  #returns percent decline
1-(predict(spp_rich_N_Mods_run$spp_rich_N_poly, test_n))/max(spp_rich_N_Mod_pred$spp_rich_N_poly$fit)
  
  #returns absolute decline
predict(spp_rich_N_Mods_run$spp_rich_N_poly, test_n)

  #Take deposition associated with incremental percent decline and calculate confidence intervals
n_pct_spp_rich_N<-c(0.08,0.86,1.70,3.50,6.65,12.8)%>% as.data.frame()
colnames(n_pct_spp_rich_N)<-"N"
n_pct_CI_spp_rich_N<-predict(spp_rich_N_Mods_run$spp_rich_N_poly, n_pct_spp_rich_N, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000)
  
  # print upper and lower bounds of confidence intervals for incremental percent decline
n_pct_CI_spp_rich_N
  # 0 % 32.79772 30.98566 34.62121
  # 5 % 31.14400 29.78724 32.55391
  # 10% 29.46393 28.26562 30.44257
  # 20% 26.21598 25.56544 26.81200
  # 50% 21.68773 20.96845 22.41112
  # 80% 17.08599 16.36759 17.73913

#Use same heuristic solving approach as above to find deposition associated with upper and lower confidence intervals of incremental precent decline
  ## Heuristically solved for X given a Y ...	34.58681
  ## 0 % (max spp_rich) Heuristically solved for X given a Y 10.41...   undef-0.08-0.94  N kg/ha/yr
  ## 5 % loss of N. Cya Heuristically solved for X given a Y 9.89...	  0.19 -0.86 -1.54 N kg/ha/yr
  ## 10% loss of N. Cya Heuristically solved for X given a Y 9.37...	  1.21 -1.70 -2.34 N kg/ha/yr
  ## 20% loss of N. Cya Heuristically solved for X given a Y 8.324995...3.15 -3.50 -3.9 N kg/ha/yr
  ## 50% loss of N. Cya Heuristically solved for X given a Y 5.203122...n/a  N kg/ha/yr
  ## 80% loss of N. Cya Heuristically solved for X given a Y 2.081249...n/a kg/ha/yr

### Scatterplots with fitted lines, confidence intervals and percent decline point estimates with and without legends and
### some with and without one of the axis labels. Uncomment the one image driver and one text label and legend (if desired).
### Label and legend positions may need to be adjusted.

      # Get coeffecients to plot in figures 
      round(coefficients(spp_rich_N_Mods_run$spp_rich_N_poly),2)
      #32.97                -2.19                 0.07 
      
      tiff(paste("spp_rich_N_Lichen_Count_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),".tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
      #jpeg(paste("spp_rich_N_Lichen_Count_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),".jpeg",sep=""))
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
      plot(spp_rich~N, data=SOURCE, xlab=expression(Nitrogen ~ (kg ~ N ~ ha^{-1} ~ y^{-1})), ylab="",  col=Area, cex.main = 1.6, cex.lab=1.7, cex.axis = 1.7)#main="90th Quantile Regresssion of Lichen \n Species Richness vs N deposition",
      title(ylab ="Lichen Species Richnesss", cex.lab=1.7, line=2.5)
      points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$clim_pred, col=alpha('grey',0.25), pch= 19)
      points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$fit, col='blue', pch= 19)
      points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$low, col='black', pch= 4, cex=0.25) 
      points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$high, col='black', pch= 4, cex=0.25)	
      #Red dot coordinates are determined by heuristically solving for incremental percent decline (right value) to determine associated deposition (left value)
      points(0.08,32.80, col='red',     pch=20,cex=2) #0% decline
      #points(0.86,31.16,col='red',     pch=20,cex=2) #5% decline
      points(1.70,29.52, col='red',     pch=20,cex=2) #10% decline
      points(3.50,26.24, col='red',     pch=20,cex=2) #20% decline
      points(6.65,5.203122, col='red', pch=20,cex=2) #50% decline
      points(12.8,2.081249, col='red', pch=20,cex=2) #80% decline
      text(max(SOURCE$N)*0.6,max(SOURCE$spp_rich)*0.65,labels=bquote(atop(.("Lichen Species richness = 32.97 -"),.("2.19*N + 0.07*")*N^2)), cex=1.5)
      #text(max(SOURCE$N)*0.6,max(SOURCE$spp_rich)*0.9,labels=bquote(atop(.("Lichen Species richness = 32.97 -"),.("2.19*N + 0.07*")*N^2)), cex=1.5)
      legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1.47)
      #legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1)
      dev.off()

################################################################################################################
#Lichen species richnesss vs sulphur deposition
################################################################################################################

  ## Input files are generated in LichenCL_QuantReg_Script.R script. I have to switch each of these sources depending on which measure of Count wanted.
SOURCE<-LichDb_sFncGrpSens_All_Abun_S

  ## Calculate correlations amongst predictors and response
cor_vars<-c("spp_rich","S","maxaug_c","mindec_c","precip_cm","continen","CMD")
select(SOURCE,cor_vars) %>% cor()
      
  #Make a list of models with different combinations of predictors
XVAR<-c("S"                                                               
 ,"poly(S, 2, raw = T)"                                                      
 ,"S+maxaug_c+mindec_c+precip_cm+CMD"         
 ,"poly(S, 2, raw = T)+maxaug_c+mindec_c+precip_cm+CMD"
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
VIF(lm(spp_rich~poly(S, 2, raw=T)+maxaug_c+mindec_c+precip_cm+CMD, data=SOURCE))

  #Give each set of predictors a short name
XVAR_names<-c("S","S_poly","S_clim","S_poly_clim","clim","S_maxaug_c","S_mindec_c","S_precip_cm","S_continen","S_CMD", "Max C Aug Temp","Min C Dec Temp","Precip_cm","Continentality","CMD");

  #Combine each set of predictors with the response, turn it into a formula and fit 90th quantile regression model
spp_rich_Mods<-paste("spp_rich ~ ",paste(XVAR),sep="")
spp_rich_Formula<-lapply(spp_rich_Mods,as.formula)
spp_rich_Mods_run<-lapply(spp_rich_Formula,function (x) rq(x, tau=0.9, data=SOURCE))
	
  ## Make null models
spp_rich_S_intercept<-rq(formula = spp_rich ~ 1, tau = 0.9, data = SOURCE, model = T)
	
  ## Grab each model fit but the intercept is hard coded 
spp_rich_Mods_fit<-mapply(function (y) 1-spp_rich_Mods_run[[y]]$rho/spp_rich_S_intercept$rho, 1:length(XVAR))
	
  ## makes a list of mod names
spp_rich_Mod_names<-lapply(XVAR_names, function (x) paste("spp_rich_",paste(x),sep=""))
  
  ##Calculate AIC for each model and combine R1 statistic
names(spp_rich_Mods_run)<-spp_rich_Mod_names
spp_rich_Mods_AIC<-lapply(spp_rich_Mods_run,AIC)

  ## Save model stats as an output file
spp_rich_Lichen_Count_vs_S_Stats<-cbind(XVAR_names,spp_rich_Mods_AIC,spp_rich_Mods_fit,paste(spp_rich_Formula))
write.csv(spp_rich_Lichen_Count_vs_S_Stats,"Count_spp_rich_S_all_models.csv")

	## Make predicted data files
	## Use a model file ... spp_rich_Mods_run as input
MyFunc2<-function (x) {as.data.frame(predict(x, SOURCE, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000))}

	##This produces a matrix with a list of 10 objects with 3 observations, each with 5322 entries and 10 columns. 
spp_rich_Mod_pred<-mapply(MyFunc2,spp_rich_Mods_run)%>%as.data.frame()

  ##Creates a pdf of all models for a given lichen group. Only used for comparing models. Not for final model plotting.
pdf(paste("spp_rich_Lichen_Count_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),".pdf",sep=""))

MyPlot<-function (x)
{
  x<-spp_rich_Mod_pred[i]
  x<-unlist(x, recursive = F, use.names =T)
  fit<-unlist(x[1], recursive = F, use.names =T)
  low<-unlist(x[2], recursive = F, use.names =T)
  high<-unlist(x[3], recursive = F, use.names =T)
  SOURCE_fPlot<-cbind(SOURCE, fit, low, high)
  SOURCE_fPlot<-subset(SOURCE_fPlot, S < 20)
  plot(spp_rich~S, data=SOURCE, xlab="CMAQ S kg/ha/yr", ylab="Count spp_rich Lichens", main="Count of  spp_rich lichen species vs S deposition")
  points(SOURCE_fPlot$S, SOURCE_fPlot$fit, col='blue', pch= 19)
  points(SOURCE_fPlot$S, SOURCE_fPlot$low, col='aquamarine4', pch= 4, cex=0.25) 
  points(SOURCE_fPlot$S, SOURCE_fPlot$high, col='aquamarine4', pch= 4, cex=0.25)	
  text(max(SOURCE$S)*0.4,max(SOURCE$spp_rich)*0.8,labels=paste("R1=",print(round(spp_rich_Mods_fit[i],2))))
  text(max(SOURCE$S)*0.4,max(SOURCE$spp_rich)*0.75,labels=paste(as.character(spp_rich_Formula[i])), cex=0.7)
  legend('topright', legend=c("Raw Data","Fitted Values", "Prediction Interval"), col=c("black","blue", "aquamarine4"), pch=16)
}	
for (i in 1:length(XVAR)) MyPlot()
dev.off()

  ##Calculate confidence intervals for percent decline in lichen metric for a bunch of deposition increments
new_s_spp_rich_S<-c(1,1.5,2,2.5,3,5,7.5,10,12.5,15,17.5,20) %>% as.data.frame()
colnames(new_s_spp_rich_S)<-"S"
100*round(1-(predict(spp_rich_Mods_run$spp_rich_S_poly, new_s_spp_rich_S))/max(SOURCE_fPlot$fit),2)
new_s_CI_spp_rich_S<-predict(spp_rich_Mods_run$spp_rich_S_poly, new_s_spp_rich_S, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000)
stuff<-as.data.frame(new_s_CI_spp_rich_S)
cbind(
   round((max(SOURCE_fPlot$fit)-(stuff$higher))/max(SOURCE_fPlot$fit)*100,0)
  ,round((max(SOURCE_fPlot$fit)-(stuff$fit))/max(SOURCE_fPlot$fit)*100,0)
  ,round((max(SOURCE_fPlot$fit)-(stuff$lower))/max(SOURCE_fPlot$fit)*100,0)
)

  ##Calculate incremental percent decline of each lichen index along the 90th quantile fitted line		  
c(
max(SOURCE_fPlot$fit) ## 0% change
,max(SOURCE_fPlot$fit)-(max(SOURCE_fPlot$fit)*0.05) ## 5% change
,max(SOURCE_fPlot$fit)-(max(SOURCE_fPlot$fit)*0.10) ## 10% change
,max(SOURCE_fPlot$fit)-(max(SOURCE_fPlot$fit)*0.20) ## 20% change
,max(SOURCE_fPlot$fit)-(max(SOURCE_fPlot$fit)*0.50) ## 50% change
,max(SOURCE_fPlot$fit)-(max(SOURCE_fPlot$fit)*0.80) ## 80% change
)

  ## Heuristically solved for X given a Y using the polynomial to find deposition associated with each incremental percent decline
test_s<-as.data.frame(0.01)
colnames(test_s)<-"S"
  
  #returns percent decline
1-(predict(spp_rich_Mods_run$spp_rich_S_poly, test_s))/max(SOURCE_fPlot$fit)
  
  #returns absolute decline
predict(spp_rich_Mods_run$spp_rich_S_poly, test_s)
 
  #Take deposition associated with incremental percent decline and calculate confidence intervals
s_pct_spp_rich_S<-c(0.21,1.51,2.90,6.00)%>% as.data.frame()
colnames(s_pct_spp_rich_S)<-"S"
s_pct_CI_spp_rich_S<-predict(spp_rich_Mods_run$spp_rich_S_poly, s_pct_spp_rich_S, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000)
  
  # print upper and lower bounds of confidence intervals for incremental percent decline
s_pct_CI_spp_rich_S
  # 0 % 29.36114 28.43805 29.96999
  # 5 % 27.88764 27.18389 28.37904
  # 10% 26.40762 25.73268 26.71669
  # 20% 23.46231 22.77852 24.00978

#Use same heuristic solving approach as above to find deposition associated with upper and lower confidence intervals of incremental precent decline
  ## 0 % (max spp_rich ) Heuristically solved for  X  S kg/ha/yr given a Y 9.08... undef-0.21-1.02 S kg/ha/yr
  ## 5 % loss of max spp_rich Heuristically solved for X S kg/ha/yr given a Y 8.62...    1.07 -1.51 -2.23 S kg/ha/yr
  ## 10% loss of max spp_rich Heuristically solved for X S kg/ha/yr given a Y 8.17...    2.6  -2.90  -3.56 S kg/ha/yr
  ## 20% loss of max spp_rich Heuristically solved for X S kg/ha/yr given a Y 7.26...    5.39 -6.00 -6.81 S kg/ha/yr
  ## 50% loss of max spp_rich Heuristically solved for X S kg/ha/yr given a Y 4.54...  n/a S kg/ha/yr
  ## 80% loss of max spp_rich Heuristically solved for X S kg/ha/yr given a Y 1.82 ... n/a S kg/ha/yr

### Scatterplots with fitted lines, confidence intervals and percent decline point estimates with and without legends and
### some with and without one of the axis labels. Uncomment the one image driver and one text label and legend (if desired).
### Label and legend positions may need to be adjusted.

      # Get coeffecients to plot in figures 
      round(coefficients(spp_rich_Mods_run$spp_rich_S_poly),2)
      #29.61                -1.18                 0.03  

      #jpeg(paste("spp_rich_Lichen_Count_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),".jpeg",sep=""))
      #tiff(paste("spp_rich_Lichen_Count_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),".tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
      #jpeg(paste("spp_rich_Lichen_Count_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.jpeg",sep=""))
      #jpeg(paste("spp_rich_Lichen_Count_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_noyaxis_nolegend.jpeg",sep=""))
      tiff(paste("spp_rich_Lichen_Count_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
      #tiff(paste("spp_rich_Lichen_Count_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_noyaxis_nolegend.tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
      
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
      plot(spp_rich~S, data=SOURCE, ylab="" ,xlab=expression(Sulfur ~(kg ~ S ~ ha^{-1} ~ y^{-1})),  col=Area, cex.main = 1.6, cex.lab=1.7, cex.axis = 1.7)#main="90th Quantile Regression of Lichen \n Species Richness vs S deposition",
      title(ylab ="Lichen Species Richnesss", cex.lab=1.7, line=2.5)
      points(SOURCE_fPlot$S, SOURCE_fPlot$clim_pred, col=alpha('grey',0.25), pch= 19)
      points(SOURCE_fPlot$S, SOURCE_fPlot$fit_1, col='blue', pch= 19)
      points(SOURCE_fPlot$S, SOURCE_fPlot$low_1, col='black', pch= 4, cex=0.25) 
      points(SOURCE_fPlot$S, SOURCE_fPlot$high_1, col='black', pch= 4, cex=0.25)	
      #Red dot coordinates are determined by heuristically solving for incremental percent decline (right value) to determine associated deposition (left value)
      points(0.21,29.361145, col='red', pch=20,cex=2) #0% decline
      #points(1.51,27.893088,col='red', pch=20,cex=2) #5% decline
      points(2.90,26.42503, col='red', pch=20,cex=2)  #10% decline
      points(6.00,23.488916,col='red', pch=20,cex=2)  #20% decline
      #points(,4.54, col='red',        pch=20,cex=2)  #50% decline
      #points(,1.82, col='red',        pch=20,cex=2)  #80% decline
  
      #text(max(SOURCE$S)*0.55,max(SOURCE$spp_rich)*0.65,labels=bquote(atop(.("Lichen Species Richness = 29.61 -"), .("1.18*S+ 0.03*")*S^2)), cex=1.5)
      text(max(SOURCE$S)*0.55,max(SOURCE$spp_rich)*0.9,labels=bquote(atop(.("Lichen Species Richness = 29.61 -"), .("1.18*S+ 0.03*")*S^2)), cex=1.5)
      #legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1)
      #legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1.47)
      dev.off()