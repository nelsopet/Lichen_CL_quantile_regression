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
#USA - Count of Sensitive vs sulphur 
################################################################################################################

	## Input files are generated in LichenCL_QuantReg_Script.R script. I have to switch each of these sources depending on which measure of Count wanted.
SOURCE<-LichDb_sFncGrpSens_All_Abun_S

  ## Calculate correlations amongst predictors and response
cor_vars<-c("Abn_Cya_sum34_max4","S","maxaug_c","mindec_c","precip_cm","continen","CMD")
select(SOURCE,cor_vars) %>% cor()

  #Make a list of models with different combinations of predictors
XVAR<-c("S"                                                               
 ,"poly(S, 2,raw=T)"                                                      
 ,"S+maxaug_c+mindec_c+precip_cm+CMD"         
 ,"poly(S, 2,raw=T)+maxaug_c+mindec_c+precip_cm+CMD"
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
 ,"CMD"
 ) 
  
  ## Variance Inflation Factor for polynomial model, which ended up being the best model 
VIF(lm(N.sensitive.~poly(S, 2, raw=T)+maxaug_c+mindec_c+precip_cm+CMD, data=SOURCE))

  #Give each set of predictors a short name  
XVAR_names<-c("S","S_poly","S_clim","S_poly_clim","clim","_maxaug_c","S_mindec_c","S_precip_cm","S_continen","S_CMD", "Max C Aug Temp","Min C Dec Temp","Precip_cm","Continentality","CMD");

  #Combine each set of predictors with the response, turn it into a formula and fit 90th quantile regression model
N.sensitive._S_Mods<-paste("N.sensitive. ~ ",paste(XVAR),sep="")
N.sensitive._S_Formula<-lapply(N.sensitive._S_Mods,as.formula)
N.sensitive._S_Mods_run<-lapply(N.sensitive._S_Formula,function (x) rq(x, tau=0.9, data=SOURCE))

	## Make null models
N.sensitive._S_intercept<-rq(formula = N.sensitive. ~ 1, tau = 0.9, data = SOURCE, model = T)

	## Grab each model fit but the intercept is hard coded 
N.sensitive._S_Mods_fit<-mapply(function (y) 1-N.sensitive._S_Mods_run[[y]]$rho/N.sensitive._S_intercept$rho, 1:length(XVAR))

	## makes a list of mod names
N.sensitive._S_Mod_names<-lapply(XVAR_names, function (x) paste("N.sensitive._",paste(x),sep=""))
names(N.sensitive._S_Mods_run)<-N.sensitive._S_Mod_names
  
  ##Calculate AIC for each model and combine R1 statistic
N.sensitive._S_Mods_AIC<-lapply(N.sensitive._S_Mods_run,AIC)
N.sensitive._S_Lichen_Count_vs_S_Stats<-cbind(XVAR_names,N.sensitive._S_Mods_AIC,N.sensitive._S_Mods_fit,paste(N.sensitive._S_Formula))
  
  ## Save model stats as an output file
write.csv(N.sensitive._S_Lichen_Count_vs_S_Stats,"Count_N.sensitive._S_all_models.csv")

	## Make predicted data files
	## Use a model file ... N.sensitive._Mods_run as input
MyFunc2<-function (x) {as.data.frame(predict(x, SOURCE, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000))}

	##This produces a matrix with a list of 10 objects with 3 observations, each with 5322 entries and 10 columns. 
N.sensitive._S_Mod_pred<-mapply(MyFunc2,N.sensitive._S_Mods_run)%>%as.data.frame()

  #Creates a pdf of all models for a given lichen group. Only used for comparing models. Not for final model plotting.
pdf(paste("N.sensitive._S_Lichen_Count_vs_S_percentiles_10thPctile","_",format(Sys.time(),"%y%d%m"),".pdf",sep=""))

MyPlot<-function (x)
{
    x<-N.sensitive._S_Mod_pred[i]
    x<-unlist(x, recursive = F, use.names =T)
    fit<-unlist(x[1], recursive = F, use.names =T)
    low<-unlist(x[2], recursive = F, use.names =T)
    high<-unlist(x[3], recursive = F, use.names =T)
    palette(c("darkorange","darkgreen"))
    plot(N.sensitive.~S, data=SOURCE, xlab="CMAQ S kg/ha/yr", ylab="Count N.sensitive. Lichens", main="Count of  N.sensitive. lichen species vs S deposition", col=Area)
    points(SOURCE$S, fit, col='blue', pch= 19)
    points(SOURCE$S, low, col='black', pch= 4, cex=0.25) 
    points(SOURCE$S, high, col='black', pch= 4, cex=0.25)	
    text(max(SOURCE$S)*0.4,max(SOURCE$N.sensitive.)*0.8,labels=paste("R1=",print(round(N.sensitive._S_Mods_fit[i],2))))
    text(max(SOURCE$S)*0.4,max(SOURCE$N.sensitive.)*0.85,labels=paste("AIC=",print(round(as.numeric(N.sensitive._S_Lichen_Count_vs_S_Stats[i,2]),1))))
    text(max(SOURCE$S)*0.4,max(SOURCE$N.sensitive.)*0.75,labels=paste(as.character(N.sensitive._S_Formula[i])), cex=0.7)
    legend('topright', legend=c("Raw Data East","Raw Data West","Fitted Values", "Prediction Interval"), col=c("darkorange","darkgreen","blue", "black"), pch=16)
}	
for (i in 1:length(XVAR)) MyPlot()
dev.off()

  #Calculate confidence intervals for percent decline in lichen metric for a bunch of deposition increments
new_s_N.sensitive._S_Mods_run<-c(1,1.5,2,2.5,3,5,7.5,10,12.5,15,17.5,20) %>% as.data.frame()
colnames(new_s_N.sensitive._S_Mods_run)<-"S"
100*round(1-(predict(N.sensitive._S_Mods_run$N.sensitive._S_poly,new_s_N.sensitive._S_Mods_run))/  max(SOURCE_S_fPlot$fit),2)
new_s_CI_N.sensitive._S_Mods_run<-predict(N.sensitive._S_Mods_run$N.sensitive._S_poly, new_s_N.sensitive._S_Mods_run, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000)
stuff<-as.data.frame(new_s_CI_N.sensitive._S_Mods_run)
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

  ##Heuristically solved for X given a Y using the polynomial to find deposition associated with each incremental percent decline
test_s<-as.data.frame(0.42)
colnames(test_s)<-"S"
  
  #returns percent decline
1-(predict(N.sensitive._S_Mods_run$N.sensitive._S_poly, test_s)/max(SOURCE_S_fPlot$fit))
  
  #returns absolute decline
predict(N.sensitive._S_Mods_run$N.sensitive._S_poly, test_s)

  #Take deposition associated with incremental percent decline and calculate confidence intervals
s_pct_N.sensitive._S<-c(0.21,0.75	,1.3,2.47	,6.66	,14.05 )%>% as.data.frame()
colnames(s_pct_N.sensitive._S)<-"S"
s_pct_CI_N.sensitive._S<-predict(N.sensitive._S_Mods_run$N.sensitive._S_poly, s_pct_N.sensitive._S, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000)
  
  # print upper and lower bounds of confidence intervals for incremental percent decline
s_pct_CI_N.sensitive._S
  #fit     lower    higher
  #1 17.277259 16.935449 17.618958
  #2 16.405792 16.092860 16.720389
  #3 15.545909 15.251000 15.833526
  #4 13.809758 13.549212 14.057165
  #5  8.630771  8.341216  8.861546
  #6  3.453643  3.240293  3.663014

#Use same heuristic solving approach as above to find deposition associated with upper and lower confidence intervals of incremental precent decline
  # 0 % loss of sensitive lichens Heuristically solved for X given a Y ...  0.005-0.21 -0.42 N kg/ha/yr
  # 5 % loss of sensitive lichens Heuristically solved for X given a Y ...	 0.55 -0.75 -0.945 N kg/ha/yr
  # 10% loss of sensitive lichens Heuristically solved for X given a Y ...	 1.12 -1.3  -1.49 N kg/ha/yr
  # 20% loss of sensitive lichens Heuristically solved for X given a Y ...  2.3  -2.47 -2.65 N kg/ha/yr
  # 50% loss of sensitive lichens Heuristically solved for X given a Y ...  6.45 -6.66 -6.95 N kg/ha/yr
  # 80% loss of sensitive lichens Heuristically solved for X given a Y ...  13.5 -14.05-14.69 N kg/ha/yr
				
### Scatterplots with fitted lines, confidence intervals and percent decline point estimates with and without legends and
### some with and without one of the axis labels. Uncomment the one image driver and one text label and legend (if desired).
### Label and legend positions may need to be adjusted.

      # Get coeffecients to plot in figures 
      round(coefficients(N.sensitive._S_Mods_run$N.sensitive._S_poly),2)
      #17.62                -1.66                 0.05 

      #tiff(paste("N.sensitive._S_Lichen_Count_vs_S_percentiles_US_max20_kgS_ha_yr","_",format(Sys.time(),"%y%d%m"),".tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
      #jpeg(paste("N.sensitive._S_Lichen_Count_vs_S_percentiles_US_max20_kgS_ha_yr",format(Sys.time(),"%y%d%m"),".jpeg",sep=""))
      #jpeg(paste("N.sensitive._S_Lichen_Count_vs_S_percentiles_US_max20_kgS_ha_yr",format(Sys.time(),"%y%d%m"),"_nolegend.jpeg",sep=""))
      tiff(paste("N.sensitive._S_Lichen_Count_vs_S_percentiles_US_max20_kgS_ha_yr","_",format(Sys.time(),"%y%d%m"),"_nolegend.tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
      
      x<-N.sensitive._S_Mod_pred[2]
      x2<-unlist(x, recursive = F, use.names =T)
      y<-N.sensitive._S_Mod_pred[4]
      clim<-unlist(y[1], recursive = F, use.names =T)
      clim_pred<-unlist(clim[1], recursive = F, use.names =T)
      fit<-unlist(x2[1], recursive = F, use.names =T)
      low<-unlist(x2[2], recursive = F, use.names =T)
      high<-unlist(x2[3], recursive = F, use.names =T)
      SOURCE_S_fPlot<-cbind(SOURCE, fit, low, high, clim_pred)
      SOURCE_S_fPlot<-subset(SOURCE_S_fPlot, S < 18)    
      palette(c("darkorange","darkgreen"))
      plot(N.sensitive.~S, data=SOURCE, xlab=expression(Sulfur ~(kg ~ S ~ ha^{-1} ~ y^{-1})), ylab="",  col=Area, cex.main = 1.6, cex.lab=1.7, cex.axis = 1.7)#main="90th Quantile Regression Count of \n Sensitive Species vs S deposition",
      title(ylab ="Count of Sensitive Species", cex.lab=1.7, line=2.5)
      points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$clim_pred, col=alpha('grey',0.25), pch= 19)
      points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$fit, col='blue', pch= 19)
      points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$low, col='black', pch= 4, cex=0.25) 
      points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$high, col='black', pch= 4, cex=0.25)	
      #Red dot coordinates are determined by heuristically solving for incremental percent decline (right value) to determine associated deposition (left value)
      points(0.21	,17.277259      ,col='red', pch=20,cex=2)
      #points(0.75	,16.413396      ,col='red', pch=20,cex=2)
      points(1.3	  ,15.549533      ,col='red', pch=20,cex=2)
      points(2.47	,13.821807      ,col='red', pch=20,cex=2)
      points(6.66	,8.638630       ,col='red', pch=20,cex=2)
      points(14.05 ,3.455452       ,col='red', pch=20,cex=2)
      #text(max(SOURCE$S)*0.60,max(SOURCE$N.sensitive.)*0.45,labels=bquote(atop(.("Count of Sensitive Species = 17.62 -"),.("1.66*S + 0.05*")*S^2)), cex=1.5)
      text(max(SOURCE$S)*0.6,max(SOURCE$N.sensitive.)*0.9,labels=bquote(atop(.("Count of Sensitive Species = 17.62 -"),.("1.66*S + 0.05*")*S^2)), cex=1.5)
      #legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1)
      #legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence\n Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1.47)
      dev.off()				
