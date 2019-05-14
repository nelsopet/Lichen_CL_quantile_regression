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
#USA - Count of Oligotroph vs Nitrogen 
#################################################################################################################

  ## Input files are generated in LichenCL_QuantReg_Script.R script. I have to switch each of these sources depending on which measure of Count wanted.
SOURCE<-LichDb_sFncGrpSens_All_Abun_N

  ## Calculate correlations amongst predictors and response
cor_vars<-c("N.oligotroph.","N","maxaug_c","mindec_c","precip_cm","continen","CMD")
select(SOURCE,cor_vars) %>% cor()

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
 ,"CMD"
 ) 

  ## Variance Inflation Factor for polynomial model, which ended up being the best model 
VIF(lm(N.oligotroph.~poly(N, 2, raw=T)+maxaug_c+mindec_c+precip_cm+CMD, data=SOURCE))

  #Give each set of predictors a short name
XVAR_names<-c("N","N_poly","N_clim","N_poly_clim","clim","N_maxaug_c","N_mindec_c","N_precip_cm","N_continen","N_CMD", "Max C Aug Temp","Min C Dec Temp","Precip_cm","Continentality","CMD");

  #Combine each set of predictors with the response, turn it into a formula and fit 90th quantile regression model
N.oligotroph._N_Mods<-paste("N.oligotroph. ~ ",paste(XVAR),sep="")
N.oligotroph._N_Formula<-lapply(N.oligotroph._N_Mods,as.formula)
N.oligotroph._N_Mods_run<-lapply(N.oligotroph._N_Formula,function (x) rq(x, tau=0.9, data=SOURCE))

	## Make null models
N.oligotroph._N_intercept<-rq(formula = N.oligotroph. ~ 1, tau = 0.9, data = SOURCE, model = T)

	## Grab each model fit but the intercept is hard coded 
N.oligotroph._N_Mods_fit<-mapply(function (y) 1-N.oligotroph._N_Mods_run[[y]]$rho/N.oligotroph._N_intercept$rho, 1:length(XVAR))

	## makes a list of mod names
N.oligotroph._N_Mod_names<-lapply(XVAR_names, function (x) paste("N.oligotroph._",paste(x),sep=""))
names(N.oligotroph._N_Mods_run)<-N.oligotroph._N_Mod_names
  
  ##Calculate AIC for each model and combine R1 statistic
N.oligotroph._N_Mods_AIC<-lapply(N.oligotroph._N_Mods_run,AIC)
N.oligotroph._N_Lichen_Count_vs_N_Stats<-cbind(XVAR_names,N.oligotroph._N_Mods_AIC,N.oligotroph._N_Mods_fit,paste(N.oligotroph._N_Formula))
  
  ## Save model stats as an output file
write.csv(N.oligotroph._N_Lichen_Count_vs_N_Stats,"output/Count_N.oligotroph._N_all_models.csv")

	## Make predicted data files
	## Use a model file ... N.oligotroph._Mods_run as input
MyFunc2<-function (x) {as.data.frame(predict(x, SOURCE, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000))}

	##This produces a matrix with a list of 10 objects with 3 observations, each with 5322 entries and 10 columns. 
N.oligotroph._N_Mod_pred<-mapply(MyFunc2,N.oligotroph._N_Mods_run)%>%as.data.frame()

  ##Creates a pdf of all models for a given lichen group. Only used for comparing models. Not for final model plotting.
pdf(paste("output/N.oligotroph._N_Lichen_Count_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),".pdf",sep=""))

MyPlot<-function (x)
{
    x<-N.oligotroph._N_Mod_pred[i]
    x<-unlist(x, recursive = F, use.names =T)
    fit<-unlist(x[1], recursive = F, use.names =T)
    low<-unlist(x[2], recursive = F, use.names =T)
    high<-unlist(x[3], recursive = F, use.names =T)
    palette(c("darkorange","darkgreen"))
    plot(N.oligotroph.~N, data=SOURCE, xlab="CMAQ N kg/ha/yr", ylab="Count N.oligotroph. Lichens", main="Count of  N.oligotroph. lichen species vs N deposition", col=Area)
    points(SOURCE$N, fit, col='blue', pch= 19)
    points(SOURCE$N, low, col='black', pch= 4, cex=0.25) 
    points(SOURCE$N, high, col='black', pch= 4, cex=0.25)	
    text(max(SOURCE$N)*0.4,max(SOURCE$N.oligotroph.)*0.8,labels=paste("R1=",print(round(N.oligotroph._N_Mods_fit[i],2))))
    text(max(SOURCE$N)*0.4,max(SOURCE$N.oligotroph.)*0.85,labels=paste("AIC=",print(round(as.numeric(N.oligotroph._N_Lichen_Count_vs_N_Stats[i,2]),1))))
    text(max(SOURCE$N)*0.4,max(SOURCE$N.oligotroph.)*0.75,labels=paste(as.character(N.oligotroph._N_Formula[i])), cex=0.7)
    legend('topright', legend=c("Raw Data East","Raw Data West","Fitted Values", "Prediction Interval"), col=c("darkorange","darkgreen","blue", "black"), pch=16)
}	
for (i in 1:length(XVAR)) MyPlot()
dev.off()
    
  ##Calculate confidence intervals for percent decline in lichen metric for a bunch of deposition increments
new_n_N.oligotroph._N<-c(1,1.5,2,2.5,3,5,7.5,10,12.5,15,17.5,20) %>% as.data.frame()
colnames(new_n_N.oligotroph._N)<-"N"
100*round(1-(predict(N.oligotroph._N_Mods_run$N.oligotroph._N_poly,new_n_N.oligotroph._N))/max(N.oligotroph._N_Mod_pred$N.oligotroph._N_poly$fit),2)
new_n_CI_N.oligotroph._N<-predict(N.oligotroph._N_Mods_run$N.oligotroph._N_poly, new_n_N.oligotroph._N, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000)
stuff<-as.data.frame(new_n_CI_N.oligotroph._N)
cbind(
   round((max(SOURCE_N_fPlot$fit)-(stuff$higher))/max(SOURCE_N_fPlot$fit)*100,0)
  ,round((max(SOURCE_N_fPlot$fit)-(stuff$fit))/max(SOURCE_N_fPlot$fit)*100,0)
  ,round((max(SOURCE_N_fPlot$fit)-(stuff$lower))/max(SOURCE_N_fPlot$fit)*100,0)
)
  
  ##Calculate incremental percent decline of each lichen index along the 90th quantile fitted line
c(
max(N.oligotroph._N_Mod_pred$N.oligotroph._N_poly$fit) ## 0% change
,max(N.oligotroph._N_Mod_pred$N.oligotroph._N_poly$fit)-(max(N.oligotroph._N_Mod_pred$N.oligotroph._N_poly$fit)*0.05) ## 5% change
,max(N.oligotroph._N_Mod_pred$N.oligotroph._N_poly$fit)-(max(N.oligotroph._N_Mod_pred$N.oligotroph._N_poly$fit)*0.10) ## 10% change
,max(N.oligotroph._N_Mod_pred$N.oligotroph._N_poly$fit)-(max(N.oligotroph._N_Mod_pred$N.oligotroph._N_poly$fit)*0.20) ## 20% change
,max(N.oligotroph._N_Mod_pred$N.oligotroph._N_poly$fit)-(max(N.oligotroph._N_Mod_pred$N.oligotroph._N_poly$fit)*0.50) ## 50% change
,max(N.oligotroph._N_Mod_pred$N.oligotroph._N_poly$fit)-(max(N.oligotroph._N_Mod_pred$N.oligotroph._N_poly$fit)*0.80) ## 80% change
)

  ## Heuristically solved for X given a Y using the polynomial to find deposition associated with each incremental percent decline
test_n<-as.data.frame(0.01)
colnames(test_n)<-"N"
  
  #returns percent decline
1-(predict(N.oligotroph._N_Mods_run$N.oligotroph._N_poly, test_n)/max(N.oligotroph._N_Mod_pred$N.oligotroph._N_poly$fit))
  
  #returns absolute decline
predict(N.oligotroph._N_Mods_run$N.oligotroph._N_poly, test_n)

  #Take deposition associated with incremental percent decline and calculate confidence intervals
n_pct_N.oligotroph._N<-c(0.08,0.81,1.56,3.11,8.29,14.85)%>% as.data.frame()
colnames(n_pct_N.oligotroph._N)<-"N"
n_pct_CI_N.oligotroph._N<-predict(N.oligotroph._N_Mods_run$N.oligotroph._N_poly, n_pct_N.oligotroph._N, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000)
  
  # print upper and lower bounds of confidence intervals for incremental percent decline
n_pct_CI_N.oligotroph._N
  #0 % 18.284188 17.669860 19.018474
  #5 % 17.367441 16.888653 17.919640
  #10% 16.446511 16.019988 16.931364
  #20% 14.610483 14.292127 15.008203
  #50%  9.132006  8.465245  9.595596
  #80%  3.646336  3.287549  4.143290

#Use same heuristic solving approach as above to find deposition associated with upper and lower confidence intervals of incremental precent decline
  ## 0 % loss of oligotrophs Heuristically solved for X given a Y ...  undf -0.08 -0.57 N kg/ha/yr
  ## 5 % loss of oligotrophs Heuristically solved for X given a Y ...	 0.38 -0.81 -1.2 N kg/ha/yr
  ## 10% loss of oligotrophs Heuristically solved for X given a Y ...	 1.17 -1.56 -1.92 N kg/ha/yr
  ## 20% loss of oligotrophs Heuristically solved for X given a Y ...  2.77 -3.11 -3.38 N kg/ha/yr
  ## 50% loss of oligotrophs Heuristically solved for X given a Y ...  7.81 -8.29 -8.99 N kg/ha/yr
  ## 80% loss of oligotrophs Heuristically solved for X given a Y ...  14.16-14.85-15.36 N kg/ha/yr
				
### Scatterplots with fitted lines, confidence intervals and percent decline point estimates with and without legends and
### some with and without one of the axis labels. Uncomment the one image driver and one text label and legend (if desired).
### Label and legend positions may need to be adjusted.

      # Get coeffecients to plot in figures 
      round(coefficients(N.oligotroph._N_Mods_run$N.oligotroph._N_poly),2)

      #tiff(paste("output/N.oligotroph._N_Lichen_Count_vs_N_percentiles_US_","_",format(Sys.time(),"%y%d%m"),".tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
      #jpeg(paste("output/N.oligotroph._N_Lichen_Count_vs_N_percentiles_US_",format(Sys.time(),"%y%d%m"),".jpeg",sep=""))
      #jpeg(paste("output/N.oligotroph._N_Lichen_Count_vs_N_percentiles_US_",format(Sys.time(),"%y%d%m"),"_nolegend.jpeg",sep=""))
       tiff(paste("output/N.oligotroph._N_Lichen_Count_vs_N_percentiles_US_","_",format(Sys.time(),"%y%d%m"),"_nolegend.tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
      
      x<-N.oligotroph._N_Mod_pred[2]
      x2<-unlist(x, recursive = F, use.names =T)
      y<-N.oligotroph._N_Mod_pred[4]
      clim<-unlist(y[1], recursive = F, use.names =T)
      clim_pred<-unlist(clim[1], recursive = F, use.names =T)
      
      fit<-unlist(x2[1], recursive = F, use.names =T)
      low<-unlist(x2[2], recursive = F, use.names =T)
      high<-unlist(x2[3], recursive = F, use.names =T)
      SOURCE_N_fPlot<-cbind(SOURCE, fit, low, high, clim_pred)
      #SOURCE_N_fPlot<-subset(SOURCE_N_fPlot, N < 12.5)    
      palette(c("darkorange","darkgreen"))
      plot(N.oligotroph.~N, data=SOURCE, xlab=expression(Nitrogen ~ (kg ~ N ~ ha^{-1} ~ y^{-1})), ylab="",  col=Area, cex.main = 1.6, cex.lab=1.7, cex.axis = 1.7)#main="90th Quantile Regresssion Count of \n Western Oligotroph Species vs N deposition",
      title(ylab ="Count of Oligotrophs", cex.lab=1.7, line=2.5)
      points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$clim_pred, col=alpha('grey',0.25), pch= 19)
      points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$fit, col='blue', pch= 19)
      points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$low, col='black', pch= 4, cex=0.25) 
      points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$high, col='black', pch= 4, cex=0.25)	
      #Red dot coordinates are determined by heuristically solving for incremental percent decline (right value) to determine associated deposition (left value)
      points(0.08,	18.284188, col='red', pch=20,cex=2) #0% decline
      #points(0.81,	17.369978, col='red', pch=20,cex=2) #5% decline
      points(1.56,	16.455769, col='red', pch=20,cex=2) #10% decline
      points(3.11,	14.627350, col='red', pch=20,cex=2) #20% decline
      points(8.29,	9.142094, col='red', pch=20,cex=2) #50% decline
      points(14.85,	3.656838, col='red',  pch=20,cex=2) #80% decline
     
      #text(max(SOURCE$N)*0.70,max(SOURCE$N.oligotroph.)*0.60,labels=bquote(atop(.("Count of Oligotrophs = 18.39 -"),.("1.27*N + 0.02*")*N^2)), cex=1.5)
      #text(max(SOURCE$N)*0.70,max(SOURCE$N.oligotroph.)*0.60,labels=bquote(atop(.("Count of Oligotrophs = 18.39 -"),.("1.27*N + 0.02*")*N^2)), cex=1.5)
      text(max(SOURCE$N)*0.65,max(SOURCE$N.oligotroph.)*0.9,labels=bquote(atop(.("Count of Oligotrophs = 18.39 -"),.("1.27*N + 0.02*")*N^2)), cex=1.5)
      #legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1)
      #legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence\n Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1.47)
      dev.off()
