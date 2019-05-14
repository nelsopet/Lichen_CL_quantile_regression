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
#Abundance of Forage Lichen vs Nitrogen 
################################################################################################################

  ## Input files are generated in LichenCL_QuantReg_Script.R script. I have to switch each of these sources depending on which measure of Abundance wanted.
SOURCE<-LichDb_sFncGrpSens_All_Abun_N

  ## Calculate correlations amongst predictors and response
cor_vars<-c("Abn_For_sum34_max4","N","maxaug_c","mindec_c","precip_cm","continen","CMD")
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
,"CMD") 

  ## Variance Inflation Factor for polynomial model, which ended up being the best model 
VIF(lm(Abn_For_sum34_max4~poly(N, 2, raw=T)+maxaug_c+mindec_c+precip_cm+CMD, data=SOURCE))

  #Give each set of predictors a short name
XVAR_names<-c("N","N_poly","N_clim","N_poly_clim","clim","N_maxaug_c","N_mindec_c","N_precip_cm","N_continen","N_CMD", "Max C Aug Temp","Min C Dec Temp","Precip_cm","Continentality","CMD");
  
  #Combine each set of predictors with the response, turn it into a formula and fit 90th quantile regression model
Abn_For_sum34_max4_N_Mods<-paste("Abn_For_sum34_max4 ~ ",paste(XVAR),sep="")
Abn_For_sum34_max4_N_Formula<-lapply(Abn_For_sum34_max4_N_Mods,as.formula)
Abn_For_sum34_max4_N_Mods_run<-lapply(Abn_For_sum34_max4_N_Formula,function (x) rq(x, tau=0.9, data=SOURCE))

	## Make null models
Abn_For_sum34_max4_N_intercept<-rq(formula = Abn_For_sum34_max4 ~ 1, tau = 0.9, data = SOURCE, model = T)

	## Grab each model fit but the intercept is hard coded 
Abn_For_sum34_max4_N_Mods_fit<-mapply(function (y) 1-Abn_For_sum34_max4_N_Mods_run[[y]]$rho/Abn_For_sum34_max4_N_intercept$rho, 1:length(XVAR))

	## makes a list of mod names
Abn_For_sum34_max4_N_Mod_names<-lapply(XVAR_names, function (x) paste("Abn_For_sum34_max4_",paste(x),sep=""))
names(Abn_For_sum34_max4_N_Mods_run)<-Abn_For_sum34_max4_N_Mod_names
  
  ##Calculate AIC for each model and combine R1 statistic
Abn_For_sum34_max4_N_Mods_AIC<-lapply(Abn_For_sum34_max4_N_Mods_run,AIC)
Abn_For_sum34_max4_N_Lichen_Abundance_vs_N_Stats<-cbind(XVAR_names,Abn_For_sum34_max4_N_Mods_AIC,Abn_For_sum34_max4_N_Mods_fit,paste(Abn_For_sum34_max4_N_Formula))
  
  ## Save model stats as an output file
write.csv(Abn_For_sum34_max4_N_Lichen_Abundance_vs_N_Stats,"Abundance_Abn_For_sum34_max4_N_all_models.csv")

	## Make predicted data files
	## Use a model file ... Abn_For_sum34_max4_Mods_run as input
MyFunc2<-function (x) {as.data.frame(predict(x, SOURCE, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000))}

	##This produces a matrix with a list of 10 objects with 3 observations, each with 5322 entries and 10 columns. 
Abn_For_sum34_max4_N_Mod_pred<-mapply(MyFunc2,Abn_For_sum34_max4_N_Mods_run)%>%as.data.frame()

  #Creates a pdf of all models for a given lichen group. Only used for comparing models. Not for final model plotting.
pdf(paste("Abn_For_sum34_max4_N_Lichen_Abundance_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),".pdf",sep=""))

MyPlot<-function (x)
{
    x<-Abn_For_sum34_max4_N_Mod_pred[i]
    x<-unlist(x, recursive = F, use.names =T)
    fit<-unlist(x[1], recursive = F, use.names =T)
    low<-unlist(x[2], recursive = F, use.names =T)
    high<-unlist(x[3], recursive = F, use.names =T)
    plot(Abn_For_sum34_max4~N, data=SOURCE, xlab="CMAQ N kg/ha/yr", ylab="Abundance Abn_For_sum34_max4 Lichens", main="Abundance of  Abn_For_sum34_max4 lichen species vs N deposition")
    points(SOURCE$N, fit, col='blue', pch= 19)
    points(SOURCE$N, low, col='aquamarine4', pch= 4, cex=0.25) 
    points(SOURCE$N, high, col='aquamarine4', pch= 4, cex=0.25)	
    text(max(SOURCE$N)*0.4,max(SOURCE$Abn_For_sum34_max4)*0.8,labels=paste("R1=",print(round(Abn_For_sum34_max4_N_Mods_fit[i],2))))
    text(max(SOURCE$N)*0.4,max(SOURCE$Abn_For_sum34_max4)*0.75,labels=paste(as.character(Abn_For_sum34_max4_N_Formula[i])), cex=0.7)
    legend('topright', legend=c("Raw Data","Fitted Values", "Prediction Interval"), col=c("black","blue", "aquamarine4"), pch=16)
}	
for (i in 1:length(XVAR)) MyPlot()
dev.off()
  
  #Calculate confidence intervals for percent decline in lichen metric for a bunch of deposition increments
new_n_Abn_For_sum34_max4<-c(1,1.5,2,2.5,3,5,7.5,10,12.5,15,17.5) %>% as.data.frame() #,20
colnames(new_n_Abn_For_sum34_max4)<-"N"
100*round(1-(predict(Abn_For_sum34_max4_N_Mods_run$Abn_For_sum34_max4_N_poly, new_n_Abn_For_sum34_max4))/max(SOURCE_N_fPlot$fit),2)
new_n_CI_Abn_For_sum34_max4<-predict(Abn_For_sum34_max4_N_Mods_run$Abn_For_sum34_max4_N_poly, new_n_Abn_For_sum34_max4, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000)
round(new_n_CI_Abn_For_sum34_max4,2)
stuff<-as.data.frame(new_n_CI_Abn_For_sum34_max4)
cbind(
  round((max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)-(stuff$higher))/max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)*100,0)
  ,round((max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)-(stuff$fit))/max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)*100,0)
  ,round((max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)-(stuff$lower))/max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)*100,0)
)

  ##Calculate incremental percent decline of each lichen index along the 90th quantile fitted line                                                                  
c(
  max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit) ## 0% change
  ,max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)-(max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)*0.05) ## 5% change
  ,max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)-(max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)*0.10) ## 10% change
  ,max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)-(max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)*0.20) ## 20% change
  ,max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)-(max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)*0.50) ## 50% change
  ,max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)-(max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)*0.80) ## 80% change
)

  ## Heuristically solved for X given a Y using the polynomial to find deposition associated with each incremental percent decline
test_n<-as.data.frame(10.4) 
colnames(test_n)<-"N"

  #returns percent decline
1-(predict(Abn_For_sum34_max4_N_Mods_run$Abn_For_sum34_max4_N_poly, test_n))/max(Abn_For_sum34_max4_N_Mod_pred$Abn_For_sum34_max4_N_poly$fit)
  #returns absolute decline
predict(Abn_For_sum34_max4_N_Mods_run$Abn_For_sum34_max4_N_poly, test_high_n) %>% as.matrix()

  #Take deposition associated with incremental percent decline and calculate confidence intervals
n_pct_Abn_For_sum34_max4<-c(0.08,0.52,1.00,1.95,5.30,10.4)%>% as.data.frame()
colnames(n_pct_Abn_For_sum34_max4)<-"N"
n_pct_CI_Abn_For_sum34_max4<-predict(Abn_For_sum34_max4_N_Mods_run$Abn_For_sum34_max4_N_poly, n_pct_Abn_For_sum34_max4, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000)
  
  # print upper and lower bounds of confidence intervals for incremental percent decline
n_pct_CI_Abn_For_sum34_max4

#Use same heuristic solving approach as above to find deposition associated with upper and lower confidence intervals of incremental precent decline
  # 0% # undf - 0.0 - 0.03   #  undf  - 0.08 - 0.385 kg N ha yr  # undf -30.08-29.02
  # 5% #(+0.03) - 0.05 - 0.08 # 0.051 - 0.52 - 0.83 kg N ha yr   # 30.18-28.57-27.51
  #10% # 0.05 - 0.1 - 0.13    # 0.55  - 1.0  - 1.27 kg N ha yr   # 28.45-27.07-26.07
  #20% # 0.16 - 0.2 - 0.22    # 1.55  - 1.95 - 2.18 kg N ha yr   # 25.18-24.06-23.23
  #50% # 0.48 - 0.5 - 0.52    # 5.08  - 5.30 - 5.52 kg N ha yr   # 15.42-15.04-14.40
  #80% # 0.78 - 0.8 - 0.82    # 9.99  - 10.4 - 11.03 kg N ha yr  #  6.49-6.02-5.30

### Scatterplots with fitted lines, confidence intervals and percent decline point estimates with and without legends and
### some with and without one of the axis labels. Uncomment the one image driver and one text label and legend (if desired).
### Label and legend positions may need to be adjusted.

      # Get coeffecients to plot in figures 
      round(coefficients(Abn_For_sum34_max4_N_Mods_run$Abn_For_sum34_max4_N_poly),2)
      # 30.36                -3.51                 0.11 

			  tiff(paste("Abn_For_sum34_max4_N_Forage Lichen_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),".tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
			  #jpeg(paste("Abn_For_sum34_max4_N_Forage Lichen_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),".jpeg",sep=""))
			  #jpeg(paste("Abn_For_sum34_max4_N_Forage Lichen_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.jpeg",sep=""))
			  #tiff(paste("Abn_For_sum34_max4_N_Forage Lichen_vs_N_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
			  
			  x<-Abn_For_sum34_max4_N_Mod_pred[2]
			  x<-unlist(x, recursive = F, use.names =T)
			  y<-Abn_For_sum34_max4_N_Mod_pred[4]
			  clim<-unlist(y[1], recursive = F, use.names =T)
			  clim_pred<-unlist(clim[1], recursive = F, use.names =T)
			  fit<-unlist(x[1], recursive = F, use.names =T)
			  low<-unlist(x[2], recursive = F, use.names =T)
			  high<-unlist(x[3], recursive = F, use.names =T)
			  SOURCE_N_fPlot<-cbind(SOURCE, fit, low, high, clim_pred)
			  SOURCE_N_fPlot<-subset(SOURCE_N_fPlot, N < 12.5)    
			  palette(c("darkorange","darkgreen"))
			  plot(Abn_For_sum34_max4~N, data=SOURCE, xlab=expression(Nitrogen ~ (kg ~ N ~ ha^{-1} ~ y^{-1})), ylab="", col=Area, cex.main = 1.6, cex.lab=1.7, cex.axis = 1.7) # main="90th Quantile Regression of \n Forage Lichen Abundance vs N deposition",
			  title(ylab ="Forage Lichen Abundance", cex.lab=1.7, line=2.5)
			  points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$clim_pred, col=alpha('grey',0.25), pch= 19)
			  points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$fit, col='blue', pch= 19)
			  points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$low, col='black', pch= 4, cex=0.25) 
			  points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$high, col='black', pch= 4, cex=0.25)	
			  #Red dot coordinates are determined by heuristically solving for incremental percent decline (right value) to determine associated deposition (left value)
			  points(0.08,30.08,col='red', pch=20,cex=2) #0% decline
			  #points(0.52,28.57,col='red',pch=20,cex=2) #5% decline
			  points(1.00,27.07,col='red', pch=20,cex=2) #10% decline
			  points(1.95,24.06,col='red', pch=20,cex=2) #20% decline
			  points(5.30,15.04,col='red', pch=20,cex=2) #50% decline
			  points(10.4,6.02 ,col='red', pch=20,cex=2) #80% decline
			  text(max(SOURCE$N)*0.6,max(SOURCE$Abn_For_sum34_max4)*0.6,labels=bquote(atop(.("Forage Lichen Abundance = 30.36 -"),.("3.51*N + 0.11*")*N^2)), cex=1.5)
			  #text(max(SOURCE$N)*0.5,max(SOURCE$Abn_For_sum34_max4)*0.9,labels=bquote(atop(.("Forage Lichen Abundance = 30.36 -"),.("3.51*N + 0.11*")*N^2)), cex=1.5)
			  #legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1)
			  
			  legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1.47)
			  dev.off()
			  
			  
################################################################################################################
#Abundance of Forage Lichens vs sulphur deposition
################################################################################################################

  ## Input files are generated in LichenCL_QuantReg_Script.R script. I have to switch each of these sources depending on which measure of Abundance wanted.
SOURCE<-LichDb_sFncGrpSens_All_Abun_S

  ## Calculate correlations amongst predictors and response
cor_vars<-c("Abn_For_sum34_max4","S","maxaug_c","mindec_c","precip_cm","continen","CMD")
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
VIF(lm(Abn_For_sum34_max4~poly(N, 2, raw=T)+maxaug_c+mindec_c+precip_cm+CMD, data=SOURCE))

  #Give each set of predictors a short name
XVAR_names<-c("S","S_poly","S_clim","S_poly_clim","clim","S_maxaug_c","S_mindec_c","S_precip_cm","S_continen","S_CMD", "Max C Aug Temp","Min C Dec Temp","Precip_cm","Continentality","CMD");

  #Combine each set of predictors with the response, turn it into a formula and fit 90th quantile regression model
Abn_For_S_sum34_max4_Mods<-paste("Abn_For_sum34_max4 ~ ",paste(XVAR),sep="")
Abn_For_S_sum34_max4_Formula<-lapply(Abn_For_S_sum34_max4_Mods,as.formula)
Abn_For_S_sum34_max4_Mods_run<-lapply(Abn_For_S_sum34_max4_Formula,function (x) rq(x, tau=0.9, data=SOURCE))
	
  ## Make null models
Abn_For_S_sum34_max4_S_intercept<-rq(formula = Abn_For_sum34_max4 ~ 1, tau = 0.9, data = SOURCE, model = T)
	
  ## Grab each model fit but the intercept is hard coded 
Abn_For_S_sum34_max4_Mods_fit<-mapply(function (y) 1-Abn_For_S_sum34_max4_Mods_run[[y]]$rho/Abn_For_S_sum34_max4_S_intercept$rho, 1:length(XVAR))
	
  ## makes a list of mod names
Abn_For_S_sum34_max4_Mod_names<-lapply(XVAR_names, function (x) paste("Abn_For_sum34_max4_",paste(x),sep=""))
names(Abn_For_S_sum34_max4_Mods_run)<-Abn_For_S_sum34_max4_Mod_names

  ##Calculate AIC for each model and combine R1 statistic
Abn_For_S_sum34_max4_Mods_AIC<-lapply(Abn_For_S_sum34_max4_Mods_run,AIC)
Abn_For_sum34_max4_LicheS_Abundance_vs_S_Stats<-cbind(XVAR_names,Abn_For_S_sum34_max4_Mods_AIC,Abn_For_S_sum34_max4_Mods_fit,paste(Abn_For_S_sum34_max4_Formula))
  
  ## Save model stats as an output file
write.csv(Abn_For_sum34_max4_LicheS_Abundance_vs_S_Stats,"Abundance_Abn_For_sum34_max4_S_all_models.csv")

	## Make predicted data files
	## Use a model file ... Abn_For_S_sum34_max4_Mods_run as input
MyFunc2<-function (x) {as.data.frame(predict(x, SOURCE, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000))}

	##This produces a matrix with a list of 10 objects with 3 observations, each with 5322 entries and 10 columns. 
Abn_For_S_sum34_max4_Mod_pred<-mapply(MyFunc2,Abn_For_S_sum34_max4_Mods_run)%>%as.data.frame()

  #Creates a pdf of all models for a given lichen group. Only used for comparing models. Not for final model plotting.        	
pdf(paste("Abn_For_sum34_max4_Lichen_Abundance_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),".pdf",sep=""))
MyPlot<-function (x)
  
{
  x<-Abn_For_S_sum34_max4_Mod_pred[i]
  x<-unlist(x, recursive = F, use.names =T)
  fit<-unlist(x[1], recursive = F, use.names =T)
  low<-unlist(x[2], recursive = F, use.names =T)
  high<-unlist(x[3], recursive = F, use.names =T)
  SOURCE_S_fPlot<-cbind(SOURCE, fit, low, high)
  SOURCE_S_fPlot<-subset(SOURCE_S_fPlot, S < 20)
  plot(Abn_For_sum34_max4~S, data=SOURCE, xlab="CMAQ S kg/ha/yr", ylab="Abundance Abn_For_sum34_max4 Lichens", main="Abundance of  Abn_For_sum34_max4 lichen species vs S deposition")
  points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$fit, col='blue', pch= 19)
  points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$low, col='aquamarine4', pch= 4, cex=0.25) 
  points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$high, col='aquamarine4', pch= 4, cex=0.25)	
  text(max(SOURCE$S)*0.4,max(SOURCE$Abn_For_sum34_max4)*0.8,labels=paste("R1=",print(round(Abn_For_S_sum34_max4_Mods_fit[i],2))))
  text(max(SOURCE$S)*0.4,max(SOURCE$Abn_For_sum34_max4)*0.75,labels=paste(as.character(Abn_For_S_sum34_max4_Formula[i])), cex=0.7)
  legend('topright', legend=c("Raw Data","Fitted Values", "Prediction Interval"), col=c("black","blue", "aquamarine4"), pch=16)
}	
for (i in 1:length(XVAR)) MyPlot()
dev.off()

  #Calculate confidence intervals for percent decline in lichen metric for a bunch of deposition increments
new_s_Abn_For_S_sum34_max4<-c(1,1.5,2,2.5,3,5,7.5,10,12.5,15,17.5,20) %>% as.data.frame()
colnames(new_s_Abn_For_S_sum34_max4)<-"S"
100*round(1-(predict(Abn_For_S_sum34_max4_Mods_run$Abn_For_sum34_max4_S_poly,new_s_Abn_For_S_sum34_max4)/max(SOURCE_S_fPlot$fit)),2)
new_s_CI_Abn_For_S_sum34_max4<-predict(Abn_For_S_sum34_max4_Mods_run$Abn_For_sum34_max4_S_poly, new_s_Abn_For_S_sum34_max4, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000)
round(new_s_CI_Abn_For_S_sum34_max4,1)
stuff<-as.data.frame(new_s_CI_Abn_For_S_sum34_max4)
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

  ## Heuristically solved for X given a Y using the polynomial to find deposition associated with each incremental percent decline
test_s<-as.data.frame(1.6)
colnames(test_s)<-"S"
  
  #returns percent decline
1-(predict(Abn_For_S_sum34_max4_Mods_run$Abn_For_sum34_max4_S_poly, test_s))/max(SOURCE_S_fPlot$fit)
  
  #returns absolute decline
predict(Abn_For_S_sum34_max4_Mods_run$Abn_For_sum34_max4_S_poly, test_s)

  #Take deposition associated with incremental percent decline and calculate confidence intervals
s_pct_Abn_For_S_sum34_max4<-c(0.21,0.80,1.40,2.60,6.85,12.70)%>% as.data.frame()
colnames(s_pct_Abn_For_S_sum34_max4)<-"S"
s_pct_CI_Abn_For_S_sum34_max4<-predict(Abn_For_S_sum34_max4_Mods_run$Abn_For_sum34_max4_S_poly, s_pct_Abn_For_S_sum34_max4, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000)
  
  # print upper and lower bounds of confidence intervals for incremental percent decline
s_pct_CI_Abn_For_S_sum34_max4
  #fit     lower    higher
  #1 23.553285 23.135071 24.475433
  #2 22.339469 21.960157 23.172020
  #3 21.138008 20.792688 21.904936
  #4 18.834700 18.506293 19.487683
  #5 11.745356 11.363380 12.214706
  #6  4.711953  4.110287  5.027694

#Use same heuristic solving approach as above to find deposition associated with upper and lower confidence intervals of incremental precent decline
  ## 0 % loss of forage lichen Heuristically solved for undef - 0.21-0.246 S kg/ha/yr given a Y  ...  S kg/ha/yr
  ## 5 % loss of forage lichen Heuristically solved for 0.41  - 0.80 S kg/ha/yr given a Y  ...  S kg/ha/yr
  ## 10% loss of forage lichen Heuristically solved for 1.0   - 1.40 -1.6 S kg/ha/yr given a Y  ...  S kg/ha/yr
  ## 20% loss of forage lichen Heuristically solved for 2.25  - 2.60 -2.78 S kg/ha/yr given a Y  ...  S kg/ha/yr
  ## 50% loss of forage lichen Heuristically solved for 6.53  - 6.85 -7.11 S kg/ha/yr given a Y  ...  S kg/ha/yr
  ## 80% loss of forage lichen Heuristically solved for 12.53 - 12.70-13.36 S kg/ha/yr given a Y  ...  S kg/ha/yr

### Scatterplots with fitted lines, confidence intervals and percent decline point estimates with and without legends and
### some with and without one of the axis labels. Uncomment the one image driver and one text label and legend (if desired).
### Label and legend positions may need to be adjusted.

        # Get coeffecients to plot in figures 
        round(coefficients(Abn_For_S_sum34_max4_Mods_run$Abn_For_sum34_max4_S_poly),2)
        #23.99                -2.10                 0.05 
        
        #tiff(paste("Abn_For_sum34_max4_S_Forage Lichen_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),".tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
        # jpeg(paste("Abn_For_sum34_max4_S_Forage Lichen_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),".jpeg",sep=""))
        # jpeg(paste("Abn_For_sum34_max4_S_Forage Lichen_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.jpeg",sep=""))
        tiff(paste("Abn_For_sum34_max4_S_Forage Lichen_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
        
        # jpeg(paste("Abn_For_sum34_max4_S_Forage Lichen_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_noyaxis_nolegend.jpeg",sep=""))
        #tiff(paste("Abn_For_sum34_max4_S_Forage Lichen_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_noyaxis_nolegend.tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
        
        #tiff(paste("Abn_For_sum34_max4_S_Forage Lichen_vs_S_percentiles","_",format(Sys.time(),"%y%d%m"),"_nolegend.tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
        
        ## This works but is hard coded for each predictor file. Have to change working directory each time
        ## Make final figure with grey dots in the background
        ## Need to add red dots for percentage change
        x<-Abn_For_S_sum34_max4_Mod_pred[2]
        y<-Abn_For_S_sum34_max4_Mod_pred[4]
        clim<-unlist(y[1], recursive = F, use.names =T)
        x<-unlist(x, recursive = F, use.names =T)
        fit<-unlist(x[1], recursive = F, use.names =T)
        low<-unlist(x[2], recursive = F, use.names =T)
        high<-unlist(x[3], recursive = F, use.names =T)
        clim_pred<-unlist(clim[1], recursive = F, use.names =T)
        SOURCE_S_fPlot<-cbind(SOURCE, fit, low, high, clim_pred)
        SOURCE_S_fPlot<-subset(SOURCE_S_fPlot, S < 20)
        palette(c("darkorange","darkgreen"))
        plot(Abn_For_sum34_max4~S, data=SOURCE, xlab=expression(Sulfur ~(kg ~ S ~ ha^{-1} ~ y^{-1})), ylab="",  col=Area, cex.main = 1.6, cex.lab=1.7, cex.axis = 1.7)#main="90th Quantile Regression of \n Forage Lichen Abundance vs S deposition",
        title(ylab ="Forage Lichen Abundance", cex.lab=1.7, line=2.5)
        points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$clim_pred, col=alpha('grey',0.25), pch= 19)
        points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$fit, col='blue', pch= 19)
        points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$low, col='black', pch= 4, cex=0.25) 
        points(SOURCE_S_fPlot$S, SOURCE_S_fPlot$high, col='black', pch= 4, cex=0.25)	
        #Red dot coordinates are determined by heuristically solving for incremental percent decline (right value) to determine associated deposition (left value)
        points(0.21, 23.553285,col='red',  pch=20,cex=2) #0% decline
        #points(0.80, 22.37562 ,col='red', pch=20,cex=2) #5% decline
        points(1.40, 21.197956 ,col='red', pch=20,cex=2) #10% decline
        points(2.60, 18.842628 ,col='red', pch=20,cex=2) #20% decline
        points(6.85, 11.776642 ,col='red', pch=20,cex=2) #50% decline
        points(12.70,4.710657 ,col='red',  pch=20,cex=2) #80% decline
        
        #text(max(SOURCE$S)*0.5,max(SOURCE$Abn_For_sum34_max4)*0.6,labels=bquote(atop(.("Forage Lichen Abundance = 23.99 -"),.("2.10*S + 0.05*")*S^2)), cex=1.5)
        text(max(SOURCE$S)*0.5,max(SOURCE$Abn_For_sum34_max4)*0.9,labels=bquote(atop(.("Forage Lichen Abundance = 23.99 -"),.("2.10*S + 0.05*")*S^2)), cex=1.5)	  
        #legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1.47)
        #legend('topright', legend=c("Raw Data East","Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkorange","darkgreen","blue", "black","grey","red"), pch=16, cex=1)
        dev.off()
