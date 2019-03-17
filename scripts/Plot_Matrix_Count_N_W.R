############################################################################################################################################
#Western USA - Count of Oligotroph vs Nitrogen 
#################################################################################################################

	## Input files are generated in Calculate Counts script. I have to switch each of these sources depending on which measure of Count wanted.
	## We decided to use the maximum Count because it was the most straightforward to explain with little to no loss in model fit.
## Instead of changing SOURCE, I should build a separate pipeline for each measure of Count (summing 3s and 4s, summing log transformed values and only max value)
	SOURCE<-LichDb_sFncGrpSens_All_Abun_N_West

XVAR<-c("N"                                                               
 ,"poly(N, 2)"
 ,"poly(N, 3)" 
 ,"poly(N, 4)" 
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
 ,"CMD"
 ) 

XVAR_names<-c("N","N_poly","N_poly3","N_poly4","N_clim","N_poly_clim","clim","N_maxaug_c","N_mindec_c","N_precip_cm","N_continen","N_CMD", "Max C Aug Temp","Min C Dec Temp","Precip_cm","Continentality","CMD");

N_Matrix_N_Mods<-paste("N_Matrix ~ ",paste(XVAR),sep="")

N_Matrix_N_Formula<-lapply(N_Matrix_N_Mods,as.formula)

N_Matrix_N_Mods_run_West<-lapply(N_Matrix_N_Formula,function (x) rq(x, tau=0.9, data=SOURCE))

	## Make null models
N_Matrix_N_intercept_West<-rq(formula = N_Matrix ~ 1, tau = 0.9, data = SOURCE, model = T)

	## Grab each model fit but the intercept is hard coded 

N_Matrix_N_Mods_fit_West<-mapply(function (y) 1-N_Matrix_N_Mods_run_West[[y]]$rho/N_Matrix_N_intercept_West$rho, 1:length(XVAR))

	## makes a list of mod names

N_Matrix_N_Mod_names<-lapply(XVAR_names, function (x) paste("N_Matrix_",paste(x),sep=""))

names(N_Matrix_N_Mods_run_West)<-N_Matrix_N_Mod_names

	## Save model stats as a spreadsheet of manuscript

N_Matrix_N_Mods_AIC<-lapply(N_Matrix_N_Mods_run_West,AIC)

N_Matrix_N_Lichen_Count_vs_N_Stats<-cbind(XVAR_names,N_Matrix_N_Mods_AIC,N_Matrix_N_Mods_fit_West,paste(N_Matrix_N_Formula))
write.csv(N_Matrix_N_Lichen_Count_vs_N_Stats,"Count_N_Matrix_N_all_models.csv")

	## Make predicted data files
	## Use a model file ... N_Matrix_Mods_run as input

MyFunc2<-function (x) {as.data.frame(predict(x, SOURCE, interval=c("confidence"), level = 0.95, type='percentile', se="boot", N=10000))}

	##This produces a matrix with a list of 10 objects with 3 observations, each with 5322 entries and 10 columns. 
N_Matrix_N_Mod_pred_West<-mapply(MyFunc2,N_Matrix_N_Mods_run_West)%>%as.data.frame()

	## This works but is hard coded for each predictor file. Have to change working directory each time

pdf(paste("N_Matrix_N_Lichen_Count_vs_N_percentiles_West","_",format(Sys.time(),"%y%d%m"),".pdf",sep=""))

MyPlot<-function (x)
    
{
    x<-N_Matrix_N_Mod_pred_West[i]
    x<-unlist(x, recursive = F, use.names =T)
    fit<-unlist(x[1], recursive = F, use.names =T)
    low<-unlist(x[2], recursive = F, use.names =T)
    high<-unlist(x[3], recursive = F, use.names =T)
    #palette(c("darkgreen","darkorange"))
    plot(N_Matrix~N, data=SOURCE, xlab="CMAQ N kg/ha/yr", ylab="Count Matrix Lichens", main="Count of Matrix lichen species \n vs N deposition in the western US", col="darkgreen")
    points(SOURCE$N, fit, col='blue', pch= 19)
    points(SOURCE$N, low, col='black', pch= 4, cex=0.25) 
    points(SOURCE$N, high, col='black', pch= 4, cex=0.25)	
		#points(0.08,5.63, col='red', pch=20)
		#points(0.44,5.35, col='red', pch=20)
		#points(0.80,5.07, col='red', pch=20)
		#points(1.58,4.51, col='red', pch=20)
		#points(4.27,2.82, col='red', pch=20)
		#points(7.95,1.13, col='red', pch=20)
    text(max(SOURCE$N)*0.7,max(SOURCE$N_Matrix)*0.35,labels=paste("R1=",print(round(N_Matrix_N_Mods_fit_West[i],2))))
    text(max(SOURCE$N)*0.7,max(SOURCE$N_Matrix)*0.3,labels=paste(as.character(N_Matrix_N_Formula[i])), cex=0.7)
    legend('bottomright', legend=c("Raw Data East","Fitted Values", "Prediction Interval"), col=c("darkorange","blue", "black"), pch=16)
 
}	

for (i in 1:length(XVAR)) MyPlot()

dev.off()
    #tiff(paste("N_Matrix_N_Lichen_Count_vs_N_percentiles_W_US_","_",format(Sys.time(),"%y%d%m"),".tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)
	  #jpeg(paste("N_Matrix_N_Lichen_Count_vs_N_percentiles_W_US_",format(Sys.time(),"%y%d%m"),".jpeg",sep=""))
	  jpeg(paste("N_Matrix_N_Lichen_Count_vs_N_percentiles_W_US_",format(Sys.time(),"%y%d%m"),"_nolegend.jpeg",sep=""))
    #tiff(paste("N_Matrix_N_Lichen_Count_vs_N_percentiles_W_US_","_",format(Sys.time(),"%y%d%m"),"_nolegend.tiff",sep=""), width=1200,  height=1200, units="px", pointsize = 24)

    x<-N_Matrix_N_Mod_pred_West[2]
    x2<-unlist(x, recursive = F, use.names =T)
	  y<-N_Matrix_N_Mod_pred_West[4]
	  clim<-unlist(y[1], recursive = F, use.names =T)
	  clim_pred<-unlist(clim[1], recursive = F, use.names =T)

    fit<-unlist(x2[1], recursive = F, use.names =T)
    low<-unlist(x2[2], recursive = F, use.names =T)
    high<-unlist(x2[3], recursive = F, use.names =T)
    SOURCE_N_fPlot<-cbind(SOURCE, fit, low, high, clim_pred)
	  SOURCE_N_fPlot<-subset(SOURCE_N_fPlot, N < 12.5)    
  	palette(c("darkorange","darkgreen"))
    plot(N_Matrix~N, data=SOURCE, xlab="CMAQ N kg/ha/yr", ylab="", main="90th Quantile Regresssion Count of \n Western Oligotroph Species vs N deposition", col=Area, cex.main = 1.6, cex.lab=1.7, cex.axis = 1.7)
    title(ylab ="Count of Oligotrophs", cex.lab=1.7, line=2.5)
    points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$clim_pred, col=alpha('grey',0.25), pch= 19)
    points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$fit, col='blue', pch= 19)
    points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$low, col='black', pch= 4, cex=0.25) 
    points(SOURCE_N_fPlot$N, SOURCE_N_fPlot$high, col='black', pch= 4, cex=0.25)	
		#points(0.08,	20.30, col='red', pch=20,cex=2)
		##points(0.46,	19.29, col='red', pch=20,cex=2)
		#points(0.85,	18.27, col='red', pch=20,cex=2)
		#points(1.67,	16.24, col='red', pch=20,cex=2)
		#points(4.60,	10.15, col='red', pch=20,cex=2)
		#points(9.45,	4.06, col='red',  pch=20,cex=2)
    #text(max(SOURCE$N)*0.8,max(SOURCE$N_Matrix)*0.7,labels=paste("R1=",print(round(N_Matrix_N_Mods_fit_West[2],2))))
    #text(max(SOURCE$N)*0.8,max(SOURCE$N_Matrix)*0.75,labels=paste(as.character(N_Matrix_N_Formula[2])), cex=1)
    #text(max(SOURCE$N)*0.8,max(SOURCE$N_Matrix)*0.65,labels=paste("AIC=",print(round(as.numeric(N_Matrix_N_Mods_AIC[2],2)))))
	
	#text(max(SOURCE$N)*0.65,max(SOURCE$N_Matrix)*0.5,labels=bquote(atop(.("Count of Oligotrophs = 13.81 -"),.("176.53*N + 60.61*")*N^2)), cex=1.5)
	#text(max(SOURCE$N)*0.65,max(SOURCE$N_Matrix)*0.9,labels=bquote(atop(.("Count of Oligotrophs = 13.81 -"),.("176.53*N + 60.61*")*N^2)), cex=1.5)
	#legend('topright', legend=c("Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence\n Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkgreen","blue", "black","grey","red"), pch=16, cex=1)
	#legend('topright', legend=c("Raw Data West", "Fitted Values Deposition only", "95% Bootstrapped Confidence\n Interval Deposition only","Fitted Values Climate + Deposition","0/10/20/50/80% loss"), col=c("darkgreen","blue", "black","grey","red"), pch=16, cex=1.47)
	dev.off()

###### Make final version of figure with grey cloud of points in the background of the climate fitted values, the raw data as open black circles
###### and the fitted polynomial line with the predicted fitted 

	five_pct_N_Matrix  <- max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)-(max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)*0.05) %>% as.data.frame()
	ten_pct_N_Matrix   <- max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)-(max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)*0.1) %>% as.data.frame()
	twenty_pct_N_Matrix<- max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)-(max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)*0.2) %>% as.data.frame()
	fifty_pct_N_Matrix <- max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)-(max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)*0.5) %>% as.data.frame()
	eighty_pct_N_Matrix<- max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)-(max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)*0.8) %>% as.data.frame()
	
	colnames(five_pct_N_Matrix)<-"N"
	colnames(ten_pct_N_Matrix)<-"N"   
	colnames(twenty_pct_N_Matrix)<-"N"
	colnames(fifty_pct_N_Matrix)<-"N"
	colnames(eighty_pct_N_Matrix)<-"N"

	
	as.data.frame(c(five_pct_N_Matrix  
	,ten_pct_N_Matrix   
	,twenty_pct_N_Matrix
	,fifty_pct_N_Matrix 
	,eighty_pct_N_Matrix))
		
		## x - 5
		## y - N_Matrix_N_Mods_run_West
		## ModRun - N_Matrix_N_poly
		## ModPred N_Matrix_N_Mod_pred_West
			
	#GetYforX<-function (x,y,ModRun,ModPred)
	#{		
	#test_n<-as.data.frame(x)
	#colnames(test_n)<-"N"
	#(predict(y$ModRun, test_n))/max(ModPred$ModRun)
	#}
	#GetYforX(5,N_Matrix_N_Mods_run_West,N_Matrix_N_poly,N_Matrix_N_Mod_pred_West)
			
			test_n<-as.data.frame(5)
			colnames(test_n)<-"N"
			1-(predict(N_Matrix_N_Mods_run_West$N_Matrix_N_poly, test_n)/max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit))
		
		## Heuristically solved for X given a Y using the polynomial 
			## Heuristically solved for X given a Y ...	34.58681
				## 0 % (max N_Matrix ) Heuristically solved for X given a Y ... 0.08 N kg/ha/yr
				## 5 % loss of N. Cya Heuristically solved for X given a Y ...	0.47 N kg/ha/yr
				## 10% loss of N. Cya Heuristically solved for X given a Y ...	 0.8 N kg/ha/yr
				## 20% loss of N. Cya Heuristically solved for X given a Y ...  1.58 N kg/ha/yr
				## 50% loss of N. Cya Heuristically solved for X given a Y ...  4.23 N kg/ha/yr
				## 80% loss of N. Cya Heuristically solved for X given a Y ...  7.9 N kg/ha/yr

				new_n_1kg	<-1.0
				new_n_1.5kg	<-1.5
				new_n_2kg	<-2.0
				new_n_2.5kg	<-2.5
				new_n_3kg	<-3.0
				new_n_5kg   <-5.0
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
				
					#new_n_1kg$N	  <-as.numeric(new_n_1kg$N)
					#new_n_1.5kg$N  <-as.numeric(new_n_1.5kg$N)
					#new_n_2kg$N	  <-as.numeric(new_n_2kg$N)
					#new_n_2.5kg$N <-as.numeric(new_n_2.5kg$N)
					#new_n_3kg$N	  <-as.numeric(new_n_3kg$N)
					#new_n_5kg$N   <-as.numeric(new_n_5kg$N)

				
			

				c(
				 1-(predict(N_Matrix_N_Mods_run_West$N_Matrix_N_poly, new_n_1kg))/max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)
				,1-(predict(N_Matrix_N_Mods_run_West$N_Matrix_N_poly,new_n_1.5kg))/max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)
				,1-(predict(N_Matrix_N_Mods_run_West$N_Matrix_N_poly,new_n_2kg))/max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)
				,1-(predict(N_Matrix_N_Mods_run_West$N_Matrix_N_poly,new_n_2.5kg))/max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)
				,1-(predict(N_Matrix_N_Mods_run_West$N_Matrix_N_poly,new_n_3kg))/max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)
				,1-(predict(N_Matrix_N_Mods_run_West$N_Matrix_N_poly,new_n_5kg))/max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)
				)
				
				c(
				max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit) ## 0% change
				,max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)-(max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)*0.05) ## 5% change
				,max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)-(max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)*0.10) ## 10% change
				,max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)-(max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)*0.20) ## 20% change
				,max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)-(max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)*0.50) ## 50% change
				,max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)-(max(N_Matrix_N_Mod_pred_West$N_Matrix_N_poly$fit)*0.80) ## 80% change
				)
				
	
