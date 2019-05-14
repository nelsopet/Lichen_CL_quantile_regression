##Final script for 90th quantile regression analyses of lichen community measures for the US based on USFS plots.
##Abundance:filter by lichen index(eg.cyanolichens),delete1’sand2’s,sumthe3’sand4’s by lichen index (eg.cyanolichen)

## Set working directory wherever you've put the files to be written
##setwd('/Users/peternelson1/Documents/UMFK/ORAU/FilesforPeter/Analysis/')

## Install required packages
#Function taken from https://nhsrcommunity.com/blog/a-simple-function-to-install-and-load-packages-in-r/
install_or_load_pack <- function(pack){
  
  create.pkg <- pack[!(pack %in% installed.packages()[, "Package"])]
  
  if (length(create.pkg))
    
    install.packages(create.pkg, dependencies = TRUE)
  
  sapply(pack, require, character.only = TRUE)
  
  #I know I should be using purr here, but this is before the Tidyverse is loaded. I know you Tidyverse trend setters will have me here.
  
}

packs<-c("plyr",
         "quantreg",
         "tidyverse",
         "scales",
         "vegan",
         "fmsb")

install_or_load_pack(packs)

## Read in data
## Lichen community data; 127001 rows of individual lichen observations, three measures of abundance,
## many columns of pollutant data through time, some plot location information. 
LichDb_raw<-read.csv(file="LichenDB9.csv",header=T)

# Climate data; 5322 rows of output from ClimateNA which are 1960-1999 normals)
ClimateNA<-read.csv(file="ForClimateNa_Normal_1961_1990Y.csv",header=T)

# 4356 rows; Five lichen community indices (N_sens, N_oligo, N_For, N_Cyan, spp_rich) and plot number
# Why is N_sens in here? That wouldn't apply to N dep. Was this included inadvertently?
Lichen_Indices_Nmods_US_toR<-read.csv(file="Lichen_Indices_Nmods_US_toR.csv", header=T)

# 5172 rows; Five lichen community indices (N_sens, N_oligo, N_For, N_Cyan, spp_rich) and plot number
# Why is N_oligo in here? That wouldn't apply to S dep. Was this included inadvertently?
Lichen_Indices_Smods_US_toR<-read.csv(file="Lichen_Indices_Smods_US_toR.csv", header=T)

# 849 rows; Four variables of Sensitivity, Pollutant, Area and Functional group with SciCode and Name
# This is from the national lichen critical load manuscript
LichFncGrp<-read.csv(file="LichenFunctionalGroups.csv",header=T)
LichFncGrp$'Functional.Group'<-as.character(LichFncGrp$'Functional.Group')
LichFncGrp$'Functional.Group'[LichFncGrp$'Functional.Group'=="hv matrix"]<-"matrix"
LichFncGrp$'Functional.Group'[LichFncGrp$'Functional.Group'=="lv matrix"]<-"matrix"
LichFncGrp$'Functional.Group'[LichFncGrp$'Functional.Group'=="small matrix"]<-"matrix"
LichFncGrp$'Functional.Group'[LichFncGrp$'Functional.Group'=="large matrix"]<-"matrix"
LichFncGrp$'Functional.Group'<-as.factor(LichFncGrp$'Functional.Group')


# 804 rows; Supplmentary table from manuscript that has various species-level data including finer levels of functional
# groups 
LichFncGrp_ms<-read.csv(file="Supplementary Table 2.csv",header=T)

# 5322 rows; Fifty seven columns, all climate variables except CMD are here; most of the functional group measures
# are also in this file as is deposition. Why then do I need the other climate files?
# This file came from the original spreadsheet Linda had been using for comparing the different indices. I took one
# tab and put it in its own spreadsheet for use in R.
NS_vs_Indices<-read.csv(file="N&SscoresVsLichenIndicesDB.csv",header=T)

## In NS_vs_Indices, make megadbid a factor so it can be joined to other tables
#NS_vs_Indices$megadbid<-as.factor(NS_vs_Indices$megadbid)
NS_vs_Indices$megadbid<-as.integer(NS_vs_Indices$megadbid)

## Add counts of larger functional groupings
NS_vs_Indices<-NS_vs_Indices %>% 
  mutate(N_Cya = N.Cya.l. + N.Cya.s.m., N_For = N.For.p. + N.For.sh., N_Matrix = N.Mtx.m.l. + N.Mtx.s.)

## For Nitrogen 
NS_vs_Indices_N_CMD_only<-inner_join(Lichen_Indices_Nmods_US_toR,ClimateNA, by= "Plot") %>% select(Plot, CMD)

## For Sulphur
NS_vs_Indices_S_CMD_only<-inner_join(Lichen_Indices_Smods_US_toR,ClimateNA, by= "Plot") %>% select(Plot, CMD)

## Make plots factors
#NS_vs_Indices_N_CMD_only$Plot<-as.factor(NS_vs_Indices_N_CMD$Plot);
#NS_vs_Indices_S_CMD_only$Plot<-as.factor(NS_vs_Indices_S_CMD$Plot);
NS_vs_Indices_N_CMD_only$Plot<-as.integer(NS_vs_Indices_N_CMD$Plot);
NS_vs_Indices_S_CMD_only$Plot<-as.integer(NS_vs_Indices_S_CMD$Plot);

## Change CMD values that are below 0 to zero
NS_vs_Indices_N_CMD_only$CMD[NS_vs_Indices_N_CMD_only$CMD<0]<-0
NS_vs_Indices_S_CMD_only$CMD[NS_vs_Indices_S_CMD_only$CMD<0]<-0

## In raw lichen community data, change NAs to 3’s and abundances >4 back to 3’s
#LichDb_raw$megadbid<-as.factor(LichDb_raw$megadbid)
LichDb_raw$megadbid<-as.integer(LichDb_raw$megadbid)
LichDb_raw$maxabun[is.na(LichDb_raw$maxabun)]<-3
LichDb_raw$maxabun[LichDb_raw$maxabun>4]<-3

##Calculate the “orginal” form of abundance, where we turn all 1’s and 2’s to 0 and sum the 3’s and 4’s
LichDb_raw$maxabun_34<-LichDb_raw$maxabun
LichDb_raw$maxabun_34[LichDb_raw$maxabun_34<3]<-0

##Calculate the second form of abundance, rescaled and summed without deleting any species
LichDb_raw$maxabun_log<-LichDb_raw$maxabun
LichDb_raw$maxabun_log[LichDb_raw$maxabun_log==1]<-0.001
LichDb_raw$maxabun_log[LichDb_raw$maxabun_log==2]<-0.01
LichDb_raw$maxabun_log[LichDb_raw$maxabun_log==3]<-0.1
LichDb_raw$maxabun_log[LichDb_raw$maxabun_log==4]<-1

##Make East and West categories consistent
LichDb_raw$us_loc<-as.character(LichDb_raw$us_loc)
LichDb_raw$us_loc[LichDb_raw$us_loc=="east"]<-"East"
LichDb_raw$us_loc[LichDb_raw$us_loc=="west"]<-"West"

##Calculate the third form of abundance, which is just the maximum abundance for each functional group in the latter section of this script
##with the pipeline where the tables are joined and grouped.
LichFncGrp_ms_unique<-LichFncGrp_ms%>%
  select(Code,Scientific.Name,Fxl.Gp)%>%
  distinct()%>%
  arrange(Code)

##The group_by section wasn’t working until I removed plyr...see thread below
##https://stackoverflow.com/questions/21653295/dplyr-issues-when-using-group-bymultiple-variables		

##Functional group abundance calculation summing 3's and 4's with a max abundance of 4
##Need to add summing 3's and 4's in subcategories of Functional groups which are then
##added together. See oringal abundance calculation script for that code, which
##I considered a mistake originally but it turned out to be a better fitting response
##combining both abundance and diversity more than summing 3's and 4's

## Combine the two different levels of functional graoups and give them different names	
Fxl.Gp_unique<-levels(LichFncGrp_ms_unique$Fxl.Gp)
Fxl.Gp_Big<-c('Cya','Cya','For','For','Mtx','Mtx')
Fxl.Gps<-cbind(Fxl.Gp_unique,Fxl.Gp_Big)%>%as.data.frame()

###############################

##Calculate abundance(sum3'sand4's) of functional groups
LichDb_wFncGrp_sum34_max8_noNA<-LichDb_raw%>%
subset(select=c(megadbid,newscicode,maxabun_34))%>%
inner_join(LichFncGrp_ms_unique,by=c("newscicode"="Code"))%>%
inner_join(Fxl.Gps,by=c("Fxl.Gp"="Fxl.Gp_unique"))%>%
group_by(.dots=c("megadbid","Fxl.Gp"))%>%
summarise(TotAbun=max(maxabun_34))%>%
as.data.frame()%>%
spread(Fxl.Gp,TotAbun)%>%
as.data.frame()

LichDb_wFncGrp_sum34_max8_noNA[is.na(LichDb_wFncGrp_sum34_max8_noNA)]<-0

colnames(LichDb_wFncGrp_sum34_max8_noNA)<-c("megadbid","Abn_Cya_l","Abn_Cya_s_m","Abn_For_p","Abn_For_sh","Abn_Mtx_m_l","Abn_Mtx_s")		

##Calculate LichDb_wFncGrp_log_final created in Calculate Abundances script.
LichDb_wFncGrp_sum34_max8_noNA$Abn_For_sum34_max8<-LichDb_wFncGrp_sum34_max8_noNA$"Abn_For_p"+LichDb_wFncGrp_sum34_max8_noNA$"Abn_For_sh"
LichDb_wFncGrp_sum34_max8_noNA$Abn_Cya_sum34_max8<-LichDb_wFncGrp_sum34_max8_noNA$"Abn_Cya_l"+LichDb_wFncGrp_sum34_max8_noNA$"Abn_Cya_s_m"
LichDb_wFncGrp_sum34_max8_noNA$Abn_Mtx_sum34_max8<-LichDb_wFncGrp_sum34_max8_noNA$"Abn_Mtx_m_l"+LichDb_wFncGrp_sum34_max8_noNA$"Abn_Mtx_s"

LichDb_wSens_sum34_noNA<-LichDb_raw%>%
subset(select=c(megadbid,newscicode,maxabun_34))%>%
inner_join(LichFncGrp_ms,by=c("newscicode"="Code"))%>%
group_by(.dots=c("megadbid","Sensitivity"))%>%
summarise(TotAbun=sum(maxabun_34))%>%
as.data.frame()%>%
spread(Sensitivity,TotAbun)%>%
as.data.frame()

LichDb_wSens_sum34_noNA[is.na(LichDb_wSens_sum34_noNA)]<-0
colnames(LichDb_wSens_sum34_noNA)<-c("megadbid","Abn_eutr","Abn_inter","Abn_meso","Abn_olig","Abn_sens","Abn_tol")

LichDb_wFncGrp_max_noNA<-LichDb_raw%>%
subset(select=c(megadbid,newscicode,maxabun_34))%>%
inner_join(LichFncGrp_ms_unique,by=c("newscicode"="Code"))%>%
inner_join(Fxl.Gps,by=c("Fxl.Gp"="Fxl.Gp_unique"))%>%
group_by(.dots=c("megadbid","Fxl.Gp_Big"))%>%
summarise(TotAbun=max(maxabun_34))%>%
as.data.frame()%>%
spread(Fxl.Gp_Big,TotAbun)				

LichDb_wFncGrp_max_noNA[is.na(LichDb_wFncGrp_max_noNA)]<-0
colnames(LichDb_wFncGrp_max_noNA)<-c("megadbid","Abn_Max_Cya","Abn_Max_For","Abn_Max_Mtx")

LichDb_wSens_max_noNA<-LichDb_raw%>%
subset(select=c(megadbid,newscicode,maxabun_34))%>%
inner_join(LichFncGrp_ms,by=c("newscicode"="Code"))%>%
group_by(.dots=c("megadbid","Sensitivity"))%>%
summarise(TotAbun=max(maxabun_34))%>%
as.data.frame()%>%
spread(Sensitivity,TotAbun)%>%
as.data.frame()

LichDb_wSens_max_noNA[is.na(LichDb_wSens_max_noNA)]<-0
colnames(LichDb_wSens_max_noNA)<-c("megadbid","Abn_Max_eutr","Abn_Max_inter","Abn_Max_meso","Abn_Max_olig","Abn_Max_sens","Abn_Max_tol")

LichDb_wFncGrp_sumlog_noNA<-LichDb_raw%>%subset(select=c(megadbid,newscicode,maxabun_log))%>%
inner_join(LichFncGrp_ms_unique,by=c("newscicode"="Code"))%>%
inner_join(Fxl.Gps,by=c("Fxl.Gp"="Fxl.Gp_unique"))%>%
group_by(.dots=c("megadbid","Fxl.Gp_Big"))%>%
summarise(TotAbun=sum(maxabun_log))%>%
as.data.frame()%>%
spread(Fxl.Gp_Big,TotAbun)	
		
LichDb_wFncGrp_sumlog_noNA[is.na(LichDb_wFncGrp_sumlog_noNA)]<-0
colnames(LichDb_wFncGrp_sumlog_noNA)<-c("megadbid","Abn_Log_Cya","Abn_Log_For","Abn_Log_Mtx")

LichDb_wSens_sumlog_noNA<-LichDb_raw%>%
subset(select=c(megadbid,newscicode,maxabun_log))%>%
inner_join(LichFncGrp_ms,by=c("newscicode"="Code"))%>%
group_by(.dots=c("megadbid","Sensitivity"))%>%
summarise(TotAbun=max(maxabun_log))%>%
as.data.frame()%>%
spread(Sensitivity,TotAbun)%>%
as.data.frame()

LichDb_wSens_sumlog_noNA[is.na(LichDb_wSens_sumlog_noNA)]<-0
colnames(LichDb_wSens_sumlog_noNA)<-c("megadbid","Abn_Log_eutr","Abn_Log_inter","Abn_Log_meso","Abn_Log_olig","Abn_Log_sens","Abn_Log_tol")

LichDb_wFncGrp_sum34_max4_noNA<-LichDb_raw%>% 
subset(select=c(megadbid,newscicode,maxabun_34))%>%
inner_join(LichFncGrp_ms_unique,by=c("newscicode"="Code"))%>%
inner_join(Fxl.Gps,by=c("Fxl.Gp"="Fxl.Gp_unique"))%>%
group_by(.dots=c("megadbid","Fxl.Gp_Big"))%>%
summarise(TotAbun=sum(maxabun_34))%>%
as.data.frame()%>%
spread(Fxl.Gp_Big,TotAbun)

LichDb_wFncGrp_sum34_max4_noNA[is.na(LichDb_wFncGrp_sum34_max4_noNA)]<-0
colnames(LichDb_wFncGrp_sum34_max4_noNA)<-c("megadbid","Abn_Cya_sum34_max4","Abn_For_sum34_max4","Abn_Mtx_sum34_max4")

LichDb_sFncGrpSens_All_Abun<-join_all(list(LichDb_wFncGrp_sum34_max8_noNA,
LichDb_wSens_sum34_noNA,
LichDb_wFncGrp_max_noNA,
LichDb_wSens_max_noNA,
LichDb_wFncGrp_sumlog_noNA,
LichDb_wSens_sumlog_noNA,
LichDb_wFncGrp_sum34_max4_noNA),by='megadbid',type='left')%>% as.data.frame()

LichDb_sFncGrpSens_All_Abun$megadbid<-as.integer(LichDb_sFncGrpSens_All_Abun$megadbid)
#LichDb_sFncGrpSens_All_Abun$megadbid<-as.factor(LichDb_sFncGrpSens_All_Abun$megadbid)

## write csv file for use outside R
LichDb_sFncGrpSens_All_Abun<-write.csv(LichDb_sFncGrpSens_All_Abun, "LichDb_sFncGrpSens_All_Abun.csv")
##############################

##Nitrogen; 4356 rows; Ninety four variables
LichDb_sFncGrpSens_All_Abun_N<-NS_vs_Indices%>%
inner_join(NS_vs_Indices_N_CMD_only,by=c("megadbid"="Plot"))%>% 
inner_join(LichDb_sFncGrpSens_All_Abun,by=c("megadbid"))

##East and West LichDb_sFncGrpSens_All_Abun for Nitrogen
LichDb_sFncGrpSens_All_Abun_N_East<-subset(LichDb_sFncGrpSens_All_Abun_N,Area=="East")
LichDb_sFncGrpSens_All_Abun_N_West<-subset(LichDb_sFncGrpSens_All_Abun_N,Area=="West")

##Sulphur; 5172 rows; Ninety four variables
LichDb_sFncGrpSens_All_Abun_S<-NS_vs_Indices%>%
inner_join(NS_vs_Indices_S_CMD_only,by=c("megadbid"="Plot"))%>%
inner_join(LichDb_sFncGrpSens_All_Abun,by=c("megadbid"))

##Remove Nitrogen to avoid using wrong dataset
#LichDb_sFncGrpSens_All_Abun_S<-LichDb_sFncGrpSens_All_Abun_S[-LichDb_sFncGrpSens_All_Abun_S$N]

##East and West LichDb_sFncGrpSens_All_Abun
LichDb_sFncGrpSens_All_Abun_S_East<-subset(LichDb_sFncGrpSens_All_Abun_S,Area=="East")
LichDb_sFncGrpSens_All_Abun_S_West<-subset(LichDb_sFncGrpSens_All_Abun_S,Area=="West")