#Set working directory
#Make sure you change this to the folder than contains all of the data, input and function files
setwd("") #Set this to the folder that contains the following three files
source('stochasticLBSPR.r')
source('YPR_SPR_RPs.r')
source('Flexible_Control_Rule_code.r')

###################################################################
### Estimating uncertainty in LBSPR through life history values ###
###################################################################
setwd("") #Reset directory to folder with input and life history files (if different than directory set above)

#Isolating March 2018 as the best sampled month
#Option 1: Use life history input file
Ex.Spp.LBSPR.LHin<-LBSPR_stochastic(
  Length.dat="Ex.Spp.csv", #.csv file of lengths
  dataType.in="raw", #or raw
  Nsamp=1000,
  LH_in=1, 
  LH.file="Ex.Spp_LH.csv",
  Linf=NA,
  Linf_CV=0.1,
  k=NA,
  k_CV=0.1,
  Cor_Linf_k=-0.9,
  M_k=NA,
  M_k_CV=0.2,
  Cor_Linf_M_k=0.2,
  M_prior=NA,
  L50=NA,
  L95=NA,
  maxF=10,
  F_step = 0.1)

#Option 2: Establish Linf, k and M distributions, then get M/k
#Use L,K, and M inputs based on Stephanus estimates
load("M_Ex.Spp.DMP")
M.in<-pdf.samples
Ex.Spp.LBSPR.LinfkM_MK<-LBSPR_stochastic(
  Length.dat="Ex.Spp.csv", #.csv file of lengths
  dataType.in="raw", #or raw
  Nsamp=1000,
  LH_in=2, 
  LH.file=NA,
  Linf=10.5,
  Linf_CV=0.1,
  k=0.77,
  k_CV=0.1,
  Cor_Linf_k=-0.9,
  M_prior=M.in,
  M_k=NA,
  M_k_CV=0.0001,
  Cor_Linf_M_k=0.2,
  L50=4,
  L95=6,
  maxF=15,
  F_step = 0.1)

#Option 3: Use L and M/k inputs
#Based on Stephanus life history estimates
Ex.Spp.LBSPR.LinfMK<-LBSPR_stochastic(
  Length.dat="Ex.Spp.csv", #.csv file of lengths
  dataType.in="raw", #or raw
  Nsamp=1000,
  LH_in=3, 
  LH.file=NA,
  Linf=10.5,
  Linf_CV=0.1,
  k=NA,
  k_CV=0.1,
  Cor_Linf_k=-0.9,
  M_prior=NA,
  M_k=2.59,
  M_k_CV=0.2,
  Cor_Linf_M_k=0.2,
  L50=4,
  L95=6,
  maxF=15,
  F_step = 0.1)

#Look at the median and 95% CI of selected values (in this case, SPR and Fishing rate F)
apply(Ex.Spp.LBSPR.LHin$SPR,2,quantile,c(0.025,0.5,0.975),na.rm=T)
apply(Ex.Spp.LBSPR.LHin$Sel50,2,quantile,c(0.025,0.5,0.975),na.rm=T)
apply(Ex.Spp.LBSPR.LHin$Sel95,2,quantile,c(0.025,0.5,0.975),na.rm=T)

###############################
### Plot outputs in barplot ###
###############################
#Prepare outputs for plotting in ggplot
LBSPRs.out.Ex.Spp<-data.frame(cbind(Ex.Spp.LBSPR.LHin$SPR,Ex.Spp.LBSPR.LinfkM_MK$SPR,Ex.Spp.LBSPR.LinfMK$SPR))
colnames(LBSPRs.out.Ex.Spp)<-c("LHrows","Linf_k_M","Lin_MK")
LBSPRs.out.df.gg<-melt(LBSPRs.out.Ex.Spp,variable.name="Data",value.name="SPR")
colnames(LBSPRs.out.df.gg)<-c("Scenario","SPR")
#Combine the 5 scenarios as equally weighted 
Combo_res<-data.frame(rbind(Ex.Spp.LBSPR.LHin$SPR,Ex.Spp.LBSPR.LinfkM_MK$SPR,Ex.Spp.LBSPR.LinfMK$SPR))
Combo_res<-data.frame("Combined",Combo_res)
colnames(Combo_res)<-c("Scenario","SPR")
LBSPRs.out.df.gg<-rbind(LBSPRs.out.df.gg,Combo_res)
#Use ggplot boxplot to visualize results. Horizontal line is the SPR target
ggplot(LBSPRs.out.df.gg,aes(Scenario,SPR, group=Scenario))+geom_boxplot(outlier.colour="red")+geom_hline(yintercept=c(0.59,0.43),col=c("blue","red"),lty=2,lwd=1.25)
