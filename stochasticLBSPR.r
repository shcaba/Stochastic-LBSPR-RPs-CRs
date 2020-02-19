library(LBSPR)
library(MASS)
library(tmvtnorm)
library(reshape)
library(reshape2)
library(ggplot2)

###################################################################
### Estimating uncertainty in LBSPR through life history values ###
###################################################################
### Monte Carlo draws of life history values###

#LH_in:
#Option 1: read in LinfMk.file with L and M/k values
#Option 2: Draw Linf, k, and M values, then calculate M/k
#Option 3: Draw Linf, k, and M/k values

LBSPR_stochastic<-function(
  Length.dat, #.csv file of lengths
  dataType.in="freq", #or raw
  Nsamp=1000,
  LH_in=1, 
  LH.file=NA,
  LH.plot = TRUE,
  RawBinWidth=2,
  Linf=NA,
  Linf_CV=0.1,
  k=NA,
  k_CV=0.1,
  Cor_Linf_k=-0.9,
  M_prior=NA,
  M_k=NA,
  M_k_CV=0.2,
  Cor_Linf_M_k=0.2,
  L50=NA,
  L95=NA,
  SPR.flag=0.99,
  Calc.RPs = TRUE,
  maxage=50,        #Maximum age for YPR analysis
  agestep=1,        #Age step for YPR
  maxF = 5,         #Max F to search for YPR
  F_step = 0.01     #F step to search for YPR
  )
{
  #Set-up output onjects, SPR and variances
  LBSPR.out_SPR<-LBSPR.out_vars<-LBSPR.out_Sel50<-LBSPR.out_Sel95<-LBSPR.out_F<-LBSPR.out_F_M<-LBSPR.out_Sel50_vars<-LBSPR.out_Sel95_vars<-LBSPR.out_F_vars<-list()
  RPs_Fmax<-RPs_SPR<-F_wted_LBSPR<-F_wted_YPR<-list() #Reference points list
  if(LH_in==1)
  {
    LH.in<-read.csv(LH.file)
    LH.samps4LBSPR<-data.frame(LH.in[sample(nrow(LH.in),size=Nsamp,replace=TRUE),])
    if(LH.plot==TRUE)
    {
      plot(LH.samps4LBSPR$Linf,LH.samps4LBSPR$M_k)
    }
    }
    
  #Option 2: Establish Linf, k and M distributions, then get M/k
  if(LH_in==2)
  {
    Covar_Linf_k<-Cor_Linf_k*(Linf_CV*Linf*k_CV*k) #caluclate covariance matrix for multivariate sampling
    Linf_k_sigma <- matrix(c((Linf*Linf_CV)^2,Covar_Linf_k,Covar_Linf_k,(k*k_CV)^2),2,2) #Set-up the variance-covariance matrix for multivariate sampling
    L50_Linf <- L50/Linf #Ratio between L50 and Linf. Want to keep this constant
    L95_Linf <- L95/Linf #Ratio between L95 and Linf. Want to keep this constant
    #Draw Nsamp values for Linf and k from mulitvariate random normal distribution
    Linf_k_samps<-rtmvnorm(Nsamp,c(Linf,k),Linf_k_sigma,lower=rep(0,length(c(Linf,k))))
    #Plot the Linf and k values to make sure random draws are consistent with median values
    if(LH.plot == TRUE)
    {
      plot(Linf_k_samps,xlab="Linf",ylab="k")
    #Point is the median values for Linf and k
      points(median(Linf_k_samps[,1]),median(Linf_k_samps[,2]),pch=21,bg="blue",col="blue",cex=1.2)
    #Input values for Linf and k
      abline(v=Linf,h=k,col="blue")
    }  
    #Sample Nsamp values of M
    M.samps<-sample(M_prior,Nsamp)
    M_k_pt<-median(M_prior)/k
    #Make object of Linf, k, and M/k values
    Linf_k_Mk<-cbind(Linf_k_samps,M.samps/Linf_k_samps[,2])
    #Check Linf vs M/k samples
    if(LH.plot == TRUE)
      {
        plot(Linf_k_Mk[,c(1,3)],xlab="Linf",ylab="M/k")
      #Point is the median values for Linf and k
        points(median(Linf_k_samps[,1]),M_k_pt,pch=21,bg="blue",col="blue",cex=1.2)
      #Input values for Linf and k
        abline(v=Linf,h=M_k_pt,col="blue")
    }
    if(dataType.in=="freq")
    {
      #Check maximum Linf value 
      #You need to do this to make sure you data file has enough length bins for all parameter draws
      #If your current data file length bins end before this max value, add more length bins to go above the max Linf value
      max(Linf_k_Mk[,1])
    }
    LH.samps4LBSPR<-data.frame(cbind(Linf_k_Mk,L50_Linf*Linf_k_Mk[,1],L95_Linf*Linf_k_Mk[,1]))
    colnames(LH.samps4LBSPR)<-c("Linf","k","M_k","L50","L95")
  }

  #Option 3: Establish Linf and M/k distributions
  if(LH_in==3)
  {
    #caluclate covariance matrix for multivariate sampling
    Covar_Linf_M_k<-Cor_Linf_M_k*(Linf_CV*Linf*M_k_CV*M_k) 
    Linf_M_k_sigma <- matrix(c((Linf*Linf_CV)^2,Covar_Linf_M_k,Covar_Linf_M_k,(M_k*M_k_CV)^2),2,2) #Set-up the variance-covariance matrix for multivariate sampling
    L50_Linf <- L50/Linf #Ratio between L50 and Linf. Want to keep this constant
    L95_Linf <- L95/Linf #Ratio between L95 and Linf. Want to keep this constant
    #Draw Nsamp values for Linf and k from mulitvariate random normal distribution
    Linf_M_k_samps<-rtmvnorm(Nsamp,mean=c(Linf,M_k),sigma=Linf_M_k_sigma,lower=rep(0,length(c(Linf,M_k))))
    if(LH.plot==TRUE)
    {
      #Plot the Linf and k values to make sure random draws are consistent with median values
      plot(Linf_M_k_samps,xlab="Linf",ylab="M/k")
      #Point is the median values for Linf and k
      points(median(Linf_M_k_samps[,1]),median(Linf_M_k_samps[,2]),pch=21,bg="blue",col="blue",cex=1.2)
      #Input values for Linf and k
      abline(v=Linf,h=M_k,col="blue")
    }
    if(dataType.in=="freq")
    {
    #Check maximum Linf value 
    #You need to do this to make sure you data file has enough length bins for all parameter draws
    #If your current data file length bins end before this max value, add more length bins to go above the max Linf value
    max(Linf_M_k_samps[,1])
    }
    LH.samps4LBSPR<-data.frame(cbind(Linf_M_k_samps,L50_Linf*Linf_M_k_samps[,1],L95_Linf*Linf_M_k_samps[,1]))
    colnames(LH.samps4LBSPR)<-c("Linf","M_k","L50","L95")
  }
  
  #Run LBSPR analysis for all parameter draws
  for(i in 1:Nsamp)
  {
    print(paste0("Sample ",i))
    #Set-up initial pars file
    Spp_pars.samp<-new("LB_pars")
    Spp_pars.samp@L_units <- "cm" 

    Spp_pars.samp@Linf <- LH.samps4LBSPR$Linf[i]
    Spp_pars.samp@MK <- LH.samps4LBSPR$M_k[i] 
    Spp_pars.samp@L50 <- LH.samps4LBSPR$L50[i] 
    Spp_pars.samp@L95 <- LH.samps4LBSPR$L95[i]
    Spp_pars.samp@CVLinf <- Linf_CV
    print(LH.samps4LBSPR[i,])

  if(dataType.in=="raw")
    {
      Spp_pars.samp@BinMin <- 0 
      Spp_pars.samp@BinMax <- max(LH.samps4LBSPR$Linf)+2 
      Spp_pars.samp@BinWidth <- RawBinWidth
    }

    #Spp_sim.samp <- LBSPRsim(Spp_pars.samp)
    Spp_lts.samp<- new("LB_lengths", LB_pars=Spp_pars.samp, file=paste0(Length.dat),dataType=dataType.in, header=TRUE)
    Spp.LBSPR.fit<-LBSPRfit(Spp_pars.samp, Spp_lts.samp,useCPP = TRUE)
    
    LBSPR.out_SPR[[i]]<-Spp.LBSPR.fit@SPR
    LBSPR.out_Sel50[[i]]<-Spp.LBSPR.fit@Ests[,1]
    LBSPR.out_Sel95[[i]]<-Spp.LBSPR.fit@Ests[,2]
    LBSPR.out_F[[i]]<-as.numeric(Spp.LBSPR.fit@Ests[,3])*LH.samps4LBSPR$M_k[i]*LH.samps4LBSPR$k[i]
    LBSPR.out_F_M[[i]]<-as.numeric(Spp.LBSPR.fit@Ests[,3])
    LBSPR.out_Sel50_vars[[i]]<-Spp.LBSPR.fit@Vars[,1]
    LBSPR.out_Sel95_vars[[i]]<-Spp.LBSPR.fit@Vars[,2]
    LBSPR.out_F_vars[[i]]<-Spp.LBSPR.fit@Vars[,3]
    
    #if(ncol(LBSPR.out_Sel50_vars[[i]])>1){Spp.LBSPR.fit@Vars[,1]}

    if(Calc.RPs == TRUE & !is.na(k))
    {
      RPs_Fmax_temp<-RPs_SPR_temp<-Fmax_wt_LBSPR_temp<-Fmax_wt_YPR_temp<-list()
      for(ii in 1:length(LBSPR.out_Sel50[[i]]))
      {
      YPR.out<-YPR_SPR_RPs(maxage=maxage,agestep=agestep,
                         Linf=LH.samps4LBSPR$Linf[i],
                         k=LH.samps4LBSPR$k[i],t0=0, 
                         Lmat50=LH.samps4LBSPR$L50[i],
                         Lmat95=LH.samps4LBSPR$L50[i],
                         Sel50=Spp.LBSPR.fit@Vars[ii,1], 
                         Sel95=Spp.LBSPR.fit@Vars[ii,2], 
                         M =LH.samps4LBSPR$M_k[i]*LH.samps4LBSPR$k[i], 
                         maxF = maxF, F_step = F_step) 
      RPs_Fmax_temp[[ii]]<-YPR.out$Fmax_SPR[1]              
      RPs_SPR_temp[[ii]]<-YPR.out$Fmax_SPR[2]         

    #Calculate weighted F values
    Fmax_wt_YPR_temp<-F_weighted(Linf=LH.samps4LBSPR$Linf[i],
                                 k=LH.samps4LBSPR$k[i],
                                 t0=0, 
                                 M =LH.samps4LBSPR$M_k[i]*LH.samps4LBSPR$k[i],
                                 Fmax=YPR.out$Fmax_SPR[1],
                                 Sel50=Spp.LBSPR.fit@Vars[ii,1], 
                                 Sel95=Spp.LBSPR.fit@Vars[ii,2], 
                                 maxage=maxage,
                                 agestep=agestep)
    Fmax_wt_LBSPR_temp<-F_weighted(Linf=LH.samps4LBSPR$Linf[i],
                                   k=LH.samps4LBSPR$k[i],
                                   t0=0, 
                                   M =LH.samps4LBSPR$M_k[i]*LH.samps4LBSPR$k[i],
                                   Fmax=LBSPR.out_F[[i]][ii],
                                   Sel50=Spp.LBSPR.fit@Vars[ii,1], 
                                   Sel95=Spp.LBSPR.fit@Vars[ii,2], 
                                   maxage=maxage,
                                   agestep=agestep)
      }
    RPs_Fmax[[i]]<-do.call(cbind,RPs_Fmax_temp)
    RPs_SPR[[i]]<-do.call(cbind,RPs_SPR_temp)
    F_wted_LBSPR[[i]]<-do.call(cbind,RPs_Fmax_temp)
    F_wted_YPR[[i]]<-do.call(cbind,RPs_SPR_temp)
    }
  }
  
  #Process data for ggplot
  LBSPRs.vals.out<-list()
  LBSPRs.vals.out[[1]]<-LH.samps4LBSPR
  LBSPRs.vals.out[[2]]<-do.call(rbind,LBSPR.out_SPR)
  LBSPRs.vals.out[[3]]<-do.call(rbind,LBSPR.out_Sel50)
  LBSPRs.vals.out[[4]]<-do.call(rbind,LBSPR.out_Sel95)
  LBSPRs.vals.out[[5]]<-do.call(rbind,LBSPR.out_F)
  LBSPRs.vals.out[[6]]<-do.call(rbind,LBSPR.out_Sel50_vars)
  LBSPRs.vals.out[[7]]<-do.call(rbind,LBSPR.out_Sel95_vars)
  LBSPRs.vals.out[[8]]<-do.call(rbind,LBSPR.out_F_vars)
  LBSPRs.vals.out[[9]]<-do.call(rbind,LBSPR.out_F_M)
  #Remove non-converged runs
  LBSPRs.vals.out[[2]][is.na(LBSPRs.vals.out[[6]])]<-NA
  LBSPRs.vals.out[[3]][is.na(LBSPRs.vals.out[[6]])]<-NA
  LBSPRs.vals.out[[4]][is.na(LBSPRs.vals.out[[6]])]<-NA
  LBSPRs.vals.out[[5]][is.na(LBSPRs.vals.out[[6]])]<-NA
  LBSPRs.vals.out[[9]][is.na(LBSPRs.vals.out[[6]])]<-NA
  #Remove runs with SPR > flagged value
  LBSPRs.vals.out[[2]][LBSPRs.vals.out[[2]]>SPR.flag]<-NA
  LBSPRs.vals.out[[3]][is.na(LBSPRs.vals.out[[2]])]<-NA
  LBSPRs.vals.out[[4]][is.na(LBSPRs.vals.out[[2]])]<-NA
  LBSPRs.vals.out[[5]][is.na(LBSPRs.vals.out[[2]])]<-NA
  LBSPRs.vals.out[[9]][is.na(LBSPRs.vals.out[[2]])]<-NA
  
  LBSPRs.vals.out[[10]]<-mapply(function(x) length(na.omit(LBSPRs.vals.out[[2]][,x])),x=1:ncol(LBSPRs.vals.out[[2]]))
  if(Calc.RPs == TRUE & !is.na(k))
  {
    LBSPRs.vals.out[[11]]<-do.call(rbind,RPs_Fmax)
    LBSPRs.vals.out[[12]]<-do.call(rbind,RPs_SPR)
    LBSPRs.vals.out[[13]]<-do.call(rbind,F_wted_LBSPR)
    LBSPRs.vals.out[[14]]<-do.call(rbind,F_wted_YPR)
    LBSPRs.vals.out[[11]][LBSPRs.vals.out[[11]][,1]==maxF,1]<-LBSPRs.vals.out[[12]][LBSPRs.vals.out[[11]][,1]==maxF,1]<-NA
    LBSPRs.vals.out[[13]][LBSPRs.vals.out[[11]][,1]==maxF,1]<-LBSPRs.vals.out[[14]][LBSPRs.vals.out[[11]][,1]==maxF,1]<-NA
  }
  
  if(Calc.RPs == TRUE & !is.na(k)){names(LBSPRs.vals.out)<-c("LH_inputs","SPR","Sel50","Sel95","F","Sel50_Var","Sel95_Var","F_Var","F_M","Kept samples","Fmax_MSY_RPs","SPR_MSY_RPs","weightedF_LBSPR","weightedF_YPR")}
  else names(LBSPRs.vals.out)<-c("LH_inputs","SPR","Sel50","Sel95","F","Sel50_Var","Sel95_Var","F_Var","F_M","Kept samples")
  return(LBSPRs.vals.out)
}


###################################################
###################################################

