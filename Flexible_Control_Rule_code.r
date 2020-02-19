####################
### Control Rule ###
####################

  SPR_Eff_CR<-function(SPR,SPR_targ,
    slope.up=c(rep(0,5),seq(0.1,0.25,0.05),0.4,seq(0.05,0.2,0.05),0.35),
    slope.dl=c(rep(0,5),seq(0.1,0.25,0.05),0.4,seq(0.05,0.2,0.05),0.35),
    shape.up=c(seq(0.1,0.3,0.05),rep(0,5),rep(0.05,5)),
    shape.dl=c(seq(0.1,0.3,0.05),rep(0,5),rep(0.05,5)),
    CRname=c(c(paste("linear",1:5)),c(paste("cubic",1:5)),c(paste("cubicpoly",1:5))))
  {
  Vt.out<-list()
  if(length(SPR)>length(SPR_targ)){SPR_targ<-rep(SPR_targ,length(SPR))}
  for(i in 1:length(SPR))
  {
    if(is.na(SPR[i]/SPR_targ[i]))
    {
      Vt.out[[i]]<-data.frame(matrix(NA,length(slope.up),2))
      colnames(Vt.out[[i]])<-c("Modifier","CR_name")
    }
    if(!is.na(SPR[i]/SPR_targ[i]))
    {
      if(SPR[i]/SPR_targ[i]>=1)
      {
        Vt<-1+(shape.up*((SPR[i]/SPR_targ[i])-1)^3+slope.up*((SPR[i]/SPR_targ[i])-1))
        Vt.out[[i]]<-data.frame(Vt,CRname)
        colnames(Vt.out[[i]])<-c("Modifier","CR_name")
      }
  
      if(SPR[i]/SPR_targ[i]<1)
      {
        Vt<-1+(shape.dl*((SPR[i]/SPR_targ[i])-1)^3+slope.dl*((SPR[i]/SPR_targ[i])-1))
        Vt.out[[i]]<-data.frame(Vt,CRname)
        colnames(Vt.out[[i]])<-c("Modifier","CR_name")
      }
    }
    # if(SPR/SPR_targ>=1 & SPRt-SPRc>=0){D=0.2}
    # if(SPR/SPR_targ>=1 & SPRt-SPRc<0){D=0}
    # if(SPR/SPR_targ<1 & SPRt-SPRc>=0){D=0}
    # if(SPR/SPR_targ<1 & SPRt-SPRc<0){D=0.2}
    # Vt_jn<-1+[shape*((SPR/SPR_targ)-1)^3+slope*((SPR/SPR_targ)-1)]*[1+D]
  }
    return(Vt.out)
  }
  ###########################################
  