####################
### Control Rule ###
####################

  SPR_Eff_CR<-function(SPR,SPR_targ,slope.up,slope.dl,shape.up,shape.dl,CRname)
  {
    if(SPR/SPR_targ>=1)
    {
      Vt<-shape.up*((SPR/SPR_targ)-1)^3+slope.up*((SPR/SPR_targ)-1)
      Vt.out<-data.frame(1+Vt,CRname)
      colnames(Vt.out)<-c("E_buffer","CR_name")
    }
    
    if(SPR/SPR_targ<1)
    {
      Vt<-shape.dl*((SPR/SPR_targ)-1)^3+slope.dl*((SPR/SPR_targ)-1)
      Vt.out<-data.frame(1+Vt,CRname)
      colnames(Vt.out)<-c("E_buffer","CR_name")
    }
    
    # if(SPR/SPR_targ>=1 & SPRt-SPRc>=0){D=0.2}
    # if(SPR/SPR_targ>=1 & SPRt-SPRc<0){D=0}
    # if(SPR/SPR_targ<1 & SPRt-SPRc>=0){D=0}
    # if(SPR/SPR_targ<1 & SPRt-SPRc<0){D=0.2}
    # Vt_jn<-1+[shape*((SPR/SPR_targ)-1)^3+slope*((SPR/SPR_targ)-1)]*[1+D]
    
    return(Vt.out)
  }
  ###########################################
  
