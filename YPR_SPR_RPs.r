    Numbers<-function(R_in,M,F_yr,Sel,ages,a_i)
        {
          M_step<-M*a_i
          F_step<-F_yr*a_i
          Z_step<-M_step+(F_step*Sel)
          Z<-M+(F_yr*Sel)
          nums<-matrix(data= ages*0, nrow= length(ages), ncol= 2) #Empty matrix for numbers at age per sex
          nums[1,1:2]<-R_in/2       #Initialize female & male columns
          for (i in 1:2)        #Loop to calculate numbers for years 1 to max-1
          {
            for(j in 2:(length(ages)-1))
            {
              nums[j,i]<-nums[j-1,i]*exp(-(Z_step[j-1]))
            }
          }
          
          for (ii in 1:2)   #Loop to calculate numbners in plus group
          {
            nums[length(ages),ii]<-(nums[length(ages)-1,ii]*exp(-(Z_step[length(ages)-1])))/(1-exp(-Z[length(ages)]))
          }
          return(nums)
        }
        
        SpawnB<-function(R_in,Maturity,Weight,M,F_yr,Sel,ages,a_i)
        {
          numbers<-Numbers(R_in,M,F_yr,Sel,ages,a_i)
          SBPRF<-sum(Maturity*Weight*numbers[,1])
          return(SBPRF)
        }
        
    F_weighted<-function(Linf,k,t0,M,Fmax,Sel50,Sel95,maxage,agestep)
      {
        R_in<-1000
        ages<-seq(0,maxage,agestep)
        Lts<-Linf*(1-exp(-k*(ages-t0)))
        sel_vec<-1/(1+exp(-(Lts-Sel50)/((Sel95-Sel50)/log(19))))
        numbers<-Numbers(R_in,M,F_yr=Fmax,Sel=sel_vec,ages,a_i=agestep)[,1]
        stand.nums<-numbers[-1]/sum(numbers[-1])
        F_Sel<-Fmax*sel_vec[-1]
        F_wted<-sum(F_Sel*stand.nums)
        return(F_wted)
      }

    YPR_SPR_RPs<-function (maxage,agestep=1,
        Linf,k,t0, 
        LW_a=0.000001,LW_b=3, 
        Lmat50,Lmat95,Sel50,Sel95, 
        M, maxF = 2, F_step = 0.01, 
        plot.YPRSPR = FALSE) 
    {
        
        R_in<-1000
        ages<-seq(0,maxage,agestep)
        Lts<-Linf*(1-exp(-k*(ages-t0)))
        mat_vec<-1/(1+exp(-(Lts-Lmat50)/((Lmat95-Lmat50)/log(19))))
        sel_vec<-1/(1+exp(-(Lts-Sel50)/((Sel95-Sel50)/log(19))))
        wt<-LW_a*Lts^LW_b
        F_search<-seq(0,maxF,F_step)
        YPR<-SPF<-SPR<-matrix(NA,length(F_search),2)
        colnames(YPR)<-c("F","SPF")
        colnames(SPF)<-c("F","SPF")
        colnames(SPF)<-c("F","SPR")
        SPF0<-SpawnB(R_in=R_in,Maturity=mat_vec,Weight=wt,M=M,F_yr=0,Sel=sel_vec,ages=ages,a_i=agestep)

    for(i in 1:length(F_search))
    {
        YPR[i,1]<-SPF[i,1]<-SPR[i,1]<-F_search[i]
        Nums_fished<-Numbers(R_in=R_in,M=M,F_yr=F_search[i],Sel=sel_vec,ages=ages,a_i=agestep)[,1]
        YPR[i,2]<-sum(Nums_fished*(1-exp(-F_search[i]*sel_vec))*wt)
        SPF[i,2]<-SpawnB(R_in=R_in,Maturity=mat_vec,Weight=wt,M=M,F_yr=F_search[i],Sel=sel_vec,ages=ages,a_i=agestep)
    }

    SPR[,2]=SPF[,2]/SPF0
    YPR_MSY<-max(YPR[,2])
    #Plot YPR and SPR
    if(plot.YPRSPR == TRUE)
    {
        plot(YPR,pch=21,col="blue",lwd=1.5,xlab="Fishing mortality",ylab="YPR")
        par(new = T)
        plot(SPR, pch=16, col="red",axes=F, xlab=NA, ylab=NA, cex=1.2,ylim=c(0,1))
        axis(side = 4)
        mtext(side = 4, line = 3, "SPR")
    }
    Fmax<-YPR[,1][which(YPR[,2] == YPR_MSY)]
    SPR_Fmax<-SPF[,2][which(YPR[,2] == YPR_MSY)]/SPF0
    Fmax_SPR<-matrix(c(Fmax,SPR_Fmax),1,2)
    colnames(Fmax_SPR)<-c("Fmax","SPR_Fmax")
    Fmax_SPR.out<-list(YPR=YPR,SPR=SPR,Fmax_SPR=Fmax_SPR)
    return(Fmax_SPR.out)
    }

