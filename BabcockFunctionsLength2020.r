# Functions for length based analysis

# To deal with virus checker problem
#debug(utils:::unpackPkgZip)

##### Froese and Binohlan functions to convert between life history parameters
#Froese, R., Binohlan, C., 2000. Empirical relationships to estimate asymptotic length,length at first maturity and length at maximum yield per recruit in fishes, with a simple method to evaluate length frequency data. J. Fish. Biol. 56, 758-773.
#Each method converts the parameter after the period to the paramteter before
Linf.Lmax = function(Lmax) 10^(0.044+0.9841*log10(Lmax))
Lopt.Linf=function(Linf) 10^(1.0421*log10(Linf)-.2742) 
Lopt.Lm=function(Lm) 10^(1.053*log10(Lm)-0.0565)
Lm.Linf=function(Linf) 10^(0.898*log10(Linf)-0.0781)
## See the file popdynJFB.xls for the rest of the functions from this paper

##### Functions to calculate Beverton invariants
# Beverton, R.J.H., 1992. Patterns of reproductive strategy parameters in some marine teleost fishes. J. Fish. Biol. 41, 137-160.Lopt.Bev=function(Linf,M,K)  3*Linf/(3+M/K)
Lopt.Bev=function(Linf,M,K)  3*Linf/(3+M/K)
Lm.Linf.bev=function(Linf) 0.66*Linf

#### Functions to calculate the Froese length-based indicator
#Froese, R., 2004. Keep it simple: three indicators to deal with overfishing. Fish. Fish 5, 86-91.
PmatFunc=function(x,Lm) {
  x=x[!is.na(x)]
  temp=length(x[x>=Lm])/length(x)
  if(is.na(Lm)|length(x)<5) temp=NA
  temp
}
PoptFunc=function(x,Lopt) {
  x=x[!is.na(x)]
  temp=length(x[x>=Lopt*0.9 & x<=Lopt*1.1])/length(x)
  if(is.na(Lopt)|length(x)<5) temp=NA
  temp
}
PmegaFunc=function(x,Lopt) {
  x=x[!is.na(x)]
  temp=length(x[x>=1.1*Lopt])/length(x)
  if(is.na(Lopt)|length(x)<5) temp=NA
  temp
}

## Z from Lbar using either Bevertion and Holt (BH.Z) or Ehrhardt and Ault (EA) 
# Beverton, R.J.H., Holt, S.J., 1957. On the Dynamics of Exploited Fish Populations.Chapman and Hall, London.
# Ehrhardt, N.M., Ault, J.S., 1992. Analysis of two length-based mortality models applied to bounded catch length frequencies. J. Am. Fish Soc. 121, 115-122.
BH.Z=function(K,Linf,Lave,Lc)     K*(Linf-Lave)/(Lave-Lc)
EA.Z=function(K,Linf,Lave,Lc,Lup,Zest,Zmin) {
  f1 = function(Z) abs(((Linf-Lup)/(Linf-Lc))^(Z/K)-((Z)*(Lc-Lave)+K*(Linf-Lave))/(Z*(Lup-Lave)+K*(Linf-Lave)))
  nlminb(Zest,f1,lower=Zmin)
  # Note that the iterative function may not converge on the first try
  # Function returns a list containing Z, the objectiv function etc.
}

# Cope and Punt function to calculate whether population is overfished
#Cope, J.M., Punt, A.E., 2009. Length-based reference points for data-limited situations:applications and restrictions. Mar. Coastal Fish. 1, 1-18.
#Output is a vector containing status (0=not overfished, 1=overfished), selectivity category (1-5), and Px/Ptarget, which is less than 1 for overfished populations
cope=function(Pmat,Popt,Pmega,Lmat,Lopt) {
  Pobj=Pmat+Popt+Pmega
  Px=Pmat
  if(Pobj<=1) {
    if(Popt+Pmega==0)  {
      selectivity=1 
      Ptarg=0.25
    } else {
      selectivity=2
      if(Lmat>0.825*Lopt) 
        Ptarg=0.25
       else  
        Ptarg=0.4
      }
     } 
  if(Pobj>1 & Pobj<2) {
    selectivity=3
    if(Lmat>0.825*Lopt)
      Ptarg=0.9
    else
      Ptarg=0.95
  }
  if(Pobj==2) {
    if(Popt<1) {
      selectivity=4 
      Px=Popt
      Ptarg=0.65
      } else 
        if(Popt==1) {
          selectivity=5 
          Px=NA
          Ptarg=NA
        } 
  }
    Px.Ptarg=Px/Ptarg
    if(is.na(Px.Ptarg)) status=-1 else status=ifelse(Px.Ptarg>=1,0,1)
    c("status"=status,"selectivity"=selectivity,"Px.Ptarg"=Px.Ptarg)
}


#Find Lower bound function. This function finds the mode of the length frequency distribution
#The Lc used to calculate Z from Lbar should be this number or higher
Lc.func=function(x) {
  #  x is length sample.Returns Lc
  z=table(x)
  z1=cumsum(z)
  z1=z1/sum(z)
  a=as.numeric(names(z))
  a1=seq(trunc(min(a))+1,max(a),by=1)
  d=loess(z1~a)
  d1=predict(d,newdata=a1)
  d2=d1[2:length(a1)]-d1[1:(length(a1)-1)]
  d3=a1[1:length(d2)][d2==max(d2)]
  plot(a1[1:length(d2)],d2)
  d3
}



# Functions to calculate M from life history parameters
# Hewitt, D.A., Lambert, D.M., Hoenig, J.M., Lipcius, R.N., Bunnell, D.B., Miller, T.J., 2007.Direct and indirect estimates of natural mortality for Chesapeake Bay blue crab.J. Am. Fish. Soc. 136, 1030-1040.
# Jensen, A.L., 1996. Beverton and Holt life history invariants result from optimal trade-off of reproduction and survival. Can. J. Fish. Aquat. Sci. 53, 820-822.
# Methods 1-8 are from Hewitt.
M1func=function(tm)  round(1.65/tm,3)
M2func=function(K)   round(1.6*K,3)    
M3func=function(tmax,K)  round(3*K/(exp(0.38*K*tmax)-1),3)
M4func=function(tmax) round(exp(1.44-0.982*log(tmax)),3)
M5func=function(Linf,K,T) round(10^(-0.0066-0.279*log10(Linf)+0.6543*log10(K)+0.4634*log10(T)),3) #Pauly
M6func=function(Winf,K,T) round(10^(-.2107-0.0824*log10(Winf)+0.6757*log10(k)+0.4627*log10(T)),3) #Pauly
M7func=function(K,tm) round(3*K/(exp(K*tm)-1),3)
M8func=function(W)  round(3*W^(-2.88),3)        
M9func=function(tmax) round(-log(0.01)/tmax,3) # With lower number than Ault
M10func=function(tmax) round(-log(0.05)/tmax,3) #Alagaraja as used by Ault, J.S., Smith, S.G., Luo, J.G., Monaco, M.E., Appeldoorn, R.S., 2008. Length-based assessment of sustainability benchmarks for coral reef fishes in Puerto Rico. Environ. Conserv. 35, 221-231.
M11func=function(tmax) round(-log(0.005)/tmax,3) # With lower bound than Ault
M12func=function(K) round(0.21+1.45*K,3) #Jenson 
M13func=function(K) round(exp(0.51+0.95*log(K)),3) #Jenson
M.Then=function(tmax) 4.899*tmax^(-.916)
M.Then.K=function(Linf,K) 4.118*K^0.73*(Linf*10)^(-0.33) #converted to cm from mm

#Function to bootstrap a length frequncy sample
boot.sample=function(x) {
  sample(x,length(x),replace=TRUE)
}

#### data summary functions
mostfreqfunc=function(x) {
  y=table(x)
  names(y)[y==max(y)][1]
}
standard.error=function(x) {
  sd(x,na.rm=T)/sqrt(length(x))
}

## Simulation function.  This function returns a simulation matrix. Each row is a sim,
## Columns are the parameters
## Takes parmean and parvar from FishLife (Thorson) so they are log transformed
MonteCarloFunc=function(nsims,sp,lengthdat,parmean,parvar,parnames,Lc,Lup,Lmax) {
 y=array(0,dim=c(nsims,26)) #Make array to hold simulation results
 dimnames(y)[[2]]=c("sp","Lmax","Lm","M","K","Loo","Lopt","Lc","Lup","Lave","Pmat","Popt","Pmega","bhz","Z","F","F.M","Lm.Lopt","Pobj","sel","cope","Px.Ptarg","SPR","SPR.FM","SL50","SL95")
 y[,"sp"]<-sp
 y[,parnames]<-exp(rmvnorm(nsims,parmean,parvar))  #Take the exponent because the original values are log transformed
 y[,"Lc"]=Lc
 y[,"Lup"]=Lup
 y[,"Lmax"]=Lmax
 a=Lopt.Bev(y[,"Loo"],y[,"M"],y[,"K"])
 b=Lopt.Lm(y[,"Lm"])
 for(i in 1:nsims)
  y[i,"Lopt"]=runif(1,min(a[i],b[i]),max(a[i],b[i]))
 # Draw length bootstrap sample
 length.boot[[sp]]=replicate(nsims,boot.sample(lengthdat))
 write.table(lengthdat,file="lengthdat.csv",col.names=FALSE,row.names=FALSE,sep=",")
 for(i in 1:nsims) {
  lensam=length.boot[[sp]][,i]
#  write.table(lensam,file="lengthdat.csv",col.names=FALSE,row.names=FALSE,sep=",")
  y[i,"Lave"]=mean(lensam[lensam>=y[i,"Lc"] & lensam<=y[i,"Lup"]],na.rm=TRUE)
  y[i,"Pmat"]=PmatFunc(lensam,y[i,"Lm"]) 
  y[i,"Popt"]=PoptFunc(lensam,y[i,"Lopt"])
  y[i,"Pmega"]=PmegaFunc(lensam,y[i,"Lopt"])
  if(y[i,"Lm"]<y[i,"Loo"]) {
    MyPars <- new("LB_pars",verbose=FALSE)
    MyPars@Species<-spfile$Species[sp]
    MyPars@Linf<-y[i,"Loo"]
    MyPars@L50<-y[i,"Lm"]
    MyPars@L95<-y[i,"Lm"]+1
    MyPars@MK<-y[i,"M"]/y[i,"K"]
    MyPars@L_units <- "cm"
    #  MyPars@BinMin=min(lengthdat)-1
    #  MyPars@BinMax=max(c(fishLH$Linf.max[sp], lengthdat))+1
    #  MyPars@BinWidth=(MyPars@BinMax-MyPars@BinMin)/50
    len1=new("LB_lengths",LB_pars=MyPars,file="lengthdat.csv",dataType="raw",verbose=FALSE)
    # plot(len1@LMids,len1@LData)
    myFit1 <- LBSPRfit(MyPars, len1,verbose=FALSE)
   y[i,"SPR"]<-myFit1@SPR
   y[i,"SPR.FM"]<-myFit1@FM
   y[i,"SL50"]<-myFit1@SL50
   y[i,"SL95"]<-myFit1@SL95
  }
 }
 y[,"bhz"]=BH.Z(y[,"K"],y[,"Loo"],y[,"Lave"],y[,"Lc"])
 Zminval=c(0.001,0.01,0.1,1,2,5)
 ## Use iterative function to calculate the Ehrhardt and Ault Z, if desired
 if(useEA) {
 for(i in 1:nsims) {
  if(y[i,"Lave"]<y[i,"Loo"] & y[i,"Lave"]>(1+y[i,"Lc"])) {
    if(y[i,"bhz"]>0 & y[i,"Loo"]>y[i,"Lup"]) {
      a=EA.Z( y[i,"K"], y[i,"Loo"], y[i,"Lave"], y[i,"Lc"], y[i,"Lup"], y[i,"bhz"],Zminval[1])
      if(a$objective[1]>1e-6)  a=EA.Z( y[i,"K"], y[i,"Loo"], y[i,"Lave"], y[i,"Lc"], y[i,"Lup"], y[i,"bhz"],Zminval[2])
      if(a$objective[1]>1e-6)  a=EA.Z( y[i,"K"], y[i,"Loo"], y[i,"Lave"], y[i,"Lc"], y[i,"Lup"], y[i,"bhz"],Zminval[3])
      if(a$objective[1]>1e-6)  a=EA.Z( y[i,"K"], y[i,"Loo"], y[i,"Lave"], y[i,"Lc"], y[i,"Lup"], y[i,"bhz"],Zminval[4])
      if(a$objective[1]>1e-6)  a=EA.Z( y[i,"K"], y[i,"Loo"], y[i,"Lave"], y[i,"Lc"], y[i,"Lup"], y[i,"bhz"],Zminval[5])
      if(a$objective[1]>1e-6)  a=EA.Z( y[i,"K"], y[i,"Loo"], y[i,"Lave"], y[i,"Lc"], y[i,"Lup"], y[i,"bhz"],Zminval[6])
      if(a$objective[1]<1e-6)  {
        y[i,"Z"]=a$par[1]
      } else  y[i,"Z"]=NA
    } else  y[i,"Z"]=NA
  } else  y[i,"Z"]=NA
}} 
else y[,"Z"]=y[,"bhz"]
  #Calculate the rest of the columns
y[,"F"]=y[,"Z"]-y[,"M"]
y[,"F"][y[,"F"]<0]=0
y[,"F.M"]=y[,"F"]/y[,"M"]
y[,"Lm.Lopt"]=y[,"Lm"]/y[,"Lopt"]
y[,"Pobj"]=(y[,"Pmat"]+y[,"Popt"]+y[,"Pmega"])
# Calculate cope status indicator
for(i in 1:nsims)  {
  temp=cope(y[i,"Pmat"],y[i,"Popt"],y[i,"Pmega"],y[i,"Lm"],y[i,"Lopt"])
  y[i,"sel"]=temp[2]
  y[i,"cope"]=temp[1]
  y[i,"Px.Ptarg"]=temp[3]
}
#Calculate the rest of the columns
y[,"F"]=y[,"Z"]-y[,"M"]
y[,"F"][y[,"F"]<0]=0
y[,"F.M"]=y[,"F"]/y[,"M"]
y[,"Lm.Lopt"]=y[,"Lm"]/y[,"Lopt"]
y[,"Pobj"]=(y[,"Pmat"]+y[,"Popt"]+y[,"Pmega"])
# Calculate cope status indicator
for(i in 1:nsims)  {
  temp=cope(y[i,"Pmat"],y[i,"Popt"],y[i,"Pmega"],y[i,"Lm"],y[i,"Lopt"])
  y[i,"sel"]=temp[2]
  y[i,"cope"]=temp[1]
  y[i,"Px.Ptarg"]=temp[3]
}
y
}

#Function leaving out parameter uncertainty
BootFunc1=function(nsims,sp,lengthdat,spfile) {
 y=array(0,dim=c(nsims,length(paramsboot))) #Make array to hold bootstrap results
 dimnames(y)[[2]]=paramsboot
 # Draw length bootstrap sample if doing
 if(!is.na(spfile$M[sp])) {
  length.boot[[sp]]=replicate(nsims,boot.sample(lengthdat))
  for(i in 1:nsims) {
   lensam=length.boot[[sp]][,i]
   y[i,"Lave"]=mean(lensam[lensam>=spfile$Lc[sp]],na.rm=TRUE)
   y[i,"Pmat"]=PmatFunc(lensam,spfile$Lm[sp]) 
   y[i,"Popt"]=PoptFunc(lensam,spfile$Lopt[sp])
   y[i,"Pmega"]=PmegaFunc(lensam,spfile$Lopt[sp])
  }
  y[,"Z"]=BH.Z(spfile$K[sp],spfile$Loo[sp],y[,"Lave"],spfile$Lc[sp])
  y[,"F"]=y[,"Z"]-spfile$M[sp]
  y[,"F"][y[,"F"]<0]=0
  y[,"F.M"]=y[,"F"]/spfile$M[sp]
  y[,"F.M"][y[,"F.M"]>4]=4
  parvalout<-data.frame(Parameter=c(paramsboot,paramslbspr),
   Value=NA,Mean=NA,SE=NA,LCI=NA,UCI=NA)
  parvalout$Value[1:4]<-c(
   mean(lengthdat[lengthdat>=spfile$Lc[sp]],na.rm=TRUE),
   PmatFunc(lengthdat,spfile$Lm[sp]), 
   PoptFunc(lengthdat,spfile$Lopt[sp]),
   PmegaFunc(lengthdat,spfile$Lopt[sp]))
  parvalout$Value[parvalout$Parameter=="Z"]=
   BH.Z(spfile$K[sp],spfile$Loo[sp],parvalout$Value[parvalout$Parameter=="Lave"],spfile$Lc[sp])
  parvalout$Value[parvalout$Parameter=="F"]=parvalout$Value[parvalout$Parameter=="Z"]-spfile$M[sp]
  parvalout$Value[parvalout$Parameter=="F"][parvalout$Value[parvalout$Parameter=="F"]<0]=0
  parvalout$Value[parvalout$Parameter=="F.M"]=parvalout$Value[parvalout$Parameter=="F"]/spfile$M[sp]
  parvalout$Value[parvalout$Parameter=="F.M" & parvalout$Value>4]=4
  for(i in 1:length(paramsboot)) {
   parvalout[parvalout$Parameter==paramsboot[i],3:6]<-c(mean(y[,paramsboot[i]],na.rm=TRUE),sd(y[,paramsboot[i]],na.rm=TRUE),quantile(y[,paramsboot[i]],c(0.025,0.975),na.rm=TRUE))
  }
  if(spfile$n.Lc[sp]<30) {
    parvalout[5:7,2:6]<-NA
  }
 } else {
  parvalout<-data.frame(Parameter=c(paramsboot,paramslbspr),
   Value=NA,Mean=NA,SE=NA,LCI=NA,UCI=NA)
 }  
  write.table(lengthdat,file="lengthdat.csv",col.names=FALSE,row.names=FALSE,sep=",")
 if(!is.na(spfile$Loo[sp]) & spfile$n[sp]>=40) {
   if((spfile$Lm[sp] < spfile$Linf[sp])) {
    MyPars <- new("LB_pars",verbose=FALSE)
    MyPars@Species<-spfile$Species[sp]
    MyPars@Linf<-spfile$Loo[sp]
    MyPars@L50<-spfile$Lm[sp]
    MyPars@L95<-spfile$Lm[sp]+1
    MyPars@MK<-spfile$M[sp]/spfile$K[sp]
    if(!is.na(spfile$MK[sp])) MyPars@MK<-spfile$MK[sp]
    MyPars@L_units <- "cm"
    len1=new("LB_lengths",LB_pars=MyPars,file="lengthdat.csv",dataType="raw",verbose=FALSE)
    myFit1 <- LBSPRfit(MyPars, len1,verbose=FALSE)
    x=sqrt(myFit1@Vars[1,])/myFit1@Ests[1,]
     if(max(x,na.rm=TRUE)<2) {
     parvalout[c(10,11,9,8),"Value"]<-myFit1@Ests[1,]
     parvalout[c(10,11,9,8),"Mean"]<-myFit1@Ests[1,]
     parvalout[c(10,11,9,8),"SE"]<-sqrt(myFit1@Vars[1,])
     parvalout[c(10,11,9,8),"LCI"]<-parvalout[c(10,11,9,8),"Mean"]-parvalout[c(10,11,9,8),"SE"]*1.96
     parvalout[c(10,11,9,8),"UCI"]<-parvalout[c(10,11,9,8),"Mean"]+parvalout[c(10,11,9,8),"SE"]*1.96
     }
   }
 }
 parvalout[parvalout$Parameter=="SPR.F.M",c("Value","Mean","LCI","UCI")][parvalout[parvalout$Parameter=="SPR.F.M",c("Value","Mean","UCI","UCI")] >4]<-4
 parvalout[parvalout$Parameter=="SPR",c("Value","Mean","LCI","UCI")][parvalout[parvalout$Parameter=="SPR",c("Value","Mean","UCI","UCI")] >1]<-1
 parvalout[parvalout$Parameter %in% c("SPR.F.M","SPR"),c("Value","Mean","LCI","UCI")][parvalout[parvalout$Parameter %in% c("SPR.F.M","SPR"),c("Value","Mean","UCI","UCI")] <0]<-0
 parvalout$sp<-rep(sp,11)
 parvalout
}


## Biodiversity

#x is count of total number of each species
simpson.index=function(x)  {
  #From Hurlbert 1971 
  D2=sum(x/sum(x)*(sum(x)-x)/(sum(x)-1))
  D1=sum(x)/(sum(x)+1)*D2
  1-D1
}

simpson.index2=function(x) 1-sum((x/sum(x))^2)

simpson.se=function(x)  {
  #From Hurlbert 1971 and Simpson 1949 
  N=sum(x)
  p=x/sum(x)
  varD2=(4*N*(N-1)*(N-2)*sum(p^3)+2*N*(N-1)*sum(p^2)-2*N*(N-1)*(2*N-3)*(sum(p^2))^2)/((N*(N-1))^2)
  # varD2=4/N*(sum(p^3)-(sum(p^2))^2)
  varD2=(N/(N-1))*varD2
  sqrt(varD2) 
}

# Pauly growth index function from Pauly et al 1998 Fishbase book. 

phi.func=function(Linf,K) log10(K) + 2*log10(Linf) 
K.phi=function(Linf,phi)  K=10^(phi-2*log10(Linf))

## Logistic selectivity

sel.func=function(L,SL50,SL95) {
  1/(1+exp(-log(19)*(L-SL50)/(SL95-SL50)))
}
plot(seq(1,100,.1),sel.func(seq(1,100,.1),38.22,47.11),type="l")
lines(c(43,43),c(0,1))

##lognormal functions
lnorm.mean=function(x1,x1e) {
 #convert lognormal to normal
 exp(x1+0.5*x1e^2)
}
lnorm.se=function(x1,x1e) {
 #convert lognormal to normal
 ((exp(x1e^2)-1)*exp(2*x1+x1e^2))^0.5
}

#Calculate variance of sum

var.sum<-function(var.x,var.y,cov.xy) var.x+var.y+2*cov(var.xy)
