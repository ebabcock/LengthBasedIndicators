#This file conducts single species length based indicator analysis, and makes the 
#summary tables and figures. Specify Belize or Guatemala

source("BabcockFunctionsLength2020.r")
library(mvtnorm) 
library(LBSPR)
library(metafor)
library(DHARMa)
library(lemon)
library(egg)

### For Belize or Guatemala
region<-c("Belize","Guatemala")[1]
if(region=="Belize") {
 Belizesimdf<-list()
 Belizelmvals<-list()
} else {
 Guatemalasimdf<-list()
 Guatemalalmvals<-list()

}
for(pars in 1:3) {
parsource<-c("Data","FL","LBBnoprior")[pars]
priorsource<-c("Data","FL","Default")[pars]
if(region=="Belize") {
 datfile<-Belize
 datfile$length<-datfile$TL2_CM
 datfile$sciname<-datfile$scinameFishbase
 if(parsource=="Data")  spfile<-BelizeAll
 if(parsource=="FL")  spfile<-BelizeAllFL
 if(parsource=="LBBnoprior")  spfile<-BelizeAllLBBnoprior
 if(parsource=="LBBFLprior")  spfile<-BelizeAllLBBFLprior
 sptodo<-which(spfile$n>=15)
 sptoplot<-sptoplotB
 datfile$Date<-as.Date(datfile$Patrol.Start.Date,format="%m/%d/%Y")
} else {
 datfile<-Guatemala
 spfile<-GuatemalaAll
 if(parsource=="Data")  spfile<-GuatemalaAll
 if(parsource=="FL")  spfile<-GuatemalaAllFL
 if(parsource=="LBBnoprior")  spfile<-GuatemalaAllLBBnoprior
 if(parsource=="LBBFLprior")  spfile<-GuatemalaAllLBBFLprior
 parfile<-Guatemalapars
 datfile$length<-datfile$Length
 datfile$sciname<-datfile$scinameFishbase
 sptodo<-which(spfile$n>=15)
 spfile$Loo<-spfile$Linf
 sptoplot<-sptoplotG
}

# Setup sims
nsims<-10000  #Should be several thousand for publication
spsum=list()
sum.mat=list()
length.boot=list()
simdf<-NULL
spfile$Lc<-rep(NA,dim(spfile)[1])
spfile$n.Lc<-rep(NA,dim(spfile)[1])
for(i in sptodo) {
  spfile$Lc[i]<-Lc.func(datfile$length[datfile$sciname==spfile$Species[i]])
  spfile$n.Lc[i]<-length(datfile$length[datfile$sciname==spfile$Species[i] & datfile$length>=spfile$Lc[i]])
  }
spfile$LoptBH<-Lopt.Bev(spfile$Loo,spfile$M,spfile$K)
spfile$LoptLm<-Lopt.Lm(spfile$Lm)
ggplot(spfile,aes(x=LoptLm,y=LoptBH))+geom_point()+geom_abline(aes(intercept=0,slope=1))
spfile$Lopt=(Lopt.Bev(spfile$Loo,spfile$M,spfile$K)+
 Lopt.Lm(spfile$Lm))/2
summary(spfile$Lopt)
paramsin<-c("Lmax","Lm","M","K","Loo","Lopt","Lc")
paramsboot<-c("Lave","Pmat","Popt","Pmega","Z","F","F.M")
paramslbspr<-c("SPR","SPR.F.M","SL50","SL95")

for(i in 1:length(sptodo))   {
 sp<-sptodo[i]
 lengthdat<-datfile$length[datfile$sciname==spfile$Species[sp]]
 temp<-BootFunc1(nsims,sp,lengthdat,spfile)
 simdf<-rbind(simdf,temp)
 print(i)
}

summary(simdf)
simdf$Species<-spfile$Species[simdf$sp]
simdf$Common<-spfile$Common[simdf$sp]
simdf$n<-spfile$n[simdf$sp]
simdf$Lmax<-spfile$Lmax[simdf$sp]
simdf$Trophic<-spfile$Trophic[simdf$sp]
x<-order(spfile$Lmax)
simdf$Common<-factor(simdf$Common,levels=spfile$Common[x])
simdf[simdf$SE>10 & !is.na(simdf$SE),]

if(region=="Belize") { y<-filter(BelizeResults$lbbres,Priors==priorsource,Region==region) 
 } else {
  y<-filter(GuatemalaResults$lbbres,Priors==priorsource,Region==region)
}
z<-data.frame(Parameter=rep(c("LBB.F.M","BB0"),each=nrow(y)),
  Value=c(y$FM.med,y$BB0.med),SE=NA,LCI=c(y$FM.lcl,y$BB0.lcl),
  Mean=c(y$FM.med,y$BB0.med),
  UCI=c(y$FM.ucl,y$BB0.ucl),Species=rep(y$Species,2)  )
z$Value[z$Parameter=="LBB.F.M" & z$Value>4]<-4
z$sp<-match(z$Species,spfile$Species)
z$Common<-spfile$Common[z$sp]
z$Lmax<-spfile$Lmax[z$sp]
z$Trophic<-spfile$Trophic[z$sp]
z$n<-spfile$n[z$sp]
summary(z)
sort(names(z))
sort(names(simdf))

x<-match(y$Species,simdf$Species)
summary(x)

indicators<-c("Pmat","Popt","Pmega","F.M","SPR.F.M","LBB.F.M", "SPR","BB0")     
simdf<-bind_rows(simdf,z)
simdf$Indicator<-factor(simdf$Parameter,levels=indicators)
summary(simdf$Indicator)

#Stats
lmlistL<-list()
lmlistT<-list()
lmvalsL<-data.frame(n=rep(NA,6),P=rep(NA,6),a=rep(NA,6),b=rep(NA,6))
lmvalsT<-data.frame(n=rep(NA,6),P=rep(NA,6),a=rep(NA,6),b=rep(NA,6))
simdf$fitT<-simdf$fitL<-
  simdf$lwrT<-simdf$lwrL<-
  simdf$uprT<-simdf$uprL<-rep(NA,dim(simdf)[1])
simdf$included<-rep(0,dim(simdf)[1])

for(i in 1:length(indicators)) {
  whichrows<-which(simdf$Indicator==indicators[i]
    & !is.na(simdf$Mean))
  if(length(whichrows)>0) {
  x<-simdf[whichrows,]
#  y<-rma(Mean ~ Lmax, SE, method="FE", data=x)
  y<-lm(Mean ~ Lmax, data=x)
  lmlistL[[i]]<-y  
  #  lmvals[i,c("P","a","b")]<-c(y$pval[2],coef(y))
  lmvalsL[i,c("P","a","b")]<-c(anova(y)[1,5],coef(y))
  z<-predict(y,interval="confidence")
  lmvalsL[i,"n"]<-length(whichrows)
#  simdf$fit[whichrows]<-z$pred
#  simdf$lwr[whichrows]<-z$ci.lb
#  simdf$upr[whichrows]<-z$ci.ub
  simdf$fitL[whichrows]<-z[,"fit"]
  simdf$lwrL[whichrows]<-z[,"lwr"]
  simdf$uprL[whichrows]<-z[,"upr"]
  simdf$included[whichrows]<-1
  y<-lm(Mean ~ Trophic, data=x)
  lmlistT[[i]]<-y  
  lmvalsT[i,c("P","a","b")]<-c(anova(y)[1,5],coef(y))
  z<-predict(y,interval="confidence")
  lmvalsT[i,"n"]<-length(whichrows)
  simdf$fitT[whichrows]<-z[,"fit"]
  simdf$lwrT[whichrows]<-z[,"lwr"]
  simdf$uprT[whichrows]<-z[,"upr"]
}}

lmvalsT$P<-round(lmvalsT$P,3)
lmvalsT$a<-round(lmvalsT$a,2)
lmvalsT$b<-round(lmvalsT$b,4)
lmvalsL$P<-round(lmvalsL$P,3)
lmvalsL$a<-round(lmvalsL$a,2)
lmvalsL$b<-round(lmvalsL$b,4)

lmvals<-cbind(data.frame(Indicator=indicators),lmvalsL,lmvalsT)
lmvals
write.csv(lmvals,paste0(region,parsource,"lmvals.csv"))

if(region=="Belize") {
  print(priorval)
  Belizesimdf[[pars]]<-simdf
  Belizelmvals[[pars]]<-lmvals
} else {
  Guatemalasimdf[[pars]]<-simdf
  Guatemalalmvals[[pars]]<-lmvals
}

}

#Pick region and parameters for figures and tables
region<-"Belize"
pars<-1
parsource<-c("Data","FL","LBBnoprior")[pars]
if(region=="Belize") {
  simdf<-Belizesimdf[[pars]]   
  lmvals<-Belizelmvals[[pars]]
  spfile<-BelizeAll
  } else {
 simdf<-Guatemalasimdf[[pars]]
 lmvals<-Guatemalalmvals[[pars]]
 spfile<-GuatemalaAll
}
#Plots
df1<-filter(simdf,Indicator %in% indicators & included==1 )
x<-match(df1$Species,spfile$Species)
summary(x)
df1$Family<-spfile$Family[x]
df1$Indicator2<-df1$Parameter
df1$Indicator2[df1$Parameter=="Pmat"]="P[mat]"
df1$Indicator2[df1$Parameter=="Popt"]="P[opt]"
df1$Indicator2[df1$Parameter=="Pmega"]="P[mega]"
df1$Indicator2[df1$Parameter=="F.M"]="F/M (BH)"
df1$Indicator2[df1$Parameter=="SPR.F.M"]="F/M (SPR)"
df1$Indicator2[df1$Parameter=="LBB.F.M"]="F/M (LBB)"
df1$Indicator2[df1$Parameter=="BB0"]="B/B[0]"
df1$Indicator2<-factor(df1$Indicator2)
summary(df1$Indicator2)
df1$Indicator2<-factor(df1$Indicator2,levels=
    levels(df1$Indicator2)[c(5,7,6,2,4,3,8,1)])
levels(df1$Indicator2)
df1$lwrL[df1$Indicator %in% indicators[lmvals$P>0.05]]=NA
df1$fitL[df1$Indicator %in% indicators[lmvals$P>0.05]]=NA
df1$uprL[df1$Indicator %in% indicators[lmvals$P>0.05]]=NA
df1$LCI[df1$LCI<0 &!is.na(df1$LCI)]<-0
df1$label2<-rep("",dim(df1)[1])
x<-match(df1$Species,sptoplot,nomatch = 0)
df1$label2[x!=0]<-(letters[1:15])[x]
table(df1$label2)
head(df1)

datText <- data.frame(label = paste0(" (",letters[1:8],")"),
  Indicator2=levels(df1$Indicator2),
  expandval=c(rep(0.2,6),rep(0.4,2)))
datText$Indicator2<-factor(datText$Indicator2,levels=levels(df1$Indicator2))
datText




 #https://andburch.github.io/ggplot_facets/
 old_plot<-ggplot(df1,aes(x=Lmax))+
  geom_point(data=filter(df1,label2==""),aes(x=Lmax,y=Value),size=1)+
  geom_text(aes(y=Value,label=label2))+
 # geom_point(aes(y=Value))+
  ggtitle("")+
  geom_hline(data=filter(df1,Parameter %in% c("Pmat","Popt","F.M","SPR.F.M","LBB.F.M")),aes(yintercept=1),lty=2)+
  geom_hline(data=filter(df1,Parameter %in% c("Pmega")),aes(yintercept=0.2),lty=2)+
  geom_hline(data=filter(df1,Parameter %in% c("SPR","BB0")),aes(yintercept=0.3),lty=2)+
  facet_rep_wrap(Indicator2~.,scales="free",
    strip.position="left",labeller=label_parsed)+
  theme_classic()+
  theme(strip.background = element_blank(),
    strip.placement="outside",axis.line=element_line(),
    panel.border=element_blank(),
    strip.text.y = element_text(size = 11),
    strip.text=element_text(hjust=0.5),
    legend.position="bottom",
    legend.title = element_blank(),
    legend.text=element_text(size=8))+
  xlab(expression(L[max]))+ylab("")+
  scale_x_continuous(breaks=seq(0,300,by=100))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  geom_text( data    = datText,
       mapping = aes(x = -Inf, y = Inf, label = label),
      hjust   = 0,  vjust   =1)
 old_plot
 new_plot <- old_plot +
 #  geom_errorbar(aes(ymin=LCI,ymax=UCI,color=Family))+
   geom_line(aes(y=fitL))+
   geom_ribbon(aes(y=fitL,ymin=lwrL,ymax=uprL),alpha=.2)


 old_plot_data <- ggplot_build(old_plot)
 new_plot_data <- ggplot_build(new_plot)

 new_plot_data$layout$panel_params <- old_plot_data$layout$panel_params
#
 if(region=="Belize") jpeg(paste0("Figure8.jpg"),height=6,width=6.5,units="in",res=300)
 if(region=="Guatemala") jpeg("Figure9.jpg",height=6,width=6.5,units="in",res=300)
 try(plot(ggplot_gtable(new_plot_data)))
 dev.off()

region
 
 names(simdf)
x<-pivot_wider(simdf,id_cols="Species",names_from="Parameter",values_from=Value) %>%
  dplyr::select(Species,Pmat,Popt,Pmega,F.M,SPR.F.M,LBB.F.M,SPR,BB0)

x$Common<-spfile$Common[match(x$Species,spfile$Species)]
x$Family<-spfile$Family[match(x$Species,spfile$Species)]

write.csv(x[order(x$Family),],file=paste0(region,"SingleIndicators.csv"))
x<-x[order(x$Pmat),]
x[1:10,c("Common","Species")]
dim(x)
x[45:47,c("Common","Species","Pmat")]


