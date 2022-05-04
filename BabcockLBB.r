library(tidyverse)
library(gridExtra)
setwd("~/Dropbox/BethGlovers/2020 finfish")

region<-c("Belize","Guatemala")[2]
for(pars in 1:3) {
parsource<-c("Data","FL","None")[pars]
if(region=="Belize") {
 datfile<-Belize
 datfile$length<-datfile$TL2_CM
 datfile$sciname<-datfile$scinameFishbase
 if(parsource=="Data")  spfile<-BelizeAll  #Use data based priors
 if(parsource=="FL")  spfile<-BelizeAllFL  # Use Fishlife priors
 if(parsource=="None")  {    #Use default LBB priors for Linf and M/K, but keep data base prior for Lm
   spfile<-BelizeAll
   spfile$Linf[]<-NA
   spfile$M[]<-NA
   spfile$K[]<-NA
   spfile$MK[]<-NA
 }   
 spfileAll<-BelizeAll
} else {
 datfile<-Guatemala
 datfile$length<-datfile$Length
 datfile$sciname<-datfile$scinameFishbase
 if(parsource=="Data")  spfile<-GuatemalaAll  #Use data based priors
 if(parsource=="FL")  spfile<-GuatemalaAllFL  # Use Fishlife priors
 if(parsource=="None")  {    #Use default LBB priors for Linf and M/K, but keep data base prior for Lm
   spfile<-GuatemalaAll
   spfile$Linf[]<-NA
   spfile$M[]<-NA
   spfile$K[]<-NA
   spfile$MK[]<-NA
 }   
 spfile$FBname<-spfile$Common
 spfileAll<-GuatemalaAll
}
sptodo<-which(spfile$n>=40)
sptodo

dat1<-datfile %>% filter(scinameFishbase %in% spfile$Species[sptodo] 
  &!is.na(length)) %>%
  mutate(FBname=spfile$FBname[match(scinameFishbase,spfile$Species)]) %>%
  mutate(Stock=FBname,Year=0,Length=round(length)*10) %>% 
  group_by(Stock,Year,Length) %>%  dplyr::summarize(CatchNo=n())
head(dat1)
length(unique(dat1$Stock))

filename<-paste0(region,"DatLBB.csv")
write.csv(dat1,filename)

id1<-data.frame(File=rep(filename,length(sptodo))) %>%
  mutate(
    Species=spfile$Species[sptodo],	
    Stock=spfile$FBname[sptodo],	
    StartYear=NA,	
    EndYear=NA,	
    Years.user=NA,	
    Year.select=NA,	
    Gears.user=NA,	
    Lcut.user=NA,	
    Lc.user=NA,	
    Lstart.user=NA,	
    Linf.user=NA,
    Linf.user=spfile$Linf[sptodo],	
    MK.user=spfile$MK[sptodo],
    mm.user=FALSE,	
    GausSel=FALSE,	
    MergeLF=FALSE,	
    Pile=0,	
    Lm50=spfile$Lm[sptodo],	
    Comment=NA,	
    Source=NA)

write.csv(id1,paste0(parsource,region,"LBBID.csv"))
}

##########################################
#Run markdown wrapper here

## Load outputs for printing
theme_set(theme_classic())
setwd("~/Dropbox/BethGlovers/2020 finfish")
lbbres<-list()
parsource<-c("Data","FL","None")
region
for(pars in 1:3) {
 priorval<-parsource[pars]
 filename<-paste0(priorval,region,"LBBouttable.csv")
 lbbres[[pars]]<-read.csv(filename)
} 
names(lbbres)<-parsource
summary(lbbres)
lbbres<-bind_rows(lbbres,.id="Priors")
 
lbbres$Priors[lbbres$Priors=="None"]<-"Default"
ggplot(filter(lbbres, !is.na(Linf.med)),aes(x=Stock,color=Priors))+
  geom_point(aes(y=Linf.med))+
  geom_errorbar(aes(ymin=Linf.lcl,ymax=Linf.ucl))+
  coord_flip()+ 
  xlab("Species")+ylab("Loo")  #+
#  ggtitle("Estimates of Loo from LBB when priors were data, \nFishLife meta-analysis and default priors")
ggplot(filter(lbbres, !is.na(MK.med)),aes(x=Stock,color=Priors))+
  geom_point(aes(y=MK.med))+
  geom_errorbar(aes(ymin=MK.lcl,ymax=MK.ucl))+
  coord_flip()+
  xlab("Species")+ylab("M/K") #+
#  ggtitle("Estimates of M/K from LBB when priors were data, \nFishLife meta-analysis and default priors")
ggplot(lbbres,aes(x=Linf.med,y=MK.med, color=Priors))+
  geom_point()+
  geom_errorbar(aes(xmin=Linf.lcl,xmax=Linf.ucl))+
  geom_errorbar(aes(ymin=MK.lcl,ymax=MK.ucl)) + 
  xlab("Loo")+
  ylab("M/K")

# Get data from other analysis to compare outputs
parsource<-c("Data","FL","LBBnoprior","LBBFLprior")
getres<-list()
for(i in 1:4) {
 filename<-paste0(region,parsource[i],"SingleIndicators.csv")
 getres[[i]]<-read.csv(filename)
} 
summary(getres)
names(getres)<-parsource
getres<-bind_rows(getres,.id="Source")
getreslong<-getres %>% pivot_longer(Pmat:SPR.F.M,names_to="Indicator",values_to="Value") %>%
  filter(Indicator %in% c("F.M","SPR","SPR.F.M")) %>%
  dplyr::select(Source,Species,Common,Indicator, Value) %>% 
  mutate(Code="All")
names(getreslong)

lbblong<-lbbres %>% dplyr::select(FM.med,BB0.med,Species,Stock,Priors) %>%
  rename(Common=Stock,Source=Priors) %>%
  pivot_longer(FM.med:BB0.med,names_to="Indicator",values_to="Value") %>%
  mutate(Code="LBB")
names(lbblong)
names(getreslong)

allres<-bind_rows(lbblong,getreslong) %>% filter(!is.na(Value))
summary(allres)
dim(allres)
table(allres$Indicator,allres$Source)

FMpars<-c("F.M","SPR.F.M","FM.med")
statpars<-c("SPR","BB0.med")

allres$Value[(allres$Indicator %in% FMpars) & (allres$Value>4)]<-4 

allres<-allres %>% mutate(Method=NA) %>%
  mutate(Method=ifelse(Indicator=="F.M",paste("Beverton/Holt",Source),Method),
         Method=ifelse(Indicator %in% c("SPR.F.M","SPR"),paste("LBSPR",Source),Method),
         Method=ifelse(Indicator %in% c("FM.med","BB0.med"),paste("LBB",Source),Method))
table(allres$Method)

ggplot(filter(allres,Indicator %in% FMpars),aes(x=Species,y=Value,color=Method))+
  geom_point()+coord_flip()+geom_hline(yintercept=1)+ylab("F/M")
ggplot(filter(allres,Indicator %in% statpars),aes(x=Species,y=Value,color=Method))+
  geom_point()+coord_flip()+geom_hline(yintercept=0.3) +ylab("Overfished status")

a<-match(allres$Species,spfileAll$Species)
summary(a)
length(unique(a))
length(a)
allres$Lmax<-spfileAll$Lmax[a]

allres$Trophic<-spfileAll$Trophic[a]

ggplot(filter(allres,Indicator %in% FMpars),aes(x=Lmax,y=Value))+
  geom_point(aes(color=Method))+
  geom_hline(yintercept=1)+
  stat_smooth(method="lm")+
  ylab("F/M")

ggplot(filter(allres,Indicator %in% statpars),aes(x=Lmax,y=Value))+
  geom_point(aes(color=Method))+
  geom_hline(yintercept=0.3)+
  stat_smooth(method="lm")+
  ylab("Overfished status")

methodstokeep<-c("LBB Data", "LBB FL" ,"LBB Default" , "Beverton/Holt Data", 
"LBSPR Data", "Beverton/Holt FL", "LBSPR FL", "Beverton/Holt LBBnoprior",
"LBSPR LBBnoprior")

sumtabF<-allres %>% filter(Indicator %in% FMpars & Method %in% methodstokeep) %>%
  mutate(Overfishing=ifelse(Value>1,1,0)) %>%
  dplyr::select(Species,Method,Overfishing) %>%
  pivot_wider(names_from = Method,values_from=Overfishing)
sumtabF$FractionF<-rowMeans(sumtabF[,-1],na.rm=TRUE)

sumtabB<-allres %>% filter(Indicator %in% statpars & Method %in% methodstokeep) %>%
  mutate(Overfished=ifelse(Value<0.3,1,0)) %>%
  dplyr::select(Species,Method,Overfished) %>%
  pivot_wider(names_from = Method,values_from=Overfished)
sumtabB$FractionB<-rowMeans(sumtabB[,-1],na.rm=TRUE)

sumtabB<-sumtabB[,sort(names(sumtabB))]
sumtabF<-sumtabF[,sort(names(sumtabF))]

sumtab<-merge(dplyr::select(sumtabF,Species,FractionF),
  dplyr::select(sumtabB,Species,FractionB),by="Species")


write.csv(sumtabB,paste0(region,"sumtabB.csv"))
write.csv(sumtabF,paste0(region,"sumtabF.csv"))


ggplot(sumtab,aes(x=FractionB,y=FractionF))+
  geom_jitter()+
  geom_abline(intercept=0,slope=1)

sumtab$Species[sumtab$FractionF>=7/8 & sumtab$FractionB>=7/8]
sumtab$Species[sumtab$FractionF<=1/8 & sumtab$FractionB<=1/8]

a<-match(sumtab$Species,spfileAll$Species)
sumtab$Lmax<-spfileAll$Lmax[a]
ggplot(pivot_longer(sumtab,FractionB:FractionF,values_to="Value",names_to="Indicator"),
  aes(x=Lmax,y=Value,color=Indicator,fill=Indicator))+
  geom_point()+
  stat_smooth(method="lm")+
  ylab("Fraction overfished or overfishing")

if(region=="Belize") {
  BelizeResults<-list(sumtab=sumtab,sumtabB=sumtabB,sumtabF=sumtabF,
    lbbres=lbbres,allres=allres,getres=getres)
}
if(region=="Guatemala") {
  GuatemalaResults<-list(sumtab=sumtab,sumtabB=sumtabB,sumtabF=sumtabF,
    lbbres=lbbres,allres=allres,getres=getres)
}

#Combine region after doing both

#M/K and Linf
lbbres<-bind_rows(list(Belize=BelizeResults$lbbres,Guatemala=GuatemalaResults$lbbres),.id="Country")
lbbres[,c("FM.med","FM.lcl","FM.ucl")][lbbres[,c("FM.med","FM.lcl","FM.ucl")]>4]<-4
summary(lbbres)
ggplot(filter(lbbres, !is.na(Linf.med)),aes(x=Stock,color=Priors))+
  geom_point(aes(y=Linf.med))+
  geom_errorbar(aes(ymin=Linf.lcl,ymax=Linf.ucl))+
  coord_flip()+ 
  facet_wrap(Country~.,nrow=2)+
  xlab("Species")+ylab("Loo")  #+
#  ggtitle("Estimates of Loo from LBB when priors were data, \nFishLife meta-analysis and default priors")
ggplot(filter(lbbres, !is.na(MK.med)),aes(x=Stock,color=Priors))+
  geom_point(aes(y=MK.med))+
  geom_errorbar(aes(ymin=MK.lcl,ymax=MK.ucl))+
  coord_flip()+
    facet_wrap(Country~.,nrow=2)+
  xlab("Species")+ylab("M/K") #+
#  ggtitle("Estimates of M/K from LBB when priors were data, \nFishLife meta-analysis and default priors")
ggplot(lbbres,aes(x=Linf.med,y=MK.med, color=Priors))+
  geom_point()+
  geom_errorbar(aes(xmin=Linf.lcl,xmax=Linf.ucl))+
  geom_errorbar(aes(ymin=MK.lcl,ymax=MK.ucl)) + 
  xlab("Loo")+
  ylab("M/K")+
  facet_wrap(Country~.,nrow=2,scales="free")


##
sumtab<-bind_rows(list(Belize=BelizeResults$sumtab,Guatemala=GuatemalaResults$sumtab),.id = "Country")
names(sumtab)
write.csv(sumtab,"senstivityS7.csv")
df2<-pivot_longer(sumtab,FractionF:FractionB,values_to="Value",names_to="Indicator")
ggplot(df2, aes(x=Lmax,y=Value))+
  geom_point(alpha=0.5)+
#  stat_smooth(method="lm")+
  ylab("Fraction overfished or overfishing")+
  facet_grid(rows = vars(Country),cols=vars(Indicator))+
  theme_bw()

df3<-df2 %>% arrange(Country,Indicator,Value)
write.csv(df3,"sensitivitySp.csv")

df3<-df2 %>% group_by(Country,Indicator) %>% 
    summarize(Total=n(),No=length(Value[Value==0])/Total,
    Yes=length(Value[Value==1])/Total,
      Uncertain=length(Value[Value>0 & Value<1])/Total    ) %>%
    mutate(Indicator=ifelse(Indicator=="FractionB","Overfished","Overfishing"))
df3
write.csv(df3,"senstivityB.csv")

df3<-df2 %>% group_by(Country,Indicator) %>% 
    summarize(SpeciesNo=toString(Species[Value==0]),
    SpeciesYes=toString(Species[Value==1]),
      Uncertain=toString(Species[Value>0 & Value<1]))   

write.csv(df3,"senstivityC.csv")
  
sumtab

ggplot(filter(df2,Lmax<150), aes(x=Lmax,y=Value,color=Indicator,fill=Indicator))+
  geom_point()+
  stat_smooth(method="lm")+
  ylab("Fraction overfished or overfishing")+
  facet_wrap(Country~.,nrow=2,scales="free")

##
head(sumtab)

sumtab %>% group_by(Country) %>% 
  summarize(Fagree=length(FractionF[FractionF>.8 | FractionF<0.2])/length(FractionF),
    Bagree=length(FractionB[FractionB>0.8 | FractionB<0.2])/length(FractionB)) 

names(sumtab)

### Table 4

allres<-bind_rows(list(Belize=BelizeResults$allres,Guatemala=GuatemalaResults$allres),.id="Country")

methodstokeep<-c("LBB Data", "LBB FL" ,"LBB Default" , "Beverton/Holt Data", 
"LBSPR Data", "Beverton/Holt FL", "LBSPR FL", "Beverton/Holt LBBnoprior",
"LBSPR LBBnoprior")


allres<-allres %>% filter(Method %in% methodstokeep) %>%
  mutate(Status=ifelse((Indicator %in% FMpars & Value>=1)|(Indicator %in% statpars & Value<=0.3),1,0 )) %>%
  mutate(Status=ifelse(is.na(Value),-1,Status)) %>%
  mutate(Source=ifelse(Source %in% c("Default","LBBnoprior"),"LBB",Source))
table(allres$Source)

sumtab<-allres %>% 
  group_by(Country,Source, Indicator) %>% 
  summarize(Mean=round(mean(Status,na.rm=TRUE),2))%>%
  pivot_wider(names_from="Source",values_from=Mean)

write.csv(sumtab,"sensitivity.csv")

sumtab$Indicator2<-rep(c("B/Bo","F/M (BH)","F/M (LBB)","SPR","F/M (SPR)"),2)
sumtablong<-pivot_longer(sumtab,Data:LBB,names_to="Parameters",values_to = "Fraction")
sumtablong$Indicator2<-factor(sumtablong$Indicator2,levels=c("F/M (BH)","F/M (SPR)","F/M (LBB)","SPR","B/Bo"))

ggplot(sumtablong,aes(x=Indicator2,y=Fraction,color=Parameters))+
  geom_point(size=2)+facet_wrap(Country~.)+
  xlab("Indicator")+ ylab("Fraction overfished or overfishing status")+
  theme_bw()


#figure out data source
x<-list(Belize=BelizeResults$sumtab,Guatemala=GuatemalaResults$sumtab)
x$Belize$Family<-BelizeAll$Family[match(x$Belize$Species,BelizeAll$Species)]
x$Guatemala$Family<-GuatemalaAll$Family[match(x$Guatemala$Species,GuatemalaAll$Species)]
x$Belize$Source<-BelizeSource$K[match(x$Belize$Species,BelizeSource$Species)]
x$Guatemala$Source<-GuatemalaSource$K[match(x$Guatemala$Species,GuatemalaSource$Species)]
x$Belize$lbb<-rep(0,nrow(x$Belize))
x$Belize$lbb[x$Belize$Species %in% lbbres$Species[lbbres$Country=="Belize" &!is.na(lbbres$FM.med) & lbbres$Priors=="Data"]]<-1
x$Guatemala$lbb<-rep(0,nrow(x$Guatemala))
x$Guatemala$lbb[x$Guatemala$Species %in% lbbres$Species[lbbres$Country=="Guatemala"&!is.na(lbbres$FM.med) & lbbres$Priors=="Data"]]<-1

sumtab<-bind_rows(x,.id = "Country") 
summary(factor(sumtab$Family))
table(sumtab$Country,sumtab$Source)
table(sumtab$Country,sumtab$lbb)

write.csv(sumtab,"Sup7.csv")

x<-sumtab[sumtab$lbb==1,]
table(x$Country,x$Source)

# Sources of data
table(BelizeSource$K[BelizeAll$n>=15])
table(GuatemalaSource$K[GuatemalaAll$n>=15])


## Methods separate
sumtabF<-bind_rows(list(Belize=BelizeResults$sumtabF,Guatemala=GuatemalaResults$sumtabF),.id = "Country")
sumtabB<-bind_rows(list(Belize=BelizeResults$sumtabB,Guatemala=GuatemalaResults$sumtabB),.id = "Country")

names(sumtabF)[2:11]<-paste("F.M",names(sumtabF)[2:11])
names(sumtabB)[2:8]<-paste("B",names(sumtabB)[2:8])
sumtabF$F.M.count<-rowSums(!is.na(sumtabF[,c(2:4,6:11)]))
sumtabB$B.count<-rowSums(!is.na(sumtabB[,c(3:8)]))
summary(sumtabB$B.count)
summary(sumtabF$F.M.count)

sumtabAll<-merge(sumtabF,sumtabB)
dim(sumtabAll)
sumtabAll$Family<-BelizeAll$Family[match(sumtabAll$Species,BelizeAll$Species)]
sumtabAll$Family[is.na(sumtabAll$Family)]<-GuatemalaAll$Family[match(sumtabAll$Species[is.na(sumtabAll$Family)],GuatemalaAll$Species)]
summary(factor(sumtabAll$Family))

write.csv(sumtabAll[,c("Country","Family","Species","F.M FractionF","F.M.count","B FractionB","B.count")],"sumtabAll.csv")


