#This file reads in the Belize data, does data cleaning and makes input files
# for further analysis and makes tables and figures. 

#Installation code
#devtools::install_github("james-thorson/FishLife")
#remotes::install_github( 'ropensci/rfishbase@fb-21.06', force=TRUE )

## Code to get rid of repeated warning
#!diagnostics off

#Load libraries
library(tidyverse)
library(ggforce)
library(gridExtra)
library( FishLife )
library(rfishbase)
library(TMB)
setwd("C:/Users/ebabcock/Dropbox/BethGlovers/2020 finfish")
source("babcockfunctionslength2020.r")

#Read in Belize data
x<-read.csv("Community Catch Data F17-J19 Fish ONLY 20190812.csv") #old
dim(x)
table(x$Year)
table(x$X..individuals)
y<-read.csv("Belize Community Catch Data 2017-2020.csv")
dim(y)
table(y$Year)
table(y$X..individuals)

Belize<-read.csv("Belize Community Catch Data 2017-2020.csv",as.is=TRUE) # Current
summary(Belize)

###Fix odd character if needed (encoding problem)
x<-grep("?",Belize$SCIENTIFIC.NAME2)
x<-grep("?",Belize$SCIENTIFIC.NAME2)
summary(x)
length(x)
Belize$SCIENTIFIC.NAME2[x]<-sub("?"," ",Belize$SCIENTIFIC.NAME2[x],fixed=TRUE)
x<-grep("?",Belize$SCIENTIFIC.NAME2)
summary(x)
length(x)

#Duplicate records for multiple individuals
Belize<-uncount(Belize,weights=X..individuals,.id="Uncount")
dim(Belize)  #19302

sort(names(Belize))
# Fix formatting of species names
Belize$SCIENTIFIC.NAME<-paste0(toupper(substring(Belize$SCIENTIFIC.NAME, 1, 1)), 
                               substring(Belize$SCIENTIFIC.NAME, 2))
Belize$SCIENTIFIC.NAME<-trimws(Belize$SCIENTIFIC.NAME)
Belize$SCIENTIFIC.NAME2<-paste0(toupper(substring(Belize$SCIENTIFIC.NAME2, 1, 1)), 
                                substring(Belize$SCIENTIFIC.NAME2, 2))
Belize$SCIENTIFIC.NAME2<-trimws(Belize$SCIENTIFIC.NAME2)
Belize$SPECIES<-trimws(Belize$SPECIES)
Belize$SPECIES2<-trimws(Belize$SPECIES2)
Belize$SPECIES2<-paste0(toupper(substring(Belize$SPECIES2, 1, 1)), 
                        tolower(substring(Belize$SPECIES2, 2)))
Belize$FAMILY<-trimws(Belize$FAMILY)
Belize$FAMILY<-paste0(toupper(substring(Belize$FAMILY, 1, 1)), 
                        tolower(substring(Belize$FAMILY, 2)))
table(Belize$FAMILY)

# Correct all variables

table(Belize$TROPHIC.GROUP)
Belize$Trophic2<-Belize$TROPHIC.GROUP
Belize$Trophic2[Belize$Trophic2=="Herbivore"]<-"Herbivore-Invertivore"
Belize$Gear<-Belize$UNKNOWN
Belize$Gear<-paste0(toupper(substring(Belize$Gear, 1, 1)), 
                        tolower(substring(Belize$Gear, 2)))
table(Belize$Gear)
Belize$Gear[Belize$Gear=="Line"]<-"Hand line"
Belize$Gear[Belize$Gear=="Other net"]<-"Beach weir"

ggplot(Belize,aes(x=TL2_CM))+geom_histogram()
summary(Belize$TL2_CM)
quantile(Belize$TL2_CM,na.rm=TRUE)
Belize$sizecat<-cut(Belize$TL2_CM,breaks=c(0,25,50,200))
table(Belize$sizecat)

table(Belize$Station)
Belize$Station[Belize$Station=="Caye Cualker"]<-"Caye Caulker"

Belize$Habitat<-rep("Coastal",dim(Belize)[1])

# Fix sailor's choice
table(Belize$SCIENTIFIC.NAME2[Belize$SPECIES2=="Sailors choice"])
Belize$SCIENTIFIC.NAME2[Belize$SPECIES2=="Sailors choice"]<-"Haemulon parra"


#### Make new column scinamefishbase with fb spellings, marine fish only
Belize$scinameFishbase<-Belize$SCIENTIFIC.NAME2
dim(Belize) #All rows 19302. 18776 without unknown gear
Belize <- filter(Belize,!scinameFishbase %in% c("Pangasius bocourti","Brycon guatemalensis","Vieja melanura","Oreochromis niloticus","Pterois volitans")) #Freshwater and lionfish
dim(Belize) #All rows 19058 without those species. 18533 without unknown gear

BelizeFam<-filter(Belize,!Gear=="Unknown")  #Correct database for analysis by family.
summary(factor(Belize$FAMILY))

#Fix odd character if needed (encoding problem)
x<-grep("?",Belize$scinameFishbase)
x<-grep("?",Belize$scinameFishbase)
summary(x)
length(x)


# Now filter Belize to species to species only
Belize$scinameFishbase[Belize$scinameFishbase=="Haemulon plumieri"]<-"Haemulon plumierii"
Belize$scinameFishbase[Belize$scinameFishbase=="Cephalopholis fulvus"]<-"Cephalopholis fulva"
Belize$scinameFishbase[Belize$scinameFishbase=="Etelis ocultus"]<-"Etelis oculatus"

Belize$scinameFishbase[grep("sp.",Belize$scinameFishbase,fixed=TRUE)]=NA  #Remove fish id'ed to genus
table(Belize$SCIENTIFIC.NAME2[is.na(Belize$scinameFishbase)])

#Take out records for fish that are not to species
Belize<-Belize[!is.na(Belize$scinameFishbase),]
dim(Belize)  #18907  
table(Belize$scinameFishbase)
BelizeGear<-filter(Belize,!Gear=="Unknown")
dim(BelizeGear) #18382

#Identify most common species
x<-sort(table(Belize$scinameFishbase))
x<-x[x>=30]
x
length(x)
Belizesp<-data.frame(Species=names(x))
for(i in 1:length(x)) {
  Belizesp$Genus[i]=unlist(strsplit(Belizesp$Species[i]," "))[1]
  Belizesp$species[i]=unlist(strsplit(Belizesp$Species[i]," "))[2]
}
x<-match(Belizesp$Species,Belize$scinameFishbase)
summary(x)
Belizesp$Common<-Belize$SPECIES2[x]
dim(Belizesp)

#FLTL conversions for common species
Belizesp$tlfl.a<-rep(NULL,length(Belizesp$Species))
Belizesp$tlfl.b<-rep(NULL,length(Belizesp$Species))
Belizesp$tlfl.r2<-rep(NULL,length(Belizesp$Species))
Belizesp$tlfl.n<-rep(NULL,length(Belizesp$Species))
g<-list()
for(i in 1:length(Belizesp$Species)) {
  x<-filter(Belize,Belize$Uncount==1 & !is.na(Belize$FL2_CM) & !is.na(Belize$TL2_CM)& Belize$scinameFishbase==Belizesp$Species[i])
  Belizesp$tlfl.n[i]<-dim(x)[1]
  if(length(which(x$FL2_CM==x$TL2_CM))== Belizesp$tlfl.n[i]) {
   Belizesp$tlfl.a[i]<-0
   Belizesp$tlfl.b[i]<-1
   Belizesp$tlfl.r2[i]<-NA
  } else   {
  y<-lm(TL2_CM~FL2_CM,data=x)
  g[[i]]<-ggplot(x,aes(x=FL2_CM,y=TL2_CM))+geom_point()+stat_smooth(method="lm")+ggtitle(paste(i,Belizesp$Species[i]))
  Belizesp$tlfl.a[i]<-round(coef(y)[1],3)
  Belizesp$tlfl.b[i]<-round(coef(y)[2],3)
  Belizesp$tlfl.r2[i]<-round(summary(y)$r.squared,2)
  }
}
write.csv(Belizesp,"FLTL.csv")

# Get fishbase data for all  
x<-sort(unique(Belize$scinameFishbase))
y<-validate_names(x)
length(x)
length(y)
y[!(y %in% x)]
x[!(x%in%y)]

BelizeSpecies<-species(x)
BelizeEcology<-ecology(x) 
BelizeAll<-dplyr::select(BelizeSpecies,Species,FBname)
BelizeSource<-expand.grid(Species=BelizeAll$Species,Lm=NA,Linf=NA,K=NA,tmax=NA,Lmax=NA)
table(BelizeSpecies$Saltwater)
BelizeGrowth<-popgrowth(x)
BelizePopChar<-popchar(x)
summary(BelizePopChar)
BelizePopChar$Lmax<-as.numeric(BelizePopChar$Lmax)
BelizeGrowth$Loo[BelizeGrowth$Species=="Melichthys niger" &BelizeGrowth$Loo==2037]<-20.37  #Typo in Fishbase, found correct value in Kavanagh 2006
BelizeMat<-maturity(x)
y<-BelizeMat %>% group_by(Species) %>%
  summarize(Lm=median(Lm,na.rm=TRUE),LengthMatMin=median(LengthMatMin,na.rm=TRUE))
z<-y$Species[is.na(y$Lm &!is.na(y$LengthMatMin))]
BelizeSource$Lm[!is.na(y$Lm)]<-"FB"
BelizeSource$Lm[is.na(y$Lm &!is.na(y$LengthMatMin))]<-"FB.LMatMin"
BelizeMat$Lm[BelizeMat$Species %in%z]<-BelizeMat$LengthMatMin[BelizeMat$Species %in%z]

#Length conversions
Belizepopll<-popll(x) 
Belizepopll$a<-as.numeric(Belizepopll$a)
Belizepopll$b<-as.numeric(Belizepopll$b)
BelizepopllFL<-Belizepopll %>% filter(Length1=="TL" & Length2=="FL") %>%
  filter(!duplicated(Species)) 
BelizepopllSL<-Belizepopll %>% filter(Length1=="TL" & Length2=="SL") %>%
  filter(!duplicated(Species)) 
summary(BelizepopllFL)
x<-merge(Belizesp,BelizepopllFL[,c("Species","a","b")],by="Species",all.x=TRUE)
dim(x)
ggplot(x,aes(x=a,y=tlfl.a))+geom_point()+geom_abline(intercept=0,slope=1)
ggplot(x,aes(x=b,y=tlfl.b))+geom_point()+geom_abline(intercept=0,slope=1)
x<-merge(BelizeAll,BelizepopllFL[,c("Species","a","b")],all.x=TRUE,by="Species")
x<-merge(x,Belizesp[,c("Species","tlfl.a","tlfl.b")],all.x=TRUE,by="Species")
head(x)
summary(x)
y<-match(Belize$scinameFishbase,x$Species)
summary(y)
z<-Belize
z$calcTLfb<-x$a[y]+x$b[y]* z$FL2_CM
z$calcTLdat<-x$tlfl.a[y]+x$tlfl.b[y]* z$FL2_CM
ggplot(z,aes(x=TL2_CM,y=calcTLfb))+geom_point()+geom_abline(intercept=0,slope=1)+
  geom_point(aes(x=calcTLdat),col="red")
ggplot(z,aes(x=calcTLfb/TL2_CM))+geom_histogram()
# It looks like conversion are fine for either fb or data calculations. 
# Will use our numbers, except where not available. 
BelizeAll<-merge(BelizeAll,BelizepopllFL[,c("Species","a","b"),],all.x=TRUE,by="Species")
x<-match(BelizeAll$Species,Belizesp$Species)
summary(x)
BelizeAll$a<-ifelse(!is.na(x),Belizesp$tlfl.a[x],BelizeAll$a)
BelizeAll$b<-ifelse(!is.na(x),Belizesp$tlfl.b[x],BelizeAll$b)

BelizepopllSL$SL.a<-BelizepopllSL$a
BelizepopllSL$SL.b<-BelizepopllSL$b
BelizeAll<-merge(BelizeAll,BelizepopllSL[,c("Species","SL.a","SL.b"),],all.x=TRUE,by="Species")
BelizeAll$SL.a[BelizeAll$Species=="Sphyraena barracuda"]<-0
BelizeAll$SL.b[BelizeAll$Species=="Sphyraena barracuda"]<-round(1/0.828,2)  # From SL~TL in Fishbase
BelizeAll$SL.a[BelizeAll$Species=="Micropogonias undulatus"]<-0
BelizeAll$SL.b[BelizeAll$Species=="Micropogonias undulatus"]<-round(1/0.8,2)  # From SL~TL in Fishbase

#a and b in BelizeAll are now FL to TL conversions, SL.a and SL.b convert from SL

# Now get Lm in TL
table(BelizeMat$Type1)
x<-match(BelizeMat$Species,BelizeAll$Species)
summary(x)
y<-BelizeAll$a[x]+BelizeAll$b[x]* BelizeMat$Lm
y[!BelizeMat$Type1=="FL"]<-NA
z<-BelizeAll$SL.a[x]+BelizeAll$SL.b[x]* BelizeMat$Lm
z[!BelizeMat$Type1=="SL"]<-NA
BelizeMat$TL.Lm<-ifelse(BelizeMat$Type1=="FL",y,BelizeMat$Lm)
BelizeMat$TL.Lm<-ifelse(BelizeMat$Type1=="SL",z,BelizeMat$TL.Lm)
BelizeMat$TL.Lm<-ifelse(is.na(BelizeMat$Type1),BelizeMat$Lm,BelizeMat$TL.Lm)
table(BelizeMat$Species,is.na(BelizeMat$TL.Lm))
BelizeMat$Lm<-BelizeMat$TL.Lm

BelizeMat<-dplyr::select(BelizeMat,Species,Type1,Lm)
summary(BelizeMat)
#Loo
x<-match(BelizeGrowth$Species,BelizeAll$Species)
y<-BelizeAll$a[x]+BelizeAll$b[x]* BelizeGrowth$Loo
y[!BelizeGrowth$Type=="FL"]<-NA
BelizeGrowth$TL.Loo<-ifelse(BelizeGrowth$Type=="FL",y,BelizeGrowth$Loo)
z<-BelizeAll$SL.a[x]+BelizeAll$SL.b[x]* BelizeGrowth$Loo
z[!BelizeGrowth$Type=="SL"]<-NA
BelizeGrowth$TL.Loo<-ifelse(BelizeGrowth$Type=="SL",z,BelizeGrowth$TL.Loo)
BelizeGrowth$TL.Loo<-ifelse(is.na(BelizeGrowth$Type),BelizeGrowth$Loo,BelizeGrowth$TL.Loo)
BelizeGrowth<-BelizeGrowth %>%  rename(Linf=TL.Loo)
#summary(factor(BelizeGrowth$Type[is.na(BelizeGrowth$TL.Loo) & !is.na(BelizeGrowth$Loo)]))
summary(factor(BelizeGrowth$Species[is.na(BelizeGrowth$TL.Loo) & !is.na(BelizeGrowth$Loo)]))
#Lmax convert to TL in popchar and Species
x<-match(BelizePopChar$Species,BelizeAll$Species)
summary(x)
y<-BelizeAll$a[x]+BelizeAll$b[x]* BelizePopChar$Lmax
y[!BelizePopChar$Type=="FL"]<-NA
z<-BelizeAll$SL.a[x]+BelizeAll$SL.b[x]* BelizePopChar$Lmax
z[!BelizePopChar$Type=="SL"]<-NA
z[BelizePopChar$Species=="Sphyraena barracuda"]<-NA #Leave out SL barracuda because doesn't convert well
BelizePopChar<- BelizePopChar %>%
  mutate(Lmax.TL=ifelse(Type=="FL",y,Lmax)) %>%
  mutate(Lmax.TL=ifelse(Type=="SL",z,Lmax.TL)) %>%
  mutate(Lmax=Lmax.TL)

# New recent data that is not in fishbase yet
vega<-read.csv("VegaCendejas.csv")
summary(vega)
names(vega)[1]<-"Species"
x<-match(vega$Species,BelizeAll$Species)
table(is.na(x))
vega<-filter(vega,!is.na(x))
x<-match(vega$Species,BelizeAll$Species)
BelizeAll$Species[x]
x<-match(BelizeAll$Species,Belize$scinameFishbase)
BelizeAll$Family<-Belize$Family[x,]
x<-match(vega$Species,BelizeAll$Species)
y<- BelizeAll$SL.a[x]+BelizeAll$SL.b[x]*vega$SLm
vega$Lm<-ifelse(is.na(vega$TLm),y,vega$TLm)
summary(vega$Lm)
y<- BelizeAll$SL.a[x]+BelizeAll$SL.b[x]*vega$SLmax
vega$Lmax<-ifelse(is.na(vega$TLmax),y,vega$TLmax)
summary(vega$Lmax)
vega$Lm<-round(vega$Lm,2)
vega  #Output with converted lengths and added to newdat
kadison<-c(Lmax.FL=134.8,Linf.FL=123.64)
round(BelizeAll$a[BelizeAll$Species=="Sphyraena barracuda"]+BelizeAll$b[BelizeAll$Species=="Sphyraena barracuda"]*kadison,2)
# Vega and Kadison are now in newdat
newdat<-read.csv("newdatBelize.csv")
names(newdat)[1]<-"Species"

BelizeGrowth<-bind_rows(list(FB=BelizeGrowth,New=dplyr::select(newdat,Species,K,Linf,tmax)),.id="id")
dim(BelizeGrowth)

BelizePopChar<-bind_rows(list(FB=dplyr::select(BelizePopChar,Species,Lmax,tmax),New=dplyr::select(newdat,Species,Lmax,tmax)),.id="id")
dim(BelizePopChar)

BelizeMat<-bind_rows(list(FB=dplyr::select(BelizeMat,Species,Lm),New=dplyr::select(newdat,Species,Lm)),.id="id")
dim(BelizeMat)

#Calculate median Lm, K, Linf, K and max tmax and Lmax
BelizeGrowthSum<-BelizeGrowth %>% group_by(Species) %>%
  summarize(Linf=median(Linf,na.rm=TRUE),
            K=median(K,na.rm=TRUE),tmax=max(tmax,na.rm=TRUE)) %>%
  mutate(tmax=ifelse(tmax<0,NA,tmax))
dim(BelizeGrowthSum)
summary(BelizeGrowthSum)

BelizePopCharSum<-BelizePopChar %>% group_by(Species) %>%
  summarize(Lmax=max(Lmax,na.rm=TRUE),tmax=max(tmax,na.rm=TRUE)) %>%
  mutate(Lmax=ifelse(Lmax<0,NA,Lmax),tmax=ifelse(tmax<0,NA,tmax))
summary(BelizePopCharSum)

BelizeMatSum<-BelizeMat %>% group_by(Species) %>%
  summarize(Lm=median(Lm,na.rm=TRUE)) %>%
  mutate(Lm=ifelse(Lm<0,NA,Lm))
summary(BelizeMatSum)
dim(BelizeMatSum)

#Get max Lmax
x<-BelizeAll$a+BelizeAll$b*BelizeSpecies$Length
x[!BelizeSpecies$LTypeMaxM=="FL"]<-NA
BelizeSpecies$Length.TL<-ifelse(BelizeSpecies$LTypeMaxM=="TL",BelizeSpecies$Length,x)
table(is.na(BelizeSpecies$Length.TL),is.na(BelizePopCharSum$Lmax))
BelizePopCharSum<-BelizePopCharSum %>% mutate(Lmax=ifelse(is.na(Lmax),
  BelizeSpecies$Length.TL,Lmax)) %>%
  mutate(Lmax=ifelse(!is.na(Lmax) & Lmax<BelizeSpecies$Length.TL,
    BelizeSpecies$Length.TL,Lmax))
summary(BelizePopCharSum[,c("Lmax","tmax")])

# Get max tmax
BelizeGrowthSum$tmax<-ifelse(is.na(BelizeGrowthSum$tmax),BelizePopCharSum$tmax,
  BelizeGrowthSum$tmax)
BelizeGrowthSum$tmax<-ifelse(!is.na(BelizeGrowthSum$tmax+BelizePopCharSum$tmax) &
    BelizePopCharSum$tmax>BelizeGrowthSum$tmax,BelizePopCharSum$tmax,BelizeGrowthSum$tmax)
summary(BelizeGrowthSum$tmax)
BelizeGrowthSum<-merge(BelizeGrowthSum,BelizeMatSum,by="Species",all=TRUE)
BelizeGrowthSum$Lmax<-BelizePopCharSum$Lmax
summary(BelizeGrowthSum)
BelizeSource$Linf[!is.na(BelizeGrowthSum$Linf)]<-"FB"
BelizeSource$K[!is.na(BelizeGrowthSum$K)]<-"FB"
BelizeSource$Lmax[!is.na(BelizeGrowthSum$Lmax)]<-"FB"
BelizeSource$tmax[!is.na(BelizeGrowthSum$tmax)]<-"FB"
BelizeSource$Lm[!is.na(BelizeGrowthSum$Lm) & is.na(BelizeSource$Lm)]<-"FB"

# Combine tables
BelizeAll<-merge(BelizeAll,BelizeGrowthSum) 
dim(BelizeAll)

#Check if Lmax in fishbase is larger than max measured length
y<-Belize %>% group_by(scinameFishbase) %>%
  summarize(LmaxDat=max(TL2_CM,na.rm=TRUE)) %>%
  rename(Species=scinameFishbase)
BelizeAll<-merge(BelizeAll,y,by="Species")
table(BelizeAll$Lmax>BelizeAll$LmaxDat)
x<-BelizeAll$Species[BelizeAll$LmaxDat > BelizeAll$Lmax *1.1 ]
df1<-filter(Belize, scinameFishbase %in% x)
df1$Lmax<-BelizeAll$Lmax[match(df1$scinameFishbase,BelizeAll$Species)]
ggplot(df1,aes(x=TL2_CM))+
  geom_histogram(fill="blue")+
  facet_wrap(scinameFishbase~.,scales="free",ncol=2) +
  geom_vline(aes(xintercept=Lmax),col="red")
BelizeSource$Lmax[BelizeAll$LmaxDat>BelizeAll$Lmax]<-"Data"
BelizeAll$Lmax<-ifelse(BelizeAll$LmaxDat>BelizeAll$Lmax,BelizeAll$LmaxDat,BelizeAll$Lmax)

#Compare to data from Stevens, M. H., S. G. Smith, and J. S. Ault. 2019. 
#Life history demographic parameter synthesis for exploited florida and caribbean 
#coral reef fishes. Fish and Fisheries 20:1196-1217.
stevens<-read.csv("stevens.csv")
summary(stevens)
names(stevens)[1]<-"Common"
x<-match(stevens$Species,BelizeAll$Species)
table(is.na(x))  
stevens$Common[is.na(x)]
stevens$Linf<-round(BelizeAll$a[x]+BelizeAll$b[x]*stevens$Linf.FL/10,2)
stevens$Common[is.na(stevens$Linf) & !is.na(x)]
table(is.na(stevens$Linf),is.na(stevens$Linf.FL))
table(is.na(stevens$Linf) & is.na(stevens$Linf.FL),is.na(x))
stevens$Common[is.na(stevens$Linf) & !is.na(stevens$Linf.FL) &!is.na(x)]
stevens$Lm<-round(BelizeAll$a[x]+BelizeAll$b[x]*stevens$Lm.FL/10,2)
stevens$Common[is.na(stevens$Lm) &!is.na(stevens$Lm.FL)& !is.na(x)]
stevens$Lmax<-round(BelizeAll$a[x]+BelizeAll$b[x]*stevens$Lmax.FL/10,2)
stevens$Common[is.na(stevens$Lmax) &!is.na(stevens$Lmax.FL)& !is.na(x)]
x<-match(BelizeAll$Species,stevens$Species)
table(is.na(x))

df1<-merge(dplyr::select(BelizeAll,Species,Linf,K,Lmax,tmax,Lm),dplyr::select(stevens,Species,Linf,K,Lmax,tmax,Lm),by="Species")
names(df1)
g1<-ggplot(df1,aes(x=Linf.x,y=Linf.y))+geom_point()+geom_abline(aes(intercept=0,slope=1))
g2<-ggplot(df1,aes(x=K.x,y=K.y))+geom_point()+geom_abline(aes(intercept=0,slope=1))
g3<-ggplot(df1,aes(x=tmax.x,y=tmax.y))+geom_point()+geom_abline(aes(intercept=0,slope=1))
g4<-ggplot(df1,aes(x=Lmax.x,y=Lmax.y))+geom_point()+geom_abline(aes(intercept=0,slope=1))
g5<-ggplot(df1,aes(x=Lm.x,y=Lm.y))+geom_point()+geom_abline(aes(intercept=0,slope=1))
grid.arrange(g1,g2,g3,g4,g5)

x<-match(BelizeAll$Species,stevens$Species)
summary(x)

BelizeAll$FBname[is.na(BelizeAll$tmax)&!is.na(stevens$tmax[x])]
BelizeSource$tmax[is.na(BelizeAll$tmax)&!is.na(stevens$tmax[x])]<-"Stevens"
BelizeAll$tmax<-ifelse(is.na(BelizeAll$tmax)&!is.na(stevens$tmax[x]),stevens$tmax[x],BelizeAll$tmax)
BelizeAll$FBname[!is.na(BelizeAll$tmax)&!is.na(stevens$tmax[x])&
    stevens$tmax[x]>BelizeAll$tmax]
BelizeSource$Stevens[!is.na(BelizeAll$tmax)&!is.na(stevens$tmax[x])&
    stevens$tmax[x]>BelizeAll$tmax]
BelizeAll$tmax<-ifelse(!is.na(BelizeAll$tmax)&!is.na(stevens$tmax[x])&
    stevens$tmax[x]>BelizeAll$tmax,stevens$tmax[x],BelizeAll$tmax)
summary(BelizeAll$tmax)
BelizeAll$FBname[is.na(BelizeAll$Lmax)&!is.na(stevens$Lmax[x])]
#No records
BelizeAll$FBname[!is.na(BelizeAll$Lmax)&!is.na(stevens$Lmax[x])&
    stevens$Lmax[x]>BelizeAll$Lmax]
BelizeSource$Lmax[!is.na(BelizeAll$Lmax)&!is.na(stevens$Lmax[x])&
    stevens$Lmax[x]>BelizeAll$Lmax]<-"Stevens"
BelizeAll$Lmax<-ifelse(!is.na(BelizeAll$Lmax)&!is.na(stevens$Lmax[x])&
    stevens$Lmax[x]>BelizeAll$Lmax,stevens$Lmax[x],BelizeAll$Lmax)
summary(BelizeAll$Lmax)
BelizeAll$FBname[is.na(BelizeAll$Linf)&!is.na(stevens$Linf[x])]
#NO records
BelizeAll$FBname[is.na(BelizeAll$Lm)&!is.na(stevens$Lm[x])]
BelizeSource$Lm[is.na(BelizeAll$Lm)&!is.na(stevens$Lm[x])]<-"Stevens"
BelizeAll$Lm<-ifelse(is.na(BelizeAll$Lm)&!is.na(stevens$Lm[x]),stevens$Lm[x],
  BelizeAll$Lm)

BelizeSource$Species[!is.na(stevens$Lm[x]) & BelizeSource$Lm=="FB.LMatMin" ]
BelizeSource$Lm[!is.na(stevens$Lm[x]) & BelizeSource$Lm=="FB.LMatMin" ]<-"Stevens"

BelizeAll$Lm[!is.na(stevens$Lm[x]) & BelizeSource$Lm=="FB.LMatMin"]<-
  stevens$Lm[!is.na(stevens$Lm[x]) & BelizeSource$Lm=="FB.LMatMin"]
summary(BelizeAll$Lm)
table(BelizeSource$Lm)

#Trophic
BelizeAll$Trophic
BelizeAll$Trophic<-ifelse(is.na(BelizeEcology$DietTroph),BelizeEcology$FoodTroph,BelizeEcology$DietTroph)
summary(BelizeAll[,c("Lm","Linf","Lmax","K","tmax","Trophic")])
BelizeAll[is.na(BelizeAll$Trophic),]
BelizeSource$Trophic<-ifelse(is.na(BelizeAll$Trophic),"fb.genus","fb")
a<-which(is.na(BelizeAll$Trophic))
for(i in a) {
  x<-species_list(Genus=local(BelizeAll$Genus[i]))
  y<-ecology(x)
  y$Trophic<-ifelse(is.na(y$DietTroph),y$FoodTroph,y$DietTroph)
  BelizeAll$Trophic[i]<-median(y$Trophic,na.rm=TRUE)
}
summary(BelizeAll$Trophic)
summary(BelizeAll[,c("Lm","Linf","Lmax","K","tmax","Trophic")])

# Common species without data
BelizeAll$Species[is.na(BelizeAll$tmax) & BelizeAll$Species %in% Belizesp$Species]
BelizeAll$FBname[is.na(BelizeAll$tmax) & BelizeAll$Species %in% Belizesp$Species]
BelizeAll$Species[is.na(BelizeAll$Linf) & BelizeAll$Species %in% Belizesp$Species]
BelizeAll$Species[is.na(BelizeAll$Lm) & BelizeAll$Species %in% Belizesp$Species]

#Sample size
x<-Belize %>% group_by(scinameFishbase) %>% 
  summarize(n=length(scinameFishbase)) %>%
  rename(Species=scinameFishbase)
BelizeAll<-merge(BelizeAll,x,by="Species")
names(BelizeAll)

#Add genus and species to All
for(i in 1:length(BelizeAll$Species)) {
  BelizeAll$Genus[i]=unlist(strsplit(BelizeAll$Species[i]," "))[1]
  BelizeAll$species[i]=unlist(strsplit(BelizeAll$Species[i]," "))[2]
}

#Select only needed columns in All
names(BelizeAll)
x<-match(BelizeAll$Species,Belize$scinameFishbase)
BelizeAll$Common<-Belize$SPECIES2[x]
summary(BelizeAll)

#Run fishlife for all
Belizepars<-list()
BelizeAll$species[BelizeAll$Species=="Diplodus argenteus"]="argenteus argenteus"
BelizeAll$predictive<-rep(1,dim(BelizeAll)[1])
for(i in 1:dim(BelizeAll)[1]) {
 x <-try(Plot_taxa( Search_species(Genus=BelizeAll$Genus[i],Species=BelizeAll$species[i])$match_taxonomy ))
 Belizepars[[i]]<-x[[1]]
 if(length(x)==5) BelizeAll$predictive[i]=0
}
summary(BelizeAll)
summary(Belizepars)

params<-c("Lm","M","K","Loo")
Belizeparmean<-NULL
Belizeparmedian<-NULL
Belizeparse<-NULL
for(i in 1:dim(BelizeAll)[1]) {
  meanvalMK<-Belizepars[[i]]$Mean_pred["M"]-Belizepars[[i]]$Mean_pred["K"]
  names(meanvalMK)<-"MK"
  varvalMK<-Belizepars[[i]]$Cov_pred["M","M"]+Belizepars[[i]]$Cov_pred["K","K"]+2*Belizepars[[i]]$Cov_pred["M","K"]
  names(varvalMK)<-"MK"
  meanval<-lnorm.mean(c(Belizepars[[i]]$Mean_pred,meanvalMK),sqrt(c(diag(Belizepars[[i]]$Cov_pred),varvalMK)))
  seval<-lnorm.se(c(Belizepars[[i]]$Mean_pred,meanvalMK),sqrt(c(diag(Belizepars[[i]]$Cov_pred),varvalMK)))
  Belizeparmean<-rbind(Belizeparmean,meanval)
  Belizeparse<-rbind(Belizeparse,seval)
  Belizeparmedian<-rbind(Belizeparmedian,exp(c(Belizepars[[i]]$Mean_pred,meanvalMK)))
}
params<-c("Lm","M","K","Loo","MK")
Belizepardf<-data.frame(cbind(Belizeparmean[,params],Belizeparse[,params]),Belizeparmedian[,params])
names(Belizepardf)
names(Belizepardf)[6:10]<-paste0(params,"se")
names(Belizepardf)[11:15]<-paste0(params,"med")

Belizepardf$Species<-BelizeAll$Species
df2<-merge(BelizeAll,Belizepardf,by="Species")
names(df2)
df2$Linf.Lmax<-ifelse(is.na(df2$Linf),Linf.Lmax(BelizeAll$Lmax),df2$Linf)
df2$Lm.Linf<-ifelse(is.na(df2$Lm.x),Lm.Linf(df2$Linf.Lmax),df2$Lm.x)
df2$M<-ifelse(is.na(df2$tmax),NA,M.Then(df2$tmax))
g1<-ggplot(df2,aes(x=Lm.x,y=Lm.y,color=factor(predictive)))+geom_point()+geom_abline(intercept=0,slope=1)+
  geom_errorbar(aes(ymin=Lm.y-2*Lmse,ymax=Lm.y+2*Lmse))+ggtitle("Lm")+ theme(legend.position = "none") 
g2<-ggplot(df2,aes(x=Lm.Linf,y=Lm.y,color=factor(predictive)))+geom_point()+geom_abline(intercept=0,slope=1)+
  geom_errorbar(aes(ymin=Lm.y-2*Lmse,ymax=Lm.y+2*Lmse))+ggtitle("Lm.Linf")+ theme(legend.position = "none") 
g3<-ggplot(df2,aes(x=Linf,y=Loo,color=factor(predictive)))+geom_point()+geom_abline(intercept=0,slope=1)+
  geom_errorbar(aes(ymin=Loo-2*Loose,ymax=Loo+2*Loose))+ggtitle("Loo")+ theme(legend.position = "none") 
g4<-ggplot(df2,aes(x=Linf.Lmax,y=Loo,color=factor(predictive)))+geom_point()+geom_abline(intercept=0,slope=1)+
  geom_errorbar(aes(ymin=Loo-2*Loose,ymax=Loo+2*Loose))+ggtitle("Linf.Lmax")+ theme(legend.position = "none") 
g5<-ggplot(df2,aes(x=K.x,y=K.y,color=factor(predictive)))+geom_point()+geom_abline(intercept=0,slope=1)+
  geom_errorbar(aes(ymin=K.y-2*Kse,ymax=K.y+2*Kse))+ggtitle("K")+ theme(legend.position = "none") 
g6<-ggplot(df2,aes(x=M.Then(tmax),y=M,color=factor(predictive)))+geom_point()+geom_abline(intercept=0,slope=1)+
  geom_errorbar(aes(ymin=M-2*Mse,ymax=M+2*Mse))+ggtitle("M")+ theme(legend.position = "none") +
  ylim(c(0,1.2))
g7<-ggplot(df2,aes(x=M/K.x,y=MK,color=factor(predictive)))+geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_errorbar(aes(ymin=MK-2*MKse,ymax=MK+2*MKse))+
  ggtitle("M/K")+ theme(legend.position = "none") +ylim(c(-1,6))
grid.arrange(g1,g2,g3,g4,g5,g6,g7,ncol=2)

##
summary(Belizepardf$MKse/Belizepardf$MK)
summary(Belizepardf$MK)
ggplot(Belizepardf,aes(x=MKmed,y=Mmed/Kmed))+geom_point()+geom_abline(intercept=0,slope=1)
ggplot(Belizepardf,aes(x=MKmed,y=MK))+geom_point()+geom_abline(intercept=0,slope=1)

#Fix too large Lmax
BelizeSource$tmax[BelizeAll$Species=="Caranx lugubris"]<-NA
BelizeAll$tmax[BelizeAll$Species=="Caranx lugubris"]
BelizeAll$tmax[BelizeAll$Species=="Caranx lugubris"]=NA #This number is from the Pacific and seems off

#See which parameters are still missing
summary(BelizeAll)
BelizeAll$Common[is.na(BelizeAll$Lm) & BelizeAll$n>30]
BelizeAll$Common[is.na(BelizeAll$tmax)  & BelizeAll$n>30]
BelizeAll$Common[is.na(BelizeAll$K) & BelizeAll$n>30]

# add all parameters to BelizeAll
BelizeSource$Linf[is.na(BelizeAll$Linf)]<-"Lmax"
BelizeAll$Linf<-ifelse(is.na(BelizeAll$Linf) ,Linf.Lmax(BelizeAll$Lmax),BelizeAll$Linf)
BelizeSource$Lm[is.na(BelizeAll$Lm)]<-"Linf"
BelizeAll$Lm<-ifelse(is.na(BelizeAll$Lm),Lm.Linf(BelizeAll$Linf),BelizeAll$Lm)
BelizeSource$M<-ifelse(is.na(BelizeAll$tmax),NA,"Then")
BelizeSource$M<-ifelse(is.na(BelizeAll$tmax) &!is.na(BelizeAll$K),
  "Then.K",BelizeSource$M)
x<-ifelse(!is.na(M.Then(BelizeAll$tmax)),M.Then(BelizeAll$tmax),M.Then.K(BelizeAll$Linf,BelizeAll$K))
summary(x)
BelizeSource$M[is.na(x)]<-"FishLife"
BelizeAll$M<-ifelse(is.na(x), Belizeparmedian[,"M"],x)
BelizeSource$K[is.na(BelizeAll$K)]<-"Fishlife"
BelizeAll$K<-ifelse(is.na(BelizeAll$K),Belizeparmedian[,"K"],BelizeAll$K)
summary(BelizeAll[,c("Lm","M","K","Linf")])

Belizepardf<-Belizepardf %>% mutate(LooCV=Loose/Loo,KCV=Kse/K,MCV=Mse/M,LmCV=Lmse/Lm)
summary(Belizepardf)

# Print source table
BelizeSource$Predictive<-BelizeAll$predictive
BelizeSource$n<-BelizeAll$n
BelizeSource$Common<-BelizeAll$Common
x<-match(BelizeAll$Species,Belize$scinameFishbase)
BelizeAll$Family<-Belize$FAMILY[x]
BelizeSource$Family<-Belize$FAMILY[x]
write.csv(BelizeSource[order(BelizeSource$Family),],"BelizeSource.csv")
  
#Add fishbase numbers to Belize
x<-match(Belize$scinameFishbase,BelizeAll$Species)
Belize$Lmax<-BelizeAll$Lmax[x]
Belize$Lm<-BelizeAll$Lm[x]
Belize$Trophic<-BelizeAll$Trophic[x]


# Check life history values
g1<-ggplot(BelizeAll,aes(x=Linf/Lmax))+geom_histogram()
g2<-ggplot(BelizeAll,aes(x=Lm/Linf))+geom_histogram()
g3<-ggplot(BelizeAll,aes(x=M/K))+geom_histogram()
g4<-ggplot(BelizeAll,aes(x=Lmax,y=M))+geom_point()
grid.arrange(g1,g2,g3,g4)

x<-filter(BelizeAll,M/K>4 | M/K <0.5)
x$M.K<-x$M/x$K
dim(x)
x
#Dolphinfish has a high growth and low M. This is okay

#Compare to what we used before
fishLH<-read.csv("FishLH2020.csv")
dim(fishLH)
x<-match(fishLH$Sci.name,BelizeAll$Species)
summary(x)
fishLH[is.na(x),]
fishLH<-fishLH[!is.na(x),]
df1<-merge(fishLH,BelizeAll[x,],by.x="Sci.name",by.y="Species")
g1<-ggplot(df1,aes(x=Linf.x,y=Linf.y))+geom_point()+geom_abline(aes(intercept=0,slope=1))
g2<-ggplot(df1,aes(x=K.x,y=K.y))+geom_point()+geom_abline(aes(intercept=0,slope=1))
g3<-ggplot(df1,aes(x=Lm.x,y=Lm.y))+geom_point()+geom_abline(aes(intercept=0,slope=1))
g4<-ggplot(df1,aes(x=Lmax.x,y=Lmax.y))+geom_point()+geom_abline(aes(intercept=0,slope=1))
g5<-ggplot(df1,aes(x=tmax.x,y=tmax.y))+geom_point()+geom_abline(aes(intercept=0,slope=1))
g6<-ggplot(df1,aes(x=M.x,y=M.y))+geom_point()+geom_abline(aes(intercept=0,slope=1))
grid.arrange(g1,g2,g3,g4,g5,g6)
df1[df1$tmax.x>df1$tmax.y &!is.na(df1$tmax.x),]
df1[df1$M.x>df1$M.y*2 &!is.na(df1$M.x),]

#Numbers look okay, but some are different from previous paper

write.csv(table(BelizeFam$Gear,BelizeFam$Station),"temp.csv")
names(BelizeFam)[names(BelizeFam)=="FAMILY"]="Family"

x<-match(BelizeAll$Species,Belize$scinameFishbase)
BelizeAll$Trophic2<-Belize$TROPHIC.GROUP[x]
write.csv(BelizeAll[order(BelizeAll$Family,BelizeAll$Species),c("Family","Species","Common","tmax","Lmax","Lm","Linf","K","M","Trophic","Trophic2")],"BelizeAll.csv")

dim(BelizeFam)

#Make additional parameter files for LBB and LBSPR
BelizeAll$MK<-BelizeAll$M/BelizeAll$K
BelizeAll$Loo<-BelizeAll$Linf
summary(BelizeAll)

BelizeAllFL<-BelizeAll
a<-match(BelizeAll$Species,Belizepardf$Species)
summary(a-1:length(a))
BelizeAllFL$Linf<-Belizepardf$Loo
BelizeAllFL$Loo<-Belizepardf$Loo
BelizeAllFL$Lm<-Belizepardf$Lm
BelizeAllFL$K<-Belizepardf$K
BelizeAllFL$M<-Belizepardf$M
BelizeAllFL$MK<-Belizepardf$MK

## length frequency
sptoplotB<-c("Caranx hippos","Seriola dumerili","Gerres cinereus","Haemulon plumierii",
"Lutjanus synagris","Ocyurus chrysurus","Lutjanus analis","Lutjanus vivanus",
"Mycteroperca bonaci","Epinephelus itajara","Epinephelus guttatus","Sphyraena barracuda",
  "Lachnolaimus maximus","Lutjanus cyanopterus","Lutjanus jocu")

label2<-paste0("(",letters[1:15],") ",sptoplotB)
Belize$label2<-label2[match(Belize$scinameFishbase,sptoplotB)]

length(sptoplotB)
gs2<-ggplot(filter(Belize,scinameFishbase %in% sptoplotB & !Gear=="Unknown"),
  aes(x=TL2_CM,col=Gear,fill=Gear))+
  geom_histogram(position = "dodge2")+  
  facet_wrap(label2 ~.,scale="free",ncol=3)+
   theme_classic()+ theme(strip.background = element_blank(),
     strip.text=element_text(hjust=0,face = "italic"),
     legend.position = "bottom",
     legend.title = element_blank())+
  ylab("Proportion")+xlab("Length (cm)")+
  geom_vline(aes(xintercept=Lm),lty=2)
gs2
ggsave("Figure4.jpg",gs2,height=9,width=6.5)
length(sptoplotB)
#BelizeAll$n<-BelizeAll$n.x
spsup<-BelizeAll$Species[BelizeAll$n>=20 & ! BelizeAll$Species %in% sptoplotB]
length(spsup)
gs3<-ggplot(filter(Belize,scinameFishbase %in% spsup & !Gear=="Unknown"),
  aes(x=TL2_CM,col=Gear,fill=Gear))+
  geom_histogram(position = "dodge2")+  
  facet_wrap(scinameFishbase ~.,scale="free",ncol=4)+
   theme_classic()+ theme(strip.background = element_blank(),
     strip.text=element_text(hjust=0,face = "italic",size=7),
     legend.position = "bottom",
     legend.title = element_blank())+
  ylab("Proportion")+xlab("Length (cm)")+
  geom_vline(aes(xintercept=Lm),lty=2)
gs3
ggsave("BelizeHistSup.jpg",gs3,height=9,width=6.5)

# #The following don't work well. 
# pdf("BelizeHistograms.pdf",height=11,width=8.5)
# for(i in 1:4) {
# print(ggplot(filter(Belize,scinameFishbase %in% sptoplotB),
#   aes(x=TL2_CM,..density..,col=Gear,fill=Gear))+
#   geom_histogram(position = "dodge2")+
#   facet_wrap_paginate(scinameFishbase ~.,scale="free",nrow=4,ncol=3,page=i)+
#    theme_classic()+ theme(strip.background = element_blank(),
#      strip.text=element_text(hjust=0,face = "italic"),
#      legend.position = "bottom")+
#   ylab("Proportion")+xlab("Length (cm)")+
#   geom_vline(aes(xintercept=Lm),lty=2))
# }
# dev.off()
# 
# pdf("BelizeHistogramsCount.pdf",height=11,width=8.5)
# for(i in 1:4) {
# print(ggplot(filter(Belize,scinameFishbase %in% sptoplotB),
#   aes(x=TL2_CM,col=Gear,fill=Gear))+
#   geom_histogram(position = "dodge2")+
#   facet_wrap_paginate(scinameFishbase ~.,scale="free",nrow=4,ncol=3,page=i)+
#    theme_classic()+ theme(strip.background = element_blank(),
#      strip.text=element_text(hjust=0,face = "italic"),
#      legend.position = "bottom")+
#   ylab("Count")+xlab("Length (cm)")+
#   geom_vline(aes(xintercept=Lm),lty=2))
# }
# dev.off()

# At this point stop and run LBB

## Read in LBB results to set up files for LBSPR
theme_set(theme_classic())
lbbnoprior<-read.csv("NoneBelizeLBBoutTable.csv")
lbbFLprior<-read.csv("FLBelizeLBBoutTable.csv")
lbbdataprior<-read.csv("DataBelizeLBBoutTable.csv")
splbbno<-lbbnoprior$Species[!is.na(lbbnoprior$Linf.med) ]
splbbfl<-lbbFLprior$Species[!is.na(lbbFLprior$Linf.med) ]
splbbdat<-lbbFLprior$Species[!is.na(lbbdataprior$Linf.med) ]
length(splbbno)
length(splbbfl)
length(splbbdat)

spllb<-splbbno
noprior<-lbbnoprior %>% mutate(Linf=Linf.med,MK=MK.med,M=NA,K=NA,Lm=Lm50) %>%
  dplyr::select(Species,Linf,MK,M,K,Lm)
flprior<-lbbFLprior %>% mutate(Linf=Linf.med,MK=MK.med,M=NA,K=NA,Lm=Lm50) %>%
  dplyr::select(Species,Linf,MK,M,K,Lm)
datprior<-lbbdataprior %>% mutate(Linf=Linf.med,MK=MK.med,M=NA,K=NA,Lm=Lm50) %>%
  dplyr::select(Species,Linf,MK,M,K,Lm)

x<-BelizeGrowth %>% mutate(Linf=Loo,M=NA,Lm=NA) %>% dplyr::select(Species,Linf,K,M,Lm)
w<-BelizeMat %>% mutate(Linf=NA,K=NA,M=NA,MK=NA) %>% dplyr::select(Species,Linf,K,M,Lm)
y<-Belizepardf %>% mutate(Linf=Loomed,MK=MKmed,M=Mmed,K=Kmed) %>%
  dplyr::select(Species,M,K,MK,Linf,Lm)
z<-BelizeAll %>% mutate(MK=M/K) %>%
  dplyr::select(Species,M,K,MK,Linf,Lm)
  
allPars<-bind_rows(list("LBB default"=noprior,"LBB FL"=flprior,"LBB data"=datprior,
  Data=x,Fishlife=y,"Data mean"=z),.id="Source")
head(allPars)


allParsFilter<-filter(allPars,Species %in% splbbno & Species %in% splbbfl)
ggplot() +
  geom_point(data=filter(allParsFilter,Source=="Data"),
    aes(x=Species,y=Linf),color="black",alpha=0.3) + 
  geom_point(data=filter(allParsFilter,Source!="Data"),
    aes(x=Species,y=Linf,color=Source),size=3) + 
  coord_flip() +
  scale_color_brewer(palette="Set1")


ggplot(filter(allParsFilter,Source!="Data"),aes(x=Species,y=MK,color=Source)) +
  geom_point(alpha=0.8) + coord_flip() +
  geom_point(data=filter(allParsFilter,Source == "Final"),size=3) + 
  geom_hline(yintercept=1.5)+
   scale_color_brewer(palette="Set1")+ylab("M/K")

summary(BelizeAll$M/BelizeAll$K)

ggplot(BelizeAll,aes(x=M/K))+geom_histogram()
ggplot(BelizeAll,aes(x=log(M/K)))+geom_histogram()

# Make LBB data files for lbspr
BelizeAllLBBnoprior<-BelizeAll
a<-match(lbbnoprior$Species,BelizeAll$Species)
summary(a)
table(is.na(a))
BelizeAllLBBnoprior$Linf[]<-NA
BelizeAllLBBnoprior$MK[]<-NA
BelizeAllLBBnoprior$M[]<-NA
BelizeAllLBBnoprior$K[]<-NA
BelizeAllLBBnoprior$Linf[a]<-lbbnoprior$Linf.med
BelizeAllLBBnoprior$Loo<-BelizeAllLBBnoprior$Linf
BelizeAllLBBnoprior$MK[a]<-lbbnoprior$MK.med
BelizeAllLBBnoprior$K[a]<-BelizeAll$K[a]
BelizeAllLBBnoprior$M<-BelizeAllLBBnoprior$MK*BelizeAllLBBnoprior$K
summary(BelizeAllLBBnoprior[,c("M","K","MK")])

BelizeAllLBBFLprior<-BelizeAll
a<-match(lbbFLprior$Species,BelizeAll$Species)
summary(a)
table(is.na(a))
BelizeAllLBBFLprior$Linf[]<-NA
BelizeAllLBBFLprior$MK[]<-NA
BelizeAllLBBFLprior$M[]<-NA
BelizeAllLBBFLprior$K[]<-NA
BelizeAllLBBFLprior$Linf[a]<-lbbFLprior$Linf.med
BelizeAllLBBFLprior$Loo<-BelizeAllLBBFLprior$Linf
BelizeAllLBBFLprior$MK[a]<-lbbFLprior$MK.med
BelizeAllLBBFLprior$K[a]<-lbbFLprior$MK.med

summary(BelizeAllLBBFLprior)
x<-BelizeAll[!is.na(BelizeAllLBBnoprior$MK),]
x$K2<-BelizeAllLBBnoprior$MK[!is.na(BelizeAllLBBnoprior$MK)]

ggplot(x,aes(x=M,y=K))+
  geom_point()+
  geom_abline(intercept=0,slope=1/1.5)+
  geom_point(aes(x=K/MK),color="red")

#Use data K, and M = K/MK    

df5<-data.frame(Data=BelizeAll$M,FL=BelizeAllFL$M)
ggplot(df5,aes(x=Data,y=FL))+geom_point()+geom_abline()
