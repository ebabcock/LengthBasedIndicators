#This code reads in the Guatemala data, does data cleaning, gets life history parameters
# and makes summary figures and tables

## Code to get rid of repeated warning
#!diagnostics off
library(tidyverse)
library(gridExtra)
library( FishLife )
library(rfishbase)
library(ggplot2)
#library(TMB)
setwd("C:/Users/ebabcock/Dropbox/BethGlovers/2020 finfish")
source("babcockfunctionslength2020.r")

#Read in Guatemala data
#old<-read.csv("20210119_DB_Pacific_Guatemala_WCS_2019_2020.csv")
Guatemala<-read.csv("20210611_DB_Pacific_Guatemala_WCS_2019_2020.csv")
#dim(old)
dim(Guatemala)
#table(old$Arte.de.pesca..fishing.gear.)
table(Guatemala$Arte.de.pesca..fishing.gear.)
  
# Fix formatting of species names and other columns
sort(table(Guatemala$Scientific.name))
sort(table(Guatemala$Common.name))
Guatemala$Scientific.name<-paste0(toupper(substring(Guatemala$Scientific.name, 1, 1)), 
                               tolower(substring(Guatemala$Scientific.name, 2)))
Guatemala$Scientific.name<-trimws(Guatemala$Scientific.name)
Guatemala$Common.name<-trimws(Guatemala$Common.name)

Guatemala$Arte.de.pesca..fishing.gear.<-trimws(Guatemala$Arte.de.pesca..fishing.gear.)
Guatemala$Arte.de.pesca..fishing.gear.<-paste0(toupper(substring(Guatemala$Arte.de.pesca..fishing.gear., 1, 1)), 
                               tolower(substring(Guatemala$Arte.de.pesca..fishing.gear., 2)))
table(Guatemala$Arte.de.pesca..fishing.gear.)
Guatemala$Gear<-Guatemala$Arte.de.pesca..fishing.gear.
Guatemala$Gear[Guatemala$Gear=="Coastal/demersal longlines"]<-"Coastal LL"
Guatemala$Gear[Guatemala$Gear=="Long line"]<-"Pelagic LL"
table(Guatemala$Gear)

format(as.Date(Guatemala$Fecha..date.[1:100],format="%m/%d/%Y"),"%m")
Guatemala$Date<-as.Date(Guatemala$Fecha..date.,format="%m/%d/%Y")
summary(Guatemala$Date)
Guatemala$Month<-format(Guatemala$Date,"%m")
Guatemala$Year<-format(Guatemala$Date,"%Y")
table(Guatemala$Month,Guatemala$Year)

Guatemala$Lugar..place.<-trimws(Guatemala$Lugar..place.)
Guatemala$Lugar..place.<-paste0(toupper(substring(Guatemala$Lugar..place., 1, 1)), 
                               tolower(substring(Guatemala$Lugar..place., 2)))
table(Guatemala$Lugar..place.)
Guatemala$Station<-Guatemala$Lugar..place.
Guatemala$Station[grep("san jose",Guatemala$Station)]<-"Puerto de San Jose"
Guatemala$Station[grep("buena vista",Guatemala$Station)]<-"Buena Vista"
Guatemala$Station[grep("Buena vista",Guatemala$Station)]<-"Buena Vista"
table(Guatemala$Station)

Guatemala$Habitat<-ifelse(Guatemala$Gear=="Pelagic LL","Offshore","Coastal")

Guatemala$Scientific.name[Guatemala$Scientific.name=="Zapterix exasperata"]="Zapteryx exasperata"
Guatemala$Scientific.name[Guatemala$Scientific.name=="Pseudobatos leucorbynchus"]="Pseudobatos leucorhynchus"
Guatemala$Scientific.name[Guatemala$Scientific.name=="Centrophomus viridis"]="Centropomus viridis"
Guatemala$Scientific.name[Guatemala$Scientific.name=="Centrophomus robalito"]="Centropomus robalito"
Guatemala$Scientific.name[Guatemala$Scientific.name=="Polidactylus opercularis"]="Polydactylus opercularis"
Guatemala$Scientific.name[Guatemala$Scientific.name=="Katsuwonus pelamises"]="Katsuwonus pelamis"

# Make new column scinamefishbase with fb spellings, only to species, marine fish
Guatemala$scinameFishbase<-Guatemala$Scientific.name
dim(Guatemala) #17628
Guatemala<-filter(Guatemala,!scinameFishbase %in% 
    c("Penaeus vannamei","Panulirus gracilis","Camaron cebra","Peneus californiensis",
      "Callinectes arcuatus","Callinectes toxotes","Camaron cebra" ,"Xiphopenaeus kroyeri",
      "Panulirus inflatus","Portunus asper","Evibacus princeps"))  #invertebrates
dim(Guatemala) #17080  #This eliminates the trawls.

GuatemalaSp<-read.csv("GuatemlaSpeciesData.csv")
GuatemalaSp$Scientific.name<-trimws(GuatemalaSp$Scientific.name)
x<-match(Guatemala$scinameFishbase,GuatemalaSp$Scientific.name)
summary(x)
sort(table(Guatemala$scinameFishbase[is.na(x)]))

Guatemala$Family<-trimws(GuatemalaSp$Family[x])
Guatemala$Trophic2<-GuatemalaSp$Trophic.group[x]
GuatemalaFam<-Guatemala

#table(c(Guatemala$Family,Belize$FAMILY))

#Remove not known to species
table(Guatemala$scinameFishbase[grep("spp",Guatemala$scinameFishbase)])
Guatemala$scinameFishbase[grep("spp",Guatemala$scinameFishbase)]=NA  #Remove fish id'ed to genus
Guatemala<-filter(Guatemala,!is.na(Guatemala$scinameFishbase))
dim(Guatemala) #15938

#Different spelling in fishbase
Guatemala$scinameFishbase[Guatemala$scinameFishbase=="Narcine vermiculatus"]="Narcine vermiculata"
Guatemala$scinameFishbase[Guatemala$scinameFishbase=="Isopisthus ramifer"]="Isopisthus remifer"
Guatemala$scinameFishbase[Guatemala$scinameFishbase=="Lobotes pacificus"]="Lobotes pacifica"
# A. laticeps is recent split from A.narinari, so use narinari life history
Guatemala$scinameFishbase[Guatemala$scinameFishbase=="Aetobatus laticeps"]="Aetobatus narinari"

# Check length types
Guatemala$Ancho.de.disco...Disk.width..cm.<-as.numeric(Guatemala$Ancho.de.disco...Disk.width..cm.)
table(Guatemala$scinameFishbase, is.na(Guatemala$Ancho.de.disco...Disk.width..cm.))
table(is.na(Guatemala$Longitud.total....cm...Total.lenght.), is.na(Guatemala$Ancho.de.disco...Disk.width..cm.))
#We have total length or disc width for nearly all records
table(Guatemala$scinameFishbase[is.na(Guatemala$Longitud.total....cm...Total.lenght.)& is.na(Guatemala$Ancho.de.disco...Disk.width..cm.)])
table(Guatemala$scinameFishbase[is.na(Guatemala$Longitud.total....cm...Total.lenght.)& is.na(Guatemala$Ancho.de.disco...Disk.width..cm.) &
    !is.na(Guatemala$Lorgitud.de.orquilla.L.cm...furcal.lenght.)])
table(Guatemala$scinameFishbase[is.na(Guatemala$Longitud.total....cm...Total.lenght.)& is.na(Guatemala$Ancho.de.disco...Disk.width..cm.) &
    is.na(Guatemala$Lorgitud.de.orquilla.L.cm...furcal.lenght.)])
table(Guatemala$scinameFishbase[is.na(Guatemala$Longitud.total....cm...Total.lenght.)& is.na(Guatemala$Ancho.de.disco...Disk.width..cm.) &
    !is.na(Guatemala$Interdorsal)])
table(is.na(Guatemala$Lorgitud.de.orquilla.L.cm...furcal.lenght.),is.na(Guatemala$Longitud.total....cm...Total.lenght.))
#Pelagic thresher and some other sharks only have interdorsal. 
Guatemala$Ltype<-rep("",dim(Guatemala)[1])
Guatemala$Ltype[!is.na(Guatemala$Longitud.total....cm...Total.lenght.)]<-"TL"
Guatemala$Ltype[Guatemala$Ltype=="" & 
    !is.na(Guatemala$Ancho.de.disco...Disk.width..cm.)]<-"WD"
Guatemala$Ltype[Guatemala$Ltype=="" & 
    !is.na(Guatemala$Lorgitud.de.orquilla.L.cm...furcal.lenght.)]<-"FL"
Guatemala$Ltype[Guatemala$Ltype=="" & 
    !is.na(Guatemala$Interdorsal)]<-"ID"
table(Guatemala$Ltype)

Guatemala[Guatemala$Ltype=="",]
write.csv((table(Guatemala$scinameFishbase,Guatemala$Ltype)),"temp.csv")
Guatemala<-Guatemala %>% 
  mutate(Length=ifelse(Ltype=="TL",Longitud.total....cm...Total.lenght.,NA))%>%
  mutate(Length=ifelse(Ltype=="WD",Ancho.de.disco...Disk.width..cm.,Length))
summary(Guatemala$Length)
Guatemala<-filter(Guatemala,!is.na(Length))
dim(Guatemala)
write.csv((table(Guatemala$scinameFishbase,Guatemala$Ltype)),"temp2.csv")

#Identify most common species
x<-sort(table(Guatemala$scinameFishbase))
x<-x[x>=30]
x
length(x)
Guatemalasp<-data.frame(Species=names(x))
x<-match(Guatemalasp$Species,Guatemala$scinameFishbase)
summary(x)
Guatemalasp$Common<-Guatemala$Common.name[x]
dim(Guatemalasp)

# Get fishbase data for all
x<-sort(unique(Guatemala$scinameFishbase))
y<-validate_names(x)
length(x)
length(y)
y[!y %in% x]
x[!(x%in%y)]

GuatemalaSpecies<-species(x)
GuatemalaEcology<-ecology(x) 
GuatemalaAll<-dplyr::select(GuatemalaSpecies,Species,FBname)
dim(GuatemalaAll)
a<-match(GuatemalaAll$Species,Guatemala$scinameFishbase)
summary(a)
GuatemalaAll$Ltype<-Guatemala$Ltype[a]
table(GuatemalaAll$Ltype)
GuatemalaAll$Species[GuatemalaAll$Ltype=="WD"]
GuatemalaSource<-expand.grid(Species=GuatemalaAll$Species,Lm=NA,Linf=NA,K=NA,tmax=NA,Lmax=NA)
table(GuatemalaSpecies$Saltwater)
list_fields("K")
list_fields("FoodTroph")
list_fields("Lm")
list_fields("tmax")
list_fields("Lmax")
GuatemalaGrowth<-popgrowth(x)
GuatemalaPopChar<-popchar(x) %>% mutate(Lmax=as.numeric(Lmax))
GuatemalaMat<-maturity(x)
#See if Lm is avilalbe or LmatMin by species
y<-GuatemalaMat %>% group_by(Species) %>%
  summarize(Lm=median(Lm,na.rm=TRUE),LengthMatMin=median(LengthMatMin,na.rm=TRUE))
z<-y$Species[is.na(y$Lm &!is.na(y$LengthMatMin))]
z
GuatemalaSource$Lm[!is.na(y$Lm)]<-"FB.Lm"
GuatemalaSource$Lm[is.na(y$Lm &!is.na(y$LengthMatMin))]<-"FB.LMatMin"
GuatemalaMat$Lm[GuatemalaMat$Species %in%z]<-GuatemalaMat$LengthMatMin[GuatemalaMat$Species %in%z]

#Length conversions
Guatemalapopll<-popll(x) 
GuatemalapopllFL<-Guatemalapopll %>% filter(Length1=="TL" & Length2=="FL") %>%
  filter(!duplicated(Species)) %>% mutate(a=as.numeric(a),b=as.numeric(b))
GuatemalapopllSL<-Guatemalapopll %>% filter(Length1=="TL" & Length2=="SL") %>%
  filter(!duplicated(Species))  %>% mutate(a=as.numeric(a),b=as.numeric(b))
GuatemalapopllWD<-Guatemalapopll %>% filter(Length1=="WD") %>%
  filter(!duplicated(Species)) 
dim(GuatemalapopllWD) #No conversions from WD
summary(GuatemalapopllFL)
a<-match(GuatemalaAll$Species,GuatemalapopllFL$Species)

x<-merge(GuatemalaAll,GuatemalapopllFL[,c("Species","a","b")],all.x=TRUE,by="Species")
head(x)
summary(x)
GuatemalaAll<-merge(GuatemalaAll,GuatemalapopllFL[,c("Species","a","b"),],all.x=TRUE,by="Species")

GuatemalapopllSL$SL.a<-GuatemalapopllSL$a
GuatemalapopllSL$SL.b<-GuatemalapopllSL$b
GuatemalaAll<-merge(GuatemalaAll,GuatemalapopllSL[,c("Species","SL.a","SL.b"),],all.x=TRUE,by="Species")
#a and b in GuatemalaAll are now FL to TL conversions, SL.a and SL.b convert from SL

# Now get Lm in TL or WD
table(GuatemalaMat$Type1)
x<-match(GuatemalaMat$Species,GuatemalaAll$Species)
summary(x)
y<-GuatemalaAll$a[x]+GuatemalaAll$b[x]* GuatemalaMat$Lm
y[!GuatemalaMat$Type1=="FL"]<-NA
z<-GuatemalaAll$SL.a[x]+GuatemalaAll$SL.b[x]* GuatemalaMat$Lm
z[!GuatemalaMat$Type1=="SL"]<-NA
GuatemalaMat$TL.Lm<-ifelse(GuatemalaMat$Type1=="FL",y,GuatemalaMat$Lm)
GuatemalaMat$TL.Lm<-ifelse(GuatemalaMat$Type1=="SL",z,GuatemalaMat$TL.Lm)
GuatemalaMat$TL.Lm<-ifelse(GuatemalaMat$Type1=="WD",GuatemalaMat$Lm,GuatemalaMat$TL.Lm)
table(is.na(GuatemalaMat$Lm),is.na(GuatemalaMat$Type1))
table(Guatemala$Family[!is.na(GuatemalaMat$Lm)&is.na(GuatemalaMat$Type1)])
GuatemalaMat$TL.Lm<-ifelse(is.na(GuatemalaMat$Type1),GuatemalaMat$Lm,GuatemalaMat$TL.Lm)
GuatemalaMat$Lm<-GuatemalaMat$TL.Lm
GuatemalaMat<-dplyr::select(GuatemalaMat,Species,Type1,Lm)

#Loo
x<-match(GuatemalaGrowth$Species,GuatemalaAll$Species)
y<-GuatemalaAll$a[x]+GuatemalaAll$b[x]* GuatemalaGrowth$Loo
y[!GuatemalaGrowth$Type=="FL"]<-NA
z<-GuatemalaAll$SL.a[x]+GuatemalaAll$SL.b[x]* GuatemalaGrowth$Loo
GuatemalaGrowth$TL.Loo<-ifelse(GuatemalaGrowth$Type=="FL",y,GuatemalaGrowth$Loo)
z[!GuatemalaGrowth$Type=="SL"]<-NA
GuatemalaGrowth$TL.Loo<-ifelse(GuatemalaGrowth$Type=="SL",z,GuatemalaGrowth$TL.Loo)
GuatemalaGrowth$TL.Loo<-ifelse(GuatemalaGrowth$Type=="WD",GuatemalaGrowth$Loo,GuatemalaGrowth$TL.Loo)
GuatemalaGrowth$TL.Loo<-ifelse(is.na(GuatemalaGrowth$Type),GuatemalaGrowth$Loo,GuatemalaGrowth$TL.Loo)
GuatemalaGrowth<-GuatemalaGrowth %>%  rename(Linf=TL.Loo)

#Lmax convert to TL in popchar and Species or use WD
x<-match(GuatemalaPopChar$Species,GuatemalaAll$Species)
summary(x)
y<-GuatemalaAll$a[x]+GuatemalaAll$b[x]* GuatemalaPopChar$Lmax
y[!GuatemalaPopChar$Type=="FL"]<-NA
z<-GuatemalaAll$SL.a[x]+GuatemalaAll$SL.b[x]* GuatemalaPopChar$Lmax
z[!GuatemalaPopChar$Type=="SL"]<-NA
GuatemalaPopChar<- GuatemalaPopChar %>%
  mutate(Lmax.TL=ifelse(Type=="FL",y,Lmax)) %>%
  mutate(Lmax.TL=ifelse(Type=="SL",z,Lmax.TL)) %>%
  mutate(Lmax.TL=ifelse(Type=="WD",Lmax,Lmax.TL)) %>%
  mutate(Lmax=Lmax.TL)
summary(GuatemalaPopChar$Lmax)

# New recent data that is not in fishbase yet
#Oshitani Linf=2.08+1.32*216.4
#Mair. Averaged across methods
# Linf<-mean(c(59.2,62.0, 59.0))
# Lm<-mean(c(33.0, 32.6))
# K<-mean(c(0.32, 0.28, 0.33))
# round(c(Linf,K,Lm),2)

newdat<-read.csv("newdatGuatemala.csv")
names(newdat)[1]<-"Species"
GuatemalaGrowth<-bind_rows(list(FB=GuatemalaGrowth,New=dplyr::select(newdat,Species,K,Linf,tmax)),.id="id")
dim(GuatemalaGrowth)
GuatemalaPopChar<-bind_rows(list(FB=dplyr::select(GuatemalaPopChar,Species,Lmax,tmax),New=dplyr::select(newdat,Species,Lmax,tmax)),.id="id")
dim(GuatemalaPopChar)
GuatemalaMat<-bind_rows(list(FB=dplyr::select(GuatemalaMat,Species,Lm),New=dplyr::select(newdat,Species,Lm)),.id="id")
dim(GuatemalaMat)

#Calculate average Lm, K, Linf, K and max tmax and Lmax
GuatemalaGrowthSum<-GuatemalaGrowth %>% group_by(Species) %>%
  summarize(Linf=median(Linf,na.rm=TRUE),
            K=median(K,na.rm=TRUE),tmax=max(tmax,na.rm=TRUE)) %>%
  mutate(tmax=ifelse(tmax<0,NA,tmax))
dim(GuatemalaGrowthSum)
summary(GuatemalaGrowthSum)

GuatemalaPopCharSum<-GuatemalaPopChar %>% group_by(Species) %>%
  summarize(Lmax=max(Lmax,na.rm=TRUE),tmax=max(tmax,na.rm=TRUE)) %>%
  mutate(Lmax=ifelse(Lmax<0,NA,Lmax),tmax=ifelse(tmax<0,NA,tmax))
summary(GuatemalaPopCharSum)

GuatemalaMatSum<-GuatemalaMat %>% group_by(Species) %>%
  summarize(Lm=median(Lm,na.rm=TRUE)) %>%
  mutate(Lm=ifelse(Lm<0,NA,Lm))
summary(GuatemalaMatSum)
dim(GuatemalaMatSum)

#Get max Lmax
x<-GuatemalaAll$a+GuatemalaAll$b*GuatemalaSpecies$Length
x[!GuatemalaSpecies$LTypeMaxM=="FL"]<-NA
GuatemalaSpecies$Length.TL<-ifelse(GuatemalaSpecies$LTypeMaxM=="FL",x,GuatemalaSpecies$Length)
table(is.na(GuatemalaSpecies$Length.TL),is.na(GuatemalaPopCharSum$Lmax))
GuatemalaPopCharSum<-GuatemalaPopCharSum %>% mutate(Lmax=ifelse(is.na(Lmax),
  GuatemalaSpecies$Length.TL,Lmax)) %>%
  mutate(Lmax=ifelse(!is.na(Lmax) & Lmax<GuatemalaSpecies$Length.TL,
    GuatemalaSpecies$Length.TL,Lmax))
summary(GuatemalaPopCharSum[,c("Lmax","tmax")])

# Get max tmax
GuatemalaGrowthSum$tmax<-ifelse(is.na(GuatemalaGrowthSum$tmax),GuatemalaPopCharSum$tmax,
  GuatemalaGrowthSum$tmax)
GuatemalaGrowthSum$tmax<-ifelse(!is.na(GuatemalaGrowthSum$tmax+GuatemalaPopCharSum$tmax) &
    GuatemalaPopCharSum$tmax>GuatemalaGrowthSum$tmax,GuatemalaPopCharSum$tmax,GuatemalaGrowthSum$tmax)
summary(GuatemalaGrowthSum$tmax)
GuatemalaGrowthSum<-merge(GuatemalaGrowthSum,GuatemalaMatSum,by="Species")
GuatemalaGrowthSum$Lmax<-GuatemalaPopCharSum$Lmax
summary(GuatemalaGrowthSum)
GuatemalaSource$Linf[!is.na(GuatemalaGrowthSum$Linf)]<-"FB"
GuatemalaSource$K[!is.na(GuatemalaGrowthSum$K)]<-"FB"
GuatemalaSource$Lmax[!is.na(GuatemalaGrowthSum$Lmax)]<-"FB"
GuatemalaSource$tmax[!is.na(GuatemalaGrowthSum$tmax)]<-"FB"
GuatemalaSource$Lm[!is.na(GuatemalaGrowthSum$Lm) & is.na(GuatemalaSource$Lm)]<-"FB"

# Combine tables
GuatemalaAll<-merge(GuatemalaAll,GuatemalaGrowthSum) 

#Check if Lmax in fishbase is larger than max measured length
y<-Guatemala %>% group_by(scinameFishbase) %>%
  summarize(LmaxDat=max(Length,na.rm=TRUE)) %>%
  rename(Species=scinameFishbase)
summary(y)
dim(y)
GuatemalaAll<-merge(GuatemalaAll,y,by="Species")
table(GuatemalaAll$Lmax>GuatemalaAll$LmaxDat)
x<-GuatemalaAll$Species[GuatemalaAll$LmaxDat > GuatemalaAll$Lmax *1.1 ]
df1<-filter(Guatemala, scinameFishbase %in% x)
df1$Lmax<-GuatemalaAll$Lmax[match(df1$scinameFishbase,GuatemalaAll$Species)]
ggplot(df1,aes(x=Length))+
  geom_histogram(fill="blue")+
  facet_wrap(scinameFishbase~.,scales="free",ncol=2) +
  geom_vline(aes(xintercept=Lmax),col="red")
GuatemalaSource$Lmax[GuatemalaAll$LmaxDat>GuatemalaAll$Lmax | is.na(GuatemalaAll$Lmax )]<-"Data"
GuatemalaAll$Lmax<-ifelse(GuatemalaAll$LmaxDat>GuatemalaAll$Lmax| is.na(GuatemalaAll$Lmax ),GuatemalaAll$LmaxDat,GuatemalaAll$Lmax)

GuatemalaAll$Lmax[GuatemalaAll$Species=="Aetobatus narinari"]<-330 #From Fishbase

#Trophic
GuatemalaAll$Trophic<-ifelse(is.na(GuatemalaEcology$DietTroph),GuatemalaEcology$FoodTroph,GuatemalaEcology$DietTroph)
GuatemalaAll$Trophic[GuatemalaAll$Species=="Ophichthus gomesii"] <- 4 #Based on similar species, from fb website
summary(GuatemalaAll[,c("Lm","Linf","Lmax","K","tmax","Trophic")])
GuatemalaAll[is.na(GuatemalaAll$Trophic),]
GuatemalaSource$Trophic<-ifelse(is.na(GuatemalaAll$Trophic),"fb.genus","fb")
a<-which(is.na(GuatemalaAll$Trophic))
for(i in a) {
  x<-species_list(Genus=GuatemalaAll$Genus[i])
  y<-ecology(x)
  y$Trophic<-ifelse(is.na(y$DietTroph),y$FoodTroph,y$DietTroph)
  GuatemalaAll$Trophic[i]<-median(y$Trophic,na.rm=TRUE)
}
summary(GuatemalaAll$Trophic)

# Common species without data
GuatemalaAll$Species[is.na(GuatemalaAll$tmax) & GuatemalaAll$Species %in% Guatemalasp$Species]
GuatemalaAll$FBname[is.na(GuatemalaAll$tmax) & GuatemalaAll$Species %in% Guatemalasp$Species]
GuatemalaAll$Species[is.na(GuatemalaAll$Linf) & GuatemalaAll$Species %in% Guatemalasp$Species]
GuatemalaAll$Species[is.na(GuatemalaAll$Lm) & GuatemalaAll$Species %in% Guatemalasp$Species]

#Sample size
x<-Guatemala %>% group_by(scinameFishbase) %>% 
  summarize(n=length(scinameFishbase)) %>%
  rename(Species=scinameFishbase)
GuatemalaAll<-merge(GuatemalaAll,x,by="Species")
names(GuatemalaAll)

#Compare to Alex's data
x<-match(GuatemalaAll$Species,Guatemala$scinameFishbase)
summary(x)
GuatemalaAll$SpeciesOld<-Guatemala$Scientific.name[x]
x<-match(GuatemalaAll$SpeciesOld,GuatemalaSp$Scientific.name)
summary(x)
df<-bind_cols(GuatemalaAll,GuatemalaSp[x,])
ggplot(df,aes(x=Lm...11,y=Lm...23))+geom_point()+geom_abline(aes(intercept=0,slope=1))

#Add genus and species to All, and change spelling to match FishBase
for(i in 1:length(GuatemalaAll$Species)) {
  GuatemalaAll$Genus[i]=unlist(strsplit(GuatemalaAll$Species[i]," "))[1]
  GuatemalaAll$species[i]=unlist(strsplit(GuatemalaAll$Species[i]," "))[2]
}

GuatemalaAll$Genus[GuatemalaAll$Species=="Hypanus longus"]<-"Dasyatis"
GuatemalaAll$species[GuatemalaAll$Species=="Hypanus longus"]<-"longa"
GuatemalaAll$Genus[GuatemalaAll$Species=="Pseudobatos leucorhynchus"]<-"Rhinobatos"
GuatemalaAll$Genus[GuatemalaAll$Species=="Styracura pacifica"]<-"Himantura"
GuatemalaAll$species[GuatemalaAll$Species=="Styracura pacifica"]<-"pacifica"
GuatemalaAll$species[GuatemalaAll$Species=="Lobotes pacifica"]<-"pacificus"
GuatemalaAll$species[GuatemalaAll$Species=="Narcine vermiculata"]<-"vermiculatus"


#Select only needed columns in All
names(GuatemalaAll)
GuatemalaAll$Common<-GuatemalaAll$FBname
summary(GuatemalaAll)

#Run fishlife for all
Guatemalapars<-list()
GuatemalaAll$species[GuatemalaAll$Species=="Diplodus argenteus"]="argenteus argenteus"
GuatemalaAll$predictive<-rep(1,dim(GuatemalaAll)[1])
for(i in 1:dim(GuatemalaAll)[1]) {
 x <-try(Plot_taxa( Search_species(Genus=GuatemalaAll$Genus[i],Species=GuatemalaAll$species[i])$match_taxonomy ))
 Guatemalapars[[i]]<-x[[1]]
 if(length(x)==5) GuatemalaAll$predictive[i]=0
}
summary(GuatemalaAll)
summary(Guatemalapars)
length(Guatemalapars)
dim(GuatemalaAll)

##
params<-c("Lm","M","K","Loo")
Guatemalaparmean<-NULL
Guatemalaparmedian<-NULL
Guatemalaparse<-NULL
for(i in 1:dim(GuatemalaAll)[1]) {
  meanvalMK<-Guatemalapars[[i]]$Mean_pred["M"]-Guatemalapars[[i]]$Mean_pred["K"]
  names(meanvalMK)<-"MK"
  varvalMK<-Guatemalapars[[i]]$Cov_pred["M","M"]+Guatemalapars[[i]]$Cov_pred["K","K"]+2*Guatemalapars[[i]]$Cov_pred["M","K"]
  names(varvalMK)<-"MK"
  meanval<-lnorm.mean(c(Guatemalapars[[i]]$Mean_pred,meanvalMK),sqrt(c(diag(Guatemalapars[[i]]$Cov_pred),varvalMK)))
  seval<-lnorm.se(c(Guatemalapars[[i]]$Mean_pred,meanvalMK),sqrt(c(diag(Guatemalapars[[i]]$Cov_pred),varvalMK)))
  Guatemalaparmean<-rbind(Guatemalaparmean,meanval)
  Guatemalaparse<-rbind(Guatemalaparse,seval)
  Guatemalaparmedian<-rbind(Guatemalaparmedian,exp(c(Guatemalapars[[i]]$Mean_pred,meanvalMK)))
}
params<-c("Lm","M","K","Loo","MK")
Guatemalapardf<-data.frame(cbind(Guatemalaparmean[,params],Guatemalaparse[,params]),Guatemalaparmedian[,params])
names(Guatemalapardf)
names(Guatemalapardf)[6:10]<-paste0(params,"se")
names(Guatemalapardf)[11:15]<-paste0(params,"med")

dim(Guatemalapardf)

Guatemalapardf$Species<-GuatemalaAll$Species
df2<-merge(GuatemalaAll,Guatemalapardf,by="Species")
df2$Linf.Lmax<-ifelse(is.na(df2$Linf),Linf.Lmax(GuatemalaAll$Lmax),df2$Linf)
df2$Lm.Linf<-ifelse(is.na(df2$Lm.x),Lm.Linf(df2$Linf.Lmax),df2$Lm.x)
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
grid.arrange(g1,g2,g3,g4,g5,g6,ncol=2)

summary(GuatemalaAll)
GuatemalaAll$Common[is.na(GuatemalaAll$Lm) & GuatemalaAll$n>30]
GuatemalaAll$Common[is.na(GuatemalaAll$tmax)  & GuatemalaAll$n>30]
GuatemalaAll$Common[is.na(GuatemalaAll$K) & GuatemalaAll$n>30]

x<-match(GuatemalaAll$SpeciesOld,GuatemalaSp$Scientific.name)
summary(x)
GuatemalaAll$Common<-GuatemalaSp$Common..English.[x]
GuatemalaAll$Family<-GuatemalaSp$Family[x]

# add all parameters to GuatemalaAll
GuatemalaSource$Linf[is.na(GuatemalaAll$Linf)]<-"Lmax"
GuatemalaAll$Linf<-ifelse(is.na(GuatemalaAll$Linf) ,Linf.Lmax(GuatemalaAll$Lmax),GuatemalaAll$Linf)
GuatemalaSource$Lm[is.na(GuatemalaAll$Lm)]<-"Linf"
GuatemalaAll$Lm<-ifelse(is.na(GuatemalaAll$Lm),Lm.Linf(GuatemalaAll$Linf),GuatemalaAll$Lm)
GuatemalaSource$M<-ifelse(is.na(GuatemalaAll$tmax),NA,"Then")
GuatemalaSource$M<-ifelse(is.na(GuatemalaAll$tmax) &!is.na(GuatemalaAll$K),"Then.K",GuatemalaSource$M)
x<-ifelse(!is.na(M.Then(GuatemalaAll$tmax)),M.Then(GuatemalaAll$tmax),M.Then.K(GuatemalaAll$Linf,GuatemalaAll$K))
summary(x)
GuatemalaSource$M[is.na(x)]<-"FishLife"
GuatemalaAll$M<-ifelse(is.na(x), Guatemalaparmedian[,"M"],x)
GuatemalaSource$K[is.na(GuatemalaAll$K)]<-"Fishlife"
GuatemalaAll$K<-ifelse(is.na(GuatemalaAll$K),Guatemalaparmedian[,"K"],GuatemalaAll$K)
summary(GuatemalaAll[,c("Lm","M","K","Linf")])

GuatemalaAll[is.na(GuatemalaAll$M),]

Guatemalapardf<-Guatemalapardf %>% mutate(LooCV=Loose/Loo,KCV=Kse/K,MCV=Mse/M,LmCV=Lmse/Lm)
summary(Guatemalapardf)

# Print source table
GuatemalaSource$Predictive<-GuatemalaAll$predictive
GuatemalaSource$n<-GuatemalaAll$n
GuatemalaSource$Common<-GuatemalaAll$Common
x<-match(GuatemalaAll$Species,Guatemala$scinameFishbase)
GuatemalaSource$Family<-GuatemalaAll$Family
write.csv(GuatemalaSource[order(GuatemalaSource$Family),],"GuatemalaSource2.csv")
  
#Add fishbase numbers to Guatemala
x<-match(Guatemala$scinameFishbase,GuatemalaAll$Species)
Guatemala$Lmax<-GuatemalaAll$Lmax[x]
Guatemala$Lm<-GuatemalaAll$Lm[x]
Guatemala$Trophic<-GuatemalaAll$Trophic[x]
z<-load_taxa()
GuatemalaAll$Family<-z$Family[match(GuatemalaAll$Species,z$Species)]
table(GuatemalaAll$Family)
Guatemala$Family<-GuatemalaAll$Family[x]
table(Guatemala$Family)

# Check life history values
g1<-ggplot(GuatemalaAll,aes(x=Linf/Lmax))+geom_histogram()
g2<-ggplot(GuatemalaAll,aes(x=Lm/Linf))+geom_histogram()
g3<-ggplot(GuatemalaAll,aes(x=K/M))+geom_histogram()
g4<-ggplot(GuatemalaAll,aes(x=Lmax,y=M))+geom_point()
grid.arrange(g1,g2,g3,g4)

x<-filter(GuatemalaAll,M/K>4 | M/K <0.5)
x$M.K<-x$M/x$K
dim(x)
x

write.csv(table(GuatemalaFam$Gear,GuatemalaFam$Station),"temp.csv")

x<-table(Guatemala$scinameFishbase,Guatemala$Gear)
x

x<-match(Guatemala$scinameFishbase,GuatemalaAll$Species)
Guatemala$Common<-GuatemalaAll$Common[x]
df1<-filter(Guatemala,scinameFishbase %in% GuatemalaAll$Species[GuatemalaAll$n>100])

x<-match(GuatemalaAll$SpeciesOld,GuatemalaSp$Scientific.name)
summary(x)
GuatemalaAll$Trophic2<-GuatemalaSp$Trophic.group[x]
GuatemalaAll$Family<-trimws(GuatemalaSp$Family[x])
write.csv(GuatemalaAll[order(GuatemalaAll$Family,GuatemalaAll$Species),c("Family","Species","Common","tmax","Lmax","Lm","Linf","K","M","Trophic","Trophic2")],"GuatemalaAll2.csv")

dim(GuatemalaFam)
dim(Guatemala)

table(Guatemala$Trophic2)
Guatemala$Trophic2[Guatemala$Trophic2=="Zooplanktivore"]<-"Invertivore"

library(ggforce)

Guatemala$sciname2<-Guatemala$scinameFishbase
Guatemala$sciname2[Guatemala$sciname2=="Aetobatus narinari"]<-"Aetobatus laticeps"
sptoplotG<-sort(c("Hypanus longus","Sphyrna lewini","Cynoscion reticulatus","Peprilus snyderi",
  "Caranx caballus","Carcharhinus falciformis","Coryphaena hippurus","Diapterus peruvianus",
  "Carcharhinus limbatus","Sphyraena ensis","Mobula thurstoni","Aetobatus laticeps",
  "Mustelus lunulatus","Lutjanus guttatus","Mobula munkiana"))
label2<-paste0("(",letters[1:15],") ",sptoplotG)
Guatemala$label2<-label2[match(Guatemala$sciname2,sptoplotG)]

gs1<-ggplot(filter(Guatemala,sciname2 %in% sptoplotG),
#  aes(x=Length,..density..,col=gear,fill=gear))+
   aes(x=Length,col=Gear,fill=Gear))+
  geom_histogram(position = "dodge2")+  
  facet_wrap(label2 ~.,scale="free",ncol=3)+
   theme_classic()+ theme(strip.background = element_blank(),
     strip.text=element_text(hjust=0,face = "italic"),
     legend.position = "bottom",
     legend.title = element_blank())+
  ylab("Count")+xlab("Length (cm)")+
  geom_vline(aes(xintercept=Lm),lty=2)
gs1
ggsave("Figure5.jpg",gs1,height=9,width=6.5)
spsup<-GuatemalaAll$Species[GuatemalaAll$n>=20 & ! GuatemalaAll$Species %in% sptoplotG]
length(spsup)

gs4<-ggplot(dplyr::filter(Guatemala,sciname2 %in% spsup),
#  aes(x=Length,..density..,col=gear,fill=gear))+
   aes(x=Length,col=Gear,fill=Gear))+
  geom_histogram(position = "dodge2")+  
  facet_wrap(sciname2 ~.,scale="free",ncol=4)+
   theme_classic()+ theme(strip.background = element_blank(),
     strip.text=element_text(hjust=0,face = "italic",size=8),
     legend.position = "bottom",
     legend.title = element_blank())+
  ylab("Count")+xlab("Length (cm)")+
  geom_vline(aes(xintercept=Lm),lty=2)
gs4
ggsave("GuatemalaHistSup1.jpg",gs4,height=9,width=6.5)

pdf("GuatemalaHistograms.pdf",height=11,width=8.5)

for(i in 1:3) {
print(ggplot(filter(Guatemala,scinameFishbase %in%
    GuatemalaAll$Species[GuatemalaAll$n>=20]),
  aes(x=Length,(..count..)/sum(..count..),col=gear,fill=gear))+
  geom_histogram(position = "dodge2")+
  facet_wrap_paginate(sciname2 ~.,scale="free",nrow=5,ncol=3,page=i)+
  theme_bw()+
  ylab("Proportion"))
}
dev.off()


#Make additional parameter files for LBB and LBSPR
GuatemalaAll$MK<-GuatemalaAll$M/GuatemalaAll$K
GuatemalaAll$Loo<-GuatemalaAll$Linf
summary(GuatemalaAll)

GuatemalaAllFL<-GuatemalaAll
a<-match(GuatemalaAll$Species,Guatemalapardf$Species)
summary(a-1:length(a))
GuatemalaAllFL$Linf<-Guatemalapardf$Loo
GuatemalaAllFL$Loo<-Guatemalapardf$Loo
GuatemalaAllFL$Lm<-Guatemalapardf$Lm
GuatemalaAllFL$K<-Guatemalapardf$K
GuatemalaAllFL$M<-Guatemalapardf$M
GuatemalaAllFL$MK<-Guatemalapardf$MK


# At this point stop and run LBB

## Read in LBB results to set up files for LBSPR
theme_set(theme_classic())
lbbnoprior<-read.csv("NoneGuatemalaLBBoutTable.csv")
lbbFLprior<-read.csv("FLGuatemalaLBBoutTable.csv")
lbbdataprior<-read.csv("DataGuatemalaLBBoutTable.csv")
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

x<-GuatemalaGrowth %>% mutate(Linf=Loo,M=NA,Lm=NA) %>% dplyr::select(Species,Linf,K,M,Lm)
w<-GuatemalaMat %>% mutate(Linf=NA,K=NA,M=NA,MK=NA) %>% dplyr::select(Species,Linf,K,M,Lm)
y<-Guatemalapardf %>% mutate(Linf=Loomed,MK=MKmed,M=Mmed,K=Kmed) %>%
  dplyr::select(Species,M,K,MK,Linf,Lm)
z<-GuatemalaAll %>% mutate(MK=M/K) %>%
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

summary(GuatemalaAll$M/GuatemalaAll$K)

ggplot(GuatemalaAll,aes(x=M/K))+geom_histogram()
ggplot(GuatemalaAll,aes(x=log(M/K)))+geom_histogram()

# Make LBB data files for lbspr
GuatemalaAllLBBnoprior<-GuatemalaAll
a<-match(lbbnoprior$Species,GuatemalaAll$Species)
summary(a)
table(is.na(a))
GuatemalaAllLBBnoprior$Linf[]<-NA
GuatemalaAllLBBnoprior$MK[]<-NA
GuatemalaAllLBBnoprior$M[]<-NA
GuatemalaAllLBBnoprior$K[]<-NA
GuatemalaAllLBBnoprior$Linf[a]<-lbbnoprior$Linf.med
GuatemalaAllLBBnoprior$Loo<-GuatemalaAllLBBnoprior$Linf
GuatemalaAllLBBnoprior$MK[a]<-lbbnoprior$MK.med
GuatemalaAllLBBnoprior$K[a]<-GuatemalaAll$K[a]
GuatemalaAllLBBnoprior$M<-GuatemalaAllLBBnoprior$MK*GuatemalaAllLBBnoprior$K
summary(GuatemalaAllLBBnoprior)


GuatemalaAllLBBFLprior<-GuatemalaAll
a<-match(lbbFLprior$Species,GuatemalaAll$Species)
summary(a)
table(is.na(a))
GuatemalaAllLBBFLprior$Linf[]<-NA
GuatemalaAllLBBFLprior$MK[]<-NA
GuatemalaAllLBBFLprior$M[]<-NA
GuatemalaAllLBBFLprior$K[]<-NA
GuatemalaAllLBBFLprior$Linf[a]<-lbbFLprior$Linf.med
GuatemalaAllLBBFLprior$Loo<-GuatemalaAllLBBFLprior$Linf
GuatemalaAllLBBFLprior$MK[a]<-lbbFLprior$MK.med



