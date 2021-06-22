#After data cleaning for Belize and Guatemala, combine with this code.

GuatemalaFam$Year<-as.numeric(GuatemalaFam$Year)
GuatemalaFam$gear<-paste0(" ",GuatemalaFam$Gear)
BelizeFam$gear<-BelizeFam$Gear
Guatemala$gear<-paste0(" ",Guatemala$Gear)
Belize$gear<-Belize$Gear
BelizeGear<-Belize[Belize$Gear!="Unknown",]
  
x<-names(GuatemalaFam)[c(59,61:68)]                          
y<-match(x,names(BelizeFam))
x[is.na(y)]
 
CombinedFam<-bind_rows(list(Guatemala=GuatemalaFam[,x],Belize=BelizeFam[,x]),.id="Source")
table(CombinedFam$Source,CombinedFam$gear)

table(CombinedFam$Source,CombinedFam$Habitat)

CombinedFam$Habitat<-factor(CombinedFam$Habitat)
write.csv(table(CombinedFam$Gear,CombinedFam$Source),"SampleSize.csv")

#

names(BelizeGear)[names(BelizeGear)=="FAMILY"]="Family"
BelizeGear$Length<-BelizeGear$TL2_CM
BelizeGear$Common<-BelizeGear$SPECIES2
x<-names(Guatemala)[c(59,61:72)]                          
y<-match(x,names(BelizeGear))
x[is.na(y)]

Guatemala$Year<-as.numeric(Guatemala$Year)
BelizeGear$Ltype<-rep("TL",dim(BelizeGear)[1])
 
Combined<-bind_rows(list(Guatemala=Guatemala[,x],Belize=BelizeGear[,x]),.id="Source")
table(Combined$Gear)

table(Combined$Source,Combined$Habitat)

Combined$Habitat<-factor(Combined$Habitat)
summary(Combined)
write.csv(table(CombinedFam$gear,CombinedFam$Habitat),"SampleSizeSpecies.csv")


Combined$Trophic2[Combined$Trophic2=="Invert-Piscivore"]="Invertivore-Piscivore"
Combined$Trophic2[Combined$Trophic2=="Herbivore-Invertivore"]="Invertivore"
table(Combined$Trophic2)

