#This code does all the multispecies analyses and multivariate indicators. 

library(vegan)
library(MASS)
library(lemon)
library(DHARMa)
library(grid)
library(ggplot2)
library(tidyverse)

region=c("Belize","Guatemala","Combined")[3]
if(region=="Belize") {
  datfile<-BelizeFam
  spfile<-BelizeAll
  datfile$length<-datfile$TL2_CM
  datfile$Species<-datfile$scinameFishbase
  datfile$gear<-datfile$Gear
  datfile<-filter(datfile,gear!="Unknown")
} 
if(region=="Guatemala") {
 datfile<-GuatemalaFam
 spfile<-GuatemalaAll
 datfile$length<-datfile$Length
 datfile$sciname<-datfile$scinameFishbase
 datfile$gear<-datfile$Gear
}
if(region=="Combined") {
 datfile<-CombinedFam
 datfile$sciname<-datfile$scinameFishbase
}
#Start analysis
if(region=="Combined") {
 x<-table(datfile$Family,datfile$Source)
 x[,1]<-x[,1]/sum(x[,1])
 x[,2]<-x[,2]/sum(x[,2])
 x<-x[x[,1]>0.01 | x[,2]>0.01,]
 x<-x[order(rowSums(x)),]
 write.csv(x[order(-rowSums(x)),],file="Figure2.csv")
 datfile$Family[!datfile$Family %in% dimnames(x)[[1]]]<-"Other"
 datfile$Family<-factor(datfile$Family, levels=c("Other",dimnames(x)[[1]]))
datfile$Type2<-datfile$Source
datfile$Type2[datfile$Type2=="Belize"]<-"(a) Belize"
datfile$Type2[datfile$Type2=="Guatemala"]<-"(b) Guatemala"
table(datfile$Type2)
df1<-datfile %>% group_by(Source,Type2,Family,Gear) %>%
  summarize(Count=n()) %>%
  mutate(Proportion=ifelse(Source=="Belize",Count/18533,Count/17080)) 

g4<-ggplot(df1,aes(y=Family,x=Proportion,fill=Gear))+
  geom_col()+ 
    xlab("Proportion")+
    facet_rep_wrap(Type2~.,ncol=3,scales="free_x")+
   theme_classic()+ theme(strip.background = element_blank(),
    panel.spacing.x = unit(-8, "mm"),
     strip.text=element_text(hjust=0))
 g4    #Fixed issue with summing to 1.0 across all facets
 ggsave("Figure2rev.jpg",g4,height=6,width=6.5)  
} else {
 x<-sort(table(datfile$Family)) 
 x
 write.csv(data.frame(x/sum(x)),paste0(region,"familycomp.csv"))
 y<-round(x/sum(x),3)
 x<-x[y>=0.01]
 x
 length(x)
 x
 datfile$Family[!datfile$Family %in% names(x)]="Other"
 datfile$Family<-factor(datfile$Family, levels=c("Other",names(x)))
}

 dim(datfile) #Total including other
 datfile<-droplevels(datfile[datfile$Family %in% names(x),])
 summary(datfile$Family)
 length(unique(datfile$Family))


### PERMANOVA and NMDS family leaving out station but including source if combined
if(region=="Combined") data1<- datfile %>%
  group_by(Year,Source,Month,gear,Family) %>% 
  summarize(Count=length(Family)) else
data1<- datfile %>%
  group_by(Year,Month,gear,Family) %>% 
  summarize(Count=length(Family))
dim(data1)
names(data1)
data2<-data1 %>% pivot_wider(id_cols=Year:gear,names_from = Family,values_from=Count,values_fn=sum)
head(data2)
dim(data2)
names(data2)
if(region=="Combined") spprop<-as.matrix(data2[,c(-1,-2,-3,-4)]) else
 spprop<-as.matrix(data2[,c(-1,-2,-3)])
spprop[is.na(spprop)]<-0
spprop<-spprop/rowSums(spprop)
summary(rowSums(spprop))
dim(spprop)
data2$Year<-factor(data2$Year)
data2$Month<-factor(data2$Month)
adonis2<-adonis(spprop~(gear+Month+Year)^2,data=data2)
adonis2
write.csv(round(adonis2$aov.tab,3),paste0(region,"permanova2.csv"))
dist2=vegdist(spprop,method="bray")
if(region=="Belize") {
 res2=metaMDS(spprop,k=3,distance="bray",trymax=4000)
 res2Belize<-res2
} 
if(region=="Guatemala") {
 res2=metaMDS(spprop,k=2,distance="bray",trymax=4000)
 res2Guatemala<-res2
}
  
res2
dev.off()
stressplot(res2,ylim=c(0,10))
data2Fam<-data2


# Make plots in ggplot
resnum<-2 #2 is family, 3 is trophic, 4 is species, 5 is size
if(resnum==2) {
  res=res2
  data2<-data2Fam
}  
if(resnum==3) {
  res=res3
  data2<-data2T
}  
if(resnum==4) {
  res=res4
  data2<-data2Sp
}  
if(resnum==5) {
  res=res5
  data2<-data2Size
}  
data2$group=data2$gear
data.scores <- as.data.frame(scores(res))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$group <-data2$group   #  add the grp variable created earlier
head(data.scores)  #look at the data
species.scores <- as.data.frame(scores(res, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data
hull.data<-NULL
groups<-unique(data2$group)
for(i in 1:length(groups)) {
 temp <- data.scores[data.scores$group == groups[i], ][chull(data.scores[data.scores$group == 
    groups[i], c("NMDS1", "NMDS2")]), ] 
 hull.data<-rbind(hull.data,temp)
} 

#
if(region=="Guatemala") {
  ghull.data<-hull.data
  gdata.scores<-data.scores
  gspecies.scores<-species.scores
}
#
if(region=="Belize") {
  bhull.data<-hull.data
  bdata.scores<-data.scores
  bspecies.scores<-species.scores
}
hull.data<-bind_rows(list(`(a) Belize`=bhull.data,`(b) Guatemala`=ghull.data),.id="Source")
data.scores<-bind_rows(list(`(a) Belize`=bdata.scores,`(b) Guatemala`=gdata.scores),.id="Source")
species.scores<-bind_rows(list(`(a) Belize`=bspecies.scores,`(b) Guatemala`=gspecies.scores),.id="Source")

#Hull plot

g1<-ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=group,group=group),alpha=0.30) + # add the convex hulls
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=group,colour=group),size=4) + # add the point markers
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  theme_classic(base_size=12) +
  theme(strip.background = element_blank(),
    panel.spacing.y = unit(8, "mm"),
     strip.text=element_text(hjust=0))+
  facet_wrap(Source~.,ncol=1,scales="free_x")+
  guides(fill=guide_legend(title="Gear"),color=guide_legend(title="Gear"),
    shape=guide_legend(title="Gear"))+
  scale_x_continuous(expand = c(.07, .07),breaks=scales::pretty_breaks(n=7))
g1
#ggsave("Fig5rev.jpg",g1,height=8,width=6.5)

# Species comp figure
region="Combined"
if(region=="Combined") {
  datfile<-Combined
  datfile$Species<-datfile$scinameFishbase
 }
table(datfile$gear)
if(region=="Combined") {
 x<-table(datfile$Species,datfile$Source)
 x[,1]<-x[,1]/sum(x[,1])
 x[,2]<-x[,2]/sum(x[,2])
 dim(x)
 write.csv(x[order(-rowSums(x)),],"Fig3.csv")
 round(x,2)
 x<-x[x[,1]>0.01 | x[,2]>0.01,]
 dim(x)
 x<-x[order(x[,1]+x[,2]),]
 datfile$Species[!datfile$Species %in% dimnames(x)[[1]]]<-"Other"
 datfile$Species[datfile$Species=="Coryphaena hippurus" & datfile$Source=="Belize"]<-"Other"
 datfile$Species<-factor(datfile$Species,levels=c("Other",dimnames(x)[[1]]))
}
datfile$Source[datfile$Source=="Guatemala"]<-"(b) Guatemala"
datfile$Source[datfile$Source=="Belize"]<-"(a) Belize"
df1<-datfile %>% group_by(Source,Species,Gear) %>%
  summarize(Count=n()) %>%
  mutate(Proportion=ifelse(Source=="(a) Belize",Count/18382,Count/15883)) 

g1<-ggplot(df1,aes(y=Species,fill=Gear,x=Proportion))+
  geom_col()+
  theme_classic(base_size=12)+
 facet_wrap(Source~.,scales="free",ncol=1)+
    theme( strip.text=element_text(hjust=0),
      strip.background = element_blank())+
  theme(axis.text.y = element_text(face="italic"))+xlab("Proportion")
g1
#ggsave("Fig3rev.jpg",g1,width=6.5,height=8)


### Ecosystem indicators (only with length and species id)
region=c("Belize","Guatemala","Combined")[3]
if(region=="Belize") {
 datfile<-filter(Belize,Gear !="Unknown")
 datfile$gear<-datfile$Gear
 datfile$length<-datfile$TL2_CM
} 
if(region=="Guatemala") {
 datfile<-Guatemala
 datfile$length<-datfile$Length
 datfile$gear<-datfile$Arte.de.pesca..fishing.gear.
}
if(region=="Combined") {
 datfile<-Combined
 datfile$length<-datfile$Length
 datfile<-filter(datfile,Gear !="Unknown")
 datfile$Species<-datfile$scinameFishbase
 datfile$gear<-datfile$Gear
 datfile$gear[datfile$Source=="Guatemala"]<-paste0(" ", datfile$gear[datfile$Source=="Guatemala"]) 
}
datfile<-filter(datfile,!is.na(datfile$length))
dim(datfile)
table(datfile$Source)
summary(datfile$Lm)
summary(datfile$length)
summary(datfile$length/datfile$Lm)
datfile$Lmat<-datfile$length/datfile$Lm
summary(datfile$Trophic)
x<-table(datfile$Trophic2,datfile$gear)
datfile$Piscivore<-ifelse(datfile$Trophic2=="Piscivore" &!is.na(datfile$Trophic2),1,0)
datfile$Invertivore<-ifelse(datfile$Trophic2=="Invertivore" &!is.na(datfile$Trophic2),1,0)
datfile$`Piscivore-Invertivore`<-ifelse(datfile$Trophic2=="Invertivore-Piscivore" &!is.na(datfile$Trophic2),1,0)
summary(datfile)
datfile<-dplyr::select(datfile,gear,Gear,Species,Habitat,Source,Station,Year, Month, Trophic,length,Lmax,Lmat,Trophic,Piscivore,Invertivore,`Piscivore-Invertivore`) 
table(datfile$Source)

datlong<-pivot_longer(datfile,cols=Trophic:`Piscivore-Invertivore`,names_to="Indicator",
   values_to="Value")
names(datfile)
names(datlong)
table(datlong$Indicator)
datlong$Indicator[datlong$Indicator=="Trophic"]<-"Trophic level"
datlong$Indicator[datlong$Indicator=="length"]<-"Length"
#datlong$Indicator<-factor(datlong$Indicator,levels=c("Length","Lmat","Lmax","Trophic level","Piscivore","Invertivore","Piscivore-Invertivore")
datlong$Indicator2<-datlong$Indicator
datlong$Indicator2[datlong$Indicator=="Lmat"]="L[mat](cm)"
datlong$Indicator2[datlong$Indicator=="Lmax"]="L[max](cm)"
datlong$Indicator2[datlong$Indicator=="Trophic level"]<-"Trophic~level"
datlong$Indicator2[datlong$Indicator=="Length"]<-"Length(cm)"
datlong$Indicator2<-factor(datlong$Indicator2)
summary(datlong$Indicator2)
datlong$Indicator2<-factor(datlong$Indicator2,
  levels=levels(datlong$Indicator2)[c(3,2,7,5,6,1,4)])
summary(datlong$Indicator2)

datlong<-bind_rows(Total=datlong,Gear=datlong,.id="Level")
datlong$gear3<-datlong$Gear
datlong$gear3[datlong$Level=="Total"]<-paste0(datlong$Source[datlong$Level=="Total"], " Total")
datlong$gear3[datlong$Source=="Guatemala"]<-paste0(" ",datlong$gear3[datlong$Source=="Guatemala"])

datlong$gear3<-factor(datlong$gear3,levels=sort(unique(datlong$gear3))[c(1,2,4,3,5,7,8,9,10,6)])
summary(datlong$gear3)

targetline<-data.frame(Indicator2=unique(datlong$Indicator2),target=c(1,rep(NA,6)))

x<-datlong %>% group_by(Source,gear3,Indicator) %>%
  summarize(mean=mean(Value),standard.error=standard.error(Value))
write.csv(x,"Fig7.csv")

#Country only
datlong %>% group_by(Source,Indicator) %>%
  summarize(mean=mean(Value),standard.error=standard.error(Value)) %>%
  arrange(Indicator)



datText <- data.frame(label = paste0("   (",letters[c(1:6)],")"),
  Indicator2=levels(datlong$Indicator2)[c(1:6)])
datText$Indicator2<-factor(datText$Indicator2,levels=levels(datlong$Indicator2))
datText<-datText[1:3,]

inds<-c("Lmat","Lmax","Trophic level")

g2<-ggplot(filter(datlong,Indicator %in% inds),
  aes(y=gear3,x=Value))+
  geom_violin(aes(fill=Source))+
 facet_rep_wrap(Indicator2~.,scales="free_x",
   strip.position="bottom",labeller=label_parsed)+
 geom_vline(data=filter(datlong,Indicator=="Lmat"),aes(xintercept=1),lty=2)+
  theme_classic()+ 
  theme(strip.background = element_blank(),
    strip.placement="outside",axis.line=element_line(),
    strip.text.x = element_text(size = 10),
    panel.border=element_blank(),
    panel.spacing.x = unit(-20, "mm"),
    panel.spacing.y = unit(1, "mm"),
    legend.position="none",
    plot.margin=unit(c(.1,1,0,.1),"cm"))+
    ylab("Gear")+xlab("")+
    guides(shape=FALSE)+
   geom_text( data    = datText,
      mapping = aes(x = -Inf, y = Inf, label = label), 
     hjust   = 0,  vjust   = 1)+
  geom_hline(yintercept=4.5,color="grey",lwd=2)
g2

jpeg("Figure7rev2.jpg",height=6,width=6.5,units="in",res=300)
g2
grid.text("Belize",x=unit(.98,"npc"),y=unit(0.75,"npc"),hjust=0.5,rot=270,gp=gpar(fontsize=10))
grid.text("Guatemala",x=unit(.98,"npc"),y=unit(0.3,"npc"),hjust=0.5,rot=270,gp=gpar(fontsize=10))
dev.off()


#Histograms

summary(datlong$Indicator2)
datlong$Indicator2<-factor(datlong$Indicator2,
  levels=levels(datlong$Indicator2)[c(7,1:6)])

datText <- data.frame(label = paste0(" (",letters[1:4],")"),
  Indicator2=levels(datlong$Indicator2)[1:4])
datText$Indicator2<-factor(datText$Indicator2,levels=levels(datlong$Indicator2)[1:4])
datText$Source=rep("Belize",4)

g4<-ggplot(filter(datlong,!Indicator %in% c("Piscivore","Invertivore","Piscivore-Invertivore")),
  aes(x=Value,y=..count../sum(..count..),fill=Source))+geom_histogram(position="dodge")+
 geom_histogram(position="dodge")  +
 facet_rep_wrap(Indicator2~.,scales="free",
   strip.position="bottom",labeller=label_parsed)+
 geom_vline(data=filter(datlong,Indicator=="Lmat"),aes(xintercept=1),lty=2)+
 theme_classic(base_size=12)+ 
 theme(strip.background = element_blank(),
    strip.placement="outside",axis.line=element_line(),
    panel.border=element_blank(),
    panel.spacing.y = unit(5, "mm"),
   strip.text.x = element_text(size = 11))+
    ylab("Proportion")+xlab("")+
    labs(fill="Country")+
     geom_text( data    = datText,
      mapping = aes(x = -Inf, y = Inf, label = label), 
     hjust   = 0,  vjust   = 1)+
   scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))

g4
#ggsave("Fig1s.jpg",g4,height=6,width=6.5)

# Various models
countryLM<-list()
indvals<-c("Length","Lmat","Lmax","Trophic",
  "Invertivore","Piscivore-Invertivore","Piscivore")
countryLM[[1]]<-lm(length~Source,data=datfile)
countryLM[[2]]<-lm(Lmat~Source,data=datfile)
countryLM[[3]]<-lm(Lmax~Source,data=datfile)
countryLM[[4]]<-lm(Trophic~Source,data=datfile)
countryLM[[5]]<-glm(Invertivore~Source,data=datfile,family="binomial")
countryLM[[6]]<-glm(`Piscivore-Invertivore`~Source,data=datfile,family="binomial")
countryLM[[7]]<-glm(Piscivore~Source,data=datfile,family="binomial")
countryLM[[8]]<-lm(log(length)~Source,data=datfile)
countryLM[[9]]<-lm(log(Lmat)~Source,data=datfile)
countryLM[[10]]<-lm(log(Lmax)~Source,data=datfile)
countryLM[[11]]<-lm(log(Trophic)~Source,data=datfile)
countryLMtab<-data.frame(Indicator=c(indvals,paste("log",indvals[1:4])))
newdat<-data.frame(Source=c("Belize","Guatemala"))
for(i in c(1:4,8:11)) {
  countryLMtab$P[i]<-anova(countryLM[[i]])["Source","Pr(>F)"]
  temp<-predict(countryLM[[i]],newdata=newdat,se.fit=TRUE)
  countryLMtab$Belize[i]<-temp$fit[1]
  countryLMtab$Guatemala[i]<-temp$fit[2]
  countryLMtab$BelizeSE[i]<-temp$se.fit[1]
  countryLMtab$GuatemalaSE[i]<-temp$se.fit[2]
}  
for(i in 5:7) {
  countryLMtab$P[i]<-anova(countryLM[[i]],test="Chi")["Source","Pr(>Chi)"]
  temp<-predict(countryLM[[i]],newdata=newdat,se.fit=TRUE,type="response")
  countryLMtab$Belize[i]<-temp$fit[1]
  countryLMtab$Guatemala[i]<-temp$fit[2]
  countryLMtab$BelizeSE[i]<-temp$se.fit[1]
  countryLMtab$GuatemalaSE[i]<-temp$se.fit[2]
}  


countryLMtab
# For large sample size, can assume normality. 



lmmods<-list()
lmmods[[1]]=stepAIC(lm(length~(Source+Year+Month+Gear)^2,data=datfile))
lmmods[[2]]=stepAIC(lm(Lmat~(Source+Year+Month+Gear)^2,data=datfile))
lmmods[[3]]=stepAIC(lm(Lmax~(Source+Year+Month+Gear)^2,data=datfile))
lmmods[[4]]=stepAIC(lm(Trophic~(Source+Year+Month+Gear)^2,data=datfile))
lmmods[[5]]=stepAIC(glm(Invertivore~(Source+Year+Month+Gear)^2,data=datfile))
lmmods[[6]]=stepAIC(glm(`Piscivore-Invertivore`~(Source+Year+Month+Gear)^2,data=datfile))
lmmods[[7]]=stepAIC(glm(Piscivore~(Source+Year+Month+Gear)^2,data=datfile))
lmmods[[8]]=stepAIC(lm(log(length)~(Source+Year+Month+Gear)^2,data=datfile))
lmmods[[9]]=stepAIC(lm(log(Lmat)~(Source+Year+Month+Gear)^2,data=datfile))
lmmods[[10]]=stepAIC(lm(log(Lmax)~(Source+Year+Month+Gear)^2,data=datfile))
lmmods[[11]]=stepAIC(lm(log(Trophic)~(Source+Year+Month+Gear)^2,data=datfile))

aovvals<-list()
for(i in 1:11) {
  if(class(lmmods[[i]])[1]=="glm") {
    aovvals[[i]]<-anova(lmmods[[i]],test="Chi")
    aovvals[[i]]$DevianceExplained=round(aovvals[[i]]$Deviance/aovvals[[i]]$`Resid. Dev`[1],2)
  }     else {
    aovvals[[i]]<-anova(lmmods[[i]])
    aovvals[[i]]$DevianceExplained=round(aovvals[[i]]$`Sum Sq`/sum(aovvals[[i]]$`Sum Sq`),2)
  }
}
aovvals  
dev.exp<-function(x) {
  sum(x$Deviance,na.rm=TRUE)/x$`Resid. Dev`[1]
}

mvtab<-data.frame(Indicator=c(indvals,paste0("log",indvals[1:4])))
for(i in 1:dim(mvtab)[1]) 
  mvtab$BestModel[i]<-deparse1(formula(lmmods[[i]]))
for(i in c(1:4,8:11)) mvtab$DevianceExplained[i]<-summary(lmmods[[i]])$r.squared
for(i in 5:7) mvtab$DevianceExplained[i]<-dev.exp(anova(lmmods[[i]],test="Chi"))
mvtab

a<-data.frame(mvtab)
dim(a)
a$DevianceExplained<-round(a$DevianceExplained,2)0
a$Belize<-paste0(round(countryLMtab$Belize,2)," (",round(countryLMtab$BelizeSE,3),")")
a$Guatemala<-paste0(round(countryLMtab$Guatemala,2)," (",round(countryLMtab$GuatemalaSE,3),")")

a

write.csv(a[1:7,],"Combinedaovtab.csv")
countryLMtab

# Check residuals 
#lmmods is stepAIC for both countries multivariate
par(mfrow=c(4,2),mar=c(4,4,1,1))
for(i in 1:11) {
 # if(class(lmmods[[i]])[1]=="glm") {
 #   x<-simulateResiduals(lmmods[[i]],plot=TRUE)
 #   mtext(indvals[i],3)
 # }
#   else 
     plot(lmmods[[i]],1:2,main=indvals[i])
}

#Combined summaries

x<-table(Combined$Length>Combined$Lm)

x<-table(Combined$Length>Combined$Lm,Combined$Source)
x/rbind(colSums(x),colSums(x))

x<-Combined %>% group_by(scinameFishbase,Source) %>%
  summarize(mature=length(Length[Length>Lm]),
    total=n()) %>%
  mutate(fraction=mature/total) %>%
  filter(total>=15) %>%
  group_by(Source) %>%
  summarize(immature=length(fraction[fraction<0.5]),total=n()) %>%
  mutate(fraction=immature/total)
x 

length(unique(GuatemalaFam$Family))
length(unique(BelizeFam$Family))
length(unique(GuatemalaAll$Family))
length(unique(BelizeAll$Family))


length(unique(Combined$scinameFishbase))
length(unique(Belize$scinameFishbase))
length(unique(Guatemala$scinameFishbase))
90+69

x<-table(Combined$scinameFishbase,Combined$Source)
dimnames(x[x[,1]>0&x[,2]>0,])


