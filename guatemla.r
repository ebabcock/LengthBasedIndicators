
## Guatemala
guatemala<-read.csv("DB_Pacific_Guatemala_WCS_2020_26_03 (Final).csv")
summary(guatemala)
guatemala$Scientific.name<-trimws(guatemala$Scientific.name)
gsp<-read.csv("Guatemala species list.csv")
gsp$Scientific.name<-trimws(gsp$Scientific.name)
x<-match(guatemala$Scientific.name,gsp$Scientific.name)
summary(x)
sort(table(guatemala$Scientific.name[is.na(x)]))
guatemala$Common.name.english<-gsp$Common..English.[x]

x<-sort(table(guatemala$Common.name.english))
x<-x[x>=40]
x

x<-sort(table(guatemala$Scientific.name))
x<-x[x>=40]
x
y<-c("Penaeus vannamei","Panulirus gracilis","Larimus spp.","Mugill spp.","Ariopsis spp1.")
x<-x[! names(x) %in% y]
length(x)
guatemalasp<-data.frame(Scientific.name=names(x))
guatemalasp$Scientific.name2<-guatemalasp$Scientific.name
guatemalasp$Scientific.name2[guatemalasp$Scientific.name2=="Isopisthus ramifer"]<-"Isopisthus remifer"
for(i in 1:length(x)) {
  guatemalasp$Genus[i]=unlist(strsplit(guatemalasp$Scientific.name2[i]," "))[1]
  guatemalasp$species[i]=unlist(strsplit(guatemalasp$Scientific.name2[i]," "))[2]
}
guatemalasp$Genus[guatemalasp$Scientific.name=="Hypanus longus"]<-"Dasyatis"
guatemalasp$species[guatemalasp$Scientific.name=="Hypanus longus"]<-"longa"
x<-match(guatemalasp$Scientific.name,guatemala$Scientific.name)
summary(x)
guatemalasp$Common<-guatemala$Common.name.english[x]

guatemalasp$Scientific.name
guafb<-popchar(guatemalasp$Scientific.name2)
x<-match(guatemalasp$Scientific.name2,guafb$Species)
summary(x)
guatemalasp$Scientific.name2[is.na(x)]
guatemalasp$Lmax<-rep(NA,dim(guatemalasp)[1])
for(i in 1:dim(guatemalasp)[1]) 
  guatemalasp$Lmax[i]<-max(guafb$Lmax[guafb$Species==guatemalasp$Scientific.name2[i]],na.rm=TRUE)
x<-match(guatemalasp$Scientific.name,gsp$Scientific.name)
summary(x)
guatemalasp$Lmax[!is.na(x)]<-gsp$Lmax[x[!is.na(x)]]
guatemalasp$Lm[!is.na(x)]<-gsp$Lm[x[!is.na(x)]]
summary(guatemalasp)
for(i in 1:dim(guatemalasp)[1]) 
  guatemalasp$n[i]<-length(guatemala$Scientific.name[guatemala$Scientific.name==
                                                       guatemalasp$Scientific.name[i] & !is.na(guatemala$Total.length.cm)])
write.csv(guatemalasp,file="guatemalasp.csv")

params<-c("Lm","M","K","Loo")
guatemalapars<-list()
for(i in 1:dim(guatemalasp)[1]) {
  guatemalapars[[i]]<-Plot_taxa( Search_species(Genus=guatemalasp$Genus[i],
                                                Species=guatemalasp$species[i])$match_taxonomy )[[1]]
}

