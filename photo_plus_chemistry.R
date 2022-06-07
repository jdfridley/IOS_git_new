#Photosynthesis data with leaf chemistry data
#JDF June 2022

library(lme4)
library(lmerTest)

dfold = "/Users/etudiant/Documents/IOS_git/IOS_data/"
hb = read.csv(paste0(dfold,"indoutHB4.csv"))
str(hb)
master = read.csv(paste0(dfold,"NSF-IOSspreadsheetJDF5.csv"))
str(master)
dat = merge(hb,master,by.x="ID",by.y="ID",all=T)
str(dat)
dat$GLI = as.numeric(as.character(dat$GLI.128))
dat$GLI[dat$GLI==1] = 100
dat$region = as.factor(substr((dat$Region),1,1))
covar = read.csv(paste0(dfold,"covariates.csv")) #covariates at the species-region level
dat = merge(dat,covar,by.x=c("species","region"),by.y=c("species","region"),all.x=T)

#inspect
sort(table(dat$sppcode))
sort(table(dat$species))
sort(table(dat$Region))
hist(dat$GLI)
  summary(dat$GLI)

#colors    
  makeTransparent<-function(someColor, alpha=150)
  {   newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)}) }
ecol = makeTransparent("darkgreen")
fcol = makeTransparent("blue")
jcol = makeTransparent("red2")
inv = rgb(0,0,0,alpha=0,max=255) #invisible  
trans = c(ecol,fcol,jcol)[as.numeric(dat$Region)]
solid = c("darkgreen","blue","red2")[as.numeric(dat$Region)]
  
#relationships

#1. Is alpha predicted by light level?
par(mar=c(4,4,1,1))
plot(log(dat$GLI[dat$alpha<1&dat$woody==1]),dat$alpha[dat$alpha<1&dat$woody==1],col=trans[dat$alpha<1&dat$woody==1],ylim=c(0,.8),pch=19,xaxt="n",ylab="Alpha",xlab="GLI")
points(log(dat$GLI[dat$alpha<1&dat$woody==1]),dat$alpha[dat$alpha<1&dat$woody==1],col=solid[dat$alpha<1&dat$woody==1],ylim=c(0,1))
abline(lsfit(log(dat$GLI[dat$alpha<1&dat$woody==1]),dat$alpha[dat$alpha<1&dat$woody==1]),lwd=2,col="gray")
xval = c(1,2,3,5,10,20,50,100)
axis(side=1,at=log(xval),labels=xval)
legend(legend=c("ENA","France","Japan"),x=-1,y=.8,bty="n",cex=.8,text.col=c("darkgreen","blue","red2"))
legend(legend=c("All GLI","GLI<100"),x=log(20),y=.8,bty="n",cex=.8,lty=c(1,2),col="gray",lwd=c(2,1))
abline(lsfit(log(dat$GLI[dat$alpha<1&dat$woody==1&dat$GLI<99]),dat$alpha[dat$alpha<1&dat$woody==1&dat$GLI<99]),lwd=1,lty=2,col="gray")

#linear model
lmer1 = lmer(alpha~GLI+(1|species),data=dat,subset=dat$alpha<1&dat$woody==1)
summary(lmer1)
  #very strong relationship when GLI=100 included
lmer2 = lmer(alpha~GLI+(1|species),data=dat,subset=dat$alpha<1&dat$woody==1&dat$GLI<99)
summary(lmer2)
  #no relationship when full sun excluded
lmer3 = lmer(alpha~GLI*region+(1|species),data=dat,subset=dat$alpha<1&dat$woody==1&dat$GLI<99)
summary(lmer3)
  #when GLI interaction included in the model, there is a marginal effect of J having lower alpha
lmer4 = lmer(alpha~GLI*region+(1|species),data=dat,subset=dat$alpha<1&dat$woody==1)
anova(lmer4)
  #ditto: regional effect is marginally sig even when GLI 100 included

#2. Is alpha predicted by chl?
par(mar=c(4,4,1,1))
plot(dat$Chlorophyll.a[dat$alpha<1&dat$woody==1],dat$alpha[dat$alpha<1&dat$woody==1],col=trans[dat$alpha<1&dat$woody==1],pch=19,ylim=c(0,.8),ylab="Alpha",xlab="Total Chl")
lmer1 = lmer(alpha~Total.Chlorophyll+(1|species),data=dat,subset=dat$alpha<1&dat$woody==1)
summary(lmer1)
  #significant but weak positive relationship, P=0.037
  #***seems strong issue with regional differences, esp. ENA

#3. Is Vcmax predicted by Rubisco? yes
par(mar=c(4,4,1,1),mfrow=c(2,1))
plot(dat$Rubisco.content[dat$woody==1],dat$Vcmax[dat$woody==1],col=trans[dat$woody==1],pch=19,ylim=c(0,150),ylab="Vcmax",xlab="Rubisco")
abline(lsfit(dat$Rubisco.content[dat$woody==1],dat$Vcmax[dat$woody==1]),col="gray",lwd=2)
mtext("Woody",side=3,line=-1)
plot(dat$Rubisco.content[dat$woody==0],dat$Vcmax[dat$woody==0],col=trans[dat$woody==0],pch=19,ylim=c(0,300),ylab="Vcmax",xlab="Rubisco")
abline(lsfit(dat$Rubisco.content[dat$woody==0],dat$Vcmax[dat$woody==0]),col="gray",lwd=2)
mtext("Herbaceous",side=3,line=-1)
lmer1 = lmer(Vcmax~Rubisco.content*region*woody+(1|species),data=dat)
anova(lmer1)
  #very strong Rubisco and woody/herb effect


#4. NUE ######################################
hist(dat$Asat)
dat$leafN = as.numeric(as.character(dat$leafN))
hist(dat$leafN)
dat$Specific.leaf.area = as.numeric(as.character(dat$Specific.leaf.area))
hist(dat$Specific.leaf.area)
  #two extreme values > 1000 cm2/g for Agrostis in Japan: real?
hist(dat$Specific.leaf.area[dat$Specific.leaf.area<1000])
dat$SLA.m2 = dat$Specific.leaf.area / 10000
dat$Narea = dat$leafN / dat$Specific.leaf.area
dat$Narea.m2 = dat$leafN / dat$SLA.m2
hist(dat$Narea)
dat$NUE = (dat$Asat * dat$SLA.m2) / (dat$leafN/100)
hist(dat$NUE[dat$NUE<100])
  #values between about 5 and 40 umol CO2 per g N per s, matches past data

par(mar=c(4,4,1,1),mfrow=c(2,1))
plot(dat$Narea.m2[dat$woody==1],dat$Asat[dat$woody==1],col=trans[dat$woody==1],pch=19,ylim=c(0,25),xlim=c(30,310),ylab="Asat",xlab="Narea (g N per m2)")
abline(lsfit(dat$Narea.m2[dat$woody==1&dat$Region=="ENA"],dat$Asat[dat$woody==1&dat$Region=="ENA"]),col="darkgreen",lwd=2)
abline(lsfit(dat$Narea.m2[dat$woody==1&dat$Region=="France"],dat$Asat[dat$woody==1&dat$Region=="France"]),col="blue",lwd=2)
abline(lsfit(dat$Narea.m2[dat$woody==1&dat$Region=="Japan"],dat$Asat[dat$woody==1&dat$Region=="Japan"]),col="red2",lwd=2)
mtext("Woody",side=3,line=-1)
plot(dat$Narea.m2[dat$woody==0],dat$Asat[dat$woody==0],col=trans[dat$woody==0],pch=19,ylim=c(0,50),xlim=c(30,310),ylab="Asat",xlab="Narea (g N per m2)")
abline(lsfit(dat$Narea.m2[dat$woody==0&dat$Region=="ENA"],dat$Asat[dat$woody==0&dat$Region=="ENA"]),col="darkgreen",lwd=2)
abline(lsfit(dat$Narea.m2[dat$woody==0&dat$Region=="France"],dat$Asat[dat$woody==0&dat$Region=="France"]),col="blue",lwd=2)
abline(lsfit(dat$Narea.m2[dat$woody==0&dat$Region=="Japan"],dat$Asat[dat$woody==0&dat$Region=="Japan"]),col="red2",lwd=2)
mtext("Herbaceous",side=3,line=-1)

#home-away contrasts: woodies from Japan to ENA
spp1 = c("Berberis thunbergii","Celastrus orbiculatus","Euonymus alatus","Lonicera japonica","Lonicera morrowii","Rosa multiflora","Viburnum dilatatum")
out1 = dat[is.element(dat$species,spp1),]
lmer1 = lmer(NUE~homeaway+(1|family/species),data=out1)
anova(lmer1)
summary(lmer1)
lmer1 = lmer(NUE~Region+(1|family/species),data=out1)
anova(lmer1)
summary(lmer1)
boxplot(NUE~homeaway,out1)
boxplot(NUE~Region,out1)
  #no homeaway or region differences in woodies from Japan to ENA

#home-away contrasts, ENA trees to France
spp2 = c("Prunus serotina","Quercus rubra","Parthenocissus sp.","Robinia pseudo-acacia","Acer negundo")
out2 = dat[is.element(dat$species,spp2),]
lmer2 = lmer(NUE~homeaway+(1|family/species),data=out2)
anova(lmer2)
summary(lmer2)
lmer2 = lmer(NUE~Region+(1|family/species),data=out2)
anova(lmer2)
summary(lmer2)
boxplot(NUE~homeaway,out2)
boxplot(NUE~Region,out2)
  #***very strong DECREASE in NUE in trees from ENA to France. Why???

#home-away contrasts, ENA herbs to Japan and France
spp3 = c("Ambrosia artemisiifolia","Bidens frondosa","Conyza canadensis","Erigeron annuus","Solidago gigantea")
out3 = dat[is.element(dat$species,spp3),]
lmer3 = lmer(NUE~homeaway+(1|species),data=out3)
anova(lmer3)
summary(lmer3)
lmer3 = lmer(NUE~Region+(1|species),data=out3)
anova(lmer3)
summary(lmer3)
boxplot(NUE~homeaway,out3)
boxplot(NUE~Region,out3)
  #no differences in ENA herbs to Japan and France, although Japan values much larger range

#home-away contrasts, France herbs to Japan and NY
spp4 = c("Anthoxanthum odoratum","Agrostis stolonifera","Chenopodium album","Daucus carota","Leucanthemum vulgaris","Plantago lanceolata","Artemisia vulgaris")
out4 = dat[is.element(dat$species,spp4),]
lmer4 = lmer(NUE~homeaway+(1|family/species),data=out4)
summary(lmer4)
anova(lmer4)
lmer4 = lmer(NUE~Region+(1|family/species),data=dat4)
summary(lmer4)
anova(lmer4)
boxplot(NUE~homeaway,out4,ylim=c(0,100))
boxplot(NUE~Region,out4,ylim=c(0,100))
  #****substantial sig increase from France to Japan, NOT France to ENA, but poor sample size in ENA (perhaps because still missing SLA?)

