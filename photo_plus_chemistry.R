#Photosynthesis data with leaf chemistry data
#JDF June 2022

library(lme4)
library(lmerTest)

dfold = "/Users/etudiant/Documents/IOS_git/IOS_data/"
hb = read.csv(paste0(dfold,"indoutHB4.csv"))
str(hb)
master = read.csv(paste0(dfold,"NSF-IOSspreadsheetJDF4.csv"))
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



