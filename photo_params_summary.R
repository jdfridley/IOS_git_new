##summary of photo params fit via HB
#JDF 6-1-22

#note: still need final version of HB fit (with more iterations)

library(lme4)
library(lmerTest)
library(multcomp)

dfold = "/Users/etudiant/Documents/IOS_git/IOS_data/"

hb = read.csv(paste0(dfold,"indoutHB3.csv"))
str(hb)
nls = read.csv(paste0(dfold,"lightcurves_nls_output.csv"))
str(nls)
dat = merge(hb,nls,by.x="ID",by.y="sample",all.x=T)
str(dat)
dat$region = substr(dat$site,1,1)
  #1 is ENA, 2 is France, 3 is Japan
dat$Region = as.factor(dat$region)
#covar = read.csv(paste0(dfold,"covariates.csv")) #covariates at the species-region level

ep = read.csv(paste0(dfold,"photo_params_ecophys.csv")) #includes covars
out = merge(ep,dat,by.x="ID",by.y="ID",all=T)
out = droplevels(out)
out$region = out$region.x
out$species = out$species.x
out$sppcode = out$sppcode.x
out$away = out$homeaway=="away"

wout = out[out$woody==T,]
wout = droplevels(wout)

hout = out[out$woody==F,]
hout = droplevels(hout)


##################################################
#general summary plot
  #parameter, growth form, region (ENA=green, France=blue, Japan=red)

makeTransparent<-function(someColor, alpha=150)
{   newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)}) }

ecol = makeTransparent("darkgreen")
fcol = makeTransparent("blue")
jcol = makeTransparent("red2")
inv = rgb(0,0,0,alpha=0,max=255) #invisible

###Asat
par(mfrow=c(2,1),mar=c(3,3,1,0))
y = wout$Asat
yrange = c(0,60)
plot(wout$sppcode[wout$region=="E"],y[wout$region=="E"],border=ecol,ylim=yrange,las=3,col=inv)
mtext("Asat",side=2,line=2.5,cex=1.5)
par(new=T)
plot(wout$sppcode[wout$region=="F"],y[wout$region=="F"],border=fcol,ylim=yrange,las=3,col=inv)
par(new=T)
plot(wout$sppcode[wout$region=="J"],y[wout$region=="J"],border=jcol,ylim=yrange,las=3,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Woody","topright",bty="n",,cex=.7)
y = hout$Asat
plot(hout$sppcode[hout$region=="E"],y[hout$region=="E"],border=ecol,ylim=yrange,las=3,cex.axis=.8,col=inv)
mtext("Asat",side=2,line=2.5,cex=1.5)
par(new=T)
plot(hout$sppcode[hout$region=="F"],y[hout$region=="F"],border=fcol,ylim=yrange,las=3,cex.axis=.8,col=inv)
par(new=T)
plot(hout$sppcode[hout$region=="J"],y[hout$region=="J"],border=jcol,ylim=yrange,las=3,cex.axis=.8,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Herbaceous","topright",bty="n",,cex=.7)

###Vcmax
par(mfrow=c(2,1),mar=c(3,3,1,0))
y = wout$Vcmax.x
yrange = c(0,300)
plot(wout$sppcode[wout$region=="E"],y[wout$region=="E"],border=ecol,ylim=yrange,las=3,col=inv)
mtext("Vcmax",side=2,line=2.5,cex=1.5)
par(new=T)
plot(wout$sppcode[wout$region=="F"],y[wout$region=="F"],border=fcol,ylim=yrange,las=3,col=inv)
par(new=T)
plot(wout$sppcode[wout$region=="J"],y[wout$region=="J"],border=jcol,ylim=yrange,las=3,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Woody","topright",bty="n",,cex=.7)
y = hout$Vcmax.x
plot(hout$sppcode[hout$region=="E"],y[hout$region=="E"],border=ecol,ylim=yrange,las=3,cex.axis=.8,col=inv)
mtext("Vcmax",side=2,line=2.5,cex=1.5)
par(new=T)
plot(hout$sppcode[hout$region=="F"],y[hout$region=="F"],border=fcol,ylim=yrange,las=3,cex.axis=.8,col=inv)
par(new=T)
plot(hout$sppcode[hout$region=="J"],y[hout$region=="J"],border=jcol,ylim=yrange,las=3,cex.axis=.8,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Herbaceous","topright",bty="n",,cex=.7)

#home-away contrasts: woodies from Japan to ENA
spp1 = c("Berberis thunbergii","Celastrus orbiculatus","Euonymus alatus","Lonicera japonica","Lonicera morrowii","Rosa multiflora","Viburnum dilatatum")
out1 = out[is.element(out$species,spp1),]
lmer1 = lmer(Vcmax.x~homeaway+(1|family/species),data=out1)
anova(lmer1)
summary(lmer1)
  #!!Vcmax increases in the away range of J to E invasive shrubs, by about 10 units
boxplot(Vcmax.x~homeaway,data=out1)


###Jmax
par(mfrow=c(2,1),mar=c(3,3,1,0))
y = wout$Jmax.x
yrange = c(0,350)
plot(wout$sppcode[wout$region=="E"],y[wout$region=="E"],border=ecol,ylim=yrange,las=3,col=inv)
mtext("Jmax",side=2,line=2.5,cex=1.5)
par(new=T)
plot(wout$sppcode[wout$region=="F"],y[wout$region=="F"],border=fcol,ylim=yrange,las=3,col=inv)
par(new=T)
plot(wout$sppcode[wout$region=="J"],y[wout$region=="J"],border=jcol,ylim=yrange,las=3,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Woody","topright",bty="n",,cex=.7)
y = hout$Jmax.x
plot(hout$sppcode[hout$region=="E"],y[hout$region=="E"],border=ecol,ylim=yrange,las=3,cex.axis=.8,col=inv)
mtext("Jmax",side=2,line=2.5,cex=1.5)
par(new=T)
plot(hout$sppcode[hout$region=="F"],y[hout$region=="F"],border=fcol,ylim=yrange,las=3,cex.axis=.8,col=inv)
par(new=T)
plot(hout$sppcode[hout$region=="J"],y[hout$region=="J"],border=jcol,ylim=yrange,las=3,cex.axis=.8,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Herbaceous","topright",bty="n",,cex=.7)

###alpha
par(mfrow=c(2,1),mar=c(3,3,1,0))
y = wout$alpha
yrange = c(0,.8)
plot(wout$sppcode[wout$region=="E"],y[wout$region=="E"],border=ecol,ylim=yrange,las=3,col=inv)
mtext("alpha",side=2,line=2.5,cex=1.5)
par(new=T)
plot(wout$sppcode[wout$region=="F"],y[wout$region=="F"],border=fcol,ylim=yrange,las=3,col=inv)
par(new=T)
plot(wout$sppcode[wout$region=="J"],y[wout$region=="J"],border=jcol,ylim=yrange,las=3,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Woody","topright",bty="n",,cex=.7)
#no herbaceous

###Rd
par(mfrow=c(2,1),mar=c(3,3,1,0))
y = wout$Rd.x
yrange = c(0,3)
plot(wout$sppcode[wout$region=="E"],y[wout$region=="E"],border=ecol,ylim=yrange,las=3,col=inv)
mtext("Rd",side=2,line=2.5,cex=1.5)
par(new=T)
plot(wout$sppcode[wout$region=="F"],y[wout$region=="F"],border=fcol,ylim=yrange,las=3,col=inv)
par(new=T)
plot(wout$sppcode[wout$region=="J"],y[wout$region=="J"],border=jcol,ylim=yrange,las=3,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Woody","topright",bty="n",,cex=.7)
y = hout$Rd.x
plot(hout$sppcode[hout$region=="E"],y[hout$region=="E"],border=ecol,ylim=yrange,las=3,cex.axis=.8,col=inv)
mtext("Rd",side=2,line=2.5,cex=1.5)
par(new=T)
plot(hout$sppcode[hout$region=="F"],y[hout$region=="F"],border=fcol,ylim=yrange,las=3,cex.axis=.8,col=inv)
par(new=T)
plot(hout$sppcode[hout$region=="J"],y[hout$region=="J"],border=jcol,ylim=yrange,las=3,cex.axis=.8,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Herbaceous","topright",bty="n",,cex=.7)


###################################################
#compare Amax/Asat calcs (should be smaller for HB)
plot(dat$Asat,dat$Amax,col=as.numeric(dat$Region)); abline(0,1)
abline(lsfit(dat$Asat,dat$Amax),lty=2,col="gray")
summary(lm(Amax~Asat,dat))
  #R2=.72, slope = .55; HB fit is sig lower than nls
  #but no extreme outliers
lmer1 = lmer(Amax~Asat*Region+(1|species),data=dat)
anova(lmer1)
  #some differences in relationship based on region, too
summary(glht(lmer1, linfct = mcp(Region = 'Tukey')))
  #post hoc Tukey suggests regional difference is between ENA and France
boxplot(log(dat$Amax/dat$Asat) ~ dat$Region,ylim=c(-1,1))
  #difference between methods is greater in ENA; curious

###################################################
#compare alpha/AQY calcs
plot(dat$alpha,dat$AQY,col=as.numeric(dat$Region),xlim=c(0,.5)); abline(0,1)
abline(lsfit(dat$alpha,dat$AQY),lty=2,col="gray")
  #no relationship of AQY and alpha... ?? AQY via nls seems far too low (median = .06; typical is >.2)

plot(out$sppcode[out$region=="E"],out$alpha[out$region=="E"],col="green2",ylim=c(0,.8),las=3)
par(new=T)
plot(out$sppcode[out$region=="F"],out$alpha[out$region=="F"],col="blue",ylim=c(0,.8),las=3)
par(new=T)
plot(out$sppcode[out$region=="J"],out$alpha[out$region=="J"],col="red2",ylim=c(0,.8),las=3)

plot(out$sppcode[out$region=="E"],out$AQY[out$region=="E"],col="green2",ylim=c(0,.2),las=3)
par(new=T)
plot(out$sppcode[out$region=="F"],out$AQY[out$region=="F"],col="blue",ylim=c(0,.2),las=3)
par(new=T)
plot(out$sppcode[out$region=="J"],out$AQY[out$region=="J"],col="red2",ylim=c(0,.2),las=3)



#Amax/Asat comparisons using both methods
lmer2 = lmer(Amax~Region+(1|species),data=dat)
summary.aov(lmer2) #no overall regional difference in Amax
lmer3 = lmer(Amax~homeaway+(1|family/species),data=out,subset=out$region=="E")
anova(lmer3)
boxplot(Amax~homeaway,out,subset=out$region=="E")

#home-away contrasts: woodies from Japan to ENA
spp1 = c("Berberis thunbergii","Celastrus orbiculatus","Euonymus alatus","Lonicera japonica","Lonicera morrowii","Rosa multiflora","Viburnum dilatatum")
out1 = out[is.element(out$species,spp1),]
lmer6 = lmer(Amax~homeaway+(1|family/species),data=out1)
anova(lmer6)
boxplot(Amax~homeaway,data=out1)
boxplot(Asat~homeaway,data=out1)

