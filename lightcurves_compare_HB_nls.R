##compare light curve parameters between HB and nls methods
#JDF 6-1-22

library(lme4)
library(lmerTest)
library(multcomp)

dfold = "/Users/etudiant/Documents/IOS_git/IOS_data/"

hb = read.csv(paste0(dfold,"indoutHB3.csv"))
str(hb)
nls = read.csv(paste0(dfold,"lightcurves_nls_output.csv"))
str(nls)
dat = merge(hb,nls,by.x="ID",by.y="sample")
str(dat)
dat$region = substr(dat$site,1,1)
  #1 is ENA, 2 is France, 3 is Japan
dat$Region = as.factor(dat$region)
#covar = read.csv(paste0(dfold,"covariates.csv")) #covariates at the species-region level

ep = read.csv(paste0(dfold,"photo_params_ecophys.csv")) #includes covars
out = merge(ep,dat,by.x="ID",by.y="ID",all=T)
#out = droplevels(out)

#licor data
dat1 = read.csv(paste0(dfold,"all_licor_data2.csv"))
str(dat1)
dat1$region = substr(dat1$site,1,1)
dat1$Date = as.Date(dat1$date,format="%m/%d/%y")
dat1$sppcode[dat1$sppcode=="amrt"] = "amart" #fix typo


##################################################
#compare curves visually for each species

sample = unique(dat1$ID)

i = 11
df = dat1[dat1$ID==sample[i],]
params = out[as.character(out$ID)==as.character(sample[i]),]

par(mar=c(1,1,1,0),mfrow=c(2,1))
##light curve
dat2 = df[df$CO2R>360&df$CO2R<440,]
plot(dat2$PARi,dat2$Photo,xlab="", ylab="", ylim=c(-2,max(dat2$Photo)+2),cex.lab=1.2,cex.axis=1.5,cex=1.5)
mtext(sample[i],side=3,line=1.5,cex=1.5)

curve((1/(2*params$theta*(summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1]-sqrt((summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1])^2-4*summary(curve.nlslrc)$coef[2,1]*summary(curve.nlslrc)$coef[4,1]*summary(curve.nlslrc)$coef[1,1]*x))-summary(curve.nlslrc)$coef[3,1],lwd=2,col="blue",add=T)


#points(PARlrc,photolrc,col=true,pch=19,cex=1.5)
points(PARlrc,photolrc,col=abs(stdres)>thresh,pch=19,cex=1.5)
mtext(expression("PPFD ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=2)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2,cex=2)
curve((1/(2*summary(curve.nlslrc)$coef[4,1]))*(summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1]-sqrt((summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1])^2-4*summary(curve.nlslrc)$coef[2,1]*summary(curve.nlslrc)$coef[4,1]*summary(curve.nlslrc)$coef[1,1]*x))-summary(curve.nlslrc)$coef[3,1],lwd=2,col="blue",add=T)
text(100,0,i)
#readline() #pause until any input




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

