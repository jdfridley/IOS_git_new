#Comparing regional protein values
#JDF 5-24-22

dat = read.csv("NSF-IOSspreadsheetJDF4.csv")
str(dat)
dat$Water.soluble.Protein.mg.cm2 = as.numeric(dat$Water.soluble.Protein.mg.cm2)

boxplot(dat$Water.soluble.Protein.mg.cm2~Region,dat,ylim=c(0,2))
boxplot(dat$SDS.soluble.Protein.mg.cm2~Region,dat,ylim=c(0,2))

boxplot(dat$Water.soluble.Protein.mg.cm2~Region,dat,ylim=c(0,2),subset=dat$Species.abbreviation=="fajap")
boxplot(dat$SDS.soluble.Protein.mg.cm2~Region,dat,ylim=c(0,2),subset=dat$Species.abbreviation=="fajap")

spp = unique(dat$Species.abbreviation)

woody = c("acneg","bethu","ceorb","cesca","euala","liben","locan","lomor","paspp","pravi","prser","rhcat","ropse","romul","viace","lojap","qurub","vidil","pains","qurob","acpse","hulup","laana")
dat$woody = rep(0,dim(dat)[1])
dat$woody[is.element(dat$Species.abbreviation,woody)] = 1

for(i in 1:length(spp)) {
  par(mfrow=c(1,2))
  boxplot(dat$Water.soluble.Protein.mg.cm2~Region,dat,ylim=c(0,2),subset=dat$Species.abbreviation==spp[i])
  boxplot(dat$SDS.soluble.Protein.mg.cm2~Region,dat,ylim=c(0,2),subset=dat$Species.abbreviation==spp[i])
  title(main=spp[i])
  readline()
}

library(lme4)
library(lmerTest)

#water sol protein
lme1 = lmer(Water.soluble.Protein.mg.cm2~Region*woody+(1|Species.abbreviation),data=dat)
summary(lme1)
anova(lme1)              
  #no regional difference
              
#SDS sol protein
lme2 = lmer(SDS.soluble.Protein.mg.cm2~Region-1+(1|Species.abbreviation),data=dat)
summary(lme2)
anova(lme2)              
fixef(lme2)
  #regional difference

#water plus SDS
lme3 = lmer(Water.soluble.Protein.mg.cm2 + SDS.soluble.Protein.mg.cm2~Region+(1|Species.abbreviation),data=dat)
summary(lme3)
anova(lme3)   
fixef(lme3)
  #no regional difference

plot(dat$Water.soluble.Protein.mg.cm2,dat$SDS.soluble.Protein.mg.cm2,col=as.numeric(as.factor(dat$woody)),xlim=c(0,2),ylim=c(0,2),cex=1.8,lwd=2)
abline(0,1)


