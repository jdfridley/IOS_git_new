##########N allocation
##2024 revision
#Fridley June 2024

################
#Datasets

#import trait and covariate data
dat = read.csv("C:\\Users\\fridley\\OneDrive - Clemson University\\academic\\projects\\IOS_FranceJapan\\NSFIOSspreadsheetJDF11.csv")
cov = read.csv("C:\\Users\\fridley\\OneDrive - Clemson University\\academic\\projects\\IOS_FranceJapan\\covariates.csv")
names(cov)[1:2] = c("Species","Region")
cov$Region = as.factor(cov$Region)
levels(cov$Region) = c("ENA","France","Japan")
dat2 = merge(dat,cov,by.x=c("Species","Region"),by.y=c("Species","Region"),all.x=T)
dat = dat2

#import photosynthesis parameters (from photo_params.R)
hb = read.csv("C:\\Users\\fridley\\OneDrive - Clemson University\\academic\\projects\\IOS_FranceJapan\\IOS_git2\\indoutHB6test.csv")
dat3 = merge(dat,hb,by.x="ID",by.y="ID",all=T)
dat3$region = as.factor(substr((dat3$Region),1,1))
dat = dat3[!is.na(dat3$region),]
dat = dat[!is.na(dat$woody),]
#dat = dat3
  #442 unique rows, but only 414 have photosynthetic data

#if posteriors from photo param model needed:
#load("jags.p_5.31.24.RData") #JAGS object saved, jags.p
#attach(jags.p)
  #Vcm.out, Jm.out, Rd.int, Asat

#make P. inserta just Parthenocissus sp. due to taxonomic questions about this group
dat$Species[dat$Species=="Parthenocissus inserta"] = "Parthenocissus sp." 

#fix dates
dat$Date = as.Date(dat$Sampling.Date,format="%m/%d/%Y")
dat$month = substr(dat$Date,7,7)

#grouping variable for nat-inv contrasts
dat$group = "native"
dat$group[dat$invader=="invader"&dat$homeaway=="away"] = "invaway"
dat$group[dat$invader=="invader"&dat$homeaway=="home"] = "invhome"

#add in predicted A values at low light for forest species
load("C:/Users/fridley/OneDrive - Clemson University/academic/projects/IOS_FranceJapan/IOS_git2/qmatout.RData")
  #data frame is 'd3'
d4 = d3[,c("ID","A.q20","A.q50","A.q100","A.q200")]
dat = merge(dat,d4,by.x="ID",by.y="ID",all.x=T)


###################
#Variables

dat$Nmass = (as.numeric(as.character(dat$leafN))/100) #g leaf N per g leaf
dat$SLA = as.numeric(dat$Specific.leaf.area)
  dat$SLA[dat$SLA>1000] = NA #two nonsensical values for Agrostis



###################
#1. PNUE modeling
  #Step 1: calculate (and statistically compare) SLA and Nmass
  #Step 2: use Asat, SLA, and Nmass to calculate PNUE
  #Step 3: model PNUE based on same process model
    #- this has to be separate run because JAGS does not allow PNUE to be defined twice
library(R2jags)
  
dfN = dat

#covariates
Nmass = dfN$Nmass*100
Asat = dfN$Asat
Alow = dfN$A.q100 #A at 100 PPFD, woodies only
group = as.numeric(as.factor(dfN$group)) #1=invaway 2=invhome 3=native
species = as.numeric(as.factor(dfN$Species)) #48 species
Nspp = max(species)
site = as.numeric(as.factor(dfN$region))
N = length(group)
woody = as.numeric(dfN$woody==1)
SLA = dfN$SLA #also log(dfN$SLA)
aqy = dfN$alpha.emp
away = as.numeric(dfN$homeaway=="away") #1 = away, 0 = home
native = as.numeric(dfN$invader=="native") #1 = native, 0 = invader
invader = as.numeric(dfN$invader=="invader") #1 = invader, 0 = native

#covariates at the species level
sppdf = dfN[duplicated(dfN$Species)==F,]
sppdf = sppdf[order(sppdf$Species),]
speciesS = c(1:dim(sppdf)[1])
woodyS = as.numeric(sppdf$woody==1)
invaderS = as.numeric(sppdf$invader=="invader")
groupS = as.numeric(as.factor(sppdf$group)) #1=invaway 2=invhome 3=native


#model and compare SLA and Nmass, calculate PNUE
mod.n <- "model
{
    #estimate Asat missing values
    for(i in 1:N) {
      asat[i] ~ dnorm(mean.a,tau.a)
      alow[i] ~ dnorm(mean.al,tau.al)
    }
    mean.a ~ dnorm(0,0.001)
    tau.a <- sigma.a^-2
    sigma.a ~ dunif(0,100) #prior for sd
    mean.al ~ dnorm(0,0.001)
    tau.al <- sigma.al^-2
    sigma.al ~ dunif(0,100) #prior for sd

  for(i in 1:N) {

    #SLA
    sla[i] ~ dnorm(mu.sla[i],tau)
    #separate base model predictors according to whether species is native or invasive
    mu.sla[i] <- ifelse(group[i] < 3,
        b0 + bW*woody[i] + bA[woody[i]+1]*away[i] + bSpp[species[i]] + bS[site[i]]
    ,
        b0 + bW*woody[i] + bSpp[species[i]] + bS[site[i]]
          #+ bI[woody[i]+1]*invader[i] #this part is now in species-level loop
    )

    #Nmass
    nmass[i] ~ dnorm(mu.nmass[i],tauN)
    #separate base model predictors according to whether species is native or invasive
    mu.nmass[i] <- ifelse(group[i] < 3,
        b0N + bWN*woody[i] + bAN[woody[i]+1]*away[i] + bSppN[species[i]] + bSN[site[i]]
    ,
        b0N + bWN*woody[i] + bSppN[species[i]] + bSN[site[i]]
    )

    #calculate PNUE, based on Asat and Narea
    narea[i] <- (nmass[i]/sla[i])*100 #save Narea for graphing results (g per m2)
    pnue[i] <- asat[i] / narea[i] #saturating light
    pnue100[i] <- alow[i] / narea[i] #at PPFD 100

}

  #priors for fixed effect slopes (groups)
  bW ~ dnorm(0, 0.001)
  b0 ~ dnorm(0, 0.001)
  bWN ~ dnorm(0, 0.001)
  b0N ~ dnorm(0, 0.001)

  for(i in 1:2) {
    bA[i] ~ dnorm(0, 0.001)
    bI[i] ~ dnorm(0, 0.001)
    bAN[i] ~ dnorm(0, 0.001)
    bIN[i] ~ dnorm(0, 0.001)
  }

  #random effect priors
  for(i in 1:3) {
    bS[i] ~ dnorm(0,tauS1) #site random effect
    bSN[i] ~ dnorm(0,tauS1N) #site random effect
  }

  #species REs
  for(i in 1:Nspp) {
    bSpp[i] ~ dnorm(mu.bSpp[i],tauSpp)
    mu.bSpp[i] <- bI[woodyS[i]+1]*invaderS[i] #invader effect must be here, else lost in species RE
  
    bSppN[i] ~ dnorm(mu.bSppN[i],tauSppN)
    mu.bSppN[i] <- bIN[woodyS[i]+1]*invaderS[i] #invader effect must be here, else lost in species RE
  }
  
  #priors, variances
  tau <- sigma^-2 
  sigma ~ dunif(0, 100) 
  tauSpp <- sigmaSpp^-2 
  sigmaSpp ~ dunif(0, 100) 
  tauS1 <- sigmaS1^-2 
  sigmaS1 ~ dunif(0, 100) 

  tauN <- sigmaN^-2 
  sigmaN ~ dunif(0, 100) 
  tauSppN <- sigmaSppN^-2 
  sigmaSppN ~ dunif(0, 100) 
  tauS1N <- sigmaS1N^-2 
  sigmaS1N ~ dunif(0, 100) 

}" #end model
write(mod.n, "model.txt")
  
#input lists for JAGS
params = c("bS","bW","bI","bSpp","b0","bA","bSN","bWN","bIN","bSppN","b0N","bAN","pnue","pnue100","narea","sla") #parameters to monitor
inits = function() list(b0N=(1.5),b0=5) #starting values of fitted parameters
input = list(N=N,sla=SLA,nmass=Nmass,asat=Asat,alow=Alow,group=group,site=site,woody=woody,species=species,Nspp=Nspp,away=away,woodyS=woodyS,invaderS=invaderS) #input data
  
#run JAGS model
jagsY <- jags(model = "model.txt",data = input,param=params, inits=inits,
               n.chains = 3, #number of separate MCMC chains 3
               n.iter =5000, #number of iterations per chain; 10000
               n.thin=5, #thinning 5
               n.burnin = 1000) #number of initial iterations to discard 2000
jagsY 
plot(jagsY)
#detach.jags()
attach.jags(jagsY)
#traceplot(jagsY)
pnue.out = colMeans(pnue)
pnue100.out = colMeans(pnue100)
Narea.out = colMeans(narea) #for plotting
sla.out = colMeans(sla) #to save for dataset
sla.sd = apply(sla,2,sd) #to save for dataset

#add estimated SLA and Narea to dat
dat$SLA2 = sla.out
dat$SLA2.se = sla.sd
dat$Narea2 = Narea.out

#calculate SLA and Narea differences by group

    #SLA
    sla.inv.away.herb = b0 + bI[,1] + bA[,1]
    sla.inv.away.woody = b0 + bI[,2] + bA[,2] + bW
    sla.inv.home.herb = b0 + bI[,1] 
    sla.inv.home.woody = b0 + bI[,2] + bW
    sla.nat.home.herb = b0  
    sla.nat.home.woody = b0 + bW
    sla.comp = data.frame(sla.inv.away.herb,sla.inv.away.woody,sla.inv.home.herb,sla.inv.home.woody,sla.nat.home.herb,sla.nat.home.woody)
    boxplot(sla.comp[,c(1,3,5)],outline=F)
    boxplot(sla.comp[,c(2,4,6)],outline=F)
    #credible intervals: for graphing
    sla.awayhome.h = sla.comp[,1]-sla.comp[,3]
    sla.awaynat.h = sla.comp[,1]-sla.comp[,5]
    sla.homenat.h = sla.comp[,3]-sla.comp[,5]
    x = sla.awayhome.h; length(x[x<0])/length(x) #P=0.515
    x = sla.awaynat.h; length(x[x<0])/length(x) #P=0.330
    x = sla.homenat.h; length(x[x<0])/length(x) #P=0.300
    sla.awayhome.w = sla.comp[,2]-sla.comp[,4]
    sla.awaynat.w = sla.comp[,2]-sla.comp[,6]
    sla.homenat.w = sla.comp[,4]-sla.comp[,6]
    x = sla.awayhome.w; length(x[x>0])/length(x) #P=0.0092**away shift to lower SLA in woodies
    x = sla.awaynat.w; length(x[x>0])/length(x) #P=0.319
    x = sla.homenat.w; length(x[x>0])/length(x) #P=0.776

    #Nmass: no differences between groups
    nmass.inv.away.herb = b0N + bIN[,1] + bAN[,1]
    nmass.inv.away.woody = b0N + bIN[,2] + bAN[,2] + bWN
    nmass.inv.home.herb = b0N + bIN[,1] 
    nmass.inv.home.woody = b0N + bIN[,2] + bWN
    nmass.nat.home.herb = b0N  
    nmass.nat.home.woody = b0N + bWN
    nmass.comp = data.frame(nmass.inv.away.herb,nmass.inv.away.woody,nmass.inv.home.herb,nmass.inv.home.woody,nmass.nat.home.herb,nmass.nat.home.woody)
    boxplot(nmass.comp[,c(1,3,5)],outline=F)
    boxplot(nmass.comp[,c(2,4,6)],outline=F)
    #credible intervals: for graphing
    nmass.awayhome.h = nmass.comp[,1]-nmass.comp[,3]
    nmass.awaynat.h = nmass.comp[,1]-nmass.comp[,5]
    nmass.homenat.h = nmass.comp[,3]-nmass.comp[,5]
    x = nmass.awayhome.h; length(x[x>0])/length(x) #P=0.463
    x = nmass.awaynat.h; length(x[x>0])/length(x) #P=0.304
    x = nmass.homenat.h; length(x[x>0])/length(x) #P=0.321
    nmass.awayhome.w = nmass.comp[,2]-nmass.comp[,4]
    nmass.awaynat.w = nmass.comp[,2]-nmass.comp[,6]
    nmass.homenat.w = nmass.comp[,4]-nmass.comp[,6]
    x = nmass.awayhome.w; length(x[x>0])/length(x) #P=0.434
    x = nmass.awaynat.w; length(x[x>0])/length(x) #P=0.413
    x = nmass.homenat.w; length(x[x>0])/length(x) #P=0.432
    

#######now run same process model for PNUE and Aarea
mod.n <- "model
{
  for(i in 1:N) {

    #PNUE
    pnue[i] ~ dnorm(mu.pnue[i],tau)
    #separate base model predictors according to whether species is native or invasive
    mu.pnue[i] <- ifelse(group[i] < 3,
        b0 + bW*woody[i] + bA[woody[i]+1]*away[i] + bSpp[species[i]] + bS[site[i]]
    ,
        b0 + bW*woody[i] + bSpp[species[i]] + bS[site[i]]
    )

    #Narea
    narea[i] ~ dnorm(mu.narea[i],tauN)
    #separate base model predictors according to whether species is native or invasive
    mu.narea[i] <- ifelse(group[i] < 3,
        b0N + bWN*woody[i] + bAN[woody[i]+1]*away[i] + bSppN[species[i]] + bSN[site[i]]
    ,
        b0N + bWN*woody[i] + bSppN[species[i]] + bSN[site[i]]
    )
}

  #priors for fixed effect slopes (groups)
  bW ~ dnorm(0, 0.001)
  b0 ~ dnorm(0, 0.001)
  bWN ~ dnorm(0, 0.001)
  b0N ~ dnorm(0, 0.001)

  for(i in 1:2) {
    bA[i] ~ dnorm(0, 0.001)
    bI[i] ~ dnorm(0, 0.001)
    bAN[i] ~ dnorm(0, 0.001)
    bIN[i] ~ dnorm(0, 0.001)
  }

  #random effect priors
  for(i in 1:3) {
    bS[i] ~ dnorm(0,tauS1) #site random effect
    bSN[i] ~ dnorm(0,tauS1N) #site random effect
  }

  #species REs
  for(i in 1:Nspp) {
    bSpp[i] ~ dnorm(mu.bSpp[i],tauSpp)
    mu.bSpp[i] <- bI[woodyS[i]+1]*invaderS[i] #invader effect must be here, else lost in species RE

    bSppN[i] ~ dnorm(mu.bSppN[i],tauSppN)
    mu.bSppN[i] <- bIN[woodyS[i]+1]*invaderS[i] #invader effect must be here, else lost in species RE
  }
  
  #priors, variances
  tau <- sigma^-2 
  sigma ~ dunif(0, 100) 
  tauSpp <- sigmaSpp^-2 
  sigmaSpp ~ dunif(0, 100) 
  tauS1 <- sigmaS1^-2 
  sigmaS1 ~ dunif(0, 100) 
  
  tauN <- sigmaN^-2 
  sigmaN ~ dunif(0, 100) 
  tauSppN <- sigmaSppN^-2 
  sigmaSppN ~ dunif(0, 100) 
  tauS1N <- sigmaS1N^-2 
  sigmaS1N ~ dunif(0, 100) 

}" #end model
write(mod.n, "model.txt")

#input lists for JAGS
params = c("bS","bW","bI","bSpp","b0","bA","bSN","bWN","bIN","bSppN","b0N","bAN") #parameters to monitor
inits = function() list(b0N=(1.5)) #starting values of fitted parameters
input = list(N=N,pnue=pnue.out,narea=Narea.out,group=group,site=site,woody=woody,species=species,Nspp=Nspp,away=away,woodyS=woodyS,invaderS=invaderS) #input data

#run JAGS model
jagsY <- jags(model = "model.txt",data = input,param=params,# inits=inits,
              n.chains = 3, #number of separate MCMC chains 3
              n.iter =5000, #number of iterations per chain; 10000
              n.thin=5, #thinning 5
              n.burnin = 1000) #number of initial iterations to discard 2000
jagsY 
plot(jagsY)
  #wow, fascinating effects here... graph
  #save posteriors for initial model to plot leaf N vs. Asat? (or summary version, like for PS)
    #want means and SEs of Narea and Asat (and A50 or A100) at the species level
      #what is the best way to derive these from posteriors?
detach.jags()
attach.jags(jagsY)
#compare groups: PNUE
pnue.inv.away.herb = b0 + bI[,1] + bA[,1]
pnue.inv.away.woody = b0 + bI[,2] + bA[,2] + bW
pnue.inv.home.herb = b0 + bI[,1] 
pnue.inv.home.woody = b0 + bI[,2] + bW
pnue.nat.home.herb = b0  
pnue.nat.home.woody = b0 + bW
pnue.comp = data.frame(pnue.inv.away.herb,pnue.inv.away.woody,pnue.inv.home.herb,pnue.inv.home.woody,pnue.nat.home.herb,pnue.nat.home.woody)
boxplot(pnue.comp[,c(1,3,5)],outline=F)
boxplot(pnue.comp[,c(2,4,6)],outline=F)

    #compare groups: Aarea
    narea.inv.away.herb = b0N + bIN[,1] + bAN[,1]
    narea.inv.away.woody = b0N + bIN[,2] + bAN[,2] + bWN
    narea.inv.home.herb = b0N + bIN[,1] 
    narea.inv.home.woody = b0N + bIN[,2] + bWN
    narea.nat.home.herb = b0N  
    narea.nat.home.woody = b0N + bWN
    narea.comp = data.frame(narea.inv.away.herb,narea.inv.away.woody,narea.inv.home.herb,narea.inv.home.woody,narea.nat.home.herb,narea.nat.home.woody)
    boxplot(narea.comp[,c(1,3,5)],outline=F)
    boxplot(narea.comp[,c(2,4,6)],outline=F)
    narea.awayhome.h = narea.comp[,1]-narea.comp[,3]
    narea.awaynat.h = narea.comp[,1]-narea.comp[,5]
    narea.homenat.h = narea.comp[,3]-narea.comp[,5]
    x = narea.awayhome.h; length(x[x<0])/length(x) #P=0.107
    x = narea.awaynat.h; length(x[x>0])/length(x) #P=0.292
    x = narea.homenat.h; length(x[x>0])/length(x) #P=0.108
    narea.awayhome.w = narea.comp[,2]-narea.comp[,4]
    narea.awaynat.w = narea.comp[,2]-narea.comp[,6]
    narea.homenat.w = narea.comp[,4]-narea.comp[,6]
    x = narea.awayhome.w; length(x[x<0])/length(x) #P=0.08*marginal shift higher in away range
    x = narea.awaynat.w; length(x[x>0])/length(x) #P=0.463
    x = narea.homenat.w; length(x[x>0])/length(x) #P=0.251


#same for PNUE100, without woody predictor
dfnew = cbind(dat,pnue100.out)
dfN = dfnew[dfnew$woody==1,]
Nmass = dfN$Nmass*100
Asat = dfN$Asat
Alow = dfN$A.q100 #A at 100 PPFD, woodies only
pnue100 = dfN$pnue100.out
group = as.numeric(as.factor(dfN$group)) #1=invaway 2=invhome 3=native
species = as.numeric(as.factor(dfN$Species)) #48 species
Nspp = max(species)
site = as.numeric(as.factor(dfN$region))
N = length(group)
woody = as.numeric(dfN$woody==1)
SLA = log(dfN$SLA)
aqy = dfN$alpha.emp
away = as.numeric(dfN$homeaway=="away") #1 = away, 0 = home
native = as.numeric(dfN$invader=="native") #1 = native, 0 = invader
invader = as.numeric(dfN$invader=="invader") #1 = invader, 0 = native

#covariates at the species level
sppdf = dfN[duplicated(dfN$Species)==F,]
sppdf = sppdf[order(sppdf$Species),]
speciesS = c(1:dim(sppdf)[1])
woodyS = as.numeric(sppdf$woody==1)
invaderS = as.numeric(sppdf$invader=="invader")
groupS = as.numeric(as.factor(sppdf$group)) #1=invaway 2=invhome 3=native


mod.n <- "model
{
  for(i in 1:N) {

    #PNUE
    pnue[i] ~ dnorm(mu.pnue[i],tau)
    #separate base model predictors according to whether species is native or invasive
    mu.pnue[i] <- ifelse(group[i] < 3,
        bA*away[i] + bSpp[species[i]] + bS[site[i]]
    ,
        bSpp[species[i]] + bS[site[i]]
    )

}

  #priors for fixed effect slopes (groups)
  #bW ~ dnorm(0, 0.001)
  #b0 ~ dnorm(0, 0.001)

  #for(i in 1:2) {
    bA ~ dnorm(0, 0.001)
    bI ~ dnorm(0, 0.001)
  #}

  #random effect priors
  for(i in 1:3) {
    bS[i] ~ dnorm(0,tauS1) #site random effect
  }

  #species REs
  for(i in 1:Nspp) {
    bSpp[i] ~ dnorm(mu.bSpp[i],tauSpp)
    mu.bSpp[i] <- bI*invaderS[i] #invader effect must be here, else lost in species RE
    #Spp.pred[i] <-  #can predict values here at the species level
  }
  
  #priors, variances
  tau <- sigma^-2 
  sigma ~ dunif(0, 10) 
  tauSpp <- sigmaSpp^-2 
  sigmaSpp ~ dunif(0, 10) 
  tauS1 <- sigmaS1^-2 
  sigmaS1 ~ dunif(0, 10) 

}" #end model
write(mod.n, "model.txt")

params = c("bS","bI","bSpp","bA","mu.pnue") #parameters to monitor
inits = function() list(bSpp=rnorm(Nspp),b0=rnorm(1),bA=rnorm(1),bI=rnorm(1)) #starting values of fitted parameters
input = list(N=N,pnue=pnue100,group=group,site=site,species=species,Nspp=Nspp,away=away,invaderS=invaderS) #input data

jagsY <- jags(model = "model.txt",data = input,param=params,inits=inits,
              n.chains = 3, #number of separate MCMC chains 3
              n.iter =10000, #number of iterations per chain; 10000
              n.thin=5, #thinning 5
              n.burnin = 5000) #number of initial iterations to discard 2000
jagsY 
plot(jagsY)
  #converges if overall intercept omitted
  #significant away increase in PNUE100 (~.12), slight nonsig increase for invaders over natives
attach.jags(jagsY)
plot(colMeans(mu.pnue),pnue100); abline(0,1)
summary(lm(pnue100~colMeans(mu.pnue))) #pseudo-R2=0.4 at ind level
#compare groups
b0 = rowMeans(mu.pnue)
pnue.inv.away.woody = b0 + bI + bA
pnue.inv.home.woody = b0 + bI
pnue.nat.home.woody = b0
pnue100.comp = data.frame(pnue.inv.away.woody,pnue.inv.home.woody,pnue.nat.home.woody)
boxplot(pnue100.comp,outline=F)

dat$Narea = Narea.out


###############Plots A-N showing differences in PNUE among groups
  #3 panel: Asat-N for herbs, Asat-N for woodies, A100-N for woodies
  #means and SEs per species-HA groups

Asat.inv.away.herb = tapply(dat$Asat[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==0],dat$Species[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==0],function(x)mean(x,na.rm=T))
Asat.inv.home.herb = tapply(dat$Asat[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==0],dat$Species[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==0],function(x)mean(x,na.rm=T))
Asat.nat.home.herb = tapply(dat$Asat[dat$invader=="native"&dat$homeaway=="home"&dat$woody==0],dat$Species[dat$invader=="native"&dat$homeaway=="home"&dat$woody==0],function(x)mean(x,na.rm=T))

Asat.inv.away.herb.sd = tapply(dat$Asat[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==0],dat$Species[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==0],function(x)sd(x,na.rm=T))
Asat.inv.home.herb.sd = tapply(dat$Asat[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==0],dat$Species[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==0],function(x)sd(x,na.rm=T))
Asat.nat.home.herb.sd = tapply(dat$Asat[dat$invader=="native"&dat$homeaway=="home"&dat$woody==0],dat$Species[dat$invader=="native"&dat$homeaway=="home"&dat$woody==0],function(x)sd(x,na.rm=T))

Asat.inv.away.herb.N = tapply(dat$Asat[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==0],dat$Species[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==0],function(x)length(x[!is.na(x)]))
Asat.inv.home.herb.N = tapply(dat$Asat[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==0],dat$Species[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==0],function(x)length(x[!is.na(x)]))
Asat.nat.home.herb.N = tapply(dat$Asat[dat$invader=="native"&dat$homeaway=="home"&dat$woody==0],dat$Species[dat$invader=="native"&dat$homeaway=="home"&dat$woody==0],function(x)length(x[!is.na(x)]))

Asat.inv.away.herb.se = Asat.inv.away.herb.sd/sqrt(Asat.inv.away.herb.N)
Asat.inv.home.herb.se = Asat.inv.home.herb.sd/sqrt(Asat.inv.home.herb.N)
Asat.nat.home.herb.se = Asat.nat.home.herb.sd/sqrt(Asat.nat.home.herb.N)

dat$Narea = Narea.out
Narea.inv.away.herb = tapply(dat$Narea[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==0],dat$Species[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==0],function(x)mean(x,na.rm=T))
Narea.inv.home.herb = tapply(dat$Narea[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==0],dat$Species[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==0],function(x)mean(x,na.rm=T))
Narea.nat.home.herb = tapply(dat$Narea[dat$invader=="native"&dat$homeaway=="home"&dat$woody==0],dat$Species[dat$invader=="native"&dat$homeaway=="home"&dat$woody==0],function(x)mean(x,na.rm=T))

Narea.inv.away.herb.sd = tapply(dat$Narea[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==0],dat$Species[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==0],function(x)sd(x,na.rm=T))
Narea.inv.home.herb.sd = tapply(dat$Narea[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==0],dat$Species[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==0],function(x)sd(x,na.rm=T))
Narea.nat.home.herb.sd = tapply(dat$Narea[dat$invader=="native"&dat$homeaway=="home"&dat$woody==0],dat$Species[dat$invader=="native"&dat$homeaway=="home"&dat$woody==0],function(x)sd(x,na.rm=T))

Narea.inv.away.herb.N = tapply(dat$Narea[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==0],dat$Species[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==0],function(x)length(x[!is.na(x)]))
Narea.inv.home.herb.N = tapply(dat$Narea[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==0],dat$Species[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==0],function(x)length(x[!is.na(x)]))
Narea.nat.home.herb.N = tapply(dat$Narea[dat$invader=="native"&dat$homeaway=="home"&dat$woody==0],dat$Species[dat$invader=="native"&dat$homeaway=="home"&dat$woody==0],function(x)length(x[!is.na(x)]))

Narea.inv.away.herb.se = Narea.inv.away.herb.sd/sqrt(Narea.inv.away.herb.N)
Narea.inv.home.herb.se = Narea.inv.home.herb.sd/sqrt(Narea.inv.home.herb.N)
Narea.nat.home.herb.se = Narea.nat.home.herb.sd/sqrt(Narea.nat.home.herb.N)


Asat.inv.away.woody = tapply(dat$Asat[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],function(x)mean(x,na.rm=T))
Asat.inv.home.woody = tapply(dat$Asat[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],function(x)mean(x,na.rm=T))
Asat.nat.home.woody = tapply(dat$Asat[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],function(x)mean(x,na.rm=T))

Asat.inv.away.woody.sd = tapply(dat$Asat[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],function(x)sd(x,na.rm=T))
Asat.inv.home.woody.sd = tapply(dat$Asat[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],function(x)sd(x,na.rm=T))
Asat.nat.home.woody.sd = tapply(dat$Asat[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],function(x)sd(x,na.rm=T))

Asat.inv.away.woody.N = tapply(dat$Asat[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],function(x)length(x[!is.na(x)]))
Asat.inv.home.woody.N = tapply(dat$Asat[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],function(x)length(x[!is.na(x)]))
Asat.nat.home.woody.N = tapply(dat$Asat[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],function(x)length(x[!is.na(x)]))

Asat.inv.away.woody.se = Asat.inv.away.woody.sd/sqrt(Asat.inv.away.woody.N)
Asat.inv.home.woody.se = Asat.inv.home.woody.sd/sqrt(Asat.inv.home.woody.N)
Asat.nat.home.woody.se = Asat.nat.home.woody.sd/sqrt(Asat.nat.home.woody.N)

Narea.inv.away.woody = tapply(dat$Narea[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],function(x)mean(x,na.rm=T))
Narea.inv.home.woody = tapply(dat$Narea[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],function(x)mean(x,na.rm=T))
Narea.nat.home.woody = tapply(dat$Narea[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],function(x)mean(x,na.rm=T))

Narea.inv.away.woody.sd = tapply(dat$Narea[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],function(x)sd(x,na.rm=T))
Narea.inv.home.woody.sd = tapply(dat$Narea[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],function(x)sd(x,na.rm=T))
Narea.nat.home.woody.sd = tapply(dat$Narea[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],function(x)sd(x,na.rm=T))

Narea.inv.away.woody.N = tapply(dat$Narea[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],function(x)length(x[!is.na(x)]))
Narea.inv.home.woody.N = tapply(dat$Narea[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],function(x)length(x[!is.na(x)]))
Narea.nat.home.woody.N = tapply(dat$Narea[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],function(x)length(x[!is.na(x)]))

Narea.inv.away.woody.se = Narea.inv.away.woody.sd/sqrt(Narea.inv.away.woody.N)
Narea.inv.home.woody.se = Narea.inv.home.woody.sd/sqrt(Narea.inv.home.woody.N)
Narea.nat.home.woody.se = Narea.nat.home.woody.sd/sqrt(Narea.nat.home.woody.N)


A100.inv.away.woody = tapply(dat$A.q100[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],function(x)mean(x,na.rm=T))
A100.inv.home.woody = tapply(dat$A.q100[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],function(x)mean(x,na.rm=T))
A100.nat.home.woody = tapply(dat$A.q100[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],function(x)mean(x,na.rm=T))

A100.inv.away.woody.sd = tapply(dat$A.q100[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],function(x)sd(x,na.rm=T))
A100.inv.home.woody.sd = tapply(dat$A.q100[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],function(x)sd(x,na.rm=T))
A100.nat.home.woody.sd = tapply(dat$A.q100[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],function(x)sd(x,na.rm=T))

A100.inv.away.woody.N = tapply(dat$A.q100[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="away"&dat$woody==1],function(x)length(x[!is.na(x)]))
A100.inv.home.woody.N = tapply(dat$A.q100[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="invader"&dat$homeaway=="home"&dat$woody==1],function(x)length(x[!is.na(x)]))
A100.nat.home.woody.N = tapply(dat$A.q100[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],dat$Species[dat$invader=="native"&dat$homeaway=="home"&dat$woody==1],function(x)length(x[!is.na(x)]))

A100.inv.away.woody.se = A100.inv.away.woody.sd/sqrt(A100.inv.away.woody.N)
A100.inv.home.woody.se = A100.inv.home.woody.sd/sqrt(A100.inv.home.woody.N)
A100.nat.home.woody.se = A100.nat.home.woody.sd/sqrt(A100.nat.home.woody.N)


#PNUE credible intervals: for graphing
pnue.awayhome.h = pnue.comp[,1]-pnue.comp[,3]
pnue.awaynat.h = pnue.comp[,1]-pnue.comp[,5]
pnue.homenat.h = pnue.comp[,3]-pnue.comp[,5]
pnue.awayhome.w = pnue.comp[,2]-pnue.comp[,4]
pnue.awaynat.w = pnue.comp[,2]-pnue.comp[,6]
pnue.homenat.w = pnue.comp[,4]-pnue.comp[,6]
pnue100.awayhome.w = pnue100.comp[,1]-pnue100.comp[,2]
pnue100.awaynat.w = pnue100.comp[,1]-pnue100.comp[,3]
pnue100.homenat.w = pnue100.comp[,2]-pnue100.comp[,3]
x = pnue.awayhome.h; length(x[x<0])/length(x) #P=0.163, marginal
x = pnue.awaynat.h; length(x[x<0])/length(x) #P=0.024*
x = pnue.homenat.h; length(x[x<0])/length(x) #P=0.074, marginal
x = pnue.awayhome.w; length(x[x<0])/length(x) #P=0.21 NS
x = pnue.awaynat.w; length(x[x<0])/length(x) #P=0.53 NS
x = pnue.homenat.w; length(x[x<0])/length(x) #P=0.68 NS
x = pnue100.awayhome.w; length(x[x<0])/length(x) #P=0.012*
x = pnue100.awaynat.w; length(x[x<0])/length(x) #P=0.095 marginal
x = pnue100.homenat.w; length(x[x<0])/length(x) #P=0.78 NS



##############graph results

library(RColorBrewer)
nat.col = brewer.pal(3,"Dark2")[1]
home.col = brewer.pal(3,"Dark2")[2]
away.col = brewer.pal(3,"Dark2")[3]

#pdf("PNUE_fig.pdf",height=6,width=10)
par(mfrow=c(2,3),oma=c(0,0,2,0),mar=c(5,6,1,1))
plot(Narea.inv.away.herb,Asat.inv.away.herb,ylim=c(8,40),xlim=c(0.4,1.7),col=away.col,pch=19,cex=1.5,
     ylab = expression(A[sat]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),
     xlab = expression(N[area]*" (g N "*m^-2*")"),
     cex.lab=1.6,cex.axis=1.4)
arrows(Narea.inv.away.herb,Asat.inv.away.herb-Asat.inv.away.herb.se,Narea.inv.away.herb,Asat.inv.away.herb+Asat.inv.away.herb.se,code=3,length=.05,angle=90,col=away.col,lwd=.5)
arrows(Narea.inv.away.herb+Narea.inv.away.herb.se,Asat.inv.away.herb,Narea.inv.away.herb-Narea.inv.away.herb.se,Asat.inv.away.herb,code=3,length=.05,angle=90,col=away.col,lwd=.5)
abline(lsfit(Narea.inv.away.herb,Asat.inv.away.herb),col=away.col,lwd=2)
points(Narea.inv.home.herb,Asat.inv.home.herb,ylim=c(0,40),xlim=c(.2,1.2),col=home.col,pch=19,cex=1.5)
arrows(Narea.inv.home.herb,Asat.inv.home.herb-Asat.inv.home.herb.se,Narea.inv.home.herb,Asat.inv.home.herb+Asat.inv.home.herb.se,code=3,length=.05,angle=90,col=home.col,lwd=.5)
arrows(Narea.inv.home.herb+Narea.inv.home.herb.se,Asat.inv.home.herb,Narea.inv.home.herb-Narea.inv.home.herb.se,Asat.inv.home.herb,code=3,length=.1,angle=90,col=home.col,lwd=.5)
abline(lsfit(Narea.inv.home.herb,Asat.inv.home.herb),col=home.col,lwd=2)
points(Narea.nat.home.herb,Asat.nat.home.herb,ylim=c(0,40),xlim=c(.2,1.2),col=nat.col,pch=19,cex=1.5)
arrows(Narea.nat.home.herb,Asat.nat.home.herb-Asat.nat.home.herb.se,Narea.nat.home.herb,Asat.nat.home.herb+Asat.nat.home.herb.se,code=3,length=.05,angle=90,col=nat.col,lwd=.5)
arrows(Narea.nat.home.herb+Narea.nat.home.herb.se,Asat.nat.home.herb,Narea.nat.home.herb-Narea.nat.home.herb.se,Asat.nat.home.herb,code=3,length=.05,angle=90,col=nat.col,lwd=.5)
abline(lsfit(Narea.nat.home.herb,Asat.nat.home.herb),col=nat.col,lwd=2)
mtext("A)",side=3,line=-1,cex=2.5,xpd=NA,adj=0,outer=T)
mtext(expression(Fields - A[sat]),side=3,line=-2,adj=.05,cex=1.2)

plot(Narea.inv.away.woody,Asat.inv.away.woody,ylim=c(3,25),xlim=c(.3,1.4),col=away.col,pch=19,cex=1.5,
     ylab = expression(A[sat]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),
     xlab = expression(N[area]*" (g N "*m^-2*")"),
     cex.lab=1.6,cex.axis=1.4)
arrows(Narea.inv.away.woody,Asat.inv.away.woody-Asat.inv.away.woody.se,Narea.inv.away.woody,Asat.inv.away.woody+Asat.inv.away.woody.se,code=3,length=.05,angle=90,col=away.col,lwd=.5)
arrows(Narea.inv.away.woody+Narea.inv.away.woody.se,Asat.inv.away.woody,Narea.inv.away.woody-Narea.inv.away.woody.se,Asat.inv.away.woody,code=3,length=.05,angle=90,col=away.col,lwd=.5)
abline(lsfit(Narea.inv.away.woody,Asat.inv.away.woody),col=away.col,lwd=2)
points(Narea.inv.home.woody,Asat.inv.home.woody,ylim=c(0,40),xlim=c(.2,1.2),col=home.col,pch=19,cex=1.5)
arrows(Narea.inv.home.woody,Asat.inv.home.woody-Asat.inv.home.woody.se,Narea.inv.home.woody,Asat.inv.home.woody+Asat.inv.home.woody.se,code=3,length=.05,angle=90,col=home.col,lwd=.5)
arrows(Narea.inv.home.woody+Narea.inv.home.woody.se,Asat.inv.home.woody,Narea.inv.home.woody-Narea.inv.home.woody.se,Asat.inv.home.woody,code=3,length=.05,angle=90,col=home.col,lwd=.5)
abline(lsfit(Narea.inv.home.woody,Asat.inv.home.woody),col=home.col,lwd=2)
points(Narea.nat.home.woody,Asat.nat.home.woody,ylim=c(0,40),xlim=c(.2,1.2),col=nat.col,pch=19,cex=1.5)
arrows(Narea.nat.home.woody,Asat.nat.home.woody-Asat.nat.home.woody.se,Narea.nat.home.woody,Asat.nat.home.woody+Asat.nat.home.woody.se,code=3,length=.05,angle=90,col=nat.col,lwd=.5)
arrows(Narea.nat.home.woody+Narea.nat.home.woody.se,Asat.nat.home.woody,Narea.nat.home.woody-Narea.nat.home.woody.se,Asat.nat.home.woody,code=3,length=.05,angle=90,col=nat.col,lwd=.5)
abline(lsfit(Narea.nat.home.woody,Asat.nat.home.woody),col=nat.col,lwd=2)
mtext(expression(Forests - A[sat]),side=3,line=-2,adj=.05,cex=1.2)

plot(Narea.inv.away.woody,A100.inv.away.woody,ylim=c(1,4.5),xlim=c(.3,1.4),col=away.col,pch=19,cex=1.5,
     ylab = expression(A[100]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),
     xlab = expression(N[area]*" (g N "*m^-2*")"),
     cex.lab=1.6,cex.axis=1.4)
arrows(Narea.inv.away.woody,A100.inv.away.woody-A100.inv.away.woody.se,Narea.inv.away.woody,A100.inv.away.woody+A100.inv.away.woody.se,code=3,length=.1,angle=90,col=away.col,lwd=.5)
arrows(Narea.inv.away.woody+Narea.inv.away.woody.se,A100.inv.away.woody,Narea.inv.away.woody-Narea.inv.away.woody.se,A100.inv.away.woody,code=3,length=.1,angle=90,col=away.col,lwd=.5)
abline(lsfit(Narea.inv.away.woody,A100.inv.away.woody),col=away.col,lwd=2)
points(Narea.inv.home.woody,A100.inv.home.woody,ylim=c(0,40),xlim=c(.2,1.2),col=home.col,pch=19,cex=1.5)
arrows(Narea.inv.home.woody,A100.inv.home.woody-A100.inv.home.woody.se,Narea.inv.home.woody,A100.inv.home.woody+A100.inv.home.woody.se,code=3,length=.1,angle=90,col=home.col,lwd=.5)
arrows(Narea.inv.home.woody+Narea.inv.home.woody.se,A100.inv.home.woody,Narea.inv.home.woody-Narea.inv.home.woody.se,A100.inv.home.woody,code=3,length=.1,angle=90,col=home.col,lwd=.5)
abline(lsfit(Narea.inv.home.woody,A100.inv.home.woody),col=home.col,lwd=2)
points(Narea.nat.home.woody,A100.nat.home.woody,ylim=c(0,40),xlim=c(.2,1.2),col=nat.col,pch=19,cex=1.5)
arrows(Narea.nat.home.woody,A100.nat.home.woody-A100.nat.home.woody.se,Narea.nat.home.woody,A100.nat.home.woody+A100.nat.home.woody.se,code=3,length=.1,angle=90,col=nat.col,lwd=.5)
arrows(Narea.nat.home.woody+Narea.nat.home.woody.se,A100.nat.home.woody,Narea.nat.home.woody-Narea.nat.home.woody.se,A100.nat.home.woody,code=3,length=.1,angle=90,col=nat.col,lwd=.5)
abline(lsfit(Narea.nat.home.woody,A100.nat.home.woody),col=nat.col,lwd=2)
mtext(expression(Forests - A[100]),side=3,line=-2,adj=.05,cex=1.2)

par(mar=c(6,8,1,2),new=F)
box1 = boxplot(pnue.comp[,c(1,3,5)],outline=F,col=c(away.col,home.col,nat.col),
        xlab="",pch=19,cex=2,cex.axis=1.4,xaxt="n",cex.lab=1.5,ylim=c(10,38),
        ylab = expression(PNUE[sat]*" ("*mu*"mol "*CO[2]*" "*g^-1*" N "*s^-1*")"))
mtext("B)",side=3,line=.5,cex=2.5,xpd=NA,adj=-.65,outer=F)
mtext(side=1,c("Invaders","Invaders","Native"),at=c(1,2,3),col=c(away.col,home.col,nat.col),line=.7,cex=.8)
mtext(side=1,c("(away)","(home)",""),at=c(1,2,3),col=c(away.col,home.col,nat.col),line=2.1,cex=.8)
mtext(expression(Fields - A[sat]),side=3,line=-2,adj=.05,cex=1.2)
text(1:3,box1$stats[5,]+2,c("a","a","b"),col=c(away.col,home.col,nat.col),cex=1.7)  
box1 = boxplot(pnue.comp[,c(2,4,6)],outline=F,col=c(away.col,home.col,nat.col),
        xlab="",pch=19,cex=2,cex.axis=1.3,xaxt="n",cex.lab=1.5,ylim=c(7,30),
        ylab = expression(PNUE[sat]*" ("*mu*"mol "*CO[2]*" "*g^-1*" N "*s^-1*")"))
mtext(side=1,c("Invaders","Invaders","Native"),at=c(1,2,3),col=c(away.col,home.col,nat.col),line=.7,cex=.8)
mtext(side=1,c("(away)","(home)",""),at=c(1,2,3),col=c(away.col,home.col,nat.col),line=2.1,cex=.8)
mtext(expression(Forests - A[sat]),side=3,line=-2,adj=.05,cex=1.2)
text(1:3,box1$stats[5,]+2,c("a","a","a"),col=c(away.col,home.col,nat.col),cex=1.7)  
box1 = boxplot(pnue100.comp,outline=F,col=c(away.col,home.col,nat.col),
        xlab="",pch=19,cex=2,cex.axis=1.3,xaxt="n",cex.lab=1.5,ylim=c(2.5,5.5),
        ylab = expression(PNUE[100]*" ("*mu*"mol "*CO[2]*" "*g^-1*" N "*s^-1*")"))
mtext(side=1,c("Invaders","Invaders","Native"),at=c(1,2,3),col=c(away.col,home.col,nat.col),line=.7,cex=.8)
mtext(side=1,c("(away)","(home)",""),at=c(1,2,3),col=c(away.col,home.col,nat.col),line=2.1,cex=.8)
mtext(expression(Forests - A[100]),side=3,line=-2,adj=.05,cex=1.2)
text(1:3,box1$stats[5,]+.2,c("a","b","ab"),col=c(away.col,home.col,nat.col),cex=1.7)  
#dev.off()


########N allocation model

###data

#pull in leaf temps for bioenergetics component
temp.dat = read.csv("C:\\Users\\fridley\\OneDrive - Clemson University\\academic\\projects\\IOS_FranceJapan\\leaftempsbyID.csv") #file of median leaf temps per sample
names(temp.dat) = c("ID","Tleaf")
dat2 = merge(dat,temp.dat,by.x="ID",by.y="ID",all.x=T) #merge leaf temps with dataset
dat = dat2

dat=dat[!is.na(dat$Nmass),] #only samples with known Nmass

N = dim(dat)[1]
SLA = dat$SLA2 #mean estimate
SLA.se = dat$SLA2.se #se of SLA estimates (those that were originally missing)
Nmass = dat$Nmass   #36 missing
Rubisco.ug = as.numeric(dat$Rubisco.content) #39 missing
WallpercN = as.numeric(dat$Cell.wall..N)/100    #g N per g cell wall, 45 missing
Wallmass = (as.numeric(dat$Cell.wall.mass) / as.numeric(dat$cell.wall.leaf.area.cm2))/1000 #cell wall mass, mg per cm2, 37 missing
DefenseN = dat$total_alkaloidN_mg_g/1000 + dat$cyanN_mg_g/1000 #g N in defense per g leaf, 23 missing
  dat$alk = dat$total_alkaloids_mg_g #mg atropine equivalents per g dry leaf mass
  dat$cyan = dat$cyan_mg_g #mg KCN equivalents per g dry leaf mass
  dat$alkN = dat$alk*.04 #4% by mass in N
  dat$cyanN = dat$cyan*.21 #21%by mass in N
  dat$defNmg = dat$alkN + dat$cyanN #mg N per g dry leaf mass
  defNperc = (dat$defNmg/1000) / Nmass
DefenseNmass = DefenseN*Nmass  
Chl.ugcm2 = as.numeric(as.character(dat$Total.Chlorophyll)) #ug chl per cm2
Chl.mmolm2 = (Chl.ugcm2 / 1000000) * 1000 * 10000 / 893.5 #mmol chl per m2 leaf, 14 missing
Tleaf = dat$Tleaf #27 missing
#calculate temp-based specific activity of electron transport, per Niinemets & Tenhunen (1997)
Jmc.temp = function(c,dHa,dHd,dS,R,Tk) {
  exp(c-(dHa/(R*Tk))) / (1 + exp(  ((dS*Tk)-dHd) / (R*Tk)) )
}
Jmc = Jmc.temp(c=14.77, dHa=24100, dHd=564150, dS=1810, R=8.314, Tk=Tleaf+273.15)
Jmax = dat$Jmax #28 missing
#Not used here, but use when relating Vcmax to Rubisco content
#calculate temp-based specific activity of Rubisco, per Niinemets & Tenhunen (1997)
  #same equation as Jmc, different constants
Vcr = Jmc.temp(c=32.9, dHa=74000, dHd=203000, dS=645, R=8.314, Tk=Tleaf+273.15)
Vcmax = dat$Vcmax #28 missing


#model N allocation 
mod.allocate <- "model
{
    #missing values
    for(i in 1:N) {
      #sla[i] ~ dnorm(sla[i],1/sla.se[i]) #sla is mostly known, esimated from prior posterior SEs
      rub[i] ~ dnorm(mean.rub,tau.rub) #fill Rubisco (area-based) missing values
      #nmass[i] ~ dnorm(0.02,tau.nmass) #fill Nmass missing values
      wallpercN[i] ~ dnorm(mean.wallN,tau.wallN) #fill wallN missing values (this is for percents)
      wallmass[i] ~ dnorm(mean.wallmass,tau.wallmass) #fill in missing wall mass values
      defenseN[i] ~ dnorm(mean.defN,tau.defN) #fill in missing defense N values
      defenseNmass[i] ~ dnorm(mean.defNm,tau.defNm) #fill in missing defense N values
      chl.mmolm2[i] ~ dnorm(mean.chl,tau.chl) #fill in missing chl values
      jmc[i] ~ dnorm(mean.jmc,tau.jmc)
      jmax[i] ~ dnorm(mean.jmax,tau.jmax)
      #vcr[i] ~ dnorm(mean.vcr,tau.vcr)
      #vcmax[i] ~ dnorm(mean.vcmax,tau.vcmax)
    }
    
    tau.rub ~ dgamma(0.01,0.01)
    #tau.nmass ~ dgamma(0.01,0.01)
    tau.wallmass ~ dgamma(0.01,0.01)
    tau.wallN ~ dgamma(0.01,0.01)
    tau.defN ~ dgamma(0.01,0.01)
    tau.defNm ~ dgamma(0.01,0.01)
    tau.chl ~ dgamma(0.01,0.01)
    tau.jmc ~ dgamma(0.01,0.01)
    tau.jmax ~ dgamma(0.01,0.01)
    tau.rubN ~ dgamma(0.01,0.01) #for latent variable rubN, below
    #tau.vcr ~ dgamma(0.01,0.01)
    #tau.vcmax ~ dgamma(0.01,0.01)
    #tau.resid ~ dgamma(0.01,0.01)

    mean.rub ~ dnorm(212,136^-2)
    mean.nmass ~ dnorm(0,0.001)
    mean.wallmass ~ dnorm(0,0.001)
    mean.wallN ~ dnorm(0.01,0.001)
    mean.defN ~ dnorm(0.01,0.001)
    mean.defNm ~ dnorm(0.01,0.001)
    mean.chl ~ dnorm(0.01,0.001)
    mean.jmc ~ dnorm(176,0.01)
    mean.jmax ~ dnorm(0.01,0.01)
    #mean.vcr ~ dnorm(0.01,0.001)
    #mean.vcmax ~ dnorm(0.01,0.001)

  for(i in 1:N) {

    #Rubisco N per unit mass
    #chemical measurement of Rubisco
      rubN[i] <- ( rub[i] * sla[i] * .16 ) / (1000000) #g N in Rubisco per g leaf
      rubNperc[i] <- rubN[i] / nmass[i]   #percent of leaf N that is in Rubisco
    #inference from Vcmax (here ignored)
      #rubNpred[i] <- vcmax[i] / (6.25*vcr[i]*(10000/sla[i]))
      #rubNperc1[i] <- rubN[i] / nmass[i]
      #rubNperc2[i] <- rubNpred[i] / nmass[i]

    #N in light harvesting, from chl content
    chl.mmolg[i] <- ( chl.mmolm2[i] * sla[i]) / 10000 #mmol chl per g leaf
    lightN[i] <- (chl.mmolg[i] * 37.1 / 1000) * 14 #g N in light harvesting per g leaf, Takashima et al. (2004)
    lightNperc[i] <- lightN[i] / nmass[i]
    
    #N in bioenergetics
    #use Niinemets & Tenhunen (1997) approach that includes temp dependence
    beN[i] <- jmax[i] * (sla[i]/10000) * (1/(8.06*jmc[i]))
    beNperc[i] <- jmax[i] * (sla[i]/10000) * (1/nmass[i]) * (1/(8.06*jmc[i]))

    #N in cell wall
    wallN[i] <- (wallmass[i] * sla[i] * wallpercN[i]) #g N in the cell wall per g leaf
    wallNperc[i] <- wallN[i] / nmass[i] #percent of leaf N that is in the cell wall
    
    #N in defensive chemistry
    defNperc[i] <- defenseN[i] #/ nmass[i] #percent of leaf N that is in alkaloids and cyanogenic glycosides
    defN[i] <- defenseNmass[i]
    
    #Residual N: left over, must all equal Nmass
    resN[i] ~ dnorm(0,.01) 
    nmass[i] ~ sum(rubN[i],lightN[i],beN[i],wallN[i],defN[i],resN[i])

}

}" #end model
write(mod.allocate, "model.txt")

#input lists for JAGS
params = c("rubNperc","wallNperc","defNperc","lightNperc","beNperc","resN") #parameters to monitor
inits = function() list(mean.nmass=.15,mean.jmc=rnorm(1,176,10)) #starting values of fitted parameters
input = list(N=N,sla=SLA,nmass=Nmass,rub=Rubisco.ug,wallmass=Wallmass,wallpercN=WallpercN,defenseN=DefenseN,defenseNmass=DefenseNmass,chl.mmolm2=Chl.mmolm2,jmc=Jmc,jmax=Jmax,vcr=Vcr,vcmax=Vcmax) #input data

#run JAGS model
jagsY <- jags(model = "model.txt",data = input,param=params, inits=inits,
              n.chains = 3, #number of separate MCMC chains 3
              n.iter =5000, #number of iterations per chain; 10000
              n.thin=5, #thinning 5
              n.burnin = 1000) #number of initial iterations to discard 2000
jagsY 
plot(jagsY)
detach.jags()
attach.jags(jagsY)
dat$rubNperc = colMeans(rubNperc)
dat$lightNperc = colMeans(lightNperc)
dat$beNperc = colMeans(beNperc)
dat$wallNperc = colMeans(wallNperc)
dat$defNperc = colMeans(defNperc)
dat$resNperc = colMeans(resN)

mean(dat$rubNperc)
mean(dat$lightNperc)
mean(dat$beNperc)
mean(dat$wallNperc)
mean(dat$defNperc)

boxplot(dat$rubNperc~dat$woody)
boxplot(dat$lightNperc~dat$woody)
boxplot(dat$beNperc~dat$woody)
boxplot(log(dat$wallNperc)~dat$woody)
boxplot(dat$defNperc~dat$woody)

#some neg residuals are due to high Rubisco N
plot(dat$rubNperc,dat$resNperc); abline(h=0,lwd=2)

#investigate relationships
plot(dat$rubNperc,dat$lightNperc,col=dat$woody+1); abline(0,1)
plot(dat$lightNperc,dat$beNperc,col=dat$woody+1); abline(0,1)
plot(dat$rubNperc+dat$lightNperc+dat$beNperc,dat$wallNperc+dat$defNperc,col=dat$woody+1,log="")
  #some evidence of negative relationship between photo and non-photo N?
  #but note these must add up to 1 (kind of)

#relationship of Vcmax to rubisco N: not bad
plot(dat$rubNperc*dat$Nmass,dat$Vcmax*dat$SLA2/10000,col=dat$woody+1)



###################################################
#######now run process model for fixed effects on N fractions

#use allocation or total N in each component? (divide by Nmass or no?)
  #if the former, need to decide on Narea or Nmass
#Y = log(dat$rubNperc) #note lognormal
group = as.numeric(as.factor(dat$group)) #1 = invaway 2 = invhome 3 = native
woody = as.numeric(dat$woody==1) #1 = woody, 0 = herbaceous
away = as.numeric(dat$homeaway=="away") #1 = away
species = as.numeric(as.factor(dat$Species))
invader = as.numeric(dat$invader=="invader") #1 = invader, 0 = native
site = as.numeric(as.factor(dat$region))
N = dim(dat)[1]

#covariates at the species level
sppdf = dat[duplicated(dat$Species)==F,]
sppdf = sppdf[order(sppdf$Species),]
speciesS = c(1:dim(sppdf)[1])
woodyS = as.numeric(sppdf$woody==1)
invaderS = as.numeric(sppdf$invader=="invader")
groupS = as.numeric(as.factor(sppdf$group)) #1=invaway 2=invhome 3=native
Nspp = dim(sppdf)[1]

#dependent variables
Ym = cbind(dat$rubNperc,dat$lightNperc,dat$beNperc,dat$wallNperc,dat$defNperc)
#Ym = Ym[,c(1:2)]
Zscale = function(x) (x - mean(x))/sd(x)
Ym = apply(Ym,2,Zscale)
K = dim(Ym)[2] #number of Y variables
R <- diag(K) # Identity matrix as a simple prior


mod.n <- "model
{

  for(i in 1:N) { #observations

    Y[i,1:K] ~ dmnorm(mu.Y[i,],Omega[,]) #multivariate normal

    mu.Y[i,1:K] <- ifelse(group[i] < 3,
          bW[1:K]*woody[i] + bA[woody[i]+1,1:K]*away[i] + bSpp[species[i],1:K] + bS[site[i],1:K] #+ b0[1:K] 
    ,
          bW[1:K]*woody[i] + bSpp[species[i],1:K] + bS[site[i],1:K] #+ b0[1:K]
    )
  }

  #Expected value of mu.Y on original (non-Z) 0-1 scale, ignoring species and site effects
  #  Yexp[i,1:K] <- ifelse(group[i] < 3,
  #        bW[1:K]*woody[i] + bA[woody[i]+1,1:K]*away[i]
  #      ,
  #        bW[1:K]*woody[i] + bI[woodyS[s]+1,j]*invader[i]
  #  )

  #Prior for the precision matrix (inverse of the covariance matrix)
  Omega[1:K, 1:K] ~ dwish(R[,], K + 1) # Prior on the precision matrix
  Sigma[1:K, 1:K] <- inverse(Omega[,]) # Covariance matrix is the inverse of the precision matrix

  for(j in 1:K) { #loop over dependent variables

    #priors for fixed effect slopes (groups)
    bW[j] ~ dnorm(0, 0.001)
    #b0[j] ~ dnorm(0, 0.001)

    for(i in 1:2) {
      bA[i,j] ~ dnorm(0, 0.001)
      bI[i,j] ~ dnorm(0, 0.001)
    }

    #site REs
    for(i in 1:3) {
      bS[i,j] ~ dnorm(0,tauS1[j])
    }

    #species REs
    for(s in 1:Nspp) {
      bSpp[s,j] ~ dnorm(mu.bSpp[s,j],tauSpp[j])
      mu.bSpp[s,j] <- bI[woodyS[s]+1,j]*invaderS[s] #invader effect must be here, else lost in species RE
    }
  
    #priors, RE variances
    #tau[j] <- sigma[j]^-2 
    #sigma[j] ~ dunif(0, 100) 
    tauSpp[j] <- sigmaSpp[j]^-2 
    sigmaSpp[j] ~ dunif(0, 100) 
    tauS1[j] <- sigmaS1[j]^-2 
    sigmaS1[j] ~ dunif(0, 100) 
}

}" #end model
write(mod.n, "model.txt")

#input lists for JAGS
params = c("bS","bW","bI","bSpp","bA","Sigma") #parameters to monitor
inits = function() list(b0=colMeans(Ym)) #starting values of fitted parameters
input = list(R=R,K=K,N=N,Y=Ym,group=group,site=site,woody=woody,species=species,Nspp=Nspp,away=away,woodyS=woodyS,invaderS=invaderS) #input data

#run JAGS model
jagsN <- jags(model = "model.txt",data = input,param=params,#inits=inits,
              n.chains = 3, #number of separate MCMC chains 3
              n.iter =2000, #number of iterations per chain; 5000
              n.thin=5, #thinning 5
              n.burnin = 500) #number of initial iterations to discard 1000
plot(jagsN)

#examine results
attach.jags(jagsN)
cov.mat = apply(Sigma,c(2,3),mean)
round(cov.mat,3) #fitted covariance matrix
  #largest variances are Nrub and Nwall

#means, 80, and 95 CIs
bW.means = apply(bW,2,mean)
bW.80 = apply(bW,2,function(x)quantile(x,c(.1,.9)))
bW.95 = apply(bW,2,function(x)quantile(x,c(.025,.975)))
bA.means = apply(bA,c(2,3),mean)
bA.80 = apply(bA,c(2,3),function(x)quantile(x,c(.1,.9)))
bA.95 = apply(bA,c(2,3),function(x)quantile(x,c(.025,.975)))
bI.means = apply(bI,c(2,3),mean)
bI.80 = apply(bI,c(2,3),function(x)quantile(x,c(.1,.9)))
bI.95 = apply(bI,c(2,3),function(x)quantile(x,c(.025,.975)))

#stack barplot of each group
library(doBy)
b = data.frame(NR=dat$rubNperc*dat$Nmass,NL=dat$lightNperc*dat$Nmass,NB=dat$beNperc*dat$Nmass,NC=dat$wallNperc*dat$Nmass,ND=dat$defNperc*dat$Nmass,group=dat$group,woody=dat$woody)
b$Nres = dat$Nmass - rowSums(b[,1:5])
bd = summaryBy(NR+NL+NB+NC+ND+Nres ~ group + woody, FUN=mean, data=b)
bdH = bd[bd$woody==0,-c(1:2)]; rownames(bdH) = c("invway","invhome","native")
bdW = bd[bd$woody==1,-c(1:2)]; rownames(bdW) = c("invway","invhome","native")

#error bars for total N
nmass.H.means = tapply(dat$Nmass[dat$woody==0]*1000,dat$group[dat$woody==0],mean)
nmass.H.sds = tapply(dat$Nmass[dat$woody==0]*1000,dat$group[dat$woody==0],sd)
nmass.H.N = tapply(dat$Nmass[dat$woody==0]*1000,dat$group[dat$woody==0],length)
nmass.H.ses = nmass.H.sds/sqrt(nmass.H.N)
nmass.W.means = tapply(dat$Nmass[dat$woody==1]*1000,dat$group[dat$woody==1],mean)
nmass.W.sds = tapply(dat$Nmass[dat$woody==1]*1000,dat$group[dat$woody==1],sd)
nmass.W.N = tapply(dat$Nmass[dat$woody==1]*1000,dat$group[dat$woody==1],length)
nmass.W.ses = nmass.H.sds/sqrt(nmass.H.N)

library(RColorBrewer)
bar.cols = brewer.pal(n = 8, name = 'Greens')[c(8,7,5,1)]
bar.cols = c(bar.cols[1:3],"burlywood4","black","darkgray")
library(pBrackets)

#pdf("Nallocation_stackedbars.pdf",height=8,width=12)

par(mar=c(4,4.5,1,7),mfrow=c(1,1),oma=c(0,0,1,1),mfrow=c(1,2))
b1 = barplot(t(as.matrix(bdH*1000)),ylab="mg N per N leaf",col=bar.cols,ylim=c(0,28),names.arg=c("Invaders","Invaders","Natives"),
             cex.lab=1.5,cex.axis=1.4,cex.names=1.5)
arrows(b1,nmass.H.means,b1,nmass.H.means+nmass.H.ses,angle=90,length=.1,lwd=3)
mtext("Rubisco***",side=4,las=1,at=4,cex=1.3,col=bar.cols[1],line=-.5)
mtext("Light harvest*",side=4,las=1,at=8.3,cex=1.3,col=bar.cols[2],line=-.5)
mtext("Bioenergetics*",side=4,las=1,at=10.7,cex=1.3,col=bar.cols[3],line=-.5)
mtext("Cell wall",side=4,las=1,at=12.8,cex=1.3,col=bar.cols[4],line=-.5)
mtext("Defense**",side=4,las=1,at=14.5,cex=1.3,col=bar.cols[5],line=-.5)
mtext("Residual",side=4,las=1,at=20,cex=1.3,col="gray45",line=-.5)
mtext("Total",side=4,las=1,at=24.5,cex=1.3,col="black",line=-.5)
mtext("A) Fields",side=3,cex=2.5,line=-1,at=.5)
mtext(side=1,c("(away)","(home)"),at=b1[1:2],line=2.3,cex=1.5)

text(b1-.3,bdH[,1]*1000-.7,c("a","b","b"),col="white",cex=1.5)
text(b1-.3,rowSums(bdH[,1:2])*1000-.7,c("a","b","b"),col="white",cex=1.5)
text(b1-.3,rowSums(bdH[,1:3])*1000-.7,c("ab","a","b"),col="white",cex=1.5)
text(b1-.3,rowSums(bdH[,1:4])*1000-.7,c("a","a","a"),col="white",cex=1.5)
text(b1-.3,rowSums(bdH[,1:5])*1000+.7,c("a","b","ab"),col="white",cex=1.5)
#brackets(4.55,11.5,4.55,3.4,xpd=NA,col=bar.cols[1],curvature=.3,h=.2,lwd=3)
brackets(5.1,11.5,5.1,3.4,xpd=NA,col=bar.cols[1],curvature=.3,h=.2,lwd=3)
text(5.4,8,"**",xpd=NA,cex=1.5,col=bar.cols[1])

b1 = barplot(t(as.matrix(bdW*1000)),ylab="",col=bar.cols,ylim=c(0,28),names.arg=c("Invaders","Invaders","Natives"),
             cex.lab=1.5,cex.axis=1.4,cex.names=1.5)
arrows(b1,nmass.W.means,b1,nmass.W.means+nmass.W.ses,angle=90,length=.1,lwd=3)
mtext("Rubisco*",side=4,las=1,at=5,cex=1.3,col=bar.cols[1],line=-.5)
mtext("Light harvest*",side=4,las=1,at=10.5,cex=1.3,col=bar.cols[2],line=-.5)
mtext("Bioenergetics*",side=4,las=1,at=13.2,cex=1.3,col=bar.cols[3],line=-.5)
mtext("Cell wall",side=4,las=1,at=17,cex=1.3,col=bar.cols[4],line=-.5)
mtext("Defense",side=4,las=1,at=20,cex=1.3,col=bar.cols[5],line=-.5)
mtext("Residual",side=4,las=1,at=22,cex=1.3,col="gray45",line=-.5)
mtext("Total",side=4,las=1,at=23.9,cex=1.3,col="black",line=-.5)
mtext("B) Forests",side=3,cex=2.5,line=-1,at=.5)
mtext(side=1,c("(away)","(home)"),at=b1[1:2],line=2.3,cex=1.5)

text(b1-.3,bdW[,1]*1000-.7,c("a","b","ab"),col="white",cex=1.5)
text(b1-.3,rowSums(bdW[,1:2])*1000-.7,c("a","ab","b"),col="white",cex=1.5)
text(b1-.3,rowSums(bdW[,1:3])*1000-.7,c("a","a","b"),col="white",cex=1.5)
text(b1-.3,rowSums(bdW[,1:4])*1000-.7,c("a","a","a"),col="white",cex=1.5)
text(b1-.3,rowSums(bdW[,1:5])*1000+.7,c("a","a","a"),col="white",cex=1.5)

#dev.off()


#plotting coefficient posteriors (not for pub?)

wcol = "darkgreen"
hcol = "goldenrod2"
off = .1 #offset

#pdf("Nallocation_effects.pdf",height=8,width=7)

par(mfrow=c(3,1),mar=c(.5,5,.5,0),oma=c(.5,.5,6,5))

plot(1:5-off,bA.means[1,],ylim=c(-.5,.7),cex=2.1,pch=19,col=hcol,xlim=c(.7,5.3),xlab="",
     xaxt="n",bty="n",ylab="Home-Away effect",cex.lab=1.6,cex.axis=1.2)
abline(h=0,lty=2,col="gray")
arrows(1:5-off,bA.80[1,1,],1:5-off,bA.80[2,1,],angle=90,code=3,length=0,lwd=5,col=hcol)
arrows(1:5-off,bA.95[1,1,],1:5-off,bA.95[2,1,],angle=90,code=3,length=0,lwd=1,col=hcol)
points(1:5+off,bA.means[2,],cex=1.7,pch=19,col=wcol)
arrows(1:5+off,bA.80[1,2,],1:5+off,bA.80[2,2,],angle=90,code=3,length=0,lwd=5,col=wcol)
arrows(1:5+off,bA.95[1,2,],1:5+off,bA.95[2,2,],angle=90,code=3,length=0,lwd=1,col=wcol)
abline(v=c(1.5,2.5,3.5,4.5),lty=1,col="darkgray")
arrows(5.5,.1,5.5,.5,xpd="n",lwd=3)
text(5.75,.3,"Away Range",xpd="n",cex=2,srt=90)

mtext(c(expression("N"[R]),expression("N"[L]),expression("N"[B]),expression("N"[W]),expression("N"[D])),
      side=3,at=1:5,cex=2,line=-1)
mtext(c("Rubisco","","","Cell Wall","Defense"),at=1:5,line=3,cex=1.3,side=3)
mtext(c("Light","Bio-"),at=c(2,3),line=4,cex=1.3,side=3)
mtext(c("harvest","energetics"),at=c(2,3),line=2,cex=1.3,side=3)

polygon(x=c(2.15,2.15,3.9,3.9),y=c(-.7,-.45,-.45,-.7),xpd="n",col="white")
mtext("Fields",col=hcol,side=1,at=2.5,cex=1.5,line=.1)
mtext("Forests",col=wcol,side=1,at=3.5,cex=1.5,line=.1)

plot(1:5-off,bI.means[1,],ylim=c(-1,1),cex=2.1,pch=19,col=hcol,xlim=c(.7,5.3),xlab="",
     xaxt="n",bty="n",ylab="Native-Invader effect",cex.lab=1.6,cex.axis=1.2)
abline(h=0,lty=2,col="gray")
arrows(1:5-off,bI.80[1,1,],1:5-off,bI.80[2,1,],angle=90,code=3,length=0,lwd=5,col=hcol)
arrows(1:5-off,bI.95[1,1,],1:5-off,bI.95[2,1,],angle=90,code=3,length=0,lwd=1,col=hcol)
points(1:5+off,bI.means[2,],cex=1.7,pch=19,col=wcol)
arrows(1:5+off,bI.80[1,2,],1:5+off,bI.80[2,2,],angle=90,code=3,length=0,lwd=5,col=wcol)
arrows(1:5+off,bI.95[1,2,],1:5+off,bI.95[2,2,],angle=90,code=3,length=0,lwd=1,col=wcol)
abline(v=c(1.5,4.5),lty=1,col="darkgray")
arrows(2.5,-1.1,2.5,.9,col="darkgray",length=0)
arrows(3.5,-1.1,3.5,.9,col="darkgray",length=0)
arrows(5.5,.1,5.5,.8,xpd="n",lwd=3)
text(5.75,.4,"Invaders",xpd="n",cex=2,srt=90)

plot(1:5,bW.means,ylim=c(-1,1),cex=2.1,pch=19,col="black",xlim=c(.7,5.3),xlab="",
     xaxt="n",bty="n",ylab="Field-Forest effect",cex.lab=1.6,cex.axis=1.2)
abline(h=0,lty=2,col="gray")
arrows(1:5,bW.80[1,],1:5,bW.80[2,],angle=90,code=3,length=0,lwd=5,col="black")
arrows(1:5,bW.95[1,],1:5,bW.95[2,],angle=90,code=3,length=0,lwd=1,col="black")
abline(v=c(1.5,2.5,3.5,4.5),lty=1,col="darkgray")
arrows(5.5,.1,5.5,.75,xpd="n",lwd=3)
text(5.75,.4,"Forests",xpd="n",cex=2,srt=90)

#dev.off()

################################
#change HB model to more simplify estimate effects of 3 groups
#herb = as.numeric(woody==0)
invaway = as.numeric(group==1)
invhome = as.numeric(group==2)
YmU = cbind(dat$rubNperc,dat$lightNperc,dat$beNperc,dat$wallNperc,dat$defNperc)
Zscale2 = function(x) (x/sd(x))
YmZ2 = apply(YmU,2,Zscale2)

mod.n2 <- "model
{

  for(i in 1:N) { #observations

    Y[i,1:K] ~ dmnorm(mu.Y[i,],Omega[,]) #multivariate normal

    mu.Y[i,1:K] <- b0[1:K] + #mean - can't use with Z scaled Y
                    bW[1:K]*woody[i] +  #herb-woody difference for natives
                    bIA[1:K]*invaway[i] +
                    bIH[1:K]*invhome[i] +
                    bIAW[1:K]*invaway[i]*woody[i] +
                    bIHW[1:K]*invhome[i]*woody[i] +
                    #bGH[group[i],1:K]*herb[i] + #group effect (1-3) for herb species
                    #bGW[group[i],1:K]*woody[i] + #group effect (1-3) for woody species
                    bSpp[species[i],1:K] + 
                    bS[site[i],1:K]
  #}

    #predictions for fixed effects
    Ypred[i,1:K] <- bW[1:K]*woody[i] +
                    bIA[1:K]*invaway[i] +
                    bIH[1:K]*invhome[i] +
                    bIAW[1:K]*invaway[i]*woody[i] +
                    bIHW[1:K]*invhome[i]*woody[i]

    }

  #Prior for the precision matrix (inverse of the covariance matrix)
  Omega[1:K, 1:K] ~ dwish(R[,], K + 1) # Prior on the precision matrix
  Sigma[1:K, 1:K] <- inverse(Omega[,]) # Covariance matrix is the inverse of the precision matrix

  for(j in 1:K) { #loop over dependent variables

    #priors for fixed effect slopes (groups)
    bW[j] ~ dnorm(0, 0.001)
    b0[j] ~ dnorm(0, 0.001)
    bIA[j] ~ dnorm(0, 0.001)
    bIH[j] ~ dnorm(0, 0.001)
    bIAW[j] ~ dnorm(0, 0.001)
    bIHW[j] ~ dnorm(0, 0.001)

    #for(i in 1:3) {
    #  bGH[i,j] ~ dnorm(0, 0.001)
    #  bGW[i,j] ~ dnorm(0, 0.001)
      #bI[i,j] ~ dnorm(0, 0.001)
    #}

    #site REs
    for(i in 1:3) {
      bS[i,j] ~ dnorm(0,tauS1[j])
    }

    #species REs
    for(s in 1:Nspp) {
      bSpp[s,j] ~ dnorm(0,tauSpp[j])
      #bSpp[s,j] ~ dnorm(mu.bSpp[s,j],tauSpp[j])
      #mu.bSpp[s,j] <- bI[woodyS[s]+1,j]*invaderS[s] #invader effect must be here, else lost in species RE
    }
  
    #priors, RE variances
    #tau[j] <- sigma[j]^-2 
    #sigma[j] ~ dunif(0, 100) 
    tauSpp[j] <- sigmaSpp[j]^-2 
    sigmaSpp[j] ~ dunif(0, 100) 
    tauS1[j] <- sigmaS1[j]^-2 
    sigmaS1[j] ~ dunif(0, 100) 
}

}" #end model
write(mod.n2, "model2.txt")

#input lists for JAGS
params = c("bS","bW","b0","bIA","bIAW","bIH","bIHW","bSpp","Sigma","Ypred") #parameters to monitor
inits = function() list(b0=colMeans(Ym)) #starting values of fitted parameters
input = list(R=R,K=K,N=N,Y=Ym,site=site,woody=woody,invaway=invaway,invhome=invhome,species=species,Nspp=Nspp) #input data

#run JAGS model
jagsN2 <- jags(model = "model2.txt",data = input,param=params,#inits=inits,
              n.chains = 3, #number of separate MCMC chains 3
              n.iter =5000, #number of iterations per chain; 5000
              n.thin=5, #thinning 5
              n.burnin = 500) #number of initial iterations to discard 1000
plot(jagsN2)
detach.jags()
attach.jags(jagsN2)

#plot pairwise group comparison distributions, from posteriors
xnames = c("Away","Home","Native")

Ypred1 = Ypred[,,1]
IA.IH.H.diff1 = apply(Ypred1,1,function(x) mean(x[group==1&woody==0]) - mean(x[group==2&woody==0]))
IA.N.H.diff1 = apply(Ypred1,1,function(x) mean(x[group==1&woody==0]) - mean(x[group==3&woody==0]))
IH.N.H.diff1 = apply(Ypred1,1,function(x) mean(x[group==2&woody==0]) - mean(x[group==3&woody==0]))
IA.IH.W.diff1 = apply(Ypred1,1,function(x) mean(x[group==1&woody==1]) - mean(x[group==2&woody==1]))
IA.N.W.diff1 = apply(Ypred1,1,function(x) mean(x[group==1&woody==1]) - mean(x[group==3&woody==1]))
IH.N.W.diff1 = apply(Ypred1,1,function(x) mean(x[group==2&woody==1]) - mean(x[group==3&woody==1]))
  #NR
  #herbs: IA greater than IH (P=.002); IA marginally greater than N (P=.108)
  #woody: IA marginally greater than IH (P=.109)

Ypred2 = Ypred[,,2]
IA.IH.H.diff2 = apply(Ypred2,1,function(x) mean(x[group==1&woody==0]) - mean(x[group==2&woody==0]))
IA.N.H.diff2 = apply(Ypred2,1,function(x) mean(x[group==1&woody==0]) - mean(x[group==3&woody==0]))
IH.N.H.diff2 = apply(Ypred2,1,function(x) mean(x[group==2&woody==0]) - mean(x[group==3&woody==0]))
IA.IH.W.diff2 = apply(Ypred2,1,function(x) mean(x[group==1&woody==1]) - mean(x[group==2&woody==1]))
IA.N.W.diff2 = apply(Ypred2,1,function(x) mean(x[group==1&woody==1]) - mean(x[group==3&woody==1]))
IH.N.W.diff2 = apply(Ypred2,1,function(x) mean(x[group==2&woody==1]) - mean(x[group==3&woody==1]))
 #NL:
  #herb: IA marginally greater than IH (P=.073), IA marginally greater than native (P=0.059)
  #wody: IA marginally greater than N (P=0.086)

Ypred3 = Ypred[,,3]
IA.IH.H.diff3 = apply(Ypred3,1,function(x) mean(x[group==1&woody==0]) - mean(x[group==2&woody==0]))
IA.N.H.diff3 = apply(Ypred3,1,function(x) mean(x[group==1&woody==0]) - mean(x[group==3&woody==0]))
IH.N.H.diff3 = apply(Ypred3,1,function(x) mean(x[group==2&woody==0]) - mean(x[group==3&woody==0]))
IA.IH.W.diff3 = apply(Ypred3,1,function(x) mean(x[group==1&woody==1]) - mean(x[group==2&woody==1]))
IA.N.W.diff3 = apply(Ypred3,1,function(x) mean(x[group==1&woody==1]) - mean(x[group==3&woody==1]))
IH.N.W.diff3 = apply(Ypred3,1,function(x) mean(x[group==2&woody==1]) - mean(x[group==3&woody==1]))
  #NB:
  #herb: IH marginally greater than native (P=0.1)
  #woody: natives are marginally greater than both IA and IH (P=0.091, 0.073)

Ypred4 = Ypred[,,4]
IA.IH.H.diff4 = apply(Ypred4,1,function(x) mean(x[group==1&woody==0]) - mean(x[group==2&woody==0]))
IA.N.H.diff4 = apply(Ypred4,1,function(x) mean(x[group==1&woody==0]) - mean(x[group==3&woody==0]))
IH.N.H.diff4 = apply(Ypred4,1,function(x) mean(x[group==2&woody==0]) - mean(x[group==3&woody==0]))
IA.IH.W.diff4 = apply(Ypred4,1,function(x) mean(x[group==1&woody==1]) - mean(x[group==2&woody==1]))
IA.N.W.diff4 = apply(Ypred4,1,function(x) mean(x[group==1&woody==1]) - mean(x[group==3&woody==1]))
IH.N.W.diff4 = apply(Ypred4,1,function(x) mean(x[group==2&woody==1]) - mean(x[group==3&woody==1]))
  #NC: no differences

Ypred5 = Ypred[,,5]
IA.IH.H.diff5 = apply(Ypred5,1,function(x) mean(x[group==1&woody==0]) - mean(x[group==2&woody==0]))
IA.N.H.diff5 = apply(Ypred5,1,function(x) mean(x[group==1&woody==0]) - mean(x[group==3&woody==0]))
IH.N.H.diff5 = apply(Ypred5,1,function(x) mean(x[group==2&woody==0]) - mean(x[group==3&woody==0]))
IA.IH.W.diff5 = apply(Ypred5,1,function(x) mean(x[group==1&woody==1]) - mean(x[group==2&woody==1]))
IA.N.W.diff5 = apply(Ypred5,1,function(x) mean(x[group==1&woody==1]) - mean(x[group==3&woody==1]))
IH.N.W.diff5 = apply(Ypred5,1,function(x) mean(x[group==2&woody==1]) - mean(x[group==3&woody==1]))
  #ND:
  #herbs: IA greater than IH (0.022)
  #woody: no differences

Ypred6 = Ypred1 + Ypred2 + Ypred3 #all photo fractions
IA.IH.H.diff6 = apply(Ypred6,1,function(x) mean(x[group==1&woody==0]) - mean(x[group==2&woody==0]))
IA.N.H.diff6 = apply(Ypred6,1,function(x) mean(x[group==1&woody==0]) - mean(x[group==3&woody==0]))
IH.N.H.diff6 = apply(Ypred6,1,function(x) mean(x[group==2&woody==0]) - mean(x[group==3&woody==0]))
IA.IH.W.diff6 = apply(Ypred6,1,function(x) mean(x[group==1&woody==1]) - mean(x[group==2&woody==1]))
IA.N.W.diff6 = apply(Ypred6,1,function(x) mean(x[group==1&woody==1]) - mean(x[group==3&woody==1]))
IH.N.W.diff6 = apply(Ypred6,1,function(x) mean(x[group==2&woody==1]) - mean(x[group==3&woody==1]))
  #N photo:
  #herbs: IA is higher than both IH and N (P=0.013, P=0.017)
  #woody: NS

#Supplemental figure: pairwise comparisons of N fractions across groups

#pdf("Nallocation_pairwise_supplement.pdf",height=8,width=12)

#NR
par(mfrow=c(2,4),mar=c(5,7,4,1),oma=c(1,2,4,1))
x = IA.IH.H.diff1; hist(x,main="Away - Home",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
mtext("Fields",side=3,adj=0,cex=2,line=-2,xpd=NA,outer=T)
mtext(expression('N'[R]*' pairwise comparisons'),xpd=NA,outer=T,adj=0,cex=2.5,line=.8)
x = IA.N.H.diff1; hist(x,main="Away - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
x = IH.N.H.diff1; hist(x,main="Home - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
boxplot(dat$rubNperc[dat$woody==0]~dat$group[dat$woody==0],ylab="N fraction",xlab="",names=xnames,main="Fields")
x = IA.IH.W.diff1; hist(x,main="Away - Home",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
mtext("Forests",side=3,adj=0,cex=2,xpd=NA,outer=T,line=-28)
x = IA.N.W.diff1; hist(x,main="Away - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
x = IH.N.W.diff1; hist(x,main="Home - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
boxplot(dat$rubNperc[dat$woody==1]~dat$group[dat$woody==1],ylab="N fraction",xlab="",names=xnames,main="Forests")

#NL
par(mfrow=c(2,4),mar=c(5,7,4,1),oma=c(1,2,4,1))
x = IA.IH.H.diff2; hist(x,main="Away - Home",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
mtext("Fields",side=3,adj=0,cex=2,line=-2,xpd=NA,outer=T)
mtext(expression('N'[L]*' pairwise comparisons'),xpd=NA,outer=T,adj=0,cex=2.5,line=.8)
x = IA.N.H.diff2; hist(x,main="Away - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
x = IH.N.H.diff2; hist(x,main="Home - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
boxplot(dat$lightNperc[dat$woody==0]~dat$group[dat$woody==0],ylab="N fraction",xlab="",names=xnames,main="Fields")
x = IA.IH.W.diff2; hist(x,main="Away - Home",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
mtext("Forests",side=3,adj=0,cex=2,xpd=NA,outer=T,line=-28)
x = IA.N.W.diff2; hist(x,main="Away - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
x = IH.N.W.diff2; hist(x,main="Home - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
boxplot(dat$lightNperc[dat$woody==1]~dat$group[dat$woody==1],ylab="N fraction",xlab="",names=xnames,main="Forests")

#NB
par(mfrow=c(2,4),mar=c(5,7,4,1),oma=c(1,2,4,1))
x = IA.IH.H.diff3; hist(x,main="Away - Home",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
mtext("Fields",side=3,adj=0,cex=2,line=-2,xpd=NA,outer=T)
mtext(expression('N'[B]*' pairwise comparisons'),xpd=NA,outer=T,adj=0,cex=2.5,line=.8)
x = IA.N.H.diff3; hist(x,main="Away - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
x = IH.N.H.diff3; hist(x,main="Home - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
boxplot(dat$beNperc[dat$woody==0]~dat$group[dat$woody==0],ylab="N fraction",xlab="",names=xnames,main="Fields")
x = IA.IH.W.diff3; hist(x,main="Away - Home",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
mtext("Forests",side=3,adj=0,cex=2,xpd=NA,outer=T,line=-28)
x = IA.N.W.diff3; hist(x,main="Away - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
x = IH.N.W.diff3; hist(x,main="Home - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
boxplot(dat$beNperc[dat$woody==1]~dat$group[dat$woody==1],ylab="N fraction",xlab="",names=xnames,main="Forests")

#NC
par(mfrow=c(2,4),mar=c(5,7,4,1),oma=c(1,2,4,1))
x = IA.IH.H.diff4; hist(x,main="Away - Home",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
mtext("Fields",side=3,adj=0,cex=2,line=-2,xpd=NA,outer=T)
mtext(expression('N'[C]*' pairwise comparisons'),xpd=NA,outer=T,adj=0,cex=2.5,line=.8)
x = IA.N.H.diff4; hist(x,main="Away - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
x = IH.N.H.diff4; hist(x,main="Home - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
boxplot(dat$wallNperc[dat$woody==0]~dat$group[dat$woody==0],ylab="N fraction",xlab="",names=xnames,main="Fields")
x = IA.IH.W.diff4; hist(x,main="Away - Home",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
mtext("Forests",side=3,adj=0,cex=2,xpd=NA,outer=T,line=-28)
x = IA.N.W.diff4; hist(x,main="Away - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
x = IH.N.W.diff4; hist(x,main="Home - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.15,line=-.5,col="red",cex=.7)
boxplot(dat$wallNperc[dat$woody==1]~dat$group[dat$woody==1],ylab="N fraction",xlab="",names=xnames,main="Forests")

#ND
par(mfrow=c(2,4),mar=c(5,7,4,1),oma=c(1,2,4,1))
x = IA.IH.H.diff5; hist(x,main="Away - Home",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.02,line=-.5,col="red",cex=.7)
mtext("Fields",side=3,adj=0,cex=2,line=-2,xpd=NA,outer=T)
mtext(expression('N'[D]*' pairwise comparisons'),xpd=NA,outer=T,adj=0,cex=2.5,line=.8)
x = IA.N.H.diff5; hist(x,main="Away - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.2,line=-.5,col="red",cex=.7)
x = IH.N.H.diff5; hist(x,main="Home - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.2,line=-.5,col="red",cex=.7)
boxplot(dat$defNperc[dat$woody==0]~dat$group[dat$woody==0],ylab="N fraction",xlab="",names=xnames,main="Fields")
x = IA.IH.W.diff5; hist(x,main="Away - Home",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.02,line=-.5,col="red",cex=.7)
mtext("Forests",side=3,adj=0,cex=2,xpd=NA,outer=T,line=-28)
x = IA.N.W.diff5; hist(x,main="Away - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.2,line=-.5,col="red",cex=.7)
x = IH.N.W.diff5; hist(x,main="Home - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.2,line=-.5,col="red",cex=.7)
boxplot(dat$defNperc[dat$woody==1]~dat$group[dat$woody==1],ylab="N fraction",xlab="",names=xnames,main="Forests",ylim=c(0,.01))

#Nphoto
par(mfrow=c(2,4),mar=c(5,7,4,1),oma=c(1,2,4,1))
x = IA.IH.H.diff6; hist(x,main="Away - Home",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.02,line=-.5,col="red",cex=.7)
mtext("Fields",side=3,adj=0,cex=2,line=-2,xpd=NA,outer=T)
mtext(expression('N'[photo]*' pairwise comparisons'),xpd=NA,outer=T,adj=0,cex=2.5,line=.8)
x = IA.N.H.diff6; hist(x,main="Away - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.2,line=-.5,col="red",cex=.7)
x = IH.N.H.diff6; hist(x,main="Home - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.2,line=-.5,col="red",cex=.7)
boxplot(dat$rubNperc[dat$woody==0]+dat$lightNperc[dat$woody==0]+dat$beNperc[dat$woody==0]~dat$group[dat$woody==0],ylab="N fraction",xlab="",names=xnames,main="Fields")
x = IA.IH.W.diff6; hist(x,main="Away - Home",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.02,line=-.5,col="red",cex=.7)
mtext("Forests",side=3,adj=0,cex=2,xpd=NA,outer=T,line=-28)
x = IA.N.W.diff6; hist(x,main="Away - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.2,line=-.5,col="red",cex=.7)
x = IH.N.W.diff6; hist(x,main="Home - Native",xlab=""); abline(v=0,col="red",lwd=3); mtext(round(length(x[x>0])/length(x),3),at=.2,line=-.5,col="red",cex=.7)
boxplot(dat$rubNperc[dat$woody==1]+dat$lightNperc[dat$woody==1]+dat$beNperc[dat$woody==1]~dat$group[dat$woody==1],ylab="N fraction",xlab="",names=xnames,main="Forests")

#dev.off()



########other trait comparisons

#stat models
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)

#SLA
sla.mod = lmer ( log(SLA2) ~group*woody + (1|species) + (1|site),dat)
summary(sla.mod)
anova(sla.mod)
  #woody P<0.03
  #post-hoc, adjust P values manually with hoch test
  t.test(log(dat$SLA2[dat$group=="invaway"&dat$woody==0]),log(dat$SLA2[dat$group=="invhome"&dat$woody==0]))
  t.test(log(dat$SLA2[dat$group=="invaway"&dat$woody==0]),log(dat$SLA2[dat$group=="native"&dat$woody==0]))
  t.test(log(dat$SLA2[dat$group=="invhome"&dat$woody==0]),log(dat$SLA2[dat$group=="native"&dat$woody==0]))
  pvec = c(.1,.5,.6)  
  p.adjust(pvec,method="hochberg")
    #NS
  t.test(log(dat$SLA2[dat$group=="invaway"&dat$woody==1]),log(dat$SLA2[dat$group=="invhome"&dat$woody==1]))
  t.test(log(dat$SLA2[dat$group=="invaway"&dat$woody==1]),log(dat$SLA2[dat$group=="native"&dat$woody==1]))
  t.test(log(dat$SLA2[dat$group=="invhome"&dat$woody==1]),log(dat$SLA2[dat$group=="native"&dat$woody==1]))
  pvec = c(..4,.6,.9)  
  p.adjust(pvec,method="hochberg")
    #NS

#Narea  
  mod = lmer ( (Narea) ~group*woody + (1|species) + (1|site),dat)
  summary(mod)
  anova(mod)
    #woody P-0.00145
  #post-hoc, adjust P values manually with hoch test
  x = dat$Narea
  t.test(x[dat$group=="invaway"&dat$woody==0],x[dat$group=="invhome"&dat$woody==0])
  t.test(x[dat$group=="invaway"&dat$woody==0],x[dat$group=="native"&dat$woody==0])
  t.test(x[dat$group=="invhome"&dat$woody==0],x[dat$group=="native"&dat$woody==0])
  pvec = c(.4,.5,.9)  
  p.adjust(pvec,method="hochberg")
  #NS
  t.test(x[dat$group=="invaway"&dat$woody==1],x[dat$group=="invhome"&dat$woody==1])
  t.test(x[dat$group=="invaway"&dat$woody==1],x[dat$group=="native"&dat$woody==1])
  t.test(x[dat$group=="invhome"&dat$woody==1],x[dat$group=="native"&dat$woody==1])
  pvec = c(.05857,.8885,.1807)  
  p.adjust(pvec,method="hochberg")
  #NS after adjustment
  
#Nmass  
  mod = lmer ( (Nmass) ~group*woody + (1|species) + (1|site),dat)
  summary(mod)
  anova(mod)
  #NS
  #post-hoc, adjust P values manually with hoch test
  x = dat$Nmass
  t.test(x[dat$group=="invaway"&dat$woody==0],x[dat$group=="invhome"&dat$woody==0])
  t.test(x[dat$group=="invaway"&dat$woody==0],x[dat$group=="native"&dat$woody==0])
  t.test(x[dat$group=="invhome"&dat$woody==0],x[dat$group=="native"&dat$woody==0])
  pvec = c(.2,.8,.9)
  p.adjust(pvec,method="hochberg")
  #NS
  t.test(x[dat$group=="invaway"&dat$woody==1],x[dat$group=="invhome"&dat$woody==1])
  t.test(x[dat$group=="invaway"&dat$woody==1],x[dat$group=="native"&dat$woody==1])
  t.test(x[dat$group=="invhome"&dat$woody==1],x[dat$group=="native"&dat$woody==1])
  pvec = c(.9,.2,.9)  
  p.adjust(pvec,method="hochberg")
  #NS 
  
dat$Rubisco.ug = Rubisco.ug  
rub.mod = lmer ( log(Rubisco.ug) ~group*woody + (1|species) + (1|site),dat)
summary(rub.mod)
anova(rub.mod)
emmeans(rub.mod, list(pairwise ~ group), adjust = "tukey")
summary(glht(rub.mod,linfct = mcp(group = "Tukey"), test = adjusted("hochberg")))
  #alternatively, could just adjust P values manually with hoch test
  t.test(log(Rubisco.ug[dat$group=="invaway"&dat$woody==0]),log(Rubisco.ug[dat$group=="invhome"&dat$woody==0]))
  t.test(log(Rubisco.ug[dat$group=="invaway"&dat$woody==0]),log(Rubisco.ug[dat$group=="native"&dat$woody==0]))
  t.test(log(Rubisco.ug[dat$group=="invhome"&dat$woody==0]),log(Rubisco.ug[dat$group=="native"&dat$woody==0]))
  pvec = c(.00011,.03126,.03476)  
  p.adjust(pvec,method="hochberg")
    #A-H: P=0.003; I-N: P=0.035; H-N: P=0.035 (a,b,c)
  t.test(log(Rubisco.ug[dat$group=="invaway"&dat$woody==1]),log(Rubisco.ug[dat$group=="invhome"&dat$woody==1]))
  t.test(log(Rubisco.ug[dat$group=="invaway"&dat$woody==1]),log(Rubisco.ug[dat$group=="native"&dat$woody==1]))
  t.test(log(Rubisco.ug[dat$group=="invhome"&dat$woody==1]),log(Rubisco.ug[dat$group=="native"&dat$woody==1]))
  pvec = c(.064,.115,.966)  
  p.adjust(pvec,method="hochberg")
    #all NS (a,a,a)

#chl
  dat$chl = as.numeric(as.character(dat$Total.Chlorophyll))
  mod = lmer ( (chl) ~group*woody + (1|species) + (1|site),dat)
    #similar resids if Y is logged or not
  summary(mod)
  anova(mod)
  #NS
  #post-hoc, adjust P values manually with hoch test
  x = dat$chl
  t.test(x[dat$group=="invaway"&dat$woody==0],x[dat$group=="invhome"&dat$woody==0])
  t.test(x[dat$group=="invaway"&dat$woody==0],x[dat$group=="native"&dat$woody==0])
  t.test(x[dat$group=="invhome"&dat$woody==0],x[dat$group=="native"&dat$woody==0])
  pvec = c(.0001156,.002659,.75444)
  p.adjust(pvec,method="hochberg")
    #A-H: P-0.0003; #A-N: P=0.0053; #H-N: P=.754 (a,b,b)
  t.test(x[dat$group=="invaway"&dat$woody==1],x[dat$group=="invhome"&dat$woody==1])
  t.test(x[dat$group=="invaway"&dat$woody==1],x[dat$group=="native"&dat$woody==1])
  t.test(x[dat$group=="invhome"&dat$woody==1],x[dat$group=="native"&dat$woody==1])
  pvec = c(.8895,.00161,.0001287)  
  p.adjust(pvec,method="hochberg")  
    #A-H: P=.8895; #A-N: P=0.00322; #H-N: P=.000386 (a,a,b)

  
###graphing traits

library(RColorBrewer)
nat.col = brewer.pal(3,"Dark2")[1]
home.col = brewer.pal(3,"Dark2")[2]
away.col = brewer.pal(3,"Dark2")[3]
boxcols = c(away.col,home.col,nat.col,away.col,home.col,nat.col)

#pdf("traits_box.pdf",heigh=9,width=6)

par(mfcol=c(3,2),mar=c(.5,2,1,2.5),oma=c(3,3,2,1))
boxplot(SLA2~group+woody,dat,at=c(1,2,3,4.5,5.5,6.5),col=boxcols,ylab="",
        ylim=c(90,680),names=rep("",6),xlab="",xaxt="n",xpd=NA)
abline(v=3.75,col="gray")
mtext(expression('SLA (g cm'^-2*')'),side=2,line=2)
mtext(c("Fields","Forests"),side=3,at=c(2,5.5),cex=1.4)
text(6,660,"Habitat*",cex=1.5)
boxplot(Narea~group+woody,dat,at=c(1,2,3,4.5,5.5,6.5),col=boxcols,ylab="",
        names=rep("",6),xlab="",xaxt="n",xpd=NA)
mtext(expression('N'[area]* ' (g m'^-2*')'),side=2,line=2)
abline(v=3.75,col="gray")
text(6,3,"Habitat*",cex=1.5)
boxplot(Nmass*100~group+woody,dat,at=c(1,2,3,4.5,5.5,6.5),col=boxcols,ylab="",
        names=rep("",6),xlab="",xaxt="n",xpd=NA)
mtext(expression('N'[mass]* ' (%)'),side=2,line=2)
abline(v=3.75,col="gray")
box1 = boxplot(Rubisco.ug~group+woody,dat,at=c(1,2,3,4.5,5.5,6.5),col=boxcols,ylab="",
        names=rep("",6),xlab="",xaxt="n",xpd=NA,ylim=c(0,700))
mtext(expression('Rubisco (ug g'^-1*')'),side=2,line=2)
abline(v=3.75,col="gray")
mtext(c("Fields","Forests"),side=3,at=c(2,5.5),cex=1.4)
text(2.5,670," Group**",cex=1.5)
text(c(1,2,3,4.5,5.5,6.5),box1$stats[5,]+40,c("a","b","c","","",""),col=boxcols,cex=1.5) 
#text(4.4,660,"NS",cex=1.5)
box1 = boxplot(chl~group+woody,dat,at=c(1,2,3,4.5,5.5,6.5),col=boxcols,ylab="",
        names=rep("",6),xlab="",xaxt="n",xpd=NA,ylim=c(0,70))
mtext(expression('chl (ug cm'^-2*')'),side=2,line=2)
abline(v=3.75,col="gray")
text(c(1,2,3,4.5,5.5,6.5),box1$stats[5,]+4,c("a","b","b","a","a","b"),col=boxcols,cex=1.5)  
text(6,68,"Habitat*",cex=1.5)
text(2.5,68," Group***",cex=1.5)
plot(1:30,1:30,col="white",xaxt="n",yaxt="n",bty="n",ylim=c(0,30),xlim=c(0,30))  
text(12,19,"Invaders (away)",cex=1.9,col=away.col,xpd=NA)
text(12,15,"Invaders (home)",cex=1.9,col=home.col,xpd=NA)
text(12,11,"Natives",cex=1.9,col=nat.col,xpd=NA)

#dev.off()




##Extras

###This creates a supplemental table of home-away differences of all species for all traits.

#create target list, woodies first then herbs
wt = table(dat$Species[dat$woody==1])
wn = names(wt[wt>2])
ht = table(dat$Species[dat$woody==0])
hn = names(ht[ht>2])
spplist = c(wn,hn)

#traits
#tr.select = c("SLA2","leafN","leafC","leafCN","tot.prot","wall","fib.tot","alk","cyan")
dat$PNUE = dat$Asat/dat$Narea
dat$PNUE100 = dat$A.q100/dat$Narea
tr.select = c("SLA","Narea","Nmass","Rubisco.ug","chl","Vcmax","Jmax","Asat","Rd","alpha","PNUE","PNUE100",
            "rubNperc","lightNperc","beNperc","wallNperc","defNperc")
colnames = c("Species",paste(rep(tr.select,each=3),rep(c("H","A","P"),times=length(tr.select))))

spptab = matrix("",nrow=length(spplist),ncol=length(colnames))
colnames(spptab) = colnames

for(i in 1:length(spplist)) {
  #select species
  sp = spplist[i]
  spptab[i,1] = sp
  
  for(j in 1:length(tr.select)) {
    #select trait
    d = na.exclude(dat[dat$Species==sp,c("Species","homeaway",tr.select[j])])
    
    #mean and s.e. of trait in home and away range
    m = tapply(d[,3],d[,2],function(x)mean(x,na.rm=T))
    se = tapply(d[,3],d[,2],function(x)sd(x,na.rm=T)/sqrt(length(x[!is.na(x)])))
    
    #t test P value
    x = table(d[,2])
    if(length(x) > 1 & min(x) > 1) p = t.test(d[,3]~d[,2])$p.value else p = "NA"
    if(p>0.05) a = "" else {if(p>0.01) a = "*" else {if(p>0.001) a = "**" else a = "***"}}
    
    #insert in output table
    spptab[i,j*3-1] = paste(signif(m["home"],3),signif(se["home"],2),sep=" \u00B1 " ) #home
    spptab[i,j*3] = paste(signif(m["away"],3),signif(se["away"],2),sep=" \u00B1 " ) #away
    spptab[i,j*3+1] = a
  }
}
#write.table(spptab,"spptabN.txt",sep=";")

#trait table, with overall mean and S.D., sample size
msdvec = rep(0,length(tr.select))
Nvec = rep(0,length(tr.select))
for(j in 1:length(tr.select)) {
  m = mean(dat[,tr.select[j]],na.rm=T)
  sd = sd(dat[,tr.select[j]],na.rm=T)
  x = dat[,tr.select[j]]
  N = length(x[!is.na(x)])
  
  msdvec[j] = paste(signif(m,3),signif(sd,2),sep=" \u00B1 " )
  Nvec[j] = N
  
}
trtab = data.frame(tr.select,msdvec,Nvec)

#alpha was only calculated for woodies, reduce to N=161
trtab[10,3] = 161
#write.table(trtab,"trtabN.txt",sep=";")

