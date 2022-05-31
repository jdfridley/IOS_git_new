#Deriving photosynthetic properties, NSF-IOS project, all regions
#light curves only, because alpha was not converging for some samples
#JDF 5-5-22

#### Libraries
library(plantecophys)
library(nlme)
library(lme4)
library(R2jags)
library(doBy)
library(RCurl)
library(bayesplot)
library(ggplot2)
library(mcmcplots)
library(gridExtra)
library(googlesheets4)

#### Dataset: on OneDrive

master = read_sheet("https://docs.google.com/spreadsheets/d/1QLL5AUeP-fHar0HRfDjDP0Mw25FvowqidiCyEbBA5og/edit#gid=0")
dat = read.csv("/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/all_licor_data2.csv")
#version 2 is the same as v1, with one licor file (Fallopia F7) added
str(dat)
dat$region = substr(dat$site,1,1)
#11592 observations, 95 cols
dat$Date = as.Date(dat$date,format="%m/%d/%y")
dat$sppcode[dat$sppcode=="amrt"] = "amart" #fix typo

#### Variable inspection

#unique IDs
length(unique(dat$ID)) #419
length(unique(dat$filename)) #421
length(unique(master$ID)) #449
setdiff(master$ID,dat$ID) #differences explained in spreadsheet: 30 files not available for various reasons
#number of samples per species, per region
Fsub = tapply(dat$ID[dat$region=="F"],dat$species[dat$region=="F"],function(x)length(unique(x)))
Esub = tapply(dat$ID[dat$region=="E"],dat$species[dat$region=="E"],function(x)length(unique(x)))
Jsub = tapply(dat$ID[dat$region=="J"],dat$species[dat$region=="J"],function(x)length(unique(x)))
Fsub; Esub; Jsub

#Anet (umol CO2 per m2 per s)
hist(dat$Photo) #values >40 exceptional
plot(dat$Area,dat$Photo) #don't seem to be caused by errors of leaf Area (except one extreme value ~100)
plot(dat$Photo)
summary(dat$Photo)
#few outliers that should be ignored
dat = dat[dat$Photo<70,]
dat = dat[dat$Photo>(-10),]

#Conductance to H2O (mol H20 per m2 per s)
hist(dat$Cond) #substantial right skew
summary(dat$Cond)

#Intercelluar CO2 conc (umol CO2 per mol)
hist(dat$Ci) #seems good
summary(dat$Ci)
dat = dat[dat$Ci>(-10),]
summary(dat$Ci)

#convert Ci to Pa units if needed
dat$Ci_Pa = dat$Ci * dat$Press/ 1000

#Transpiration rate (mm H2O per m2 per s)
hist(dat$Trmmol) #ok

#VPD based on leaf temp (kPa)
hist(dat$VpdL) #ok

#Leaf area (cm2)
hist(dat$Area) #a few values ~4: extremely high

#Stomatal Ratio (0-1; 1 is all stomata on underside of leaf)
table(dat$StmRat) #which are 0? does it make a difference? don't think so...
table(dat$StmRat,dat$region)

#Boundary layer conductance (mol per m2 per s)
hist(dat$BLCond) #are values > 3 due to smaller leaf area?
plot(dat$Area,dat$BLCond) #yes, linear decrease for some

#Air temp
hist(dat$Tair)
quantile(prob=c(.025,.975),dat$Tair) #95% within 24 and 33C 
boxplot(Tair~region,dat) #significantly cooler air temp in France

#Leaf temp
hist(dat$Tleaf)
quantile(prob=c(.025,.975),dat$Tleaf) #95% within 24 and 33; need to consider in curve fitting?
length(dat$Tleaf[dat$Tleaf>28&dat$Tleaf<30])/length(dat$Tleaf)
boxplot(Tleaf~region,dat) #significantly cooler air temp in France
hist(dat$Tleaf[dat$region=="E"])  
#nearly all ENA Tleaf at 30
hist(dat$Tleaf[dat$region=="J"])  
#nearly all Japan Tleaf at 30
hist(dat$Tleaf[dat$region=="F"])  
plot(dat$Date[dat$region=="F"],dat$Tleaf[dat$region=="F"])
#25C temps used in summer 2020 in France, presumably when necessary to reduce VPD due to very dry conditions

#ignoring TBlk, CO2R, CO2S, H2OR, H2OS, RH_R, RH_S

#Flow rate
hist(dat$Flow) #constant at 500 umol per s

#PAR on leaf (umol photos per m2 per s)
hist(dat$PARi) #ok

#ambient PAR (umol photos per m2 per s)
hist(dat$PARo) #majority of samples in shade or heavy cloud

#Atm pressure (kPa)
hist(dat$Press) #std atm is 101.3; some relatively low values must be higher elevation

#Stability status of licor
table(dat$Status) #nearly all 111115; 42 are 111135; 67 are 111125, 5 are 111105 - these are CO2 mixer flags
#11 are 111215: the 2 is a flow control flag

#Labels:
table(dat$filename)
table(dat$date)
hist(dat$Date,breaks=20)
length(table(dat$site)) #122 sites: 39 ENA, 57 France, 26 Japan (8 less than actual named sites due to missing files) 
table(dat$ID)
table(dat$species); length(table(dat$species))  #49 species!
table(dat$sppcode); length(table(dat$sppcode))    
table(dat$species,dat$site) #number of rows per species per site
sort(colSums(table(dat$species,dat$site)>0)) #number of species per site; about 1/3 have just 1, max is 12
#so, need site RE?

#### Temperature adjusted coefficients (Mason's code)

# Constants published in Sharkey et al (2007) Plant Cell Env 30: 1035-1040 
# Measured using transgenic tobacco (ASSUMED to be similar across higher plants)
# Ci units in Pa; Sharkey et al (2007) recommend partial pressures
# **Be sure units are correct for your input data** (Ci is in Pa or ppm?)

R=0.008314 #(kJ mol^-1 K^-1)
dat$Kc=exp(35.9774-80.99/(R*(dat$Tleaf+273.15))) #Michaelis-Menten constant for Rubisco for O2 (Pa)
dat$Ko=exp(12.3772-23.72/(R*(dat$Tleaf+273.15))) #Michaelis-Menten constant for Rubisco for CO2 (kPa)
dat$GammaStar=exp(11.187-24.46/(R*(dat$Tleaf+273.15))) #Photorespiration compensation point (Pa)
dat$O=21 #oxygen (O2) partial pressure (kPa)  

#### Which species are A-Ci only, which are both A-Ci and A-q?

#by species, take q values, estimate variance or range
spp.list = unique(dat$species)
out = matrix(0,nrow=length(spp.list),ncol=16)
rownames(out) = spp.list
colnames(out) = seq(0,1500,length=16)
for(i in 1:length(spp.list)) {
  ds = dat[dat$species==spp.list[i],]
  out[i,] = table(cut(ds$PARi,seq(-50,1600,by=100)))
}

#if more than 3 readings at PAR=300, light curves were attempted
spp.forest = spp.list[out[,4]>3]  #21 out of 49 species
spp.field = spp.list[out[,4]<4]   #28 out of 49 species


##create limited dataset for light-curve fitting, using only CO2R values around 400 ppm
dat2 = dat[dat$Ci>360&dat$Ci<440,]




############################
##HB via JAGS, all species, A-q curves only

spp1 = "Rosa multiflora"
spp2 = "Euonymus alatus"
#spp3 = "Artemisia vulgaris"
df = dat2[dat2$species==spp1|dat2$species==spp2,]
#df = dat2
N = dim(df)[1]
ind = as.numeric(as.factor(df$filename)) #grouping vector (individual)
N.indiv = max(ind)
spp = as.numeric(as.factor(df$species))
N.spp = max(spp)
ind.spp = as.numeric(apply(table(ind,spp)>0,1,function(x)names(which(x==T)))) #lists species to which an individual belongs, numerically

mod.photo <- "model
{
    #Priors
    alpha.int ~ dnorm(0.24,100) #strongish prior, low variation across species (Feng & Dietze 2013)
    theta.int ~ dunif(0,1) #allowable range of theta

    #for Rd treat as fixed (unpooled) effect not random (difficult to estimate with >0 constraint)
    for(i in 1:N.indiv) {
      Rd.int[i] ~ dnorm(0,1)   
    }
    
    ind.tau.alpha <- ind.sigma.Jmax^-2 
    ind.sigma.alpha ~ dunif(0, 100) 
    spp.tau.alpha <- spp.sigma.Jmax^-2 
    spp.sigma.alpha ~ dunif(0, 100)

    ind.tau.theta <- ind.sigma.Jmax^-2 
    ind.sigma.theta ~ dunif(0, 100) 
    spp.tau.theta <- spp.sigma.Jmax^-2 
    spp.sigma.theta ~ dunif(0, 100)

    #residual error
    tau <- sigma^-2 #coverts sd to precision
    sigma ~ dunif(0, 100)  #uniform prior for standard deviation

    #NOTE: site-level RE ignored (for now)

    for(i in 1:N) { #lowest level (N=N)
        
        Anet[i] ~ dnorm(mu[i],tau) #residual error
        
        #mu[i] <- min(mu.v[i],mu.j[i]) - Rd.int[ind[i]] #minimum of RuBP and Rubisco limitation; TPU limitation ignored

        mu[i] <- (alpha[i]*q[i] + Amax[i] - sqrt (  (alpha[i]*q[i] + Amax[i])^2 - 4*theta[i]*alpha[i]*q[i]*Amax[i] )) / (2*theta[i]) - Rd.int[ind[i]]
        
        # Photo params that have individual pooling nested within species pooling
        alpha[i] <- alpha.int + b0.ind.alpha[ind[i]] #here b0.ind.alpha is informed by b0.spp.alpha
        theta[i] <- theta.int + b0.ind.theta[ind[i]] 
        Amax[i] <- Amax.int + b0.ind.Amax[ind[i]] 

        #J light dependency according to Tenhunen et al 1976
	      #J[i]<-(alpha[i]*q[i]/(sqrt(1+(alpha[i]*alpha[i]*q[i]*q[i])/(Jmax[i]*Jmax[i]))))
	      
    }

	  #random intercept for individual effects, nested within species
      for(i in 1:N.indiv) {
        b0.ind.Vcmax[i] ~ dnorm(b0.spp.Vcmax[ind.spp[i]],ind.tau.Vcmax) #individual value sampled from species' distribution
        b0.ind.Jmax[i] ~ dnorm(b0.spp.Jmax[ind.spp[i]],ind.tau.Jmax) #individual value sampled from species' distribution
        b0.ind.alpha[i] ~ dnorm(b0.spp.alpha[ind.spp[i]],ind.tau.alpha) #individual value sampled from species' distribution
        b0.ind.theta[i] ~ dnorm(b0.spp.theta[ind.spp[i]],ind.tau.theta) #individual value sampled from species' distribution
        
        #output monitoring
        alpha.out[i] <- alpha.int + b0.ind.alpha[i]
        theta.out[i] <- theta.int + b0.ind.theta[i]
      }  

    #random intercept for species effects
      for(i in 1:N.spp) {
        b0.spp.alpha[i] ~ dnorm(0,spp.tau.alpha)
        b0.spp.theta[i] ~ dnorm(0,spp.tau.theta)
        
        #monitor spp-level values
        alpha.spp[i] <- alpha.int + b0.spp.alpha[i]
        theta.spp[i] <- theta.int + b0.spp.theta[i]
      }


}" #end model
  write(mod.photo, "model.txt")
  
  #input lists for JAGS
  params = c("Asat","Vcm.out","Jm.out","alpha.out","Rd.int","sigma","Vcm.spp","Jm.spp","alpha.spp","Vcmax.int","Jmax.int","alpha.int") #parameters to monitor
  inits = function() list(Vcmax.int=rnorm(1,100,10),Jmax.int=rnorm(1,75,10),alpha.int=runif(1),Rd.int=rnorm(N.indiv)) #starting values of fitted parameters
  input = list(N=N,Anet=df$Photo,Ci_Pa=df$Ci_Pa,q=df$PARi,GammaStar=df$GammaStar,Kc=dat$Kc,Ko=dat$Ko,O=df$O,ind=ind,N.indiv=N.indiv,ind.spp=ind.spp,N.spp=N.spp) #input data
  
  #run JAGS model
  jags.p <- jags(model = "model.txt",data = input,inits=inits,param=params,
                 n.chains = 3, #number of separate MCMC chains
                 n.iter =3000, #number of iterations per chain
                 n.thin=3, #thinning
                 n.burnin = 500) #number of initial iterations to discard
  
  jags.p
  
  