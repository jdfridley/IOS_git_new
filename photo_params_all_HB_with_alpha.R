###HB model of photosynthesis, A-Ci and A-q curves
#forest species only
#JDF April 2022


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
#dat = read.csv("/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/all_licor_data.csv")
dat = read.csv("all_licor_data.csv")
dat$region = substr(dat$site,1,1)
#11592 observations, 95 cols
dat$Date = as.Date(dat$date)#,format="%m/%d/%y")
dat$sppcode[dat$sppcode=="amrt"] = "amart" #fix typo

#### Variable inspection

#unique IDs
length(unique(dat$ID)) #419
length(unique(dat$filename)) #420
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


#### HB via JAGS, A-Ci and A-q curves estimated simultaneously for forest species

spp.forest

#spp1 = "Rosa multiflora"
#spp2 = "Lonicera morrowii"
#df = dat[dat$species==spp1|dat$species==spp2,]
df = dat[is.element(dat$species,spp.forest),]
#df = df[df$PARi>800,]
df = df[df$Ci>=0,]
N = dim(df)[1]
ind = as.numeric(as.factor(df$ID)) #grouping vector (individual)
N.indiv = max(ind)

#color codes for region
test = unique(cbind(as.numeric(ind),df$region))
test = as.data.frame(test); test$V1 = as.numeric(test$V1)
reg = test[order(test[,1]),2]
cols = rep("blue",length(reg))
cols[reg=="E"] = "darkgreen"
cols[reg=="J"] = "red2"

mod.photo <- "model
{
    #Priors
    Vcmax.int ~ dnorm(75,0.001)T(0,) #weak
    Jmax.int ~ dnorm(75,0.001)T(0,) #weak
    #J.int ~ dnorm(100,0.001) #very weak; J is for those plants without A-q curves
    alpha.int ~ dnorm(0.24,100) #strong prior, low variation across species (Feng & Dietze 2013)

    #for Rd, which can't be negative, treat as fixed (unpooled) effect not random (difficult to estimate with >0 constraint)
    for(i in 1:N.indiv) {
      Rd[i] ~ dnorm(0,1)T(0,) #note cannot take on negative values with T(0,)
    }
    
    #individual and species-level variance in photo params
    ind.tau.Vcmax <- ind.sigma.Vcmax^-2 
    ind.sigma.Vcmax ~ dunif(0, 100)
    #spp.tau.Vcmax <- spp.sigma.Vcmax^-2 
    #spp.sigma.Vcmax ~ dunif(0, 100)
    
    ind.tau.Jmax <- ind.sigma.Jmax^-2 
    ind.sigma.Jmax ~ dunif(0, 100) 
    #spp.tau.Jmax <- spp.sigma.Jmax^-2 
    #spp.sigma.Jmax ~ dunif(0, 100)

    #ind.tau.J <- ind.sigma.J^-2 
    #ind.sigma.J ~ dunif(0, 100) 
    #spp.tau.J <- spp.sigma.J^-2 
    #spp.sigma.J ~ dunif(0, 100)

    ind.tau.alpha <- ind.sigma.alpha^-2 
    ind.sigma.alpha ~ dunif(0, 100) 
    #spp.tau.alpha <- spp.sigma.alpha^-2 
    #spp.sigma.alpha ~ dunif(0, 100)

    #level 1 variance (error)
    tau <- sigma^-2 #coverts sd to precision
    sigma ~ dunif(0, 100)  #uniform prior for standard deviation

    #NOTE need to add site-level variation (non-nested)

    for(i in 1:N) {
        
        Anet[i] ~ dnorm(mu[i],tau)
        
        mu[i] <- min(mu.v[i],mu.j[i]) - Rd[ind[i]] #minimum of RuBP and Rubisco limitation; TPU limitation ignored
        
        Vcmax[i] <- Vcmax.int + b0.ind.Vcmax[ind[i]]

        Jmax[i] <- Jmax.int + b0.ind.Jmax[ind[i]]
        
        alpha[i] <- alpha.int + b0.ind.alpha[ind[i]]

        #J[i] <- J.int + b0.ind.J[ind[i]]

        # ---RuBisCO limited portion---
        mu.v[i] <- (Vcmax[i]*(Ci_Pa[i]-GammaStar[i]))/(Ci_Pa[i]+(Kc[i]*(1+(O[i]/Ko[i]))))
  
        # ---RUBP limited portion---
        #mu.j[i] <- ((Jmax[i]*(Ci_Pa[i]-GammaStar[i]))/((4*Ci_Pa[i])+(8*GammaStar[i])))
	      mu.j[i] <- ((J[i]*(Ci_Pa[i]-GammaStar[i]))/((4*Ci_Pa[i])+(8*GammaStar[i]))) 

        #J light dependency according to Tenhunen et al 1976
	      J[i]<-(alpha[i]*q[i]/(sqrt(1+(alpha[i]*alpha[i]*q[i]*q[i])/(Jmax[i]*Jmax[i]))))
  
    }

	  #random intercept for individual effect on Vcmax
      for(i in 1:N.indiv) {
        b0.ind.Vcmax[i] ~ dnorm(0,ind.tau.Vcmax) 
        b0.ind.Jmax[i] ~ dnorm(0,ind.tau.Jmax) 
        b0.ind.alpha[i] ~ dnorm(0,ind.tau.alpha)
        #b0.ind.J[i] ~ dnorm(0,ind.tau.J)
        
        #output monitoring
        Vcm.out[i] = Vcmax.int + b0.ind.Vcmax[i]
        Jm.out[i] = Jmax.int + b0.ind.Jmax[i]
        alpha.out[i] = alpha.int + b0.ind.alpha[i]
        #J.out[i] = J.int + b0.ind.J[i]
        #calculate Asat at 40 Ci_Pa and saturating light
        Asat[i] <-min(((Jm.out[i]*(40-4.275))/((4*40)+(8*4.275))),Vcm.out[i]*(40-4.275)/(40+(40.49*(1+(O/27.84))))) 

        #hypothesis testing

      }  

}" #end model
  write(mod.photo, "model.txt")
  
  #input lists for JAGS
  params = c("Asat","Vcm.out","Jm.out","alpha.out","Rd","sigma") #parameters to monitor
  inits = function() list(Vcmax.int=rnorm(1,100,10),Jmax.int=rnorm(1,75,10),alpha.int=runif(1),Rd=rnorm(N.indiv)) #starting values of fitted parameters
  input = list(N=N,Anet=df$Photo,Ci_Pa=df$Ci_Pa,q=df$PARi,GammaStar=df$GammaStar,Kc=dat$Kc,Ko=dat$Ko,O=df$O,ind=ind,N.indiv=N.indiv,region=as.numeric(as.factor(reg))) #input data
  
  #run JAGS model
  jags.p <- jags(model = "model.txt",data = input,inits=inits,param=params,
                 n.chains = 3, #number of separate MCMC chains
                 n.iter =5000, #number of iterations per chain
                 n.thin=3, #thinning
                 n.burnin = 1000) #number of initial iterations to discard
  
  jags.p
  
  
  #ind-level output
  detach.jags()
  attach.jags(jags.p)
  
  ind.out = summaryBy(filename~filename+date+site+species+sppcode+ID,df)
  ind.out = ind.out[,-7]
  ind.out = cbind(ind.out,data.frame(Asat=apply(Asat,2,mean),Vcmax=apply(Vcm.out,2,mean),Jmax=apply(Jm.out,2,mean),Rd=apply(Rd,2,median),alpha=apply(alpha.out,2,mean), Asat.se=apply(Asat,2,sd), Vcmax.se=apply(Vcm.out,2,sd),Jmax.se=apply(Jm.out,2,sd),Rd.se=apply(Rd,2,sd),alpha.se=apply(alpha.out,2,sd)))
  
  str(ind.out)
  