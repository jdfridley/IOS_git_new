#Deriving photosynthetic properties, NSF-IOS project, version 1, all regions
#JDF 3-30-22

#to do list
#1 reliable output for all individs, A-Ci curves
#3 summary posteriors of major contrasts (home-away, native-nonnative, etc)
  #fixed effect (eg of region) or random (with ind nested within?)
#4 integrate other data (N sources, light, etc)

#### Data folder (on Amiens iMac)
dfold = "/Users/fridley/Documents/IOS/IOS_FranceJapan/"

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

#master = read_sheet("https://docs.google.com/spreadsheets/d/1QLL5AUeP-fHar0HRfDjDP0Mw25FvowqidiCyEbBA5og/edit#gid=0")
master = read.csv(paste0(dfold,"NSF-IOSspreadsheetJDF4.csv"))
#dat = read.csv("/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/all_licor_data2.csv")
dat = read.csv(paste0(dfold,"all_licor_data2.csv"))
  #version 2 is the same as v1, with one licor file (Fallopia F7) added
str(dat)
dat$region = substr(dat$site,1,1)
  #11592 observations, 95 cols
dat$Date = as.Date(dat$date,format="%m/%d/%y")
dat$sppcode[dat$sppcode=="amrt"] = "amart" #fix typo
dat$ID[dat$filename=="2021-09-13-bifro-dewitt.xlsx"] = "bifro-E37-1" #fix typo

#### Variable inspection

  #unique IDs
  length(unique(dat$ID)) #419
    #some IDs only indicate separate licor files that should be combined under the same ID:
    dat$ID[dat$ID=="ropse-F51-2"] = "ropse-F51-1"
    dat$ID[dat$ID=="acneg-F37-2"] = "acneg-F37-1"
    dat$ID[dat$ID=="acpse-F41-2"] = "acpse-F41-1"
    dat$ID[dat$ID=="paspp-F28-2"] = "paspp-F28-1"
    
    #some dates need to be altered for duplicate samples on different days:
    dat$date[dat$ID=="acneg-F37-1"&dat$date=="6/10/20"] = "6/5/20"
    dat$date[dat$ID=="acpse-F41-1"&dat$date=="7/29/20"] = "7/28/20"

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
  dat$Press[is.na(dat$Press)] = 100 #two samples missing Press in ENA; nearly sea level
  
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
  

############################
  ##HB via JAGS, all species, combined A-Ci and A-q curves (similar structure to Heberling & Fridley
    #but with ind effects nested within species)

  spp1 = "Lonicera canadensis"
  spp2 = "Euonymus alatus"
  df = dat[dat$species==spp1|dat$species==spp2,]
  df = dat
  df = df[df$Ci>=0,]
  N = dim(df)[1]
  #ind = as.numeric(as.factor(df$filename)) #grouping vector (individual)
  ind = as.numeric(as.factor(df$ID)) #change to ID for duplicate licor files per ind
  N.indiv = max(ind)
  spp = as.numeric(as.factor(df$species))
  N.spp = max(spp)
  ind.spp = as.numeric(apply(table(ind,spp)>0,1,function(x)names(which(x==T)))) #lists species to which an individual belongs, numerically
  
  mod.photo <- "model
{
    #Priors
    Vcmax.int ~ dnorm(75,0.001)T(0,) #very weak
    Jmax.int ~ dnorm(100,0.001)T(0,) #very weak
    alpha.int ~ dnorm(0.24,100) #strong prior, low variation across species (Feng & Dietze 2013)

    #for Rd, which can't be negative, treat as fixed (unpooled) effect not random (difficult to estimate with >0 constraint)
    for(i in 1:N.indiv) {
      Rd.int[i] ~ dnorm(0,1)T(0,) #note cannot take on negative values with T(0,)
    }
    
    #individual and species-level variance in photo params
    ind.tau.Vcmax <- ind.sigma.Vcmax^-2 
    ind.sigma.Vcmax ~ dunif(0, 100)
    spp.tau.Vcmax <- spp.sigma.Vcmax^-2 
    spp.sigma.Vcmax ~ dunif(0, 100)
    
    ind.tau.Jmax <- ind.sigma.Jmax^-2 
    ind.sigma.Jmax ~ dunif(0, 100) 
    spp.tau.Jmax <- spp.sigma.Jmax^-2 
    spp.sigma.Jmax ~ dunif(0, 100)

    ind.tau.alpha <- ind.sigma.alpha^-2 
    ind.sigma.alpha ~ dunif(0, 100) 
    spp.tau.alpha <- spp.sigma.alpha^-2 
    spp.sigma.alpha ~ dunif(0, 100)

    #residual error
    tau <- sigma^-2 #coverts sd to precision
    sigma ~ dunif(0, 100)  #uniform prior for standard deviation

    #NOTE: site-level RE ignored (for now)

    for(i in 1:N) { #lowest level (N=N)
        
        Anet[i] ~ dnorm(mu[i],tau) #residual error
        
        mu[i] <- min(mu.v[i],mu.j[i]) - Rd.int[ind[i]] #minimum of RuBP and Rubisco limitation; TPU limitation ignored
        
        # Photo params that have individual pooling nested within species pooling
        Vcmax[i] <- Vcmax.int + b0.ind.Vcmax[ind[i]]  #here b0.ind.Vcmax is informed by b0.spp.Vcmax
        Jmax[i] <- Jmax.int + b0.ind.Jmax[ind[i]] #here b0.ind.Jmax is informed by b0.spp.Jmax
        alpha[i] <- alpha.int + b0.ind.alpha[ind[i]] #here b0.ind.alpha is informed by b0.spp.alpha

        # ---RuBisCO limited portion---
        mu.v[i] <- (Vcmax[i]*(Ci_Pa[i]-GammaStar[i]))/(Ci_Pa[i]+(Kc[i]*(1+(O[i]/Ko[i]))))
  
        # ---RUBP limited portion---
        #mu.j[i] <- ((Jmax[i]*(Ci_Pa[i]-GammaStar[i]))/((4*Ci_Pa[i])+(8*GammaStar[i])))
	      mu.j[i] <- ((J[i]*(Ci_Pa[i]-GammaStar[i]))/((4*Ci_Pa[i])+(8*GammaStar[i]))) 

        #J light dependency according to Tenhunen et al 1976
	      J[i]<-(alpha[i]*q[i]/(sqrt(1+(alpha[i]*alpha[i]*q[i]*q[i])/(Jmax[i]*Jmax[i]))))
	      
    }

	  #random intercept for individual effects, nested within species
      for(i in 1:N.indiv) {
        b0.ind.Vcmax[i] ~ dnorm(b0.spp.Vcmax[ind.spp[i]],ind.tau.Vcmax) #individual value sampled from species' distribution
        b0.ind.Jmax[i] ~ dnorm(b0.spp.Jmax[ind.spp[i]],ind.tau.Jmax) #individual value sampled from species' distribution
        b0.ind.alpha[i] ~ dnorm(b0.spp.alpha[ind.spp[i]],ind.tau.alpha) #individual value sampled from species' distribution
        
        #output monitoring
        Vcm.out[i] <- Vcmax.int + b0.ind.Vcmax[i] #includes spp constraint
        Jm.out[i] <- Jmax.int + b0.ind.Jmax[i]
        alpha.out[i] <- alpha.int + b0.ind.alpha[i]
        #calculate Asat at 40 Ci_Pa and saturating light
        Asat[i] <-min(((Jm.out[i]*(40-4.275))/((4*40)+(8*4.275))),Vcm.out[i]*(40-4.275)/(40+(40.49*(1+(O/27.84))))) 
        #calculate Ac-Aj Ci transition point (value of Ci_Pa of equal limitation by Rubisco and RUBP)
        Citr[i] <- ((Vcm.out[i]/Jm.out[i])*8*GammaStar[i] - (Kc[i]*(1+(O[i]/Ko[i])))) / (1-(4*(Vcm.out[i]/Jm.out[i])))
      }  

    #random intercept for species effects
      for(i in 1:N.spp) {
        b0.spp.Vcmax[i] ~ dnorm(0,spp.tau.Vcmax)
        b0.spp.Jmax[i] ~ dnorm(0,spp.tau.Jmax)
        b0.spp.alpha[i] ~ dnorm(0,spp.tau.alpha)
        
        #monitor spp-level values
        Vcm.spp[i] <- Vcmax.int + b0.spp.Vcmax[i]
        Jm.spp[i] <- Jmax.int + b0.spp.Jmax[i]
        alpha.spp[i] <- alpha.int + b0.spp.alpha[i]
      }

}" #end model
  write(mod.photo, "model.txt")
  
  #input lists for JAGS
  params = c("Asat","Vcm.out","Jm.out","alpha.out","Rd.int","Citr","sigma","Vcm.spp","Jm.spp","alpha.spp","Vcmax.int","Jmax.int","alpha.int") #parameters to monitor
  inits = function() list(Vcmax.int=rnorm(1,100,10),Jmax.int=rnorm(1,75,10),alpha.int=runif(1),Rd.int=rnorm(N.indiv)) #starting values of fitted parameters
  input = list(N=N,Anet=df$Photo,Ci_Pa=df$Ci_Pa,q=df$PARi,GammaStar=df$GammaStar,Kc=dat$Kc,Ko=dat$Ko,O=df$O,ind=ind,N.indiv=N.indiv,ind.spp=ind.spp,N.spp=N.spp) #input data
  
  #run JAGS model
  jags.p <- jags(model = "model.txt",data = input,inits=inits,param=params,
                 n.chains = 3, #number of separate MCMC chains
                 n.iter =100, #number of iterations per chain
                 n.thin=3, #thinning
                 n.burnin = 20) #number of initial iterations to discard
  
  jags.p
  
  #examine Rhats (8th col of output)
  hist(jags.p$BUGSoutput$summary[,8])
  summary(jags.p$BUGSoutput$summary[,8])
  #save(jags.p,file="jags.p_5.31.22.RData") #JAGS object saved
    #some convergence problems
  
  #create output spreadsheets
  detach.jags()
  #attach(jags.p,overwrite=T)
  attach.jags(jags.p)
  
  #species-level: note summaries of Rd are based on the median
  spp.out = data.frame(species=levels(as.factor(df$species)),Vcmax=apply(Vcm.spp,2,mean),Jmax=apply(Jm.spp,2,mean),Rd=tapply(apply(Rd.int,2,median),ind.spp,mean),
                       Vcmax.se=apply(Vcm.spp,2,sd),Jmax.se=apply(Jm.spp,2,sd),Rd.se=tapply(apply(Rd.int,2,median),ind.spp,sd))
  woody = c(0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,0,1,0) #7 woodies
  invader = c(1,1,1,1,0,1,0,1,0,1,0,1,0,0,0,0,1,0,0,0,1,0,1,0,0,0)
  plot(spp.out$Vcmax,spp.out$Jmax,col=woody+1,pch=19,cex=1.5) #woodies have lowest Jmax, Vcmax (except Cryptotaenia japonica) values
  plot(spp.out$Vcmax,spp.out$Jmax,col=invader+1,pch=19,cex=1.5) #invaders are higher 
  
  #ind-level
  #ind.out = summaryBy(filename~filename+ID+date+site+species+sppcode,df)
  ind.out = summaryBy(ID~ID+date+site+species+sppcode,df)
  ind.out = ind.out[,-6]
  ind.out = cbind(ind.out,data.frame(Asat=apply(Asat,2,mean),Vcmax=apply(Vcm.out,2,mean),Jmax=apply(Jm.out,2,mean),Rd=apply(Rd.int,2,median),alpha=apply(alpha.out,2,mean), Asat.se=apply(Asat,2,sd), Vcmax.se=apply(Vcm.out,2,sd),Jmax.se=apply(Jm.out,2,sd),Rd.se=apply(Rd.int,2,sd),alpha.se=apply(alpha.out,2,sd)))
  
  #same as csv file
  #write.csv(ind.out,"indoutHB3.csv")

  #combine with ecophys derived A-Ci dataset
  #epset = read.csv("photo_params_ecophys.csv")
  #out = merge(epset,ind.out,by.x="ID",by.y="ID",all=T)
  #dim(out)
  #head(out)
  #write.csv(out,"photo_params_HBandEcoPhys.csv")
    #'x' values are ecophys derived, 'y' values are HB derived

  out = read.csv("photo_params_HBandEcoPhys.csv") #load in dataset with ecophys and former HB (A-Ci only) values
  
  test = merge(out,ind.out,by.x="ID",by.y="ID")
  #Vcmax and Jmax values using A-Ci only ('y' vars) in HB vs. using A-Ci and A-q HB curves (as above)
  plot(test$Vcmax.y,test$Vcmax,col=out$woody+2); abline(0,1) #woodies are green
  plot(test$Jmax.y,test$Jmax,col=out$woody+2); abline(0,1) #woodies are green
  plot(test$Rd.y,test$Rd,col=out$woody+2); abline(0,1) #woodies are green
  hist(test$alpha[test$alpha<1&test$alpha>0])
  hist(test$Asat)
  
  #compare HB (combined A-Ci and A-q) to ecophys
  plot(test$Vcmax.x,test$Vcmax,col=out$woody+2); abline(0,1) #woodies are green
  plot(test$Jmax.x,test$Jmax,col=out$woody+2,xlim=c(0,300)); abline(0,1) #woodies are green
  plot(test$Rd.x,test$Rd,col=out$woody+2,xlim=c(-10,10)); abline(0,1) #woodies are green
    #good; as expected, HB values are less extreme


##Quick summaries
ind.out$region = substr(ind.out$site,1,1)  
sort(tapply(ind.out$Asat,ind.out$species,mean))
sort(tapply(ind.out$Vcmax,ind.out$species,mean))
sort(tapply(ind.out$Jmax,ind.out$species,mean))
sort(tapply(ind.out$alpha,ind.out$species,mean)) #not working, need new approach
sort(tapply(ind.out$Rd,ind.out$species,mean))
plot(tapply(ind.out$Rd,ind.out$species,mean),tapply(ind.out$Asat,ind.out$species,mean),cex=1.5,pch=19)
  
  