#Deriving photosynthetic properties, NSF-IOS project, version 1, all regions
#JDF 3-30-22

#to do list
#1 reliable output for all individs, A-Ci curves
#3 summary posteriors of major contrasts (home-away, native-nonnative, etc)
  #fixed effect (eg of region) or random (with ind nested within?)
#4 integrate other data (N sources, light, etc)

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
dat = read.csv("/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/all_licor_data.csv")
str(dat)
dat$region = substr(dat$site,1,1)
  #11592 observations, 95 cols
dat$Date = as.Date(dat$date,format="%m/%d/%y")
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
  
  
#### Example A-Ci analysis 1

#loop over species
  spp = unique(dat$species)
  par(mfrow=c(3,3),oma=c(1.5,1.5,3,1),mar=c(3,3,0,0))
  
  for(i in 2:5) {  
  eg = dat[dat$species==spp[i],]
  table(eg$filename)
  #eg = eg[eg$filename==levels(as.factor(eg$filename))[1],]
  eg = eg[eg$PARi>1010,] #only saturaing light
    
  #par(mar=c(3,3,0,0),oma=c(1.5,1.5,3,1))
  plot(eg$Ci_Pa,eg$Photo,ylab="", xlab="",cex.lab=1.2,cex.axis=1.5,cex=2)
  mtext(expression("Intercellular "*CO[2]*" Pressure (Pa)"),side=1,line=3.3,cex=1.5)
  mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)
  points(eg$Ci_Pa,eg$Photo,col=as.numeric(as.factor(eg$filename)),pch=19,cex=2)
  mtext(spp[i],line=-20,at=100,cex=2)
  
  # ---RuBisCO limited portion---
  #(Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko)))))-Rd
  
  # ---RUBP limited portion---
  #((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar)))-Rd
  
  # Simultaneous estimation method described by Dubois et al. 2007 New Phyt 176:402-414
  # Could change optimization algorithm (default here is Gauss-Newton)
  # Could also do a "grid search" if estimates are sensitive to starting values
  
    #Version 1: nls (all individuals fit together, no REs), using Mason's code
  
    aci.fit <- nls(Photo~ifelse(((Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko)))))<((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar))),((Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko))))),((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar))))-Rd,start=list(Vcmax=50,J=100,Rd=0.5),data=eg) 
  
    summary(aci.fit)
  
    Vcmax<-summary(aci.fit)$coef[1,1]
    J<-summary(aci.fit)$coef[2,1]
    Rd<-summary(aci.fit)$coef[3,1]
    nls.param = c(Vcmax,J,Rd)
  
    par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
    plot(eg$Ci_Pa,eg$Photo,ylab="", xlab="",cex.lab=1.2,cex.axis=1.5,cex=2)
    mtext(expression("Intercellular "*CO[2]*" Pressure (Pa)"),side=1,line=3.3,cex=1.5)
    mtext(expression("Net photosynthetic rate (umol "* CO[2]* m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)
    curve(ifelse(((Vcmax*(x+mean(eg$GammaStar)))/(x+(mean(eg$Kc)*(1+(O/mean(eg$Ko))))))<((J*(x+mean(eg$GammaStar)))/((4*x)+(8*mean(eg$GammaStar)))),((Vcmax*(x+mean(eg$GammaStar)))/(x+(mean(eg$Kc)*(1+(O/mean(eg$Ko)))))),((J*(x+mean(eg$GammaStar)))/((4*x)+(8*mean(eg$GammaStar)))))-Rd,add=T)

    #Version 2: include random intercepts for individuals using nlme
    
    eg$O = O
    aci.fit2<-nlme(model=Photo~ifelse(((Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko)))))<((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar))),((Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko))))),((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar))))-Rd,start=c(Vcmax=nls.param[1],J=nls.param[2],Rd=nls.param[3]),data=eg,groups=~filename,fixed=Vcmax+J+Rd~1,random=Vcmax+J+Rd~1,control=list(msMaxIter=1000)) 
  
    summary(aci.fit2)  
  
    Vcmax<-aci.fit2$coef$fixed[1]
    J<-aci.fit2$coef$fixed[2]
    Rd<-aci.fit2$coef$fixed[3]
    
    curve(ifelse(((Vcmax*(x+mean(eg$GammaStar)))/(x+(mean(eg$Kc)*(1+(O/mean(eg$Ko))))))<((J*(x+mean(eg$GammaStar)))/((4*x)+(8*mean(eg$GammaStar)))),((Vcmax*(x+mean(eg$GammaStar)))/(x+(mean(eg$Kc)*(1+(O/mean(eg$Ko)))))),((J*(x+mean(eg$GammaStar)))/((4*x)+(8*mean(eg$GammaStar)))))-Rd,add=T,col="red")
 
    par(add=F)
}    
  
  #many curves aren't working
  #variance requires light levels: need to estimate params individually; use nlme predict?
  
  #fit all data simultaneously       

    aci.fit2 <- nlme(
    model=Photo~ifelse(((Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko)))))<((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar))),((Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko))))),((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar))))-Rd,
    start=c(Vcmax=50,J=100,Rd=.5),
    data=dat[200:400,],
    #groups=~species,
    fixed=Vcmax+J+Rd~1,
    random=Vcmax+J+Rd~1|species/filename,
    control=list(msMaxIter=1000)  )

  summary(aci.fit2) 
  
  #breaks easily
  
  ##Version 3: ecophys package: works well; note BETH is light curves only

  spp = unique(dat$species)
  for(i in 1:length(spp)) {
  df = dat[dat$species==spp[i],]
  df = df[df$PARi>1010,]
  df = df[df$Ci>=0,]
  
  if(spp[i]=="Berberis thunbergii") next #only light curves for this species
  
  f = fitacis(df,
             varnames = list(ALEAF = "Photo", Tleaf = "Tleaf", Ci = "Ci", PPFD = "PARi", Rd = "Rd", Patm = "Press"),
             Tcorrect=F,
             Patm = mean(dat$Press),
             group = "filename",
             fitmethod = "bilinear" #interestingly, bilinear works and default doesn't
            )
  par(mar=c(5,5,5,5),mfrow=c(1,1))
  plot(f,how="oneplot")
  title(spp[i])
  summary(f)
  coef(f)
  readline()
}  
    
  ##Version 4: mixed effect version with nlmer (not working)

  df = dat[dat$species=="Chenopodium album",]
  df = df[df$PARi>1010,]
  df = df[df$Ci>=0,]
  
  start.df = c(Vcmax=50,J=100,Rd=.5)
  
  aci.fit3 <- nlmer(Photo~(ifelse(((Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko)))))<((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar))),((Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko))))),((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar))))-Rd ) ~ 0 
                    + Vcmax+J+Rd + 
                      (0 + Vcmax+J+Rd | filename),
            start=start.df,
            data=df, nAGQ=0,
            verbose=T)
    #not working
  
  #### Version 5: HB via JAGS, target species only
  
  ##forest species: A-q and A-Ci curves
  spp.forest
  
  spp = "Rosa multiflora"
  df = dat[dat$species==spp,]
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
                n.iter =3000, #number of iterations per chain
                n.thin=3, #thinning
                n.burnin = 500) #number of initial iterations to discard
  
  jags.p

  #proper caterpillar plotting
  attach.jags(jags.p)
  quartz()
  par(mfrow=c(1,5),mar=c(5,5,3,1))
  #Rd
  mat95 = apply(Rd,2,function(x)quantile(x,c(.025,.975)))
  mat80 = apply(Rd,2,function(x)quantile(x,c(.1,.9)))
  med = apply(Rd,2,median)
  range = c(min(mat95),max(mat95))
  plot(med,c(1:length(med)),xlim=range,pch=19,cex=1.5,ylab="",col=cols,xlab="Rd",cex.axis=1.3,cex.lab=1.5)
  segments(mat95[1,],c(1:length(med)),mat95[2,],c(1:length(med)),col=cols)
  segments(mat80[1,],c(1:length(med)),mat80[2,],c(1:length(med)),lwd=2,col=cols)
  #Vcmax
  mat95 = apply(Vcm.out,2,function(x)quantile(x,c(.025,.975)))
  mat80 = apply(Vcm.out,2,function(x)quantile(x,c(.1,.9)))
  med = apply(Vcm.out,2,median)
  range = c(min(mat95),max(mat95))
  plot(med,c(1:length(med)),xlim=range,pch=19,cex=1.5,ylab="",col=cols,xlab="Vcmax",cex.axis=1.3,cex.lab=1.5)
  segments(mat95[1,],c(1:length(med)),mat95[2,],c(1:length(med)),col=cols)
  segments(mat80[1,],c(1:length(med)),mat80[2,],c(1:length(med)),lwd=2,col=cols)
  #Jmax
  mat95 = apply(Jm.out,2,function(x)quantile(x,c(.025,.975)))
  mat80 = apply(Jm.out,2,function(x)quantile(x,c(.1,.9)))
  med = apply(Jm.out,2,median)
  range = c(min(mat95),max(mat95))
  plot(med,c(1:length(med)),xlim=range,pch=19,cex=1.5,ylab="",col=cols,xlab="Jmax",cex.axis=1.3,cex.lab=1.5)
  segments(mat95[1,],c(1:length(med)),mat95[2,],c(1:length(med)),col=cols)
  segments(mat80[1,],c(1:length(med)),mat80[2,],c(1:length(med)),lwd=2,col=cols)
  #alpha
  mat95 = apply(alpha.out,2,function(x)quantile(x,c(.025,.975)))
  mat80 = apply(alpha.out,2,function(x)quantile(x,c(.1,.9)))
  med = apply(alpha.out,2,median)
  range = c(min(mat95),max(mat95))
  plot(med,c(1:length(med)),xlim=range,pch=19,cex=1.5,ylab="",col=cols,xlab="Alpha",cex.axis=1.3,cex.lab=1.5)
  segments(mat95[1,],c(1:length(med)),mat95[2,],c(1:length(med)),col=cols)
  segments(mat80[1,],c(1:length(med)),mat80[2,],c(1:length(med)),lwd=2,col=cols)
  #Asat
  mat95 = apply(Asat,2,function(x)quantile(x,c(.025,.975)))
  mat80 = apply(Asat,2,function(x)quantile(x,c(.1,.9)))
  med = apply(Asat,2,median)
  range = c(min(mat95),max(mat95))
  plot(med,c(1:length(med)),xlim=range,pch=19,cex=1.5,ylab="",col=cols,xlab="Asat",cex.axis=1.3,cex.lab=1.5)
  segments(mat95[1,],c(1:length(med)),mat95[2,],c(1:length(med)),col=cols)
  segments(mat80[1,],c(1:length(med)),mat80[2,],c(1:length(med)),lwd=2,col=cols)
  
  
  
  #compare to nls
  aci.fit <- nls(Photo~ifelse(((Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko)))))<((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar))),((Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko))))),((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar))))-Rd,start=list(Vcmax=50,J=100,Rd=0.5),data=df) 
  summary(aci.fit)
  
  #compare to plantecophys
  f = fitacis(df,
              varnames = list(ALEAF = "Photo", Tleaf = "Tleaf", Ci = "Ci", PPFD = "PARi", Rd = "Rd", Patm = "Press"),
              Tcorrect=F,
              Patm = mean(dat$Press),
              group = "ID",
              fitmethod = "bilinear" #interestingly, bilinear works and default doesn't
  )
  par(mar=c(5,5,5,5),mfrow=c(1,1))
  plot(f,how="oneplot"); title(spp)
  coef(f)
  
  #compare results
  attach.jags(jags.p)
  par(mar=c(.3,.3,.1,.1),mfrow=c(3,N.indiv))
  for(i in 1:N.indiv) {
    hist(Vcm.out[,i],main=paste0("Vcmax",i),xlim=c(40,180)); abline(v=coef(f)[i,2],col="red") }
  for(i in 1:N.indiv) {
    hist(Jm.out[,i],main=paste0("Jmax",i),xlim=c(60,240)); abline(v=coef(f)[i,3],col="red") }
  for(i in 1:N.indiv) {
    hist(Rd[,i],main=paste0("Rd",i),xlim=c(0,2)); abline(v=coef(f)[i,4],col="red") }
  
  #plot comparison of HB (pooled) results with plantecophys (unpooled)
  dfP = df
  dfP$Vcmax = apply(Vcm.out,2,median)[ind]
  dfP$Jmax = apply(Jm.out,2,median)[ind]
  dfP$Rd = apply(Rd,2,median)[ind]
  
  pred.func = function(Vcmax,J,Rd,Ci_Pa,GammaStar,Ko,Kc,O) ifelse(((Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko)))))<((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar))),((Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko))))),((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar))))-Rd
  #need to change by adding alpha etc
  predY = pred.func(dfP$Vcmax,dfP$Jmax,dfP$Rd,dfP$Ci_Pa,dfP$GammaStar,dfP$Ko,dfP$Kc,dfP$O)  

  quartz()
  par(mar=c(5,5,5,5),mfrow=c(1,1))
  plot(f,how="oneplot"); title(spp)
  points(dfP$Ci,predY,col=as.numeric(as.factor(df$region)),pch=19,cex=1)
  

  
  ##other graphing
  out = as.mcmc(jags.p)
  #Plot with caterplot
  par(mfrow=c(1,3))
  caterplot(out,regex="Vcm.out",col=cols,reorder=F)
  caterplot(out,regex="Jm.out",col=cols,reorder=F)
  caterplot(out,regex="Rd",col=cols,reorder=F)
  
  #Plot posteriors via bayesplot  
  posterior = as.matrix(out)
  par(mfrow=c(3,1),mar=c(5,5,1,1))
  color_scheme_set("red")
  p1 = mcmc_areas(posterior,
                  regex_pars = c("Vcm.out"),
                  prob=.8,
                  rhat=c(1,1,1,1,1,2,2,2))
  p2 = mcmc_areas(posterior,
                  regex_pars = c("Jm.out"),
                  prob=.8)
  p3 = mcmc_areas(posterior,
                  regex_pars = c("Rd"),
                  prob=.8)
  grid.arrange(p1,p2,p3,nrow=1)
  
  
  

  ##Version 6: HB via JAGS, all species
  
  #spp1 = "Bidens frondosa"
  #spp2 = "Chenopodium album"
  #spp3 = "Conyza canadensis"
  #df = dat[dat$species==spp1|dat$species==spp2|dat$species==spp3,]
  #df = dat[is.element(dat$species,unique(dat$species)[1:4]),]
  df = dat
  df = df[df$PARi>1010,]
  df = df[df$Ci>=0,]
  N = dim(df)[1]
  ind = as.numeric(as.factor(df$filename)) #grouping vector (individual)
  N.indiv = max(ind)
  spp = as.numeric(as.factor(df$species))
  N.spp = max(spp)
  ind.spp = as.numeric(apply(table(ind,spp)>0,1,function(x)names(which(x==T)))) #lists species to which an individual belongs, numerically
  
  mod.photo <- "model
{
    #Priors
    Vcmax.int ~ dnorm(75,0.001) #very weak
    Jmax.int ~ dnorm(100,0.001) #very weak
    
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

    #residual error
    tau <- sigma^-2 #coverts sd to precision
    sigma ~ dunif(0, 100)  #uniform prior for standard deviation

    #NOTE: add site-level RE
    #NOTE need to add light limitation component for forest species

    for(i in 1:N) { #lowest level (N=N)
        
        Anet[i] ~ dnorm(mu[i],tau) #residual error
        
        mu[i] <- min(mu.v[i],mu.j[i]) - Rd.int[ind[i]] #minimum of RuBP and Rubisco limitation; TPU limitation ignored
        
        # Photo params that have individual pooling nested within species pooling
        Vcmax[i] <- Vcmax.int + b0.ind.Vcmax[ind[i]]  #here b0.ind.Vcmax is informed by b0.spp.Vcmax
        Jmax[i] <- Jmax.int + b0.ind.Jmax[ind[i]] #here b0.ind.Jmax is informed by b0.spp.Jmax

        # ---RuBisCO limited portion---
        mu.v[i] <- (Vcmax[i]*(Ci_Pa[i]-GammaStar[i]))/(Ci_Pa[i]+(Kc[i]*(1+(O[i]/Ko[i]))))
  
        # ---RUBP limited portion---
        mu.j[i] <- ((Jmax[i]*(Ci_Pa[i]-GammaStar[i]))/((4*Ci_Pa[i])+(8*GammaStar[i])))

    }

	  #random intercept for individual effects, nested within species
      for(i in 1:N.indiv) {
        b0.ind.Vcmax[i] ~ dnorm(b0.spp.Vcmax[ind.spp[i]],ind.tau.Vcmax) #individual value sampled from species' distribution
        b0.ind.Jmax[i] ~ dnorm(b0.spp.Jmax[ind.spp[i]],ind.tau.Jmax) #individual value sampled from species' distribution
        
        #output monitoring
        Vcm.out[i] <- Vcmax.int + b0.ind.Vcmax[i] #includes spp constraint
        Jm.out[i] <- Jmax.int + b0.ind.Jmax[i]
      }  

    #random intercept for species effects
      for(i in 1:N.spp) {
        b0.spp.Vcmax[i] ~ dnorm(0,spp.tau.Vcmax)
        b0.spp.Jmax[i] ~ dnorm(0,spp.tau.Jmax)
        
        #monitor spp-level values
        Vcm.spp[i] <- Vcmax.int + b0.spp.Vcmax[i]
        Jm.spp[i] <- Jmax.int + b0.spp.Jmax[i]

      }

}" #end model
  write(mod.photo, "model.txt")
  
  #input lists for JAGS
  params = c("Vcm.out","Jm.out","Rd.int","sigma","Vcm.spp","Jm.spp","Vcmax.int","Jmax.int") #parameters to monitor
  inits = function() list(Vcmax.int=rnorm(1),Jmax.int=rnorm(1),Rd.int=rlnorm(N.indiv)) #starting values of fitted parameters
  input = list(N=N,Anet=df$Photo,Ci_Pa=df$Ci_Pa,GammaStar=df$GammaStar,Kc=dat$Kc,Ko=dat$Ko,O=df$O,ind=ind,N.indiv=N.indiv,ind.spp=ind.spp,N.spp=N.spp) #input data
  
  #run JAGS model
  jags.p <- jags(model = "model.txt",data = input,inits=inits,param=params,
                 n.chains = 3, #number of separate MCMC chains
                 n.iter =3000, #number of iterations per chain
                 n.thin=3, #thinning
                 n.burnin = 500) #number of initial iterations to discard
  
  jags.p
  
  
  
  #compare to nls
  aci.fit <- nls(Photo~ifelse(((Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko)))))<((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar))),((Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko))))),((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar))))-Rd,start=list(Vcmax=50,J=100,Rd=0.5),data=df) 
  summary(aci.fit)
  
  #compare to plantecophys
  f = fitacis(df,
              varnames = list(ALEAF = "Photo", Tleaf = "Tleaf", Ci = "Ci", PPFD = "PARi", Rd = "Rd", Patm = "Press"),
              Tcorrect=F,
              Patm = mean(dat$Press),
              group = "filename",
              fitmethod = "bilinear" #interestingly, bilinear works and default doesn't
  )
  par(mar=c(5,5,5,5),mfrow=c(1,1))
  plot(f,how="oneplot")
  coef(f)
  
  #compare results: IGNORED, comparison below
  #attach.jags(jags.p)
  #par(mar=c(3,3,1,1),mfrow=c(3,N.indiv))
  #for(i in 1:N.indiv) {
  #  hist(Vcm.out[,i],main=paste0("Vcmax",i),xlim=c(40,180)); abline(v=coef(f)[i,2],col="red") }
  #for(i in 1:N.indiv) {
  #  hist(Jm.out[,i],main=paste0("Jmax",i),xlim=c(60,240)); abline(v=coef(f)[i,3],col="red") }
  #for(i in 1:N.indiv) {
  #  hist(Rd[,i],main=paste0("Rd",i),xlim=c(0,2)); abline(v=coef(f)[i,4],col="red") }
  
  #plot comparison of HB (pooled) results with plantecophy (unpooled)
  dfP = df
  dfP$Vcmax = apply(Vcm.out,2,median)[ind]
  dfP$Jmax = apply(Jm.out,2,median)[ind]
  dfP$Rd = apply(Rd,2,median)[ind]
  
  pred.func = function(Vcmax,J,Rd,Ci_Pa,GammaStar,Ko,Kc,O) ifelse(((Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko)))))<((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar))),((Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko))))),((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar))))-Rd
  predY = pred.func(dfP$Vcmax,dfP$Jmax,dfP$Rd,dfP$Ci_Pa,dfP$GammaStar,dfP$Ko,dfP$Kc,dfP$O)  
  
  quartz()
  par(mar=c(5,5,5,5),mfrow=c(1,1))
  plot(f,how="oneplot")
  points(dfP$Ci,predY,col="blue",pch=19,cex=1)
  
  #create output spreadsheets
  detach.jags()
  attach(jags.p,overwrite=T)
  
  #species-level: note summaries of Rd are based on the median
  spp.out = data.frame(species=levels(as.factor(df$species)),Vcmax=apply(Vcm.spp,2,mean),Jmax=apply(Jm.spp,2,mean),Rd=tapply(apply(Rd.int,2,median),ind.spp,mean),
                       Vcmax.se=apply(Vcm.spp,2,sd),Jmax.se=apply(Jm.spp,2,sd),Rd.se=tapply(apply(Rd.int,2,median),ind.spp,sd))
  woody = c(0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,0,1,0) #7 woodies
  invader = c(1,1,1,1,0,1,0,1,0,1,0,1,0,0,0,0,1,0,0,0,1,0,1,0,0,0)
  plot(spp.out$Vcmax,spp.out$Jmax,col=woody+1,pch=19,cex=1.5) #woodies have lowest Jmax, Vcmax (except Cryptotaenia japonica) values
  plot(spp.out$Vcmax,spp.out$Jmax,col=invader+1,pch=19,cex=1.5) #invaders are higher 
  
  #ind-level
  ind.out = summaryBy(filename~filename+date+site+species+sppcode,df)
  ind.out = ind.out[,-6]
  ind.out = cbind(ind.out,data.frame(Vcmax=apply(Vcm.out,2,mean),Jmax=apply(Jm.out,2,mean),Rd=apply(Rd.int,2,median),
                       Vcmax.se=apply(Vcm.out,2,sd),Jmax.se=apply(Jm.out,2,sd),Rd.se=apply(Rd.int,2,sd)))
  #if want to add plantecophys params for comparison
  fout = coef(f); names(fout) = paste0(names(fout),".f")
  ind.out2 = cbind(ind.out,fout)
  
  #plot comparison values
  quartz()
  par(mar=c(5,5,2,2),mfrow=c(2,2))
  plot(ind.out2$Vcmax,ind.out2$Vcmax.f,col=as.numeric(as.factor(ind.out2$species)),pch=19,cex=1.5); abline(0,1,lwd=2,lty=2)
  arrows(ind.out2$Vcmax,ind.out2$Vcmax.f-ind.out2$Vcmax_SE.f,ind.out2$Vcmax,ind.out2$Vcmax.f+ind.out2$Vcmax_SE.f,length=.1,angle=90,code=3,col="gray")
  arrows(ind.out2$Vcmax-ind.out2$Vcmax.se,ind.out2$Vcmax.f,ind.out2$Vcmax+ind.out2$Vcmax.se,ind.out2$Vcmax.f,length=.1,angle=90,code=3,col="gray")
  points(ind.out2$Vcmax,ind.out2$Vcmax.f,col=as.numeric(as.factor(ind.out2$species)),pch=19,cex=1.5); abline(0,1,lwd=2,lty=2)
  plot(ind.out2$Jmax,ind.out2$Jmax.f,col=as.numeric(as.factor(ind.out2$species)),pch=19,cex=1.5); abline(0,1,lwd=2,lty=2)
  arrows(ind.out2$Jmax,ind.out2$Jmax.f-ind.out2$Jmax_SE.f,ind.out2$Jmax,ind.out2$Jmax.f+ind.out2$Jmax_SE.f,length=.1,angle=90,code=3,col="gray")
  arrows(ind.out2$Jmax-ind.out2$Jmax.se,ind.out2$Jmax.f,ind.out2$Jmax+ind.out2$Jmax.se,ind.out2$Jmax.f,length=.1,angle=90,code=3,col="gray")
  points(ind.out2$Jmax,ind.out2$Jmax.f,col=as.numeric(as.factor(ind.out2$species)),pch=19,cex=1.5); abline(0,1,lwd=2,lty=2)
  plot(ind.out2$Rd,ind.out2$Rd.f,col=as.numeric(as.factor(ind.out2$species)),pch=19,cex=1.5); abline(0,1,lwd=2,lty=2)
  arrows(ind.out2$Rd,ind.out2$Rd.f-ind.out2$Rd_SE.f,ind.out2$Rd,ind.out2$Rd.f+ind.out2$Rd_SE.f,length=.1,angle=90,code=3,col="gray")
  arrows(ind.out2$Rd-ind.out2$Rd.se,ind.out2$Rd.f,ind.out2$Rd+ind.out2$Rd.se,ind.out2$Rd.f,length=.1,angle=90,code=3,col="gray")
  points(ind.out2$Rd,ind.out2$Rd.f,col=as.numeric(as.factor(ind.out2$species)),pch=19,cex=1.5); abline(0,1,lwd=2,lty=2)
  