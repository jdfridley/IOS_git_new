###########Deriving photosynthetic properties, all regions
#J. Fridley 2024, final version 6-4-24

#### Data directory 
dfold = "C:\\Users\\fridley\\OneDrive - Clemson University\\academic\\projects\\IOS_FranceJapan\\IOS_git2\\IOS_data\\"

#### Libraries
library(R2jags)
library(doBy)
library(RCurl)
library(bayesplot)
library(ggplot2)
library(mcmcplots)
library(gridExtra)

#### Data
master = read.csv("C:\\Users\\fridley\\OneDrive - Clemson University\\academic\\projects\\IOS_FranceJapan\\NSFIOSspreadsheetJDF11.csv")
dat = read.csv(paste0(dfold,"all_licor_data2.csv")) 
str(dat)
dat$region = substr(dat$site,1,1)
  #11592 observations, 95 cols
dat$Date = as.Date(dat$date,format="%m/%d/%y")
dat$sppcode[dat$sppcode=="amrt"] = "amart" #fix typo
dat$ID[dat$filename=="2021-09-13-bifro-dewitt.xlsx"] = "bifro-E37-1" #fix typo

#### Variable inspection

  #unique IDs
  length(unique(dat$ID)) #416
    #some IDs only indicate separate licor files that should be combined under the same ID:
    dat$ID[dat$ID=="ropse-F51-2"] = "ropse-F51-1"
    dat$ID[dat$ID=="acneg-F37-2"] = "acneg-F37-1"
    dat$ID[dat$ID=="acpse-F41-2"] = "acpse-F41-1"
    dat$ID[dat$ID=="paspp-F28-2"] = "paspp-F28-1"
    
    #some dates need to be altered for duplicate samples on different days:
    dat$date[dat$ID=="acneg-F37-1"&dat$date=="6/10/20"] = "6/5/20"
    dat$date[dat$ID=="acpse-F41-1"&dat$date=="7/29/20"] = "7/28/20"

  length(unique(dat$filename)) #421
  length(unique(master$ID)) #433

  #Anet (umol CO2 per m2 per s)
  summary(dat$Photo)
    #omit a handful of extreme readings
    dat = dat[dat$Photo<70,]
    dat = dat[dat$Photo>(-10),]
  
  #Intercelluar CO2 conc (umol CO2 per mol)
  dat = dat[dat$Ci>(-10),]
  dat$Press[is.na(dat$Press)] = 100 #two samples missing Press in ENA; nearly sea level
  
  #convert Ci to Pa units
  dat$Ci_Pa = dat$Ci * dat$Press/ 1000
  
  #O2 partial pressure
  dat$O = 21 * (dat$Press/100)
  
#### Temperature adjusted coefficients
  
  # Constants published in Sharkey et al (2007) Plant Cell Env 30: 1035-1040 
  R=0.008314 #(kJ mol^-1 K^-1)
  dat$Kc=exp(35.9774-80.99/(R*(dat$Tleaf+273.15))) #Michaelis-Menten constant for Rubisco for O2 (Pa)
  dat$Ko=exp(12.3772-23.72/(R*(dat$Tleaf+273.15))) #Michaelis-Menten constant for Rubisco for CO2 (kPa)
  dat$GammaStar=exp(11.187-24.46/(R*(dat$Tleaf+273.15))) #Photorespiration compensation point (Pa)

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
  

#### empirical estimates of AQY for A-q inds, rather than latent variables

  #Variation in AQY is  too small to estimate with the below FvCB approach unless more q measurements
  #were taken in the light-limited region of the curve; AQY is also based on absorbed light as q rather than 
  #q emitted from the LED, so priors based on published data would be inaccurate with an integrating sphere
  #simpler and more accurate approach is to just use linear A-q fit from light-limited region (<100 PPFD)
  
  #by ID, take q values, estimate variance or range
  id.list = unique(dat$ID)
  out = matrix(0,nrow=length(id.list),ncol=16)
  rownames(out) = id.list
  colnames(out) = seq(0,1500,length=16)
  for(i in 1:length(id.list)) {
    ds = dat[dat$ID==id.list[i],]
    out[i,] = table(cut(ds$PARi,seq(-50,1600,by=100)))
  }
  qID = id.list[out[,2]>1]
  q2 = dat[is.element(dat$ID,qID),]
  q2 = q2[q2$PARi<110,] #limit to light-limited portion
  table(q2$species) #21 species total
  
  #loop through curves
  ids = unique(q2$ID)
  alphaout = rep(0,length(ids))
  for(i in 1:length(ids)) {
    y = q2$Photo[q2$ID==ids[i]]
    #if(length(y)<3) next
    x = q2$PARi[q2$ID==ids[i]]
    plot(x,y)
    abline(lsfit(x,y))
    alphaout[i] = lsfit(x,y)$coef[2]
  }
  hist(alphaout)
  mean(alphaout); sd(alphaout)
    #mean = 0.043, sd = 0.0104  
  names(alphaout) = ids
  
  #insert empirical alphas into dat
  dat$alpha.emp = NA
  for(i in 1:length(alphaout)) {
    a = alphaout[i]; id = names(alphaout[i])
    dat$alpha.emp[dat$ID==id] = a
  }
  
############################
######A-Ci params with HB via JAGS
  
  df = dat
  df = df[df$Ci>=0,]
  df = df[df$PARi>950,] #only use data from A-Ci curves
  N = dim(df)[1]
  ind = as.numeric(as.factor(df$ID)) #change to ID for duplicate licor files per ind
  N.indiv = max(ind)
  spp = as.numeric(as.factor(df$species))
  N.spp = max(spp)
  ind.spp = as.numeric(apply(table(ind,spp)>0,1,function(x)names(which(x==T)))) #lists species to which an individual belongs, numerically
  reg = as.numeric(as.factor(df$region)) #1 = ENA 2 = France 3 = Japan
  
  #examine curves
  #plot(df$Ci_Pa[df$region=="E"],df$Photo[df$region=="E"],col="darkgreen")
  #plot(df$Ci_Pa[df$region=="F"],df$Photo[df$region=="F"],col="purple")
  #plot(df$Ci_Pa[df$region=="J"],df$Photo[df$region=="J"],col="orange")
  
  
  mod.photo <- "model
{

    #Priors
    Vcmax.int ~ dnorm(75,0.001)T(0,) #very weak
    Jmax.int ~ dnorm(100,0.001)T(0,) #very weak
    alpha.int ~ dnorm(0.052,100000) #strong prior, low variation across species (Skillman 2008 JEB)

    #for Rd, which can't be negative, treat as fixed (unpooled) effect not random (difficult to estimate with >0 constraint)
    for(i in 1:N.indiv) {
      Rd.int[i] ~ dnorm(0,1)T(0,) #note cannot take on negative values with T(0,)
    }
    
    #individual and species-level variance in photo params
    ind.tau.Vcmax <- ind.sigma.Vcmax^-2 
    ind.sigma.Vcmax ~ dunif(0, 100)
    spp.tau.Vcmax <- spp.sigma.Vcmax^-2 
    spp.sigma.Vcmax ~ dunif(0, 100)
    reg.tau.Vcmax <- reg.sigma.Vcmax^-2 
    reg.sigma.Vcmax ~ dunif(0, 100)
    
    ind.tau.Jmax <- ind.sigma.Jmax^-2 
    ind.sigma.Jmax ~ dunif(0, 100) 
    spp.tau.Jmax <- spp.sigma.Jmax^-2 
    spp.sigma.Jmax ~ dunif(0, 100)
    reg.tau.Jmax <- reg.sigma.Jmax^-2 
    reg.sigma.Jmax ~ dunif(0, 100)

    ind.tau.alpha <- ind.sigma.alpha^-2 
    ind.sigma.alpha ~ dunif(0, 100000) #originally prec of 100
    spp.tau.alpha <- spp.sigma.alpha^-2 
    spp.sigma.alpha ~ dunif(0, 100000) #originally prec of 100
    reg.tau.alpha <- reg.sigma.alpha^-2 
    reg.sigma.alpha ~ dunif(0, 100000) #originally prec of 100

    #residual error
    tau <- sigma^-2 #coverts sd to precision
    sigma ~ dunif(0, 100)  #uniform prior for standard deviation

    for(i in 1:N) { #lowest level (N=N)
        
        Anet[i] ~ dnorm(mu[i],tau) #residual error
        
        mu[i] <- min(mu.v[i],mu.j[i]) - Rd.int[ind[i]] #minimum of RuBP and Rubisco limitation; TPU limitation ignored
        
        # Photo params that have individual pooling nested within species pooling, plus regional effect
        Vcmax[i] <- Vcmax.int + b0.ind.Vcmax[ind[i]] + b0.reg.Vcmax[reg[i]]  #here b0.ind.Vcmax is informed by b0.spp.Vcmax
        Jmax[i] <- Jmax.int + b0.ind.Jmax[ind[i]] + b0.reg.Jmax[reg[i]] #here b0.ind.Jmax is informed by b0.spp.Jmax
        alpha[i] <- alpha.int + b0.ind.alpha[ind[i]] + b0.reg.alpha[reg[i]] #here b0.ind.alpha is informed by b0.spp.alpha

        # ---RuBisCO limited portion---
        mu.v[i] <- (Vcmax[i]*(Ci_Pa[i]-GammaStar[i]))/(Ci_Pa[i]+(Kc[i]*(1+(O[i]/Ko[i]))))
  
        # ---RUBP limited portion---
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
      
    #random intercept for region (E-F-J) effects
      for(i in 1:3) {
        b0.reg.Vcmax[i] ~ dnorm(0,reg.tau.Vcmax)
        b0.reg.Jmax[i] ~ dnorm(0,reg.tau.Jmax)
        b0.reg.alpha[i] ~ dnorm(0,reg.tau.alpha)
      }

}" #end model
  write(mod.photo, "model.txt")
  
  #input lists for JAGS
  params = c("Asat","Vcm.out","Jm.out","alpha.out","Rd.int","sigma","Vcmax.int","Jmax.int","alpha.int","b0.reg.Vcmax","b0.reg.Jmax","b0.reg.alpha") #parameters to monitor
  inits = function() list(Vcmax.int=rnorm(1,100,10),Jmax.int=rnorm(1,75,10),Rd.int=rnorm(N.indiv)) #starting values of fitted parameters
  input = list(N=N,Anet=df$Photo,Ci_Pa=df$Ci_Pa,q=df$PARi,GammaStar=df$GammaStar,Kc=dat$Kc,Ko=dat$Ko,O=df$O,ind=ind,N.indiv=N.indiv,ind.spp=ind.spp,N.spp=N.spp,reg=reg) #input data
  
  #run JAGS model
  jags.p <- jags(model = "model.txt",data = input,inits=inits,param=params,
                 n.chains = 3, #number of separate MCMC chains 3
                 n.iter =10000, #number of iterations per chain; 10000
                 n.thin=5, #thinning 5
                 n.burnin = 2000) #number of initial iterations to discard 2000
  
  jags.p
  
  #examine Rhats (8th col of output)
  hist(jags.p$BUGSoutput$summary[,8])
  summary(jags.p$BUGSoutput$summary[,8])
  #save(jags.p,file="jags.p_5.31.24.RData") #JAGS object saved

  #create output spreadsheets
  detach.jags()
  attach.jags(jags.p)
  
  #ind-level summary
  ind.out = summaryBy(ID~ID+date+site+species+sppcode+alpha.emp,df)
  ind.out = cbind(ind.out,data.frame(Asat=apply(Asat,2,mean),Vcmax=apply(Vcm.out,2,mean),Jmax=apply(Jm.out,2,mean),Rd=apply(Rd.int,2,median),alpha=apply(alpha.out,2,mean),Asat.se=apply(Asat,2,sd), Vcmax.se=apply(Vcm.out,2,sd),Jmax.se=apply(Jm.out,2,sd),Rd.se=apply(Rd.int,2,sd)),alpha.se=apply(alpha.int,2,sd))

  #save as csv file
  #write.csv(ind.out,"indoutHB.csv")


##################################################
###Photosynthesis comparison, natives and invaders (home and away)

  #######calculate predicted A-Ci curves from posteriors
  load("jags.p_5.31.24.RData") #JAGS object saved, jags.p
  
  attach.jags(jags.p)
  
  CiN = 20
  Ci = seq(0,40,length=CiN)
  aci.mat = matrix(0,nrow=dim(Vcm.out)[2],ncol=length(Ci))
  
  #fill A-Ci matrix
  for(i in 1:dim(aci.mat)[1]) {
    #calculate An for given Ci for all permutations, take mean
    for(j in 1:length(Ci)) {
      C = Ci[j]
      A.J = ((Jm.out[,i]*(C-4.275))/((4*C)+(8*4.275)))
      A.V = Vcm.out[,i]*(C-4.275)/(C+(40.49*(1+(21/27.84))))
      An = apply(cbind(A.J,A.V),1,min) - Rd.int[,i]
      aci.mat[i,j] = mean(An)
    }
  }
  colnames(aci.mat) = paste0("Ci.",c(1:CiN))
  
  #ind-level
  ind.out = summaryBy(ID~ID+date+site+species+sppcode+alpha.emp,df)
  ind.out = cbind(ind.out,as.data.frame(aci.mat))
  ind.out$Region = substr(ind.out$site,1,1)
  cov = read.csv("C:\\Users\\fridley\\OneDrive - Clemson University\\academic\\projects\\IOS_FranceJapan\\covariates.csv")
  names(cov)[1:2] = c("species","Region")
  dc = merge(cov,ind.out,by=c("species","Region"))
  
  aci.mat = dc[,15:34]
  
  #take means and ses, A-Ci curves
  aci.inv.away.woody = apply(aci.mat[dc$invader=="invader"&dc$homeaway=="away"&dc$woody==1,],2,function(x)mean(x,na.rm=T))
  aci.nat.home.woody = apply(aci.mat[dc$invader=="native"&dc$homeaway=="home"&dc$woody==1,],2,function(x)mean(x,na.rm=T))
  aci.inv.home.woody = apply(aci.mat[dc$invader=="invader"&dc$homeaway=="home"&dc$woody==1,],2,function(x)mean(x,na.rm=T))
  
  aci.inv.away.woody.sd = apply(aci.mat[dc$invader=="invader"&dc$homeaway=="away"&dc$woody==1,],2,function(x)sd(x,na.rm=T))
  aci.nat.home.woody.sd = apply(aci.mat[dc$invader=="native"&dc$homeaway=="home"&dc$woody==1,],2,function(x)sd(x,na.rm=T))
  aci.inv.home.woody.sd = apply(aci.mat[dc$invader=="invader"&dc$homeaway=="home"&dc$woody==1,],2,function(x)sd(x,na.rm=T))
  
  aci.inv.away.woody.n = apply(aci.mat[dc$invader=="invader"&dc$homeaway=="away"&dc$woody==1,],2,function(x)length(x[!is.na(x)]))
  aci.nat.home.woody.n = apply(aci.mat[dc$invader=="native"&dc$homeaway=="home"&dc$woody==1,],2,function(x)length(x[!is.na(x)]))
  aci.inv.home.woody.n = apply(aci.mat[dc$invader=="invader"&dc$homeaway=="home"&dc$woody==1,],2,function(x)length(x[!is.na(x)]))
  
  aci.inv.away.woody.se = aci.inv.away.woody.sd/sqrt(aci.inv.away.woody.n)
  aci.nat.home.woody.se = aci.nat.home.woody.sd/sqrt(aci.nat.home.woody.n)
  aci.inv.home.woody.se = aci.inv.home.woody.sd/sqrt(aci.inv.home.woody.n)
  
  aci.inv.away.herb = apply(aci.mat[dc$invader=="invader"&dc$homeaway=="away"&dc$woody==0,],2,function(x)mean(x,na.rm=T))
  aci.nat.home.herb = apply(aci.mat[dc$invader=="native"&dc$homeaway=="home"&dc$woody==0,],2,function(x)mean(x,na.rm=T))
  aci.inv.home.herb = apply(aci.mat[dc$invader=="invader"&dc$homeaway=="home"&dc$woody==0,],2,function(x)mean(x,na.rm=T))
  
  aci.inv.away.herb.sd = apply(aci.mat[dc$invader=="invader"&dc$homeaway=="away"&dc$woody==0,],2,function(x)sd(x,na.rm=T))
  aci.nat.home.herb.sd = apply(aci.mat[dc$invader=="native"&dc$homeaway=="home"&dc$woody==0,],2,function(x)sd(x,na.rm=T))
  aci.inv.home.herb.sd = apply(aci.mat[dc$invader=="invader"&dc$homeaway=="home"&dc$woody==0,],2,function(x)sd(x,na.rm=T))
  
  aci.inv.away.herb.n = apply(aci.mat[dc$invader=="invader"&dc$homeaway=="away"&dc$woody==0,],2,function(x)length(x[!is.na(x)]))
  aci.nat.home.herb.n = apply(aci.mat[dc$invader=="native"&dc$homeaway=="home"&dc$woody==0,],2,function(x)length(x[!is.na(x)]))
  aci.inv.home.herb.n = apply(aci.mat[dc$invader=="invader"&dc$homeaway=="home"&dc$woody==0,],2,function(x)length(x[!is.na(x)]))
  
  aci.inv.away.herb.se = aci.inv.away.woody.sd/sqrt(aci.inv.away.woody.n)
  aci.nat.home.herb.se = aci.nat.home.woody.sd/sqrt(aci.nat.home.woody.n)
  aci.inv.home.herb.se = aci.inv.home.woody.sd/sqrt(aci.inv.home.woody.n)
  
  #plot A-Ci curves
  library(RColorBrewer)
  nat.col = brewer.pal(3,"Dark2")[1]
  home.col = brewer.pal(3,"Dark2")[2]
  away.col = brewer.pal(3,"Dark2")[3]
  
  par(mar=c(5,5,1,1))
  plot(Ci*10,aci.inv.away.woody,ylim=c(0,25),type="l",col=away.col,lwd=3,cex.lab=1.4,cex.axis=1.4,
       xlab = expression('C'[i]*' (ppm)'),
       ylab = expression(A[n]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"))
  lines(Ci*10,aci.inv.away.woody+aci.inv.away.woody.se,ylim=c(0,25),type="l",col=away.col,lty=2)
  lines(Ci*10,aci.inv.away.woody-aci.inv.away.woody.se,ylim=c(0,25),type="l",col=away.col,lty=2)
  lines(Ci*10,aci.nat.home.woody,col=nat.col,lwd=3)
  lines(Ci*10,aci.nat.home.woody+aci.nat.home.woody.se,ylim=c(0,25),type="l",col=nat.col,lty=2)
  lines(Ci*10,aci.nat.home.woody-aci.nat.home.woody.se,ylim=c(0,25),type="l",col=nat.col,lty=2)
  lines(Ci*10,aci.inv.home.woody,col=home.col,lwd=3)
  lines(Ci*10,aci.inv.home.woody+aci.inv.home.woody.se,ylim=c(0,25),type="l",col=home.col,lty=2)
  lines(Ci*10,aci.inv.home.woody-aci.inv.home.woody.se,ylim=c(0,25),type="l",col=home.col,lty=2)
  lines(Ci*10,aci.inv.away.herb,col=away.col,lwd=3)
  lines(Ci*10,aci.inv.away.herb+aci.inv.away.herb.se,ylim=c(0,25),type="l",col=away.col,lty=2)
  lines(Ci*10,aci.inv.away.herb-aci.inv.away.herb.se,ylim=c(0,25),type="l",col=away.col,lty=2)
  lines(Ci*10,aci.nat.home.herb,col=nat.col,lwd=3)
  lines(Ci*10,aci.nat.home.herb+aci.nat.home.herb.se,ylim=c(0,25),type="l",col=nat.col,lty=2)
  lines(Ci*10,aci.nat.home.herb-aci.nat.home.herb.se,ylim=c(0,25),type="l",col=nat.col,lty=2)
  lines(Ci*10,aci.inv.home.herb,col=home.col,lwd=3)
  lines(Ci*10,aci.inv.home.herb+aci.inv.home.herb.se,ylim=c(0,25),type="l",col=home.col,lty=2)
  lines(Ci*10,aci.inv.home.herb-aci.inv.home.herb.se,ylim=c(0,25),type="l",col=home.col,lty=2)
  text(240,22.5,"Fields",cex=2)
  text(300,5,"Forests",cex=2)
  legend(x=-50,y=26,c("Invaders (away)","Invaders (home)","Natives"),bty="n",text.col=c(away.col,home.col,nat.col),
         cex=1.2,y.intersp=.75)
  
#####A-q curves
  
cov = read.csv("C:\\Users\\fridley\\OneDrive - Clemson University\\academic\\projects\\IOS_FranceJapan\\covariates.csv")
names(cov)[1:2] = c("species","Region")
hb = read.csv("C:\\Users\\fridley\\OneDrive - Clemson University\\academic\\projects\\IOS_FranceJapan\\IOS_git2\\indoutHB.csv")
hb$Region = substr(hb$site,1,1)
d2 = merge(cov,hb,by=c("species","Region"),all=T)

d2$group = "native"
d2$group[d2$invader=="invader"&d2$homeaway=="away"] = "invaway"
d2$group[d2$invader=="invader"&d2$homeaway=="home"] = "invhome"

d3 = d2[!is.na(d2$alpha.emp),]
boxplot(alpha.emp ~ group, d3,las=2,xlab="")
anova(lm(alpha.emp ~ group,d3))
summary(lm(alpha.emp ~ group,d3))

#A-q response curves; fit theta parameter of nonrectangular hyperbola
q = seq(0,2000,by=100)
qcurve = function(Am,AQY,Rd,q,theta) {
  (1/(2*theta))*(AQY*q+Am-sqrt((AQY*q+Am)^2-4*AQY*theta*Am*q))-Rd
}
plot(q,qcurve(Am=20,AQY=0.052,Rd=.5,theta=.5,q=q),ylim=c(0,25))
lines(q,qcurve(Am=21,AQY=0.06,Rd=.5,theta=.5,q=q))
lines(q,qcurve(Am=20,AQY=0.052,Rd=.5,theta=.9,q=q))

Amg = tapply(d3$Asat,d3$group,mean)
amg = tapply(d3$alpha.emp,d3$group,mean)
rmg = tapply(d3$Rd,d3$group,mean)

plot(q,qcurve(Am=Amg[1],AQY=amg[1],Rd=rmg[1],theta=.5,q=q),ylim=c(0,15))
lines(q,qcurve(Am=Amg[2],AQY=amg[2],Rd=rmg[2],theta=.5,q=q))
lines(q,qcurve(Am=Amg[3],AQY=amg[3],Rd=rmg[3],theta=.5,q=q))

q = seq(0,2000,by=10)
q.mat = matrix(0,ncol=length(q),nrow=dim(d3)[1])

#fill A-q matrix
for(i in 1:dim(q.mat)[1]) {
  alpha = d3$alpha.emp[i]
  Rd = d3$Rd[i]
  Asat = d3$Asat[i]
  for(j in 1:length(q)) {
    q.mat[i,] = qcurve(Am=Asat,AQY=alpha,Rd=Rd,q=q,theta=.5)
  }
}

plot(q,q.mat[1,],ylim=c(-20,45))    
for(i in 2:dim(q.mat)[1]) {
  points(q,q.mat[i,])
}

#take means and ses, A-q curves (woody only)
q.inv.away.woody = apply(q.mat[d3$invader=="invader"&d3$homeaway=="away",],2,function(x)mean(x,na.rm=T))
q.nat.home.woody = apply(q.mat[d3$invader=="native"&d3$homeaway=="home",],2,function(x)mean(x,na.rm=T))
q.inv.home.woody = apply(q.mat[d3$invader=="invader"&d3$homeaway=="home",],2,function(x)mean(x,na.rm=T))

q.inv.away.woody.sd = apply(q.mat[d3$invader=="invader"&d3$homeaway=="away",],2,function(x)sd(x,na.rm=T))
q.nat.home.woody.sd = apply(q.mat[d3$invader=="native"&d3$homeaway=="home",],2,function(x)sd(x,na.rm=T))
q.inv.home.woody.sd = apply(q.mat[d3$invader=="invader"&d3$homeaway=="home",],2,function(x)sd(x,na.rm=T))

q.inv.away.woody.n = apply(q.mat[d3$invader=="invader"&d3$homeaway=="away",],2,function(x)length(x[!is.na(x)]))
q.nat.home.woody.n = apply(q.mat[d3$invader=="native"&d3$homeaway=="home",],2,function(x)length(x[!is.na(x)]))
q.inv.home.woody.n = apply(q.mat[d3$invader=="invader"&d3$homeaway=="home",],2,function(x)length(x[!is.na(x)]))

q.inv.away.woody.se = q.inv.away.woody.sd/sqrt(q.inv.away.woody.n)
q.nat.home.woody.se = q.nat.home.woody.sd/sqrt(q.nat.home.woody.n)
q.inv.home.woody.se = q.inv.home.woody.sd/sqrt(q.inv.home.woody.n)

#save predicted A-q values for woody species for later PNUE analysis
d3$A.q20 = q.mat[,q==20]
d3$A.q50 = q.mat[,q==50]
d3$A.q100 = q.mat[,q==100]
d3$A.q200 = q.mat[,q==200]
#save(d3,file="qmatout.RData")


#plot A-q curves
par(mar=c(5,5,1,1))
plot(q,q.inv.away.woody,ylim=c(0,13),type="l",col=away.col,lwd=3,cex.lab=1.4,cex.axis=1.4,
     xlab = expression("PPFD ("*mu*"mol photons "*m^-2*s^-1*")"),
     ylab = expression(A[n]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"))
lines(q,q.inv.away.woody+q.inv.away.woody.se,ylim=c(0,25),type="l",col=away.col,lty=2)
lines(q,q.inv.away.woody-q.inv.away.woody.se,ylim=c(0,25),type="l",col=away.col,lty=2)
lines(q,q.nat.home.woody,col=nat.col,lwd=3)
lines(q,q.nat.home.woody+q.nat.home.woody.se,ylim=c(0,25),type="l",col=nat.col,lty=2)
lines(q,q.nat.home.woody-q.nat.home.woody.se,ylim=c(0,25),type="l",col=nat.col,lty=2)
lines(q,q.inv.home.woody,col=home.col,lwd=3)
lines(q,q.inv.home.woody+q.inv.home.woody.se,ylim=c(0,25),type="l",col=home.col,lty=2)
lines(q,q.inv.home.woody-q.inv.home.woody.se,ylim=c(0,25),type="l",col=home.col,lty=2)
text(350,12,"Forests",cex=2)
legend(x=2300,y=3,c("Invaders (away)","Invaders (home)","Natives"),bty="n",text.col=c(away.col,home.col,nat.col),
       cex=1.2,y.intersp=.75,xjust=1)

########################
#2-panel plot for publication, A-Ci and A-q

#pdf("Photo_curves.pdf",height=6,width=10)
par(mfrow=c(1,2),mar=c(5,5,1,1),oma=c(1,1,1,1))
plot(Ci*10,aci.inv.away.woody,ylim=c(0,25),type="l",col=away.col,lwd=3,cex.lab=1.4,cex.axis=1.4,
     xlab = expression('C'[i]*' (ppm)'),
     ylab = expression(A[n]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"))
lines(Ci*10,aci.inv.away.woody+aci.inv.away.woody.se,ylim=c(0,25),type="l",col=away.col,lty=2)
lines(Ci*10,aci.inv.away.woody-aci.inv.away.woody.se,ylim=c(0,25),type="l",col=away.col,lty=2)
lines(Ci*10,aci.nat.home.woody,col=nat.col,lwd=3)
lines(Ci*10,aci.nat.home.woody+aci.nat.home.woody.se,ylim=c(0,25),type="l",col=nat.col,lty=2)
lines(Ci*10,aci.nat.home.woody-aci.nat.home.woody.se,ylim=c(0,25),type="l",col=nat.col,lty=2)
lines(Ci*10,aci.inv.home.woody,col=home.col,lwd=3)
lines(Ci*10,aci.inv.home.woody+aci.inv.home.woody.se,ylim=c(0,25),type="l",col=home.col,lty=2)
lines(Ci*10,aci.inv.home.woody-aci.inv.home.woody.se,ylim=c(0,25),type="l",col=home.col,lty=2)
lines(Ci*10,aci.inv.away.herb,col=away.col,lwd=3)
lines(Ci*10,aci.inv.away.herb+aci.inv.away.herb.se,ylim=c(0,25),type="l",col=away.col,lty=2)
lines(Ci*10,aci.inv.away.herb-aci.inv.away.herb.se,ylim=c(0,25),type="l",col=away.col,lty=2)
lines(Ci*10,aci.nat.home.herb,col=nat.col,lwd=3)
lines(Ci*10,aci.nat.home.herb+aci.nat.home.herb.se,ylim=c(0,25),type="l",col=nat.col,lty=2)
lines(Ci*10,aci.nat.home.herb-aci.nat.home.herb.se,ylim=c(0,25),type="l",col=nat.col,lty=2)
lines(Ci*10,aci.inv.home.herb,col=home.col,lwd=3)
lines(Ci*10,aci.inv.home.herb+aci.inv.home.herb.se,ylim=c(0,25),type="l",col=home.col,lty=2)
lines(Ci*10,aci.inv.home.herb-aci.inv.home.herb.se,ylim=c(0,25),type="l",col=home.col,lty=2)
text(210,22.5,"Fields",cex=2)
text(300,5,"Forests",cex=2)
mtext("A)",side=3,line=-2,cex=2.5,xpd=NA,adj=0,outer=T)

plot(q,q.inv.away.woody,ylim=c(0,13),type="l",col=away.col,lwd=3,cex.lab=1.4,cex.axis=1.4,
     xlab = expression("PPFD ("*mu*"mol photons "*m^-2*s^-1*")"),
     ylab = expression(A[n]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"))
lines(q,q.inv.away.woody+q.inv.away.woody.se,ylim=c(0,25),type="l",col=away.col,lty=2)
lines(q,q.inv.away.woody-q.inv.away.woody.se,ylim=c(0,25),type="l",col=away.col,lty=2)
lines(q,q.nat.home.woody,col=nat.col,lwd=3)
lines(q,q.nat.home.woody+q.nat.home.woody.se,ylim=c(0,25),type="l",col=nat.col,lty=2)
lines(q,q.nat.home.woody-q.nat.home.woody.se,ylim=c(0,25),type="l",col=nat.col,lty=2)
lines(q,q.inv.home.woody,col=home.col,lwd=3)
lines(q,q.inv.home.woody+q.inv.home.woody.se,ylim=c(0,25),type="l",col=home.col,lty=2)
lines(q,q.inv.home.woody-q.inv.home.woody.se,ylim=c(0,25),type="l",col=home.col,lty=2)
text(550,11.5,"Forests",cex=2)
legend("bottomright",c("Invaders (away)","Invaders (home)","Natives"),bty="n",text.col=c(away.col,home.col,nat.col),
       cex=1.4,y.intersp=.95,xjust=1)
mtext("B)",side=3,line=-1,cex=2.5,xpd=NA,adj=0,outer=F,at=-700)
#dev.off()


############################################
#####Stats, photosynthetic params

#for all parameters fit by HB model (Vcmax, Jmax, Asat, Rd), query posteriors of differences between three groups
rm(Asat)
attach.jags(jags.p)

dc$group = "native"
dc$group[dc$invader=="invader"&dc$homeaway=="away"] = "invaway"
dc$group[dc$invader=="invader"&dc$homeaway=="home"] = "invhome"
  #order below is invaway, invhome, native

#Fields

Vcm.out1 = Vcm.out[,dc$woody==0]
Jm.out1 = Jm.out[,dc$woody==0]
Asat1 = Asat[,dc$woody==0]
Rd.int1 = Rd.int[,dc$woody==0]
dc1 = dc[dc$woody==0,]

#for each permutation, calculate mean of each group, then examine distribution of group differences
Vcm.mat1 = matrix(0,nrow=dim(Vcm.out1)[1],ncol=3)
Jm.mat1 = matrix(0,nrow=dim(Vcm.out1)[1],ncol=3)
As.mat1 = matrix(0,nrow=dim(Vcm.out1)[1],ncol=3)
Rd.mat1 = matrix(0,nrow=dim(Vcm.out1)[1],ncol=3)
for(i in 1:dim(Vcm.out)[1]) {
  Vcm.mat1[i,] = tapply(Vcm.out1[i,],dc1$group,function(x)mean(x,na.rm=T))
  Jm.mat1[i,] =  tapply(Jm.out1[i,],dc1$group,function(x)mean(x,na.rm=T))
  As.mat1[i,] = tapply(Asat1[i,],dc1$group,function(x)mean(x,na.rm=T))
  Rd.mat1[i,] = tapply(Rd.int1[i,],dc1$group,function(x)mean(x,na.rm=T))
}
Vcm.awayhome = Vcm.mat1[,1] - Vcm.mat1[,2]
Vcm.awaynat = Vcm.mat1[,1] - Vcm.mat1[,3]
Vcm.homenat = Vcm.mat1[,2] - Vcm.mat1[,3]
Jm.awayhome = Jm.mat1[,1] - Jm.mat1[,2]
Jm.awaynat = Jm.mat1[,1] - Jm.mat1[,3]
Jm.homenat = Jm.mat1[,2] - Jm.mat1[,3]
As.awayhome = As.mat1[,1] - As.mat1[,2]
As.awaynat = As.mat1[,1] - As.mat1[,3]
As.homenat = As.mat1[,2] - As.mat1[,3]
Rd.awayhome = Rd.mat1[,1] - Rd.mat1[,2]
Rd.awaynat = Rd.mat1[,1] - Rd.mat1[,3]
Rd.homenat = Rd.mat1[,2] - Rd.mat1[,3]

x = Vcm.awayhome; length(x[x<0])/length(x) #P=0 ***
x = Vcm.awaynat; length(x[x<0])/length(x) #P=0 ***
x = Vcm.homenat; length(x[x<0])/length(x) #P=0.376 NS
x = Jm.awayhome; length(x[x<0])/length(x) #P=0.357 NS
x = Jm.awaynat; length(x[x<0])/length(x) #P=0.004 **
x = Jm.homenat; length(x[x<0])/length(x) #P=0.031 *
x = As.awayhome; length(x[x<0])/length(x) #P=0.106 NS
x = As.awaynat; length(x[x<0])/length(x) #P=0 ***
x = As.homenat; length(x[x<0])/length(x) #P=0.0196 *
x = Rd.awayhome; length(x[x<0])/length(x) #P=0.468 NS
x = Rd.awaynat; length(x[x<0])/length(x) #P=0.262 NS
x = Rd.homenat; length(x[x<0])/length(x) #P=0.291 NS

par(mfrow=c(4,3),mar=c(5,5,1,1),oma=c(1,1,1,1))
hist(Vcm.awayhome); hist(Vcm.awaynat); hist(Vcm.homenat)
hist(Jm.awayhome); hist(Jm.awaynat); hist(Jm.homenat)
hist(As.awayhome); hist(As.awaynat); hist(As.homenat)
hist(Rd.awayhome); hist(Rd.awaynat); hist(Rd.homenat)

#Forests

Vcm.out2 = Vcm.out[,dc$woody==1]
Jm.out2 = Jm.out[,dc$woody==1]
Asat2 = Asat[,dc$woody==1]
Rd.int2 = Rd.int[,dc$woody==1]
dc2 = dc[dc$woody==1,]

#for each permutation, calculate mean of each group, then examine distribution of group differences
Vcm.mat2 = matrix(0,nrow=dim(Vcm.out2)[1],ncol=3)
Jm.mat2 = matrix(0,nrow=dim(Vcm.out2)[1],ncol=3)
As.mat2 = matrix(0,nrow=dim(Vcm.out2)[1],ncol=3)
Rd.mat2 = matrix(0,nrow=dim(Vcm.out2)[1],ncol=3)
for(i in 1:dim(Vcm.out)[1]) {
  Vcm.mat2[i,] = tapply(Vcm.out2[i,],dc2$group,function(x)mean(x,na.rm=T))
  Jm.mat2[i,] =  tapply(Jm.out2[i,],dc2$group,function(x)mean(x,na.rm=T))
  As.mat2[i,] = tapply(Asat2[i,],dc2$group,function(x)mean(x,na.rm=T))
  Rd.mat2[i,] = tapply(Rd.int2[i,],dc2$group,function(x)mean(x,na.rm=T))
}
Vcm.awayhome = Vcm.mat2[,1] - Vcm.mat2[,2]
Vcm.awaynat = Vcm.mat2[,1] - Vcm.mat2[,3]
Vcm.homenat = Vcm.mat2[,2] - Vcm.mat2[,3]
Jm.awayhome = Jm.mat2[,1] - Jm.mat2[,2]
Jm.awaynat = Jm.mat2[,1] - Jm.mat2[,3]
Jm.homenat = Jm.mat2[,2] - Jm.mat2[,3]
As.awayhome = As.mat2[,1] - As.mat2[,2]
As.awaynat = As.mat2[,1] - As.mat2[,3]
As.homenat = As.mat2[,2] - As.mat2[,3]
Rd.awayhome = Rd.mat2[,1] - Rd.mat2[,2]
Rd.awaynat = Rd.mat2[,1] - Rd.mat2[,3]
Rd.homenat = Rd.mat2[,2] - Rd.mat2[,3]

x = Vcm.awayhome; length(x[x<0])/length(x) #P=0.082, marginal
x = Vcm.awaynat; length(x[x<0])/length(x) #P=0.197 NS
x = Vcm.homenat; length(x[x<0])/length(x) #P=0.524 NS
x = Jm.awayhome; length(x[x<0])/length(x) #P=0.105 NS
x = Jm.awaynat; length(x[x<0])/length(x) #P=0.782 NS
x = Jm.homenat; length(x[x<0])/length(x) #P=0.907 NS; ie .093, marginal
x = As.awayhome; length(x[x<0])/length(x) #P=0.093 NS
x = As.awaynat; length(x[x<0])/length(x) #P=0.456 NS
x = As.homenat; length(x[x<0])/length(x) #P=0.717 NS
x = Rd.awayhome; length(x[x<0])/length(x) #P=0.868 NS
x = Rd.awaynat; length(x[x<0])/length(x) #P=0.783 NS
x = Rd.homenat; length(x[x<0])/length(x) #P=0.441 NS

par(mfrow=c(4,3),mar=c(5,5,1,1),oma=c(1,1,1,1))
hist(Vcm.awayhome); hist(Vcm.awaynat); hist(Vcm.homenat)
hist(Jm.awayhome); hist(Jm.awaynat); hist(Jm.homenat)
hist(As.awayhome); hist(As.awaynat); hist(As.homenat)
hist(Rd.awayhome); hist(Rd.awaynat); hist(Rd.homenat)

#for AQY, compare empirical fits
library(lme4)
library(lmerTest)
library(emmeans)
lme1 = lmer(alpha.emp ~ group + (1|Region) + (1|sppcode),d3)
summary(lme1)
emmeans(lme1, list(pairwise ~ group), adjust = "tukey")
  #away - native comparison: away inv higher by .012, t=-4.1, P<0.001
  #away - home: NS (P=.29)
  #home - nat: home higher by .0085; t=2.863, P=0.0189

AQY.mean = tapply(d3$alpha.emp,d3$group,function(x)mean(x,na.rm=T))
AQY.q80l = tapply(d3$alpha.emp,d3$group,function(x)quantile(x,.1))
AQY.q80u = tapply(d3$alpha.emp,d3$group,function(x)quantile(x,.9))
AQY.q95l = tapply(d3$alpha.emp,d3$group,function(x)quantile(x,.025))
AQY.q95u = tapply(d3$alpha.emp,d3$group,function(x)quantile(x,.975))


##############################
#Plot photo param comparisons

#pdf("Photo_params_fig_rev.pdf",height=6,width=10)
par(mfcol=c(2,5),mar=c(.5,3,1,.5),oma=c(1,10,5,1))
box1 = boxplot(Vcm.mat1[,1],Vcm.mat1[,2],Vcm.mat1[,3],xlab="",ylab="",pch=19,col=c(away.col,home.col,nat.col),
                 cex=2,cex.axis=1.3,ylim=c(70,160),xaxt="n")
  mtext(expression(V['c,max']),cex=1.7,line=1,xpd=F)
  mtext("Fields",side=2,cex=1.8,las=2,outer=T,at=.75,line=.75)
  mtext("Forests",side=2,cex=1.8,las=2,outer=T,at=.28,line=.75)
  text(1:3,box1$stats[5,]+8,c("a","b","b"),col=c(away.col,home.col,nat.col),cex=1.7)  
box1 = boxplot(Vcm.mat2[,1],Vcm.mat2[,2],Vcm.mat2[,3],xlab="",ylab="",pch=19,col=c(away.col,home.col,nat.col),
          cex=2,cex.axis=1.3,ylim=c(30,130),xaxt="n")  
  text(1:3,box1$stats[5,]+8,c("a","a","a"),col=c(away.col,home.col,nat.col),cex=1.7)  
box1 = boxplot(Jm.mat1[,1],Jm.mat1[,2],Jm.mat1[,3],xlab="",ylab="",pch=19,col=c(away.col,home.col,nat.col),
                 cex=2,cex.axis=1.3,ylim=c(90,170),xaxt="n")
  mtext(expression(J['max']),cex=1.7,line=1,xpd=F)
  text(1:3,box1$stats[5,]+8,c("a","a","b"),col=c(away.col,home.col,nat.col),cex=1.7)  
box1 = boxplot(Jm.mat2[,1],Jm.mat2[,2],Jm.mat2[,3],xlab="",ylab="",pch=19,col=c(away.col,home.col,nat.col),
                 cex=2,cex.axis=1.3,ylim=c(30,130),xaxt="n",outline=F)  
  text(1:3,box1$stats[5,]+8,c("ab","a","b"),col=c(away.col,home.col,nat.col),cex=1.7) 
box1 = boxplot(As.mat1[,1],As.mat1[,2],As.mat1[,3],xlab="",ylab="",pch=19,col=c(away.col,home.col,nat.col),
                 cex=2,cex.axis=1.3,ylim=c(17,30),xaxt="n")
  mtext(expression(A['sat']),cex=1.7,line=1,xpd=F)
  text(1:3,box1$stats[5,]+1,c("a","a","b"),col=c(away.col,home.col,nat.col),cex=1.7)  
box1 = boxplot(As.mat2[,1],As.mat2[,2],As.mat2[,3],xlab="",ylab="",pch=19,col=c(away.col,home.col,nat.col),
                 cex=2,cex.axis=1.3,ylim=c(5,20),xaxt="n",outline=F)  
  text(1:3,box1$stats[5,]+1,c("a","a","a"),col=c(away.col,home.col,nat.col),cex=1.7) 
box1 = boxplot(Rd.mat1[,1],Rd.mat1[,2],Rd.mat1[,3],xlab="",ylab="",pch=19,col=c(away.col,home.col,nat.col),
                 cex=2,cex.axis=1.3,ylim=c(.6,1.1),xaxt="n",outline=F)
  mtext(expression(R['d']),cex=1.7,line=1,xpd=F)
  text(1:3,box1$stats[5,]+.05,c("a","a","a"),col=c(away.col,home.col,nat.col),cex=1.7)  
box1 = boxplot(Rd.mat2[,1],Rd.mat2[,2],Rd.mat2[,3],xlab="",ylab="",pch=19,col=c(away.col,home.col,nat.col),
                 cex=2,cex.axis=1.3,ylim=c(.5,1.2),xaxt="n",outline=F)  
  text(1:3,box1$stats[5,]+.05,c("a","a","a"),col=c(away.col,home.col,nat.col),cex=1.7) 
plot(1:30,1:30,col="white",xaxt="n",yaxt="n",bty="n",ylim=c(0,30),xlim=c(0,30))  
  text(12,19,"Invaders (away)",cex=1.9,col=away.col,xpd=NA)
  text(12,15,"Invaders (home)",cex=1.9,col=home.col,xpd=NA)
  text(12,11,"Natives",cex=1.9,col=nat.col,xpd=NA)
box1 = boxplot(alpha.emp~group,data=d3,xlab="",ylab="",pch=19,col=c(away.col,home.col,nat.col),
                 cex=2,cex.axis=1.3,ylim=c(.017,.08),xaxt="n",outline=F)
  mtext("AQY",cex=1.7,line=1,xpd=F)
  text(1:3,box1$stats[5,]+.005,c("a","a","b"),col=c(away.col,home.col,nat.col),cex=1.7) 
#dev.off()  
  

##Alternative plotting of photo params, showing distribution of sample rather than posterior mean of each group
  #for parallel method to fig 3 (empirical data)
  
  pdf("Photo_params_fig_empirical.pdf",height=9,width=6)
  
  library(RColorBrewer)
  nat.col = brewer.pal(3,"Dark2")[1]
  home.col = brewer.pal(3,"Dark2")[2]
  away.col = brewer.pal(3,"Dark2")[3]
  boxcols = c(away.col,home.col,nat.col,away.col,home.col,nat.col)
  labs = rep(c("Inv-Away","Inv-Home","Native"),2)
  
  
  par(mfcol=c(3,2),mar=c(.5,2,1,2.5),oma=c(5,3,2,1))
  box1 = boxplot(Vcmax~group+woody,dat,at=c(1,2,3,4.5,5.5,6.5),col=boxcols,ylab="",
        names=rep("",6),xlab="",xaxt="n",xpd=NA,ylim=c(30,220))
  abline(v=3.75,col="gray")
  mtext(expression(V[cmax]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2)
  mtext(c("Fields","Forests"),side=3,at=c(2,5.5),cex=1.4)
  text(6,660,"Habitat*",cex=1.5)
  text(c(1,2,3,4.5,5.5,6.5),box1$stats[5,]+20,c("a","b","b","a","a","a"),col=boxcols,cex=1.5)  
  box1 = boxplot(Jmax~group+woody,dat,at=c(1,2,3,4.5,5.5,6.5),col=boxcols,ylab="",
                 names=rep("",6),xlab="",xaxt="n",xpd=NA,ylim=c(0,250))
  abline(v=3.75,col="gray")
  mtext(expression(J[max]*" ("*mu*"mol "*e^-" "*m^-2*s^-1*")"),side=2,line=2)
  text(6,660,"Habitat*",cex=1.5)
  text(c(1,2,3,4.5,5.5,6.5),box1$stats[5,]+20,c("a","a","b","ab","a","b"),col=boxcols,cex=1.5)  
  box1 = boxplot(Asat~group+woody,dat,at=c(1,2,3,4.5,5.5,6.5),col=boxcols,ylab="",
                 names=rep("",6),xlab="",xaxt="n",xpd=NA,ylim=c(0,50))
  abline(v=3.75,col="gray")
  text(x,-3,labs,srt=45,xpd=NA,adj=1,col=boxcols,cex=1.5)
  mtext(expression(A[sat]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2)
  text(c(1,2,3,4.5,5.5,6.5),box1$stats[5,]+2,c("a","a","b","a","a","a"),col=boxcols,cex=1.5)  
  box1 = boxplot(Rd~group+woody,dat,at=c(1,2,3,4.5,5.5,6.5),col=boxcols,ylab="",
                 names=rep("",6),xlab="",xaxt="n",xpd=NA,ylim=c(0,2))
  abline(v=3.75,col="gray")
  mtext(expression(R[d]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2)
  mtext(c("Fields","Forests"),side=3,at=c(2,5.5),cex=1.4)
  text(c(1,2,3,4.5,5.5,6.5),box1$stats[5,]+.2,c("a","a","a","a","a","a"),col=boxcols,cex=1.5)  
  box1 = boxplot(alpha.emp~group,d3,at=c(4.5,5.5,6.5),col=boxcols,ylab="",
                 names=rep("",3),xlab="",xaxt="n",xpd=NA,xlim=c(.5,7),ylim=c(.02,.07))
  abline(v=3.75,col="gray")
  text(c(1,2,3,4.5,5.5,6.5),box1$stats[5,]+.005,c("","","","a","a","b"),col=boxcols,cex=1.5)  
  mtext(expression(AQY*" (mol "*e^-" "*mol*" "*quanta^-")"),side=2,line=2)
  #mtext(expression(AQY*" (mol "*e^-" "*mol*" "*quanta^-*")"),side=2,line=2)
  text(x,0.015,labs,srt=45,xpd=NA,adj=1,col=boxcols,cex=1.5)
  #plot(1:30,1:30,col="white",xaxt="n",yaxt="n",bty="n",ylim=c(0,30),xlim=c(0,30))  
  #text(12,19,"Invaders (away)",cex=1.9,col=away.col,xpd=NA)
  #text(12,15,"Invaders (home)",cex=1.9,col=home.col,xpd=NA)
  #text(12,11,"Natives",cex=1.9,col=nat.col,xpd=NA)
  dev.off()
  
  
  #pdf("Photo_params_fig_rev2.pdf",height=8,width=8)
    #compare to q80 and q95 of posteriors
  x = c(1,2,3,4.5,5.5,6.5)
  labs = rep(c("Inv-Away","Inv-Home","Native"),2)
  par(mfcol=c(2,2),mar=c(.5,2,1,2.5),oma=c(5,3,2,1))
  m1 = colMeans(Vcm.mat1)
  q951 = apply(Vcm.mat1,2,function(x)quantile(x,c(.025,.975)))
  q801 = apply(Vcm.mat1,2,function(x)quantile(x,c(.1,.9)))
  m2 = colMeans(Vcm.mat2)
  q952 = apply(Vcm.mat2,2,function(x)quantile(x,c(.025,.975)))
  q802 = apply(Vcm.mat2,2,function(x)quantile(x,c(.1,.9)))
  plot(x,c(m1,m2),xlab="",xaxt="n",ylab="",
       col=boxcols,cex=2,pch=19,xlim=c(.5,7),ylim=c(50,145))
  arrows(x,c(q951[1,],q952[1,]),x,c(q951[2,],q952[2,]),length=0,code=3,angle=90,col=boxcols,lwd=2)
  arrows(x,c(q801[1,],q802[1,]),x,c(q801[2,],q802[2,]),length=0,code=3,angle=90,col=boxcols,lwd=4)
  abline(v=3.75,col="gray")
  mtext(expression(V[cmax]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2)
  mtext(c("Fields","Forests"),side=3,at=c(2,5.5),cex=1.4,line=.7)
  text(c(1,2,3,4.5,5.5,6.5),c(q951[2,],q952[2,])+5,c("a","b","b","a","a","a"),col=boxcols,cex=1.2)  
  
  m1 = colMeans(Jm.mat1)
  q951 = apply(Jm.mat1,2,function(x)quantile(x,c(.025,.975)))
  q801 = apply(Jm.mat1,2,function(x)quantile(x,c(.1,.9)))
  m2 = colMeans(Jm.mat2)
  q952 = apply(Jm.mat2,2,function(x)quantile(x,c(.025,.975)))
  q802 = apply(Jm.mat2,2,function(x)quantile(x,c(.1,.9)))
  plot(x,c(m1,m2),ylim=c(45,155),xlab="",xaxt="n",ylab="",
       col=boxcols,cex=2,pch=19,xlim=c(.5,7))
  arrows(x,c(q951[1,],q952[1,]),x,c(q951[2,],q952[2,]),length=0,code=3,angle=90,col=boxcols,lwd=2)
  arrows(x,c(q801[1,],q802[1,]),x,c(q801[2,],q802[2,]),length=0,code=3,angle=90,col=boxcols,lwd=4)
  abline(v=3.75,col="gray")
  mtext(expression(J[max]*" ("*mu*"mol "*e^-" "*m^-2*s^-1*")"),side=2,line=2)
  text(c(1,2,3,4.5,5.5,6.5),c(q951[2,],q952[2,])+5,c("a","a","b","ab","a","b"),col=boxcols,cex=1.2)  
  text(x,37,labs,srt=45,xpd=NA,adj=1,col=boxcols,cex=1.5)
  
  m1 = colMeans(As.mat1)
  q951 = apply(As.mat1,2,function(x)quantile(x,c(.025,.975)))
  q801 = apply(As.mat1,2,function(x)quantile(x,c(.1,.9)))
  m2 = colMeans(As.mat2)
  q952 = apply(As.mat2,2,function(x)quantile(x,c(.025,.975)))
  q802 = apply(As.mat2,2,function(x)quantile(x,c(.1,.9)))
  plot(x,c(m1,m2),ylim=c(8,28),xlab="",xaxt="n",ylab="",
       col=boxcols,cex=2,pch=19,xlim=c(.5,7))
  arrows(x,c(q951[1,],q952[1,]),x,c(q951[2,],q952[2,]),length=0,code=3,angle=90,col=boxcols,lwd=2)
  arrows(x,c(q801[1,],q802[1,]),x,c(q801[2,],q802[2,]),length=0,code=3,angle=90,col=boxcols,lwd=4)
  abline(v=3.75,col="gray")
  mtext(expression(A[sat]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2)
  mtext(c("Fields","Forests"),side=3,at=c(2,5.5),cex=1.4,line=.7)
  text(c(1,2,3,4.5,5.5,6.5),c(q951[2,],q952[2,])+1,c("a","a","b","a","a","a"),col=boxcols,cex=1.2)  
  
  m1 = colMeans(Rd.mat1)
  q951 = apply(Rd.mat1,2,function(x)quantile(x,c(.025,.975)))
  q801 = apply(Rd.mat1,2,function(x)quantile(x,c(.1,.9)))
  m2 = colMeans(Rd.mat2)
  q952 = apply(Rd.mat2,2,function(x)quantile(x,c(.025,.975)))
  q802 = apply(Rd.mat2,2,function(x)quantile(x,c(.1,.9)))
  plot(x,c(m1,m2),ylim=c(.6,1.1),xlab="",xaxt="n",ylab="",
       col=boxcols,cex=2,pch=19,xlim=c(.5,7))
  arrows(x,c(q951[1,],q952[1,]),x,c(q951[2,],q952[2,]),length=0,code=3,angle=90,col=boxcols,lwd=2)
  arrows(x,c(q801[1,],q802[1,]),x,c(q801[2,],q802[2,]),length=0,code=3,angle=90,col=boxcols,lwd=4)
  abline(v=3.75,col="gray")
  mtext(expression(R[d]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2)
  text(c(1,2,3,4.5,5.5,6.5),c(q951[2,],q952[2,])+.02,c("a","a","a","a","a","a"),col=boxcols,cex=1.2)  
  text(x,.56,labs,srt=45,xpd=NA,adj=1,col=boxcols,cex=1.5)
  
  #dev.off()
  