#Deriving photosynthetic properties, NSF-IOS project, version 2, all regions
#JDF 4-4-22
#this version uses the plantecophys program to create a comparison dataset of photo params

#### Libraries
library(plantecophys)
library(nlme)
library(lme4)
library(doBy)
library(RCurl)
library(ggplot2)
library(googlesheets4)
library(lmerTest)

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


##########  Plant Ecophys R package analysis, looped over species, A-Ci only 

spp = unique(dat$species)
coef.out = NULL

for(i in 1:length(spp)) {
  df = dat[dat$species==spp[i],]
  df = df[df$PARi>800,]
  df = df[df$Ci>=0,]
  if(spp[i]=="Fallopia japonica") df = df[df$ID!="fajap-F31-1",] #bad curve
  if(spp[i]=="Berberis thunbergii") next #only light curves for this species
  if(spp[i]=="Parthenocissus inserta") next #bad curve
  if(spp[i]=="Laburnum anagyroides") df = df[df$ID!="laana-F35-1",] #bad curve
  if(spp[i]=="Laburnum anagyroides") df = df[df$ID!="laana-F43-2",] #bad curve
  
    
  plot(df$Ci_Pa,df$Photo)
  
  f = fitacis(df,
              varnames = list(ALEAF = "Photo", Tleaf = "Tleaf", Ci = "Ci", PPFD = "PARi", Rd = "Rd", Patm = "Press"),
              Tcorrect=F,
              Patm = mean(dat$Press),
              group = "ID",
              fitmethod = "bilinear" #interestingly, bilinear works and default doesn't
  )
  par(mar=c(5,5,5,5),mfrow=c(1,1))
  plot(f,how="oneplot")
  title(spp[i])
  regout = unique(data.frame(df$ID,df$region,df$species,df$sppcode,df$Date))
  regout = regout[order(regout$df.ID),]
  names(regout) = c("ID","region","species","sppcode","Date")
  coef1 = cbind(coef(f),regout)
  coef.out = rbind(coef.out,coef1)
  #readline()
}  
coef.out

### Partition variance in photo params

covar = read.csv("covariates.csv") #covariates at the species-region level
out = merge(coef.out,covar,all.x=T)
#write.csv(out,"photo_params_ecophys.csv")

#Acer negundo has unclear role in this study: omit
out = out[out$species!="Acer negundo",]

#Vcmax
lme1 = lme(fixed=Vcmax~1,random=~1|family/species/region,data=out)
VarCorr(lme1)
  #largest component of variation is within-region (and within-species); species and family split

lme1 = lme(fixed=Jmax~1,random=~1|family/species/region,data=out)
VarCorr(lme1)
  #no region-level variance; most within-region, lots between species, some between family

lme1 = lme(fixed=Rd~1,random=~1|family/species/region,data=out)
VarCorr(lme1)
  #nearly all residual; some between families

#fixed effects for Vcmax

lme2 = lme(fixed=Vcmax~woody,random=~1|family/species/region,data=out)
summary(lme2)
  #much much less in woodies

lme3 = lmer(Vcmax~woody*region+(1|family/species),data=out)
summary(lme3)
anova(lme3)
  #interaction of woody and region is marginal

#ENA: invader vs. native
lme4 = lmer(Vcmax~invader+(1|family/species),data=out,subset=out$region=="E"&out$woody==1)
summary(lme4)
anova(lme4)
  #greater invader Vcmax compared to native woodies in ENA

lme5 = lmer(Vcmax~invader+(1|family/species),data=out,subset=out$region=="E"&out$woody==0)
summary(lme5)
anova(lme5)
  #not true for herb spp

#France: invader vs. native
lme4 = lmer(Vcmax~invader+(1|family/species),data=out,subset=out$region=="F"&out$woody==1)
summary(lme4)
anova(lme4)
  #no N-I diff for woodies

lme5 = lmer(Vcmax~invader+(1|family/species),data=out,subset=out$region=="F"&out$woody==0)
summary(lme5)
anova(lme5)
  #no N-I diff for herbs

#Japan: invader vs. native
lme4 = lmer(Vcmax~invader+(1|family/species),data=out,subset=out$region=="J"&out$woody==1)
summary(lme4)
anova(lme4)
  #no N-I diff for woodies

lme5 = lmer(Vcmax~invader+(1|family/species),data=out,subset=out$region=="J"&out$woody==0)
summary(lme5)
anova(lme5)
  #big difference for herbs: natives much lower Vcmax than invaders

  #N-I contrast summary:
  #larger Vcmax max for ENA woody invaders compared to natives, and Japan herbaceous invaders compared to natives
  #no differences in France

#home-away contrasts: woodies from Japan to ENA
spp1 = c("Berberis thunbergii","Celastrus orbiculatus","Euonymus alatus","Lonicera japonica","Lonicera morrowii","Rosa multiflora","Viburnum dilatatum")
out1 = out[is.element(out$species,spp1),]
lme6 = lmer(Vcmax~homeaway+(1|family/species),data=out1)
anova(lme6)
summary(lme6)
  #!!Vcmax increases in the away range of J to E invasive shrubs
boxplot(Vcmax~homeaway,data=out1)

#home-away contrasts, ENA trees to France
spp2 = c("Prunus serotina","Quercus rubra","Parthenocissus sp.","Robinia pseudo-acacia","Acer negundo")
out2 = out[is.element(out$species,spp2),]
lme7 = lmer(Vcmax~homeaway+(1|species),data=out2) #all species from different families
anova(lme7)
summary(lme7)
  #no evidence for shift in Vcmax for ENA trees to France
boxplot(Vcmax~homeaway,data=out2)
dotplot(coef(lme7))

#home-away contrasts, ENA herbs to Japan and France
spp3 = c("Ambrosia artemisiifolia","Bidens frondosa","Conyza canadensis","Erigeron annuus","Solidago gigantea")
out3 = out[is.element(out$species,spp3),]
lme7 = lmer(Vcmax~homeaway+(1|species),data=out3) #they are all asters
anova(lme7)
summary(lme7)
  #no evidence for shift in Vcmax for ENA herbs to France and Japan

#home-away contrasts, France herbs to Japan and NY
spp4 = c("Anthoxanthum odoratum","Agrostis stolonifera","Chenopodium album","Daucus carota","Leucanthemum vulgaris","Plantago lanceolata","Artemisia vulgaris")
out4 = out[is.element(out$species,spp4),]
lme7 = lmer(Vcmax~homeaway+(1|family/species),data=out4)
summary(lme7)
anova(lme7)
boxplot(Vcmax~homeaway,data=out4)
dotplot(coef(lme7))
  #!!Vcmax increases in the away range of France to J/E invasive herbaceous spp

  ##summary of home-away Vcmax contrasts:
    #Group A: Japan shrubs to ENA: large increase in away range
    #Group B: ENA herbs to Japan and France: no shift
    #Group C: France herbs to Japan and ENA: large increase in away range
    #Group D: ENA trees to France: no shift

#Story so far with Vcmax (carboxylation)
    #Group A: Japan shrub invaders increase Vcmax in away range in ENA, making it higher than natives
    #Group B: invasive herbs in Japan have higher Vcmax than natives; those from ENA don't shift in away range, but those from France do
    #Group C: no differences between herbs in 

#5 separate scenarios:
out$regionNrange = paste0(out$region,out$nativerange) #combination of sampled region and native range

par(mfrow=c(3,2),mar=c(5,5,1,1))

#1. Woody species in ENA from Japan
invaway = out$Vcmax[out$region=="E"&out$invader=="invader"&out$homeaway=="away"&out$regionNrange=="EJ"&out$woody==1]
invhome = out$Vcmax[out$region=="J"&out$invader=="invader"&out$homeaway=="home"&out$regionNrange=="JJ"&out$woody==1]
nathome = out$Vcmax[out$region=="E"&out$invader=="native"&out$regionNrange=="EE"&out$woody==1]
boxplot(invaway,invhome,nathome,main="Japan->ENA woody",names=c("Inv-Away","Inv-Home","Native"),ylim=c(0,120))

#2. Woody species in France from ENA
invaway = out[out$region=="F"&out$invader=="invader"&out$homeaway=="away"&out$regionNrange=="FE"&out$woody==1,]
invhome = out[out$region=="E"&out$invader=="invader"&out$homeaway=="home"&out$regionNrange=="EE"&out$woody==1,]
nathome = out[out$region=="F"&out$invader=="native"&out$regionNrange=="FF"&out$woody==1,]
boxplot(invaway$Vcmax,invhome$Vcmax,nathome$Vcmax,main="ENA->France woody",names=c("Inv-Away","Inv-Home","Native"),ylim=c(0,120))

#3. Herbaceous species from ENA to France, exclude grasses
invaway = out[out$region=="F"&out$invader=="invader"&out$homeaway=="away"&out$regionNrange=="FE"&out$woody==0,]
invhome = out[out$region=="E"&out$invader=="invader"&out$homeaway=="home"&out$regionNrange=="EE"&out$woody==0,]
nathome = out[out$region=="F"&out$invader=="native"&out$regionNrange=="FF"&out$woody==0&out$grass==0,]
boxplot(invaway$Vcmax,invhome$Vcmax,nathome$Vcmax,main="ENA->France herbaceous",names=c("Inv-Away","Inv-Home","Native"),ylim=c(0,200))

#4. Herbaceous species from ENA to Japan
invaway = out[out$region=="J"&out$invader=="invader"&out$homeaway=="away"&out$regionNrange=="JE"&out$woody==0,]
invhome = out[out$region=="E"&out$invader=="invader"&out$homeaway=="home"&out$regionNrange=="EE"&out$woody==0,]
nathome = out[out$region=="J"&out$invader=="native"&out$regionNrange=="JJ"&out$woody==0&out$grass==0,]
boxplot(invaway$Vcmax,invhome$Vcmax,nathome$Vcmax,main="ENA->Japan herbaceous",names=c("Inv-Away","Inv-Home","Native"),ylim=c(0,200))

#5. Herbaceous species from France to ENA, exclude grasses
invaway = out[out$region=="E"&out$invader=="invader"&out$homeaway=="away"&out$regionNrange=="EF"&out$woody==0,]
invhome = out[out$region=="F"&out$invader=="invader"&out$homeaway=="home"&out$regionNrange=="FF"&out$woody==0&out$grass==0,]
nathome = out[out$region=="E"&out$regionNrange=="EE"&out$woody==0,]
boxplot(invaway$Vcmax,invhome$Vcmax,nathome$Vcmax,main="France->ENA herbaceous",names=c("Inv-Away","Inv-Home","Native"),ylim=c(0,200))
  #missing some key species in ENA: Chenopodium album, Daucus carota; had to use invader home range values for 'nathome'
  #thus comparison a bit of a mess

#6. Herbaceous species from France to Japan
invaway = out[out$region=="J"&out$invader=="invader"&out$homeaway=="away"&out$regionNrange=="JF"&out$woody==0,]
  invaway = invaway[invaway$species!="Chenopodium album",]
invhome = out[out$region=="F"&out$invader=="invader"&out$homeaway=="home"&out$regionNrange=="FF"&out$woody==0,]
  invhome = invhome[invhome$species!="Artemisia vulgaris"&invhome$species!="Chenopodium album"&invhome$species!="Daucus carota",]
nathome = out[out$region=="J"&out$invader=="native"&out$regionNrange=="JJ"&out$woody==0&out$grass==0,]
boxplot(invaway$Vcmax,invhome$Vcmax,nathome$Vcmax,main="France->Japan herbaceous",names=c("Inv-Away","Inv-Home","Native"),ylim=c(0,200))


#next steps:
#plow forward with HB code to get all photo params, refer to this for analysis

