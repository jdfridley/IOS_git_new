#NSF project, light curves from nls, for Julie's analysis

dat = read.csv("all_licor_data2.csv")
#version 2 is the same as v1, with one licor file (Fallopia F7) added
str(dat)
dat$region = substr(dat$site,1,1)
dat$Date = as.Date(dat$date,format="%m/%d/%y")
dat$sppcode[dat$sppcode=="amrt"] = "amart" #fix typo

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
dat2 = dat[dat$CO2R>360&dat$CO2R<440,]
dat2 = dat2[is.element(dat2$species,spp.forest),]
str(dat2)
sample = unique(dat2$ID)

#empty vectors for output
AQY = NULL
Amax = NULL
Rd = NULL
theta = NULL

#loop over each sample
for(i in 1:length(sample)) {

  #skip bad curves (should be investigated, some have only a few poor values)
  bad = c(14,20,24,36,55,61,68,79,88,127,130,131,152)
  if(is.element(i,bad)) next
    
  df = dat2[dat2$ID==sample[i],]

  #Mason's A-q code 
  PARlrc<-df$PARi #PAR (aka PPFD or Q)
  photolrc<-df$Photo #net photosynthetic rate (Anet)
  curvelrc<-data.frame(PARlrc,photolrc)

  #fit curve
  curve.nlslrc = nls(photolrc ~ (1/(2*theta))*(AQY*PARlrc+Am-sqrt((AQY*PARlrc+Am)^2-4*AQY*theta*Am*PARlrc))-Rd,start=list(Am=(max(photolrc)-min(photolrc)),AQY=0.05,Rd=-min(photolrc),theta=1),control=list(warnOnly=T,tol=.1)) 
  summary(curve.nlslrc) #summary of model fit; coef 1 is Amax, 2 is AQY, 3 is Rd, 4 is theta

  #check fit
  par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
  plot(PARlrc,photolrc,xlab="", ylab="", ylim=c(-2,max(photolrc)+2),cex.lab=1.2,cex.axis=1.5,cex=2)
  mtext(expression("PPFD ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=2)
  mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2,cex=2)
  curve((1/(2*summary(curve.nlslrc)$coef[4,1]))*(summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1]-sqrt((summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1])^2-4*summary(curve.nlslrc)$coef[2,1]*summary(curve.nlslrc)$coef[4,1]*summary(curve.nlslrc)$coef[1,1]*x))-summary(curve.nlslrc)$coef[3,1],lwd=2,col="blue",add=T)
  readline() #pause until any input
  
  #save params
  Amax[i] = coef(curve.nlslrc)[1]
  AQY[i] = coef(curve.nlslrc)[2]
  Rd[i] = coef(curve.nlslrc)[3]
  theta[i] = coef(curve.nlslrc)[4]
}

#combine into data frame
out = data.frame(sample,Amax,AQY,Rd,theta)

#note some bad model fits, should exclude:
out = out[out$AQY<1,]
out = out[out$Amax<50,]
