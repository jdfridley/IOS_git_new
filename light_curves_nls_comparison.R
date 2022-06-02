#NSF project, light curves from nls, for Julie's analysis

library(nlstools)

dfold = "/Users/etudiant/Documents/IOS_git/IOS_data/"

dat = read.csv(paste0(dfold,"all_licor_data2.csv"))
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
  bad = c(14,20,24,36,55,61,4,7)
    #14 is acneg-F37-1, no light curve
    #20 is laana-F35-1, no light curve
    #24 is laana-F43-1, ditto
    #36 is paspp-F28-1, no PAR < 800
    #55 is acpse-F41-1, no light curve
    #61 is ropse-F51-2, no light curve
    #4 is qurob-F18-1, bad light curve
    #7 is qurub-F31-1, bad light curve
    if(is.element(i,bad)) next
    
  df = dat2[dat2$ID==sample[i],]

  #some curves need manual point omission
  if(sample[i]=="acneg-E10-1") {
    df = df[(df$PARi>800&df$Photo<10)==F,]
  }
  
  if(sample[i]=="bethu-E7-1") {
    df = df[(df$PARi>800&df$Photo<6)==F,]
  }
  
  if(sample[i]=="liben-E3-1") {
    df = df[(df$PARi>800&df$Photo<6)==F,]
  }
  
  if(sample[i]=="romul-E11-1") {
    df = df[(df$PARi>800&df$Photo<10)==F,]
  }
  
  if(sample[i]=="bethu-E6-1") {
    df = df[(df$PARi>800&df$Photo<10)==F,]
  }
  
  if(sample[i]=="paspp-E1-1") {
    df = df[(df$PARi>800&df$Photo<10)==F,]
  }
  
  if(sample[i]=="ceorb-J22-1") {
    df = df[c(1,3,5,7,9,11,13,15,17,19,21),]
  }
  
  if(sample[i]=="ropse-E17-1") {
    df = df[(df$PARi>800&df$Photo>9)==F,]
  }
  
  if(sample[i]=="qurub-F33-1") {
    df = df[(df$PARi>1400)==F,]
  }
  
  if(sample[i]=="laana-F10-1") {
    df = df[(df$Photo>4)==F,]
  }
  
  if(sample[i]=="euala-E3-1") {
    df = df[(df$PARi>200&df$Photo<1)==F,]
  }
  
  if(sample[i]=="romul-E3-1") {
    df = df[(df$PARi>800&df$Photo<4)==F,]
  }
  
  if(sample[i]=="liben-E3-1") {
    df = df[(df$PARi>800&df$Photo<6)==F,]
  }
  
  if(sample[i]=="ropse-E20-1") {
    df = df[(df$Photo>13)==F,]
  }
  
  #Mason's A-q code 
  PARlrc<-df$PARi #PAR (aka PPFD or Q)
  photolrc<-df$Photo #net photosynthetic rate (Anet)
  curvelrc<-data.frame(PARlrc,photolrc)

  #fit curve
  curve.nlslrc = nls(photolrc ~ (1/(2*theta))*(AQY*PARlrc+Am-sqrt((AQY*PARlrc+Am)^2-4*AQY*theta*Am*PARlrc))-Rd,start=list(Am=(max(photolrc)-min(photolrc)),AQY=0.05,Rd=-min(photolrc),theta=1),control=list(warnOnly=T,tol=.1)) 
  summary(curve.nlslrc) #summary of model fit; coef 1 is Amax, 2 is AQY, 3 is Rd, 4 is theta

  res = nlsResiduals(curve.nlslrc) #extracts std residuals from nls fit
  stdres = res$resi2[,2]  #standardized residuals
  
  #refit without std res > 2
  thresh = 2
  thresh2 = c(32,64,65,66,67,69,70,71,73,74,85,86,100,122,125,128,132,158,165,167)
  if(is.element(i,thresh2)) thresh = 1 #some curves need lower threshold for std res
  skip = c(31,151) #a few curves break with fewer points
  if(is.element(i,skip)==F) {
    curve.nlslrc = nls(photolrc ~ (1/(2*theta))*(AQY*PARlrc+Am-sqrt((AQY*PARlrc+Am)^2-4*AQY*theta*Am*PARlrc))-Rd,start=list(Am=(max(photolrc)-min(photolrc)),AQY=0.05,Rd=-min(photolrc),theta=1),control=list(warnOnly=T,tol=.1),subset=abs(stdres)<thresh) 
  }
  
  #cut based on time
  #groups = as.numeric(cut(df$FTime,2)) #two groups
  #true = groups == median(groups) #most common group
  
  #refit without smaller group
  #curve.nlslrc = nls(photolrc ~ (1/(2*theta))*(AQY*PARlrc+Am-sqrt((AQY*PARlrc+Am)^2-4*AQY*theta*Am*PARlrc))-Rd,start=list(Am=(max(photolrc)-min(photolrc)),AQY=0.05,Rd=-min(photolrc),theta=1),control=list(warnOnly=T,tol=.1),subset=true)

  #check fit
  par(mar=c(3,3,3,0))
  plot(PARlrc,photolrc,xlab="", ylab="", ylim=c(-2,max(photolrc)+2),cex.lab=1.2,cex.axis=1.5,cex=2)
  title(main=sample[i],line=2)
  #points(PARlrc,photolrc,col=true,pch=19,cex=1.5)
  points(PARlrc,photolrc,col=abs(stdres)>thresh,pch=19,cex=1.5)
  mtext(expression("PPFD ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=2)
  mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2,cex=2)
  curve((1/(2*summary(curve.nlslrc)$coef[4,1]))*(summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1]-sqrt((summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1])^2-4*summary(curve.nlslrc)$coef[2,1]*summary(curve.nlslrc)$coef[4,1]*summary(curve.nlslrc)$coef[1,1]*x))-summary(curve.nlslrc)$coef[3,1],lwd=2,col="blue",add=T)
  text(100,0,i)
  #readline() #pause until any input
  
  #save params
  Amax[i] = coef(curve.nlslrc)[1]
  AQY[i] = coef(curve.nlslrc)[2]
  Rd[i] = coef(curve.nlslrc)[3]
  theta[i] = coef(curve.nlslrc)[4]
}

#combine into data frame
out = data.frame(sample,Amax,AQY,Rd,theta)

#examine, make a list of bad samples
hist(out$Amax) #good
hist(out$AQY) #a few a bit high, but ok
hist(out$Rd) #it is what it is

plot(out$Amax,out$AQY)
plot(out$AQY,out$Rd)


bad.samples = as.character(out$sample[is.na(out$Amax)]) #initial set
summary(out$Amax)
out$Amax[out$Amax>40]
out$sample[out$Amax>40]
which(out$Amax>40)
bad.samples = c(bad.samples,as.character(out$sample[out$Amax>40&!is.na(out$Amax)]))

summary(out$AQY)
out$sample[out$AQY>.35]
which(out$AQY>.35)
bad.samples = unique(c(bad.samples,as.character(out$sample[out$AQY>1&!is.na(out$AQY)])))

summary(out$Rd)
out$sample[out$Rd>5]
  #no new bad samples

#dataset is good: export
#write.csv(out,"lightcurves_nls_output.csv")


