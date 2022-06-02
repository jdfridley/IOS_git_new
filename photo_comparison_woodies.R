#look for artifacts across regions in licor data
#JDF 6-2-22

library(nlstools)

dfold = "/Users/etudiant/Documents/IOS_git/IOS_data/"
dat = read.csv(paste0(dfold,"all_licor_data2.csv"))
#version 2 is the same as v1, with one licor file (Fallopia F7) added
str(dat)
dat$region = substr(dat$site,1,1)
#11592 observations, 95 cols
dat$Date = as.Date(dat$date,format="%m/%d/%y")
dat$sppcode[dat$sppcode=="amrt"] = "amart" #fix typo
dat$ID[dat$filename=="2021-09-13-bifro-dewitt.xlsx"] = "bifro-E37-1" #fix typo

#convert Ci to Pa units if needed
dat$Ci_Pa = dat$Ci * dat$Press/ 1000

dat2 = dat[dat$CO2R>360&dat$CO2R<440,]

spp.select = "Prunus serotina"
df = dat2[dat2$species==spp.select,]
str(df)

prse.alpha = c(.205,.224,.204,.275,.272,.405,.389,.279,.329,.332)

#plot all light curves separately, indicate whether E or F
sample = as.character(na.exclude(unique(as.character(df$ID))))
i=1
par(mfcol=c(5,2),mar=c(3,3,1,1))
for(i in 1:length(sample)) {
  dfs = df[df$ID==sample[i],]
  plot(dfs$PARi,dfs$Photo,ylim=c(-2,13),xlim=c(0,1600))
  abline(h=0,lty=2,col="gray")
  title(main=sample[i])
  
  #if(i==5) next
  
  #if(sample[i]=="paspp-E1-1") {
  #  dfs = dfs[(dfs$PARi>800&dfs$Photo<10)==F,]
  #}
  
  #Mason's A-q code 
  #PARlrc<-dfs$PARi #PAR (aka PPFD or Q)
  #photolrc<-dfs$Photo #net photosynthetic rate (Anet)
  curvelrc<-na.exclude(data.frame(dfs$PARi,dfs$Photo))
  PARlrc <- curvelrc$dfs.PARi
  photolrc <- curvelrc$dfs.Photo

  
  #fit curve
  curve.nlslrc = nls(photolrc ~ (1/(2*theta))*(AQY*PARlrc+Am-sqrt((AQY*PARlrc+Am)^2-4*AQY*theta*Am*PARlrc))-Rd,start=list(Am=(max(photolrc)-min(photolrc)),AQY=0.05,Rd=-min(photolrc),theta=1),control=list(warnOnly=T,tol=.1)) 
  summary(curve.nlslrc) #summary of model fit; coef 1 is Amax, 2 is AQY, 3 is Rd, 4 is theta
  
  res = nlsResiduals(curve.nlslrc) #extracts std residuals from nls fit
  stdres = res$resi2[,2]  #standardized residuals
  
  #refit without std res > 2
  thresh = 1
  curve.nlslrc = nls(photolrc ~ (1/(2*theta))*(AQY*PARlrc+Am-sqrt((AQY*PARlrc+Am)^2-4*AQY*theta*Am*PARlrc))-Rd,start=list(Am=(max(photolrc)-min(photolrc)),AQY=0.05,Rd=-min(photolrc),theta=1),control=list(warnOnly=T,tol=.1),subset=abs(stdres)<thresh) 
  curve((1/(2*summary(curve.nlslrc)$coef[4,1]))*(summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1]-sqrt((summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1])^2-4*summary(curve.nlslrc)$coef[2,1]*summary(curve.nlslrc)$coef[4,1]*summary(curve.nlslrc)$coef[1,1]*x))-summary(curve.nlslrc)$coef[3,1],lwd=2,col="blue",add=T)
  mtext(paste0("alpha=",prse.alpha[i]),side=1,line=-2,at=1400,cex=.7)
}

