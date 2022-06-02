##summary of photo params fit via HB
#JDF 6-1-22

#note: still need final version of HB fit (with more iterations)

library(lme4)
library(lmerTest)
library(multcomp)
library(lsmeans)

dfold = "/Users/etudiant/Documents/IOS_git/IOS_data/"

hb = read.csv(paste0(dfold,"indoutHB3.csv"))
str(hb)
nls = read.csv(paste0(dfold,"lightcurves_nls_output.csv"))
str(nls)
dat = merge(hb,nls,by.x="ID",by.y="sample",all.x=T)
str(dat)
dat$region = substr(dat$site,1,1)
  #1 is ENA, 2 is France, 3 is Japan
dat$Region = as.factor(dat$region)
#covar = read.csv(paste0(dfold,"covariates.csv")) #covariates at the species-region level

ep = read.csv(paste0(dfold,"photo_params_ecophys.csv")) #includes covars
out = merge(ep,dat,by.x="ID",by.y="ID",all=T)
out = droplevels(out)
out$region = out$region.x
out$species = out$species.x
out$sppcode = out$sppcode.x
out$away = out$homeaway=="away"

wout = out[out$woody==T,]
wout = droplevels(wout)

hout = out[out$woody==F,]
hout = droplevels(hout)


##################################################
#general summary plot
  #parameter, growth form, region (ENA=green, France=blue, Japan=red)

makeTransparent<-function(someColor, alpha=150)
{   newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)}) }

ecol = makeTransparent("darkgreen")
fcol = makeTransparent("blue")
jcol = makeTransparent("red2")
inv = rgb(0,0,0,alpha=0,max=255) #invisible

###Asat
par(mfrow=c(2,1),mar=c(3,3,1,0))
y = wout$Asat
yrange = c(0,60)
plot(wout$sppcode[wout$region=="E"],y[wout$region=="E"],border=ecol,ylim=yrange,las=3,col=inv)
mtext("Asat",side=2,line=2.5,cex=1.5)
par(new=T)
plot(wout$sppcode[wout$region=="F"],y[wout$region=="F"],border=fcol,ylim=yrange,las=3,col=inv)
par(new=T)
plot(wout$sppcode[wout$region=="J"],y[wout$region=="J"],border=jcol,ylim=yrange,las=3,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Woody","topright",bty="n",,cex=.7)
y = hout$Asat
plot(hout$sppcode[hout$region=="E"],y[hout$region=="E"],border=ecol,ylim=yrange,las=3,cex.axis=.8,col=inv)
mtext("Asat",side=2,line=2.5,cex=1.5)
par(new=T)
plot(hout$sppcode[hout$region=="F"],y[hout$region=="F"],border=fcol,ylim=yrange,las=3,cex.axis=.8,col=inv)
par(new=T)
plot(hout$sppcode[hout$region=="J"],y[hout$region=="J"],border=jcol,ylim=yrange,las=3,cex.axis=.8,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Herbaceous","topright",bty="n",,cex=.7)

#home-away contrasts: woodies from Japan to ENA
spp1 = c("Berberis thunbergii","Celastrus orbiculatus","Euonymus alatus","Lonicera japonica","Lonicera morrowii","Rosa multiflora","Viburnum dilatatum")
out1 = out[is.element(out$species,spp1),]
lmer1 = lmer(Asat~homeaway+(1|family/species),data=out1)
anova(lmer1)
summary(lmer1)
  #no change in Asat from home to away
boxplot(Asat~homeaway,data=out1)

#home-away contrasts, ENA trees to France
spp2 = c("Prunus serotina","Quercus rubra","Parthenocissus sp.","Robinia pseudo-acacia","Acer negundo")
out2 = out[is.element(out$species,spp2),]
lmer2 = lmer(Asat~homeaway+(1|species),data=out2) #all species from different families
anova(lmer2)
summary(lmer2)
  #no evidence for shift in Asat for ENA trees to France

#home-away contrasts, ENA herbs to Japan and France
spp3 = c("Ambrosia artemisiifolia","Bidens frondosa","Conyza canadensis","Erigeron annuus","Solidago gigantea")
out3 = out[is.element(out$species,spp3),]
lmer3 = lmer(Asat~homeaway+(1|species),data=out3) #they are all asters
anova(lmer3)
summary(lmer3)
  #no change in Asat

#home-away contrasts, France herbs to Japan and NY
spp4 = c("Anthoxanthum odoratum","Agrostis stolonifera","Chenopodium album","Daucus carota","Leucanthemum vulgaris","Plantago lanceolata","Artemisia vulgaris")
out4 = out[is.element(out$species,spp4),]
lmer4 = lmer(Asat~homeaway+(1|family/species),data=out4)
summary(lmer4)
anova(lmer4)
  #Asat is higher in away range, by 4.5 units (big shift!)
boxplot(Asat~homeaway,data=out4)
boxplot(Asat~region,data=out4)
lmer4b = lmer(Asat~region+(1|family/species),data=out4)
lsmeans(lmer4b,pairwise~region) 
  #sig difference is between France and Japan, although effect size is the same for Japan and ENA difference with France
glht(lmer4b,linfct=mcp(region="Tukey"))

###Vcmax
par(mfrow=c(2,1),mar=c(3,3,1,0))
y = wout$Vcmax.x
yrange = c(0,300)
plot(wout$sppcode[wout$region=="E"],y[wout$region=="E"],border=ecol,ylim=yrange,las=3,col=inv)
mtext("Vcmax",side=2,line=2.5,cex=1.5)
par(new=T)
plot(wout$sppcode[wout$region=="F"],y[wout$region=="F"],border=fcol,ylim=yrange,las=3,col=inv)
par(new=T)
plot(wout$sppcode[wout$region=="J"],y[wout$region=="J"],border=jcol,ylim=yrange,las=3,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Woody","topright",bty="n",,cex=.7)
y = hout$Vcmax.x
plot(hout$sppcode[hout$region=="E"],y[hout$region=="E"],border=ecol,ylim=yrange,las=3,cex.axis=.8,col=inv)
mtext("Vcmax",side=2,line=2.5,cex=1.5)
par(new=T)
plot(hout$sppcode[hout$region=="F"],y[hout$region=="F"],border=fcol,ylim=yrange,las=3,cex.axis=.8,col=inv)
par(new=T)
plot(hout$sppcode[hout$region=="J"],y[hout$region=="J"],border=jcol,ylim=yrange,las=3,cex.axis=.8,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Herbaceous","topright",bty="n",,cex=.7)

#home-away contrasts: woodies from Japan to ENA
spp1 = c("Berberis thunbergii","Celastrus orbiculatus","Euonymus alatus","Lonicera japonica","Lonicera morrowii","Rosa multiflora","Viburnum dilatatum")
out1 = out[is.element(out$species,spp1),]
lmer1 = lmer(Vcmax.x~homeaway+(1|family/species),data=out1)
anova(lmer1)
summary(lmer1)
  #!Vcmax increases in the away range of J to E invasive shrubs, by about 10 units
boxplot(Vcmax.x~homeaway,data=out1)

#home-away contrasts, ENA trees to France
spp2 = c("Prunus serotina","Quercus rubra","Parthenocissus sp.","Robinia pseudo-acacia","Acer negundo")
out2 = out[is.element(out$species,spp2),]
lmer2 = lmer(Vcmax.x~homeaway+(1|species),data=out2) #all species from different families
anova(lmer2)
summary(lmer2)
  #marginal evidence for decline in Vcmax for ENA trees to France (higher away)
boxplot(Vcmax.x~homeaway,data=out2)

#home-away contrasts, ENA herbs to Japan and France
spp3 = c("Ambrosia artemisiifolia","Bidens frondosa","Conyza canadensis","Erigeron annuus","Solidago gigantea")
out3 = out[is.element(out$species,spp3),]
lmer3 = lmer(Vcmax.x~homeaway+(1|species),data=out3) #they are all asters
anova(lmer3)
summary(lmer3)
  #no change in Vcmax

#home-away contrasts, France herbs to Japan and NY
spp4 = c("Anthoxanthum odoratum","Agrostis stolonifera","Chenopodium album","Daucus carota","Leucanthemum vulgaris","Plantago lanceolata","Artemisia vulgaris")
out4 = out[is.element(out$species,spp4),]
lmer4 = lmer(Vcmax.x~homeaway+(1|family/species),data=out4)
summary(lmer4)
anova(lmer4)
  #Vcmax is much higher in away range, by 32 units (big shift!)
boxplot(Vcmax.x~homeaway,data=out4)


###Jmax
par(mfrow=c(2,1),mar=c(3,3,1,0))
y = wout$Jmax.x
yrange = c(0,350)
plot(wout$sppcode[wout$region=="E"],y[wout$region=="E"],border=ecol,ylim=yrange,las=3,col=inv)
mtext("Jmax",side=2,line=2.5,cex=1.5)
par(new=T)
plot(wout$sppcode[wout$region=="F"],y[wout$region=="F"],border=fcol,ylim=yrange,las=3,col=inv)
par(new=T)
plot(wout$sppcode[wout$region=="J"],y[wout$region=="J"],border=jcol,ylim=yrange,las=3,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Woody","topright",bty="n",,cex=.7)
y = hout$Jmax.x
plot(hout$sppcode[hout$region=="E"],y[hout$region=="E"],border=ecol,ylim=yrange,las=3,cex.axis=.8,col=inv)
mtext("Jmax",side=2,line=2.5,cex=1.5)
par(new=T)
plot(hout$sppcode[hout$region=="F"],y[hout$region=="F"],border=fcol,ylim=yrange,las=3,cex.axis=.8,col=inv)
par(new=T)
plot(hout$sppcode[hout$region=="J"],y[hout$region=="J"],border=jcol,ylim=yrange,las=3,cex.axis=.8,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Herbaceous","topright",bty="n",,cex=.7)

#home-away contrasts: woodies from Japan to ENA
spp1 = c("Berberis thunbergii","Celastrus orbiculatus","Euonymus alatus","Lonicera japonica","Lonicera morrowii","Rosa multiflora","Viburnum dilatatum")
out1 = out[is.element(out$species,spp1),]
lmer1 = lmer(Jmax.x~homeaway+(1|family/species),data=out1)
anova(lmer1)
summary(lmer1)
  #no change in Jmax

#home-away contrasts, ENA trees to France
spp2 = c("Prunus serotina","Quercus rubra","Parthenocissus sp.","Robinia pseudo-acacia","Acer negundo")
out2 = out[is.element(out$species,spp2),]
lmer2 = lmer(Jmax.x~homeaway+(1|species),data=out2) #all species from different families
anova(lmer2)
summary(lmer2)
  #no evidence for change in Asat for ENA trees to France (higher away)

#home-away contrasts, ENA herbs to Japan and France
spp3 = c("Ambrosia artemisiifolia","Bidens frondosa","Conyza canadensis","Erigeron annuus","Solidago gigantea")
out3 = out[is.element(out$species,spp3),]
lmer3 = lmer(Jmax.x~homeaway+(1|species),data=out3) #they are all asters
anova(lmer3)
summary(lmer3)
  #no change in Jmax

#home-away contrasts, France herbs to Japan and NY
spp4 = c("Anthoxanthum odoratum","Agrostis stolonifera","Chenopodium album","Daucus carota","Leucanthemum vulgaris","Plantago lanceolata","Artemisia vulgaris")
out4 = out[is.element(out$species,spp4),]
lmer4 = lmer(Jmax.x~homeaway+(1|family/species),data=out4)
summary(lmer4)
anova(lmer4)
  #Jmax is much higher in away range, by 39 units (big shift!)
boxplot(Jmax.x~homeaway,data=out4)

###alpha
par(mfrow=c(2,1),mar=c(3,3,1,0))
y = wout$alpha
yrange = c(0,.8)
plot(wout$sppcode[wout$region=="E"],y[wout$region=="E"],border=ecol,ylim=yrange,las=3,col=inv)
mtext("alpha",side=2,line=2.5,cex=1.5)
par(new=T)
plot(wout$sppcode[wout$region=="F"],y[wout$region=="F"],border=fcol,ylim=yrange,las=3,col=inv)
par(new=T)
plot(wout$sppcode[wout$region=="J"],y[wout$region=="J"],border=jcol,ylim=yrange,las=3,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Woody","topright",bty="n",,cex=.7)
#no herbaceous

#home-away contrasts: woodies from Japan to ENA
spp1 = c("Berberis thunbergii","Celastrus orbiculatus","Euonymus alatus","Lonicera japonica","Lonicera morrowii","Rosa multiflora","Viburnum dilatatum")
out1 = out[is.element(out$species,spp1),]
lmer1 = lmer(alpha~homeaway+(1|family/species),data=out1)
anova(lmer1)
summary(lmer1)
  #alpha increases in the away range of J to E invasive shrubs, by about 0.03
boxplot(alpha~homeaway,data=out1)
boxplot(alpha~region,data=out1)
lsmeans(lmer4b,pairwise~region) 

#home-away contrasts, ENA trees to France
spp2 = c("Prunus serotina","Quercus rubra","Parthenocissus sp.","Robinia pseudo-acacia","Acer negundo")
out2 = out[is.element(out$species,spp2),]
lmer2 = lmer(alpha~homeaway+(1|species),data=out2) #all species from different families
anova(lmer2)
summary(lmer2)
  #much higher alpha in home range for trees to France
boxplot(alpha~homeaway,data=out2)
lmer2b = lmer(alpha~region+(1|species),data=out2) #all species from different families
boxplot(alpha~region,data=out2)

#(no alpha contrasts for herbs)

###Rd##################################
par(mfrow=c(2,1),mar=c(3,3,1,0))
y = wout$Rd.x
yrange = c(0,3)
plot(wout$sppcode[wout$region=="E"],y[wout$region=="E"],border=ecol,ylim=yrange,las=3,col=inv)
mtext("Rd",side=2,line=2.5,cex=1.5)
par(new=T)
plot(wout$sppcode[wout$region=="F"],y[wout$region=="F"],border=fcol,ylim=yrange,las=3,col=inv)
par(new=T)
plot(wout$sppcode[wout$region=="J"],y[wout$region=="J"],border=jcol,ylim=yrange,las=3,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Woody","topright",bty="n",,cex=.7)
y = hout$Rd.x
plot(hout$sppcode[hout$region=="E"],y[hout$region=="E"],border=ecol,ylim=yrange,las=3,cex.axis=.8,col=inv)
mtext("Rd",side=2,line=2.5,cex=1.5)
par(new=T)
plot(hout$sppcode[hout$region=="F"],y[hout$region=="F"],border=fcol,ylim=yrange,las=3,cex.axis=.8,col=inv)
par(new=T)
plot(hout$sppcode[hout$region=="J"],y[hout$region=="J"],border=jcol,ylim=yrange,las=3,cex.axis=.8,col=inv)
legend(legend=c("ENA","France","Japan"),"topleft",bty="n",,cex=.5,text.col=c(ecol,fcol,jcol))
legend(legend="Herbaceous","topright",bty="n",,cex=.7)

#home-away contrasts: woodies from Japan to ENA
spp1 = c("Berberis thunbergii","Celastrus orbiculatus","Euonymus alatus","Lonicera japonica","Lonicera morrowii","Rosa multiflora","Viburnum dilatatum")
out1 = out[is.element(out$species,spp1),]
lmer1 = lmer(Rd.x~homeaway+(1|family/species),data=out1)
anova(lmer1)
summary(lmer1)
  #no change in Rd

#home-away contrasts, ENA trees to France
spp2 = c("Prunus serotina","Quercus rubra","Parthenocissus sp.","Robinia pseudo-acacia","Acer negundo")
out2 = out[is.element(out$species,spp2),]
lmer2 = lmer(Rd.x~homeaway+(1|species),data=out2) #all species from different families
anova(lmer2)
summary(lmer2)
  #no change in Rd

#home-away contrasts, ENA herbs to Japan and France
spp3 = c("Ambrosia artemisiifolia","Bidens frondosa","Conyza canadensis","Erigeron annuus","Solidago gigantea")
out3 = out[is.element(out$species,spp3),]
lmer3 = lmer(Rd.x~homeaway+(1|species),data=out3) #they are all asters
anova(lmer3)
summary(lmer3)
  #no change in Rd

#home-away contrasts, France herbs to Japan and NY
spp4 = c("Anthoxanthum odoratum","Agrostis stolonifera","Chenopodium album","Daucus carota","Leucanthemum vulgaris","Plantago lanceolata","Artemisia vulgaris")
out4 = out[is.element(out$species,spp4),]
lmer4 = lmer(Rd.x~homeaway+(1|family/species),data=out4)
summary(lmer4)
anova(lmer4)
  #Rd is marginally higher in home range, by .25 units
boxplot(Rd.x~homeaway,data=out4)

