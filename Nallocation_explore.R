#examine N partitioning data, exploration 1
#JDF 5-20-22

library(googlesheets4)

#current master
master = read_sheet("https://docs.google.com/spreadsheets/d/1e_a7xrEdfgMd5gFor2pkEwf-wVwuLapQJF8Xzp512EE/edit#gid=0")
dat = as.data.frame(master)
str(dat)
names(dat) = make.names(names(dat))
names(dat)

#water and SDS protein
x = as.numeric(unlist(dat$Water.soluble.Protein[dat$Region=="ENA"]))
plot(density(x),col="red2")
y = na.omit(as.numeric(unlist(dat$Water.soluble.Protein[dat$Region=="France"])))
lines(density(y),col="blue")
  #ENA is higher, particulary see Fallopia


#N allocation, different components, note protein is 16% N

#1. Rubisco; reported microg per cm2

dat$Rub.g = as.numeric((dat$Rubisco.content)) * as.numeric(as.character(dat$Specific.leaf.area))  #ug Rubisco per g leaf
dat$RubN = dat$Rub.g / 16 #ug N in Rubisco per g leaf
dat$Rub.perc = (dat$RubN / 1000000) / (as.numeric(as.character(dat$leafN))/100) #Rubisco g per g leaf N
hist(dat$Rub.perc)
summary(dat$Rub.perc) #just a few outliers; median is 12%

#2. Chl (total); assuming Chl is 6.27% of chl on a mass basis
  #Chla = C55 H72 O5 N4 Mg = molecular weight of 893.5
dat$Chl.g = as.numeric(as.character(dat$Total.Chlorophyll)) * as.numeric(as.character(dat$Specific.leaf.area))  #ug Chl per g leaf
dat$ChlN = dat$Chl.g / (1/.0627) #ug N in Chl per g leaf
dat$Chl.perc = (dat$ChlN / 1000000) / (as.numeric(as.character(dat$leafN))/100) #Rubisco g per g leaf N
hist(dat$Chl.perc)
summary(dat$Chl.perc)
  #fairly low percentage, median is 1.4%, nearly all <5%
  #this is nearly exactly what is expected (Evans 1989)
boxplot(dat$Chl.perc~dat$Region)

#3. Water soluble N
ws.ena.perc = (as.numeric(as.character(dat$Water.soluble.Protein[dat$Region=="ENA"]))/(16*1000)) / (as.numeric(as.character(dat$leafN[dat$Region=="ENA"]))/100)
hist(ws.ena.perc)
  #seemingly too low
ws.japan.perc = (as.numeric(as.character(dat$Water.soluble.Protein[dat$Region=="Japan"]))/(16*1000)) / (as.numeric(as.character(dat$leafN[dat$Region=="Japan"]))/100)
hist(ws.japan.perc)
  #seems much more reasonable
summary(ws.japan.perc)
  #median is 20%

#4. SDS soluble N
sds.japan.perc = (as.numeric(as.character(dat$SDS.soluble.Protein[dat$Region=="Japan"]))/(16*1000)) / (as.numeric(as.character(dat$leafN[dat$Region=="Japan"]))/100)
hist(sds.japan.perc)
summary(sds.japan.perc)
hist(sds.japan.perc + ws.japan.perc)
  #a few outliers but reasonable

#5. Cell wall N
#% of leaf N allocated to cell wall = [cell wall mass per leaf area] x [%N of cell wall] / [leaf Narea]
dat$wallN.area = as.numeric(as.character(dat$Cell.wall..N)) * as.numeric(as.character(dat$Cell.wall.mass))/1000
  #however, don't have area of leaf disc removed... need it to convert to cell wall N mass per g leaf N


