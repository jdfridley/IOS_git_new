#create master licor dataset
#NSF IOS project
#Fridley Jan 2022

#library(data.table)
library(readxl)
library(stringr)
library(googlesheets4)

#sample master downloaded from google
master = read_sheet("https://docs.google.com/spreadsheets/d/1QLL5AUeP-fHar0HRfDjDP0Mw25FvowqidiCyEbBA5og/edit#gid=0")
#master = read.csv("/Users/fridley/Documents/academic/projects/IOS_FranceJapan/NSF-IOS spreadsheet_Jan24.csv")
str(master)
unique(master$Species)
  #need to fix latin name spellings and capitalization
unique(master$Sampling.Date)
  #need to update with months and days for some
unique(master$Site.final)
master$Date = sapply(master$Sampling.Date,function(x)as.character(x))

#folder of licor files
folder.local = "/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/LICOR_data files_Japan/"
files = list.files(folder.local)

#output
out = NULL
i = 1

#loop over each file, examine, and add to 'out'
for(i in 1:length(files)) {

#sample information to include in licor dataset
date = master$Date[master$LICOR.file.in.hand==files[i]]
site = master$Site.final[master$LICOR.file.in.hand==files[i]]
species = master$Species[master$LICOR.file.in.hand==files[i]]
code = master$Species.abbreviation[master$LICOR.file.in.hand==files[i]]
ID = master$ID[master$LICOR.file.in.hand==files[i]]

if(str_sub(files[i],start=-3)  == "csv") { #csv files
dat = read.table(paste0(folder.local,files[i]),header=F,blank.lines.skip=F,sep="\t")
start = which(substr(dat[,1],1,3)=="Obs") - 1
dat = read.csv(paste0(folder.local,files[i]),header=T,skip=start)
dat = dat[(!is.na(dat$Time)),] #remove rows of no data
dat$filename = files[i] #append column of filename
dat$date = date
dat$site = site
dat$species = species
dat$sppcode = code
dat$ID = ID
par(mfrow=c(1,2))
plot(dat$PARi,dat$Photo,pch=19,cex=1.3,main="A-q")
plot(dat$Ci,dat$Photo,pch=19,cex=1.3,main="A-Ci")
mtext(files[i],at=-300,line=0)
out = rbind(out,dat)
#readline() #uncomment if you want to inspect each curve upon import
}
  
if(str_sub(files[i],start=-4) == "xlsx" | str_sub(files[i],start=-4) == "xls" ) { #excel files
dat = read_excel(paste0(folder.local,files[i]),col_names=FALSE)
start = which(dat[,1]=="Obs") - 1
dat = as.data.frame(read_excel(paste0(folder.local,files[i]),col_names=TRUE,skip=start))
dat = dat[(!is.na(dat$Time)),] #remove rows of no data
dat$filename = files[i] #append column of filename
dat$date = date
dat$site = site
dat$species = species
dat$sppcode = code
dat$ID = ID
par(mfrow=c(1,2))
plot(dat$PARi,dat$Photo,pch=19,cex=1.3,main="A-q")
plot(dat$Ci,dat$Photo,pch=19,cex=1.3,main="A-Ci")
mtext(files[i],at=-300,line=0)
out = rbind(out,dat)
#readline() #uncomment if you want to inspect each curve upon import
}
}

#write.csv(out,file="/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/Japan_licor_final.csv")
