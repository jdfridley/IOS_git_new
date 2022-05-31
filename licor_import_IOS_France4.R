#create master licor dataset
#NSF IOS project
#Fridley Jan 2022

#library(data.table)
library(readxl)
library(stringr)
library(googlesheets4)
library(dplyr)

#sample master downloaded from google
master = read_sheet("https://docs.google.com/spreadsheets/d/1QLL5AUeP-fHar0HRfDjDP0Mw25FvowqidiCyEbBA5og/edit#gid=0")
#master = read.csv("/Users/fridley/Documents/academic/projects/IOS_FranceJapan/NSF-IOS spreadsheet_Jan24.csv")
str(master)
unique(master$Species)
unique((master$Sampling.Date))
unique(master$Site.final)
master$Date = sapply(master$Sampling.Date,function(x)as.character(x))


#folder of licor files ***use the current France folder on as-fridley (correct files and filenames)
folder.local = "/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/LICOR_data_files_France_Marchupdate/" 
#folder.local = "C:\\jdfhome\\NSFjapanfrance\\LICOR\\LICOR data files_France\\" #on PC
files = list.files(folder.local)
length(files)

#compare files to spreadsheet filenames
labels = master$LICOR.file.in.handFrance[master$Region=="France"]
length(labels)  
setdiff(labels,files)
  #all good


#output
out = NULL
i = 1

#loop over each file, examine, and add to 'out'
for(i in i:length(files)) {

#sample information to include in licor dataset
date = master$Date[master$LICOR.file.in.handFrance==files[i]]
site = master$Site.final[master$LICOR.file.in.handFrance==files[i]]
species = master$Species[master$LICOR.file.in.handFrance==files[i]]
code = master$Species.abbreviation[master$LICOR.file.in.handFrance==files[i]]
ID = master$ID[master$LICOR.file.in.handFrance==files[i]]

#all xlsx files
  if(str_sub(files[i],start=-4) == "xlsx") {
  dat = read_excel(paste0(folder.local,files[i]),col_names=FALSE)
  start = which(dat[,1]=="Obs") - 1
  dat = as.data.frame(read_excel(paste0(folder.local,files[i]),col_names=TRUE,skip=start))
  dat = dat[-1,]
  dat = dat[(!is.na(dat$FTime)),] #remove rows of no data
  dat[,c(1,3:56)] = apply(dat[,c(1,3:56)], 2, function(x) as.numeric(as.character(x))) #change character to numeric (all except HHMMSS)
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
  if(class(dat[,57]) == "character" & dim(dat)[2]>65) {
    dat[,c(57:83)] = apply(dat[,c(57:83)], 2, function(x) as.numeric(as.character(x))) #change character to numeric (all except HHMMSS)
  }
  out = bind_rows(out,dat) #from dplyr
  #readline() #uncomment if you want to inspect each curve upon import
}

#text files: all have 36 columns
  if(str_sub(files[i],start=-4) != "xlsx") {
  dat = read.table(paste0(folder.local,files[i]),header=F,blank.lines.skip=F,sep="\t")
  start = which(substr(dat[,1],1,3)=="Obs") - 1
  dat = read.csv(paste0(folder.local,files[i]),header=T,skip=start,sep="\t")
  dat = dat[(!is.na(dat$FTime)),] #remove rows of no data
  dat$filename = files[i] #append column of filename
  dat$date = date
  dat$site = site
  dat$species = species
  dat$sppcode = code
  dat$ID = ID
  dat$Obs = as.numeric(dat$Obs)
  par(mfrow=c(1,2))
  plot(dat$PARi,dat$Photo,pch=19,cex=1.3,main="A-q")
  plot(dat$Ci,dat$Photo,pch=19,cex=1.3,main="A-Ci")
  mtext(files[i],at=-300,line=0)
  out = bind_rows(out,dat) #from dplyr
  #readline() #uncomment if you want to inspect each curve upon import
} 
}

length(table(out$filename))
str(out)
  #167 files
hist(out$Photo)
hist(out$Cond)
hist(out$PARi)
hist(out$Ci)
plot(out$Ci,out$Photo,col=as.numeric(as.factor(out$filename)))

#write.csv(out,file="/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/France_licor_final.csv")
#write.csv(out,"France_licor3.csv")

##combine licor datafiles from all three regions

france = read.csv("/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/France_licor_final.csv")
ena = read.csv("/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/ENA_licor2.csv")
japan = read.csv("/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/Japan_licor_final2.csv")

all = bind_rows(france,ena,japan) #from dplyr

#write.csv(all,file="/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/all_licor_data.csv")

