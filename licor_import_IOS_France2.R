#create master licor dataset
#NSF IOS project
#Fridley Feb 2022

#France files: 3 types of licor files 
#note I could not get .xls files from France to import; used raw output text files instead

#library(data.table)
library(readxl)
library(stringr)
library(lubridate)
library(nlme)

#sample master downloaded from google
master = read.csv("/Users/fridley/Documents/academic/projects/IOS_FranceJapan/NSF-IOS spreadsheet_Feb2.csv")
str(master)
unique(master$Species)
  #need to fix latin name spellings and capitalization
sort(table(master$Species))
unique(master$Sampling.Date)
  #need to update with months and days for somegit
unique(master$Site.final)

#folder of licor files
folder.local = "/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/LICOR_data_files_France_rev/"
files = list.files(folder.local)

#create desired list of filenames without .xls
fnames = master$LICOR.file.in.hand[master$Region=="France"]
for(i in 1:length(fnames)) {if(str_sub(fnames[i],start=-4)==".xls" fnames[i] = }

#compare spreadsheet labels with filenames

#output
out = NULL #2019 text files
out2 = NULL #2020-21 text files
out3 = NULL #2019 csv files: not working, may need to re-save by hand due to '.' and ',' decimal delimiters

#loop over each file, examine, and add to 'out'
for(i in 1:length(files)) {

#sample information to include in licor dataset: needs work
#date = master$Sampling.Date[master$LICOR.file.in.hand==files[i]]
#site = master$Site.final[master$LICOR.file.in.hand==files[i]]
#species = master$Species[master$LICOR.file.in.hand==files[i]]
#code = master$Species.abbreviation[master$LICOR.file.in.hand==files[i]]

#2019 xls files (unique format)
if(substr(files[i],1,4) == "2019") {
dat = read.table(paste0(folder.local,files[i]),header=FALSE,sep="\t")
start = which(dat$V1=="$STARTOFDATA$") 
dat = read.table(paste0(folder.local,files[i]),header=TRUE,skip=start,sep="\t",fill=T)
dat = dat[(!is.na(dat$FTime)),] #remove rows of no data
dat$filename = files[i] #append column of filename
#dat$date = date
#dat$site = site
#dat$species = species
#dat$sppcode = code
par(mfrow=c(1,2))
plot(dat$PARi,dat$Photo,pch=19,cex=1.3,main="A-q")
plot(dat$Ci,dat$Photo,pch=19,cex=1.3,main="A-Ci")
mtext(files[i],at=-300,line=0)
out = rbind(out,dat)
#readline() #uncomment if you want to inspect each curve upon import
}

#2019 csv files: INGORE
  if(str_sub(files[i],start=-3)  == "csv") { next #csv files 
    dat = read.delim2(paste0(folder.local,files[i]),header=T)
    start = which(substr(dat[,1],1,3)=="Obs") - 1
    dat = read.csv(paste0(folder.local,files[i]),header=T,skip=start)
    dat = dat[(!is.na(dat$Time)),] #remove rows of no data
    dat$filename = files[i] #append column of filename
    #dat$date = date
    #dat$site = site
    #dat$species = species
    #dat$sppcode = code
    par(mfrow=c(1,2))
    plot(dat$PARi,dat$Photo,pch=19,cex=1.3,main="A-q")
    plot(dat$Ci,dat$Photo,pch=19,cex=1.3,main="A-Ci")
    mtext(files[i],at=-300,line=0)
    out = rbind(out,dat)
    #readline() #uncomment if you want to inspect each curve upon import
  }  
  
#2020/2021 xls files
if(substr(files[i],1,4) == "2020" | substr(files[i],1,4) == "2021") {
    dat = read.table(paste0(folder.local,files[i]),header=FALSE,sep="\t")
    start = which(dat$V1=="$STARTOFDATA$") 
    dat = read.table(paste0(folder.local,files[i]),header=TRUE,skip=start,sep="\t",fill=T)
    dat = dat[(!is.na(dat$FTime)),] #remove rows of no data
    dat$filename = files[i] #append column of filename
    #dat$date = date
    #dat$site = site
    #dat$species = species
    #dat$sppcode = code
    par(mfrow=c(1,2))
    plot(dat$PARi,dat$Photo,pch=19,cex=1.3,main="A-q")
    plot(dat$Ci,dat$Photo,pch=19,cex=1.3,main="A-Ci")
    mtext(files[i],at=-300,line=0)
    out2 = rbind(out2,dat)
    #readline() #uncomment if you want to inspect each curve upon import
  }  

}

#write.csv(out,file="/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/France_licor1.csv")
#write.csv(out2,file="/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/France_licor2.csv")
