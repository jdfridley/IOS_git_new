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
unique(master$Sampling.Date)
unique(master$Site.final)
master$Date = sapply(master$Sampling.Date,function(x)as.character(x))


#folder of licor files
folder.local = "/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/LICOR data files_ENA/" #downloaded after Robert's fixes, 2-8-22
files = list.files(folder.local)
length(files)

#compare files to spreadsheet filenames
ena.labels = master$LICOR.file.in.handENA[master$Region=="ENA"]
  #length = 150
#xls files saved as xlsx files; add 'x' to end of spreadsheet label names
#for(i in 1:length(ena.labels)) {
#  if(str_sub(ena.labels[i],start=-3)  == "xls") ena.labels[i] = paste0(ena.labels[i],"x")
#}
setdiff(ena.labels,files)
#two problem files:
  #"210809_s31_sogi.xlsx"  is empty (appears to be a duplicate, ignore)
  #"190806_s12_romul.xlsx" is corrupted, ignore for now
  #earlier: these names were fixed on 2-8-22: "210808_s27_bifro.xls"          "combined_190625_S6_Bethu.xlsx" "combined_190625_S6_ceorb.xlsx"

#remove 2 6800 files: incompatible merge with 6400 file format (add manually)
files = files[-which(files=="2021-09-13-bifro-dewitt.xlsx" | files == "2021-09-17_amart dewitt.xlsx" )]

#file output types:
#1. 2019 csv files: standard csv with no headers
#2. 2021 text files (2): 6800 output, tab delim, investigate
#3. 2000 text files (8): 6400 output, tab delim, std header info
#the rest are xlsx files (JDF updated)
  #xlsx files can be imported successfully
  #xls files cannot unless saved again as xlsx files: see
    #https://github.com/tidyverse/readxl/issues/598


#output
out = NULL
i=1

#loop over each file, examine, and add to 'out'
for(i in i:length(files)) {

#sample information to include in licor dataset
date = master$Date[master$LICOR.file.in.handENA==files[i]]
site = master$Site.final[master$LICOR.file.in.handENA==files[i]]
species = master$Species[master$LICOR.file.in.handENA==files[i]]
code = master$Species.abbreviation[master$LICOR.file.in.handENA==files[i]]
ID = master$ID[master$LICOR.file.in.handENA==files[i]]

if(str_sub(files[i],start=-3)  == "csv") { #csv files
dat = read.table(paste0(folder.local,files[i]),header=F,blank.lines.skip=F,sep="\t")
start = which(substr(dat[,1],1,3)=="Obs") - 1
dat = read.csv(paste0(folder.local,files[i]),header=T,skip=start)
dat = dat[(!is.na(dat$FTime)),] #remove rows of no data
dat$filename = files[i] #append column of filename
dat$date = date
dat$site = site
dat$species = species
dat$sppcode = code
dat$ID = ID
names(dat) = names(out)
par(mfrow=c(1,2))
plot(dat$PARi,dat$Photo,pch=19,cex=1.3,main="A-q")
plot(dat$Ci,dat$Photo,pch=19,cex=1.3,main="A-Ci")
mtext(files[i],at=-300,line=0)
out = rbind(out,dat)
#readline() #uncomment if you want to inspect each curve upon import
}
  
if(gregexpr(".",files[i],fixed=T)<0) { #tab delim text files
  dat = read.table(paste0(folder.local,files[i]),header=F,blank.lines.skip=F,sep="\t")
  start = which(substr(dat[,1],1,3)=="Obs") - 1
  dat = read.csv(paste0(folder.local,files[i]),header=T,skip=start,sep="\t")
  dat = dat[(!is.na(dat$FTime)),] #remove rows of no data
  dat[,41:58] = NA  #add 18 extra rows (missing in this version of the licor files)
  dat$filename = files[i] #append column of filename
  dat$date = date
  dat$site = site
  dat$species = species
  dat$sppcode = code
  dat$ID = ID
    names(dat) = names(out)
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
dat = dat[-1,]
dat = dat[(!is.na(dat$FTime)),] #remove rows of no data
dat[,c(1,3:58)] = apply(dat[,c(1,3:58)], 2, function(x) as.numeric(as.character(x))) #change character to numeric (all except HHMMSS)
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

#write.csv(out,file="/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/ENA_licor2.csv")
