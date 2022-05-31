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
unique(master$Sampling.Date)
unique(master$Site.final)


#folder of licor files ***use the current France folder on as-fridley (correct files and filenames)
folder.local = "/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/LICOR_data_files_France_Marchupdate/" 
#folder.local = "C:\\jdfhome\\NSFjapanfrance\\LICOR\\LICOR data files_France\\" #on PC
files = list.files(folder.local)
length(files)

#compare files to spreadsheet filenames
labels = master$LICOR.file.in.handFrance[master$Region=="France"]
  #length = 150
#xls files saved as xlsx files; add 'x' to end of spreadsheet label names
#for(i in 1:length(ena.labels)) {
#  if(str_sub(ena.labels[i],start=-3)  == "xls") ena.labels[i] = paste0(ena.labels[i],"x")
#}
setdiff(labels,files)
  #all good

#file output types:
#1. 9 .xls files (do separately)
#2. tab-delim text files

#complicated set of file formats:
#2019: 15 text files have 58 columns, 7 xlsx files have 83 columns
#2020: all text files have 36 columns, 2 xlsx files have 56 columns
#2021: all text files have 36 columns

i = 1 #needed for an unknown reason

#output
out = NULL
out19 = NULL
out19xls = NULL
out20xls = NULL

#loop over each file, examine, and add to 'out'
for(i in i:length(files)) {

#sample information to include in licor dataset
date = master$Sampling.Date[master$LICOR.file.in.handENA==files[i]]
site = master$Site.final[master$LICOR.file.in.handENA==files[i]]
species = master$Species[master$LICOR.file.in.handENA==files[i]]
code = master$Species.abbreviation[master$LICOR.file.in.handENA==files[i]]

print(i)

#1. 2020 and 2021 text files: all have 36 columns
if(gregexpr(".",files[i],fixed=T)<0 & substr(files[i],1,4)!="2019") { #tab delim text files from 2020, 2021
  dat = read.table(paste0(folder.local,files[i]),header=F,blank.lines.skip=F,sep="\t")
  start = which(substr(dat[,1],1,3)=="Obs") - 1
  dat = read.csv(paste0(folder.local,files[i]),header=T,skip=start,sep="\t")
  dat = dat[(!is.na(dat$FTime)),] #remove rows of no data
  dat$filename = files[i] #append column of filename
  dat$date = date
  dat$site = site
  dat$species = species
  dat$sppcode = code
  par(mfrow=c(1,2))
  plot(dat$PARi,dat$Photo,pch=19,cex=1.3,main="A-q")
  plot(dat$Ci,dat$Photo,pch=19,cex=1.3,main="A-Ci")
  mtext(files[i],at=-300,line=0)
  out = rbind(out,dat)
  readline() #uncomment if you want to inspect each curve upon import
}

print(i)

#2. 2019 text files: 58 columns
if(gregexpr(".",files[i],fixed=T)<0 & substr(files[i],1,4)=="2019") { #tab delim text files from 2019
  dat = read.table(paste0(folder.local,files[i]),header=F,blank.lines.skip=F,sep="\t")
  start = which(substr(dat[,1],1,3)=="Obs") - 1
  dat = read.csv(paste0(folder.local,files[i]),header=T,skip=start,sep="\t")
  dat = dat[(!is.na(dat$FTime)),] #remove rows of no data
  dat$filename = files[i] #append column of filename
  dat$date = date
  dat$site = site
  dat$species = species
  dat$sppcode = code
  par(mfrow=c(1,2))
  plot(dat$PARi,dat$Photo,pch=19,cex=1.3,main="A-q")
  plot(dat$Ci,dat$Photo,pch=19,cex=1.3,main="A-Ci")
  mtext(files[i],at=-300,line=0)
  out19 = rbind(out19,dat)
  readline() #uncomment if you want to inspect each curve upon import
}

print(i)

#3. 2019 excel files
if(str_sub(files[i],start=-4) == "xlsx" & substr(files[i],1,4)=="2019") { #excel files of 2019
  dat = read_excel(paste0(folder.local,files[i]),col_names=FALSE)
  start = which(dat[,1]=="Obs") - 1
  dat = as.data.frame(read_excel(paste0(folder.local,files[i]),col_names=TRUE,skip=start))
  dat = dat[-1,]
  dat = dat[(!is.na(dat$FTime)),] #remove rows of no data
  dat[,c(1,3:83)] = apply(dat[,c(1,3:83)], 2, function(x) as.numeric(as.character(x))) #change character to numeric (all except HHMMSS)
  dat$filename = files[i] #append column of filename
  dat$date = date
  dat$site = site
  dat$species = species
  dat$sppcode = code
  par(mfrow=c(1,2))
  plot(dat$PARi,dat$Photo,pch=19,cex=1.3,main="A-q")
  plot(dat$Ci,dat$Photo,pch=19,cex=1.3,main="A-Ci")
  mtext(files[i],at=-300,line=0)
  out19xls = rbind(out19xls,dat)
  readline() #uncomment if you want to inspect each curve upon import
}

print(i)

#4. 2020 excel files
if(str_sub(files[i],start=-4) == "xlsx" & substr(files[i],1,4)=="2020") { #excel files of 2019
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
  par(mfrow=c(1,2))
  plot(dat$PARi,dat$Photo,pch=19,cex=1.3,main="A-q")
  plot(dat$Ci,dat$Photo,pch=19,cex=1.3,main="A-Ci")
  mtext(files[i],at=-300,line=0)
  out20xls = rbind(out20xls,dat)
  readline() #uncomment if you want to inspect each curve upon import
}

print(i)

} #close files loop

dim(out)
dim(out19)
dim(out19xls)
dim(out20xls)

#combine 4 out files, getting the columns in the right places
out1 = bind_rows(out,out19) #from dplyr
out1[,c(1,3:36,42:63)] = apply(out1[,c(1,3:36,42:63)], 2, function(x) as.numeric(as.character(x))) #change character to numeric (all except HHMMSS)
out2 = bind_rows(out19xls,out20xls)
out.final = bind_rows(out1,out2)
str(out.final)
  #164 files
out.final2 = out.final[!is.na(out.final$Photo),] #still some NA rows; omit


#write.csv(out,file="/Users/fridley/Documents/academic/projects/IOS_FranceJapan/licor_files/France_licor2.csv")
#write.csv(out.final2,"France_licor2.csv")
