##Setting directories and loading packages

source(config.R)
library(data.table,lib.loc = packages.dir)
library(MCMCpack,lib.loc = packages.dir)
library(plyr,lib.loc = packages.dir)
library(MASS,lib.loc = packages.dir)
library(ggplot2,lib.loc = packages.dir)
library(synthpop,lib.loc = packages.dir)

##Loading the two files that were downloaded directly from the ACS
file.1 = paste(data.dir,"ss16pusa.csv",sep="/")
file.2 = paste(data.dir,"ss16pusb.csv",sep = "/")
df.1=fread(file.1)
df.2=fread(file.2)

##Merging these files together
df=rbind(df.1,df.2)
dim(df)
str(df)

df=as.data.frame(df)

##Variables of interest are citizenship, race, sex, age, income, and geography
##Extracting these characteristics
names = c("CIT","AGEP","RACWHT","SEX","ST","PINCP","PWGTP")

df.subset=df[,names]

dim(df.subset)

head(df.subset)


##Now, recoding one by one

##Citizenship
##1 if U.S. citizen, 0 if not
df.subset$CIT = ifelse(df.subset$CIT==5,0,1)


##Age
##1 if age is equal to or over 18, 0 otherwise
df.subset$AGEP = ifelse(df.subset$AGEP >= 18, 1, 0)


##RACWHT is already coded to be 0/1
##0 is nonwhite/1 is white

##SEX is already coded to be 1/2 -> 1 for males, 2 for females
##Just subtracting 1 so its consistent with 0 for males, 1 for females
df.subset$SEX=df.subset$SEX-1

##ST
##Dropping Puerto Rico first
df.subset=df.subset[df.subset$ST!=72,]
##Now, we are recoding based on Northeast/Midwest/South/West
##1=northeast/2=midwest/3=south/4=west/5 = pacific
ne = c(9,23,25,33,34,36,42,44,50)
mw = c(17,18,19,20,26,27,29,31,38,39,46,55)
south = c(1,5,10,11,12,13,21,22,24,28,37,40,45,47,48,51,54)
west = c(4,6,8,16,30,32,35,41,49,53,56)
pacific = c(2,15)

df.subset$ST = ifelse(df.subset$ST %in% ne, 1, ifelse(df.subset$ST %in% mw, 2, ifelse(df.subset$ST %in% south, 3, ifelse(df.subset$ST %in% west, 4, 5))))

##PINCP
##We are going to do above/below the poverty line
##This depends on the total number of people in the household
##Right now, we are just going to assume $11,880 which was the FPL for 1 person
df.subset$PINCP = ifelse(df.subset$PINCP > 11880, 1, 0)

##Now, deleting NAs
df.subset = na.omit(df.subset)

##And saving the data frame
fwrite(df.subset, paste(data.dir,"ACSsubset.csv",sep = "/"))
