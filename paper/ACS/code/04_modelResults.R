source(config.R)

library(data.table,lib.loc = packages.dir)
library(MCMCpack,lib.loc = packages.dir)
library(plyr,lib.loc = packages.dir)
library(MASS,lib.loc = packages.dir)
library(ggplot2,lib.loc = packages.dir)
library(synthpop,lib.loc = packages.dir)
library(tidyr,lib.loc = packages.dir)

##Reading in the model data
ACS.df = fread(paste(data.dir,"ACSsubset.csv",sep="/"))

##We use a sample size of 10,000, but can be changed here
samp.size = 10000

##Seed for reproducability
set.seed(1)

##Subsetting the data so it has a sample size = samp.size
samps = sample(1:nrow(ACS.df),samp.size,replace=FALSE)
ACS.df =  ACS.df[samps,]

##Subsetting the variables that we want
ACS.df.tab = as.data.frame(table(ACS.df[,c(1,2,3,4,6)]))
##These are the true data counts
trueCounts = plyr::count(ACS.df.tab)$Freq
##Since different functions aggregate the data differently, this is reordering the counts data to match the output of later functions
trueData = plyr::count(ACS.df.tab)[c(1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31,2,18,10,26,6,22,14,30,4,20,12,28,8,24,16,32),-7]

##Fitting the full (saturated) log linear model on the orginal data
log.mod = loglm(Freq~CIT+AGEP+RACWHT+SEX+PINCP+CIT*AGEP+CIT*RACWHT+CIT*SEX+CIT*PINCP+AGEP*RACWHT+AGEP*SEX+AGEP*PINCP+RACWHT*SEX+RACWHT*PINCP+SEX*PINCP,data=trueData)
##Using a stepwise approach to select important terms
step.log.mod = step(log.mod)
##Renaming the terms
names = c("CIT","AGEP","RACWHT","SEX","PINCP","CIT:AGEP","CIT:RACWHT","CIT:SEX","CIT:PINCP","AGEP:RACWHT","AGEP:SEX","AGEP:PINCP","RACWHT:SEX","RACWHT:PINCP","SEX:PINCP")

##Creating a list of terms selected by the procedure
selected.real = lapply(step.log.mod,attributes)$term$term.labels
selected.real = as.integer(names %in% selected.real)

##A vector of coefficients for the saturated log linear model
coefs = glm(Freq~.,data=ACS.df.tab,family = "poisson")$coefficients
S=c(plyr::count(ACS.df,c("CIT","AGEP"))$freq,plyr::count(ACS.df,c("CIT","RACWHT"))$freq,plyr::count(ACS.df,c("CIT","SEX"))$freq,plyr::count(ACS.df,c("CIT","PINCP"))$freq,plyr::count(ACS.df,c("AGEP","RACWHT"))$freq,plyr::count(ACS.df,c("AGEP","SEX"))$freq,plyr::count(ACS.df,c("AGEP","PINCP"))$freq,plyr::count(ACS.df,c("RACWHT","SEX"))$freq,plyr::count(ACS.df,c("RACWHT","PINCP"))$freq,plyr::count(ACS.df,c("SEX","PINCP"))$freq)

####Graphs using the underlying non-noisy model
N.runs = 10000

set.seed(2020)
P=TwoWay_LCM(S,G=5,samp.size,N.runs=N.runs,prop.weights)
probs.reduced.inds = seq(5000,10000,by=5)
P.full = P$Full[probs.reduced.inds,]
P.marg = P$Marginal[probs.reduced.inds,]

dat.sum = data.frame(Predicted = colMeans(P.full), True = trueData$Freq/samp.size)

dodge <- position_dodge(width = 0)
ggplot(dat.sum, aes(x=True, y=Predicted)) + 
  geom_point() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=16),legend.position = "none") + 
  geom_abline(slope = 1, intercept = 0) + xlab("True Probability") + ylab("Estimated Probability")

ggsave(paste(results.dir,"fullProbs_noPrivacy.pdf",sep="\\"))

dat.sum = data.frame(Predicted = colMeans(P.marg), True = S/samp.size)

dodge <- position_dodge(width = 0)
ggplot(dat.sum, aes(x=True, y=Predicted)) + 
  geom_point() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=16),legend.position = "none") + 
  geom_abline(slope = 1, intercept = 0) + xlab("True Probability") + ylab("Estimated Probability")

ggsave(paste(results.dir,"margProbs_noPrivacy.pdf",sep="\\"))


##Now creating graphs of full and marginal table probabilities
epsilon = c(.25,.5,1)
ep.names = c("0_25","0_50","1_00")
E = -(epsilon/10)/2

##The number of  repeat simulations
kk = 5
Ps.full = array(NA, c(kk, 32,length(epsilon)))
Ps.marg = array(NA, c(kk, 40,length(epsilon)))
Ps.noisy = array(NA, c(kk, 40,length(epsilon)))

for(q in 1:length(epsilon)){
  for(i in 1:kk){
    set.seed(i)
    errs = rgeom(40,1-exp(E[q])) - rgeom(40,1-exp(E[q]))
    
    noisyS = S + errs
    noisyS = ifelse(noisyS<0,0,noisyS)
    prop.weights = noisyS
    P=TwoWay_NoisyLCM(noisyS,G=5,samp.size,N.runs=10000,E[q],prop.weights)
    probs.reduced.inds = seq(5000,10000,by=5)
    P.full = P$Full[probs.reduced.inds,]
    P.marg = P$Marginal[probs.reduced.inds,]
    
    Ps.full[i,,q] = colMeans(P.full)
    Ps.marg[i,,q] = colMeans(P.marg)
    
    write.csv(P.full,file.name = paste(results.dir, "/LCMFullE",ep.names[q],"_",i,".csv",sep=""))
    write.csv(P.marg,file.name = paste(results.dir, "/LCMMargE",ep.names[q],"_",i,".csv",sep=""))

    print(i)
  }
  print(q)
}

###Creating the graphs
dat.sum = data.frame(Probs = c(colMeans(Ps.full[,,1]),colMeans(Ps.full[,,2]),colMeans(Ps.full[,,3])),sd = c(apply(Ps.full[,,1], 2, sd),apply(Ps.full[,,2], 2, sd),apply(Ps.full[,,3], 2, sd)),Epsilon = c(rep(0.25,32),rep(0.50,32),rep(1.00,32)),True.Probs = rep(trueData$Freq/samp.size,3))
dat.sum$Epsilon = factor(dat.sum$Epsilon)
ggplot(dat.sum, aes(x=True.Probs, y=Probs, colour=Epsilon)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=16)) + xlab("True Probability") + ylab("Estimated Probability")
ggsave(paste(results.dir,"fullCombine.pdf",sep="\\"))

dat.sum = data.frame(Probs = c(colMeans(Ps.marg[,,1]),colMeans(Ps.marg[,,2]),colMeans(Ps.marg[,,3])),sd = c(apply(Ps.marg[,,1], 2, sd),apply(Ps.marg[,,2], 2, sd),apply(Ps.marg[,,3], 2, sd)),Epsilon = c(rep(0.25,40),rep(0.50,40),rep(1.00,40)),True.Probs = rep(S/samp.size,3))
dat.sum$Epsilon = factor(dat.sum$Epsilon)
ggplot(dat.sum, aes(x=True.Probs, y=Probs, colour=Epsilon)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=16)) + xlab("True Probability") + ylab("Estimated Probability")
ggsave(paste(results.dir,"margCombine.pdf",sep="\\"))



##Now coverage probabilities for the marginal counts
epsilon = c(.25,.5,1)
E = -(epsilon/10)/2

ep_name = c("E0_25","E0_50", "E1_00")

##The number of  repeat simulations
kk = 100

P.mat.marg = matrix(NA,nrow=kk,ncol=40)
CI.length.marg =  matrix(NA,nrow=kk,ncol=40)

for(j in 1:length(E)){
  for(i in 1:kk){
    set.seed(i*i)
    
    ##Subsetting the data so it has a sample size = samp.size
    samps = sample(1:nrow(ACS.df),samp.size,replace=FALSE)
    ACS.df =  ACS.df[samps,]
    
    S=c(plyr::count(ACS.df,c("CIT","AGEP"))$freq,plyr::count(ACS.df,c("CIT","RACWHT"))$freq,plyr::count(ACS.df,c("CIT","SEX"))$freq,plyr::count(ACS.df,c("CIT","PINCP"))$freq,plyr::count(ACS.df,c("AGEP","RACWHT"))$freq,plyr::count(ACS.df,c("AGEP","SEX"))$freq,plyr::count(ACS.df,c("AGEP","PINCP"))$freq,plyr::count(ACS.df,c("RACWHT","SEX"))$freq,plyr::count(ACS.df,c("RACWHT","PINCP"))$freq,plyr::count(ACS.df,c("SEX","PINCP"))$freq)
    
    errs = rgeom(40,1-exp(E[j])) - rgeom(40,1-exp(E[j]))
    
    noisyS = S + errs
    noisyS = ifelse(noisyS<0,0,noisyS)
    prop.weights = noisyS
    
    P=TwoWay_NoisyLCM(noisyS,G=5,samp.size,N.runs=15000,E[j],prop.weights)
    probs.reduced.inds = seq(5000,15000,by=10)
    
    ###calculating coverage probabilities
    tmp.P = P$Counts[,probs.reduced.inds]
    ##Converting to coda object
    chain = as.mcmc(t(tmp.P))
    
    MCMC.sums = summary(chain)
    se.mcmc = MCMC.sums$statistics[,2]
    ci.lower = MCMC.sums$statistics[,1]- qnorm(.975)*se.mcmc
    ci.upper = MCMC.sums$statistics[,1] + qnorm(.975)*se.mcmc
    
    P.mat.marg[i,] = ifelse((ci.lower <= S) & (ci.upper >= S),1,0)
    CI.length.marg[i,] = ci.upper-ci.lower
    
    file.name = paste(results.dir,paste("BayesLCM_Full_",ep_name,"_",i,".csv",sep=""),sep="\\")
    write.csv(P$Full,file = file.name)
    
    file.name = paste(results.dir,paste("BayesLCM_Marginal_",ep_name,"_",i,".csv",sep=""),sep="\\")
    write.csv(P$Marginal,file = file.name)
    
  }
}

file.name = paste(results.dir,paste("BayesLCM_covProb_",ep_name,".csv",sep=""),sep="\\")
write.csv(P.mat.marg,file = file.name)

file.name = paste(results.dir,paste("BayesLCM_CIlength_",ep_name,".csv",sep=""),sep="\\")
write.csv(CI.length.marg,file = file.name)


##Comparison graphs
##Epslion = .25


dfpriv = read.csv(paste(python.dir,"\\PrivBayes\\DataSynthesizer\\out\\correlated_attribute_mode\\privBayes_E0_0.csv",sep=""))[,-1]

dfpriv = as.data.frame(table(dfpriv))$Freq/samp.size


dflcm = read.csv(paste(results.dir, "\\BayesLCM_Full_E0_25_1.csv",sep=""))[,-1]
dflcm = colMeans(dflcm)

dfgm = read.csv(paste(python.dir,"\\Graphical Models Example\\GM_0_25_1.csv",sep=""),header=FALSE)/samp.size
dfgm = unlist(dfgm)

data.map = data.frame(True = c(rep((trueData$Freq/samp.size),2),trueCounts/samp.size),Method = c(rep("Proposed Method",32),rep("PrivBayes",32),rep("GM",32)),Syn = c(dflcm,dfpriv,dfgm))

data.map$Method = factor(data.map$Method)
ggplot(data=data.map,aes(x=True, y=Syn, group=Method,color=Method)) +
  geom_point()+
  geom_abline(slope=1, intercept=0) + xlab("True Probability") + ylab("Estimated Probability") + theme(axis.text=element_text(size=16),
                                                                                                                 axis.title=element_text(size=16),legend.text=element_text(size=16),legend.title=element_text(size=16),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(paste(results.dir,"compareE0_25.pdf",sep="/"))

##Epsilon = .50

dfpriv = read.csv(paste(python.dir,"\\PrivBayes\\DataSynthesizer\\out\\correlated_attribute_mode\\privBayes_E1_0.csv",sep=""))[,-1]

dfpriv = as.data.frame(table(dfpriv))$Freq/samp.size


dflcm = read.csv(paste(results.dir, "\\BayesLCM_Full_E0_50_1.csv",sep=""))[,-1]
dflcm = colMeans(dflcm)

dfgm = read.csv(paste(python.dir,"\\Graphical Models Example\\GM_0_50_1.csv",sep=""),header=FALSE)/samp.size
dfgm = unlist(dfgm)

data.map = data.frame(True = c(rep((trueData$Freq/samp.size),2),trueCounts/samp.size),Method = c(rep("Proposed Method",32),rep("PrivBayes",32),rep("GM",32)),Syn = c(dflcm,dfpriv,dfgm))

data.map$Method = factor(data.map$Method)
ggplot(data=data.map,aes(x=True, y=Syn, group=Method,color=Method)) +
  geom_point()+
  geom_abline(slope=1, intercept=0) + xlab("True Probability") + ylab("Estimated Probability") + theme(axis.text=element_text(size=16),
                                                                                                       axis.title=element_text(size=16),legend.text=element_text(size=16),legend.title=element_text(size=16),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(paste(results.dir,"compareE0_50.pdf",sep="/"))

###Epsilon - 1.00

dfpriv = read.csv(paste(python.dir,"\\PrivBayes\\DataSynthesizer\\out\\correlated_attribute_mode\\privBayes_E1_0.csv",sep=""))[,-1]

dfpriv = as.data.frame(table(dfpriv))$Freq/samp.size


dflcm = read.csv(paste(results.dir, "\\BayesLCM_Full_E1_00_1.csv",sep=""))[,-1]
dflcm = colMeans(dflcm)

dfgm = read.csv(paste(python.dir,"\\Graphical Models Example\\GM_1_00_1.csv",sep=""),header=FALSE)/samp.size
dfgm = unlist(dfgm)

data.map = data.frame(True = c(rep((trueData$Freq/samp.size),2),trueCounts/samp.size),Method = c(rep("Proposed Method",32),rep("PrivBayes",32),rep("GM",32)),Syn = c(dflcm,dfpriv,dfgm))

data.map$Method = factor(data.map$Method)
ggplot(data=data.map,aes(x=True, y=Syn, group=Method,color=Method)) +
  geom_point()+
  geom_abline(slope=1, intercept=0) + xlab("True Probability") + ylab("Estimated Probability") + theme(axis.text=element_text(size=16),
                                                                                                       axis.title=element_text(size=16),legend.text=element_text(size=16),legend.title=element_text(size=16),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(paste(results.dir,"compareE1_00.pdf",sep="/"))


###Inference results
###Use two models, w/ and w/o measurement error component
epsilon = c(0.25,.5,1)

privBayes.model1 = array(NA, c(100, 11, length(epsilon)))
privBayes.model2 = array(NA, c(100, 11, length(epsilon)))
LCM.model1 = array(NA, c(100, 11, length(epsilon)))
LCM.model2 = array(NA, c(100, 11, length(epsilon)))
GM.model1 = array(NA, c(100, 11, length(epsilon)))
GM.model2 = array(NA, c(100, 11, length(epsilon)))

privBayes.model1.se = array(NA, c(100, 11, length(epsilon)))
privBayes.model2.se = array(NA, c(100, 11, length(epsilon)))
LCM.model1.se = array(NA, c(100, 11, length(epsilon)))
LCM.model2.se = array(NA, c(100, 11, length(epsilon)))
GM.model1.se = array(NA, c(100, 11, length(epsilon)))
GM.model2.se = array(NA, c(100, 11, length(epsilon)))


##True Betas
mod.true.1 = glm(CIT ~  PINCP +AGEP + RACWHT +SEX+ AGEP*SEX + AGEP*RACWHT + AGEP*PINCP + PINCP*RACWHT + PINCP*SEX + RACWHT*SEX , data=ACS.df,family = "binomial")
mod.true.2 = glm(PINCP ~  CIT +AGEP + RACWHT +SEX+ AGEP*SEX + AGEP*RACWHT + AGEP*CIT + CIT*RACWHT + CIT*SEX + RACWHT*SEX, data=ACS.df,family = "binomial")


for(i in 1:length(epsilon)){
  visit = seq(1,5,by=1)
  names(visit) = c("CIT","AGEP","RACWHT","SEX","PINCP")
  method = rep("ipf",5)
  names(method) = c("CIT","AGEP","RACWHT","SEX","PINCP")
  
  for(q in 1:100){
    a = i-1
    b=q-1
    
    ##PrivBayes
    file.name = paste(python.dir,"\\PrivBayes\\DataSynthesizer\\out\\correlated_attribute_mode\\privBayes_E",a,"_",b,".csv",sep="")
    dfpriv = read.csv(file.name)[,-1]
    names(dfpriv) = c("CIT","AGEP","RACWHT","SEX","PINCP")
    
    fitt = glm(CIT ~  PINCP +AGEP + RACWHT +SEX+ AGEP*SEX + AGEP*RACWHT + AGEP*PINCP + PINCP*RACWHT + PINCP*SEX + RACWHT*SEX , data=dfpriv,family = "binomial")
    privBayes.model1[q,,i] = fitt$coefficients
    privBayes.model1.se[q,,i] = summary(fitt)$coefficients[,2]
    fitt = glm(PINCP ~  CIT +AGEP + RACWHT +SEX+ AGEP*SEX + AGEP*RACWHT + AGEP*CIT + CIT*RACWHT + CIT*SEX + RACWHT*SEX , data=dfpriv,family = "binomial")
    privBayes.model2[q,,i] = fitt$coefficients
    privBayes.model2.se[q,,i] = summary(fitt)$coefficients[,2]
    
    ##GM
    ep.names = c("0_25","0_50","1_00")
    file.name = paste(python.dir,"\\Graphical Models Example\\GM_",ep.names[i],"_",q,".csv",sep="")
    dfgm = unlist(read.csv(file.name,header=FALSE))
    df.gm = plyr::count(ACS.df.tab)
    df.gm$Freq = dfgm
    
    dfgm = uncount(df.gm,Freq)
    
    fitt = glm(CIT ~  PINCP +AGEP + RACWHT +SEX+ AGEP*SEX + AGEP*RACWHT + AGEP*PINCP + PINCP*RACWHT + PINCP*SEX + RACWHT*SEX , data=dfgm,family = "binomial")
    GM.model1[q,,i] = fitt$coefficients
    hold = data.frame(nam = names(summary(fitt)$coefficients[,2]), vals = summary(fitt)$coefficients[,2])
    merg = data.frame(nam = names(fitt$coefficients))
    tmp.merge = merge(hold,merg,by="nam",all=TRUE)
    row.names(tmp.merge) = tmp.merge$nam
    tmp.merge = tmp.merge[names(fitt$coefficients),]
    
    GM.model1.se[q,,i] = tmp.merge$vals
    fitt = glm(PINCP ~  CIT +AGEP + RACWHT +SEX+ AGEP*SEX + AGEP*RACWHT + AGEP*CIT + CIT*RACWHT + CIT*SEX + RACWHT*SEX , data=dfgm,family = "binomial")
    GM.model2[q,,i] = fitt$coefficients
    
    hold = data.frame(nam = names(summary(fitt)$coefficients[,2]), vals = summary(fitt)$coefficients[,2])
    merg = data.frame(nam = names(fitt$coefficients))
    tmp.merge = merge(hold,merg,by="nam",all=TRUE)
    row.names(tmp.merge) = tmp.merge$nam
    tmp.merge = tmp.merge[names(fitt$coefficients),]
    
    GM.model2.se[q,,i] = tmp.merge$vals
    
    
    ##LCM
    dflcm = read.csv(paste(results.dir,"\\BayesLCM_E",i,"_",q,".csv",sep=""))[,-1]
    ##Sample some data sets
    inds = sample(1:nrow(dflcm),50,replace=FALSE)
    dflcm = dflcm[inds,]
    dflcm = dflcm*samp.size
    lcmSyns = list()
    for(r in 1:50){
      df.lcm = ACS.df.tab
      df.lcm$Freq = round(unlist(dflcm[r,]))
      
      lcmSyns[[r]] = uncount(df.lcm,Freq)
    }#end of for
    Xsyn = list(syn = lcmSyns, m = 50,method=method,proper=FALSE,visit.sequence=visit)
    class(Xsyn) <- "synds"
    
    fitt = glm.synds(CIT ~  PINCP +AGEP + RACWHT +SEX+ AGEP*SEX + AGEP*RACWHT + AGEP*PINCP + PINCP*RACWHT + PINCP*SEX + RACWHT*SEX , data=Xsyn,family = "binomial")
    LCM.model1[q,,i] = fitt$mcoefavg
    LCM.model1.se[q,,i] = fitt$mvaravg
    fitt = glm.synds(PINCP ~  CIT +AGEP + RACWHT +SEX+ AGEP*SEX + AGEP*RACWHT + AGEP*CIT + CIT*RACWHT + CIT*SEX + RACWHT*SEX , data=Xsyn,family = "binomial")
    LCM.model2[q,,i] = fitt$mcoefavg
    LCM.model2.se[q,,i] = fitt$mvaravg
    print(q)
  }##End of for
}##End of for

names.mod1 = c("(Intercept)","PINCP","AGEP","RACWHT","SEX","AGEP:SEX","AGEP:RACWHT" ,"PINCP:AGEP","PINCP:RACWHT","PINCP:SEX","RACWHT:SEX") 
names.mod2 = c("(Intercept)", "CIT","AGEP","RACWHT","SEX" ,"AGEP:SEX","AGEP:RACWHT","CIT:AGEP","CIT:RACWHT","CIT:SEX","RACWHT:SEX" )

##Creating each model data frames for each epsilon (0.25)
df1 = as.data.frame(cbind(names.mod1,paste(round(colMeans(LCM.model1[,,1]),2), " (", round(colMeans(LCM.model1.se[,,1]),2), ")",sep="" ),paste(round(colMeans(privBayes.model1[,,1]),2), " (", round(colMeans(privBayes.model1.se[,,1]),2), ")",sep="" ),paste(round(colMeans(GM.model1[,,1], na.rm=TRUE),2), " (", round(colMeans(GM.model1.se[,,1], na.rm=TRUE),2), ")",sep="" )))
names(df1) = c("names","V1","V2","V3")
df1$true.mod1 = paste(round(summary(mod.true.1)$coefficients[,1],2), " (",round(summary(mod.true.1)$coefficients[,2],2) ,")", sep="")

df2  = as.data.frame(cbind(names.mod2,paste(round(colMeans(LCM.model2[,,1]),2), " (", round(colMeans(LCM.model2.se[,,1]),2), ")",sep="" ),paste(round(colMeans(privBayes.model2[,,1]),2), " (", round(colMeans(privBayes.model2.se[,,1]),2), ")",sep="" ),paste(round(colMeans(GM.model2[,,1]),2), " (", round(colMeans(GM.model2.se[,,1]),2), ")",sep="" )))

names(df2) = c("names","V4","V5","V6")
df2$true.mod2 = paste(round(summary(mod.true.2)$coefficients[,1],2), " (",round(summary(mod.true.2)$coefficients[,2],2) ,")", sep="")


##Creating each model data frames for each epsilon (0.50)
df1 = as.data.frame(cbind(names.mod1,paste(round(colMeans(LCM.model1[,,2]),2), " (", round(colMeans(LCM.model1.se[,,2]),2), ")",sep="" ),paste(round(colMeans(privBayes.model1[,,2]),2), " (", round(colMeans(privBayes.model1.se[,,2]),2), ")",sep="" ),paste(round(colMeans(GM.model1[,,2], na.rm=TRUE),2), " (", round(colMeans(GM.model1.se[,,2], na.rm=TRUE),2), ")",sep="" )))
names(df1) = c("names","V1","V2","V3")
df1$true.mod1 = paste(round(summary(mod.true.1)$coefficients[,1],2), " (",round(summary(mod.true.1)$coefficients[,2],2) ,")", sep="")

df2  = as.data.frame(cbind(names.mod2,paste(round(colMeans(LCM.model2[,,2]),2), " (", round(colMeans(LCM.model2.se[,,2]),2), ")",sep="" ),paste(round(colMeans(privBayes.model2[,,2]),2), " (", round(colMeans(privBayes.model2.se[,,2]),2), ")",sep="" ),paste(round(colMeans(GM.model2[,,2]),2), " (", round(colMeans(GM.model2.se[,,2]),2), ")",sep="" )))

names(df2) = c("names","V4","V5","V6")
df2$true.mod2 = paste(round(summary(mod.true.2)$coefficients[,1],2), " (",round(summary(mod.true.2)$coefficients[,2],2) ,")", sep="")


##Creating each model data frames for each epsilon (1.00)
df1 = as.data.frame(cbind(names.mod1,paste(round(colMeans(LCM.model1[,,3],na.rm=TRUE),2), " (", round(colMeans(LCM.model1.se[,,3],na.rm=TRUE),2), ")",sep="" ),paste(round(colMeans(privBayes.model1[,,3],na.rm=TRUE),2), " (", round(colMeans(privBayes.model1.se[,,3],na.rm=TRUE),2), ")",sep="" ),paste(round(colMeans(GM.model1[,,3], na.rm=TRUE),2), " (", round(colMeans(GM.model1.se[,,3], na.rm=TRUE),2), ")",sep="" )))
names(df1) = c("names","V1","V2","V3")
df1$true.mod1 = paste(round(summary(mod.true.1)$coefficients[,1],2), " (",round(summary(mod.true.1)$coefficients[,2],2) ,")", sep="")

df2  = as.data.frame(cbind(names.mod2,paste(round(colMeans(LCM.model2[,,3],na.rm=TRUE),2), " (", round(colMeans(LCM.model2.se[,,3],na.rm=TRUE),2), ")",sep="" ),paste(round(colMeans(privBayes.model2[,,3],na.rm=TRUE),2), " (", round(colMeans(privBayes.model2.se[,,3],na.rm=TRUE),2), ")",sep="" ),paste(round(colMeans(GM.model2[,,3],na.rm=TRUE),2), " (", round(colMeans(GM.model2.se[,,3],na.rm=TRUE),2), ")",sep="" )))

names(df2) = c("names","V4","V5","V6")
df2$true.mod2 = paste(round(summary(mod.true.2)$coefficients[,1],2), " (",round(summary(mod.true.2)$coefficients[,2],2) ,")", sep="")


##Repeating the stepwise selection
selected.terms = matrix(NA,nrow = 100, ncol = 15)
selected.terms.GM = array(NA, c(100, 15, length(epsilon)))
selected.terms.LCM = array(NA, c(100, 15, length(epsilon)))
selected.terms.PB = array(NA, c(100, 15, length(epsilon)))

for(q in 1:100){
  set.seed(i)
  ##Subsetting the data so it has a sample size = samp.size
  samps = sample(1:nrow(ACS.df),samp.size,replace=FALSE)
  ACS.df =  ACS.df[samps,]
  
  ##Subsetting the variables that we want
  ACS.df.tab = as.data.frame(table(ACS.df[,c(1,2,3,4,6)]))
  ##Since different functions aggregate the data differently, this is reordering the counts data to match the output of later functions
  trueData = plyr::count(ACS.df.tab)[c(1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31,2,18,10,26,6,22,14,30,4,20,12,28,8,24,16,32),-7]
  
  ##Fitting the full (saturated) log linear model on the orginal data
  log.mod = loglm(Freq~CIT+AGEP+RACWHT+SEX+PINCP+CIT*AGEP+CIT*RACWHT+CIT*SEX+CIT*PINCP+AGEP*RACWHT+AGEP*SEX+AGEP*PINCP+RACWHT*SEX+RACWHT*PINCP+SEX*PINCP,data=trueData)
  ##Using a stepwise approach to select important terms
  step.log.mod = step(log.mod)
  ##Renaming the terms
  names = c("CIT","AGEP","RACWHT","SEX","PINCP","CIT:AGEP","CIT:RACWHT","CIT:SEX","CIT:PINCP","AGEP:RACWHT","AGEP:SEX","AGEP:PINCP","RACWHT:SEX","RACWHT:PINCP","SEX:PINCP")
  
  ##Creating a list of terms selected by the procedure
  selected.real = lapply(step.log.mod,attributes)$term$term.labels
  selected.terms[i,] = as.integer(names %in% selected.real)
  
}

names = c("CIT","AGEP","RACWHT","SEX","PINCP","CIT:AGEP","CIT:RACWHT","CIT:SEX","CIT:PINCP","AGEP:RACWHT","AGEP:SEX","AGEP:PINCP","RACWHT:SEX","RACWHT:PINCP","SEX:PINCP")


for(i in 1:length(epsilon)){
  for(q in 1:100){
    a = i-1
    b=q-1
    file.name = paste(python.dir,"\\PrivBayes\\DataSynthesizer\\out\\correlated_attribute_mode\\privBayes_E",a,"_",b,".csv",sep="")
    dfpriv = read.csv(file.name)[,-1]
    names(dfpriv) = c("CIT","AGEP","RACWHT","SEX","PINCP")
    
    dfpriv = as.data.frame(table(dfpriv))
    
    log.mod = loglm(Freq~CIT+AGEP+RACWHT+SEX+PINCP+CIT*AGEP+CIT*RACWHT+CIT*SEX+CIT*PINCP+AGEP*RACWHT+AGEP*SEX+AGEP*PINCP+RACWHT*SEX+RACWHT*PINCP+SEX*PINCP,data=dfpriv)
    ##Using a stepwise approach to select important terms
    step.log.mod = step(log.mod)

    ##Creating a list of terms selected by the procedure
    selected.real = lapply(step.log.mod,attributes)$term$term.labels
    selected.terms.PB[q,,i] = as.integer(names %in% selected.real)
    
    ##GM
    ep.names = c("0_25","0_50","1_00")
    file.name = paste(python.dir,"\\Graphical Models Example\\GM_",ep.names[i],"_",q,".csv",sep="")
    dfgm = unlist(read.csv(file.name,header=FALSE))
    df.gm = plyr::count(ACS.df.tab)
    df.gm$Freq = dfgm
    
    log.mod = loglm(Freq~CIT+AGEP+RACWHT+SEX+PINCP+CIT*AGEP+CIT*RACWHT+CIT*SEX+CIT*PINCP+AGEP*RACWHT+AGEP*SEX+AGEP*PINCP+RACWHT*SEX+RACWHT*PINCP+SEX*PINCP,data=df.gm)
    ##Using a stepwise approach to select important terms
    step.log.mod = step(log.mod)
    
    ##Creating a list of terms selected by the procedure
    selected.real = lapply(step.log.mod,attributes)$term$term.labels
    selected.terms.GM[q,,i] = as.integer(names %in% selected.real)
    
    ##LCM
    dflcm = read.csv(paste(results.dir,"\\BayesLCM_E",i,"_",q,".csv",sep=""))[,-1]
    dflcm = colMeans(dflcm)*samp.size
    df.lcm = trueData
    df.lcm$Freq = round(dflcm)
      
    log.mod = loglm(Freq~CIT+AGEP+RACWHT+SEX+PINCP+CIT*AGEP+CIT*RACWHT+CIT*SEX+CIT*PINCP+AGEP*RACWHT+AGEP*SEX+AGEP*PINCP+RACWHT*SEX+RACWHT*PINCP+SEX*PINCP,data=df.lcm)
    ##Using a stepwise approach to select important terms
    step.log.mod = step(log.mod)
    
    ##Creating a list of terms selected by the procedure
    selected.real = lapply(step.log.mod,attributes)$term$term.labels
    selected.terms.LCM[q,,i] = as.integer(names %in% selected.real)
    
  }
}

