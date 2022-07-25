library(data.table)
library(plyr)
library(dplyr)
library(MCMCpack)
library(foreign)
library(matrixStats)
library(Rfast)
library(foreach)
library(doParallel)
library(parallel)
library(ggplot2)

base.dir = getwd()

data.dir = paste(base.dir,"/data",sep="")
results.dir = paste(base.dir,"/out",sep="")
packages.dir = paste(base.dir,"/packages",sep="")



df = read.dta(paste(data.dir, "da09681-0006.dta", sep="/"))

print("Reading in the data...")

df.samp =  df[,c("SCN_15_A_Y04", "SCN_15_B_Y04", "SCN_15_C_Y04", "SCN_15_D_Y04", "SCN_15_E_Y04", "SCN_15_F_Y04", "SCN_15_G_Y04", "SCN_15_H_Y04", "SCN_15_I_Y04", "SCN_17_A_Y04","SCN_17_B_Y04", "SCN_17_C_Y04", "SCN_17_D_Y04", "SCN_17_E_Y04", "SCN_17_F_Y04", "SCN_17_G_Y04")]
df.samp = df.samp[which(df.samp$SCN_15_A_Y04 %in% c("Yes", "No")),]
df.samp = df.samp[which(df.samp$SCN_15_B_Y04 %in% c("Yes", "No")),]
df.samp = df.samp[which(df.samp$SCN_15_C_Y04 %in% c("Yes", "No")),]
df.samp = df.samp[which(df.samp$SCN_15_D_Y04 %in% c("Yes", "No")),]
df.samp = df.samp[which(df.samp$SCN_15_E_Y04 %in% c("Yes", "No")),]
df.samp = df.samp[which(df.samp$SCN_15_F_Y04 %in% c("Yes", "No")),]
df.samp = df.samp[which(df.samp$SCN_15_G_Y04 %in% c("Yes", "No")),]
df.samp = df.samp[which(df.samp$SCN_15_H_Y04 %in% c("Yes", "No")),]
df.samp = df.samp[which(df.samp$SCN_15_I_Y04 %in% c("Yes", "No")),]
df.samp = df.samp[which(df.samp$SCN_17_A_Y04 %in% c("Yes", "No")),]
df.samp = df.samp[which(df.samp$SCN_17_B_Y04 %in% c("Yes", "No")),]
df.samp = df.samp[which(df.samp$SCN_17_C_Y04 %in% c("Yes", "No")),]
df.samp = df.samp[which(df.samp$SCN_17_D_Y04 %in% c("Yes", "No")),]
df.samp = df.samp[which(df.samp$SCN_17_E_Y04 %in% c("Yes", "No")),]
df.samp = df.samp[which(df.samp$SCN_17_F_Y04 %in% c("Yes", "No")),]
df.samp = df.samp[which(df.samp$SCN_17_G_Y04 %in% c("Yes", "No")),]


##Two way margins 

##Construct the matrix
##Construct the margins -- 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
freq = c()
keep = c()
P = ncol(df.samp)
comb.mat = matrix(ncol = P)

pairs = fread("privbayes_pairs_e5_0.csv")
for(i in 1:(P-1)){
  for(j in (i+1):P){
    counts = plyr::count(df.samp[,c(i,j)])
    freq = c(freq, counts$freq)
    tmp.total = max(apply(pairs, MARGIN = 1, FUN = function(vec) sum(vec %in% names(counts)[1:2])))
    tmp.keep = ifelse(tmp.total == 2, TRUE, FALSE)
    keep = c(keep, rep(tmp.keep,4))
    tmp.mat = matrix(NA, nrow = 4, ncol = P)
    tmp.mat[,i] = counts[,1]
    tmp.mat[,j] = counts[,2]
    comb.mat = rbind(comb.mat, tmp.mat)
  }
}

comb.mat = comb.mat[-1,]
comb.mat = comb.mat[keep,]
freq = freq[keep]
comb.mat = comb.mat -5

P.vec = freq/nrow(df.samp)

load("5_iter100.Rdata")
samps1 = samps$marg_probs
load("5_iter100b.Rdata")
samps2 = samps$marg_probs
load("5_iter100c.Rdata")
samps3 = samps$marg_probs

marg_probs_5 = rbind(samps1,samps2,samps3) 

load("2_iter100.Rdata")
samps1 = samps$marg_probs
load("2_iter100b.Rdata")
samps2 = samps$marg_probs
load("2_iter100c.Rdata")
samps3 = samps$marg_probs

marg_probs_2 = rbind(samps1,samps2,samps3)

load("1_iter100.Rdata")
samps1 = samps$marg_probs
load("1_iter100b.Rdata")
samps2 = samps$marg_probs
load("1_iter100c.Rdata")
samps3 = samps$marg_probs

marg_probs_1 = rbind(samps1,samps2,samps3)

load("0.5_iter100.Rdata")
samps1 = samps$marg_probs
load("0.5_iter100b.Rdata")
samps2 = samps$marg_probs
load("0.5_iter100c.Rdata")
samps3 = samps$marg_probs

marg_probs_0.5 = rbind(samps1,samps2,samps3)

set.seed(100)
samps = sample(c(1:29),3)
qq[(4*samps[1]-3):(4*samps[1])]
qq[(4*samps[2]-3):(4*samps[2])]
qq[(4*samps[3]-3):(4*samps[3])]

P.vec1 = P.vec[(4*samps[1]-3):(4*samps[1])]
P.vec2 = P.vec[(4*samps[2]-3):(4*samps[2])]
P.vec3 = P.vec[(4*samps[3]-3):(4*samps[3])]

ci.lower = apply(na.omit(marg_probs_2[-c(1:2000),]), 2, quantile, prob = 0.025)
ci.upper = apply(na.omit(marg_probs_2[-c(1:2000),]), 2, quantile, prob = 0.975)

ci.lower[(4*samps[1]-3):(4*samps[1])]
ci.upper[(4*samps[1]-3):(4*samps[1])]

ci.lower[(4*samps[2]-3):(4*samps[2])]
ci.upper[(4*samps[2]-3):(4*samps[2])]

ci.lower[(4*samps[3]-3):(4*samps[3])]
ci.upper[(4*samps[3]-3):(4*samps[3])]

Pnoisy[(4*samps[1]-3):(4*samps[1])]
Pnoisy[(4*samps[2]-3):(4*samps[2])]
Pnoisy[(4*samps[3]-3):(4*samps[3])]



avgProb_5 = apply(marg_probs_5[-c(1:2000),],2,mean, na.rm = TRUE)
avgProb_2 = apply(marg_probs_2[-c(1:2000),],2,mean, na.rm = TRUE)
avgProb_1 = apply(marg_probs_1[-c(1:2000),],2,mean, na.rm = TRUE)
avgProb_0.5 = apply(marg_probs_0.5[-c(1:2000),],2,mean, na.rm = TRUE)

Q = length(P.vec)
graph.df = data.frame("Estimated" = c(avgProb_0.5, avgProb_1, avgProb_2, avgProb_5), "True" = rep(P.vec,4), "Epsilon" = c(rep("0.5",Q),rep("1",Q),rep("2",Q), rep("5",Q)))

graph.df = graph.df %>% group_by(Epsilon) %>% arrange(Epsilon, True) %>%
  ungroup() %>%
  mutate(LogDiff = log(Estimated/True)) %>%
  mutate(Marginal = rep(1:Q,4))

ggplot(graph.df, aes(x = Marginal, y = LogDiff, color = Epsilon, shape = Epsilon)) +
  geom_point(size = 2) +
  ylab("Log Difference of Estimated & True") + 
  xlab("Marginal Count") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=14)) +
  geom_hline(yintercept = 0, linetype = 'dashed')


ggsave("nltcs_graphic.pdf", height = 4, width = 5, units = "in")



avgProb_5 = apply(marg_probs_5[-c(1:2000),],2,mean, na.rm = TRUE)
avgProb_2 = apply(marg_probs_2[-c(1:2000),],2,mean, na.rm = TRUE)
avgProb_1 = apply(marg_probs_1[-c(1:2000),],2,mean, na.rm = TRUE)
avgProb_0.5 = apply(marg_probs_0.5[-c(1:2000),],2,mean, na.rm = TRUE)

lowerProb_5 = apply(marg_probs_5[-c(1:2000),],2,quantile, prob = 0.025, na.rm = TRUE)
lowerProb_2 = apply(marg_probs_2[-c(1:2000),],2,quantile, prob = 0.025, na.rm = TRUE)
lowerProb_1 = apply(marg_probs_1[-c(1:2000),],2,quantile, prob = 0.025, na.rm = TRUE)
lowerProb_0.5 = apply(marg_probs_0.5[-c(1:2000),],2,quantile, prob = 0.025, na.rm = TRUE)

upperProb_5 = apply(marg_probs_5[-c(1:2000),],2,quantile, prob = 0.975, na.rm = TRUE)
upperProb_2 = apply(marg_probs_2[-c(1:2000),],2,quantile, prob = 0.975, na.rm = TRUE)
upperProb_1 = apply(marg_probs_1[-c(1:2000),],2,quantile, prob = 0.975, na.rm = TRUE)
upperProb_0.5 = apply(marg_probs_0.5[-c(1:2000),],2,quantile, prob = 0.975, na.rm = TRUE)

Q = length(P.vec)
graph.df = data.frame("Estimated" = c(avgProb_0.5, avgProb_1, avgProb_2, avgProb_5),"Lower" = c(lowerProb_0.5, lowerProb_1, lowerProb_2, lowerProb_5), "Upper" = c(upperProb_0.5, upperProb_1, upperProb_2, upperProb_5), "True" = rep(P.vec,4), "Epsilon" = c(rep("0.5",Q),rep("1",Q),rep("2",Q), rep("5",Q)))

ggplot(graph.df, aes(x = True, y = Estimated, color = Epsilon, shape = Epsilon)) +
  geom_point(size = 2) +
  ylab("Estimated Probability") + 
  xlab("True Probability") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=14)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

ggplot(graph.df, aes(x = True, y = Estimated, color = Epsilon, shape = Epsilon)) +
  geom_point(size = 2) +
  ylab("Estimated Probability") + 
  xlab("True Probability") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=14)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.05) 


###Synthetic data example
fullTab.probs <- function(mix.probs, class.probs, vals, sample = NULL){
  ##Vals: a P dim vector
  if(length(vals) != dim(class.probs)[4]){
    return("Evaluated combination is of the wrong length!!")
  }
  
  ##Now, need to find the prob of each combination
  ##Easy way to do it: first find the prob of the combination for each latent class
  ##Then, multiply by the mixing probabilities
  
  if(is.null(sample)){
    ##Take the mean and compute
    mix.avg = apply(mix.probs, 2, "mean")
    class.avg = list()
    for(i in 1:length(class.probs)){
      class.avg[[i]] = apply(class.probs[[i]], c(2,3), "mean")
    }
    
    prob.vec = rep(1, length(mix.avg))
    for(i in 1:length(vals)){
      prob.vec = prob.vec * class.avg[[vals[i]]][,i]
    }
    prob = sum(mix.avg * prob.vec)
  } else{
    prob.vec = rep(1, ncol(mix.probs))
    for(i in 1:length(vals)){
      prob.vec = prob.vec * class.probs[as.numeric(vals[i]), sample,,i]
    }
    prob = sum(mix.probs[sample,] * prob.vec)
  }
  
  return(prob)
}#End of function

##For each, randomly pick 20 runs between 2,000 and 12,000
d = 20

set.seed(2022)
###Epsilon = 2
data_pts = sample(2000:12000,d)

full_probs_e2 = matrix(NA,nrow = 2^16, ncol = d)
full_mat = expand.grid("V1" = c(1,2), "V2" = c(1,2), "V3" = c(1,2), "V4" = c(1,2), "V5" = c(1,2), "V6" = c(1,2), "V7" = c(1,2), "V8" = c(1,2), "V9" = c(1,2), "V10" = c(1,2), "V11" = c(1,2), "V12" = c(1,2), "V13" = c(1,2), "V14" = c(1,2), "V15" = c(1,2), "V16" = c(1,2))

for(i in 1:d){
  if(data_pts[i] < 4000){
    load("2_iter100.Rdata")
    for(z in 1:nrow(full_mat)){
      full_probs_e2[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = data_pts[i])
      if(z %% 1000 == 0) print(z)
    }
  } else if(data_pts[i] > 8000){
    load("2_iter100c.Rdata")
    for(z in 1:nrow(full_mat)){
      full_probs_e2[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = data_pts[i]-8000)
      if(z %% 1000 == 0) print(z)
    }
  } else{
    load("2_iter100b.Rdata")
    for(z in 1:nrow(full_mat)){
      full_probs_e2[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = data_pts[i]-4000)
      if(z %% 1000 == 0) print(z)
    }
  }
  print(paste0("The ", i, "th data set finished."))
}

est_2 = matrix(NA, nrow = d, ncol = 12)
var_2 = matrix(NA, nrow = d, ncol = 12)
set.seed(100)
samps_used = c(13,5,23)
for(j in 1:d){
  tmp_freq = round(nrow(df.samp) * full_probs_e2[,j])
  tmp_df_short = cbind(full_mat, tmp_freq) %>% filter(tmp_freq > 0)
  tmp_df = tmp_df_short[rep(row.names(tmp_df_short), tmp_df_short$tmp_freq), 1:16]
  names(tmp_df) = names(df.samp)
  
  for(q in 1:3){
    start = 4*q-3
    df_sub = tmp_df[,names(tmp_df) %in% c(pairs[samps_used[q],])]
    est_2[j,start:(start+3)] = as.data.frame(table(df_sub))$Freq/nrow(df_sub)
    
    hold_11 = ifelse(df_sub[,1] == 1 & df_sub[,2] == 1, 1, 0)
    hold_21 = ifelse(df_sub[,1] == 2 & df_sub[,2] == 1, 1, 0)
    hold_12 = ifelse(df_sub[,1] == 1 & df_sub[,2] == 2, 1, 0)
    hold_22 = ifelse(df_sub[,1] == 2 & df_sub[,2] == 2, 1, 0)
    
    var_2[j,start] = mean((hold_11-est_2[j,start])^2)/length(hold_11)
    var_2[j,(start+1)] = mean((hold_21-est_2[j,(start+1)])^2)/length(hold_21)
    var_2[j,(start+2)] = mean((hold_12-est_2[j,(start+2)])^2)/length(hold_12)
    var_2[j,(start+3)] = mean((hold_22-est_2[j,(start+3)])^2)/length(hold_22)
  }
  
}
within2 = apply(est_2, 2,var)/d

b2 = colMeans(est_2)
varFinal_2 = within2 + colMeans(var_2)


###Epsilon = 1
set.seed(200)
data_pts = sample(2000:12000,d)

full_probs_e1 = matrix(NA,nrow = 2^16, ncol = d)
full_mat = expand.grid("V1" = c(1,2), "V2" = c(1,2), "V3" = c(1,2), "V4" = c(1,2), "V5" = c(1,2), "V6" = c(1,2), "V7" = c(1,2), "V8" = c(1,2), "V9" = c(1,2), "V10" = c(1,2), "V11" = c(1,2), "V12" = c(1,2), "V13" = c(1,2), "V14" = c(1,2), "V15" = c(1,2), "V16" = c(1,2))

for(i in 1:d){
  if(data_pts[i] < 4000){
    load("1_iter100.Rdata")
    for(z in 1:nrow(full_mat)){
      full_probs_e1[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = data_pts[i])
      if(z %% 1000 == 0) print(z)
    }
  } else if(data_pts[i] > 8000){
    load("1_iter100c.Rdata")
    for(z in 1:nrow(full_mat)){
      full_probs_e1[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = data_pts[i]-8000)
      if(z %% 1000 == 0) print(z)
    }
  } else{
    load("1_iter100b.Rdata")
    for(z in 1:nrow(full_mat)){
      full_probs_e1[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = data_pts[i]-4000)
      if(z %% 1000 == 0) print(z)
    }
  }
  print(paste0("The ", i, "th data set finished."))
}

est_1 = matrix(NA, nrow = d, ncol = 12)
var_1 = matrix(NA, nrow = d, ncol = 12)
set.seed(100)
samps_used = c(13,5,23)
for(j in 1:d){
  tmp_freq = round(nrow(df.samp) * full_probs_e1[,j])
  tmp_df_short = cbind(full_mat, tmp_freq) %>% filter(tmp_freq > 0)
  tmp_df = tmp_df_short[rep(row.names(tmp_df_short), tmp_df_short$tmp_freq), 1:16]
  names(tmp_df) = names(df.samp)
  
  for(q in 1:3){
    start = 4*q-3
    df_sub = tmp_df[,names(tmp_df) %in% c(pairs[samps_used[q],])]
    est_1[j,start:(start+3)] = as.data.frame(table(df_sub))$Freq/nrow(df_sub)
    
    hold_11 = ifelse(df_sub[,1] == 1 & df_sub[,2] == 1, 1, 0)
    hold_21 = ifelse(df_sub[,1] == 2 & df_sub[,2] == 1, 1, 0)
    hold_12 = ifelse(df_sub[,1] == 1 & df_sub[,2] == 2, 1, 0)
    hold_22 = ifelse(df_sub[,1] == 2 & df_sub[,2] == 2, 1, 0)
    
    var_1[j,start] = mean((hold_11-est_1[j,start])^2)/length(hold_11)
    var_1[j,(start+1)] = mean((hold_21-est_1[j,(start+1)])^2)/length(hold_21)
    var_1[j,(start+2)] = mean((hold_12-est_1[j,(start+2)])^2)/length(hold_12)
    var_1[j,(start+3)] = mean((hold_22-est_1[j,(start+3)])^2)/length(hold_22)
  }
  
}
within1 = apply(est_1, 2,var)/d

b1 = colMeans(est_1)
varFinal_1 = within1 + colMeans(var_1)


###Epsilon = 0.5
set.seed(300)
data_pts = sample(2000:12000,d)

full_probs_e0_5 = matrix(NA,nrow = 2^16, ncol = d)
full_mat = expand.grid("V1" = c(1,2), "V2" = c(1,2), "V3" = c(1,2), "V4" = c(1,2), "V5" = c(1,2), "V6" = c(1,2), "V7" = c(1,2), "V8" = c(1,2), "V9" = c(1,2), "V10" = c(1,2), "V11" = c(1,2), "V12" = c(1,2), "V13" = c(1,2), "V14" = c(1,2), "V15" = c(1,2), "V16" = c(1,2))

for(i in 1:d){
  if(data_pts[i] < 4000){
    load("0.5_iter100.Rdata")
    for(z in 1:nrow(full_mat)){
      full_probs_e0_5[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = data_pts[i])
      if(z %% 1000 == 0) print(z)
    }
  } else if(data_pts[i] > 8000){
    load("0.5_iter100c.Rdata")
    for(z in 1:nrow(full_mat)){
      full_probs_e0_5[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = data_pts[i]-8000)
      if(z %% 1000 == 0) print(z)
    }
  } else{
    load("0.5_iter100b.Rdata")
    for(z in 1:nrow(full_mat)){
      full_probs_e0_5[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = data_pts[i]-4000)
      if(z %% 1000 == 0) print(z)
    }
  }
  print(paste0("The ", i, "th data set finished."))
}

est_0_5 = matrix(NA, nrow = d, ncol = 12)
var_0_5 = matrix(NA, nrow = d, ncol = 12)
set.seed(100)
samps_used = c(13,5,23)
for(j in 1:d){
  tmp_freq = round(nrow(df.samp) * full_probs_e0_5[,j])
  tmp_df_short = cbind(full_mat, tmp_freq) %>% filter(tmp_freq > 0)
  tmp_df = tmp_df_short[rep(row.names(tmp_df_short), tmp_df_short$tmp_freq), 1:16]
  names(tmp_df) = names(df.samp)
  
  for(q in 1:3){
    start = 4*q-3
    df_sub = tmp_df[,names(tmp_df) %in% c(pairs[samps_used[q],])]
    est_0_5[j,start:(start+3)] = as.data.frame(table(df_sub))$Freq/nrow(df_sub)
    
    hold_11 = ifelse(df_sub[,1] == 1 & df_sub[,2] == 1, 1, 0)
    hold_21 = ifelse(df_sub[,1] == 2 & df_sub[,2] == 1, 1, 0)
    hold_12 = ifelse(df_sub[,1] == 1 & df_sub[,2] == 2, 1, 0)
    hold_22 = ifelse(df_sub[,1] == 2 & df_sub[,2] == 2, 1, 0)
    
    var_0_5[j,start] = mean((hold_11-est_0_5[j,start])^2)/length(hold_11)
    var_0_5[j,(start+1)] = mean((hold_21-est_0_5[j,(start+1)])^2)/length(hold_21)
    var_0_5[j,(start+2)] = mean((hold_12-est_0_5[j,(start+2)])^2)/length(hold_12)
    var_0_5[j,(start+3)] = mean((hold_22-est_0_5[j,(start+3)])^2)/length(hold_22)
  }
  
}
within0_5 = apply(est_0_5, 2,var)/d

b0_5 = colMeans(est_0_5)
varFinal_0_5 = within0_5 + colMeans(var_0_5)

r_0_5 = (1/d)*within0_5/varFinal_0_5
dof_0_5 = (d-1)*(1+(1/r_0_5))^2

ci_upper0_5 = b0_5 + qt(.975,dof_0_5)*sqrt(varFinal_0_5)
ci_lower0_5 = b0_5 - qt(.975,dof_0_5)*sqrt(varFinal_0_5)

r_1 = (1/d)*within1/varFinal_1
dof_1 = (d-1)*(1+(1/r_1))^2

ci_upper1 = b1 + qt(.975,dof_1)*sqrt(varFinal_1)
ci_lower1 = b1 - qt(.975,dof_1)*sqrt(varFinal_1)

r_2 = (1/d)*within2/varFinal_2
dof_2 = (d-1)*(1+(1/r_2))^2
ci_upper2 = b2 + qt(.975,dof_2)*sqrt(varFinal_2)
ci_lower2 = b2 - qt(.975,dof_2)*sqrt(varFinal_2)

ci_0_5 = paste0("(", round(ci_lower0_5,3),", ",round(ci_upper0_5,3), ")")
ci_1 = paste0("(", round(ci_lower1,3),", ",round(ci_upper1,3), ")")
ci_2 = paste0("(", round(ci_lower2,3),", ",round(ci_upper2,3), ")")

library(stargazer)
stargazer(cbind(round(b0_5,3), ci_0_5, round(b1,3), ci_1, round(b2,3), ci_2))



###Calculating them all, just in case

full_probs_e2 = matrix(NA,nrow = 2^16, ncol = 10000)
full_mat = expand.grid("V1" = c(1,2), "V2" = c(1,2), "V3" = c(1,2), "V4" = c(1,2), "V5" = c(1,2), "V6" = c(1,2), "V7" = c(1,2), "V8" = c(1,2), "V9" = c(1,2), "V10" = c(1,2), "V11" = c(1,2), "V12" = c(1,2), "V13" = c(1,2), "V14" = c(1,2), "V15" = c(1,2), "V16" = c(1,2))

load("2_iter100.Rdata")
for(i in 2000:4000){
    for(z in 1:nrow(full_mat)){
      full_probs_e2[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = i)
      if(z %% 1000 == 0) print(z)
    }
  print(paste0("The ", i, "th data set finished."))
}

load("2_iter100b.Rdata")
for(i in 4002:8000){
  for(z in 1:nrow(full_mat)){
    samp_val = i - 4000
    full_probs_e2[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = samp_val)
    if(z %% 1000 == 0) print(z)
  }
  print(paste0("The ", i, "th data set finished."))
}

load("2_iter100c.Rdata")
for(i in 8002:12000){
  for(z in 1:nrow(full_mat)){
    samp_val = i - 8000
    full_probs_e2[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = samp_val)
    if(z %% 1000 == 0) print(z)
  }
  print(paste0("The ", i, "th data set finished."))
}

fwrite(full_probs_e2, "fullProbsE2.csv")
full_probs_e1 = matrix(NA,nrow = 2^16, ncol = 10000)
full_mat = expand.grid("V1" = c(1,2), "V2" = c(1,2), "V3" = c(1,2), "V4" = c(1,2), "V5" = c(1,2), "V6" = c(1,2), "V7" = c(1,2), "V8" = c(1,2), "V9" = c(1,2), "V10" = c(1,2), "V11" = c(1,2), "V12" = c(1,2), "V13" = c(1,2), "V14" = c(1,2), "V15" = c(1,2), "V16" = c(1,2))

load("1_iter100.Rdata")
for(i in 2000:4000){
  for(z in 1:nrow(full_mat)){
    full_probs_e1[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = i)
    if(z %% 1000 == 0) print(z)
  }
  print(paste0("The ", i, "th data set finished."))
}

load("1_iter100b.Rdata")
for(i in 4002:8000){
  for(z in 1:nrow(full_mat)){
    samp_val = i - 4000
    full_probs_e1[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = samp_val)
    if(z %% 1000 == 0) print(z)
  }
  print(paste0("The ", i, "th data set finished."))
}

load("1_iter100c.Rdata")
for(i in 8002:12000){
  for(z in 1:nrow(full_mat)){
    samp_val = i - 8000
    full_probs_e1[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = samp_val)
    if(z %% 1000 == 0) print(z)
  }
  print(paste0("The ", i, "th data set finished."))
}

fwrite(full_probs_e1, "fullProbsE1.csv")

full_probs_e0_5 = matrix(NA,nrow = 2^16, ncol = 10000)

load("0.5_iter100.Rdata")
for(i in 2000:4000){
  for(z in 1:nrow(full_mat)){
    full_probs_e0_5[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = i)
    if(z %% 1000 == 0) print(z)
  }
  print(paste0("The ", i, "th data set finished."))
}

load("0.5_iter100b.Rdata")
for(i in 4002:8000){
  for(z in 1:nrow(full_mat)){
    samp_val = i - 4000
    full_probs_e0_5[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = samp_val)
    if(z %% 1000 == 0) print(z)
  }
  print(paste0("The ", i, "th data set finished."))
}

load("0.5_iter100c.Rdata")
for(i in 8002:12000){
  for(z in 1:nrow(full_mat)){
    samp_val = i - 8000
    full_probs_e0_5[z,i] = fullTab.probs(samps$Pi, samps$Psi, full_mat[z,], sample = samp_val)
    if(z %% 1000 == 0) print(z)
  }
  print(paste0("The ", i, "th data set finished."))
}

fwrite(full_probs_e0_5, "fullProbsE0_5.csv")
