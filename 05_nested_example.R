##Coding for the nested example

source(config.R)
library(MCMCpack,lib.loc = packages.dir)




Summary.Probs <- function(G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3){
  summary.probs = rep(0,5)
  sum1 = rep(0,G)
  sum3a = matrix(0,nrow=9,ncol=G)
  sum3b = rep(0,9)
  sum4 = rep(0,G)
  sum5a = matrix(0,nrow=9,ncol=G)
  sum5b = rep(0,9)
  sum5c = rep(0,9)
  for(a in 1:G){ ##Looping over all the latent classes (household level)
    summary.probs[2] = summary.probs[2] + pi[a]*lambda1[1,a]
    for(b in 1:M){##Looping over all the latent classes (individual level)
      sum1[a] = sum1[a] + w[a,b]*psi2[[a]][1,b]
      sum4[a]= sum4[a] + w[a,b]*psi3[[a]][1,b]
    }
    summary.probs[1] = summary.probs[1] + pi[a]*lambda2[1,a]*sum1[a]^2
    summary.probs[4] = summary.probs[4] + 2*pi[a]*sum4[a] - (pi[a] *sum4[a]^2)
  }
  for(q in 1:9){
    for(a in 1:G){
      for(b in 1:M){
        sum3a[q,a] = sum3a[q,a] + w[a,b]*psi2[[a]][q,b]
        sum5a[q,a] = sum5a[q,a] + w[a,b]*psi2[[a]][q,b]*psi3[[a]][1,b]
      }
      sum3b[q] = sum3b[q] + pi[a]*lambda2[q,a] * sum3a[q,a]^2
      sum5b[q] = sum5b[q] + pi[a]*lambda2[q,a]* sum5a[q,a]
      sum5c[q] = sum5c[q] + pi[a] *lambda2[q,a]*sum5a[q,a]^2
    }
    summary.probs[3] = summary.probs[3] + sum3b[q]
    summary.probs[5] = summary.probs[5] +  2*sum5b[q] - sum5c[q]
  }
  return(summary.probs)
}

likelihood <- function(counts,samp.size, G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3){
  probs = Summary.Probs(G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3)
  lik = sum(dbinom(counts,samp.size, probs,log=T))
  return(lik)
}#end of function

lap.lik <- function(data,prop,E=-1/6){
  val = 0
  alpha = exp(E)
  for(i in 1:length(data)){
    val = val + log(1-alpha) - log(1+alpha) + abs(data[i] - prop[i])*log(alpha)
  }
  return(val)
}#end of function

TwoWay_NoisyLCM <- function(counts, G, M, samp.size, N.runs,E,prop.weights){
  ##Initializing vectors to hold the class membership probabilities
  ##Each is a matrix-- "c*" represents variable * where .1 or .2 represent the levels
  ##Each column represents a latent class
  sum.probs = matrix(NA,nrow=N.runs,ncol=5)
  
  pi = rep(NA,G)
  w = matrix(NA, nrow = G, ncol = M)
  
  lambda1 = matrix(NA,nrow = 2, ncol = G)
  lambda2 = matrix(NA,nrow = 9, ncol = G)

  psi1 = lapply(1:G, matrix, data= NA, nrow=2, ncol=M)
  psi2 = lapply(1:G, matrix, data= NA, nrow=9, ncol=M)
  psi3 = lapply(1:G, matrix, data= NA, nrow=12, ncol=M)
  
  S = counts
  
  
  ##Now, setting inital values
  pi = rdirichlet(1,alpha = rep(1,G))
  
  for(a in 1:G){
    w[a,] = rdirichlet(1,alpha = rep(1,M))
    lambda1[,a] = rdirichlet(1,alpha = c(1,1))
    lambda2[,a] = rdirichlet(1,alpha = 5*c(9:1))
    for(b in 1:M){
      psi1[[a]][,b] = rdirichlet(1,alpha = c(1,1))
      psi2[[a]][,b] = rdirichlet(1,alpha = 5*c(9:1))
      psi3[[a]][,b] = rdirichlet(1,alpha = 5*c(12:1))
    }
  }
  ##Now, beging the M-H algorithm
  
  for(j in 2:N.runs){
    k = j-1
    
    ##Mixing proportions
    prop.mix = rdirichlet(1,alpha = 5*pi+.5)
    
    
    prop.post = likelihood(S,samp.size,G, M, prop.mix, lambda1, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(prop.mix,alpha = rep(1/G,G)))
    past.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(pi,alpha = rep(1/G,G)))
    ppropos = log(ddirichlet(pi,alpha = 5*pi+.5)) - log(ddirichlet(prop.mix,alpha = 5*pi+.5))
    
    r = exp(prop.post - past.post + ppropos)
    r = ifelse(is.na(r),-Inf,r)
    
    
    if(runif(1)<r){
      pi = prop.mix
    }else{
      pi = pi
    }
    
    ##Now, for lambda1
    for(a in 1:G){
      tmp.prop = rdirichlet(1,alpha = 5*lambda1[,a]+.5)
      lambda1.prop = lambda1
      lambda1.prop[,a] = tmp.prop
      
      prop.post = likelihood(S,samp.size,G, M, pi, lambda1.prop, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(tmp.prop,alpha = rep(1/2,2)))
      past.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(lambda1[,a],alpha = rep(1/2,2)))
      ppropos = log(ddirichlet(lambda1[,a],alpha = 5*lambda1[,a]+.5)) - log(ddirichlet(tmp.prop,alpha = 5*lambda1[,a]+.5))
      
      r = exp(prop.post - past.post + ppropos)
      r = ifelse(is.na(r),-Inf,r)
      
      
      if(runif(1)<r){
        lambda1 = lambda1.prop
      }
    }
    
    ##Lambda2
    for(a in 1:G){
      tmp.prop = rdirichlet(1,5*lambda2[,a]+.5)
      lambda2.prop = lambda2
      lambda2.prop[,a] = tmp.prop
      
      prop.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2.prop, w, psi1, psi2, psi3) + log(ddirichlet(tmp.prop,alpha = rep(1/2,9)))
      past.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(lambda2[,a],alpha = rep(1/2,9)))
      ppropos = log(ddirichlet(lambda2[,a],alpha = 5*lambda2[,a]+.5)) - log(ddirichlet(tmp.prop,alpha = 5*lambda2[,a]+.5))
      
      r = exp(prop.post - past.post + ppropos)
      r = ifelse(is.na(r),-Inf,r)
      
      
      if(runif(1)<r){
        lambda2 = lambda2.prop
      }
    }

  ##The mixing parameter w
    for(a in 1:G){
      tmp.prop = rdirichlet(1,alpha = 5*w[a,]+.5)
      w.prop = w
      w.prop[a,] = tmp.prop
      
      prop.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w.prop, psi1, psi2, psi3) + log(ddirichlet(tmp.prop,alpha = rep(1/2,M)))
      past.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(w[a,],alpha = rep(1/2,M)))
      ppropos = log(ddirichlet(w[a,],alpha = 5*w[a,]+.5)) - log(ddirichlet(tmp.prop,alpha = 5*w[a,]+.5))
      
      r = exp(prop.post - past.post + ppropos)
      r = ifelse(is.na(r),-Inf,r)
      
      
      if(runif(1)<r){
        w = w.prop
      }
    }
    
    ##Now for psi1
    for(a in 1:G){
      for(b in 1:M){
        tmp.prop = rdirichlet(1,alpha = 5*psi1[[a]][,b]+.5)
        psi1.prop = psi1
        psi1.prop[[a]][,b] = tmp.prop
        
        prop.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1.prop, psi2, psi3) + log(ddirichlet(tmp.prop,alpha = rep(1/2,2)))
        past.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(psi1[[a]][,b],alpha = rep(1/2,2)))
        ppropos = log(ddirichlet(psi1[[a]][,b],alpha = 5*psi1[[a]][,b]+.5)) - log(ddirichlet(tmp.prop,alpha = 5*psi1[[a]][,b]+.5))
        
        r = exp(prop.post - past.post + ppropos)
        r = ifelse(is.na(r),-Inf,r)
        
        
        if(runif(1)<r){
          psi1 = psi1.prop
        }
      }
    }
    
    ##Now for psi2
    for(a in 1:G){
      for(b in 1:M){
        tmp.prop = rdirichlet(1,alpha = 5*psi2[[a]][,b]+.5)
        psi2.prop = psi2
        psi2.prop[[a]][,b] = tmp.prop
        
        prop.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2.prop, psi3) + log(ddirichlet(tmp.prop,alpha = rep(1/2,9)))
        past.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(psi2[[a]][,b],alpha = rep(1/2,9)))
        ppropos = log(ddirichlet(psi2[[a]][,b],alpha = 5*psi2[[a]][,b]+.5)) - log(ddirichlet(tmp.prop,alpha = 5*psi2[[a]][,b]+.5))
        
        r = exp(prop.post - past.post + ppropos)
        r = ifelse(is.na(r),-Inf,r)
        
        
        if(runif(1)<r){
          psi2 = psi2.prop
        }
      }
    }
    
    ##Now for psi3
    for(a in 1:G){
      for(b in 1:M){
        tmp.prop = rdirichlet(1,alpha = 5*psi3[[a]][,b]+.5)
        psi3.prop = psi3
        psi3.prop[[a]][,b] = tmp.prop
        
        prop.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3.prop) + log(ddirichlet(tmp.prop,alpha = rep(1/2,12)))
        past.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(psi3[[a]][,b],alpha = rep(1/2,12)))
        ppropos = log(ddirichlet(psi3[[a]][,b],alpha = 5*psi3[[a]][,b]+.5)) - log(ddirichlet(tmp.prop,alpha = 5*psi3[[a]][,b]+.5))
        
        r = exp(prop.post - past.post + ppropos)
        r = ifelse(is.na(r),-Inf,r)
        
        
        if(runif(1)<r){
          psi3 = psi3.prop
        }
      }
    }
    
    ##Now, the measurement error component
    S.tmp = S + rgeom(5,1-1/exp(E)) - rgeom(5,1-1/exp(E))
    
    prop.post = lap.lik(counts,S.tmp,E=-1/E) #+ dmultinom(M.prop,N,prob=rep(.25,4),log=TRUE)
    past.post = lap.lik(counts,S,E=-1/E) #+ dmultinom(M[(b:c),k],N,prob=rep(.25,4),log=TRUE)
      
      r = exp(prop.post - past.post)
      
      if(runif(1)<r){
        S = S.tmp
      }
      
    sum.probs[j,] = Summary.Probs(G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3)
    print(j)
    }#end of for
  return(sum.probs)
}

TwoWay_LCM <- function(counts, G, M, samp.size, N.runs,E,prop.weights){
  ##Initializing vectors to hold the class membership probabilities
  ##Each is a matrix-- "c*" represents variable * where .1 or .2 represent the levels
  ##Each column represents a latent class
  sum.probs = matrix(NA,nrow=N.runs,ncol=5)
  
  pi = rep(NA,G)
  w = matrix(NA, nrow = G, ncol = M)
  
  lambda1 = matrix(NA,nrow = 2, ncol = G)
  lambda2 = matrix(NA,nrow = 9, ncol = G)
  
  psi1 = lapply(1:G, matrix, data= NA, nrow=2, ncol=M)
  psi2 = lapply(1:G, matrix, data= NA, nrow=9, ncol=M)
  psi3 = lapply(1:G, matrix, data= NA, nrow=12, ncol=M)
  
  S = counts
  
  
  ##Now, setting inital values
  pi = rdirichlet(1,alpha = rep(1,G))
  
  for(a in 1:G){
    w[a,] = rdirichlet(1,alpha = rep(1,M))
    lambda1[,a] = rdirichlet(1,alpha = c(1,1))
    lambda2[,a] = rdirichlet(1,alpha = 5*c(9:1))
    for(b in 1:M){
      psi1[[a]][,b] = rdirichlet(1,alpha = c(1,1))
      psi2[[a]][,b] = rdirichlet(1,alpha = 5*c(9:1))
      psi3[[a]][,b] = rdirichlet(1,alpha = 5*c(12:1))
    }
  }
  ##Now, beging the M-H algorithm
  
  for(j in 2:N.runs){
    k = j-1
    
    ##Mixing proportions
    prop.mix = rdirichlet(1,alpha = 5*pi+.5)
    
    
    prop.post = likelihood(S,samp.size,G, M, prop.mix, lambda1, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(prop.mix,alpha = rep(1/G,G)))
    past.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(pi,alpha = rep(1/G,G)))
    ppropos = log(ddirichlet(pi,alpha = 5*pi+.5)) - log(ddirichlet(prop.mix,alpha = 5*pi+.5))
    
    r = exp(prop.post - past.post + ppropos)
    r = ifelse(is.na(r),-Inf,r)
    
    
    if(runif(1)<r){
      pi = prop.mix
    }else{
      pi = pi
    }
    
    ##Now, for lambda1
    for(a in 1:G){
      tmp.prop = rdirichlet(1,alpha = 5*lambda1[,a]+.5)
      lambda1.prop = lambda1
      lambda1.prop[,a] = tmp.prop
      
      prop.post = likelihood(S,samp.size,G, M, pi, lambda1.prop, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(tmp.prop,alpha = rep(1/2,2)))
      past.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(lambda1[,a],alpha = rep(1/2,2)))
      ppropos = log(ddirichlet(lambda1[,a],alpha = 5*lambda1[,a]+.5)) - log(ddirichlet(tmp.prop,alpha = 5*lambda1[,a]+.5))
      
      r = exp(prop.post - past.post + ppropos)
      r = ifelse(is.na(r),-Inf,r)
      
      
      if(runif(1)<r){
        lambda1 = lambda1.prop
      }
    }
    
    ##Lambda2
    for(a in 1:G){
      tmp.prop = rdirichlet(1,5*lambda2[,a]+.5)
      lambda2.prop = lambda2
      lambda2.prop[,a] = tmp.prop
      
      prop.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2.prop, w, psi1, psi2, psi3) + log(ddirichlet(tmp.prop,alpha = rep(1/2,9)))
      past.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(lambda2[,a],alpha = rep(1/2,9)))
      ppropos = log(ddirichlet(lambda2[,a],alpha = 5*lambda2[,a]+.5)) - log(ddirichlet(tmp.prop,alpha = 5*lambda2[,a]+.5))
      
      r = exp(prop.post - past.post + ppropos)
      r = ifelse(is.na(r),-Inf,r)
      
      
      if(runif(1)<r){
        lambda2 = lambda2.prop
      }
    }
    
    ##The mixing parameter w
    for(a in 1:G){
      tmp.prop = rdirichlet(1,alpha = 5*w[a,]+.5)
      w.prop = w
      w.prop[a,] = tmp.prop
      
      prop.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w.prop, psi1, psi2, psi3) + log(ddirichlet(tmp.prop,alpha = rep(1/2,M)))
      past.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(w[a,],alpha = rep(1/2,M)))
      ppropos = log(ddirichlet(w[a,],alpha = 5*w[a,]+.5)) - log(ddirichlet(tmp.prop,alpha = 5*w[a,]+.5))
      
      r = exp(prop.post - past.post + ppropos)
      r = ifelse(is.na(r),-Inf,r)
      
      
      if(runif(1)<r){
        w = w.prop
      }
    }
    
    ##Now for psi1
    for(a in 1:G){
      for(b in 1:M){
        tmp.prop = rdirichlet(1,alpha = 5*psi1[[a]][,b]+.5)
        psi1.prop = psi1
        psi1.prop[[a]][,b] = tmp.prop
        
        prop.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1.prop, psi2, psi3) + log(ddirichlet(tmp.prop,alpha = rep(1/2,2)))
        past.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(psi1[[a]][,b],alpha = rep(1/2,2)))
        ppropos = log(ddirichlet(psi1[[a]][,b],alpha = 5*psi1[[a]][,b]+.5)) - log(ddirichlet(tmp.prop,alpha = 5*psi1[[a]][,b]+.5))
        
        r = exp(prop.post - past.post + ppropos)
        r = ifelse(is.na(r),-Inf,r)
        
        
        if(runif(1)<r){
          psi1 = psi1.prop
        }
      }
    }
    
    ##Now for psi2
    for(a in 1:G){
      for(b in 1:M){
        tmp.prop = rdirichlet(1,alpha = 5*psi2[[a]][,b]+.5)
        psi2.prop = psi2
        psi2.prop[[a]][,b] = tmp.prop
        
        prop.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2.prop, psi3) + log(ddirichlet(tmp.prop,alpha = rep(1/2,9)))
        past.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(psi2[[a]][,b],alpha = rep(1/2,9)))
        ppropos = log(ddirichlet(psi2[[a]][,b],alpha = 5*psi2[[a]][,b]+.5)) - log(ddirichlet(tmp.prop,alpha = 5*psi2[[a]][,b]+.5))
        
        r = exp(prop.post - past.post + ppropos)
        r = ifelse(is.na(r),-Inf,r)
        
        
        if(runif(1)<r){
          psi2 = psi2.prop
        }
      }
    }
    
    ##Now for psi3
    for(a in 1:G){
      for(b in 1:M){
        tmp.prop = rdirichlet(1,alpha = 5*psi3[[a]][,b]+.5)
        psi3.prop = psi3
        psi3.prop[[a]][,b] = tmp.prop
        
        prop.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3.prop) + log(ddirichlet(tmp.prop,alpha = rep(1/2,12)))
        past.post = likelihood(S,samp.size,G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3) + log(ddirichlet(psi3[[a]][,b],alpha = rep(1/2,12)))
        ppropos = log(ddirichlet(psi3[[a]][,b],alpha = 5*psi3[[a]][,b]+.5)) - log(ddirichlet(tmp.prop,alpha = 5*psi3[[a]][,b]+.5))
        
        r = exp(prop.post - past.post + ppropos)
        r = ifelse(is.na(r),-Inf,r)
        
        
        if(runif(1)<r){
          psi3 = psi3.prop
        }
      }
    }
    
    
    sum.probs[j,] = Summary.Probs(G, M, pi, lambda1, lambda2, w, psi1, psi2, psi3)
    print(j)
  }#end of for
  return(sum.probs)
}

setwd("C:\\Users\\m_pis\\Box Sync\\Summer Project at Duke (Pistner)\\Master Code Files")
##Reading in the data
df = read.csv("origdata.txt", sep=" ")
S = rep(NA,5)

##Subselect households of size 3
tabs = table(df$HHIndex)
vals = labels(tabs[tabs==3])
vals = as.numeric(unlist(vals))
df.subset = df[df$HHIndex %in% vals,]

##Convert to wide format
df.tmp = reshape(df.subset,idvar="HHIndex",timevar = "WithinHHIndex",direction = "wide")

##Everyone is white
df.tmp$same = ifelse((df.tmp$Race.1 == 1) & (df.tmp$Race.2 == 1)  & (df.tmp$Race.3 == 1),1,0)
S[1] = sum(df.tmp$same)

##Ownership
S[2] = sum(ifelse(df.tmp$Owner.1==1,1,0))

##Everyone has the same race
df.tmp$same = ifelse((df.tmp$Race.1 == df.tmp$Race.2) & (df.tmp$Race.1 == df.tmp$Race.3),1,0)
S[3] = sum(df.tmp$same)

##Spouse present
df.tmp$SP = ifelse((df.tmp$Relate.2 == 2) | (df.tmp$Relate.3 == 2),1,0)
S[4] = sum(df.tmp$SP)

##Same race CP
df.tmp$SRCP = ifelse(((df.tmp$Relate.2 == 2) & (df.tmp$Race.1 == df.tmp$Race.2)) | ((df.tmp$Relate.3 == 2) & (df.tmp$Race.1 == df.tmp$Race.3)),1,0)
S[5] = sum(df.tmp$SRCP)

####Measuring without noise

set.seed(1)
probs = TwoWay_NoisyLCM(S, G=2, M=2, samp.size, N.runs=20000,E)
probs.reduced.inds = seq(10000,20000,by=5)
probs.reduced = probs[probs.reduced.inds,]

dat.sum = data.frame(Predicted = colMeans(probs.reduced),True = S/samp.size)

dodge <- position_dodge(width = 0)
ggplot(dat.sum, aes(x=True, y=Predicted)) + 
  geom_point() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=16),legend.position = "none") + 
  geom_abline(slope = 1, intercept = 0) + xlab("True Probability") + ylab("Estimated Probability") +scale_x_continuous(limits = c(.6, 1)) +scale_y_continuous(limits = c(.6, 1))

ggsave(paste(results.dir,"nested_noNoise.pdf",sep="/"))

###Now, adding privacy

E =  c(0.01,.1,1)/5
dat.sum = matrix(NA,nrow = length(E)*5,ncol = 3)

set.seed(1)
for(q in 1:length(E)){
  errs = rgeom(5,1-1/exp(E[q])) - rgeom(5,1-1/exp(E[q]))
  noisyS = S + errs
  noisyS = ifelse(noisyS<0,0,noisyS)
  noisyS = ifelse(noisyS>samp.size,samp.size,noisyS)
  counts = noisyS
  
  probs = TwoWay_NoisyLCM(counts, G=2, M=2, samp.size, N.runs=20000,E[q])
  
  probs.reduced.inds = seq(10000,20000,by=5)
  probs.reduced = probs[probs.reduced.inds,]
  dat.sum[(5*q-4):(5*q),1] = colMeans(probs.reduced)
  dat.sum[(5*q-4):(5*q),3] = rep(E[q]*5,5)
}

dat.sum[,2] = rep(S/samp.size,3)
dat.sum = data.frame(dat.sum)
names(dat.sum) = c("Predicted","True","Epsilon")
dat.sum$Epsilon = factor(dat.sum$Epsilon)

dodge <- position_dodge(width = 0)
ggplot(dat.sum, aes(x=True, y=Predicted, color = Epsilon)) + 
  geom_point() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=16)) + 
  geom_abline(slope = 1, intercept = 0) +  xlab("True Probability")+ ylab("Estimated Probability") +scale_x_continuous(limits = c(.6, 1)) +scale_y_continuous(limits = c(.6, 1))

ggsave(paste(results.dir,"nested_privacy.pdf",sep="/"))

