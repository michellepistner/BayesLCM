library(data.table)
library(plyr)
library(MCMCpack)
library(foreign)
library(matrixStats)
library(Rfast)
library(parallel)

print("Loaded the libraries...")
getProbs = function(comb, cols, samp, G.used, class.probs){
  class.probs[cbind(comb,samp,G.used, cols)]
}

row.prob.vec <- function(combs.full, mix.probs, class.probs, P, sample, G){
  iter.vec = rep(1:P, nrow(combs.full))
  samp = rep(sample, nrow(combs.full)*ncol(combs.full))
  combs = c(t(combs.full))
  tmp.probs =  sapply(1:G, FUN = getProbs, comb = combs, col = iter.vec, samp = samp, class.probs = class.probs)
  tmp.probs = rowprods(matrix(tmp.probs, ncol = P, byrow = TRUE))
  prob = sum(tmp.probs * rep(mix.probs[sample,],each = length(tmp.probs)/G))
  return(prob)
}


row.prob <- function(combs, mix.probs, class.probs, P, sample, G){
    ##We need a fixed value for val1 and val2 and iterate over the combinations
    iter.vec = 1:P
    samp = rep(sample, length(combs))
    tmp.probs = rowProds(t(sapply(1:G, FUN = getProbs, comb = combs, col = iter.vec, samp = samp)))
    ##tmp.probs = apply(apply(cbind(combs,iter.vec),1,FUN = function(info.vec, class.probs, sample){class.probs[info.vec[1], sample, , info.vec[2]]}, class.probs = class.probs, sample = sample), 1, prod)
    prob = sum(mix.probs[sample,]*tmp.probs)
    return(prob)
}


twoWay.probs.vectorized <- function(counts.vec, mix.probs, class.probs, index1=NULL, index2=NULL, val1=NULL, val2 = NULL, sample = 1, G){
  if(is.null(index1)){
    inds = which(!is.na(counts.vec))
    index1 = inds[1]
    index2 = inds[2]
    val1 = counts.vec[inds[1]]
    val2 = counts.vec[inds[2]]
  }
  
  ##Find all possible combinations
  P = dim(class.probs)[4]
  
  ###First, we need to find all combinations to sum over
  p.opt = seq(1,P)[!(seq(1,P) %in% c(index1,index2))]
  

  K = dim(class.probs)[1]
  
  tmp.list = rep(list(seq(1,K)), length(p.opt))
  combs = expand.grid(tmp.list)
  
  names(combs) = p.opt
  vec1 = rep(val1, nrow(combs))
  vec2 = rep(val2, nrow(combs))
  
  ##Splice in these vectors at the correct place
  combs.full = matrix(NA, nrow = nrow(combs), ncol = ncol(combs) + 2)
  combs.full[,index1] = vec1
  combs.full[,index2] = vec2
  tmp.counter = 1
  for(i in 1:ncol(combs.full)){
    if(is.na(combs.full[1,i])){
      combs.full[,i] = combs[,tmp.counter]
      tmp.counter = tmp.counter + 1}
  }
  combs.full = data.table(combs.full)
  prob = row.prob.vec(combs.full, mix.probs, class.probs, P, sample, G)

  return(prob)
}#end of function


fullTab.probs <- function(mix.probs, class.probs, vals, sample = NULL){
  ##Vals: a P dim vector
  if(length(vals) != dim(class.probs[[1]])[3]){
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
      prob.vec = prob.vec * class.probs[[vals[i]]][sample,,i]
    }
    prob = sum(mix.probs[sample,] * prob.vec)
  }
  
  return(prob)
}#End of function


lap.lik <- function(data,prop,E=-1/6){
  val = 0
  alpha = exp(E)
  for(i in 1:length(data)){
    val = val + abs(data[i] - prop[i])*log(alpha)
  }
  return(val)
}#end of function

dT2Lap <- function(noisyS,S,samp.size,E){
  alpha = exp(E)
  pmf = 0
  for(i in 1:length(noisyS)){
    if(noisyS[i] == 0){
      pmf = pmf + log(alpha^abs(noisyS[i] - S[i])/(1+alpha))
    }else if(noisyS[i] == samp.size){
      pmf = pmf + log(alpha^(samp.size-S[i])/(1+alpha))
    }else{
      pmf = pmf + log((1-alpha)*alpha^abs(noisyS[i] - S[i])/(1+alpha))
    }
  }
  return(pmf)
}

likelihood <- function(counts,samp.size, G, sample, counts.mat, mix.probs, comp.probs, index = NULL){
  probs = rep(NA, length(counts))
  
  if(is.null(index)){
    for(y in 1:nrow(counts.mat)){
      inds = which(!is.na(counts.mat[y,]))
      val1 = counts.mat[y, inds[1]]
      val2 = counts.mat[y, inds[2]]
      probs[y] = twoWay.probs(mix.probs, comp.probs, index1 = inds[1], index2 = inds[2], val1, val2, sample = sample)
    }##End of for
  } else{
    inds = which(!is.na(counts.mat[index,]))
    val1 = counts.mat[index, inds[1]]
    val2 = counts.mat[index, inds[2]]
    probs = twoWay.probs(mix.probs, comp.probs, index1 = inds[1], index2 = inds[2], val1, val2, sample = sample)
    
  }
  lik = sum(dbinom(counts,samp.size, probs,log=T))
  return(lik)
}#end of function

likelihood.vec <- function(counts,samp.size, G, sample, counts.mat, mix.probs, comp.probs, index = NULL, par = FALSE, clust = cl){
  probs = rep(NA, length(counts))
  
  if(is.null(index)){
    if(par == TRUE){
      probs = parApply(cl = clust, counts.mat,MARGIN = 1,FUN = twoWay.probs.vectorized, mix.probs= mix.probs, class.probs = comp.probs, sample = sample, G=G)
    } else{
      for(y in 1:nrow(counts.mat)){
        inds = which(!is.na(counts.mat[y,]))
        val1 = counts.mat[y, inds[1]]
        val2 = counts.mat[y, inds[2]]
        probs[y] = twoWay.probs.vectorized(counts.vec = NULL, mix.probs, comp.probs, index1 = inds[1], index2 = inds[2], val1, val2, sample = sample, G=G)
      }##End of for
    }
  } else{
    inds = which(!is.na(counts.mat[index,]))
    val1 = counts.mat[index, inds[1]]
    val2 = counts.mat[index, inds[2]]
    probs = twoWay.probs.vectorized(counts.vec = NULL, mix.probs, comp.probs, index1 = inds[1], index2 = inds[2], val1, val2, sample = sample, G=G)
    
  }
  lik = sum(dbinom(counts,samp.size, probs,log=T))
  return(lik)
}#end of function


###Sampling function
mcmc.sampler <- function(counts, counts.mat, eps = 1, P, G = 5, nsamples = 2000, samp.size, E.var = NULL){
  
  ###1.) Create noisy counts
  Q = length(counts)/(length(unique(!is.na(c(comb.mat))))^2)
  E = -(eps/Q)/2
  if(is.null(E.var)){E.var = E}
  errs = rgeom(length(counts),1-exp(E)) - rgeom(length(counts),1-exp(E))
  noisyM = counts + errs
  
  
  ###Initialize the starting values
  ##Mixing probabilities
  mix.probs = matrix(NA, nrow=nsamples, ncol = G)
  mix.probs[1,] = rep(1/G,G) ##Starting value: equal weights
  
  ###Counts
  M = matrix(NA, nrow = nsamples, ncol = length(noisyM))
  M[1,] = noisyM + rgeom(length(counts),1-exp(E)) - rgeom(length(counts),1-exp(E))
  M[1,] = ifelse(M[1,] < 0, 0, M[1,])
  M[1,] = ifelse(M[1,] > samp.size, samp.size, M[1,])
  
  ###Psi
  ###Currently only implemented for two binary variables
  comp.probs = list()
  comp.probs[[1]] = array(NA, dim = c(nsamples, G, P))
  
  ##Random start
  start.val = cbind(comb.mat, noisyM/samp.size)
  tmp.mat = matrix(NA, nrow = G, ncol = P)
  for(i in 1:ncol(tmp.mat)){
    tmp.mat[,i] = mean(start.val[which(start.val[,i] == 1),ncol(tmp.mat) + 1])
  }
  comp.probs[[1]][1,,] = tmp.mat
  comp.probs[[2]] = 1 - comp.probs[[1]]
  
  ##Marginal probabilities
  marg.prob = matrix(NA, nsamples, length(noisyM))
  
  ##Full probabilities
  full.prob = matrix(NA, nsamples, 2^P)
  
  tmp.list = rep(list(seq(1,2)), P)
  full.mat = expand.grid(tmp.list)
  full.mat = as.matrix(full.mat)
  accept = matrix(0,nrow=G, ncol=P)
  accept.pi = 0
  accept.M = rep(0, length(noisyM))
  
  
  for(i in 2:nsamples){
    
    ##Sample for M
    ##Proposal
    for(q in 1:length(noisyM)){
      M.prop = M[i,]
      M.prop[q] = noisyM[q] + rgeom(1,1-exp(E.var)) - rgeom(1,1-exp(E.var))
      if(q != length(noisyM)){M.prop[(q+1):length(noisyM)] = M[(i-1),(q+1):length(noisyM)]}
      
      if(M.prop[q]< 0 | M.prop[q] > samp.size){
        r = -Inf
      }else{
        prop.post = lap.lik(M.prop[q],noisyM[q], E=E) + likelihood(M.prop[q], samp.size, G, i-1, counts.mat, mix.probs, comp.probs, index = q)
        past.post = lap.lik(M[(i-1),q], noisyM[q], E=E) + likelihood(M[(i-1),q],samp.size, G, i-1, counts.mat, mix.probs, comp.probs, index = q)
        ppropos = dT2Lap(M[(i-1),q],M.prop[q],samp.size,E) - dT2Lap(M.prop[q],M[(i-1),q],samp.size,E)
        r = exp(prop.post - past.post + ppropos)
      }
      if(runif(1)<r){
        M[i,q] = M.prop[q]
        accept.M[q] = accept.M[q] + 1
      }else{
        M[i,q] = M[(i-1),q]
      }
      print(q)
    }
    
    
    ##Sample for mixing probabilities
    
    prop.mix = mix.probs
    prop.mix[(i-1),] = mix.probs[(i-1),] + rnorm(G,0,sd = 0.005)
    prop.mix[(i-1),] = prop.mix[(i-1),]/sum(prop.mix[(i-1),])
    
    if(min(prop.mix[(i-1),]) < 0){
      r = -Inf
    }else{
      prop.post = likelihood(M[i,], samp.size, G, sample = i-1, counts.mat, mix.probs, comp.probs) + log(ddirichlet(mix.probs[(i-1),], alpha = rep(1/G,G)))
      past.post = likelihood(M[i,], samp.size, G, sample = i-1, counts.mat, prop.mix, comp.probs) + log(ddirichlet(prop.mix[(i-1),], alpha = rep(1/G,G)))
      ##ppropos = log(ddirichlet(mix.probs[(i-1),], alpha = prop.mix[(i-1),])) - log(ddirichlet(prop.mix[(i-1),], alpha = mix.probs[(i-1),]))
      ppropos = 0
      r = exp(prop.post - past.post + ppropos)
    }
    
    
    if(runif(1)<r){
      mix.probs[i,] = prop.mix[(i-1),]
      accept.pi = accept.pi + 1
    }else{
      mix.probs[i,] = mix.probs[(i-1),]
    }
    
    ##For each latent class and variable, sample the component probs
    for(k in 1:P){
      for(h in 1:G){
        
        ##tmp.prop.all = rdirichlet(1, alpha = 500*c(comp.probs[[1]][1,h,k], comp.probs[[2]][1,h,k]) + 1)#c(comp.probs[[1]][1,h,k], comp.probs[[2]][1,h,k]) + rnorm(2,0,0.1)
        tmp.prop.all = c(comp.probs[[1]][(i-1),h,k], comp.probs[[2]][(i-1),h,k]) + rnorm(2,0,0.01)
        tmp.prop.all = tmp.prop.all / sum(tmp.prop.all)
        if(min(tmp.prop.all) < 0 | max(tmp.prop.all) > 1){
          r = -Inf
        } else{
          tmp.prop = tmp.prop.all[1]
          prop.probs = list()
          prop.probs[[1]] = comp.probs[[1]][(i-1):i,,]
          prop.probs[[1]][2,,] = ifelse(is.na(prop.probs[[1]][2,,]), prop.probs[[1]][1,,], prop.probs[[1]][2,,])
          prop.probs[[2]] = 1-prop.probs[[1]]
          
          prop.tmp.probs = prop.probs
          prop.tmp.probs[[1]][2,h,k] = tmp.prop
          prop.tmp.probs[[2]] = 1 - prop.tmp.probs[[1]]
          
          mix.probs.hold = rbind(mix.probs[1,], mix.probs[i,])
          prop.post = likelihood(M[i,],samp.size,G, sample = 2, counts.mat, mix.probs.hold, prop.tmp.probs) + log(ddirichlet(c(tmp.prop, 1-tmp.prop),alpha = rep(.5,2)))
          past.post = likelihood(M[i,],samp.size,G, sample = 2, counts.mat, mix.probs.hold, prop.probs) + log(ddirichlet(c(comp.probs[[1]][i-1,h,k], comp.probs[[2]][i-1,h,k]),alpha = rep(.5,2)))
          ##ppropos = log(ddirichlet(c(comp.probs[[1]][i-1,h,k], comp.probs[[2]][i-1,h,k]),alpha = tmp.prop.all)) - log(ddirichlet(tmp.prop.all,alpha = c(comp.probs[[1]][i-1,h,k], comp.probs[[2]][i-1,h,k])))
          ppropos = 0
          r = exp(prop.post - past.post + ppropos)
        }
        if(runif(1)<r){
          comp.probs[[1]][i,h,k] = tmp.prop
          accept[h,k] = accept[h,k] + 1
        }else{
          comp.probs[[1]][i,h,k] = comp.probs[[1]][i-1,h,k]
        }
        comp.probs[[2]] = 1 - comp.probs[[1]]
      }##End of for
    }##End of for
    
    
    ###Calculate marginal probability for sample
    for(y in 1:nrow(counts.mat)){
      inds = which(!is.na(counts.mat[y,]))
      val1 = counts.mat[y, inds[1]]
      val2 = counts.mat[y, inds[2]]
      marg.prob[i,y] = twoWay.probs(mix.probs, comp.probs, index1 = inds[1], index2 = inds[2], val1, val2, sample = i)
    }##End of for
    
    ###Calculate full probability for sample
    #for(z in 1:nrow(full.mat)){
    #  full.prob[i,z] = fullTab.probs(mix.probs, comp.probs, vals = full.mat[z,], sample = i)
    #}
    print(i)
  }##End of for
  return(list("marg_probs" = marg.prob, "full_probs" = full.prob, "M" = M, "Psi" = comp.probs, "Pi" = mix.probs, "accept_rate_pi" = accept.pi, "accept_rate_psi" = accept, "accept_rate_M" = accept.M))
}##End of function

###Sampling function
mcmc.sampler.v2 <- function(counts, counts.mat, eps = 1, P, G = 5, nsamples = 2000, samp.size, E.var = NULL){
  
  ###1.) Create noisy counts
  Q = length(counts)/(length(unique(!is.na(c(comb.mat))))^2)
  E = -(eps/Q)/2
  if(is.null(E.var)){E.var = E}
  errs = rgeom(length(counts),1-exp(E)) - rgeom(length(counts),1-exp(E))
  noisyM = counts + errs
  
  
  ###Initialize the starting values
  ##Mixing probabilities
  mix.probs = matrix(NA, nrow=nsamples, ncol = G)
  mix.probs[1,] = rep(1/G,G) ##Starting value: equal weights
  
  ###Counts
  M = matrix(NA, nrow = nsamples, ncol = length(noisyM))
  M[1,] = noisyM + rgeom(length(counts),1-exp(E)) - rgeom(length(counts),1-exp(E))
  M[1,] = ifelse(M[1,] < 0, 0, M[1,])
  M[1,] = ifelse(M[1,] > samp.size, samp.size, M[1,])
  
  ###Psi
  ###Currently only implemented for two binary variables
  comp.probs = array(NA, dim = c(2,nsamples,G,P))

  ##Random start
  start.val = cbind(comb.mat, noisyM/samp.size)
  tmp.mat = matrix(NA, nrow = G, ncol = P)
  for(i in 1:ncol(tmp.mat)){
    tmp.mat[,i] = mean(start.val[which(start.val[,i] == 1),ncol(tmp.mat) + 1])
  }
  comp.probs[1,1,,] = tmp.mat
  comp.probs[2,,,] = 1 - comp.probs[1,,,]
  
  ##Marginal probabilities
  marg.prob = matrix(NA, nsamples, length(noisyM))
  
  ##Full probabilities
  full.prob = matrix(NA, nsamples, 2^P)
  
  tmp.list = rep(list(seq(1,2)), P)
  full.mat = expand.grid(tmp.list)
  full.mat = as.matrix(full.mat)
  accept = matrix(0,nrow=G, ncol=P)
  accept.pi = 0
  accept.M = rep(0, length(noisyM))
  
  
  for(i in 2:nsamples){
    
    ##Sample for M
    ##Proposal
    for(q in 1:length(noisyM)){
      M.prop = M[i,]
      M.prop[q] = noisyM[q] + rgeom(1,1-exp(E.var)) - rgeom(1,1-exp(E.var))
      if(q != length(noisyM)){M.prop[(q+1):length(noisyM)] = M[(i-1),(q+1):length(noisyM)]}
      
      if(M.prop[q]< 0 | M.prop[q] > samp.size){
        r = -Inf
      }else{
        prop.post = lap.lik(M.prop[q],noisyM[q], E=E) + likelihood.vec(M.prop[q], samp.size, G, i-1, counts.mat, mix.probs, comp.probs, index = q)
        if(q == 1){
          past.post = lap.lik(M[(i-1),q], noisyM[q], E=E) + likelihood.vec(M[(i-1),q],samp.size, G, i-1, counts.mat, mix.probs, comp.probs, index = q)
        }
        ppropos = dT2Lap(M[(i-1),q],M.prop[q],samp.size,E) - dT2Lap(M.prop[q],M[(i-1),q],samp.size,E)
        r = exp(prop.post - past.post + ppropos)
      }
      if(runif(1)<r){
        M[i,q] = M.prop[q]
        accept.M[q] = accept.M[q] + 1
        past.post = prop.post
      }else{
        M[i,q] = M[(i-1),q]
        past.post = past.post
      }
      print(q)
    }

    
    ##Sample for mixing probabilities
    
    prop.mix = mix.probs
    prop.mix[(i-1),] = mix.probs[(i-1),] + rnorm(G,0,sd = 0.005)
    prop.mix[(i-1),] = prop.mix[(i-1),]/sum(prop.mix[(i-1),])
    
    if(min(prop.mix[(i-1),]) < 0){
      r = -Inf
    }else{
      prop.post = likelihood.vec(M[i,], samp.size, G, sample = i-1, counts.mat, mix.probs, comp.probs, par = TRUE) + log(ddirichlet(mix.probs[(i-1),], alpha = rep(1/G,G)))
      past.post = likelihood.vec(M[i,], samp.size, G, sample = i-1, counts.mat, prop.mix, comp.probs) + log(ddirichlet(prop.mix[(i-1),], alpha = rep(1/G,G)))
      ##ppropos = log(ddirichlet(mix.probs[(i-1),], alpha = prop.mix[(i-1),])) - log(ddirichlet(prop.mix[(i-1),], alpha = mix.probs[(i-1),]))
      ppropos = 0
      r = exp(prop.post - past.post + ppropos)
    }

    
    
    if(runif(1)<r){
      mix.probs[i,] = prop.mix[(i-1),]
      accept.pi = accept.pi + 1
      past.post = prop.post
    }else{
      mix.probs[i,] = mix.probs[(i-1),]
      past.post = past.post
    }
    
    ##For each latent class and variable, sample the component probs
    for(k in 1:P){
      for(h in 1:G){
        
        ##tmp.prop.all = rdirichlet(1, alpha = 500*c(comp.probs[[1]][1,h,k], comp.probs[[2]][1,h,k]) + 1)#c(comp.probs[[1]][1,h,k], comp.probs[[2]][1,h,k]) + rnorm(2,0,0.1)
        tmp.prop.all = c(comp.probs[1,(i-1),h,k], comp.probs[2,(i-1),h,k]) + rnorm(2,0,0.01)
        tmp.prop.all = tmp.prop.all / sum(tmp.prop.all)
        if(min(tmp.prop.all) < 0 | max(tmp.prop.all) > 1){
          r = -Inf
        } else{
          tmp.prop = tmp.prop.all[1]
          prop.probs = comp.probs[,(i-1):i,,]
          prop.probs[1,2,,] = ifelse(is.na(prop.probs[1,2,,]), prop.probs[1,1,,], prop.probs[1,2,,])
          prop.probs[2,,,] = 1-prop.probs[1,,,]
          
          prop.tmp.probs = prop.probs
          prop.tmp.probs[1,2,h,k] = tmp.prop
          prop.tmp.probs[2,,,] = 1 - prop.tmp.probs[1,,,]
          
          mix.probs.hold = rbind(mix.probs[1,], mix.probs[i,])
          start = Sys.time()
          prop.post = likelihood.vec(M[i,],samp.size,G, sample = 2, counts.mat, mix.probs = mix.probs.hold, comp.probs =prop.tmp.probs, par = TRUE) + log(ddirichlet(c(tmp.prop, 1-tmp.prop),alpha = rep(.5,2)))
          end = Sys.time()
          ##past.post = likelihood.vec(M[i,],samp.size,G, sample = 2, counts.mat, mix.probs = mix.probs.hold, comp.probs = prop.probs) + log(ddirichlet(c(comp.probs[1,i-1,h,k], comp.probs[2,i-1,h,k]),alpha = rep(.5,2)))
          ##ppropos = log(ddirichlet(c(comp.probs[[1]][i-1,h,k], comp.probs[[2]][i-1,h,k]),alpha = tmp.prop.all)) - log(ddirichlet(tmp.prop.all,alpha = c(comp.probs[[1]][i-1,h,k], comp.probs[[2]][i-1,h,k])))
          ppropos = 0
          r = exp(prop.post - past.post + ppropos)
        }
        if(runif(1)<r){
          comp.probs[1,i,h,k] = tmp.prop
          accept[h,k] = accept[h,k] + 1
          past.post = prop.post
        }else{
          comp.probs[1,i,h,k] = comp.probs[1,i-1,h,k]
          past.post = past.post
        }
        comp.probs[2,,,] = 1 - comp.probs[1,,,]
      }##End of for
    }##End of for
    
    
    ###Calculate marginal probability for sample
    marg.prob[i,] = parApply(cl = cl, counts.mat,MARGIN = 1,FUN = twoWay.probs.vectorized, mix.probs= mix.probs, class.probs = comp.probs, sample = i, G=G)

    
    ###Calculate full probability for sample
    #for(z in 1:nrow(full.mat)){
    #  full.prob[i,z] = fullTab.probs(mix.probs, comp.probs, vals = full.mat[z,], sample = i)
    #}
    print(i)
  }##End of for
  return(list("marg_probs" = marg.prob, "M" = M, "Psi" = comp.probs, "Pi" = mix.probs, "accept_rate_pi" = accept.pi, "accept_rate_psi" = accept, "accept_rate_M" = accept.M))
}##End of function

mcmc.sampler.v3 <- function(counts, counts.mat, eps = 1, P, G = 5, nsamples = 2000, samp.size, E.var = NULL, results.dir, .methodPi = "normal"){
  
  ###1.) Create noisy counts
  Q = length(counts)/(length(unique(!is.na(c(comb.mat))))^2)
  E = -(eps/Q)/2
  if(is.null(E.var)){E.var = E}
  errs = rgeom(length(counts),1-exp(E)) - rgeom(length(counts),1-exp(E))
  noisyM = counts + errs
  
  
  ###Initialize the starting values
  ##Mixing probabilities
  load(paste0(results.dir, "/1_iter100b.Rdata"))
  mix.probs = matrix(NA, nrow=nsamples, ncol = G)
  mix.probs[1,] = samps$Pi[4000,] ##rdirichlet(1,rep(1/G,G)) ##Starting value: equal weights
  
  ###Counts
  M = matrix(NA, nrow = nsamples, ncol = length(noisyM))
  M[1,] = noisyM + rgeom(length(counts),1-exp(E)) - rgeom(length(counts),1-exp(E))
  M[1,] = ifelse(M[1,] < 0, 0, M[1,])
  M[1,] = ifelse(M[1,] > samp.size, samp.size, M[1,])
  
  ###Psi
  ###Currently only implemented for two binary variables
  comp.probs = array(NA, dim = c(2,nsamples,G,P))
  
  ##Random start
  #start.val = cbind(comb.mat, noisyM/samp.size, noisyM)
  #tmp.mat = matrix(NA, nrow = G, ncol = P)
  #for(i in 1:ncol(tmp.mat)){
  #  nrow.starts = nrow(start.val[which(start.val[,i] == 1),])
  #  tmp.mat[,i] = sum(start.val[which(start.val[,i] == 1),ncol(tmp.mat) + 2])/(samp.size*nrow.starts/2)
  #
  #}
  comp.probs[1,1,,] = samps$Psi[1,4000,,]
  comp.probs[2,,,] = 1 - comp.probs[1,,,]
  ##Marginal probabilities
  marg.prob = matrix(NA, nsamples, length(noisyM))
  
  ##Full probabilities
  full.prob = matrix(NA, nsamples, 2^P)
  
  tmp.list = rep(list(seq(1,2)), P)
  full.mat = expand.grid(tmp.list)
  full.mat = as.matrix(full.mat)
  accept = rep(0,G)
  accept.pi = 0
  accept.M = rep(0, length(noisyM))
  
  
  for(i in 2:nsamples){
    
    ##Sample for M
    ##Proposal
    for(q in 1:length(noisyM)){
      M.prop = M[(i-1),]
      M.prop[q] = noisyM[q] + rgeom(1,1-exp(E.var)) - rgeom(1,1-exp(E.var))
      if(q != length(noisyM)){M.prop[(q+1):length(noisyM)] = M[(i-1),(q+1):length(noisyM)]}
      
      if(M.prop[q]< 0 | M.prop[q] > samp.size){
        r = -Inf
      }else{
        bin.lik.new = likelihood.vec(M.prop[q], samp.size, G, i-1, counts.mat, mix.probs, comp.probs, index = q)
        prop.post = lap.lik(M.prop[q],noisyM[q], E=E) + bin.lik.new
        bin.lik= likelihood.vec(M[(i-1),q],samp.size, G, i-1, counts.mat, mix.probs, comp.probs, index = q)

        past.post = lap.lik(M[(i-1),q], noisyM[q], E=E) + bin.lik
        ppropos = dT2Lap(M[(i-1),q],M.prop[q],samp.size,E) - dT2Lap(M.prop[q],M[(i-1),q],samp.size,E)
        r = exp(prop.post - past.post + ppropos)
      }
      if(runif(1)<r){
        M[i,q] = M.prop[q]
        accept.M[q] = accept.M[q] + 1
      } else{
        M[i,q] = M[(i-1),q]
      }
    }
    
    
    ##Sample for mixing probabilities
    if(.methodPi == "dirichlet"){
      prop.mix = mix.probs
      prop.mix[(i-1),] = rdirichlet(1, alpha = 50 * mix.probs[(i-1),] + .5) ###50 * mix.probs[(i-1),] + 1
      
      if(min(prop.mix[(i-1),]) < 0){
        r = -Inf
      }else{
        prop.post = likelihood.vec(M[i,], samp.size, G, sample = i-1, counts.mat, prop.mix, comp.probs, par = TRUE) + log(ddirichlet(prop.mix[(i-1),], alpha = rep(1/G,G)))
        past.post = likelihood.vec(M[i,], samp.size, G, sample = i-1, counts.mat, mix.probs, comp.probs, par = TRUE) + log(ddirichlet(mix.probs[(i-1),], alpha = rep(1/G,G)))
        ppropos = log(ddirichlet(mix.probs[(i-1),], alpha = 50 * prop.mix[(i-1),] + .5)) - log(ddirichlet(prop.mix[(i-1),], alpha = 50 * mix.probs[(i-1),] +.5))
        ##ppropos = log(ddirichlet(mix.probs[(i-1),], alpha = rep(1/G, G))) - log(ddirichlet(prop.mix[(i-1),], alpha = rep(1/G,G)))
        
        r = exp(prop.post - past.post + ppropos)
      }
      
      
      if(runif(1)<r){
        mix.probs[i,] = prop.mix[(i-1),]
        accept.pi = accept.pi + 1
      }else{
        mix.probs[i,] = mix.probs[(i-1),]
      }
    } else{
      prop.mix = mix.probs
      eta = mix.probs[(i-1),] + rnorm(G,0,sd = 0.01)
      prop.mix[(i-1),] = eta
      prop.mix[(i-1),] = prop.mix[(i-1),]/sum(prop.mix[(i-1),])
      
      if(min(prop.mix[(i-1),]) < 0){
        r = -Inf
      }else{
        prop.post = likelihood.vec(M[i,], samp.size, G, sample = i-1, counts.mat, mix.probs, comp.probs, cl=cl) + sum(dgamma(mix.probs[(i-1),], shape = 1/G, scale = 1, log = TRUE))
        past.post = likelihood.vec(M[i,], samp.size, G, sample = i-1, counts.mat, prop.mix, comp.probs, cl=cl) + sum(dgamma(prop.mix[(i-1),], shape = 1/G, scale = 1, log = TRUE))
        ppropos = 0
        r = exp(prop.post - past.post + ppropos)
      }
      
      
      if(runif(1)<r){
        mix.probs[i,] = prop.mix[(i-1),]
        accept.pi = accept.pi + 1
      }else{
        mix.probs[i,] = mix.probs[(i-1),]
      }
      
    }
    
    ##For each latent class and variable, sample the component probs
    for(h in 1:G){
      tmp.prop.all = matrix(NA, nrow = 2, ncol = P)
      for(k in 1:P){
	  tmp.prop.all[,k] = rdirichlet(1, alpha = 50 * c(comp.probs[1,(i-1),h,k], comp.probs[2,(i-1),h,k])+ .25)
      }

      ##tmp.prop.all = rdirichlet(1, alpha = 5*c(comp.probs[[1]][1,h,k], comp.probs[[2]][1,h,k]) + 1)#c(comp.probs[[1]][1,h,k], comp.probs[[2]][1,h,k]) + rnorm(2,0,0.1)
      #tmp.prop.all = matrix(c(comp.probs[1,(i-1),h,], comp.probs[2,(i-1),h,]), ncol = P, byrow=TRUE) + rnorm(2*P,0,0.01)
      #tmp.prop.all = apply(tmp.prop.all,2,FUN = function(vec){vec/sum(vec)})
      if(min(tmp.prop.all) < 0 | max(tmp.prop.all) > 1){
        r = -Inf
      } else{
        prop.probs = comp.probs[,(i-1):i,,]
        prop.probs[1,2,,] = ifelse(is.na(prop.probs[1,2,,]), prop.probs[1,1,,], prop.probs[1,2,,])
        prop.probs[2,,,] = 1-prop.probs[1,,,]
        
        prop.tmp.probs = prop.probs
        prop.tmp.probs[1,2,h,] = tmp.prop.all[1,]
        prop.tmp.probs[2,,,] = 1 - prop.tmp.probs[1,,,]
        
        mix.probs.hold = rbind(mix.probs[1,], mix.probs[i,])
        
        if(h == 1){
          bin.lik = likelihood.vec(M[i,],samp.size,G, sample = 2, counts.mat, mix.probs = mix.probs.hold, comp.probs =prop.probs, par = TRUE)
        }
        bin.lik.new = likelihood.vec(M[i,],samp.size,G, sample = 2, counts.mat, mix.probs = mix.probs.hold, comp.probs =prop.tmp.probs, par = TRUE)
        prop.post = bin.lik.new + sum(log(ddirichlet(t(tmp.prop.all),alpha = rep(.5,2))))
        past.post = bin.lik + sum(log(ddirichlet(cbind(comp.probs[1,i-1,h,], comp.probs[2,i-1,h,]),alpha = rep(.5,2))))
	ppropos = 0
	for(k in 1:P){
        	ppropos =  ppropos + log(ddirichlet(c(comp.probs[1,(i-1),h,k], comp.probs[2,(i-1),h,k]),alpha = 50 *tmp.prop.all[,k] + .25)) - log(ddirichlet(tmp.prop.all[,k],alpha = 50 *c(comp.probs[1,(i-1),h,k], comp.probs[2,(i-1),h,k]) + .25))
	}
        r = exp(prop.post - past.post + ppropos)
      }
      if(runif(1)<r){
        comp.probs[1,i,h,] = tmp.prop.all[1,]
        accept[h] = accept[h] + 1
        bin.lik = bin.lik.new
      }else{
        comp.probs[1,i,h,] = comp.probs[1,i-1,h,]
      }
      comp.probs[2,,,] = 1 - comp.probs[1,,,]
    }##End of for
    
    
    ###Calculate marginal probability for sample
    marg.prob[i,] = parApply(cl = cl, counts.mat,MARGIN = 1,FUN = twoWay.probs.vectorized, mix.probs= mix.probs, class.probs = comp.probs, sample = i, G=G)
    
    ###Calculate full probability for sample
    #for(z in 1:nrow(full.mat)){
    #  full.prob[i,z] = fullTab.probs(mix.probs, comp.probs, vals = full.mat[z,], sample = i)
    #}
    print(i)
  }##End of for
  return(list("marg_probs" = marg.prob, "full_probs" = full.prob, "M" = M, "Psi" = comp.probs, "Pi" = mix.probs, "accept_rate_pi" = accept.pi, "accept_rate_psi" = accept, "accept_rate_M" = accept.M))
}##End of function


print("Read in the functions...")

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
##Run the algorithm

eps = 1*.75
reps = 2
cl = makeCluster(detectCores()-1)
clusterExport(cl, varlist = c("twoWay.probs.vectorized", "row.prob.vec","getProbs", "data.table", "rowprods", "P"))

for(zz in 100){
  set.seed(zz)
  start = Sys.time()
  samps = mcmc.sampler.v3(freq, comb.mat, P=P,nsamples = 150, eps = eps, samp.size = nrow(df.samp), G=7, results.dir = results.dir)
  end = Sys.time()
  print(end-start)
  save(samps, file = paste0(results.dir, "/", eps, "_iter",zz,"test.Rdata"))
}
stopCluster(cl)
