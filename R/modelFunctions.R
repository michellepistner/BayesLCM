#' Calculate the likelihood for the discrete Laplace distribution
#'
#' @param data A vector of counts
#' @param prop A vector representing the mean of the counts
#' @param E The transformed epsilon to match the geometric distribution
#'
#' @return The value of the likelihood evaulated at `data`
lap.lik <- function(data,prop,E=-1/6){
  val = 0
  alpha = exp(E)
  for(i in 1:length(data)){
    val = val + abs(data[i] - prop[i])*log(alpha)
  }
  return(val)
}#end of function

#' Density function for the two-sided geometric distribution
#'
#' @param noisyS A count vector representing the new draw
#' @param S A count vector denoting the mean of the proposal distribution
#' @param samp.size Sample size
#' @param E The transformed epsilon to match the geometric distribution
#' @return The density function evaluated at `noisyS`
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

#' The binomial likelihood function for a specific count and probability vector
#'
#' @param counts A vector of noisy counts
#' @param samp.size Sample size
#' @param G The number of latent classes
#' @param sample What sample in the chain should the likelihood be evaulated for?
#' @param counts.mat A matrix encoding what counts correspond to what variables
#' @param mix.probs A # of samples * G matrix of mixing probabilities
#' @param comp.probs An array of the component probabilities for each latent class and run
#' @param index Defaults to NULL. If a specific value, will only calculate the likelihood corresponding to the count in that index.
#' @param par Boolean. Do we parallelize the code?
#' @param clust If paralleized, what cluster do we use
#' @return The evaluated likelihood
likelihood.vec <- function(counts,samp.size, G, sample, counts.mat, mix.probs, comp.probs, index = NULL, par = FALSE, clust = cl){
  probs = rep(NA, length(counts))

  if(is.null(index)){
    if(par == TRUE){
      probs = parallel::parApply(cl = clust, counts.mat,MARGIN = 1,FUN = twoWay.probs.vectorized, mix.probs= mix.probs, class.probs = comp.probs, sample = sample, G=G)
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
  lik = sum(stats::dbinom(counts,samp.size, probs,log=T))
  return(lik)
}#end of function


#' The main sampling function for the Bayesian Latent Class approach
#'
#' @param counts A vector of the true counts
#' @param counts.mat A matrix encoding what counts correspond to what variables
#' @param eps The value for epsilon
#' @param P The number of variables
#' @param G The number of latent classes
#' @param nsamples The number of samples to be drawn
#' @param samp.size Sample size
#' @param E.var Defaults to NULL. A value can be supplied to change the tuning parameter when sampling for \eqn{\Tilde{M}}
#' @param .calculateFullTabProbs A Boolean. Do we calculate the full table probabilities? It can be expensive to do so when dimension or the number of latent classes is large.
#' @param .PiTuningParam A small constant used to tune the proposal distribution of \eqn{\pi}.
#' @param .PsiTuningParam A small constant used to tune the proposal distribution of \eqn{\psi}.
#' @param cl Cluster used if parallelized
#' @return A list containing samples of \eqn{\Tilde{M}}, \eqn{\pi}, and \eqn{\psi} as well as marginal probabilities and full table probabilities (if desired). Also contains acceptance rates for all of the parameters.
#' @export
mcmc.sampler <- function(counts, counts.mat, eps = 1, P, G = 5, nsamples = 2000, samp.size, E.var = NULL, .calculateFullTabProbs = FALSE, .PiTuningParam = 0.5, .PsiTuningParam = 10, .methodPi = "normal", cl = NULL){
  par = ifelse(is.null(cl), FALSE, TRUE)
  ###1.) Create noisy counts
  Q = length(counts)/(length(unique(!is.na(c(counts.mat))))^2)
  E = -(eps/Q)/2
  if(is.null(E.var)){E.var = E}
  errs = stats::rgeom(length(counts),1-exp(E)) - stats::rgeom(length(counts),1-exp(E))
  noisyM = counts + errs


  ###Initialize the starting values
  ##Mixing probabilities
  mix.probs = matrix(NA, nrow=nsamples, ncol = G)
  mix.probs[1,] = rep(1/G,G) ##Starting value: equal weights

  ###Counts
  M = matrix(NA, nrow = nsamples, ncol = length(noisyM))
  M[1,] = noisyM + stats::rgeom(length(counts),1-exp(E)) - stats::rgeom(length(counts),1-exp(E))
  M[1,] = ifelse(M[1,] < 0, 0, M[1,])
  M[1,] = ifelse(M[1,] > samp.size, samp.size, M[1,])

  ###Psi
  ###Currently only implemented for two binary variables
  comp.probs = array(NA, dim = c(2,nsamples,G,P))

  ##Random start
  start.val = cbind(counts.mat, noisyM/samp.size, noisyM)
  tmp.mat = matrix(NA, nrow = G, ncol = P)
  for(i in 1:ncol(tmp.mat)){
    nrow.starts = nrow(start.val[which(start.val[,i] == 1),])
    tmp.mat[,i] = sum(start.val[which(start.val[,i] == 1),ncol(tmp.mat) + 2])/(samp.size*nrow.starts/2)

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
  accept = rep(0,G)
  accept.pi = 0
  accept.M = rep(0, length(noisyM))

  if(par){
    marg.prob[1,] = parallel::parApply(cl = cl, counts.mat,MARGIN = 1,FUN = twoWay.probs.vectorized, mix.probs= mix.probs, class.probs = comp.probs, sample = 1, G=G)
  } else{
    marg.prob[1,] = apply(counts.mat,MARGIN = 1,FUN = twoWay.probs.vectorized, mix.probs= mix.probs, class.probs = comp.probs, sample = 1, G=G)

  }

  for(i in 2:nsamples){

    ##Sample for M
    ##Proposal
    for(q in 1:length(noisyM)){
      M.prop = M[(i-1),]
      M.prop[q] = noisyM[q] + stats::rgeom(1,1-exp(E.var)) - stats::rgeom(1,1-exp(E.var))
      if(q != length(noisyM)){M.prop[(q+1):length(noisyM)] = M[(i-1),(q+1):length(noisyM)]}

      if(M.prop[q]< 0 | M.prop[q] > samp.size){
        r = -Inf
      }else{
        bin.lik.new = likelihood.vec(M.prop[q], samp.size, G, i-1, counts.mat, mix.probs, comp.probs, index = q, clust = cl, par = par)
        prop.post = lap.lik(M.prop[q],noisyM[q], E=E) + bin.lik.new
        bin.lik= likelihood.vec(M[(i-1),q],samp.size, G, i-1, counts.mat, mix.probs, comp.probs, index = q, clust = cl, par = par)

        past.post = lap.lik(M[(i-1),q], noisyM[q], E=E) + bin.lik
        ppropos = dT2Lap(M[(i-1),q],M.prop[q],samp.size,E) - dT2Lap(M.prop[q],M[(i-1),q],samp.size,E)
        r = exp(prop.post - past.post + ppropos)
      }
      if(stats::runif(1)<r){
        M[i,q] = M.prop[q]
        accept.M[q] = accept.M[q] + 1
      } else{
        M[i,q] = M[(i-1),q]
      }
    }


    ##Sample for mixing probabilities
    if(.methodPi == "dirichlet"){
      prop.mix = mix.probs
      prop.mix[(i-1),] = rdirichlet(1, alpha = .PiTuningParam * mix.probs[(i-1),] + .5)

      if(min(prop.mix[(i-1),]) < 0){
        r = -Inf
      }else{
        prop.post = likelihood.vec(M[i,], samp.size, G, sample = i-1, counts.mat, prop.mix, comp.probs, par = par) + log(ddirichlet(prop.mix[(i-1),], alpha = rep(1/G,G)))
        past.post = likelihood.vec(M[i,], samp.size, G, sample = i-1, counts.mat, mix.probs, comp.probs, par = par) + log(ddirichlet(mix.probs[(i-1),], alpha = rep(1/G,G)))
        ppropos = log(ddirichlet(mix.probs[(i-1),], alpha = .PiTuningParam * prop.mix[(i-1),] + .5)) - log(ddirichlet(prop.mix[(i-1),], alpha = .PiTuningParam * mix.probs[(i-1),] +.5))
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
      eta = mix.probs[(i-1),] + rnorm(G,0,sd = .PiTuningParam)
      prop.mix[(i-1),] = eta
      prop.mix[(i-1),] = prop.mix[(i-1),]/sum(prop.mix[(i-1),])

      if(min(prop.mix[(i-1),]) < 0){
        r = -Inf
      }else{
        prop.post = likelihood.vec(M[i,], samp.size, G, sample = i-1, counts.mat, mix.probs, comp.probs, cl=cl, par = par) + sum(dgamma(mix.probs[(i-1),], shape = 1/G, scale = 1, log = TRUE))
        past.post = likelihood.vec(M[i,], samp.size, G, sample = i-1, counts.mat, prop.mix, comp.probs, cl=cl, par = par) + sum(dgamma(prop.mix[(i-1),], shape = 1/G, scale = 1, log = TRUE))
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
        tmp.prop.all[,k] = MCMCpack::rdirichlet(1, alpha = .PsiTuningParam * c(comp.probs[1,(i-1),h,k], comp.probs[2,(i-1),h,k])+ 1)
      }

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
          bin.lik = likelihood.vec(M[i,],samp.size,G, sample = 2, counts.mat, mix.probs = mix.probs.hold, comp.probs =prop.probs, par = par, clust = cl)
        }
        bin.lik.new = likelihood.vec(M[i,],samp.size,G, sample = 2, counts.mat, mix.probs = mix.probs.hold, comp.probs =prop.tmp.probs, par = par, clust = cl)
        prop.post = bin.lik.new + sum(log(MCMCpack::ddirichlet(t(tmp.prop.all),alpha = rep(.5,2))))
        past.post = bin.lik + sum(log(MCMCpack::ddirichlet(cbind(comp.probs[1,i-1,h,], comp.probs[2,i-1,h,]),alpha = rep(.5,2))))
        ppropos = 0
        for(k in 1:P){
          ppropos =  ppropos + log(MCMCpack::ddirichlet(c(comp.probs[1,(i-1),h,k], comp.probs[2,(i-1),h,k]),alpha = .PsiTuningParam *tmp.prop.all[,k] + 1)) - log(MCMCpack::ddirichlet(tmp.prop.all[,k],alpha = .PsiTuningParam *c(comp.probs[1,(i-1),h,k], comp.probs[2,(i-1),h,k]) + 1))
        }
        r = exp(prop.post - past.post + ppropos)
      }
      if(stats::runif(1)<r){
        comp.probs[1,i,h,] = tmp.prop.all[1,]
        accept[h] = accept[h] + 1
        bin.lik = bin.lik.new
      }else{
        comp.probs[1,i,h,] = comp.probs[1,i-1,h,]
      }
      comp.probs[2,,,] = 1 - comp.probs[1,,,]
    }##End of for


    ###Calculate marginal probability for sample
    if(par){
      marg.prob[i,] = parallel::parApply(cl = cl, counts.mat,MARGIN = 1,FUN = twoWay.probs.vectorized, mix.probs= mix.probs, class.probs = comp.probs, sample = i, G=G)
    } else{
      marg.prob[i,] = apply(counts.mat,MARGIN = 1,FUN = twoWay.probs.vectorized, mix.probs= mix.probs, class.probs = comp.probs, sample = i, G=G)
    }

    ###Calculate full probability for sample
    if(.calculateFullTabProbs){
      for(z in 1:nrow(full.mat)){
        full.prob[i,z] = fullTab.probs(mix.probs, comp.probs, vals = full.mat[z,], sample = i)
      }
    }
    if(i %% 100 == 0) print(i)
  }##End of for
  return(list("marg_probs" = marg.prob, "full_probs" = full.prob, "M" = M, "Psi" = comp.probs, "Pi" = mix.probs, "accept_rate_pi" = accept.pi, "accept_rate_psi" = accept, "accept_rate_M" = accept.M, "epsilon" = eps))
}##End of function
