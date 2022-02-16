#' Helper function to calculate two-way probabilities
#' @keywords internal
getProbs = function(comb, cols, samp, G.used, class.probs){
  class.probs[cbind(comb,samp,G.used, cols)]
}


#' Helper function to calculate two-way probabilities
#'
#' @param combs.full A matrix of all combinations of the other variables that we need to sum over. Note: for two-way probabilities, two variables remain fixed, must iterate over all combinations of the others.
#' @param mix.probs Matrix of mixing probability
#' @param class.probs Array of class probability for latent classes
#' @param P The dimension of the table
#' @param sample What sample should the marginal probability be calculated for?
#' @param G The number of latent classes
#' @return A component of the marginal probability
#' @keywords internal
row.prob.vec <- function(combs.full, mix.probs, class.probs, P, sample, G){

  iter.vec = rep(1:P, nrow(combs.full))
  samp = rep(sample, nrow(combs.full)*ncol(combs.full))
  combs = c(t(combs.full))
  tmp.probs =  sapply(1:G, FUN = getProbs, comb = combs, col = iter.vec, samp = samp, class.probs = class.probs)
  tmp.probs = Rfast::rowprods(matrix(tmp.probs, ncol = P, byrow = TRUE))
  prob = sum(tmp.probs * rep(mix.probs[sample,],each = length(tmp.probs)/G))
  return(prob)
}

#' Vectorized function to count two-way marginal probabilities for a full vector or single count
#'
#' @param counts.vec A vector signifying chick combination is desired. Should only have two enteries and the rest should be `NA`
#' @param mix.probs Matrix of mixing probability
#' @param class.probs Array of class probability for latent classes
#' @param index1 Index of first marginal cell
#' @param index2 Index of second marginal cell. Must be different from index1.
#' @param val1 Value corresponding to the first marginal cell.
#' @param val2 Value corresponding to the second marginal cell.
#' @param sample What sample should the marginal probability be calculated for?
#' @param G The number of latent classes
#' @return A marginal probability for a cell
#' @keywords internal
#' @export
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
  combs.full = data.table::data.table(combs.full)
  prob = row.prob.vec(combs.full, mix.probs, class.probs, P, sample, G)

  return(prob)
}#end of function

#' Calculates full table probabilities
#'
#' @param mix.probs Matrix of mixing probability
#' @param class.probs Array of class probability for latent classes
#' @param vals The full table cell that should be calculated
#' @param sample What sample should the full cell probability be calculated for? `NULL` defaults to 1.
#' @return A full table probability for a given cell
#' @keywords internal
#' @export
fullTab.probs <- function(mix.probs, class.probs, vals, sample = NULL){
  ##Vals: a P dim vector
  if(length(vals) != dim(class.probs)[4]){
    return("Evaluated combination is of the wrong length!!")
  }
  ##Now, need to find the prob of each combination
  ##Easy way to do it: first find the prob of the combination for each latent class
  ##Then, multiply by the mixing probabilities
  
  if(is.null(sample)){
    sample = 1
  }
  
  prob.vec = rep(1, ncol(mix.probs))
  for(i in 1:length(vals)){
    prob.vec = prob.vec * class.probs[as.numeric(vals[i]), sample,,i]
  }
  prob = sum(mix.probs[sample,] * prob.vec)
  
  return(prob)
}#End of function
