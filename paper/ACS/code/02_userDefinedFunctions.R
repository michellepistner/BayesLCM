##This package specifies all necessary user-defined functions
##See comments on each function for more specific details
#########################################################################################

##Specifies full cell probabilities according to the latent class assignment
##Recall P_cell = sum_(k in 1:G) v_k prod(j=1:5) P_j
##where G = number of latent classes
##v_k is the mixing probability
##P_j represents the probability that a certain cell equals zero or one
########################################################################################
source(config.R)
Full.Probs <- function(G,pi, psi1, psi2, psi3, psi4, psi5){
  full.probs = rep(0,32)
  
  for(a in 1:G){ ##Looping over all the latent classes
    ##00000
    full.probs[1] = full.probs[1] + pi[a]*psi1[1,a]*psi2[1,a]*psi3[1,a]*psi4[1,a]*psi5[1,a]
    ##10000
    full.probs[2] = full.probs[2] + pi[a]*psi1[2,a]*psi2[1,a]*psi3[1,a]*psi4[1,a]*psi5[1,a]
    ##01000
    full.probs[3] = full.probs[3] + pi[a]*psi1[1,a]*psi2[2,a]*psi3[1,a]*psi4[1,a]*psi5[1,a]
    ##11000
    full.probs[4] = full.probs[4] + pi[a]*psi1[2,a]*psi2[2,a]*psi3[1,a]*psi4[1,a]*psi5[1,a]
    ##00100
    full.probs[5] = full.probs[5] + pi[a]*psi1[1,a]*psi2[1,a]*psi3[2,a]*psi4[1,a]*psi5[1,a]
    ##10100
    full.probs[6] = full.probs[6] + pi[a]*psi1[2,a]*psi2[1,a]*psi3[2,a]*psi4[1,a]*psi5[1,a]
    ##01100
    full.probs[7] = full.probs[7] + pi[a]*psi1[1,a]*psi2[2,a]*psi3[2,a]*psi4[1,a]*psi5[1,a]
    ##11100
    full.probs[8] = full.probs[8] + pi[a]*psi1[2,a]*psi2[2,a]*psi3[2,a]*psi4[1,a]*psi5[1,a]
    ##00010
    full.probs[9] = full.probs[9] + pi[a]*psi1[1,a]*psi2[1,a]*psi3[1,a]*psi4[2,a]*psi5[1,a]
    ##10010
    full.probs[10] = full.probs[10] + pi[a]*psi1[2,a]*psi2[1,a]*psi3[1,a]*psi4[2,a]*psi5[1,a]
    ##01010
    full.probs[11] = full.probs[11] + pi[a]*psi1[1,a]*psi2[2,a]*psi3[1,a]*psi4[2,a]*psi5[1,a]
    ##11010
    full.probs[12] = full.probs[12] + pi[a]*psi1[2,a]*psi2[2,a]*psi3[1,a]*psi4[2,a]*psi5[1,a]
    ##00110
    full.probs[13] = full.probs[13] + pi[a]*psi1[1,a]*psi2[1,a]*psi3[2,a]*psi4[2,a]*psi5[1,a]
    ##10110
    full.probs[14] = full.probs[14] + pi[a]*psi1[2,a]*psi2[1,a]*psi3[2,a]*psi4[2,a]*psi5[1,a]
    ##01110
    full.probs[15] = full.probs[15] + pi[a]*psi1[1,a]*psi2[2,a]*psi3[2,a]*psi4[2,a]*psi5[1,a]
    ##11110
    full.probs[16] = full.probs[16] + pi[a]*psi1[2,a]*psi2[2,a]*psi3[2,a]*psi4[2,a]*psi5[1,a]
    ##00001
    full.probs[17] = full.probs[17] + pi[a]*psi1[1,a]*psi2[1,a]*psi3[1,a]*psi4[1,a]*psi5[2,a]
    ##10001
    full.probs[18] = full.probs[18] + pi[a]*psi1[2,a]*psi2[1,a]*psi3[1,a]*psi4[1,a]*psi5[2,a]
    ##01001
    full.probs[19] = full.probs[19] + pi[a]*psi1[1,a]*psi2[2,a]*psi3[1,a]*psi4[1,a]*psi5[2,a]
    ##11001
    full.probs[20] = full.probs[20] + pi[a]*psi1[2,a]*psi2[2,a]*psi3[1,a]*psi4[1,a]*psi5[2,a]
    ##00101
    full.probs[21] = full.probs[21] + pi[a]*psi1[1,a]*psi2[1,a]*psi3[2,a]*psi4[1,a]*psi5[2,a]
    ##10101
    full.probs[22] = full.probs[22] + pi[a]*psi1[2,a]*psi2[1,a]*psi3[2,a]*psi4[1,a]*psi5[2,a]
    ##01101
    full.probs[23] = full.probs[23] + pi[a]*psi1[1,a]*psi2[2,a]*psi3[2,a]*psi4[1,a]*psi5[2,a]
    ##11101
    full.probs[24] = full.probs[24] + pi[a]*psi1[2,a]*psi2[2,a]*psi3[2,a]*psi4[1,a]*psi5[2,a]
    ##00011
    full.probs[25] = full.probs[25] + pi[a]*psi1[1,a]*psi2[1,a]*psi3[1,a]*psi4[2,a]*psi5[2,a]
    ##10011
    full.probs[26] = full.probs[26] + pi[a]*psi1[2,a]*psi2[1,a]*psi3[1,a]*psi4[2,a]*psi5[2,a]
    ##01011
    full.probs[27] = full.probs[27] + pi[a]*psi1[1,a]*psi2[2,a]*psi3[1,a]*psi4[2,a]*psi5[2,a]
    ##11011
    full.probs[28] = full.probs[28] + pi[a]*psi1[2,a]*psi2[2,a]*psi3[1,a]*psi4[2,a]*psi5[2,a]
    ##00111
    full.probs[29] = full.probs[29] + pi[a]*psi1[1,a]*psi2[1,a]*psi3[2,a]*psi4[2,a]*psi5[2,a]
    ##10111
    full.probs[30] = full.probs[30] + pi[a]*psi1[2,a]*psi2[1,a]*psi3[2,a]*psi4[2,a]*psi5[2,a]
    ##01111
    full.probs[31] = full.probs[31] + pi[a]*psi1[1,a]*psi2[2,a]*psi3[2,a]*psi4[2,a]*psi5[2,a]
    ##11111
    full.probs[32] = full.probs[32] + pi[a]*psi1[2,a]*psi2[2,a]*psi3[2,a]*psi4[2,a]*psi5[2,a]
  }
  return(full.probs)
}


##Similar to the "Full.Probs# function, but for all two-way marginal probabilities instead
##Sums over the correct full table probabilities to obtain the correct marginal probability

Marginal.Probs <- function(G,pi, psi1, psi2, psi3, psi4, psi5){
  marg.probs = rep(0,40)
  
  for(a in 1:G){ ##Looping over all the latent classes
    ##00...
    marg.probs[1] = marg.probs[1] + pi[a]*psi1[1,a]*psi2[1,a]*(psi3[1,a]*psi4[1,a]*psi5[1,a] + psi3[1,a]*psi4[1,a]*psi5[2,a] + psi3[1,a]*psi4[2,a]*psi5[1,a] + psi3[2,a]*psi4[1,a]*psi5[1,a] + psi3[1,a]*psi4[2,a]*psi5[2,a] + psi3[2,a]*psi4[1,a]*psi5[2,a] + psi3[2,a]*psi4[2,a]*psi5[1,a] + psi3[2,a]*psi4[2,a]*psi5[2,a])
    ##01...
    marg.probs[2] = marg.probs[2] + pi[a]*psi1[1,a]*psi2[2,a]*(psi3[1,a]*psi4[1,a]*psi5[1,a] + psi3[1,a]*psi4[1,a]*psi5[2,a] + psi3[1,a]*psi4[2,a]*psi5[1,a] + psi3[2,a]*psi4[1,a]*psi5[1,a] + psi3[1,a]*psi4[2,a]*psi5[2,a] + psi3[2,a]*psi4[1,a]*psi5[2,a] + psi3[2,a]*psi4[2,a]*psi5[1,a] + psi3[2,a]*psi4[2,a]*psi5[2,a])
    ##10...
    marg.probs[3] = marg.probs[3] + pi[a]*psi1[2,a]*psi2[1,a]*(psi3[1,a]*psi4[1,a]*psi5[1,a] + psi3[1,a]*psi4[1,a]*psi5[2,a] + psi3[1,a]*psi4[2,a]*psi5[1,a] + psi3[2,a]*psi4[1,a]*psi5[1,a] + psi3[1,a]*psi4[2,a]*psi5[2,a] + psi3[2,a]*psi4[1,a]*psi5[2,a] + psi3[2,a]*psi4[2,a]*psi5[1,a] + psi3[2,a]*psi4[2,a]*psi5[2,a])
    ##11...
    marg.probs[4] = marg.probs[4] + pi[a]*psi1[2,a]*psi2[2,a]*(psi3[1,a]*psi4[1,a]*psi5[1,a] + psi3[1,a]*psi4[1,a]*psi5[2,a] + psi3[1,a]*psi4[2,a]*psi5[1,a] + psi3[2,a]*psi4[1,a]*psi5[1,a] + psi3[1,a]*psi4[2,a]*psi5[2,a] + psi3[2,a]*psi4[1,a]*psi5[2,a] + psi3[2,a]*psi4[2,a]*psi5[1,a] + psi3[2,a]*psi4[2,a]*psi5[2,a])
    
    ##0.0..
    marg.probs[5] = marg.probs[5] + pi[a]*psi1[1,a]*psi3[1,a]*(psi2[1,a]*psi4[1,a]*psi5[1,a] + psi2[1,a]*psi4[1,a]*psi5[2,a] + psi2[1,a]*psi4[2,a]*psi5[1,a] + psi2[2,a]*psi4[1,a]*psi5[1,a] + psi2[1,a]*psi4[2,a]*psi5[2,a] + psi2[2,a]*psi4[1,a]*psi5[2,a] + psi2[2,a]*psi4[2,a]*psi5[1,a] + psi2[2,a]*psi4[2,a]*psi5[2,a])
    ##0.1..
    marg.probs[6] = marg.probs[6] + pi[a]*psi1[1,a]*psi3[2,a]*(psi2[1,a]*psi4[1,a]*psi5[1,a] + psi2[1,a]*psi4[1,a]*psi5[2,a] + psi2[1,a]*psi4[2,a]*psi5[1,a] + psi2[2,a]*psi4[1,a]*psi5[1,a] + psi2[1,a]*psi4[2,a]*psi5[2,a] + psi2[2,a]*psi4[1,a]*psi5[2,a] + psi2[2,a]*psi4[2,a]*psi5[1,a] + psi2[2,a]*psi4[2,a]*psi5[2,a])
    ##1.0..
    marg.probs[7] = marg.probs[7] + pi[a]*psi1[2,a]*psi3[1,a]*(psi2[1,a]*psi4[1,a]*psi5[1,a] + psi2[1,a]*psi4[1,a]*psi5[2,a] + psi2[1,a]*psi4[2,a]*psi5[1,a] + psi2[2,a]*psi4[1,a]*psi5[1,a] + psi2[1,a]*psi4[2,a]*psi5[2,a] + psi2[2,a]*psi4[1,a]*psi5[2,a] + psi2[2,a]*psi4[2,a]*psi5[1,a] + psi2[2,a]*psi4[2,a]*psi5[2,a])
    ##1.1..
    marg.probs[8] = marg.probs[8] + pi[a]*psi1[2,a]*psi3[2,a]*(psi2[1,a]*psi4[1,a]*psi5[1,a] + psi2[1,a]*psi4[1,a]*psi5[2,a] + psi2[1,a]*psi4[2,a]*psi5[1,a] + psi2[2,a]*psi4[1,a]*psi5[1,a] + psi2[1,a]*psi4[2,a]*psi5[2,a] + psi2[2,a]*psi4[1,a]*psi5[2,a] + psi2[2,a]*psi4[2,a]*psi5[1,a] + psi2[2,a]*psi4[2,a]*psi5[2,a])
    
    ##0..0.
    marg.probs[9] = marg.probs[9] + pi[a]*psi1[1,a]*psi4[1,a]*(psi2[1,a]*psi3[1,a]*psi5[1,a] + psi2[1,a]*psi3[1,a]*psi5[2,a] + psi2[1,a]*psi3[2,a]*psi5[1,a] + psi2[2,a]*psi3[1,a]*psi5[1,a] + psi2[1,a]*psi3[2,a]*psi5[2,a] + psi2[2,a]*psi3[1,a]*psi5[2,a] + psi2[2,a]*psi3[2,a]*psi5[1,a] + psi2[2,a]*psi3[2,a]*psi5[2,a])
    ##0..1.
    marg.probs[10] = marg.probs[10] + pi[a]*psi1[1,a]*psi4[2,a]*(psi2[1,a]*psi3[1,a]*psi5[1,a] + psi2[1,a]*psi3[1,a]*psi5[2,a] + psi2[1,a]*psi3[2,a]*psi5[1,a] + psi2[2,a]*psi3[1,a]*psi5[1,a] + psi2[1,a]*psi3[2,a]*psi5[2,a] + psi2[2,a]*psi3[1,a]*psi5[2,a] + psi2[2,a]*psi3[2,a]*psi5[1,a] + psi2[2,a]*psi3[2,a]*psi5[2,a])
    ##1..0.
    marg.probs[11] = marg.probs[11] + pi[a]*psi1[2,a]*psi4[1,a]*(psi2[1,a]*psi3[1,a]*psi5[1,a] + psi2[1,a]*psi3[1,a]*psi5[2,a] + psi2[1,a]*psi3[2,a]*psi5[1,a] + psi2[2,a]*psi3[1,a]*psi5[1,a] + psi2[1,a]*psi3[2,a]*psi5[2,a] + psi2[2,a]*psi3[1,a]*psi5[2,a] + psi2[2,a]*psi3[2,a]*psi5[1,a] + psi2[2,a]*psi3[2,a]*psi5[2,a])
    ##1..1.
    marg.probs[12] = marg.probs[12] + pi[a]*psi1[2,a]*psi4[2,a]*(psi2[1,a]*psi3[1,a]*psi5[1,a] + psi2[1,a]*psi3[1,a]*psi5[2,a] + psi2[1,a]*psi3[2,a]*psi5[1,a] + psi2[2,a]*psi3[1,a]*psi5[1,a] + psi2[1,a]*psi3[2,a]*psi5[2,a] + psi2[2,a]*psi3[1,a]*psi5[2,a] + psi2[2,a]*psi3[2,a]*psi5[1,a] + psi2[2,a]*psi3[2,a]*psi5[2,a])
    
    ##0...0
    marg.probs[13] = marg.probs[13] + pi[a]*psi1[1,a]*psi5[1,a]*(psi2[1,a]*psi3[1,a]*psi4[1,a] + psi2[1,a]*psi3[1,a]*psi4[2,a] + psi2[1,a]*psi3[2,a]*psi4[1,a] + psi2[2,a]*psi3[1,a]*psi4[1,a] + psi2[1,a]*psi3[2,a]*psi4[2,a] + psi2[2,a]*psi3[1,a]*psi4[2,a] + psi2[2,a]*psi3[2,a]*psi4[1,a] + psi2[2,a]*psi3[2,a]*psi4[2,a])
    ##0...1
    marg.probs[14] = marg.probs[14] + pi[a]*psi1[1,a]*psi5[2,a]*(psi2[1,a]*psi3[1,a]*psi4[1,a] + psi2[1,a]*psi3[1,a]*psi4[2,a] + psi2[1,a]*psi3[2,a]*psi4[1,a] + psi2[2,a]*psi3[1,a]*psi4[1,a] + psi2[1,a]*psi3[2,a]*psi4[2,a] + psi2[2,a]*psi3[1,a]*psi4[2,a] + psi2[2,a]*psi3[2,a]*psi4[1,a] + psi2[2,a]*psi3[2,a]*psi4[2,a])
    ##1...0
    marg.probs[15] = marg.probs[15] + pi[a]*psi1[2,a]*psi5[1,a]*(psi2[1,a]*psi3[1,a]*psi4[1,a] + psi2[1,a]*psi3[1,a]*psi4[2,a] + psi2[1,a]*psi3[2,a]*psi4[1,a] + psi2[2,a]*psi3[1,a]*psi4[1,a] + psi2[1,a]*psi3[2,a]*psi4[2,a] + psi2[2,a]*psi3[1,a]*psi4[2,a] + psi2[2,a]*psi3[2,a]*psi4[1,a] + psi2[2,a]*psi3[2,a]*psi4[2,a])
    ##1...1
    marg.probs[16] = marg.probs[16] + pi[a]*psi1[2,a]*psi5[2,a]*(psi2[1,a]*psi3[1,a]*psi4[1,a] + psi2[1,a]*psi3[1,a]*psi4[2,a] + psi2[1,a]*psi3[2,a]*psi4[1,a] + psi2[2,a]*psi3[1,a]*psi4[1,a] + psi2[1,a]*psi3[2,a]*psi4[2,a] + psi2[2,a]*psi3[1,a]*psi4[2,a] + psi2[2,a]*psi3[2,a]*psi4[1,a] + psi2[2,a]*psi3[2,a]*psi4[2,a])
    
    ##.00..
    marg.probs[17] = marg.probs[17] + pi[a]*psi2[1,a]*psi3[1,a]*(psi1[1,a]*psi4[1,a]*psi5[1,a] + psi1[1,a]*psi4[1,a]*psi5[2,a] + psi1[1,a]*psi4[2,a]*psi5[1,a] + psi1[2,a]*psi4[1,a]*psi5[1,a] + psi1[1,a]*psi4[2,a]*psi5[2,a] + psi1[2,a]*psi4[1,a]*psi5[2,a] + psi1[2,a]*psi4[2,a]*psi5[1,a] + psi1[2,a]*psi4[2,a]*psi5[2,a])
    ##.01..
    marg.probs[18] = marg.probs[18] + pi[a]*psi2[1,a]*psi3[2,a]*(psi1[1,a]*psi4[1,a]*psi5[1,a] + psi1[1,a]*psi4[1,a]*psi5[2,a] + psi1[1,a]*psi4[2,a]*psi5[1,a] + psi1[2,a]*psi4[1,a]*psi5[1,a] + psi1[1,a]*psi4[2,a]*psi5[2,a] + psi1[2,a]*psi4[1,a]*psi5[2,a] + psi1[2,a]*psi4[2,a]*psi5[1,a] + psi1[2,a]*psi4[2,a]*psi5[2,a])
    ##.10..
    marg.probs[19] = marg.probs[19] + pi[a]*psi2[2,a]*psi3[1,a]*(psi1[1,a]*psi4[1,a]*psi5[1,a] + psi1[1,a]*psi4[1,a]*psi5[2,a] + psi1[1,a]*psi4[2,a]*psi5[1,a] + psi1[2,a]*psi4[1,a]*psi5[1,a] + psi1[1,a]*psi4[2,a]*psi5[2,a] + psi1[2,a]*psi4[1,a]*psi5[2,a] + psi1[2,a]*psi4[2,a]*psi5[1,a] + psi1[2,a]*psi4[2,a]*psi5[2,a])
    ##.11..
    marg.probs[20] = marg.probs[20] + pi[a]*psi2[2,a]*psi3[2,a]*(psi1[1,a]*psi4[1,a]*psi5[1,a] + psi1[1,a]*psi4[1,a]*psi5[2,a] + psi1[1,a]*psi4[2,a]*psi5[1,a] + psi1[2,a]*psi4[1,a]*psi5[1,a] + psi1[1,a]*psi4[2,a]*psi5[2,a] + psi1[2,a]*psi4[1,a]*psi5[2,a] + psi1[2,a]*psi4[2,a]*psi5[1,a] + psi1[2,a]*psi4[2,a]*psi5[2,a])
    
    ##.0.0.
    marg.probs[21] = marg.probs[21] + pi[a]*psi2[1,a]*psi4[1,a]*(psi1[1,a]*psi3[1,a]*psi5[1,a] + psi1[1,a]*psi3[1,a]*psi5[2,a] + psi1[1,a]*psi3[2,a]*psi5[1,a] + psi1[2,a]*psi3[1,a]*psi5[1,a] + psi1[1,a]*psi3[2,a]*psi5[2,a] + psi1[2,a]*psi3[1,a]*psi5[2,a] + psi1[2,a]*psi3[2,a]*psi5[1,a] + psi1[2,a]*psi3[2,a]*psi5[2,a])
    ##.0.1.
    marg.probs[22] = marg.probs[22] + pi[a]*psi2[1,a]*psi4[2,a]*(psi1[1,a]*psi3[1,a]*psi5[1,a] + psi1[1,a]*psi3[1,a]*psi5[2,a] + psi1[1,a]*psi3[2,a]*psi5[1,a] + psi1[2,a]*psi3[1,a]*psi5[1,a] + psi1[1,a]*psi3[2,a]*psi5[2,a] + psi1[2,a]*psi3[1,a]*psi5[2,a] + psi1[2,a]*psi3[2,a]*psi5[1,a] + psi1[2,a]*psi3[2,a]*psi5[2,a])
    ##.1.0.
    marg.probs[23] = marg.probs[23] + pi[a]*psi2[2,a]*psi4[1,a]*(psi1[1,a]*psi3[1,a]*psi5[1,a] + psi1[1,a]*psi3[1,a]*psi5[2,a] + psi1[1,a]*psi3[2,a]*psi5[1,a] + psi1[2,a]*psi3[1,a]*psi5[1,a] + psi1[1,a]*psi3[2,a]*psi5[2,a] + psi1[2,a]*psi3[1,a]*psi5[2,a] + psi1[2,a]*psi3[2,a]*psi5[1,a] + psi1[2,a]*psi3[2,a]*psi5[2,a])
    ##.1.1.
    marg.probs[24] = marg.probs[24] + pi[a]*psi2[2,a]*psi4[2,a]*(psi1[1,a]*psi3[1,a]*psi5[1,a] + psi1[1,a]*psi3[1,a]*psi5[2,a] + psi1[1,a]*psi3[2,a]*psi5[1,a] + psi1[2,a]*psi3[1,a]*psi5[1,a] + psi1[1,a]*psi3[2,a]*psi5[2,a] + psi1[2,a]*psi3[1,a]*psi5[2,a] + psi1[2,a]*psi3[2,a]*psi5[1,a] + psi1[2,a]*psi3[2,a]*psi5[2,a])
    
    ##.0..0
    marg.probs[25] = marg.probs[25] + pi[a]*psi2[1,a]*psi5[1,a]*(psi1[1,a]*psi3[1,a]*psi4[1,a] + psi1[1,a]*psi3[1,a]*psi4[2,a] + psi1[1,a]*psi3[2,a]*psi4[1,a] + psi1[2,a]*psi3[1,a]*psi4[1,a] + psi1[1,a]*psi3[2,a]*psi4[2,a] + psi1[2,a]*psi3[1,a]*psi4[2,a] + psi1[2,a]*psi3[2,a]*psi4[1,a] + psi1[2,a]*psi3[2,a]*psi4[2,a])
    ##.0..1
    marg.probs[26] = marg.probs[26] + pi[a]*psi2[1,a]*psi5[2,a]*(psi1[1,a]*psi3[1,a]*psi4[1,a] + psi1[1,a]*psi3[1,a]*psi4[2,a] + psi1[1,a]*psi3[2,a]*psi4[1,a] + psi1[2,a]*psi3[1,a]*psi4[1,a] + psi1[1,a]*psi3[2,a]*psi4[2,a] + psi1[2,a]*psi3[1,a]*psi4[2,a] + psi1[2,a]*psi3[2,a]*psi4[1,a] + psi1[2,a]*psi3[2,a]*psi4[2,a])
    ##.1..0
    marg.probs[27] = marg.probs[27] + pi[a]*psi2[2,a]*psi5[1,a]*(psi1[1,a]*psi3[1,a]*psi4[1,a] + psi1[1,a]*psi3[1,a]*psi4[2,a] + psi1[1,a]*psi3[2,a]*psi4[1,a] + psi1[2,a]*psi3[1,a]*psi4[1,a] + psi1[1,a]*psi3[2,a]*psi4[2,a] + psi1[2,a]*psi3[1,a]*psi4[2,a] + psi1[2,a]*psi3[2,a]*psi4[1,a] + psi1[2,a]*psi3[2,a]*psi4[2,a])
    ##.1..1
    marg.probs[28] = marg.probs[28] + pi[a]*psi2[2,a]*psi5[2,a]*(psi1[1,a]*psi3[1,a]*psi4[1,a] + psi1[1,a]*psi3[1,a]*psi4[2,a] + psi1[1,a]*psi3[2,a]*psi4[1,a] + psi1[2,a]*psi3[1,a]*psi4[1,a] + psi1[1,a]*psi3[2,a]*psi4[2,a] + psi1[2,a]*psi3[1,a]*psi4[2,a] + psi1[2,a]*psi3[2,a]*psi4[1,a] + psi1[2,a]*psi3[2,a]*psi4[2,a])
    
    ##..00.
    marg.probs[29] = marg.probs[29] + pi[a]*psi3[1,a]*psi4[1,a]*(psi1[1,a]*psi3[1,a]*psi5[1,a] + psi1[1,a]*psi3[1,a]*psi5[2,a] + psi1[1,a]*psi3[2,a]*psi5[1,a] + psi1[2,a]*psi3[1,a]*psi5[1,a] + psi1[1,a]*psi3[2,a]*psi5[2,a] + psi1[2,a]*psi3[1,a]*psi5[2,a] + psi1[2,a]*psi3[2,a]*psi5[1,a] + psi1[2,a]*psi3[2,a]*psi5[2,a])
    ##..01.
    marg.probs[30] = marg.probs[30] + pi[a]*psi3[1,a]*psi4[2,a]*(psi1[1,a]*psi3[1,a]*psi5[1,a] + psi1[1,a]*psi3[1,a]*psi5[2,a] + psi1[1,a]*psi3[2,a]*psi5[1,a] + psi1[2,a]*psi3[1,a]*psi5[1,a] + psi1[1,a]*psi3[2,a]*psi5[2,a] + psi1[2,a]*psi3[1,a]*psi5[2,a] + psi1[2,a]*psi3[2,a]*psi5[1,a] + psi1[2,a]*psi3[2,a]*psi5[2,a])
    ##..10.
    marg.probs[31] = marg.probs[31] + pi[a]*psi3[2,a]*psi4[1,a]*(psi1[1,a]*psi3[1,a]*psi5[1,a] + psi1[1,a]*psi3[1,a]*psi5[2,a] + psi1[1,a]*psi3[2,a]*psi5[1,a] + psi1[2,a]*psi3[1,a]*psi5[1,a] + psi1[1,a]*psi3[2,a]*psi5[2,a] + psi1[2,a]*psi3[1,a]*psi5[2,a] + psi1[2,a]*psi3[2,a]*psi5[1,a] + psi1[2,a]*psi3[2,a]*psi5[2,a])
    ##..11.
    marg.probs[32] = marg.probs[32] + pi[a]*psi3[2,a]*psi4[2,a]*(psi1[1,a]*psi3[1,a]*psi5[1,a] + psi1[1,a]*psi3[1,a]*psi5[2,a] + psi1[1,a]*psi3[2,a]*psi5[1,a] + psi1[2,a]*psi3[1,a]*psi5[1,a] + psi1[1,a]*psi3[2,a]*psi5[2,a] + psi1[2,a]*psi3[1,a]*psi5[2,a] + psi1[2,a]*psi3[2,a]*psi5[1,a] + psi1[2,a]*psi3[2,a]*psi5[2,a])
    
    ##..0.0
    marg.probs[33] = marg.probs[33] + pi[a]*psi3[1,a]*psi5[1,a]*(psi1[1,a]*psi3[1,a]*psi4[1,a] + psi1[1,a]*psi3[1,a]*psi4[2,a] + psi1[1,a]*psi3[2,a]*psi4[1,a] + psi1[2,a]*psi3[1,a]*psi4[1,a] + psi1[1,a]*psi3[2,a]*psi4[2,a] + psi1[2,a]*psi3[1,a]*psi4[2,a] + psi1[2,a]*psi3[2,a]*psi4[1,a] + psi1[2,a]*psi3[2,a]*psi4[2,a])
    ##..0.1
    marg.probs[34] = marg.probs[34] + pi[a]*psi3[1,a]*psi5[2,a]*(psi1[1,a]*psi3[1,a]*psi4[1,a] + psi1[1,a]*psi3[1,a]*psi4[2,a] + psi1[1,a]*psi3[2,a]*psi4[1,a] + psi1[2,a]*psi3[1,a]*psi4[1,a] + psi1[1,a]*psi3[2,a]*psi4[2,a] + psi1[2,a]*psi3[1,a]*psi4[2,a] + psi1[2,a]*psi3[2,a]*psi4[1,a] + psi1[2,a]*psi3[2,a]*psi4[2,a])
    ##..1.0
    marg.probs[35] = marg.probs[35] + pi[a]*psi3[2,a]*psi5[1,a]*(psi1[1,a]*psi3[1,a]*psi4[1,a] + psi1[1,a]*psi3[1,a]*psi4[2,a] + psi1[1,a]*psi3[2,a]*psi4[1,a] + psi1[2,a]*psi3[1,a]*psi4[1,a] + psi1[1,a]*psi3[2,a]*psi4[2,a] + psi1[2,a]*psi3[1,a]*psi4[2,a] + psi1[2,a]*psi3[2,a]*psi4[1,a] + psi1[2,a]*psi3[2,a]*psi4[2,a])
    ##..1.1
    marg.probs[36] = marg.probs[36] + pi[a]*psi3[2,a]*psi5[2,a]*(psi1[1,a]*psi3[1,a]*psi4[1,a] + psi1[1,a]*psi3[1,a]*psi4[2,a] + psi1[1,a]*psi3[2,a]*psi4[1,a] + psi1[2,a]*psi3[1,a]*psi4[1,a] + psi1[1,a]*psi3[2,a]*psi4[2,a] + psi1[2,a]*psi3[1,a]*psi4[2,a] + psi1[2,a]*psi3[2,a]*psi4[1,a] + psi1[2,a]*psi3[2,a]*psi4[2,a])
    
    ##...00
    marg.probs[37] = marg.probs[37] + pi[a]*psi4[1,a]*psi5[1,a]*(psi1[1,a]*psi2[1,a]*psi3[1,a] + psi1[1,a]*psi2[1,a]*psi3[2,a] + psi1[1,a]*psi2[2,a]*psi3[1,a] + psi1[2,a]*psi2[1,a]*psi3[1,a] + psi1[1,a]*psi2[2,a]*psi3[2,a] + psi1[2,a]*psi2[1,a]*psi3[2,a] + psi1[2,a]*psi2[2,a]*psi3[1,a] + psi1[2,a]*psi2[2,a]*psi3[2,a])
    ##...01
    marg.probs[38] = marg.probs[38] + pi[a]*psi4[1,a]*psi5[2,a]*(psi1[1,a]*psi2[1,a]*psi3[1,a] + psi1[1,a]*psi2[1,a]*psi3[2,a] + psi1[1,a]*psi2[2,a]*psi3[1,a] + psi1[2,a]*psi2[1,a]*psi3[1,a] + psi1[1,a]*psi2[2,a]*psi3[2,a] + psi1[2,a]*psi2[1,a]*psi3[2,a] + psi1[2,a]*psi2[2,a]*psi3[1,a] + psi1[2,a]*psi2[2,a]*psi3[2,a])
    ##...10
    marg.probs[39] = marg.probs[39] + pi[a]*psi4[2,a]*psi5[1,a]*(psi1[1,a]*psi2[1,a]*psi3[1,a] + psi1[1,a]*psi2[1,a]*psi3[2,a] + psi1[1,a]*psi2[2,a]*psi3[1,a] + psi1[2,a]*psi2[1,a]*psi3[1,a] + psi1[1,a]*psi2[2,a]*psi3[2,a] + psi1[2,a]*psi2[1,a]*psi3[2,a] + psi1[2,a]*psi2[2,a]*psi3[1,a] + psi1[2,a]*psi2[2,a]*psi3[2,a])
    ##...11
    marg.probs[40] = marg.probs[40] + pi[a]*psi4[2,a]*psi5[2,a]*(psi1[1,a]*psi2[1,a]*psi3[1,a] + psi1[1,a]*psi2[1,a]*psi3[2,a] + psi1[1,a]*psi2[2,a]*psi3[1,a] + psi1[2,a]*psi2[1,a]*psi3[1,a] + psi1[1,a]*psi2[2,a]*psi3[2,a] + psi1[2,a]*psi2[1,a]*psi3[2,a] + psi1[2,a]*psi2[2,a]*psi3[1,a] + psi1[2,a]*psi2[2,a]*psi3[2,a])
  }
  return(marg.probs)
}

##Likelihood function for marginal counts
##Recall M_i ~ Bin(N,P) indep.
likelihood <- function(counts,samp.size, G, pi, psi1, psi2, psi3, psi4, psi5){
  probs = Marginal.Probs(G, pi, psi1, psi2, psi3, psi4, psi5)
  lik = sum(dbinom(counts,samp.size, probs,log=T))
  return(lik)
}#end of function

ind.likelihood <- function(counts,samp.size, G, pi, psi1, psi2, psi3, psi4, psi5,b){
  probs = Marginal.Probs(G, pi, psi1, psi2, psi3, psi4, psi5)[b]
  lik = dbinom(counts,samp.size, probs,log=T)
  return(lik)
}#end of function


##Likelihood defined for the measurement error part of the model
##Recall that the noisy count is drawn from a discrete Laplace distribution
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
  if(noisyS == 0){
    pmf = alpha^S/(1+alpha)
  }else if(noisyS == samp.size){
    pmf = alpha^(samp.size-S)/(1+alpha)
  }else{
    pmf = (1-alpha)*alpha^S/(1+alpha)
  }
  return(log(pmf))
}


###This function defines the Metropolis-Hastings algorithm for the model
##"counts"= noisy counts as input, "G" = # of latent classes
## "samp.size" = overall data size, "N.runs" = number of runs
## "E" = alpha as defined by the double Laplace distribution (function of epsilon)
## "prop.weights" = proposal weights, if something different from equal is desired
TwoWay_NoisyLCM <- function(counts, G, samp.size, N.runs,E,prop.weights){
  ##Initializing vectors to hold the class membership probabilities
  ##Each is a matrix-- "c*" represents variable * where .1 or .2 represent the levels
  ##Each column represents a latent class
  psi.c1.1 = matrix(NA,nrow = N.runs, ncol = G)
  psi.c1.2 = matrix(NA,nrow = N.runs, ncol = G)
  
  psi.c2.1 = matrix(NA,nrow = N.runs, ncol = G)
  psi.c2.2 = matrix(NA,nrow = N.runs, ncol = G)
  
  psi.c3.1 = matrix(NA,nrow = N.runs, ncol = G)
  psi.c3.2 = matrix(NA,nrow = N.runs, ncol = G)
  
  psi.c4.1 = matrix(NA,nrow = N.runs, ncol = G)
  psi.c4.2 = matrix(NA,nrow = N.runs, ncol = G)
  
  psi.c5.1 = matrix(NA,nrow = N.runs, ncol = G)
  psi.c5.2 = matrix(NA,nrow = N.runs, ncol = G)
  
  pi.runs = matrix(NA, nrow = N.runs, ncol = G)
  
  M = matrix(NA, ncol=N.runs, nrow=40)
  
  P.marg = matrix(NA, nrow = N.runs, ncol = 40)
  
  proposal.weights.1 = c(sum(counts[1:2]),sum(counts[3:4]))
  proposal.weights.2 = c(sum(counts[c(1,3)]),sum(counts[c(2,4)]))
  proposal.weights.3 = c(sum(counts[c(5,7)]),sum(counts[c(6,8)]))
  proposal.weights.4 = c(sum(counts[c(9,11)]),sum(counts[c(10,12)]))
  proposal.weights.5 = c(sum(counts[c(13,15)]),sum(counts[c(14,16)]))
  
  ##Startings values
  ##Just assuming equal membership across latent classes
  pi.runs[1,] = rep(1/G,G)
  psi.c1.1[1,] = rep(.5,G)
  psi.c1.2[1,] = rep(.5,G)
  
  psi.c2.1[1,] = rep(.5,G)
  psi.c2.2[1,] = rep(.5,G)
  
  psi.c3.1[1,] = rep(.5,G)
  psi.c3.2[1,] = rep(.5,G)
  
  psi.c4.1[1,] = rep(.5,G)
  psi.c4.2[1,] = rep(.5,G)
  
  psi.c5.1[1,] = rep(.5,G)
  psi.c5.2[1,] = rep(.5,G)
  
  M[,1] = counts
  
  psi1 = rbind(psi.c1.1[1,],psi.c1.2[1,])
  psi2 = rbind(psi.c2.1[1,],psi.c2.2[1,])
  psi3 = rbind(psi.c3.1[1,],psi.c3.2[1,])
  psi4 = rbind(psi.c4.1[1,],psi.c4.2[1,])
  psi5 = rbind(psi.c5.1[1,],psi.c5.2[1,])
  pi = pi.runs[1,]
  
  accept1 = matrix(0,G)
  accept2 = matrix(0,G)
  accept3 = matrix(0,G)
  accept4 = matrix(0,G)
  accept5 = matrix(0,G)
  accept.pi = 0
  accept.counts = rep(0,40)
  
  ##Now, beging the M-H algorithm
  
  for(j in 2:N.runs){
    k = j-1
    
    psi1 = rbind(psi.c1.1[k,],psi.c1.2[k,])
    psi2 = rbind(psi.c2.1[k,],psi.c2.2[k,])
    psi3 = rbind(psi.c3.1[k,],psi.c3.2[k,])
    psi4 = rbind(psi.c4.1[k,],psi.c4.2[k,])
    psi5 = rbind(psi.c5.1[k,],psi.c5.2[k,])
    pi = pi.runs[k,]
    
    M.tmp = M[,k]
    for(b in 1:40){
      
      samp.E = -0.01
      M.prop = counts[b] + rgeom(1,1-exp(samp.E)) - rgeom(1,1-exp(samp.E))
      if(M.prop < 0){
        M.prop = 1}
      if(M.prop > samp.size){
        M.prop = samp.size  
        }

        prop.post = lap.lik(counts[b],M.prop,E=E) + ind.likelihood(M.prop,samp.size,G,pi,psi1,psi2,psi3,psi4,psi5,b)
        past.post = lap.lik(counts[b],M[b,k],E=E) + ind.likelihood(M[b,k],samp.size,G,pi,psi1,psi2,psi3,psi4,psi5,b)
        ppropos = dT2Lap(M[b,k],M.prop,samp.size,samp.E) - dT2Lap(M.prop,M[b,k],samp.size,samp.E)

                
        r = exp(prop.post - past.post +ppropos)
        
        
        if(runif(1)<r){
          M[b,j] = M.prop
          M.tmp[b] = M.prop
          accept.counts[b] = accept.counts[b] + 1
        }else{
          M[b,j] = M[b,k]
        }
    }#end of for
    
    ##Mixing proportions
    prop.mix = rnorm(G,pi.runs[k,]*100,sd=.5)##rdirichlet(1,alpha = pi.runs[k,]+1)
    prop.mix = prop.mix/sum(prop.mix)
    if(min(prop.mix) < 0){
      r = -Inf
    }else{
    psi1 = rbind(psi.c1.1[k,],psi.c1.2[k,])
    psi2 = rbind(psi.c2.1[k,],psi.c2.2[k,])
    psi3 = rbind(psi.c3.1[k,],psi.c3.2[k,])
    psi4 = rbind(psi.c4.1[k,],psi.c4.2[k,])
    psi5 = rbind(psi.c5.1[k,],psi.c5.2[k,])
    
    prop.post = likelihood(M[,k],samp.size,G,prop.mix,psi1,psi2,psi3,psi4,psi5) + log(ddirichlet(prop.mix,alpha = rep(1,G)))
    past.post = likelihood(M[,k],samp.size,G,pi.runs[k,],psi1,psi2,psi3,psi4,psi5) + log(ddirichlet(pi.runs[k,],alpha = rep(1,G)))
    ##ppropos = log(ddirichlet(pi.runs[k,],alpha = prop.mix)) - log(ddirichlet(prop.mix,alpha = pi.runs[k,]))
    ##ppropos = log(ddirichlet(pi.runs[k,], alpha=rep(1,G))) - log(ddirichlet(prop.mix, alpha=rep(1,G)))
    
    r = exp(prop.post - past.post)
    }
    if(runif(1)<r){
      pi.runs[j,] = prop.mix
      accept.pi = accept.pi+1
    }else{
      pi.runs[j,] = pi.runs[k,]
    }
    
    ##Now, for psi1
    pi = pi.runs[j,]
    psi1 = rbind(psi.c1.1[k,],psi.c1.2[k,])
    psi2 = rbind(psi.c2.1[k,],psi.c2.2[k,])
    psi3 = rbind(psi.c3.1[k,],psi.c3.2[k,])
    psi4 = rbind(psi.c4.1[k,],psi.c4.2[k,])
    psi5 = rbind(psi.c5.1[k,],psi.c5.2[k,])
    for(a in 1:G){
      #tmp.prop = rnorm(2,psi1[,a]*100,sd=5)
      #tmp.prop = tmp.prop/sum(tmp.prop)
      
      tmp.prop = rdirichlet(1,alpha = 2*psi5[,a]+.2)
      psi1.prop = psi1
      psi1.prop[,a] = tmp.prop

      prop.post = likelihood(M[,k],samp.size,G,pi,psi1.prop,psi2,psi3,psi4,psi5) + log(ddirichlet(tmp.prop,alpha = c(1,5)))
      past.post = likelihood(M[,k],samp.size,G,pi,psi1,psi2,psi3,psi4,psi5) + log(ddirichlet(psi1[,a],alpha = c(1,5)))
      #ppropos = log(ddirichlet(psi1[,a],alpha = tmp.prop)) - log(ddirichlet(tmp.prop,alpha = psi1[,a]))
      ppropos = log(ddirichlet(psi1[,a], 2*psi1[,a]+.2)) - log(ddirichlet(tmp.prop, 2*psi1[,a]+.2))
      
      r = exp(prop.post - past.post + ppropos)
      
      if(runif(1)<r){
        psi.c1.1[j,a] = tmp.prop[1]
        psi.c1.2[j,a] = tmp.prop[2]
        psi1[,a] = tmp.prop
        accept1[a] = accept1[a]+1
      }else{
        psi.c1.1[j,a] = psi.c1.1[k,a]
        psi.c1.2[j,a] = psi.c1.2[k,a]
      }
    }
    
    pi = pi.runs[j,]
    psi1 = rbind(psi.c1.1[j,],psi.c1.2[j,])
    psi2 = rbind(psi.c2.1[k,],psi.c2.2[k,])
    psi3 = rbind(psi.c3.1[k,],psi.c3.2[k,])
    psi4 = rbind(psi.c4.1[k,],psi.c4.2[k,])
    psi5 = rbind(psi.c5.1[k,],psi.c5.2[k,])
    for(a in 1:G){
      #tmp.prop = rnorm(2,psi2[,a]*100,sd=5)
      #tmp.prop = tmp.prop/sum(tmp.prop)
      tmp.prop = rdirichlet(1,alpha = 2*psi2[,a]+.2)
      psi2.prop = psi2
      psi2.prop[,a] = tmp.prop
      
      prop.post = likelihood(M[,k],samp.size,G,pi,psi1,psi2.prop,psi3,psi4,psi5) + log(ddirichlet(tmp.prop,alpha = c(1,5)))
      past.post = likelihood(M[,k],samp.size,G,pi,psi1,psi2,psi3,psi4,psi5) + log(ddirichlet(psi2[,a],alpha = c(1,5)))
      ##ppropos = log(ddirichlet(psi2[,a],alpha = tmp.prop)) - log(ddirichlet(tmp.prop,alpha = psi2[,a]))
      ppropos = log(ddirichlet(psi2[,a], 2*psi2[,a]+.2)) - log(ddirichlet(tmp.prop, 2*psi2[,a]+.2))
      
      r = exp(prop.post - past.post + ppropos)
      
      if(runif(1)<r){
        psi.c2.1[j,a] = tmp.prop[1]
        psi.c2.2[j,a] = tmp.prop[2]
        psi2[,a] = tmp.prop
        accept2[a] = accept2[a]+1
      }else{
        psi.c2.1[j,a] = psi.c2.1[k,a]
        psi.c2.2[j,a] = psi.c2.2[k,a]
      }
    }
    
    pi = pi.runs[j,]
    psi1 = rbind(psi.c1.1[j,],psi.c1.2[j,])
    psi2 = rbind(psi.c2.1[j,],psi.c2.2[j,])
    psi3 = rbind(psi.c3.1[k,],psi.c3.2[k,])
    psi4 = rbind(psi.c4.1[k,],psi.c4.2[k,])
    psi5 = rbind(psi.c5.1[k,],psi.c5.2[k,])
    for(a in 1:G){
      #tmp.prop = rnorm(2,psi3[,a]*100,sd=5)
      #tmp.prop = tmp.prop/sum(tmp.prop)
      tmp.prop = rdirichlet(1,alpha = 2*psi3[,a]+.2)
      psi3.prop = psi3
      psi3.prop[,a] = tmp.prop
      
      prop.post = likelihood(M[,k],samp.size,G,pi,psi1,psi2,psi3.prop,psi4,psi5) + log(ddirichlet(tmp.prop,alpha = c(1,5)))
      past.post = likelihood(M[,k],samp.size,G,pi,psi1,psi2,psi3,psi4,psi5) + log(ddirichlet(psi3[,a],alpha = c(1,5)))
      ##ppropos = log(ddirichlet(psi3[,a],alpha = tmp.prop)) - log(ddirichlet(tmp.prop,alpha = psi3[,a]))
      ppropos = log(ddirichlet(psi3[,a], 2*psi3[,a]+.2)) - log(ddirichlet(tmp.prop,2*psi3[,a]+.2))
      
      r = exp(prop.post - past.post + ppropos)
      
      if(runif(1)<r){
        psi.c3.1[j,a] = tmp.prop[1]
        psi.c3.2[j,a] = tmp.prop[2]
        psi3[,a] = tmp.prop
        accept3[a] = accept3[a]+1
      }else{
        psi.c3.1[j,a] = psi.c3.1[k,a]
        psi.c3.2[j,a] = psi.c3.2[k,a]
      }
    }
    
    pi = pi.runs[j,]
    psi1 = rbind(psi.c1.1[j,],psi.c1.2[j,])
    psi2 = rbind(psi.c2.1[j,],psi.c2.2[j,])
    psi3 = rbind(psi.c3.1[j,],psi.c3.2[j,])
    psi4 = rbind(psi.c4.1[k,],psi.c4.2[k,])
    psi5 = rbind(psi.c5.1[k,],psi.c5.2[k,])
    for(a in 1:G){
      #tmp.prop = rnorm(2,psi4[,a]*100,sd=5)
      #tmp.prop = tmp.prop/sum(tmp.prop)
      tmp.prop = rdirichlet(1,alpha = 2*psi4[,a]+.2)
      psi4.prop = psi4
      psi4.prop[,a] = tmp.prop
      
      prop.post = likelihood(M[,k],samp.size,G,pi,psi1,psi2,psi3,psi4.prop,psi5) + log(ddirichlet(tmp.prop,alpha = rep(1,2)))
      past.post = likelihood(M[,k],samp.size,G,pi,psi1,psi2,psi3,psi4,psi5) + log(ddirichlet(psi4[,a],alpha = rep(1,2)))
      ##ppropos = log(ddirichlet(psi4[,a],alpha = tmp.prop)) - log(ddirichlet(tmp.prop,alpha = psi4[,a]))
      ppropos = log(ddirichlet(psi4[,a], 2*psi4[,a]+.2)) - log(ddirichlet(tmp.prop, 2*psi4[,a]+.2))
      
      r = exp(prop.post - past.post + ppropos)
      
      if(runif(1)<r){
        psi.c4.1[j,a] = tmp.prop[1]
        psi.c4.2[j,a] = tmp.prop[2]
        psi4[,a] = tmp.prop
        accept4[a] = accept4[a]+1
      }else{
        psi.c4.1[j,a] = psi.c4.1[k,a]
        psi.c4.2[j,a] = psi.c4.2[k,a]
      }
    }
    
    pi = pi.runs[j,]
    psi1 = rbind(psi.c1.1[j,],psi.c1.2[j,])
    psi2 = rbind(psi.c2.1[j,],psi.c2.2[j,])
    psi3 = rbind(psi.c3.1[j,],psi.c3.2[j,])
    psi4 = rbind(psi.c4.1[j,],psi.c4.2[j,])
    psi5 = rbind(psi.c5.1[k,],psi.c5.2[k,])
    for(a in 1:G){
      
      tmp.prop = rdirichlet(1,alpha = 2*psi5[,a]+.2)
      #tmp.prop = rnorm(2,psi5[,a]*100,sd=5)
      #tmp.prop = tmp.prop/sum(tmp.prop)
      psi5.prop = psi5
      psi5.prop[,a] = tmp.prop
      
      prop.post = likelihood(M[,k],samp.size,G,pi,psi1,psi2,psi3,psi4,psi5.prop) + log(ddirichlet(tmp.prop,alpha = rep(1,2)))
      past.post = likelihood(M[,k],samp.size,G,pi,psi1,psi2,psi3,psi4,psi5) + log(ddirichlet(psi5[,a],alpha = rep(1,2)))
      ##ppropos = log(ddirichlet(psi5[,a],alpha = 2* proposal.weights.5/samp.size)) - log(ddirichlet(tmp.prop,alpha = 2* proposal.weights.5/samp.size))
      ppropos = log(ddirichlet(psi5[,a], 2*psi5[,a]+.2)) - log(ddirichlet(tmp.prop, 2*psi5[,a]+.2))
      
      r = exp(prop.post - past.post + ppropos)
      
      if(runif(1)<r){
        psi.c5.1[j,a] = tmp.prop[1]
        psi.c5.2[j,a] = tmp.prop[2]
        psi5[,a] = tmp.prop
        accept5[a] = accept5[a]+1
      }else{
        psi.c5.1[j,a] = psi.c5.1[k,a]
        psi.c5.2[j,a] = psi.c5.2[k,a]
      }
    }
  
    

    psi1 = rbind(psi.c1.1[j,],psi.c1.2[j,])
    psi2 = rbind(psi.c2.1[j,],psi.c2.2[j,])
    psi3 = rbind(psi.c3.1[j,],psi.c3.2[j,])
    psi4 = rbind(psi.c4.1[j,],psi.c4.2[j,])
    psi5 = rbind(psi.c5.1[j,],psi.c5.2[j,])
    pi = pi.runs[j,]
    
    
    P.marg[j,] = Marginal.Probs(G,pi,psi1,psi2,psi3,psi4,psi5)
    if((j %% 100) == 0){print(j)}
  }
  
  P.full = matrix(NA, nrow = N.runs, ncol = 32)
  for(i in 1:N.runs){
    psi1 = rbind(psi.c1.1[i,],psi.c1.2[i,])
    psi2 = rbind(psi.c2.1[i,],psi.c2.2[i,])
    psi3 = rbind(psi.c3.1[i,],psi.c3.2[i,])
    psi4 = rbind(psi.c4.1[i,],psi.c4.2[i,])
    psi5 = rbind(psi.c5.1[i,],psi.c5.2[i,])
    pi = pi.runs[i,]
    
    
    P.full[i,] = Full.Probs(G,pi,psi1,psi2,psi3,psi4,psi5)
  }
  Acc.Probs = list(accept.counts,accept.pi,accept1,accept2,accept3,accept4,accept5)
  return.list <- list("Full" = P.full, "Marginal" = P.marg, "Counts" = M,"Acc.Probs" = Acc.Probs)
  return(return.list)
}


TwoWay_LCM <- function(counts, G, samp.size, N.runs,prop.weights){
  ##Initializing vectors to hold the class membership probabilities
  ##Each is a matrix-- "c*" represents variable * where .1 or .2 represent the levels
  ##Each column represents a latent class
  psi.c1.1 = matrix(NA,nrow = N.runs, ncol = G)
  psi.c1.2 = matrix(NA,nrow = N.runs, ncol = G)
  
  psi.c2.1 = matrix(NA,nrow = N.runs, ncol = G)
  psi.c2.2 = matrix(NA,nrow = N.runs, ncol = G)
  
  psi.c3.1 = matrix(NA,nrow = N.runs, ncol = G)
  psi.c3.2 = matrix(NA,nrow = N.runs, ncol = G)
  
  psi.c4.1 = matrix(NA,nrow = N.runs, ncol = G)
  psi.c4.2 = matrix(NA,nrow = N.runs, ncol = G)
  
  psi.c5.1 = matrix(NA,nrow = N.runs, ncol = G)
  psi.c5.2 = matrix(NA,nrow = N.runs, ncol = G)
  
  pi.runs = matrix(NA, nrow = N.runs, ncol = G)
  
  ##Startings values
  ##Just assuming equal membership across latent classes
  pi.runs[1,] = rep(1/G,G)
  psi.c1.1[1,] = rep(.5,G)
  psi.c1.2[1,] = rep(.5,G)
  
  psi.c2.1[1,] = rep(.5,G)
  psi.c2.2[1,] = rep(.5,G)
  
  psi.c3.1[1,] = rep(.5,G)
  psi.c3.2[1,] = rep(.5,G)
  
  psi.c4.1[1,] = rep(.5,G)
  psi.c4.2[1,] = rep(.5,G)
  
  psi.c5.1[1,] = rep(.5,G)
  psi.c5.2[1,] = rep(.5,G)
  
  M = counts

  ##Now, beging the M-H algorithm
  
  for(j in 2:N.runs){
    k = j-1
    
    ##Mixing proportions
    prop.mix = rdirichlet(1,alpha = 5*pi.runs[k,]+.5)
    
    psi1 = rbind(psi.c1.1[k,],psi.c1.2[k,])
    psi2 = rbind(psi.c2.1[k,],psi.c2.2[k,])
    psi3 = rbind(psi.c3.1[k,],psi.c3.2[k,])
    psi4 = rbind(psi.c4.1[k,],psi.c4.2[k,])
    psi5 = rbind(psi.c5.1[k,],psi.c5.2[k,])
    
    prop.post = likelihood(M,samp.size,G,prop.mix,psi1,psi2,psi3,psi4,psi5) + log(ddirichlet(prop.mix,alpha = rep(1/G,G)))
    past.post = likelihood(M,samp.size,G,pi.runs[k,],psi1,psi2,psi3,psi4,psi5) + log(ddirichlet(prop.mix,alpha = rep(1/G,G)))
    ppropos = log(ddirichlet(pi.runs[k,],alpha = prop.mix)) - log(ddirichlet(prop.mix,alpha = pi.runs[k,]))
    
    r = exp(prop.post - past.post + ppropos)
    
    if(runif(1)<r){
      pi.runs[j,] = prop.mix
    }else{
      pi.runs[j,] = pi.runs[k,]
    }
    
    ##Now, for psi1
    pi = pi.runs[j,]
    psi1 = rbind(psi.c1.1[k,],psi.c1.2[k,])
    psi2 = rbind(psi.c2.1[k,],psi.c2.2[k,])
    psi3 = rbind(psi.c3.1[k,],psi.c3.2[k,])
    psi4 = rbind(psi.c4.1[k,],psi.c4.2[k,])
    psi5 = rbind(psi.c5.1[k,],psi.c5.2[k,])
    for(a in 1:G){
      
      tmp.prop = rdirichlet(1,alpha = 5*c(psi.c1.1[k,a],psi.c1.2[k,a])+.5)
      psi1.prop = psi1
      psi1.prop[,a] = tmp.prop
      
      prop.post = likelihood(M,samp.size,G,pi,psi1.prop,psi2,psi3,psi4,psi5) + log(ddirichlet(tmp.prop,alpha = rep(1/2,2)))
      past.post = likelihood(M,samp.size,G,pi,psi1,psi2,psi3,psi4,psi5) + log(ddirichlet(psi1[,a],alpha = rep(1/2,2)))
      ppropos = log(ddirichlet(psi1[,a],alpha = tmp.prop)) - log(ddirichlet(tmp.prop,alpha = psi1[,a]))
      
      r = exp(prop.post - past.post + ppropos)
      
      if(runif(1)<r){
        psi.c1.1[j,a] = tmp.prop[1]
        psi.c1.2[j,a] = tmp.prop[2]
        psi1[,a] = tmp.prop
      }else{
        psi.c1.1[j,a] = psi.c1.1[k,a]
        psi.c1.2[j,a] = psi.c1.2[k,a]
      }
    }
    
    pi = pi.runs[j,]
    psi1 = rbind(psi.c1.1[j,],psi.c1.2[j,])
    psi2 = rbind(psi.c2.1[k,],psi.c2.2[k,])
    psi3 = rbind(psi.c3.1[k,],psi.c3.2[k,])
    psi4 = rbind(psi.c4.1[k,],psi.c4.2[k,])
    psi5 = rbind(psi.c5.1[k,],psi.c5.2[k,])
    for(a in 1:G){
      
      tmp.prop = rdirichlet(1,alpha = 5*c(psi.c2.1[k,a],psi.c2.2[k,a])+.5)
      psi2.prop = psi2
      psi2.prop[,a] = tmp.prop
      
      prop.post = likelihood(M,samp.size,G,pi,psi1,psi2.prop,psi3,psi4,psi5) + log(ddirichlet(tmp.prop,alpha = rep(1/2,2)))
      past.post = likelihood(M,samp.size,G,pi,psi1,psi2,psi3,psi4,psi5) + log(ddirichlet(psi2[,a],alpha = rep(1/2,2)))
      ppropos = log(ddirichlet(psi2[,a],alpha = tmp.prop)) - log(ddirichlet(tmp.prop,alpha = psi2[,a]))
      
      r = exp(prop.post - past.post + ppropos)
      
      if(runif(1)<r){
        psi.c2.1[j,a] = tmp.prop[1]
        psi.c2.2[j,a] = tmp.prop[2]
        psi2[,a] = tmp.prop
      }else{
        psi.c2.1[j,a] = psi.c2.1[k,a]
        psi.c2.2[j,a] = psi.c2.2[k,a]
      }
    }
    
    pi = pi.runs[j,]
    psi1 = rbind(psi.c1.1[j,],psi.c1.2[j,])
    psi2 = rbind(psi.c2.1[j,],psi.c2.2[j,])
    psi3 = rbind(psi.c3.1[k,],psi.c3.2[k,])
    psi4 = rbind(psi.c4.1[k,],psi.c4.2[k,])
    psi5 = rbind(psi.c5.1[k,],psi.c5.2[k,])
    for(a in 1:G){
      
      tmp.prop = rdirichlet(1,alpha = 5*c(psi.c3.1[k,a],psi.c3.2[k,a])+.5)
      psi3.prop = psi3
      psi3.prop[,a] = tmp.prop
      
      prop.post = likelihood(M,samp.size,G,pi,psi1,psi2,psi3.prop,psi4,psi5) + log(ddirichlet(tmp.prop,alpha = rep(1/2,2)))
      past.post = likelihood(M,samp.size,G,pi,psi1,psi2,psi3,psi4,psi5) + log(ddirichlet(psi3[,a],alpha = rep(1/2,2)))
      ppropos = log(ddirichlet(psi3[,a],alpha = tmp.prop)) - log(ddirichlet(tmp.prop,alpha = psi3[,a]))
      
      r = exp(prop.post - past.post + ppropos)
      
      if(runif(1)<r){
        psi.c3.1[j,a] = tmp.prop[1]
        psi.c3.2[j,a] = tmp.prop[2]
        psi3[,a] = tmp.prop
      }else{
        psi.c3.1[j,a] = psi.c3.1[k,a]
        psi.c3.2[j,a] = psi.c3.2[k,a]
      }
    }
    
    pi = pi.runs[j,]
    psi1 = rbind(psi.c1.1[j,],psi.c1.2[j,])
    psi2 = rbind(psi.c2.1[j,],psi.c2.2[j,])
    psi3 = rbind(psi.c3.1[j,],psi.c3.2[j,])
    psi4 = rbind(psi.c4.1[k,],psi.c4.2[k,])
    psi5 = rbind(psi.c5.1[k,],psi.c5.2[k,])
    for(a in 1:G){
      
      tmp.prop = rdirichlet(1,alpha = 5*c(psi.c4.1[k,a],psi.c4.2[k,a])+.5)
      psi4.prop = psi4
      psi4.prop[,a] = tmp.prop
      
      prop.post = likelihood(M,samp.size,G,pi,psi1,psi2,psi3,psi4.prop,psi5) + log(ddirichlet(tmp.prop,alpha = rep(1/2,2)))
      past.post = likelihood(M,samp.size,G,pi,psi1,psi2,psi3,psi4,psi5) + log(ddirichlet(psi4[,a],alpha = rep(1/2,2)))
      ppropos = log(ddirichlet(psi4[,a],alpha = tmp.prop)) - log(ddirichlet(tmp.prop,alpha = psi4[,a]))
      
      r = exp(prop.post - past.post + ppropos)
      
      if(runif(1)<r){
        psi.c4.1[j,a] = tmp.prop[1]
        psi.c4.2[j,a] = tmp.prop[2]
        psi4[,a] = tmp.prop
      }else{
        psi.c4.1[j,a] = psi.c4.1[k,a]
        psi.c4.2[j,a] = psi.c4.2[k,a]
      }
    }
    
    pi = pi.runs[j,]
    psi1 = rbind(psi.c1.1[j,],psi.c1.2[j,])
    psi2 = rbind(psi.c2.1[j,],psi.c2.2[j,])
    psi3 = rbind(psi.c3.1[j,],psi.c3.2[j,])
    psi4 = rbind(psi.c4.1[j,],psi.c4.2[j,])
    psi5 = rbind(psi.c5.1[k,],psi.c5.2[k,])
    for(a in 1:G){
      
      tmp.prop = rdirichlet(1,alpha = 5*c(psi.c5.1[k,a],psi.c5.2[k,a])+.5)
      psi5.prop = psi5
      psi5.prop[,a] = tmp.prop
      
      prop.post = likelihood(M,samp.size,G,pi,psi1,psi2,psi3,psi4,psi5.prop) + log(ddirichlet(tmp.prop,alpha = rep(1/2,2)))
      past.post = likelihood(M,samp.size,G,pi,psi1,psi2,psi3,psi4,psi5) + log(ddirichlet(psi5[,a],alpha = rep(1/2,2)))
      ppropos = log(ddirichlet(psi5[,a],alpha = tmp.prop)) - log(ddirichlet(tmp.prop,alpha = psi5[,a]))
      
      r = exp(prop.post - past.post + ppropos)
      
      if(runif(1)<r){
        psi.c5.1[j,a] = tmp.prop[1]
        psi.c5.2[j,a] = tmp.prop[2]
        psi5[,a] = tmp.prop
      }else{
        psi.c5.1[j,a] = psi.c5.1[k,a]
        psi.c5.2[j,a] = psi.c5.2[k,a]
      }
    }

    print(j)}
  
  P.marg = matrix(NA, nrow = N.runs, ncol = 40)
  for(i in 1:N.runs){
    psi1 = rbind(psi.c1.1[i,],psi.c1.2[i,])
    psi2 = rbind(psi.c2.1[i,],psi.c2.2[i,])
    psi3 = rbind(psi.c3.1[i,],psi.c3.2[i,])
    psi4 = rbind(psi.c4.1[i,],psi.c4.2[i,])
    psi5 = rbind(psi.c5.1[i,],psi.c5.2[i,])
    pi = pi.runs[i,]
    
    
    P.marg[i,] = Marginal.Probs(G,pi,psi1,psi2,psi3,psi4,psi5)
  }
  
  P.full = matrix(NA, nrow = N.runs, ncol = 32)
  for(i in 1:N.runs){
    psi1 = rbind(psi.c1.1[i,],psi.c1.2[i,])
    psi2 = rbind(psi.c2.1[i,],psi.c2.2[i,])
    psi3 = rbind(psi.c3.1[i,],psi.c3.2[i,])
    psi4 = rbind(psi.c4.1[i,],psi.c4.2[i,])
    psi5 = rbind(psi.c5.1[i,],psi.c5.2[i,])
    pi = pi.runs[i,]
    
    
    P.full[i,] = Full.Probs(G,pi,psi1,psi2,psi3,psi4,psi5)
  }
  return.list <- list("Full" = P.full, "Marginal" = P.marg, "Counts" = M,Psi1 = psi1, Psi2 = psi2, Psi3 = psi3, Psi4 = psi4, Psi5 = psi5, Pi = pi)
  return(return.list)
}
