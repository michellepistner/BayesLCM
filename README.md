  <!-- badges: start -->
  [![R-CMD-check](https://github.com/michellepistner/BayesLCM/workflows/R-CMD-check/badge.svg)](https://github.com/michellepistner/BayesLCM/actions)
  <!-- badges: end -->

# The `PrivLCM` package

Welcome to the `PrivLCM` package. Here, we implement the approach of [Nixon et. al (2022)](https://arxiv.org/abs/2201.10545). 

Here, we provide the code to fit our differentially-private Bayesian latent class model for synthetic data creation or posterior inference. The main steps of our approach are:

1. Subselect a set of marginal counts and add noise using the [Geometric mechanism for differential privacy](https://arxiv.org/abs/0811.2841).
2. Post-process these counts to be able to estimate probabilities for any table cell. This is done using ideas from composite likelihood and the [Dirichlet Process Mixture of Product Multinomials model](https://pubmed.ncbi.nlm.nih.gov/23606777/). 

Please see the small vignette (*example.Rmd*) to get started!
