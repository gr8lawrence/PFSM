## This set of functions are used for initiating random matrices
## Assume M = GC + E, we need a function to get G, C and M based on the G and C and a corresponding E

# Required packages -------------------------------------------------------

library(extraDistr)
## functions to generate a random C
## uniform distribution is not "varied" enough
## use a Dirichlet distributin
# get_C <- function(n_cell_types, n_subjects) {
#   
#   # get a random positive matrix and then divide each element by row sum:
#   C = matrix(rdirichlet(n=n_subjects, alpha=rep(1, n_cell_types)), 
#              nrow=n_cell_types, 
#              ncol=n_subjects,
#              byrow=FALSE)
# 
#   return(C)
# 
# }


# Helper functions for submatrices ----------------------------------------

## functions for simulating G
## markers
sim_marker_G <- function(n_total_genes, n_cell_types, n_markers, 
                         n_new_genes, n_good_cell_types) {
  
  if (n_markers + n_new_genes != n_total_genes) {
    stop("Sum of known genes and new genes has to equal the number of total genes.")
  }
  
  # simulate G
  G = matrix(0L, n_total_genes, n_cell_types)
  G[seq_len(n_markers), seq_len(n_good_cell_types)] = sim_marker_GEP(n_markers, n_good_cell_types)
  G[-seq_len(n_markers), -seq_len(n_good_cell_types)] = sample(1:(2 * n_markers), n_new_genes)
  
  return(G)
  
}

sim_marker_GEP <- function(n_genes, n_cell_types) {
  
  G = matrix(0L, n_genes, n_cell_types)
  G = t(apply(G, 1, function(x) choose_one(x)))
  
  # the gene expression of each row is different: reduce data redundancy
  return(diag(sample(1:(2 * n_genes), n_genes)) %*% G)
  
}

## random expression
sim_random_GEP <- function(n_genes, n_cell_types, method="chisq") {
  
  if (method == "chisq") {
    G = matrix(rchisq(n_genes * n_cell_types, df=2), n_genes, n_cell_types)
  }
  
  return(G)
}

## helper function to choose one 1 per each input vector of 0s
choose_one <- function(vec_0) {
  
  len = length(vec_0)
  ind_1 = sample(1:len, 1)
  vec_0[ind_1] = 1
  
  return(vec_0)
  
}


# Generating matrices -----------------------------------------------------

## functions to generate a random C (old version)
get_C <- function(n_cell_types, n_subjects) {

  # get a random positive matrix and then divide each element by row sum:
  C = matrix(runif(n_cell_types * n_subjects),
             nrow=n_cell_types,
             ncol=n_subjects)
  C = apply(C, 2, function(x) x/sum(x))

  return(C)

}

## functions to get a random G
## generate a mixture of markers and non-markers
get_G <- function(n_genes, n_subjects, n_cell_types, n_markers, n_good_cell_types, marker_ratio=0.5) {
  
  # a mixture of markers and random
  if (!is.numeric(marker_ratio)) {
    stop("Marker ratio needs to be a number.")
  } else if (marker_ratio < 0) {
    stop("Marker ratio needs to be between 0 and 1.")
  } else if (marker_ratio > 1) {
    stop("Marker ratio needs to be between 0 and 1.")
  }
  
  # two other counts 
  n_new_genes = n_genes - n_markers
  n_bad_cell_types = n_cell_types - n_good_cell_types
  
  # get the numbers of known and unknown genes for each type
  n_known_marker_genes = floor(n_markers * marker_ratio)
  n_known_random_genes = n_markers - n_known_marker_genes
  n_unknown_marker_genes = floor(n_new_genes * marker_ratio)
  n_unknown_random_genes = n_new_genes - n_unknown_marker_genes
  
  n_total_marker_genes = n_known_marker_genes + n_unknown_marker_genes
  n_total_random_genes = n_known_random_genes + n_unknown_random_genes
  
  # simulate the marker and random pieces
  marker_piece = sim_marker_G(n_total_marker_genes, n_cell_types, n_known_marker_genes,
                              n_unknown_marker_genes, n_good_cell_types)
  random_piece = sim_random_GEP(n_total_random_genes, n_cell_types)
  
  # combine the pieces
  G_sim = rbind(
    marker_piece[seq_len(n_known_marker_genes), ],
    random_piece[seq_len(n_known_random_genes), ],
    marker_piece[n_known_marker_genes + seq_len(n_unknown_marker_genes), ],
    random_piece[n_known_random_genes + seq_len(n_unknown_random_genes), ]
  )
  
  return(G_sim)
  
}