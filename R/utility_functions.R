## the following functions are for the utilities of our algorithm
## return the objective function calculated from a certain update
PSMF_obj_function <- function(M, G_hat, C_hat, G_0, Delta, alpha, xi, beta) {
  Delta_c = 1 - Delta
  return(
    1/2 * norm(M - G_hat %*% C_hat, type="F") + 
      1/2 * alpha * norm(Delta * G_hat - G_0, type="F") +
      1/2 * xi * norm(Delta_c * G_hat, type="F") +
      1/2 * beta * norm(colSums(C_hat) - 1, type="2")
    )
}

## return a delta matrix according to 
get_Delta <- function(n_cell_types, n_good_cell_types, n_genes, n_known_genes) {
  
  Delta = rbind(
    matrix(
      rep(c(rep(1, n_good_cell_types), rep(0, n_cell_types - n_good_cell_types)), n_known_genes),
      n_known_genes,
      n_cell_types,
      byrow=TRUE 
    ),
    matrix(
      0L,
      n_genes - n_known_genes,
      n_cell_types
    )
  )
  
  return(Delta)
  
}


## here are some functions that are useful for calculating matrices
# matrix difference
matrix_diff <- function(A, B) {
  return(A - B)
}

# get the norm of matrix differences
# A is a list
matrix_diff_norm <- function(A, B, i) {
  return(A %>% pluck(i) %>% matrix_diff(B) %>% norm(., type="F"))
}

# get the max element difference
get_max_element_diff <- function(A, B, i) {
  return(A %>% pluck(i) %>% max_element_abs_diff(B))
}

# get the max column sum difference from 1
get_max_column_sum_diff_from_one <- function(A, i) {
  return(A %>% pluck(i) %>% max_column_sum_diff_from_one)
}

# max absolute diffence of elements
max_element_abs_diff <- function(A, B) {
  return(max(abs(A - B)))
}

# max column sum difference from 1
max_column_sum_diff_from_one <- function(A) {
  return(max(abs(colSums(A) - 1)))
}
