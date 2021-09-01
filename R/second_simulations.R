## ----setup, include=FALSE-----------------------------------------------------------------------------
rm(list=ls())
pacman::p_load(tidyverse, stringr, purrr, furrr)
source("./core_algorithm_functions.R")
source("./matrix_initiation_functions.R")


## ----get.matrices-------------------------------------------------------------------------------------
## set a seed
set.seed(100)

## constants
n_genes <- 20
n_cell_types <- 5
n_subjects <- 50
n_known_genes <- 12
n_good_cell_types <- 4

## generate a G
G_true <- get_G(n_genes, n_subjects, n_cell_types, n_known_genes, n_good_cell_types)

## generate a C
C_true <- get_C(n_cell_types, n_subjects)

## generate an error matrix
E <- matrix(runif(n_genes * n_subjects), n_genes, n_subjects)

## generate an M
M <- G_true %*% C_true + E  

## get G0 from G
G0 <- matrix(0L, n_genes, n_cell_types)
G0[1:n_known_genes, 1:n_good_cell_types] <- G_true[1:n_known_genes, 1:n_good_cell_types]

# initial values
G_init <- get_G(n_genes, n_subjects, n_cell_types, n_known_genes, n_good_cell_types)
C_init <- get_C(n_cell_types, n_subjects)


## -----------------------------------------------------------------------------------------------------
## here are some functions that are useful 
# max absolute diffence of elements
max_element_abs_diff <- function(A, B) {
  return(max(abs(A - B)))
}

# max column sum difference from 1
max_column_sum_diff_from_one <- function(A) {
  return(max(abs(colSums(A) - 1)))
}

# max_column_sum_diff_from_one(matrix(rnorm(9), 3, 3))
# max_column_sum_diff_from_one(matrix(1/3, 3, 3))


## -----------------------------------------------------------------------------------------------------
alpha_v = 2^seq(2, 20, 2)
xi_v = 2^seq(-4, 4, 1)
beta_v = 2^seq(2, 20, 2)

## make a grid of parameters
tune_grid <- expand_grid(
  alpha = alpha_v,
  xi = xi_v,
  beta = beta_v
  )

tune_grid


## -----------------------------------------------------------------------------------------------------
results <- tune_grid %>%
  filter(row_number() <= 20) %>% 
  mutate(
    solve_direct = pmap(list(a = alpha, x = xi, b = beta),
    function(a, x, b) {
      PSMF_solve(
        M = M,
        G_0 = G0,
        G_init = G_init,
        C_init = C_init,
        n_markers = n_known_genes,
        n_good_cell_types = n_good_cell_types,
        alpha = a,
        xi = x,
        beta = b,
        eps = 1e-5,
        method = 'direct'
      )
    }),
    solve_auxiliary = pmap(list(a = alpha, x = xi, b = beta),
    function(a, x, b) {
      PSMF_solve(
        M = M,
        G_0 = G0,
        G_init = G_init,
        C_init = C_init,
        n_markers = n_known_genes,
        n_good_cell_types = n_good_cell_types,
        alpha = a,
        xi = x,
        beta = b,
        eps = 1e-5,
        method = 'auxiliary'
      )
    }),
    direct_G_res = map(solve_direct, ~. %>% pluck(1) %>% norm(. - G_true, type="F")),
    auxiliary_G_res = map(solve_auxiliary, ~. %>% pluck(1) %>% norm(. - G_true, type="F")),
    direct_C_res = map(solve_direct, ~. %>% pluck(2) %>% norm(. - C_true, type="F")),
    auxiliary_C_res = map(solve_auxiliary, ~. %>% pluck(2) %>% norm(. - C_true, type="F"))
  )

## -----------------------------------------------------------------------------------------------------
all_results <- results %>% 
  mutate(
    direct_G_max_difference = map(solve_direct, ~. %>% pluck(1) %>% max_element_abs_diff(G_true)),
    direct_C_max_difference = map(solve_direct, ~. %>% pluck(1) %>% max_element_abs_diff(C_true)),
    direct_C_row_difference = map(solve_direct, ~. %>% pluck(2) %>% max_column_sum_diff_from_one)
  ) 

## -----------------------------------------------------------------------------------------------------

