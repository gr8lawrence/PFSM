pacman::p_load(tidyverse, stringr, purrr, furrr)
source("../core_algorithm_functions.R")
source("../matrix_initiation_functions.R")
source("../utility_functions.R")

# this file is 
# some parameters that are fixed
n_genes <- 50
n_cell_types <- 5
n_known_genes <- 30
n_good_cell_types <- 4

n = c(50, 100, 200, 500, 1000)
m = 100 # 100 simulations for each sample size
k = 1:5 # five initial values



