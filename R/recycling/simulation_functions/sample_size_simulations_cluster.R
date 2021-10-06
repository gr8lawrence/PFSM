pacman::p_load(tidyverse, stringr, purrr, furrr, doParallel, patchwork)
#setwd("/Users/gr8lawrence/Documents/Dissertation/conv_opt_code/cleaned_algorithm/PFSM/R")
source("../core_algorithm_functions.R")
source("../matrix_initiation_functions.R")
source("../utility_functions.R")
plan(multisession, workers=4)
# set.seed(100)
furrr_options(seed=TRUE)

# this file is for simulation using different sample sizes
# some parameters that are fixed
n_genes <- 50
n_cell_types <- 5
n_known_genes <- 30
n_good_cell_types <- 4

# 312,500 * 2 = 625,000 runs of the algorithm total
n_seq = c(50, 100, 200, 500, 1000) # five different sample sizes
m_seq = 1:100 # 100 simulations for each sample size
k_seq = 1:5 # five initial values
method_seq = c("direct", "auxiliary") # two methods

# tune grid (values of hyperparameters determined by experiment simulations)
alpha_v = 2^seq(8, 16, 2)
xi_v = 2^seq(4, 8, 1)
beta_v = 2^seq(12, 20, 2)

## result data frame
total_df <- expand_grid(
  n_subjects = n_seq,
  m = m_seq,
  k = k_seq,
  a = alpha_v,
  x = xi_v,
  b = beta_v,
  method = method_seq
  )

dim(total_df) # 625,000 rows

## constants
Delta <- rbind(
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

## simulation function
run_simulation <- function(M, G0, G_init, C_init, a, x, b, n_subjects, method) {
  
  # record the solutions
  sol = PSMF_solve(
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
    method = method
  )
  
  G_hat = sol$G_hat
  C_hat = sol$C_hat
  res = sol$res[length(sol$res)]
  obj_val = PSMF_obj_function(M, 
                              G_hat, 
                              C_hat, 
                              G0, 
                              Delta, 
                              a, 
                              x, 
                              b)
  
  # record the values
  return(
    list(
      obj = obj_val,
      G_hat = G_hat,
      C_hat = C_hat,
      res = res
    ) 
  )
  
}

true_matrix_df <- expand_grid(n_subjects = n_seq,
                              m = m_seq) %>% 
  mutate(
    G_true = future_map(.x=n_subjects, ~get_G(n_genes, .x, n_cell_types, n_known_genes, n_good_cell_types),
                        .options = furrr_options(seed=123)),
    C_true = future_map(.x=n_subjects, ~get_C(n_cell_types, .x),
                        .options = furrr_options(seed=123)),
    G0 = future_map(.x=G_true, ~Delta * .x),
    E = future_map(.x=n_subjects, 
                   function(.x) {matrix(runif(n_genes * .x), n_genes, .x)},
                   .options=furrr_options(seed=123)),
    M = future_pmap(.l=list(G = G_true, C = C_true, E = E), 
                    function(G, C, E) {
                      G %*% C + E 
                    })
  )

all_df <- total_df %>% 
  left_join(true_matrix_df, by=c("n_subjects", "m")) %>% 
  # group_by(n_subjects, m, k, method) %>% 
  # filter(row_number() <= 2) %>% 
  # ungroup() %>% 
  mutate(
    G_init = future_map(.x=n_subjects, ~get_G(n_genes, .x, n_cell_types, n_known_genes, n_good_cell_types),
                        .options = furrr_options(seed=123)),
    C_init = future_map(.x=n_subjects, ~get_C(n_cell_types, .x),
                        .options = furrr_options(seed=123)),
    sol = future_pmap(
      .l=list(M, G0, G_init, C_init, a, x, b, n_subjects, method),
      run_simulation
    )
  ) 

rm(true_matrix_df, total_df)

best_obj_df <- all_df %>% 
  mutate(
    obj = future_map_dbl(.x=sol, ~ .x$obj)
  ) %>% 
  group_by(n_subjects, m, method) %>% 
  arrange(obj) %>% 
  filter(row_number() == 1) %>% 
  mutate(
    C_res = future_map2_dbl(.x=sol, .y=C_true, ~ norm(matrix_diff(.x %>% pluck(3), .y), type="F"))
  )

best_M_res_df <- all_df %>% 
  mutate(
    M_res = future_map_dbl(.x=sol, ~ .x$res)
  ) %>% 
  group_by(n_subjects, m, method) %>% 
  arrange(M_res) %>% 
  filter(row_number() == 1) %>% 
  mutate(
    C_res = future_map2_dbl(.x=sol, .y=C_true, ~ norm(matrix_diff(.x %>% pluck(3), .y), type="F"))
  )

best_C_res_df <- all_df %>% 
  mutate(
    C_res = future_map2_dbl(.x=sol, .y=C_true, ~ norm(matrix_diff(.x %>% pluck(3), .y), type="F"))
  ) %>% 
  group_by(n_subjects, m, method) %>% 
  arrange(C_res) %>% 
  filter(row_number() == 1)

plot_simulation_boxplots <- function(df) {
  
  df %>% 
    ggplot(aes(x=n_subjects, y=C_res, group=n_subjects, col=n_subjects)) +
    geom_boxplot() + 
    facet_wrap(.~method, scales="free") +
    theme(legend.position="right")
  
}

pdf("sample_size_simulation_plot.pdf")
p1 <- plot_simulation_boxplots(best_obj_df)
p2 <- plot_simulation_boxplots(best_M_res_df)
p3 <- plot_simulation_boxplots(best_C_res_df)

p1 / p2 / p3
dev.off()
                                       
                                       