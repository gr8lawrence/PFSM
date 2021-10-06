pacman::p_load(tidyverse, stringr, purrr, furrr, doParallel, patchwork)
setwd("/Users/gr8lawrence/Documents/Dissertation/conv_opt_code/cleaned_algorithm/PFSM/R")
source("./core_algorithm_functions.R")
source("./matrix_initiation_functions.R")
source("./utility_functions.R")
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
# n_seq = c(50, 100, 200, 500, 1000) # five different sample sizes
n_seq = c(50, 100) 
# m_seq = 1:100 # 100 simulations for each sample size
m_seq = 1:5
k_seq = 1:5 # five initial values
method_seq = c("direct", "auxiliary") # two methods

# tune grid (values of hyperparameters determined by experiment simulations)
alpha_v = 2^seq(8, 16, 2)
xi_v = 2^seq(4, 8, 1)
beta_v = 2^seq(12, 20, 2)

# testing only
alpha_v = 2^seq(8, 10, 2)
xi_v = 2^seq(4, 5, 1)
beta_v = 2^seq(12, 16, 2)

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
  group_by(n_subjects, m, k, method) %>% 
  filter(row_number() <= 2) %>% 
  ungroup() %>% 
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

## simulation data
# direct method
# foreach(i = m_seq, .combine = 'rbind') %do% {
#   function(i) {
#     
#     # result data frame
#     
#     # generate five Gs, Cs, and Es, each for one sample size, and then get Ms
#     G_ls <- lapply(n_ls, function(i) get_G(n_genes, i , n_cell_types, n_known_genes, n_good_cell_types))
#     C_ls <- lapply(n_ls, function(i) get_C(n_cell_types, i))    
#     E_ls <- lapply(n_ls, function(i) {matrix(runif(n_genes * i), n_genes, i)})
#     M_ls <- sapply(1:5, function(i) G_ls[[i]] %*% C_ls[[i]] + E_ls[[i]])
#     
#     # pick the best parameter for 5 diiferent initial values
#     for (j in 1:length(G_ls)) {
#       
#       res_df = tibble(
#         obj = rep(NaN, length(k_seq)),
#         G_res = rep(NaN, length(k_seq)),
#         C_res = rep(NaN, length(k_seq)),
#         M_res = rep(NaN, length(k_seq))
#       )
#       
#       tune_grid = expand_grid(
#         alpha = alpha_v,
#         xi = xi_v,
#         beta_v = beta_v
#       ) %>% mutate(
#         sol = pmap(.l=list(a = alpha, x = xi, b = beta),
#                    function(a, x, b){
#                      params = c(a, x, b)
#                      run_simulation(params, n_seq[i], "direct")
#                    }),
#         obj = map(.x=sol, ~.x %>% pluck(2)), # objective functions
#         est = map(.x=sol, ~.x %>% pluck(1)) # all estimates
#       ) %>% 
#         arrange(obj) %>% 
#         filter(row_number() == 1)
#     
#       
#     }
# }
# 
# # auxiliary method

                                       
                                       
                                       
                                       