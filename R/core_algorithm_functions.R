## the follwing functions are for iterative solving the optimization problem
## paramaeters:
## G0: partially known signature matrix (constructed separately)
## G_init: initial value for G
## alpha, xi, beta: cross-validated hyperparameters for penalization

source("./core_update_functions.R")

PSMF_solve <- function(M, G0, G_init, C_init=NULL, n_markers, n_good_cell_types, alpha, xi, beta, max_iter=1000, epsilon, method) {
  
  # parameters
  n_genes = nrow(M)
  n_subjects = ncol(M)
  n_cell_types = nrow(G_0)
  
  # initiate the counter
  n_iter = 0
  
  # initiate matrices
  G_old = G_init
  C_old = ifelse(method == 'auxiliary', C_init, matrix(1e-16, n_cell_types, n_subjects))
  
  # initiate the relative change in norm to infinity
  C_change = Inf
  G_change = Inf
    
  # start the algorithm
  while(n_iter < max_iter & min(C_change, G_change) >= eps) {
    
    if (n_iter %% 2 == 0) {
      # if n_iter is even, update C
      C_new = update_C(M=M, 
                       G_old=G_old, 
                       C_old=C_old, 
                       beta=beta, 
                       method=method)
      
      # get the relavitve change in norm for C
      C_change = norm(C_new - C_old, type="F")/norm(C_old, type="F")
      
    } else if (n_iter %% 2 == 0) {
      # if n_iter is odd, update G
      G_new = update_G(G_old=G_old, 
                       M=M, 
                       C_new=C_new, 
                       G_0=G_0, 
                       n_markers=n_markers, 
                       n_good_cell_types=n_good_cell_types, 
                       alpha=alpha, 
                       xi=xi)
      
      # get the relavitve change in norm for G
      G_change = norm(G_new - G_old, type="F")/norm(G_old, type="F")
      
    }
    
  }
  
}
