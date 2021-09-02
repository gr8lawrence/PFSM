## the follwing functions are for iterative solving the optimization problem
## paramaeters:
## G_0: partially known signature matrix (constructed separately)
## G_init: initial value for G
## alpha, xi, beta: cross-validated hyperparameters for penalization

source("./core_update_functions_new.R")

PSMF_solve <- function(M, G_0, G_init, C_init=NULL, n_markers, n_good_cell_types, alpha, xi, beta, max_iter=1e5, eps, method) {
  
  # parameters
  n_genes = nrow(M)
  n_subjects = ncol(M)
  n_cell_types = ncol(G_0)
  
  # initiate the counter
  n_iter = 0
  
  # inititate a residual vector
  res <- vector()
  
  # initiate matrices
  G_old = G_init
  if (method == 'auxiliary') {
    C_old = C_init
  } else {
    C_old = matrix(1e-8, n_cell_types, n_subjects)
  }  
  
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
      
      # replace the old C with new C
      C_old = C_new
      
      # calculate the residual
      res[n_iter + 1] = norm(M - G_old %*% C_old, type="F")
      
    } else if (n_iter %% 2 == 1) {
      # if n_iter is odd, update G
      G_new = update_G(G_old=G_old, 
                       M=M, 
                       C_new=C_old, 
                       G_0=G_0, 
                       n_markers=n_markers, 
                       n_good_cell_types=n_good_cell_types, 
                       alpha=alpha, 
                       xi=xi,
                       method=method)
      
      # get the relavitve change in norm for G
      G_change = norm(G_new - G_old, type="F")/norm(G_old, type="F")
      
      # replace the old G with new G
      G_old = G_new
      
      # calculate the residual
      res[n_iter + 1] = norm(M - G_old %*% C_old, type="F")
      
    }
    
    n_iter = n_iter + 1
    
  }
  
  #print(paste0("The algorithm terminated after finishing ", n_iter + 1, " iterations."))
  
  return(list(G_hat = G_old, C_hat = C_old, n_iter = n_iter, res_vec = res))
  
}
