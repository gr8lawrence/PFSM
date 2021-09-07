## the follwing functions are for iterative solving the optimization problem
## paramaeters:
## G_0: partially known signature matrix (constructed separately)
## G_init: initial value for G
## alpha, xi, beta: cross-validated hyperparameters for penalization

source("./core_update_functions_new.R")
source("./utility_functions.R")

PSMF_solve_test <- function(M, G_0, G_init, C_init=NULL, n_markers, n_good_cell_types, 
                            alpha, xi, beta, max_iter=1e5, eps, method, Delta, G_true, C_true) {
  
  # parameters
  n_genes = nrow(M)
  n_subjects = ncol(M)
  n_cell_types = ncol(G_0)
  
  # initiate the counter
  n_iter = 0
  
  # inititate a residual vector
  obj <- vector()
  M_res <- vector()
  G_res <- vector()
  C_res <- vector()
  G_fixed_res <- vector()
  C_colsum_res <- vector()
  
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
      
      # calculate the residuals
      obj[n_iter + 1] = PSMF_obj_function(M, G_old, C_old, G_0, Delta, alpha, xi, beta)
      M_res[n_iter + 1] = norm(M - G_old %*% C_old, type="F")
      C_res[n_iter + 1] = norm(C_old - C_true, type="F")
      C_colsum_res[n_iter + 1] = norm(colSums(C_old) - 1, type="2")
      
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
      
      # calculate the residuals
      obj[n_iter + 1] = PSMF_obj_function(M, G_old, C_old, G_0, Delta, alpha, xi, beta)
      M_res[n_iter + 1] = norm(M - G_old %*% C_old, type="F")
      G_res[n_iter + 1] = norm(G_old - G_true, type="F")
      G_fixed_res[n_iter + 1] = norm(Delta * G_old - G_0, type="F")
      
    }
    
    n_iter = n_iter + 1
    
  }
  
  #print(paste0("The algorithm terminated after finishing ", n_iter + 1, " iterations."))
  
  return(list(G_hat = G_old, C_hat = C_old, n_iter = n_iter, obj=obj, M_res_vec = M_res,
              G_res_vec = G_res, C_res_vec = C_res, C_colsum_res_vec = C_colsum_res,
              G_fixed_res_vec = G_fixed_res))
  
}
