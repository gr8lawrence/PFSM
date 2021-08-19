## the following functions are for updating C
## solving C directly by minimizing the objective function

# main update function
update_C <- function(M, G_old, C_old, beta, method) {
  
  if (method == 'direct') {
    C_new = update_C_direct(
      M=M,
      G_old=G_old,
      beta=beta
    )
  } else if (method == 'auxiliary') {
    C_new = update_C_auxiliary(
      C_old=C_old,
      M=M, 
      G_old=G_old, 
      beta=beta
    )
  } else {
    stop("Method must be either 'direct' or 'auxiliary'.")
  }
  
  return(C_new)
  
}

# two methods of update
update_C_direct <- function(M, G_old, beta) {
  
  # dimension 
  n_cell_types = ncol(G_old) 
  n_subjects = ncol(M)

  # update C directly
  C_new = solve(t(G_old) %*% G_old + 
                  beta * matrix(1, n_cell_types, n_cell_types)) %*%
    (t(G_old) %*% M + beta * matrix(1, n_cell_types, n_subjects))
  
  return(C_new)
  
}

## solving C by minimizing an auxiliary function
# standalone function to be called from outside
update_C_auxiliary <- function(C_old, M, G_old, beta) {
  
  # dimension
  n_cell_types = nrow(C_old) 
  n_subjects = ncol(C_old)
  
  # constants
  num_matrix = t(G_old) %*% M + matrix(1L, n_cell_types, n_subjects)
  denom_matrix = t(G_old) %*% G_old %*% C_old + matrix(1L, n_cell_types, n_subjects)
  
  # updates (using the element-wise products)
  C_new = C_old * num_matrix/denom_matrix 
  
  return(C_new)
  
}

## the following are for updating G
# main update function
update_G <- function(G_old, M, C_new, G_0, n_markers, n_good_cell_types, alpha, xi, method) {
  
  if (method == 'direct') {
    G_new = update_G_direct(
      M=M, 
      C_new=C_new, 
      G_0=G_0, 
      alpha=alpha, 
      xi=xi, 
      n_markers=n_markers, 
      n_good_cell_types=n_good_cell_types
    )
  } else if (method == 'auxiliary') {
    G_new = update_G_auxiliary(
      G_old=G_old, 
      M=M, 
      C_new=C_new, 
      G_0=G_0, 
      n_markers=n_markers, 
      n_good_cell_types=n_good_cell_types, 
      alpha=alpha, 
      xi=xi
    )
  } else {
    stop("Method must be either 'direct' or 'auxiliary'.")
  }
  
  return(G_new)
  
}


update_G_direct <- function(M, C_new, G_0, alpha, xi, n_markers, n_good_cell_types) {
  
  # dimension
  n_genes = nrow(M)
  n_new_genes = n_genes - n_markers
  n_cell_types = nrow(C_new)
  n_bad_cell_types = n_cell_types - n_good_cell_types
  
  # constants
  I_p = diag(1, n_cell_types, n_cell_types)
  V = diag(c(rep(1, n_good_cell_types), rep(0, n_bad_cell_types)))
  V_c = I_p - V
  CCt = C_new %*% t(C_new) 
  
  # update G directly
  G_new = matrix(0L, n_genes, n_cell_types)
  G_new[1:n_markers, ] = (M[1:n_markers, ] %*% t(C_new) + alpha * G_0[1:n_markers, ] %*% V) %*% solve(CCt + alpha * V + xi * V_c)
  G_new[(n_markers + 1):n_genes, ] = (M[(n_markers + 1):n_genes, ] %*% t(C_new)) %*% solve(CCt + xi * I_p) 

  return(G_new)
  
}

# matrix level update for G
update_G_auxiliary <- function(G_old, M, C_new, G_0, n_markers, n_good_cell_types, alpha, xi) {
  
  # dimension
  n_genes = nrow(M)
  n_new_genes = n_genes - n_markers
  n_cell_types = nrow(C_new)
  n_bad_cell_types = n_cell_types - n_good_cell_types
  
  # constants
  I_p = diag(1, n_cell_types, n_cell_types)
  V = diag(c(rep(1, n_good_cell_types), rep(0, n_bad_cell_types)))
  V_c = I_p - V
  CCt = C_new %*% t(C_new) 
  G_old_known = G_old[1:n_markers, ]
  G_old_new_genes = G_old[(n_markers + 1):n_genes, ]
  
  # update G (part-by-part)
  G_new = matrix(0L, n_genes, n_cell_types)
  # for the part with known markers
  num_matrix_known = M[1:n_markers, ] %*% t(C_new) + alpha * G_0[1:n_markers, ] %*% V
  denom_matrix_known = G_old_known %*% (CCt + alpha * V + xi * V_c)
  G_new[1:n_markers, ] = G_old_known * num_matrix_known/denom_matrix_known
  # for the part with unknown markers
  num_matrix_new_genes = M[(n_markers + 1):n_genes, ] %*% t(C_new) 
  denom_matrix_new_genes = G_old_new_genes %*% (CCt + alpha * V + xi * V_c)
  G_new[(n_markers + 1):n_genes, ] = G_old_new_genes * num_matrix_new_genes/denom_matrix_new_genes
  
  return(G_new)
  
}