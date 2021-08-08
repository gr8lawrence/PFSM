## the following functions are for updating C
## solving C directly by minimizing the objective function

# main update function
update_C <- function(M, G_old, C_old, beta, method) {
  
  if (method == 'direct') {
    C_new = update_C_direct(
      M=M,
      G_old=G_old,
      beta=beta,
      eta=eta
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
update_C_direct <- function(M, G_old, beta, eta) {
  
  # dimension 
  n_cell_types = nrow(G_old) 
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
  
  C_new = matrix(0, nrow(C_old), ncol(C_old)) 
  
  for (i in 1:ncol(C_old)) {
    
    C_new[, i] = update_C_auxiliary_vec(c_old=C_old[, i], 
                                        m=M[, i], 
                                        G_old=G_old, 
                                        beta=beta)
    
  }
  
  return(C_new)
  
}

# update for each column of C
update_C_auxiliary_vec <- function(c_old, m, G_old, beta) {
  
  # dimension
  n_cell_types = length(c_old)
  
  # constants
  Gtm = t(G_old) %*% m
  GtGc = t(G) %*% G %*% c_old
  beta_c_old_sum = beta * sum(c_old)
  
  # update C
  c_new = vector(mode="numeric", length=n_cell_types)
  
  for (k in 1:n_cell_types) {
    c_new[k] = (Gtm[k] + beta)/(GtGc[k] + beta * beta_c_old_sum)
  }
  
  return(c_new)
  
}

## the following are for updating G
# main update function
update_G <- function(M, G_old, C_old, beta, method) {
  
  if (method == 'direct') {
   
  } else if (method == 'auxiliary') {
   
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
  G_new = matrix(0, n_genes, n_cell_types)
  
  for (j in 1:n_genes) {
    if (j <= n_markers) {
    
      G_new[j, ] = (M[j, ] %*% t(C_new) + alpha * G_0[j, ] %*% V) %*%
        solve(CCt + alpha * V + xi * V_c)
      
    } else {
      
      G_new[j, ] = (M[j, ] %*% t(C_new)) %*% solve(CCt + alpha * I_p) 
      
    }
  }
  
  return(G_new)
  
}

# matrix level update for G
update_G_auxiliary <- function(G_old, M, C_new, G_0, n_markers, n_good_cell_types) {
  
}

# vector level update for G
update_G_auxiliary_vec <- function(g_old, m_j, C_new, gamma_j, V, alpha, xi) {
  
  # dimension
  n_cell_types = length(g_old) 
  
  # constants
  V_c = diag(1, n_cell_types, n_cell_types) - V
  numerator = m_j %*% t(C_new) + alpha * gamma_j %*% V
  denominator = g_old %*% (C_new %*% t(C_new) + alpha * V + xi * V_c)
  
  # update G
  g_new = vector(mode="numeric", length=n_cell_types)
  for (r in 1:n_cell_types) {
    
    g_new[r] = g_old[r] * (numerator[r]/denominator[r])
    
  }
  
  return(g_new)
  
}