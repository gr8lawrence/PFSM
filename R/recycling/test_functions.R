## this file is used to write functions that test all the updates functions
## on small samples only!!
## do not export this file's functions or use them on large matrices

## for direct solutions, we build two functions for vector- (row- or column-wise) updates
update_G_direct_test <- function(M, C_new, G_0, alpha, xi, n_markers, n_good_cell_types) {
  
  # dimension
  n_genes = nrow(M)
  n_new_genes = n_genes - n_markers
  n_cell_types = nrow(C_new)
  n_bad_cell_types = n_cell_types - n_good_cell_types
  
  # constants
  I_p = diag(rep(1, n_cell_types))
  V = diag(c(rep(1, n_good_cell_types), rep(0, n_bad_cell_types)))
  V_c = I_p - V
  CCt = C_new %*% t(C_new) 
  
  # update
  G_new = matrix(NaN, n_genes, n_cell_types)
  
  for (j in 1:n_genes) {
    
    if (j > n_markers) {
      V = diag(rep(0, n_cell_types))
      V_c = I_p
    }
    
    G_new[j, ] = (M[j, ] %*% t(C_new)  + alpha * G_0[j, ] %*% V) %*% solve(CCt + alpha * V + xi * V_c) 
    
  }
  
  return(G_new)
  
}

update_C_direct_test <- function(M, G_old, beta) {
  
  # dimension 
  n_cell_types = ncol(G_old) 
  n_subjects = ncol(M)
  
  # constant
  GtG = t(G_old) %*% G_old 
  
  # update
  C_new = matrix(NaN, n_cell_types, n_subjects)
  
  for (i in 1:n_subjects) {
    C_new[, i] = solve(GtG + beta) %*% (t(G) %*% M[, i] + beta) 
  }
  
  return(C_new)
  
}

## for auxiliary solutions, we build two functions for element-wise updates (small matrices only)
update_G_auxiliary_test <- function(G_old, M, C_new, G_0, n_markers, n_good_cell_types, alpha, xi) {
  
  
  # dimension
  n_genes = nrow(M)
  n_new_genes = n_genes - n_markers
  n_cell_types = nrow(C_new)
  n_bad_cell_types = n_cell_types - n_good_cell_types
  
  # constants
  CCt = C_new %*% t(C_new) 
  I_p = diag(rep(1, n_cell_types))
  V = diag(c(rep(1, n_good_cell_types), rep(0, n_bad_cell_types)))
  V_c = I_p - V
  
  # update G
  G_new = matrix(NaN, n_genes, n_cell_types)
  
  for (j in 1:n_genes) {
    
    if (j > n_markers) {
      V = diag(rep(0, n_cell_types))
      V_c = I_p
    }
    
    m_j = M[j, ]
    g0_j = G_0[j, ]
    numerator = m_j %*% t(C_new) + alpha * g0_j %*% V
    denominator = G_old[j, ] %*% (CCt + alpha * V + xi * V_c)  
    
    for (r in 1:n_cell_types) {
      
      G_new[j, r] = G_old[j, r] * numerator[r]/denominator[r]
      
    }
  }
  
  return(G_new)
  
}

update_C_auxiliary_test <- function(C_old, M, G_old, beta) {
  
  # dimensions 
  n_genes = ncol(G_old)
  n_cell_types = nrow(C_old) 
  n_subjects = ncol(C_old)
  
  # constants
  GtM = t(G_old) %*% M
  GtGC = t(G_old) %*% G_old %*% C  
  
  # build a new G matrix
  C_new <- matrix(NaN, n_cell_types, n_subjects)
  
  # constants
  for (k in 1:n_cell_types) {
    for (i in 1:n_subjects) {
      C_new[k, i] = C_old[k, i] * (GtM[k, i] + beta)/(GtGC[k, i] + beta)
    }
  }
  
  return(C_new)
    
} 

## wrapper functions for easy calling
update_G_test <- function(G_old, M, C_new, G_0, n_markers, n_good_cell_types, alpha, xi, method) {
  
  if (method == 'direct') {
    G_new = update_G_direct_test(
      M=M, 
      C_new=C_new, 
      G_0=G_0, 
      alpha=alpha, 
      xi=xi, 
      n_markers=n_markers, 
      n_good_cell_types=n_good_cell_types
    )
  } else if (method == 'auxiliary') {
    G_new = update_G_auxiliary_test(
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


update_C_test <- function(M, G_old, C_old, beta, method) {
  
  if (method == 'direct') {
    C_new = update_C_direct_test(
      M=M,
      G_old=G_old,
      beta=beta
    )
  } else if (method == 'auxiliary') {
    C_new = update_C_auxiliary_test(
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