## the following functions are for the utilities of our algorithm
## return the objective function calculated from a certain update
PSMF_obj_function <- function(M, G_hat, C_hat, G_0, Delta, alpha, xi, beta) {
  Delta_c = 1 - Delta
  return(
    1/2 * norm(M - G_hat %*% C_hat, type="F") + 
      1/2 * alpha * norm(Delta * G - G_0, type="F") +
      1/2 * xi * norm(Delta_c * G, type="F") +
      1/2 * beta * norm(colSums(C_hat) - 1, type="F")
    )
}


