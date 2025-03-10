#' Simulate data for a Cormack-Jolly-Seber Survival and Detection Model
#' @author Ryan N. Kinzer
#' @param n
#' @param phi
#' @param p
#'
#' @return
#' @export
#'
#' @examples
simulate_cjs_data <- function(n, phi, p) {
  # n_marks: Number of marked individuals
  # phi1, phi2, phi3: survival probabilities between time steps
  # p2, p3, p4: detection probabilities at recapture events
  
  # unknown state and detection observations
  z_mat <- matrix(NA,nrow = n, ncol = length(phi) + 1)
  z_mat[,1] <- rep(1, length(n))
  
  y_mat <- matrix(NA, nrow = n, ncol = length(p) + 1)
  y_mat[,1] <- rep(1, length(n))
  
  for(i in 1:length(phi)){
    z_mat[,i+1] <- rbinom(n, 1, phi[i]) * z_mat[,i]
    y_mat[,i+1] <- rbinom(n, 1, p[i]) * z_mat[,i+1]
  }
  
  ch <- apply(y_mat, 1, paste, collapse = "")
  data <- data.frame(ch = ch, freq = 1)
  
  return(data)
}
