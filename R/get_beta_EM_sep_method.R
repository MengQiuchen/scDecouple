#' Perform EM-based deconvolution using the scDecouple method.
#'
#' This function performs EM-based deconvolution to estimate the cellular response using the scDecouple method.
#'
#' @param y A matrix of observed data, cell-by-PC, from the perturbation group.
#' @param number.component The number of components (clusters) to estimate.
#' @param em A list of component means for each cluster.
#' @param es A list of component covariance matrices for each cluster.
#' @param el A vector of component proportions for each cluster.
#' @param beta_init The initial estimate of the cellular response.
#' @param max.step The maximum number of EM iterations. Default is 100.
#' @param fix_es Logical, indicating whether to fix the component covariance matrices. Default is FALSE.
#'
#' @return A list containing the following results:
#'   \itemize{
#'     \item beta: The estimated cellular response.
#'     \item es: A list of updated component covariance matrices.
#'     \item el: A vector of updated component proportions.
#'     \item em: A list of updated component means.
#'     \item Q_after: The log-likelihood value after the EM iterations.
#'     \item ey: A matrix of component assignment probabilities for each observation.
#'   }
#'
#' @import mixtools
#'
#' @examples
#' # Example usage:
#' result <- get_beta_EM_sep_method(y, number.component, em, es, el, beta_init, max.step, fix_es)
#'
#' @export
get_beta_EM_sep_method <- function(y, number.component, em, es, el, beta_init,
                                   max.step = 100, fix_es = FALSE) {
  
  library(mixtools)
  
  if (!length(dim(y))) { # one-dim
    number.sample <- length(y)
    dim <- 1
  } else {
    number.sample <- dim(y)[1]
    dim <- dim(y)[2]
  }
  
  Q_before <- 0
  beta <- beta_init
  ey <- matrix(NA, nrow = number.sample, ncol = number.component)
  dy <- matrix(NA, nrow = number.sample, ncol = number.component)
  
  for (k in 1:max.step) {
    
    for (j in 1:number.component) {
      if (dim == 1) {
        dy[, j] <- (dnorm(data.matrix(y), as.matrix(em[[j]] +
                                                      beta), es[[j]])) * el[j]
      } else {
        dy[, j] <- (dmvnorm(data.matrix(y), as.matrix(em[[j]] +
                                                        beta), es[[j]])) * el[j]
      }
      
    }
    for (i in 1:number.sample) {
      sy <- sum(dy[i, ])
      if (sy == 0) {
        sy <- 1e-308
      }
      ey[i, ] <- dy[i, ] / sy
    }
    em_sum <- 0
    
    flag = TRUE
    for (j in 1:number.component) {
      el[j] <- sum(ey[, j]) / number.sample
      
      if (!sum(ey[, j])) {
        flag = FALSE
        print('WARNING: NOT CONVERGENT!')
        break
      }
      if (dim == 1) {
        em_sum <- em_sum + sum(ey[, j] * y) / sum(ey[, j])
      } else {
        em_sum <- em_sum + colSums(ey[, j] * y) / sum(ey[, j])
      }
      
    }
    if (!flag) {
      break
    }
    beta <- em_sum / number.component - colMeans(matrix(unlist(em),
                                                      ncol = dim, byrow = TRUE))
    if (!fix_es) {
      for (j in 1:number.component) {
        es_sum <- matrix(0, nrow = dim, ncol = dim)
        for (i in 1:number.sample) {
          if (dim == 1) {
            es_sum <- es_sum + ey[i, j] * (t(y[i] - t(as.matrix(em[[j]] +
                                                                  beta))) %*% (y[i] - t(as.matrix(em[[j]] +
                                                                                                    beta))))
          } else {
            es_sum <- es_sum + ey[i, j] * (t(y[i, ] - t(as.matrix(em[[j]] +
                                                                    beta))) %*% (y[i, ] - t(as.matrix(em[[j]] +
                                                                                                        beta))))
          }
          
        }
        es[[j]] <- es_sum / sum(ey[, j])
      }
    }
    Q_after <- sum(log(rowSums(dy)))
    if (abs(Q_after - Q_before) < 1e-06) {
      break
    } else {
      Q_before <- Q_after
    }
  }
  
  result <- list(beta = beta, es = es, el = el, em = em, Q_after = Q_after,
                 ey = ey)
  
  return(result)
}
