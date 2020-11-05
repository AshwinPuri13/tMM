#' Complete data log likelihood of a t mixture model.
#'
#' Returns the complete data likelihood of a t mixture model with data x and parameters lambda, mu, sigma, nu.
#'
#' @param x Matrix of data.
#' @param lambda Vector of mixing probabilities. Must sum to 1.
#' @param mu List of mean vectors.
#' @param sigma List of covariance matricies.
#' @param nu Vector of degrees of freedom. Must be greater than 0.
#'
#' @return The complete data loglikelihood (numeric).
#'
#' @importFrom mvtnorm dmvt
#'
#' @export
mvtmixll <-
  function (x,
            lambda = NULL,
            mu = NULL,
            sigma = NULL,
            nu = NULL) {
    n <- nrow(x)
    G <- length(lambda)
    sum(sapply(1:n, function (i) {
      log(sum(sapply(1:G, function (g) {
        lambda[g] * mvtnorm::dmvt(
          x = x[i,],
          delta = mu[[g]],
          sigma = sigma[[g]],
          df = nu[g],
          log = FALSE
        )
      })))
    }))
  }
