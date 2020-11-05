#' Sample from a t mixture model.
#'
#' Sample from a t mixture model with parameters lambda, mu, sigma, nu.
#'
#' @param n Number of points to sample.
#' @param lambda Vector of mixing probabilities. Must sum to 1.
#' @param mu List of mean vectors.
#' @param sigma List of covariance matricies.
#' @param nu Vector of degrees of freedom. Must be greater than 0.
#'
#' @return An n x p matrix in which each row is an independently generated realization from the specified t mixture model.
#'
#' @importFrom mvtnorm rmvt
#'
#' @export
rmvtmix <-
  function(n = 1,
           lambda = NULL,
           mu = NULL,
           sigma = NULL,
           nu = NULL) {
    G <- length(lambda)
    t(sapply(1:n, function(i) {
      Zi <- sample(x = 1:G,
                   size = 1 ,
                   prob = lambda)
      mvtnorm::rmvt(
        n = 1,
        delta = mu[[Zi]],
        sigma = sigma[[Zi]],
        df = nu[Zi]
      )
    }))
  }
