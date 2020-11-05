#' Generate EM initial values for a t mixture model.
#'
#' Generate EM initial values for a t mixture model with G components.
#'
#' @param x Matrix of data.
#' @param nu0 Vector of degrees of freedom. Must be greater than 0. If nu0 is specified then values for nu will not be generated.
#' @param G Number of components. Must be greater than 1.
#' @param initialization Initialization method. Must be either 'none', random' or 'kmeans'. If initialization != 'none' and nu0 has not been specified, each element of nu0 is randomally sampled from Unif(0,100).
#'
#' @return Returns a list with elements:
#' \itemize{
#'   \item lambda0 - Vector of mixing probabilities.
#'   \item mu0 - List of mean vectors.
#'   \item sigma0 - List of covariance matricies.
#'   \item nu0 - Vector of degrees of freedom.
#' }
#
#' @export
mvtmixInitializeEM <-
  function(x,
           nu0 = NULL,
           G = NULL,
           initialization = 'kmeans') {
    cov2 <- function(x) {
      #this is defined to deal with the division by n-1 when n=1
      p <- ncol(x)
      if (nrow(x) == 1) {
        matrix(rep(0, p * p), nrow = p)
      } else {
        cov(x)
      }
    }

    x <- data.matrix(x)
    n <- nrow(x)

    if (initialization == 'kmeans') {
      km <- kmeans(x = x, centers = G)
      mu0 <- lapply(1:G, function(g)
        km$centers[g, ])
      sigma0 <- lapply(1:G, function(g)
        cov2(x[km$cluster == g, ]))
      lambda0 <- km$size / n

      if (is.null(nu0)) {
        nu0 <- runif(n = G, min = 0, max = 100)
      }

    }

    if (initialization == 'random') {
      clusters <- sample(x = 1:G ,
                         size = n,
                         replace = TRUE)
      mu0 <- lapply(1:G, function(g)
        colMeans(x[clusters == g, ]))
      sigma0 <- lapply(1:G, function(g)
        cov2(x[clusters == g, ]))
      lambda0 <- sapply(1:G, function(g)
        mean(clusters == g))

      if (is.null(nu0)) {
        nu0 <- runif(n = G, min = 0, max = 100)
      }

    }

    list(
      lambda0 = lambda0,
      mu0 = mu0,
      sigma0 = sigma0,
      nu0 = nu0
    )
  }
