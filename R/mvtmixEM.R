#' Train a t mixture model.
#'
#' Train a t mixture model on data x.
#'
#' @param x Matrix of data.
#' @param lambda0 Initial EM value for lambda, the vector of mixing probabilities. Must sum to 1.
#' @param mu0 Initial EM value for mu, the list of mean vectors.
#' @param sigma0 Initial EM value for sigma, the list of covariance matricies.
#' @param nu0 Initial EM value for nu, the vector of degrees of freedom. Must be greater than 0. If nu0 is specified and initialization != 'none' then values for nu will not be generated.
#' @param G Number of components. Must be an integer greater than 1. Must specify exactly one of G or G_range.
#' @param G_range Vector of G values to test (each value must be an integer greater than 1). The best value will be chosen via a selection criterion. Must specify exactly one of G or G_range.
#' @param tolerance: EM algorithm convergence tolerance.
#' @param max_iterations: EM algorithm max iterations.
#' @param optimize_nu: Logical value indicating whether to optimize for nu. If false, nu0 must be specified.
#' @param EM_variant: Variant of the EM algorithm to use for model fitting. Must be one of 'EM', 'SEM', 'SAEM'.
#' @param delta_seq: Vector of delta values for the SAEM variant of EM. Must be a decreasing sequence with each value <= 1.
#' @param initialization: Initialization method. Must be one of 'none', 'random', or 'kmeans'. If not 'none', overwrites lambda0, mu0, and sigma0 and, if nu0 has not specified, each element of nu0 is randomally sampled from Unif(0,100).
#' @param sc: Selection criteria to choose optimal number of components if G_range has been specified. Must be one of 'aic', 'bic', or 'icl'.
#' @param stoppage: Stoppage method for EM algorithm. Must be one of 'likelihood' (i.e. convergence is determined on the loglikelihood) or 'aitken' (i.e. convergence is determined on the aitkens acceleration values).
#'
#' @return Returns a list of class tMM with elements:
#' \itemize{
#'   \item lambda - Vector of mixing probabilities.
#'   \item mu - List of mean vectors.
#'   \item sigma - List of covariance matricies.
#'   \item nu - Vector of degrees of freedom.
#'   \item G -  Number of components. If G_range is specified then G is the value in G_range that yields the optimal sc value.
#'   \item DG - Number of parameters estimated.
#'   \item AIC - AIC value.
#'   \item BIC - BIC value.
#'   \item ICL - ICL value.
#'   \item loglik - Vector of the log likelihood values at each iteration.
#'   \item logLG - Final log likelihood value.
#'   \item initial - List of the lambda, mu, sigma and nu used to initialize EM.
#' }
#'
#' @importFrom mvtnorm dmvt
#'
#' @export
mvtmixEM <-
  function(x,
           lambda0 = NULL,
           mu0 = NULL,
           sigma0 = NULL,
           nu0 = NULL,
           G = NULL,
           G_range = NULL,
           tolerance = 1e-8,
           max_iterations = 1000,
           optimize_nu = TRUE,
           EM_variant = 'EM',
           delta_seq = (max_iterations:1) / max_iterations,
           initialization = 'kmeans',
           sc = 'BIC',
           stoppage = 'aitken') {
    x <- data.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

    if ((!is.null(G) &
         !is.null(G_range)) |
        (is.null(G) &
         is.null(G_range)))
      stop('Must specify exactly one of G or G_range')
    if (is.null(initialization) &
        length(G_range) > 1)
      stop('Must specify the exact number of clusters (G) if starting points are included')
    if (!is.null(G))
      G_range = G

    models <- lapply(G_range, function(G) {
      if (initialization != 'none') {
        initial <-
          mvtmixInitializeEM(
            x = x,
            nu0 = nu0,
            G = G,
            initialization = initialization
          )

        mu0 <- initial$mu0
        sigma0 <- initial$sigma0
        nu0 <- initial$nu0
        lambda0 <- initial$lambda0
      }

      loglik <- vector(length = max_iterations + 1)
      a <- vector(length = max_iterations)
      la <- rep(Inf, max_iterations)

      loglik[1] <- mvtmixll(
        x = x,
        mu = mu0,
        sigma = sigma0,
        nu = nu0,
        lambda = lambda0
      )

      mu <- mu0
      sigma <- sigma0
      nu <- nu0
      lambda <- lambda0

      if (EM_variant == 'SAEM') {
        E <-
          mvtmixEstep(
            x = x,
            mu = mu,
            sigma = sigma,
            nu = nu,
            lambda = lambda,
            G = G,
            n = n,
            p = p,
            EM_variant = 'EM'
          )

        M <-
          mvtmixMstep(
            x = x,
            G = G,
            n = n,
            p = p,
            E = E,
            optimize_nu = optimize_nu,
            EM_variant = 'EM'
          )

        muEM <- M$mu
        sigmaEM <- M$sigma
        nuEM <- M$nu
        lambdaEM <- M$lambda

        E <-
          mvtmixEstep(
            x = x,
            mu = mu,
            sigma = sigma,
            nu = nu,
            lambda = lambda,
            G = G,
            n = n,
            p = p,
            EM_variant = 'SEM'
          )

        M <-
          mvtmixMstep(
            x = x,
            G = G,
            n = n,
            p = p,
            E = E,
            optimize_nu = optimize_nu,
            EM_variant = 'SEM'
          )

        muSEM <- M$mu
        sigmaSEM <- M$sigma
        nuSEM <- M$nu
        lambdaSEM <- M$lambda

        mu <- lapply(1:G, function(g) {
          ((1 - delta_seq[1]) * muEM[[g]]) + (delta_seq[1] * muSEM[[g]])
        })
        sigma <- lapply(1:G, function(g) {
          ((1 - delta_seq[1]) * sigmaEM[[g]]) + (delta_seq[1] * sigmaSEM[[g]])
        })
        #rounding as weighted average may not yield an integer
        nu <-
          round(((1 - delta_seq[1]) * nuEM) + (delta_seq[1] * nuSEM))
        lambda <-
          ((1 - delta_seq[1]) * lambdaEM) + (delta_seq[1] * lambdaSEM)

        loglik[2] <-
          mvtmixll(
            x = x,
            mu = mu,
            sigma = sigma,
            nu = nu,
            lambda = lambda
          )

        i <- 3

        while (i <= max_iterations) {
          E <-
            mvtmixEstep(
              x = x,
              mu = mu,
              sigma = sigma,
              nu = nu,
              lambda = lambda,
              G = G,
              n = n,
              p = p,
              EM_variant = 'EM'
            )

          M <-
            mvtmixMstep(
              x = x,
              G = G,
              n = n,
              p = p,
              E = E,
              optimize_nu = optimize_nu,
              EM_variant = 'EM'
            )

          muEM <- M$mu
          sigmaEM <- M$sigma
          nuEM <- M$nu
          lambdaEM <- M$lambda

          E <-
            mvtmixEstep(
              x = x,
              mu = mu,
              sigma = sigma,
              nu = nu,
              lambda = lambda,
              G = G,
              n = n,
              p = p,
              EM_variant = 'SEM'
            )

          M <-
            mvtmixMstep(
              x = x,
              G = G,
              n = n,
              p = p,
              E = E,
              optimize_nu = optimize_nu,
              EM_variant = 'SEM'
            )

          muSEM <- M$mu
          sigmaSEM <- M$sigma
          nuSEM <- M$nu
          lambdaSEM <- M$lambda

          mu <- lapply(1:G, function(g) {
            ((1 - delta_seq[i - 1]) * muEM[[g]]) + (delta_seq[i - 1] * muSEM[[g]])
          })

          sigma <- lapply(1:G, function(g) {
            ((1 - delta_seq[i - 1]) * sigmaEM[[g]]) + (delta_seq[i - 1] * sigmaSEM[[g]])
          })

          #rounding as weighted average may not yield an integer
          nu <-
            round(((1 - delta_seq[i - 1]) * nuEM) + (delta_seq[i - 1] * nuSEM))

          lambda <-
            ((1 - delta_seq[i - 1]) * lambdaEM) + (delta_seq[i - 1] * lambdaSEM)

          loglik[i] <-
            mvtmixll(
              x = x,
              mu = mu,
              sigma = sigma,
              nu = nu,
              lambda = lambda
            )

          if (stoppage == 'aitken') {
            a[i - 1] <-
              (loglik[i] - loglik[i - 1]) / (loglik[i - 1] - loglik[i - 2])
            la[i] <-
              loglik[i - 1] + (loglik[i] - loglik[i - 1]) / (1 - a[i - 1])

            if (abs(la[i] - la[i - 1]) < tolerance) {
              break
            }
          }

          else if (stoppage == 'likelihood') {
            if (abs(loglik[i] - loglik[i - 1]) < tolerance) {
              break
            }
          }
          i <- i + 1
        }

        loglik <- loglik[1:(i - 1)]

      } else {
        E <-
          mvtmixEstep(
            x = x,
            mu = mu,
            sigma = sigma,
            nu = nu,
            lambda = lambda,
            G = G,
            n = n,
            p = p,
            EM_variant = EM_variant
          )

        M <-
          mvtmixMstep(
            x = x,
            G = G,
            n = n,
            p = p,
            E = E,
            optimize_nu = optimize_nu,
            EM_variant = EM_variant
          )

        mu <- M$mu
        sigma <- M$sigma
        nu <- M$nu
        lambda <- M$lambda

        loglik[2] <- mvtmixll(
          x = x,
          mu = mu,
          sigma = sigma,
          nu = nu,
          lambda = lambda
        )

        i <- 3

        while (i <= max_iterations) {
          E <-
            mvtmixEstep(
              x = x,
              mu = mu,
              sigma = sigma,
              nu = nu,
              lambda = lambda,
              G = G,
              n = n,
              p = p,
              EM_variant = EM_variant
            )

          M = mvtmixMstep(
            x = x,
            G = G,
            n = n,
            p = p,
            E = E,
            optimize_nu = optimize_nu,
            EM_variant = EM_variant
          )

          mu <- M$mu
          sigma <- M$sigma
          nu <- M$nu
          lambda <- M$lambda

          loglik[i] <-
            mvtmixll(
              x = x,
              mu = mu,
              sigma = sigma,
              nu = nu,
              lambda = lambda
            )

          if (stoppage == 'aitken') {
            a[i - 1] <-
              (loglik[i] - loglik[i - 1]) / (loglik[i - 1] - loglik[i - 2])
            la[i] <-
              loglik[i - 1] + (loglik[i] - loglik[i - 1]) / (1 - a[i - 1])

            if (abs(la[i] - la[i - 1]) < tolerance)
              break
          }

          else if (stoppage == 'likelihood') {
            if (abs(loglik[i] - loglik[i - 1]) < tolerance)
              break
          }
          i <- i + 1
        }
        loglik <- loglik[1:(i - 1)]
      }

      #D is the number of parameters, G - 1 lambdas (the last one is then uniquely determined), G mus, G sigmas, G nus
      D <- (G - 1) + G + G + G
      ENT <- -sum(E$EZ_X * log(E$EZ_X))

      #BIC here is a pseudo-Bayesian criterion, since it does not depend on the prior distribution of theta
      AIC <- -2 * loglik[length(loglik)] + (2 * D)
      BIC <- -2 * loglik[length(loglik)] + (D * log(n))
      ICL <- BIC + 2 * ENT

      model <-
        list(
          lambda = lambda,
          mu = mu,
          sigma = sigma,
          nu = nu,
          G = G,
          DG = D,
          AIC = AIC,
          BIC = BIC,
          ICL = ICL,
          loglik = loglik,
          logLG = loglik[length(loglik)],
          initial = initial
        )

      class(model) <- 'tMM'
      model
    })

    ics <- list(
      AIC = sapply(1:length(models), function(i)
        models[[i]]$AIC),
      BIC = sapply(1:length(models), function(i)
        models[[i]]$BIC),
      ICL = sapply(1:length(models), function(i)
        models[[i]]$ICL),
      logLG = sapply(1:length(models), function(i)
        models[[i]]$logLG),
      DG = sapply(1:length(models), function(i)
        models[[i]]$DG)
    )

    #p_opt = 2 * coef(lm(logLGs ~ DGs))[2]
    #SH = logLGs + (p_opt * DGs)
    return(models[[which.min(ics[[sc]])]])
  }

#E-step
mvtmixEstep <-
  function (x,
            lambda,
            mu,
            sigma,
            nu,
            G,
            n,
            p,
            EM_variant = 'EM') {
    #Expected value of Z Given X; as a matrix
    EZ_X <- t(sapply(1:n, function (i) {
      sapply(1:G, function (g) {
        lambda[g] * mvtnorm::dmvt(
          x = x[i,],
          delta = mu[[g]],
          sigma = sigma[[g]],
          df = nu[g],
          log = FALSE
        ) /
          sum(sapply(1:G, function (j) {
            lambda[j] * mvtnorm::dmvt(
              x = x[i,],
              delta = mu[[j]],
              sigma = sigma[[j]],
              df = nu[j],
              log = FALSE
            )
          }))
      })
    }))

    if (EM_variant == 'SEM') {
      #k is sampled
      k <- sapply(1:n, function(i) {
        sample(x = 1:length(EZ_X[i,]),
               size = 1,
               prob = EZ_X[i,])
      })

      #Expected value of U given X and Z; as a vector
      EU_XZ <- sapply(1:n, function (i) {
        (nu[k[i]] + p) / (nu[k[i]] + t(x[i, ] - mu[[k[i]]]) %*%
                            solve(sigma[[k[i]]]) %*% (x[i, ] - mu[[k[i]]]))
      })
    } else {
      #Expected value of U given X and Z; as a vector
      EU_XZ <- t(sapply(1:n, function (i) {
        sapply(1:G, function(g) {
          (nu[g] + p) / (nu[g] + t(x[i,] - mu[[g]]) %*%
                           solve(sigma[[g]]) %*% (x[i,] - mu[[g]]))
        })
      }))
    }
    #Returning all expectations in a list
    list(EZ_X = EZ_X, EU_XZ = EU_XZ)
  }

#M-step
mvtmixMstep <-
  function (x,
            G,
            n,
            p,
            E,
            optimize_nu = TRUE,
            EM_variant = 'EM') {
    EZ_X <- E$EZ_X
    EU_XZ <- E$EU_XZ

    #lambda parameter update
    lambda <- sapply(1:G, function (g) {
      mean(EZ_X[, g])
    })

    if (EM_variant == 'SEM') {
      #mu parameter updates
      mu <- lapply(1:G, function (g) {
        rowSums(sapply(1:n, function (i) {
          EZ_X[, g][i] * EU_XZ[i] * x[i,]
        })) /
          sum(EZ_X[, g] * EU_XZ)
      })

      #sigma parameter updates
      sigma <- lapply(1:G, function (g) {
        rowSums(sapply(1:n, function (i) {
          EZ_X[, g][i] * EU_XZ[i] * (x[i,] - mu[[g]]) %*% t(x[i,] - mu[[g]])
        }, simplify = "array"), dims = 2) /
          sum(EZ_X[, g])
      })

      if (optimize_nu == TRUE) {
        nu <- sapply(1:G, function(g) {
          suppressWarnings(uniroot(
            f = function(nu) {
              -digamma(0.5 * nu) + log(0.5 * nu) + 1 + (sum(EZ_X[, g] * (log(EU_XZ) - EU_XZ)) / sum(EZ_X[, g])) +
                digamma(0.5 * (nu + p)) - log(0.5 * (nu + p))
            },
            lower = 1,
            upper = 100,
            extendInt = 'yes'
          ))$root
        })
      }
    } else {
      #mu parameter updates
      mu <- lapply(1:G, function (g) {
        rowSums(sapply(1:n, function (i) {
          EZ_X[, g][i] * EU_XZ[, g][i] * x[i,]
        })) /
          sum(EZ_X[, g] * EU_XZ[, g])
      })

      #sigma parameter updates
      sigma <- lapply(1:G, function (g) {
        rowSums(sapply(1:n, function (i) {
          EZ_X[, g][i] * EU_XZ[, g][i] * (x[i,] - mu[[g]]) %*% t(x[i,] - mu[[g]])
        }, simplify = "array"), dims = 2) /
          sum(EZ_X[, g])
      })

      #nu parameter updates
      if (optimize_nu == TRUE) {
        nu <- sapply(1:G, function(g) {
          suppressWarnings(uniroot(
            f = function(nu) {
              -digamma(0.5 * nu) + log(0.5 * nu) + 1 + (sum(EZ_X[, g] * (log(EU_XZ[, g]) - EU_XZ[, g])) / sum(EZ_X[, g])) +
                digamma(0.5 * (nu + p)) - log(0.5 * (nu + p))
            },
            lower = 1,
            upper = 100,
            extendInt = 'yes'
          ))$root
        })
      }
    }
    list(
      mu = mu,
      sigma = sigma,
      nu = nu,
      lambda = lambda
    )
  }
