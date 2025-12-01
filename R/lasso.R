#' LASSO Regression via Coordinate Descent Algorithm (CDA)
#'
#' This function computes the LASSO regression coefficient estimates using
#' the coordinate descent algorithm (CDA). Optionally, the predictor matrix
#' can be standardized before estimation. The function also computes the
#' Bayesian Information Criterion (BIC) based on the fitted model.
#'
#' @param X A numeric matrix of dimension \code{n × p} representing the predictor variables.
#' @param y A numeric response vector of length \code{n}.
#' @param lambda A positive tuning parameter for the LASSO penalty.
#' @param max_iter Maximum number of iterations for the CDA algorithm (default: 1000).
#' @param tol Convergence tolerance. The algorithm stops when the maximum
#'   coordinate-wise change in coefficients is below this threshold (default: 1e-5).
#' @param scaling Logical; if \code{TRUE}, the columns of \code{X} are standardized
#'   and \code{y} is centered before estimation. If \code{FALSE}, estimation is
#'   performed on the original unscaled data (default: FALSE for BIC evaluation).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{beta}}{A numeric vector of estimated LASSO regression coefficients.}
#'   \item{\code{BIC}}{The Bayesian Information Criterion value for the fitted model.}
#' }
#'
#' @description
#' The algorithm updates each coefficient in turn using the soft-thresholding operator.
#' After convergence or reaching \code{max_iter}, the fitted values are obtained
#' and the BIC is computed as:
#' \deqn{ \mathrm{BIC} = n \log(\mathrm{RSS}/n) + df \log n, }
#' where \code{df} is the number of nonzero coefficients.
#'
#'
lasso_cda = function(X, y, lambda, max_iter = 1000, tol = 1e-5, scaling = TRUE)
{
  n = nrow(X)
  p = ncol(X)

  if (scaling) {
    X = scale(X)
    y = y - mean(y)
  }

  beta = rep(0, p)
  beta_old = beta

  r = y

  for (iter in 1:max_iter) {
    max_change = 0

    for (j in 1:p) {
      X_j = X[, j]
      X_j_norm2 = sum(X_j^2)
      rho = sum(X_j * (r + X_j * beta[j]))

      beta_old_j = beta[j]
      beta[j] = ST(rho / X_j_norm2, lambda / X_j_norm2)

      r = r + X_j * (beta_old_j - beta[j])

      max_change = max(max_change, abs(beta[j] - beta_old_j))
    }

    if (max_change < tol) break
  }

  y_hat = X %*% beta
  RSS = sum((y - y_hat)^2)
  sigma2 = RSS / n
  df = sum(beta != 0)

  BIC = n * log(sigma2) + df * log(n)
  return(list(beta = beta, BIC = BIC))
}
