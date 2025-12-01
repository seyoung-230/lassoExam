# lassoExam

This package implements LASSO regression using the coordinate descent algorithm (CDA).
It also computes the Bayesian Information Criterion (BIC) for model selection.

## Usage Example

```r
n = 200
p = 10
set.seed(1)
X = matrix(rnorm(n * p), n, p)
beta = c(2, 1, 0, 0, 0, 0, 0, 0, 0, 0)
y = X %*% beta_true + rnorm(n, 0, 0.1)

lambda_list = c(0.01, 1, 3, 5, 10)
results = lapply(lambda_list, function(lam) {
  lasso_cda_BIC(X, y, lambda = lam, scaling = TRUE)
})

results

```
