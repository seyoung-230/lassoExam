# 1
## (a)
# 변수선택을 가능하게 하는 penalty term/ penalty의 크기를 조절해서 모델 복잡도를 결정

## (b)
ST = function(a, lambda){
  if (a > lambda) return(a - lambda)
  if (a < -lambda) return(a + lambda)
  return(0)
}

ST(5,3)
ST(-2,3)

## (c)
l2_R = function(v) {
  s = 0
  for (i in 1:length(v)) {
    s = s + v[i] * v[i]
  }
  return(s)
}

library(Rcpp)

cppFunction('
double l2_cpp(NumericVector v) {
  int n = v.size();
  double s = 0;
  for (int i = 0; i < n; i++) {
    s += v[i] * v[i];
  }
  return s;
}
')

set.seed(1)
v = runif(100000)

library(microbenchmark)

microbenchmark(
  R_version = l2_R(v),
  Rcpp_version = l2_cpp(v),
  times = 50
)

## (d)
lasso_cda_r = function(X, y, lambda, max_iter = 1000, tol = 1e-5, scaling = TRUE)
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
  
  return(beta)
}

## (e)
lasso_cda_BIC = function(X, y, lambda, max_iter = 1000, tol = 1e-5, scaling = TRUE)
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

## (g)
n = 200
p = 10
set.seed(1)

X = matrix(rnorm(n * p), n, p)
beta_true = c(2, 1, 0, 0, 0, 0, 0, 0, 0, 0)
y = X %*% beta_true + rnorm(n, 0, 0.1)

lambda_list = c(1, 5, 10)
lambda_list = c(0.01, 1, 3, 5, 10)

result_beta = list()
result_BIC = numeric(length(lambda_list))

for (i in 1:length(lambda_list)) {
  lam = lambda_list[i]
  
  out = lasso_cda_BIC(X, y, lambda = lam, scaling = TRUE)
  
  result_beta[[i]] = out$beta
  result_BIC[i] = out$BIC
}


result_beta
result_BIC

lambda_list[which.min(result_BIC)]
