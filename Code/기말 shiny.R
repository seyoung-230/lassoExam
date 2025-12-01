library(shiny)

# ----------------------------------------------------------
# LASSO CDA 함수 (BIC 계산 포함)
# ----------------------------------------------------------

ST = function(a, lambda) {
  if (a > lambda) return(a - lambda)
  if (a < -lambda) return(a + lambda)
  return(0)
}

lasso_cda = function(X, y, lambda, max_iter = 1000, tol = 1e-5, scaling = TRUE)
{
  n = nrow(X)
  p = ncol(X)
  
  if (scaling) {
    X = scale(X)
    y = y - mean(y)
  }
  
  beta = rep(0, p)
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

# ----------------------------------------------------------
# (g)의 시뮬레이션 데이터를 미리 생성
# ----------------------------------------------------------

set.seed(1)
n = 200
p = 10
X = matrix(rnorm(n * p), n, p)
beta_true = c(2, 1, rep(0, 8))
y = X %*% beta_true + rnorm(n, 0, 0.1)

# ----------------------------------------------------------
# UI (문제에서 지정한 그대로)
# ----------------------------------------------------------

ui = fluidPage(
  
  titlePanel("LASSO Coefficient Plot (CDA Method)"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        "lambda",
        "Select Lambda (λ):",
        min = 1e-10, max = 1, value = 1
      )
    ),
    
    mainPanel(
      plotOutput("coef_plot")   
    )
  )
)

# ----------------------------------------------------------
# SERVER
# ----------------------------------------------------------

server = function(input, output) {
  
  output$coef_plot = renderPlot({
    
    # λ 가져오기
    lam = input$lambda
    
    # LASSO CDA 실행 (문제 조건: scaling = TRUE)
    beta_hat = lasso_cda(X, y, lambda = lam, scaling = TRUE)
    
    # 변수 이름
    var_names = paste0("X", 1:p)
    
    # 막대그래프 그리기
    barplot(
      beta_hat,
      names.arg = var_names,
      col = "#0072B2",
      main = paste("LASSO Coefficients (lambda =", lam, ")"),
      ylab = "Estimated Beta",
      xlab = "Predictor Variables",
      las = 2
    )
  })
}

# ----------------------------------------------------------
# 앱 실행
# ----------------------------------------------------------

shinyApp(ui = ui, server = server)
