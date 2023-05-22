library(splines2)
library(quadprog)
library(boot)
library(SemiPar)
library(plotrix)
standardize <- function(x, y) {
     sx <- scale(x); sy <- scale(y)
     x_bar <- attr(sx, "scaled:center"); x_sd <- attr(sx, "scaled:scale")
     y_bar <- attr(sy, "scaled:center"); y_sd <- attr(sy, "scaled:scale")
     list(sx = sx, x_bar = x_bar, x_sd = x_sd,
                      sy = sy, y_bar = y_bar, y_sd = y_sd)
 }
qp_solve <- function(standardized, A_mat, x_mat, ...) {
     d_vec <- with(standardized, crossprod(sx, sy)) 
     D_mat <- crossprod(standardized$sx) 
     sbeta <- quadprog::solve.QP(D_mat, d_vec, A_mat, ...)$solution
     beta0 <- with(standardized,
                                        as.numeric(crossprod(x_bar * y_sd / x_sd, sbeta)) + y_bar)
     beta1 <- with(standardized, y_sd / x_sd * sbeta)
     list(x_mat = x_mat,
                      coefficients = c("(Intercept)" = beta0, beta1),
                      fitted = as.numeric(x_mat %*% beta1) + beta0)
 }
mono_fit <- function(x, y, monotone = c("increasing", "decreasing"), ...) {
     monotone <- match.arg(monotone)
     x_mat <- iSpline(x, ..., intercept = TRUE) 
     standardized <- standardize(x_mat, y) 
     A_mat <- diag(ncol(standardized$sx)) 
     if (monotone == "decreasing") A_mat <- -A_mat  
     qp_solve(standardized, A_mat, x_mat)
 }
data(age.income, package = "SemiPar")
dfs <- c(6)
mfits <- lapply(dfs, function(d) {
     with(age.income, mono_fit(age, log.income, df = d))
 })



#bootstrap sampling 100 times -- get 100 βs


n_bootstrap <- 1000
sample_size <- 10
dfs <- c(6)
mfits_list <- vector("list", n_bootstrap)

fit_model <- function(data, indices, dfs) {
  boot_sample <- data[indices, ]
  
  mfits <- lapply(dfs, function(d) {
    with(boot_sample, mono_fit(age, log.income, df = d))
  })
  
  coefficients <- sapply(mfits, function(mfit) mfit$coefficients)
  return(coefficients)
}

boot_results <- boot(data = age.income, statistic = fit_model, R = n_bootstrap, dfs = dfs)

coefficients <- boot_results$coefficients

conf_intervals <- boot.ci(boot_results, type = "bca")

for (i in 1:nrow(conf_intervals$basic)) {
  coef_name <- names(age.income)[i]  
  coef_est <- coefficients[i]  
  conf_interval <- conf_intervals$basic[i, , drop = FALSE]  
  
 
  cat("Coefficient:", coef_name, "\n")
  cat("Estimate:", coef_est, "\n")
  cat("Confidence Interval:", conf_interval, "\n")
  cat("\n")
}
