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

new_x_points = seq(range(age.income$age)[1],range(age.income$age)[2],length = 10)
n <- nrow(age.income)
n_bootstrap <- 1000
ci_result = list()

for(j in seq_along(new_x_points)){
  
  newx = new_x_points[j]
  boot_means <- rep(NA, n_bootstrap)
  
  for (i in 1:n_bootstrap) {
    boot_indices <- sample(n, replace = TRUE)  
    boot_sample <- age.income[boot_indices, ]  
    mfitss <- lapply(dfs, function(d) {
      with(boot_sample, mono_fit(age, log.income, df = d))
    })
    
    new_x = predict(mfitss[[1]]$x_mat,newx = newx)
    new_fit = mfitss[[1]]$coefficients[1]+sum(mfitss[[1]]$coefficients[-1]*new_x)
    boot_means[i] = as.numeric(new_fit)
  }
  
  CI = quantile(boot_means,probs = c(0.025,0.975))
  ci_result[[j]] = data.frame(x = newx,lower = CI[1],upper = CI[2])
  
}

ci_result = do.call('rbind',ci_result)
rownames(ci_result) = NULL
ci_result

for (h in seq_along(new_x_points)){
  
  newx = new_x_points[h]
  boot_meanss <- rep(NA, n_bootstrap)
  
  for (i in 1:n_bootstrap) {
    boot_indices <- sample(n, replace = TRUE)  
    boot_sample <- age.income[boot_indices, ]  
    mfits <- lapply(dfs, function(d) {
      with(boot_sample, mono_fit(age, log.income, df = d))
    })
    
    new_x1 = predict(mfits[[1]]$x_mat, newx = newx)
    new_fit1 = mfits[[1]]$coefficients[1] + sum(mfits[[1]]$coefficients[-1] * new_x1)
    boot_meanss[i] = as.numeric(new_fit1)
  }
  
}

result_matrix <- matrix(boot_meanss, nrow = n_bootstrap, ncol = length(new_x_points), byrow = TRUE)
answer <- matrix(NA, nrow = n_bootstrap, ncol = length(new_x_points))

for (j in 1:length(new_x_points)) {
  for (i in 1:n_bootstrap) {
    value <- result_matrix[i, j]
    in_interval <- value >= ci_result$lower[j] & value <= ci_result$upper[j]
    answer[i, j] <- in_interval
  }
}
print(answer)
count_result <- table(answer)
true_count <- count_result["TRUE"]
false_count <- count_result["FALSE"]
true_count
false_count
