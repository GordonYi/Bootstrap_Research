library(splines2)
library(quadprog)
library(boot)
library(SemiPar)
library(plotrix)
library(tidyverse)
library(ggthemes)

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
                -as.numeric(crossprod(x_bar * y_sd / x_sd, sbeta)) + y_bar)
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
dfs <- c(6, 10, 14)
mfits <- lapply(dfs, function(d) {
  with(age.income, mono_fit(age, log.income, df = d))
})
curv_fit <- function(x, y, monotone = c("none", "increasing", "decreasing"),
                     curvature = c("convex", "concave"), ...) {
  monotone <- match.arg(monotone); curvature <- match.arg(curvature)
  min_x <- min(x); range_x <- max(x) - min_x
  x_mat <- cSpline(x, ..., intercept = TRUE) 
  if (monotone != "none") {
    right_side <- (monotone == "decreasing" && curvature == "convex") ||
      (monotone == "increasing" && curvature == "concave")
    side_knot <- knots(x_mat, "boundary")[right_side + 1]
    dx_mat <- cbind(1 / range_x, deriv(predict(x_mat, side_knot)))
  }
  x_mat <- cbind("(Linear)" = (x - min_x) / range_x, x_mat)
  standardized <- standardize(x_mat, y)
  A_mat <- rbind(0, diag(ncol(x_mat) - 1L))
  if (curvature == "concave") A_mat <- - A_mat
  if (monotone == "decreasing") dx_mat <- - dx_mat
  if (monotone != "none") {
    for (j in seq_len(ncol(dx_mat)))
      dx_mat[, j] <- dx_mat[, j] / standardized$x_sd[j]
    A_mat <- cbind(A_mat, as.numeric(dx_mat))
  }
  qp_solve(standardized, A_mat, x_mat)
}
cfits <- lapply(dfs, function(d) {
  with(age.income, curv_fit(age, log.income, "none", "concave", df = d))
})

ggplot() +
  geom_point(data = age.income,
             mapping = aes(x = age,y = log.income)) +
  geom_point(data = data.frame(age = age.income$age,
                               log.income = mfits[[1]]$fitted),
             mapping = aes(x = age,y = log.income),
             color = 'red') +
  theme_bw()

dfs <- c(6)
mfits <- lapply(dfs, function(d) {
  with(age.income, mono_fit(age, log.income, df = d))
})

new_x_points = seq(range(age.income$age)[1],range(age.income$age)[2],length = 10)
fit.log.income = rep(NA,length(new_x_points))
for (h in seq_along(new_x_points)){
  
  newx = new_x_points[h]
  new_x = predict(mfits[[1]]$x_mat, newx = newx)
  new_fit = mfits[[1]]$coefficients[1] + sum(mfits[[1]]$coefficients[-1] * new_x)
  fit.log.income[h] = new_fit
  
}

truth = data.frame(age = new_x_points,
                   fit.log.income = fit.log.income,
                   lower = NA,
                   upper = NA)

n <- nrow(age.income)
n_bootstrap <- 1000

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
  truth$lower[j] = CI[1]
  truth$upper[j] = CI[2]
  
}

truth$cover = truth$fit.log.income>truth$lower & truth$fit.log.income<truth$upper
truth

library(tidyverse)
ggplot() +
  geom_point(data = age.income,
             mapping = aes(x = age,y = log.income)) +
  geom_point(data = data.frame(age = age.income$age,
                               log.income = mfits[[1]]$fitted),
             mapping = aes(x = age,y = log.income),
             color = 'red') +
  geom_point(data = truth,
             mapping = aes(x = age,y = fit.log.income),
             color = 'blue') +
  geom_errorbar(data = truth,
                mapping = aes(x = age,ymin = lower,ymax = upper),
                width = 1.5) +
  theme_bw()

confidence_interval <- data.frame(age = new_x_points,
                                  lower = truth$lower,
                                  upper = truth$upper)


ggplot() +
  geom_point(data = age.income, aes(x = age, y = log.income), color = 'black', alpha = 0.6) +
  geom_line(data = data.frame(age = age.income$age, log.income = mfits[[1]]$fitted), 
            aes(x = age, y = log.income), color = 'red', size = 1.2) +
  geom_smooth(data = age.income, aes(x = age, y = log.income), method = "loess",
              color = 'orange', size = 1.5) +
  geom_point(data = truth,
             mapping = aes(x = age,y = fit.log.income),
             color = 'purple') +
  geom_errorbar(data = truth,
                mapping = aes(x = age,ymin = lower,ymax = upper),
                width = 1.5) +
  labs(x = "Age", y = "log (Income)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
  annotate("text", x = max(age.income$age), y = max(age.income$log.income),
         label = "df = 6", hjust = 1, vjust = 1, color = "black", size = 12)