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
dfs <- c(6)
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





# y = 0.25*(x-0.9)


y <- function(x) {
  return(0.25*(x-0.9))
}
x <- seq(0, 30, length.out = 1000)
y_values <- y(x)

noise <- rnorm(length(x), mean = 0, sd = 1)  
y_values_with_noise <- y_values + noise

plot(x, y_values_with_noise, type = "p", xlab = "x", ylab = "y(x) + noise", main = "Simulated points for the known function", col = "blue", pch = 16)

dataset <- data.frame(x = x, y = y_values_with_noise)
head(dataset)

dfs <- c(6, 10, 15)
mfitsss <- lapply(dfs, function(d) {
  with(dataset, mono_fit(x, y, df = d))
})


ggplot() +
  geom_point(data = dataset, aes(x = x, y = y), color = 'pink', alpha = 0.6) +
  geom_line(data = data.frame(x = dataset$x, y = mfitsss[[1]]$fitted), 
            aes(x = x, y = y), color = 'red', size = 1.2) +
  labs(x = "x", y = "y") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  geom_abline(slope = 0.25, intercept = -0.25*0.9, linetype = "dashed",
              aes(color = "y(x) = 0.25 * (x - 0.9)")) +
  annotate("text", x = max(dataset$x), y = max(dataset$y),
           label = "df = 6", hjust = 1, vjust = 1, color = "black", size = 5)



new_x_points <- seq(1,30,length = 10)
n_bootstrap <- 1000
iter = 100
n = nrow(dataset)
result = list()

for (j in seq_along(new_x_points)) {
  newx <- new_x_points[j]
  true_y = y(newx)
  boot_means <- rep(NA, n_bootstrap)
  
  ci = list()
  
  for(k in 1:iter){
    
    fit_value_list = list()
    
    for (i in 1:n_bootstrap) {
      boot_indices <- sample(n, replace = TRUE)
      boot_sample <- dataset[boot_indices, ]
      
      mfitsss <- lapply(dfs, function(d) {
        with(boot_sample, mono_fit(x, y, df = d))
      })
      
      fit_list = list()
      
      for(m in 1:length(mfitsss)){
        
        new_x <- predict(mfitsss[[m]]$x_mat, newx = newx)
        new_fit <- mfitsss[[m]]$coefficients[1] + sum(mfitsss[[m]]$coefficients[-1] * new_x)
        
        fit_list[[m]] = new_fit
        
      }
      
      y_fit = do.call('c',fit_list)
      fit_value_list[[i]] = y_fit
      
    }
    
    fit_value_list = do.call('rbind',fit_value_list)
    colnames(fit_value_list) = str_c('df_',dfs)
    ci[[k]] = apply(X = fit_value_list,MARGIN = 2,FUN = quantile,probs = c(0.025, 0.975))
    
  }
  
  coverage_list = list()
  
  for(i in 1:length(ci)){
    
    temp_coverage = rep(NA,ncol(ci[[i]]))
    
    for(k in 1:ncol(ci[[i]])){
      
      temp_ci = ci[[i]][,k]
      test1 = temp_ci[1] <= true_y
      test2 = temp_ci[2] >= true_y
      temp_coverage[k] = all(test1,test2)
      
    }
    
    coverage_list[[i]] = temp_coverage
    
  }
  
  coverage_list = do.call('rbind',coverage_list)
  colnames(coverage_list) = str_c('df_',dfs)
  coverage = colMeans(coverage_list) %>% t() %>% as_tibble()
  result[[j]] = tibble(x = newx,
                       y = true_y) %>% 
    bind_cols(coverage)
}

result = do.call('rbind',result)
result


#################################################################################
#################################################################################
#################################################################################
cilrep = function(tis, fun, df = 6, conf.level = .95) {

  ## generate data: x and y given f(t), noise variance
  x = seq(0, 30, length.out = 1000)
  y_values = fun(x)
  noise = rnorm(length(x), mean = 0, sd = 1)  
  y_values_with_noise = y_values + noise
  dataset = data.frame(x = x, y = y_values_with_noise)
  
  ## fit f(x) given (x, y) and estimate theta (point estimate)
  mfitsss = with(dataset, mono_fit(x, y, df = df))
  theta = mfitsss$coefficients
  
  ## bootstrap 1000 times, and obtain theta b, b = 1, \ldots, B
  B = 1000
  x_mat = list()
  boot_theta = list()
  
  for(i in 1:B){
    id_boot = sample(x = 1:nrow(dataset),size = nrow(dataset),replace = T)
    boot_data = dataset[id_boot,]
    boot_mfitss = with(boot_data, mono_fit(x, y, df = df))
    x_mat[[i]] = boot_mfitss$x_mat
    boot_theta[[i]] = boot_mfitsss$coefficients
  }
  
  ## construct CI for each f(t), t in tis (length lt)## return ciFts # alt x 2 matrix
  ## vector of TRUE/FALSE of length lt
  coverage = rep(NA,length(tis))
  
  for(i in seq_along(tis)){
    
    fit_value = rep(NA,B)
    
    for(j in 1:B){
      new_t = predict(x_mat[[j]],newx = tis[i])
      fit_value[j] = sum(boot_theta[[j]] * c(1,new_t))
    }
    
    ci_fit = quantile(fit_value,probs = c((1-conf.level)/2,1-(1-conf.level)/2))
    true_fit = fun(tis[i])
    coverage[i] = true_fit>=ci_fit[1] & true_fit<=ci_fit[2]
    
  }
  
  return(coverage)
  
}


tis = seq(1,30,length = 10)
fun = function(x) {
  return(0.25*(x-0.9))
}
nrep = 100
sim = replicate(nrep, cilrep(tis,fun,df = 6))
colMeans(sim)






