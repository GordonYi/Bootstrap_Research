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

x <- seq(0, 30, by = 0.1)  
x
y_values <- y(x)
y_values

num_points <- length(x) 
num_points
num_noise <- 1  
noise <- rnorm(num_points, mean = 0, sd = 1)  
noise
y_values_with_noise <- y_values + noise

plot(x, y_values_with_noise, type = "p", xlab = "x", ylab = "y(x) + noise", main = "Simulated points for the known function", col = "blue", pch = 16)

count <- 0

for (i in 1:length(y_values_with_noise)) {
  count <- count + 1
}

cat(count, "\n")
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



new_x_points <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
n_bootstrap <- 1000
truth2 <- data.frame(x = new_x_points, y = NA, lower = NA, upper = NA)

for (j in seq_along(new_x_points)) {
  newx <- new_x_points[j]
  boot_means <- rep(NA, n_bootstrap)
  
  for (i in 1:n_bootstrap) {
    boot_indices <- sample(n, replace = TRUE)  
    boot_sample <- dataset[boot_indices, ]  
    mfitsss <- lapply(dfs, function(d) {
      with(boot_sample, mono_fit(x, y, df = 6))
    })
    
    new_x <- predict(mfitsss[[1]]$x_mat, newx = newx)
    new_fit <- mfitsss[[1]]$coefficients[1] + sum(mfitsss[[1]]$coefficients[-1] * new_x)
    boot_means[i] <- as.numeric(new_fit)
  }
  
  CI <- quantile(boot_means, probs = c(0.025, 0.975))
  truth2$lower[j] <- CI[1]
  truth2$upper[j] <- CI[2]
  truth2$y[j] <- new_fit
}

truth2



set.seed(1)  

x_points <- 1:10  
num_points <- length(x_points) 
num_noise <- 100  
noise_mean <- 0  
noise_sd <- 0.1  

y_values_with_noise <- numeric(num_points * num_noise)

for (i in 1:num_points) {
  x <- x_points[i]
  y <- 0.25 * (x - 0.9)
  noise <- rnorm(num_noise, mean = noise_mean, sd = noise_sd)
  y_with_noise <- y + noise
  start_index <- (i - 1) * num_noise + 1
  end_index <- i * num_noise
  y_values_with_noise[start_index:end_index] <- y_with_noise
}

y_values_with_noise



x_points <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

ci_data <- data.frame(x = truth2$x,
                      y = truth2$y,
                      lower = truth2$lower,
                      upper = truth2$upper)

ggplot() +
  geom_point(data = ci_data, aes(x = x, y = y), color = 'black', alpha = 0.6) +
  geom_line(data = data.frame(x = ci_data$x, y = ci_data$y), 
            aes(x = x, y = y), color = 'red', size = 1.2) +
  geom_errorbar(data = ci_data, aes(x = x, ymin = lower, ymax = upper),
                color = 'blue', width = 0.2) +
  labs(x = "x", y = "y") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))




x_points <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

y_points <- 0.25 * (x_points - 0.9)

set.seed(1)  
num_noise <- 100
noise <- rnorm(length(x_points) * num_noise, mean = 0, sd = 0.1)
y_values_with_noise <- rep(y_points, each = num_noise) + noise

data <- data.frame(x = rep(x_points, each = num_noise),
                   y_with_noise = y_values_with_noise)

ci_data <- data.frame(x = truth2$x,
                      y = truth2$y,
                      lower = truth2$lower,
                      upper = truth2$upper)

ggplot() +
  geom_point(data = data, aes(x = x, y = y_with_noise), color = 'blue', alpha = 0.6) +
  geom_line(data = data.frame(x = x_points, y = y_points), 
            aes(x = x, y = y), color = 'red', size = 1.2) +
  geom_line(data = ci_data, aes(x = x, y = y), color = 'green', size = 1.2) +
  geom_errorbar(data = ci_data, aes(x = x, ymin = lower, ymax = upper),
                color = 'green', width = 0.2) +
  labs(x = "x", y = "y", title = "y = 0.25(x)-0.9 with Confidence Intervals") +
  scale_color_manual(values = c("red", "green"),
                     labels = c("True Function", "Confidence Intervals")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))



group_size <- 100
num_groups <- length(x_points)
coverage <- rep(0, num_groups)

for (i in 1:num_groups) {
  start_index <- (i - 1) * group_size + 1
  end_index <- i * group_size
  y_values_group <- y_values_with_noise[start_index:end_index]
  coverage[i] <- mean(y_values_group >= ci_data$lower[i] & y_values_group <= ci_data$upper[i]) * 100
}

truth2$coverage <- coverage

truth2



# y = 0.02(x^3 +1.2)
y <- function(x) {
  return(0.02 * (x^3 - 1.2))
}

x <- seq(0, 30, by = 0.1)
y_values <- y(x)

num_points <- length(x)
num_noise <- 1
noise <- rnorm(num_points, mean = 0, sd = 1)
y_values_with_noise <- y_values + noise

plot(x, y_values_with_noise, type = "p", xlab = "x", ylab = "y(x) + noise", main = "Simulated points for the known function", col = "blue", pch = 16)

count <- length(y_values_with_noise)
cat(count, "\n")

dataset <- data.frame(x = x, y = y_values_with_noise)

dfs <- c(6, 10, 15)
mfitsss <- lapply(dfs, function(d) {
  with(dataset, mono_fit(x, y, df = d))
})


new_x_points <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
n_bootstrap <- 1000
truth3 <- data.frame(x = new_x_points, y = NA, lower = NA, upper = NA)

for (j in seq_along(new_x_points)) {
  newx <- new_x_points[j]
  boot_means <- rep(NA, n_bootstrap)
  for (i in 1:n_bootstrap) {
    boot_indices <- sample(n, replace = TRUE)  
    boot_sample <- dataset[boot_indices, ]  
    mfitsss <- lapply(dfs, function(d) {
      with(boot_sample, mono_fit(x, y, df = 6))
    })
    
    new_x <- predict(mfitsss[[1]]$x_mat, newx = newx)
    new_fit <- mfitsss[[1]]$coefficients[1] + sum(mfitsss[[1]]$coefficients[-1] * new_x)
    boot_means[i] <- as.numeric(new_fit)
  }
  
  CI <- quantile(boot_means, probs = c(0.025, 0.975))
  truth3$lower[j] <- CI[1]
  truth3$upper[j] <- CI[2]
  truth3$y[j] <- new_fit
}

truth3

x_points <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

y_points <- 0.02 * (x_points^3 +1.2)

set.seed(1)  
num_noise <- 100
noise <- rnorm(length(x_points) * num_noise, mean = 0, sd = 0.1)
y_values_with_noise <- rep(y_points, each = num_noise) + noise

data <- data.frame(x = rep(x_points, each = num_noise),
                   y_with_noise = y_values_with_noise)

ci_data <- data.frame(x = truth3$x,
                      y = truth3$y,
                      lower = truth3$lower,
                      upper = truth3$upper)

ggplot() +
  geom_point(data = data, aes(x = x, y = y_with_noise), color = 'blue', alpha = 0.6) +
  geom_line(data = data.frame(x = x_points, y = y_points), 
            aes(x = x, y = y), color = 'red', size = 1.2) +
  geom_line(data = ci_data, aes(x = x, y = y), color = 'green', size = 1.2) +
  geom_errorbar(data = ci_data, aes(x = x, ymin = lower, ymax = upper),
                color = 'green', width = 0.2) +
  labs(x = "x", y = "y", title = "y = 0.02(x+1.2) with Confidence Intervals") +
  scale_color_manual(values = c("red", "green"),
                     labels = c("True Function", "Confidence Intervals")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))



group_size <- 100
num_groups <- length(x_points)
coverage <- rep(0, num_groups)

for (i in 1:num_groups) {
  start_index <- (i - 1) * group_size + 1
  end_index <- i * group_size
  y_values_group <- y_values_with_noise[start_index:end_index]
  coverage[i] <- mean(y_values_group >= ci_data$lower[i] & y_values_group <= ci_data$upper[i]) * 100
}

truth3$coverage <- coverage
truth3






# y = 0.1 * (5 * pnorm((x - 2) / 0.3)) + 1)
y <- function(x) {
  return( 0.1 * (5 * pnorm((x - 2) / 0.3)) + 1)
}

x <- seq(0, 20, by = 0.05)
y_values <- y(x)

num_points <- length(x)
num_noise <- 100
noise <- rnorm(num_points, mean = 0, sd = 0.1)
y_values_with_noise <- y_values + noise

plot(x, y_values_with_noise, type = "p", xlab = "x", ylab = "y(x) + noise", main = "Simulated points for the known function", col = "blue", pch = 16)

count <- length(y_values_with_noise)
cat(count, "\n")

dataset3 <- data.frame(x = x, y = y_values_with_noise)

dfs <- c(6, 10, 15)  
mfitssss <- lapply(dfs, function(d) {
  with(dataset3, mono_fit(x, y, df = d))
})

ggplot() +
  geom_point(data = dataset3, aes(x = x, y = y_values_with_noise), color = 'blue', alpha = 0.6) +
  geom_line(data = data.frame(x = dataset3$x, y = mfitssss[[1]]$fitted), 
            aes(x = x, y = y), color = 'red', size = 1.2) +
  labs(x = "x", y = "y") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  annotate("text", x = max(dataset3$x), y = max(dataset3$y),
           label = paste("df =", dfs[1]), hjust = 1, vjust = 1, color = "black", size = 5)

new_x_points <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
n_bootstrap <- 1000
truth4 <- data.frame(x = new_x_points, y = NA, lower = NA, upper = NA)

for (j in seq_along(new_x_points)) {
  newx <- new_x_points[j]
  boot_means <- rep(NA, n_bootstrap)
  
  for (i in 1:n_bootstrap) {
    boot_indices <- sample(num_points, replace = TRUE)  
    boot_sample <- dataset3[boot_indices, ]  
    mfitssss <- lapply(dfs, function(d) {
      with(boot_sample, mono_fit(x, y, df = d))
    })
    
    new_x <- predict(mfitssss[[1]]$x_mat, newx = newx)
    new_fit <- mfitssss[[1]]$coefficients[1] + sum(mfitssss[[1]]$coefficients[-1] * new_x)
    boot_means[i] <- as.numeric(new_fit)
  }
  
  CI <- quantile(boot_means, probs = c(0.025, 0.975))
  truth4$lower[j] <- CI[1]
  truth4$upper[j] <- CI[2]
  truth4$y[j] <- new_fit
}

truth4

x_points <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

y_points <- 0.1 * (5 * pnorm((x_points - 2) / 0.3)) + 1

set.seed(1)  
num_noise <- 100
noise <- rnorm(length(x_points) * num_noise, mean = 0, sd = 0.1)
y_values_with_noise <- rep(y_points, each = num_noise) + noise

data <- data.frame(x = rep(x_points, each = num_noise),
                   y_with_noise = y_values_with_noise)

ci_data <- data.frame(x = truth4$x,
                      y = truth4$y,
                      lower = truth4$lower,
                      upper = truth4$upper)

ggplot() +
  geom_point(data = data, aes(x = x, y = y_with_noise), color = 'blue', alpha = 0.6) +
  geom_line(data = data.frame(x = x_points, y = y_points), 
            aes(x = x, y = y), color = 'red', size = 1.2) +
  geom_line(data = ci_data, aes(x = x, y = y), color = 'green', size = 1.2) +
  geom_errorbar(data = ci_data, aes(x = x, ymin = lower, ymax = upper),
                color = 'green', width = 0.2) +
  labs(x = "x", y = "y", title = "y = 0.1(5Φ(x − 2)/0.3) + 1) with Confidence Intervals") +
  scale_color_manual(values = c("red", "green"),
                     labels = c("True Function", "Confidence Intervals")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

group_size <- 100
num_groups <- length(x_points)
coverage <- rep(0, num_groups)

for (i in 1:num_groups) {
  start_index <- (i - 1) * group_size + 1
  end_index <- i * group_size
  y_values_group <- y_values_with_noise[start_index:end_index]
  coverage[i] <- mean(y_values_group >= ci_data$lower[i] & y_values_group <= ci_data$upper[i]) * 100
}

truth4$coverage <- coverage
truth4

# We can see that in some situations the 1000 times bootstrap CI had poor performance in certain points -- however, if you increase the bootstrap sampling from 1000 to 2000
# the behavior will be much better

# It seems bootstrap CI have a poor performance in function 3 -- still finding out the reason for that
