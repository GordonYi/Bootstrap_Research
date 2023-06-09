# Bootstrap_Research

This is the repository for my first research project!! About growth curve model and bootstrap!

This is the meeting minute notes for the research, and it is writting in flashback, so the latest note is on the top.

This is the RMarkdown version of the meeting note, if you want to see it in .txt file, please see MeetingNote.txt in the repo

## Research Note:

### Report for meeting 6

7/17

#### Progress:

1. Did the I-Spline regression for the age.income data

2. Did the bootstrap sampling and the 95% CI for the dataset

#### What are we doing:

I dont know what the truth it

In a simulation study, I will generate the data and my true values will be known:

1. define a monotone function -- smooth, monotone function , use this function as a true μ(t) -- the expectation of y(t), use this to generate the dataset 1000 times, and then
add some noise

y = ax + bx^2 + c

generate some t whatever it value is

use normal distribution for the noise

Repeat the process 1000 times, the percentage of the truth value

y = μ(t) + noise

I need to know t first, and then i know μ(t),

and I generate y by setting the noise to normal(some variance of normal, variance is unchange),

my goal is to estimate μ(t)

put my bootstrap inside this setting to check is that real a 95% CI

suppose I have 5 time points, at every time point we want to check what is the empirical coverage of those CI,

at least 95% time inside 1000 times are inside the CI

look at what CI really is


### Report for meeting 5

7/3

#### Progress:

1. Rewrite the code to correct the error

2. Clean the Repo!

#### Question：

I think I did those things right and I think the next step is to do the simulation study, and I am thinking about that

#### What are we doing:

Fit the Splines regression for the data in Semipar

Do the bootstrap sampling for n = 1000 times, find out the 95% CI for 10 time points inside the regression

Things need to solve: do a simulation study, find out wheather those simulated points are dropped inside the confidence intervals

Do it 1000 times, find out how many of those points are inside the CI

### Report for meeting 4

6/22

#### Progress:

Fit the Splines regression for the data in Semipar

Do the bootstrap sampling for n = 1000 times, find out the 95% CI for 10 time points inside the regression

Things need to solve: do a simulation study, find out wheather those simulated points are dropped inside the confidence intervals

Do it 1000 times, find out how many of those points are inside the CI

#### Questions:

##### 1. About regression fitting: I think the code I found from the resource is somehow not appropriate, but I don't know how to fix that:

I use the code from the paper you gave me to do the splines regression fitting:(This part of the code is copied from the paper)

```{r,include = FASLE}
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
```
If we draw a line of this regression we can see that the regression is above the top of the dataset:

```{r,include = FASLE}
ggplot() +
  geom_point(data = age.income,
             mapping = aes(x = age,y = log.income)) +
  geom_point(data = data.frame(age = age.income$age,
                               log.income = mfits[[1]]$fitted),
             mapping = aes(x = age,y = log.income),
             color = 'red') +
  theme_bw()
```
We can see that the fitting is not correct.

If we change the dfs from 6 to 4, it would be better, but still not very good. If the df goes higher, the line will be much more higher.

In my code: 

line 47:

```{r,include = FASLE}
new_x = predict(mfits[[1]]$x_mat, newx = newx)
```
This line wants to convert x to x_mat, for fit the splines regression

line 48:

```{r,include = FASLE}
new_fit = mfits[[1]]$coefficients[1] + sum(mfits[[1]]$coefficients[-1] * new_x)
```

This line wants to fit the regression by hand because it seems we can not use the predict function directly -- the function is stuck inside a chunck of codes

the manually fitted result is the same as the fitted result generated by the code line 7 to line 40,

to check that, we can do:

```{r,include = FASLE}
mfits[[1]]$fitted[1] = 15.60871
fit.log.income[1] = 15.60871
```

Which is the same.

We can plot the function generated by the code out:

```{r,include = FASLE}
ggplot() +
  geom_point(data = age.income,
             mapping = aes(x = age,y = log.income)) +
  geom_point(data = data.frame(age = age.income$age,
                               log.income = mfits[[1]]$fitted),
             mapping = aes(x = age,y = log.income),
             color = 'red') +
  theme_bw()
```

While we can see that the fitted regression is not inside the data points, so I think the code inside the paper is wrong

##### 2. How to do simulation study

I am sorry about that, but I can't find a good resource online to do a simulation study,

I don't understand what is "give out a signal to the original data to do the simulation study"

Please give me a resource to do that

祝老师端午节安康！

### Report for meeting 3

6/5

#### Progress:
1.Use dataset age.income from package SemiPar,

2.fit the ispline monotone regression for the data, get x_mat, coefficient, fitted value,

3.Resampling by bootstrap method n_bootstrap = 1000 times, for each time, fit the ispline monotone regression for the bootstrap data, get x_mat, coefficient, fitted value, 

4.Find out the 10th, 20th, 30th... 90th percentiles of the dataset, fit f(10th), f(20th)... for n_bootstrap = 1000 times

5.Calculate the Confidence Interval for the fitted value

6.Find the fitted value generated by the regression by the real data

7.Check wheather the data is fitted inside the CI or not

       age         income    lower    upper   cover   SE generated by 1000 bootstrap resampling
1  21.00000       15.34254 13.84109 15.73719  TRUE    0.6699820
2  25.88889       17.18788 15.08217 17.95271  TRUE    0.8516536
3  30.77778       17.33796 15.22696 18.18516  TRUE    0.8718291
4  35.66667       17.34434 15.29922 18.14767  TRUE    0.8359431
5  40.55556       17.34571 15.30806 18.17145  TRUE    0.8719740
6  45.44444       17.34579 15.33284 18.27347  TRUE    0.8557299
7  50.33333       17.34579 15.31484 18.14822  TRUE    0.8721111
8  55.22222       17.34579 15.28231 18.15750  TRUE    0.8618083
9  60.11111       17.34579 15.27811 18.26073  TRUE    0.8716605
10 65.00000       17.34579 15.31106 18.21312  TRUE    0.8671026

Trying to estimate standard error in those 10 points

Should I repeat this process for 1000 times, to see wheather in >95% situations those CIs are covering the true value?

I tried to repeat it 1000 times, everytime the 95% CI (generated by the bootstrap) will cover the true value generated by the true data function

Is it prove that bootstrap is a good way to do estimation in this situation?

cover_true_count
1000 1000 1000 1000 1000 1000 1000 1000 1000 1000
cover_false_count
0 0 0 0 0 0 0 0 0 0
                           
In 1000 simulation studies, every bootstrap CIs covers the true value generated by the iSpline function

Which proves bootstrap is a good way to do estimation in this iSpline situation

footnote: 

age is generated by the 10th, 20th, 30th... percentile of the dataset

income is estimated by the iSpline function from the paper, which is a true value that generated by the iSpline regression

Lower and Upper boundaries are the 95% CI created by 1000 times bootstrap fitted regressions

All the true points generated by the original iSpline function are covered inside the 95% CI

SE are generated by the 1000 bootstrap resamplings

#### Question:

1. I estimate the CI by bootstrapping the data for 1000 times, and I can see that in every f(t) I choose, the CI of the bootstrap covers the real value, 

which might indicate that bootstrap is a good method to predict f(t), and these intervals are truely 95% intervals, but this is not about Standard Error..? 

Or can I say that since every true point falls within its corresponding confidence interval (CI), you can tentatively conclude that the Bootstrap standard error (SE) is reliable?

2. 

“do a simulation study, need to have a model, if the fitted model is the true model, generated r = 1000 dataset from the distribution”

I think I might do this wrong, let me think about that

--- I tried repeat the thing I did for 1000 times, all the CIs are covering the true value

Some interesting findings:

CIs are longer on the bottom side and shorter for the upper side

### Report for meeting 2

5/22

The thinking is wrong

I need to do a bootstrap sampling, and I get 205 values,

and I fit a function by the values, 

and I get the coefficients,

and I use this coefficients to get f(20),f(30).. I need 10 values,

and I need to use those values get CIs

#### Questions:
1. What is quadratic programming? (not that important)
2. Did the bootstrap out and make those CIs

### Report for meeting 1

5/8

#### Progress:
Some functions must be shape-restricted because of natural intuitions and prior knowledge. Some shape-restricted functions examples: Monotone Regressions, Convex Regressions... 

e.g. growth curves in population health, dose-response models in phase I clinical trials, and utility functions in economics

mainly about non-decreasing curve

A non-decreasing curve can be fitted by a linear combination of I-splines with an additional intercept where the coefficients of the I-spline basis functions are constrained to be nonnegative

We want to do least square to fit the I-splines regression, which contains several intervals, which means several betas, since we want the regression to be monotonically increasing,

we want the beta in every intervals be non-negative (I'm not very sure about this, but I guess this is the reason)

And while we fitted the regression, we want to test that can bootstrap do a good job in estimating Standard Error with the shape-restricted regression:

So first we extract the fitted value from our regression, and do a bootstrap,
than we use the formula of Standard Error to see the real standard error, and compare them.
According to the dataset given by the paper, I did a I-splines regression(with df = 6), and apply bootstrap with it, find out the value of SE it given is very close to the SE value given by the formula.

#### Questions:
1. Due to time limitation, I'm still not very farmiliar with the mathematical details behind those methods, which will be a hinder of the paper writting, will learn more by reading more
2. I am not very sure about what is the SE I need to deal with -- I did a bootstrap for the fitted value given by the non-negative I-spline function, but not sure while is it the right thing to do
3. And shape-restricted regression should be
4. **quardratic programming
What is that?
I do need to know how this works, MLE will not apply when the estimator is on the boundary**
#### Things to do: 
We want to construct CI for the f(t) over all t 

did the bootstrap on the original data,

get the bootstrap copy of the data,

coefficients are on the boundary so we can't get SE directly,

take 10 values equally spaced from the lowest to the highest,

(How to get that?)

get f^(t) by the function,

we have 10 time points

suppose the range is from 0 to 10 -- i look at the time estimation from 0 to 10

repeat this process B = 100 times,

100 by 10 matrix,

got 1000 values for each quantity,

do a 95% CI for them, and I got 1 interval,

is it really a 95% interval?

do a simulation study, need to have a model, if the fitted model is the true model, generated r = 1000 dataset from the distribution, 

See whether 95% of the intervals cover the truth value

To see whether in what situation this interval make sense and some situations are not.

In this setting, you will check how many of the bootstrap CI really covers the real value.

Set this up in the repo, every meeting record will be inside it.

### Main method to research

bootstrap -- resampling approach -- estimate the uncertainty

1. find estimated theta
2. use boostrap to estimate se of theta
3. calculate the expectation of estimate of se of theta

1. how to use bootstrap to find uncertainty
2. learn how to use splines 2 in R via reading the paper

### Things to do next time:
1. set up a repo via github √
https://github.com/GordonYi/Bootstrap_Research

2. Know how to set up a simulation via bootstrap √
se_boot <- rep(NA, 1000)
for (i in 1:1000) {
      boot_sample <- sample(fitted_value, replace = TRUE)
      fitted_boot <- mean(boot_sample)
      se_boot[i] <- sd(boot_sample) / sqrt(length(boot_sample))
      }
 mean_se_boot <- mean(se_boot)
 
 Bootstrap_SE <- mean_se_boot
 

3.  Write down the outline of the research √
