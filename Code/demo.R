
## for each given dataset and sequence of times, return CI for
## f(t) at the given time points

ci1rep <- function(tis, df = 6, conf.level = .95) {
    ## generate data: x and y given f(t), noise variance
    
    ## fit f(x) given (x, y) and estimate theta (point estimate)

    ## bootstrap 1000 times, and obtain theta_b, b = 1, \ldots, B
    ## Put them in a B by df matrix
 
    ## construct CI for each f(t), t in tis (length lt)
    ## return ciFts # alt x 2 matrix

    ## vector of TRUE/FALSE of length lt
    
}


nrep <- 1000
sim <- replicate(nrep, ci1rep(tis))


