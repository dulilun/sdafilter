
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sdafilter

<!-- badges: start -->

<!-- badges: end -->

The goal of sdafilter is to provide an R package for our new paper on
arXiv; see the link: <https://arxiv.org/abs/2002.11992>

## Installation

You can install the released version of sdafilter from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("sdafilter")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dulilun/sdafilter")
```

## Example 1: one sample t test

This is a basic example which shows you how to solve a simple problem:

``` r
library(sdafilter)
n = 50
p = 100
dat = matrix(rnorm(n*p), nrow=n)
mu = rep(0, p)
mu[1:as.integer(0.1*p)]=0.3
dat = dat+rep(1, n)%*%t(mu)
alpha = 0.2
out = SDA_M(dat, alpha, diag(p))
#> Warning in min(t[which(Ta <= alpha)]): no non-missing arguments to min;
#> returning Inf

#> Warning in min(t[which(Ta <= alpha)]): no non-missing arguments to min;
#> returning Inf

#> Warning in min(t[which(Ta <= alpha)]): no non-missing arguments to min;
#> returning Inf

#> Warning in min(t[which(Ta <= alpha)]): no non-missing arguments to min;
#> returning Inf

#> Warning in min(t[which(Ta <= alpha)]): no non-missing arguments to min;
#> returning Inf
print(out)
#> [1] "no rejection"

p = 100
n = 30
dat_I = matrix(rnorm(n*p),nrow = n)
mu = rep(0, p)
mu[1:10] = 1.5
dat_I = dat_I = rep(1, n)%*%t(mu)
dat_II = matrix(rnorm(n*p), nrow = n)
out = SDA_2S(dat_I, dat_II, alpha=0.05)
#> Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 9%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 19%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 30%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 40%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 50%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 60%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 70%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 80%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 90%
#> Conducting the graphical lasso (glasso)....done.                                          
#> Conducting extended Bayesian information criterion (ebic) selection....done
print(out)
#>  [1]  4  6  1  5  7  9  2 10  8  3 81
```

## Example 2: two sample t test

``` r
p = 100
n = 30
dat_I = matrix(rnorm(n*p),nrow = n)
mu = rep(0, p)
mu[1:10] = 1.5
dat_I = dat_I = rep(1, n)%*%t(mu)
dat_II = matrix(rnorm(n*p), nrow = n)
out = SDA_2S(dat_I, dat_II, alpha=0.05)
#> Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 9%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 19%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 30%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 40%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 50%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 60%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 70%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 80%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 90%
#> Conducting the graphical lasso (glasso)....done.                                          
#> Conducting extended Bayesian information criterion (ebic) selection....done
print(out)
#>  [1]  5  8  9  2 10  6  7  3  1  4 23 85
```
