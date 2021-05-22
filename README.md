
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

## Example 1: multiple testing: one sample t test

This is a basic example which shows you how to solve a simple problem:

``` r
library(sdafilter)
n = 50
p = 100

rho = 0.8
Sig = matrix(rho, p, p)
diag(Sig) = 1

dat <- MASS::mvrnorm(n, rep(0, p), Sig)
mu = rep(0, p)
mu[1:as.integer(0.1*p)]=0.5
dat = dat+rep(1, n)%*%t(mu)

alpha = 0.2
out = SDA_M(dat, alpha)
#> Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 9%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 19%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 30%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 40%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 50%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 60%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 70%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 80%Conducting the graphical lasso (glasso) wtih lossless screening....in progress: 90%
#> Conducting the graphical lasso (glasso)....done.                                          
#> Conducting Subsampling....in progress:5% Conducting Subsampling....in progress:10% Conducting Subsampling....in progress:15% Conducting Subsampling....in progress:20% Conducting Subsampling....in progress:25% Conducting Subsampling....in progress:30% Conducting Subsampling....in progress:35% Conducting Subsampling....in progress:40% Conducting Subsampling....in progress:45% Conducting Subsampling....in progress:50% Conducting Subsampling....in progress:55% Conducting Subsampling....in progress:60% Conducting Subsampling....in progress:65% Conducting Subsampling....in progress:70% Conducting Subsampling....in progress:75% Conducting Subsampling....in progress:80% Conducting Subsampling....in progress:85% Conducting Subsampling....in progress:90% Conducting Subsampling....in progress:95% Conducting Subsampling....in progress:100% Conducting Subsampling....done.                  
print(out)
#>  [1]  1  2  3  4  5  6  7  8  9 10 55 93
```

## Example 2: multiple testing: two sample t test

``` r
p = 100
n = 30

dat_I = matrix(rnorm(n*p),nrow = n)
mu = rep(0, p)
mu[1:10] = 1.5
dat_I = dat_I = rep(1, n)%*%t(mu)

dat_II = matrix(rnorm(n*p), nrow = n)

Sigma_I = diag(p)
Sigma_II = diag(p)

out = SDA_2S(dat_I, dat_II, alpha=0.05, Sigma_I, Sigma_II)
print(out)
#>  [1]  4  3  5  7  2  1  6  8  9 10
```

## Example 3: new function for one-sample t test, with various test statistics for the first sample provided

``` r
n = 50
p = 100

rho = 0.8
Sig = matrix(rho, p, p)
diag(Sig) = 1

dat <- MASS::mvrnorm(n, rep(0, p), Sig)
mu = rep(0, p)
mu[1:as.integer(0.1*p)]=0.5
dat = dat+rep(1, n)%*%t(mu)
alpha = 0.2
# kwd = {'lasso', 'de_lasso', 'innovate', 'pfa'}
out = SDA_robust(dat, alpha, solve(Sig), kwd='innovate')
print(out)
#>  [1]  1  2  3  4  5  6  7  8  9 10 78
```
