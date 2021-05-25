#' Symmetrized Data Aggregation for one sample t-test; other options for calculating the test
#' statistics from the first sample
#'
#' Other commomly used test statistics for the first sample are allowed in this function.
#'
#'
#' @param dat a n by p data matrix
#' @param alpha the FDR level
#' @param Omega the inverse covariance matrix; if it is missing, it will be estimated
#' by the glasso package.
#' @param nonsparse if it is TRUE, the covariance matrix will be estimated by the POET
#' package; otherwise it will be fitted by glasso by default.
#' @param stable if it is TRUE, the sample will be randomly splitted \eqn{B=10} times for stability
#' performance; otherwise, only single sample splitting is used.
#' @param kwd various methods for calculating the test statistics from the first sample
#' @param scale if it is TRUE, the test statistic from the first sample will be standardized.
#'
#' @details We provide other commonly used test statistics for the first sample. These include
#' the debiased lasso, innovated transformation, and the factor-adjusted test statistics.
#'
#' @return the indices of the hypotheses rejected
#' @export
#' @examples
#' n = 50
#' p = 100
#' rho = 0.8
#' Sig = matrix(rho, p, p)
#' diag(Sig) = 1
#' dat <- MASS::mvrnorm(n, rep(0, p), Sig)
#' mu = rep(0, p)
#' mu[1:as.integer(0.1*p)]=0.5
#' dat = dat+rep(1, n)%*%t(mu)
#' alpha = 0.2
#' out = SDA_M(dat, alpha, solve(Sig), kwd='innovate')
#' print(out)
#'
SDA_M <- function(dat,
                  alpha,
                  Omega,
                  nonsparse = FALSE,
                  stable=TRUE,
                  kwd=c('lasso', 'de_lasso', 'innovate', 'pfa'),
                  scale=TRUE){
  # The core functions for the SDA method
  # The inverse of Sigma
  INV <- function(Sigma){
    EgSig0 = eigen(Sigma)
    EgVec = EgSig0$vectors
    Omega = EgVec%*%diag( 1/EgSig0$values )%*%t(EgVec)
    return( Omega )
  }

  # Omega^{1/2}
  Sqrt <- function(Omega_1){
    EgSig0 = eigen(Omega_1)
    EgVec = EgSig0$vectors
    Gamma = EgVec%*%diag( sqrt(EgSig0$values) )%*%t(EgVec)
    return( Gamma )
  }

  # Glasso for estimation the inverse covariance matrix
  Omega_est <- function(dat, lambda){

    p = ncol(dat)
    n = nrow(dat)
    Q = diag(n)-1/n*rep(1, n)%*%t(rep(1, n))
    Sigma_est = t(dat)%*%Q%*%dat/(n-1)
    # glasso packages
    GG = glasso::glasso(Sigma_est, lambda)
    #Omega_1 = GG$wi
    Sigma_1 = GG$w
    Sigma_1 = GG$w-diag(p)*lambda # bias correction on the diagonal
    Omega_1 = INV(Sigma_1)

    return(Omega_1)
  }

  # POET for Sigma/Omega
  POET_sig <- function(dat, K_hat){
    n = nrow(dat)
    Q = diag(n)-1/n*rep(1, n)%*%t(rep(1, n))
    Y = t(dat)%*%Q
    PP = POET::POET(Y, K_hat, 0.5, thres='soft', matrix = 'cor')
    Sigma_1 = PP$SigmaY
    Omega_1 = INV(Sigma_1)
    return(Omega_1)
  }


  # LASSO for screening
  Model_Select <- function(Y, X){
    # Method I: AIC
    fit <- glmnet::glmnet(X, Y, family = "gaussian")
    k <- fit$df
    AIC <- stats::deviance(fit)+2*k
    i_min <- which.min(AIC)
    lambda_select <- fit$lambda[i_min]
    fit_AIC = glmnet::glmnet(X, Y, family = "gaussian", lambda = lambda_select)
    w1 = fit_AIC$beta[,1]
    #w1 <- coef(fit, s = lambda_select, exact=TRUE)[-1]
    sv = as.vector( which(w1!=0) )

    # method II: upper bound p/3
    k= ceiling(p/3)
    if( length(sv)>k ){
      wv = which(fit$df==max(fit$df[fit$df<k]))[1]
      w1 = fit$beta[, wv]
    }
    return(w1)
  }


  # AIC algorithm by ourselves
  AIC_sda <- function(Y, X, w1){
    p <- length(w1)
    cut <- sort( abs(w1) )
    aic_v <- sapply(1:(p-1), function(i){
      sv = (1:p)[ abs(w1) > cut[i] ]
      fit <- stats::lm(Y~X[, sv]-1)
      aic <- stats::deviance(fit)+2*length(sv)
      return(aic)
    })

    i_s <- which.min(aic_v)
    sv = (1:p)[ abs(w1) > cut[i_s] ]

    k= ceiling(p/3)
    if( length(sv)>k ){
      sv = (1:p)[ abs(w1) > cut[2*k] ]
    }

    return(sv)
  }


  # generate W_j: one splitting
  W_value <- function(Gamma, dat, index, kwd, scale){
    n=nrow(dat)
    p=ncol(dat)
    n_1 = length(index)
    n_2 = n-n_1
    dat1 = apply(dat[index,], 2, mean)
    dat2 = apply(dat[-index,],2, mean)


    X= sqrt(n_1)*Gamma
    Xt = sqrt(n_2)*Gamma
    Y = as.vector(sqrt(n_1)*dat1%*%Gamma)
    Yt = as.vector(sqrt(n_2)*dat2%*%Gamma)


    # various options for T1
    if(kwd == 'lasso'){
      w1 = Model_Select(Y, X)
      sv = which( w1!=0 )

    }else if(kwd == 'innovate'){
      w1 = X%*%Y/sqrt(n_1)
      sv = AIC_sda(Y, X, w1)

    }else if(kwd == 'pfa'){
      Sigma <- INV(Gamma%*%Gamma)
      Test <- as.vector( sqrt(n_1)*dat1 )
      PFA <- pfa::pfa.test(Test, Sigma=Sigma, reg="L2", plot="none")
      aP <- PFA$adjPvalue[, 1]
      aP[PFA$adjPvalue[, 2]] <- aP
      w1 = -stats::qnorm(aP/2)
      w1 = w1*sign(Test)
      sv = AIC_sda(Y, X, w1)

    }else if(kwd == 'de_lasso'){
      # required by the debias lasso
      X_s = scale(X, center = TRUE, scale = FALSE)
      fit <- glmnet::glmnet(X_s, Y, standardize=FALSE, family = "gaussian")
      k <- fit$df
      AIC <- stats::deviance(fit)+2*k
      i_min = which.min(AIC)
      lambda_select = fit$lambda[i_min]
      # extract coef for a given lambda
      fit_AIC = glmnet::glmnet(X_s, Y, family = "gaussian", lambda = lambda_select)
      beta_hat = fit_AIC$beta[,1]
      #beta_hat = coef(fit, x=X_s, y=Y, s=lambda_select, exact=TRUE)[-1]

      # use the R package: selectiveInference
      # compute fixed lambda test statistics
      sigma <- selectiveInference::estimateSigma(X_s, Y, intercept = TRUE, standardize=FALSE)$sigmahat
      n = length(Y)
      out = selectiveInference::fixedLassoInf(X_s, Y, beta_hat, lambda_select*n, sigma=sigma)

      w1 = rep(0, length(beta_hat))
      w1[out$vars] = out$coef0
      sv = which(w1 != 0)
    }else{
      stop('Incorrect keyword')
    }

    ####refit with Z2 part#####
    if ( length(sv)>0&&length(sv)<p ){
      bt2 = stats::lm(Yt~Xt[,sv]-1)$coefficients
      w2 = rep(0,p)
      w2[sv]=bt2
      sigma_w = rep(1, p)
      DIAG = diag( solve( t(Xt[,sv])%*%Xt[,sv] ) )

      if(scale){
        sd = DIAG
      }else{
        sd = DIAG^(1/2)
      }

      sigma_w[sv] = sd
      W = w1*w2/sigma_w
    }else{
      W = rep(0, p)
    }
    return(W)
  }

  # FDR control based mirror
  W_det <- function(Wj, alpha, options){
    t = sort(abs(Wj))

    if(options=='+'){
      Ta = sapply(t,function(x){(1+length(Wj[Wj<=(-x)]))/max(1,length(Wj[Wj>=x]))})
    }else{
      Ta = sapply(t,function(x){(length(Wj[Wj<=(-x)]))/max(1,length(Wj[Wj>=x]))})
    }

    bestlam = min(t[which(Ta<=alpha)])
    det=which(Wj>=bestlam)
    aa = rep(0, length(Wj))
    aa[det] = 1
    return(aa)
  }

  ##########################################
  # create K=10 sample splits
  # multiple splitting
  p = ncol(dat)
  n = nrow(dat)
  m = as.integer(2/3*n)

  # stable option: single or multiple splitting
  if(stable){
    K=10
  }else{K=1}
  #---------------------------------

  # Need to estimate Omega
  if ( missing(Omega) ){

    if(nonsparse){
      # poet
      n = nrow(dat)
      Q = diag(n)-1/n*rep(1, n)%*%t(rep(1, n))
      Y = t(dat)%*%Q
      K_poet = POET::POETKhat(Y)$K1BN

      Index_K = lapply(1:K, function(x){sample(1:n, m)})
      Omega_K = lapply(1:K, function(x){POET_sig(dat[Index_K[[x]], ], K_poet)})
      Gamma_K = lapply(1:K, function(x){Sqrt(Omega_K[[x]])})
    }else{
      # glasso
      # lambda: tuning parameter for glasso
      out.glasso = huge::huge(dat, method = "glasso")
      out.select = huge::huge.select(out.glasso, criterion = "stars")
      lambda = out.select$opt.lambda
      #
      Index_K = lapply(1:K, function(x){sample(1:n, m)})
      Omega_K = lapply(1:K, function(x){Omega_est(dat[Index_K[[x]], ], lambda)})
      Gamma_K = lapply(1:K, function(x){Sqrt(Omega_K[[x]])})
    }

  }else{
    Index_K = lapply(1:K, function(x){sample(1:n, m)})
    Gamma = Sqrt(Omega)
    Gamma_K = lapply(1:K, function(x){Gamma})
  }


  Wa = sapply(1:K, function(x){ W_value(Gamma_K[[x]], dat, Index_K[[x]], kwd, scale)})
  deta = sapply(1:K, function(x){ W_det(Wa[, x], alpha, '-') })
  deta = as.matrix(deta)
  mdet = apply(deta, 1, sum)
  det1 = which(mdet>=as.integer(0.5*K)+1)
  # the majority vote
  aa=rep(0, p)
  if(length(det1)>0) aa[det1]=1

  # pick the best one
  k_hat = which.min(sapply(1:K,function(x){sum(abs(aa-deta[,x]))}))
  det1=(1:p)[deta[, k_hat]==1]
  if (length(det1)==0){
    out = 'no rejection'
  }else{
    out = det1
  }
  return(out)
}
