#' Symmetrized Data Aggregation for two-sample t-test
#'
#' @param dat_I a \eqn{n_1} by \eqn{p} data matrix, the first part of data
#' @param dat_II a \eqn{n_2} by \eqn{p} data matrix, the second part of data
#' @param alpha the FDR level
#' @param Sigma_I the covariance matrix of sample 1; if it is missing, it will be estimated
#' by the glasso package.
#' @param Sigma_II the covariance matrix of sample 2; if it is missing, it will be estimated
#' by the glasso package.
#' @param stable If it is TRUE, the sample will be randomly splitted \eqn{B=10} times for stability
#' performance; otherwise, only single sample splitting is used.
#'
#' @return the indices of the hypotheses rejected
#' @export
#' @examples
#' p = 100
#' n = 30
#' dat_I = matrix(rnorm(n*p),nrow = n)
#' mu = rep(0, p)
#' mu[1:10] = 1.5
#' dat_I = dat_I = rep(1, n)%*%t(mu)
#'
#' dat_II = matrix(rnorm(n*p), nrow = n)
#' Sigma_I = diag(p)
#' Sigma_II = diag(p)
#' out = SDA_2S(dat_I, dat_II, alpha=0.05, Sigma_I, Sigma_II)
#' print(out)
#'
SDA_2S <- function(dat_I, dat_II, alpha, Sigma_I, Sigma_II, stable=TRUE){
  # SDA method with two sample t-tests
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

  # de-mean process
  de_mean = function(dat){
    n = nrow(dat)
    out = (diag(n)-1/n*rep(1, n)%*%t(rep(1, n)))%*%dat
    return(out) # n by p
  }

  # re-construction of two data
  dat_RE <- function(dat1, dat2){
    n_1 = nrow(dat1)
    n_2 = nrow(dat2)
    dat1 = de_mean(dat1)
    dat2 = de_mean(dat2)

    dat = NULL
    for(i in 1:n_1){
      for(j in 1:n_2){
        V =  dat1[i, ]-sqrt(n_1/n_2)*dat2[j, ]
        dat = rbind(dat, V)
      }
    }
    dat = as.matrix(dat) # n_1*n_2 by p
    return(dat)
  }

  # Sample covariance matrix
  Sample_COV <- function(dat){
    n = nrow(dat)
    Q = diag(n)-1/n*rep(1, n)%*%t(rep(1, n))
    out = t(dat)%*%Q%*%dat/(n-1)
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

  # LASSO for screening
  Model_Select <- function(Y, X){
    # Method I: AIC
    fit <- glmnet::glmnet(X, Y, family = "gaussian")
    k <- fit$df
    AIC <- stats::deviance(fit)+2*k
    i_min = which.min(AIC)
    lambda_select = fit$lambda[i_min]
    fit_AIC = glmnet::glmnet(X, Y, family = "gaussian", lambda = lambda_select)
    w1 = fit_AIC$beta[,1]
    sv = as.vector( which(w1!=0) )

    # method II: upper bound p/3
    k= ceiling(p/3)
    if( length(sv)>k ){
      wv = which(fit$df==max(fit$df[fit$df<k]))[1]
      sv = which(fit$beta[, wv]!=0)
      w1 = fit$beta[, wv]
    }

    return( list(index_S=sv, w1 = w1))
  }

  # mean difference of two sample: Y and X
  Y_X <- function(Gamma, dat1, dat2){
    # Gamma = Omega^(-1/2)
    n_1 = nrow(dat1)
    n_2 = nrow(dat2)
    n = n_1*n_2/(n_1+n_2)

    Z_1 = as.vector( apply(dat1, 2, mean) )
    Z_2 = as.vector( apply(dat2, 2, mean) )
    Y = as.vector( sqrt(n)*Gamma%*%(Z_1-Z_2) )
    X = as.vector(Y)
    X = sqrt(n)*Gamma
    return( list(Y=Y, X=X))
  }

  # generate W_j: splitting
  MF_S<- function(Gamma, dat_I, index_I, dat_II, index_II){
    # Gamma: only from the first sample
    # from the first part of data
    out = Y_X(Gamma, dat_I[index_I, ], dat_II[index_II, ])
    X = out$X
    Y = out$Y
    # from the second sample
    out = Y_X(Gamma, dat_I[-index_I, ], dat_II[-index_II, ])
    Xt = out$X
    Yt = out$Y
    #
    MS = Model_Select(Y, X)
    sv = MS$index_S
    w1 = MS$w1
    ####refit with Z1 part#####
    if ( length(sv)>0&&length(sv)<p ){
      bt2 = stats::lm(Yt~Xt[,sv]-1)$coefficients
      w2 = rep(0,p)
      w2[sv]=bt2
      sigma_w = rep(1, p)
      DIAG = diag( solve( t(Xt[,sv])%*%Xt[,sv] ) )
      sigma_w[sv] = DIAG
      Wj = w1*w2/sigma_w
    }else{
      Wj = rep(0, p)
    }
    return(Wj)
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

  # Method: multiple splitting
  Mirror_MM_new <- function(Gamma_K, dat_I, index_I_K, dat_II, index_II_K, alpha){

    p = ncol(dat_I)
    K = length(index_I_K)
    Wa = sapply(1:K, function(x){ MF_S(Gamma_K[[x]], dat_I, index_I_K[[x]], dat_II, index_II_K[[x]]) })
    deta = sapply(1:K,function(x){ W_det(Wa[, x], alpha, '-') })
    mdet = apply(deta, 1, sum)
    det1 = which(mdet>=as.integer(0.5*K)+1)

    aa=rep(0, p)
    if(length(det1)>0) aa[det1]=1

    k_hat = which.min(sapply(1:K,function(x){sum(abs(aa-deta[,x]))}))
    det1=(1:p)[deta[, k_hat]==1]
    # rank the output
    W_pick = Wa[, k_hat]
    det1 = det1[order(W_pick[det1], decreasing = TRUE)]
    return(  sort(det1) )
  }

  ##########################################
  # Main file
  n_1 = dim(dat_I)[1]
  n_2 = dim(dat_II)[1]
  n = n_1+n_2
  p = dim(dat_I)[2]

  # stable option: single or multiple splitting
  if(stable){
    K=10
  }else{K=1}
  #---------------------------------

  if( missing(Sigma_I)||missing(Sigma_II) ){
    # fix the tuning parameter
    index_1 = sample(1:n_1, floor(2/3*n_1))
    index_2 = sample(1:n_2, floor(2/3*n_2))
    dat1 = dat_I[index_1, ]
    dat2 = dat_II[index_2, ]
    # estimate the pooling covariance matrix: tuning parameter
    dat_combine = dat_RE(dat1, dat2)
    # lambda: tuning parameter for glasso
    out.glasso = huge::huge(dat_combine, method = "glasso")
    out.select = huge::huge.select(out.glasso, criterion = "ebic")
    lambda = out.select$opt.lambda

    index_I_K = lapply(1:K, function(x){sample(1:n_1, floor(2/3*n_1))})
    index_II_K = lapply(1:K, function(x){sample(1:n_2, floor(2/3*n_2))})
    Gamma_K = lapply(1:K, function(x){
      index_1 = index_I_K[[x]]
      index_2 = index_II_K[[x]]
      dat1 = dat_I[index_1, ]
      dat2 = dat_II[index_2, ]
      dat12 = dat_RE(dat1, dat2)
      Sigma_MM = Sample_COV(dat12) # only use as input
      Omega_hat = Omega_est(Sigma_MM, lambda)
      Omega_hat = (1+length(index_1)/length(index_2))*Omega_hat
      Gamma_hat = Sqrt(Omega_hat)
      return(Gamma_hat)
    })
  }else{
    Omega <- INV(n_2/n*Sigma_I+n_1/n*Sigma_II)
    index_I_K = lapply(1:K, function(x){sample(1:n_1, floor(2/3*n_1))})
    index_II_K = lapply(1:K, function(x){sample(1:n_2, floor(2/3*n_2))})
    Gamma = Sqrt(Omega)
    Gamma_K = lapply(1:K, function(x){Gamma})
  }


  det_01 = Mirror_MM_new(Gamma_K, dat_I, index_I_K, dat_II, index_II_K, alpha)



  # output
  if (length(det_01)==0){
    out = 'no rejection'
  }else{
    out = det_01
  }
  return(out)
}
