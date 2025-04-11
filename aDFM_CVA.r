# # # # # # # # # # # # # # # # #
### routines for the estimation of aDFM models ###
# # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # 
### system definition           ####
# # # # # # # # # # # # # # # # # # 

# specify integers
N = 200
T= 500
r = 3
n = 2
q = 2

# specify system

A = matrix(0,q,q)
A[1,1]=.5
A[2,2]= -0.8

# I(1)
A[2,2]=1

C = diag(q)

K = diag(q) + matrix(rnorm(n*q),n,q)

# ensure that system is invertible
lam = 1
ev <- eigen(A-lam*K%*%C)
while (max(abs(ev$values))>0.99){
  lam = lam/2
  ev <- eigen(A-lam*K%*%C)
}

K <- lam*K

Cext <- rbind(C,matrix(rnorm((r-n)*n),ncol=n))
Dext <- rbind(diag(q),matrix(rnorm((r-q)*q),r-q,q))

# generate loading matrix
LambdaN <- matrix(rnorm(N*r),ncol=r)
LambdaN[1:r,1:r] <- diag(r)


# # # # # # # # # # # # # # # # # 
### simulate data            ####
# # # # # # # # # # # # # # # # # 

# simulate factors 
xt <- matrix(rnorm(n),n,1)
Ft <- matrix(0,r,T)

for (t in 1:T){
  u <- rnorm(q)
  Ft[,t] <- Cext %*% xt + Dext %*% u
  xt <- A %*% xt + K%*%u
}


# simulate idiosyncratic terms as white noise
sigma <- 0.1
xitN <- matrix(rnorm(N*T)*sigma,ncol=T)

# combine the two parts 
ytN <- LambdaN %*% Ft + xitN

# # # # # # # # # # # # # # # # 
### estimation             ####
# # # # # # # # # # # # # # # # 

# PCA
# calculate covariance matrix 
GammaTN <- ytN %*% t(ytN)/(N*T)

# EVs
ev <- eigen(GammaTN)
plot(log(ev$values))

# loadings 
hatLambdaN <- ev$vectors[,1:r]
hatLambdaN <- hatLambdaN %*% solve(hatLambdaN[1:r,1:r])

plot(LambdaN,hatLambdaN)

# PCs
hatFt <- solve( t(hatLambdaN) %*% hatLambdaN) %*% t(hatLambdaN) %*% ytN
plot(Ft,hatFt)

# syst_convert: definition
# converts an unrestricted state space system to a system 
# such that the observability matrix starts with an identity matrix
# (identifiable system representation).
conv_syst <- function(A, C, K){
  
  # check parameters
  if(!is.matrix(A) || !is.matrix(C) || !is.matrix(K)){
    stop("A, C and K must be matrices")
  }
  
  n <- dim(A)[1]  # number of states
  s <- dim(C)[1]  # number of observations
  
  # calculate observability.
  O <- matrix(0, s*n, n)  # initialize observability matrix
  cur <- diag(n)  # initialize identity matrix of size n
  for (j in 1:n){  # iterate over number of states
    O[(j-1)*s+c(1:s), ] <- C %*% cur  # calculate current observability block
    cur <- cur %*% A  # update current state transition
  }
  
  Trafo <- O[1:n,]  # extract transformation matrix
  iTrafo <- solve(Trafo)  # calculate inverse of transformation matrix
  
  At <- Trafo %*% A %*% iTrafo  # transform F matrix
  Ct <- C %*% iTrafo  # transform H matrix
  Kt <- Trafo %*% K  # transform K matrix
  
  # return transformed system matrices
  return(list(A=At, C=Ct, K=Kt))
}

# definition of mlag: provides matrix with rows containing lagged values
# mlag produces a matrix, whose rows are lags of a multivariate process
# This is used in VAR estimation and in the CVA subspace method.
mlag <- function(y, k) {
  
  if(!is.matrix(y)) {
    stop("y must be a matrix")
  }
  
  s <- dim(y)[1] # number of rows (observations) in y
  T <- dim(y)[2] # number of columns (time points) in y
  
  if (k<=0) {
    stop(sprintf("Minimum lag length must be 1, but k=%s is smaller.", k))
  }
  if (k>=T) {
    stop(sprintf("Maximum lag length must be smaller than number of data points, but k=%s is not smaller than #time points: %s.", k, T))
  }
  
  Yk <- matrix(0, s * (k + 1), T - k) # initialize a matrix to store the lags
  Yk[1:s, ] <- y[, (k + 1):T] # copy the non-lagged values to the first s rows
  
  for (j in 1:k) { # loop over each lag
    Yk[(j * s) + 1:s, ] <- y[, (k - j) + 1:(T - k)] # fill in the lagged values
  }
  
  return(Yk) # return the matrix with lags
}


# CVA: definition 
# canonical variate analysis: special kind of subspace algorithm.
# Has been proposed by W. Larimore in 1983, adapted to the case of 
# aDFM to include q<r (less dynamic common factors than static common factors)
CVA_adfm <- function(y,f,p,n,q,we){

  # check parameters
  if (!is.matrix(y)) {
    stop("y must be a matrix")
  }
  if (f <= 0 || p <= 0 || n <= 0 || q <= 0) {
    stop("f, p, n and q must be positive integers")
  }
  kmax <- f+p  # Calculate maximum lags
  s <- dim(y)[1]  # Number of variables in the system
  Teff <- dim(y)[2]-f-p  # Effective sample size after accounting for lags
  
  # set up matrices
  Yk <- mlag(y, kmax)  # Set up matrices for lagged observations
  Yf <- matrix(Yk[1:(s*f), ], ncol=Teff)  # Extract future observations matrix
  Yp <- matrix(Yk[s*f + c(1:(s*p)), ], ncol=Teff)  # Extract past observations matrix
  Hfp <- Yf %*% t(Yp)  # Cross-covariance matrix between future and past observations
  Wf <- diag(f*s)
  if (we == 1){
    Wf <- solve(t(chol(Yf %*% t(Yf))))  # Cholesky decomposition of future covariance matrix
  }
  if (we ==2){
    
  }
  Wp <- chol(Yp %*% t(Yp))  # Cholesky decomposition of past covariance matrix
  
  beta <- Wf %*% Hfp %*% solve(Wp)  # Calculate beta using the inverse of the Cholesky factors
  
  # SVD, singular value decomposition
  usv <- svd(beta)  # Singular value decomposition of beta
  plot((usv$d))
  Khat <- t(usv$v[, 1:n]) %*% solve(t(Wp))  # Estimate state transition matrix Khat
  xhat <- Khat %*% Yp  # Estimate state sequence
  
  # estimate system matrices
  yeff <- matrix(Yk[(f-1)*s + (1:s), ], ncol=Teff)  # Effective observations
  Chat <- yeff %*% t(xhat) %*% solve(xhat %*% t(xhat))  # Estimate observation matrix Hhat
  res <- yeff[, 1:(Teff-1)] - Chat %*% xhat[, 1:(Teff-1)]  # Residuals of the state estimation
  uhat <- res
  
  SigmaT <- res %*% t(res)/Teff  # Covariance matrix of residuals
  
  # reduce innovations, if q<r
  if (q<r){
    Dhat <- rbind(diag(q),SigmaT[(q+1):r,1:q] %*% solve(SigmaT[1:q,1:q]))
    Ddagger <- solve( t(Dhat) %*% Dhat) %*% t(Dhat)
    uhat <- res[1:q,] 
  } else {
    uhat <- res
  }
  
  OmegaT <- uhat %*% t(uhat)/Teff 
  
  xh <- matrix(xhat[, 1:(Teff-1)], ncol=Teff-1)  # Lagged state matrices
  xh1 <- matrix(xhat[, 2:(Teff)], ncol=Teff-1) # Transition the state using the lagged state matrix xh
  uhat <- uhat[,1:(Teff-1)]
  
  Xext <- rbind(xh,uhat)
  AKhat <- xh1 %*% t(Xext) %*% solve(Xext %*% t(Xext))  # Estimate state transition matrix Fhat
  Ahat <- AKhat[,1:n]
  Khat <- AKhat[,(n+1):(n+q)]
  
  # Convert to canonical form
  syst <- conv_syst(A=Ahat, C=Chat, K=Khat)
  
  # account for Dhat 
  Tr <- Dhat[1:q,]
  Dhat <- Dhat %*% solve(Tr)
  
  OmegaT <- Tr %*% OmegaT %*% t(Tr)
  Khat <- syst$K %*% Tr 
  
  # Return system matrices and covariance
  return(list(Ahat=syst$A, Chat=syst$C, Khat=Khat, Dhat = Dhat, Omega=OmegaT, Sigma = SigmaT))
}



## CVA: application 
estimates <- CVA_adfm(y=hatFt,f=5,p=5,n=2,q=2,we = 1)


# # # # # # # # # # # # # # # # # 
### compare results.         ####                 
# # # # # # # # # # # # # # # # #
(A - estimates$Ahat)
(Cext - estimates$Chat)
(Dext - estimates$Dhat)
(K - estimates$Khat)

