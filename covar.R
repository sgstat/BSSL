library(mvtnorm)
library(splines)
library(tmvmixnorm)
library(statmod)
library(quadprog)
library(matrixcalc)
library(Matrix)
library(matrixStats)
library(dplyr)

# Read the dataset
data <- read.csv("hcc_cov.csv", header = TRUE, stringsAsFactors = FALSE)
data$Gender <- as.factor(data$Gender)
Zs = model.matrix(~ Gender - 1, data = data)

data = data[-24]
data = cbind(Zs, data)

# Ensure the label column is a factor
data$Class <- as.factor(data$Class)


scaledX = function(X){
  X = as.matrix(X)
  means = apply(X, 2, mean, na.rm=TRUE)
  sds = apply(X, 2, sd)
  ncols = ncol(X)
  do.call(cbind,
          lapply(1:ncols, function(z){
            # print(X[, z])
            # ifelse((min(X[, z]) >= 0 && max(X[, z]) <= 1), X[, z], pnorm(X[, z], mean=means[z], sd=sds[z]))
            pnorm(as.numeric(X[, z]), mean=means[z], sd=sds[z])
          }))
}

# args = commandArgs(trailingOnly = TRUE)
# argslen = length(args)
# if(argslen > 1) stop('Error: Too Many Arguments')
# if(argslen < 1) stop('Error: Takes One Arguments p, n.total and n.label')
# 
# args1 = dict[[as.numeric(args)]]


set.seed(999)
seed=sample(1:1000, 10)
mcc=fpr=fnr=rep(NA, 10)
N = nrow(data)
n = 100

# Hamiltonian Monte Carlo
HMC_exact <- function(F, g, M, mu_r, cov, L, initial_X) {
  
  # Implementation of the algorithm described in http://arxiv.org/abs/1208.4118
  # Author: Ari Pakman
  # Translated into R
  
  # Returns samples from a d-dimensional Gaussian with m constraints given by F * X + g > 0
  # If cov == TRUE
  # then M is the covariance and the mean is mu = mu_r
  # if cov == FALSE
  # then M is the precision matrix and the log-density is -1/2 * X' * M * X + r' * X
  
  # Input
  # F: m x d matrix
  # g: m x 1 vector
  # M: d x d matrix, must be symmetric and positive definite
  # mu_r: d x 1 vector
  # cov: see explanation above
  # L: number of samples desired
  # initial_X: d x 1 vector. Must satisfy the constraint.
  
  # Output
  # Xs: d x L matrix, each column is a sample
  # bounce_count: number of times the particle bounced
  
  # Go to a whitened frame
  
  if(is.null(nrow(g)))
    g = matrix(g, ncol = 1)
  if(is.null(nrow(initial_X)))
    initial_X = matrix(initial_X, ncol = 1)
  
  m = nrow(g)
  if (nrow(F) != m) {
    cat("error\n")
    return(NULL)
  }
  
  if (cov) {
    mu = mu_r
    g = g + F %*% mu
    R = chol(M)
    F = F %*% t(R)
    initial_X = initial_X - mu
    initial_X = solve(t(R), initial_X)
  } else {
    r = mu_r
    R = chol(M)      # M = R' * R
    mu = solve(R, solve(t(R), r))
    g = g + F %*% mu
    F = F %*% solve(R)  # This is the time-consuming step in this code section
    initial_X = initial_X - mu
    initial_X = R %*% initial_X
  }
  
  d = nrow(initial_X)
  Xs = matrix(NA, nrow = d, ncol = L)
  bounce_count = 0
  nearzero = 1e+04 * .Machine$double.eps
  
  # Verify that initial_X is feasible
  c = F %*% initial_X + g
  if (any(c < 0)) {
    cat("error: inconsistent initial condition\n")
    return(NULL)
  }
  
  # Unsparcify the linear constraints
  g = as.vector(g)
  F2 = rowSums(F^2)  # Squared norm of the rows of F, needed for reflecting the velocity
  F = as.matrix(F)
  Ft = t(F)
  
  # Sampling loop
  last_X = initial_X
  Xs[, 1] = initial_X[, 1]
  
  i = 2
  while (i <= L) {
    stop = 0
    j = 0
    V0 = rnorm(d)  # Initial velocity
    X = last_X
    
    T = pi / 2  # Total time the particle will move
    tt = 0  # Records how much time the particle already moved
    
    while (TRUE) {
      a = V0
      a = Re(a)
      b = X
      
      fa = F %*% a
      fb = F %*% b
      
      U = sqrt(fa^2 + fb^2)
      phi = atan2(-fa, fb)  # -pi < phi < +pi
      
      pn = abs(g / U) <= 1  # These are the walls that may be hit
      
      # Find the first time constraint becomes zero
      if (any(pn)) {
        inds = which(pn)
        
        phn = phi[pn]
        
        t1 = -phn + acos(-g[pn] / U[pn])  # Time at which coordinates hit the walls
        
        # If there was a previous reflection (j > 0) and there is a potential reflection at the sample plane
        # Make sure that a new reflection at j is not found because of numerical error
        if (j > 0) {
          if (pn[j] == 1) {
            cs = cumsum(pn)
            indj = cs[j]
            tt1 = t1[indj]
            if (abs(tt1) < nearzero || abs(tt1 - 2 * pi) < nearzero) {
              t1[indj] = Inf
            }
          }
        }
        
        mt = min(t1)
        m_ind = which.min(t1)
        
        # Find the reflection plane
        j = inds[m_ind]  # j is an index in the full vector of dim-m, not in the restricted vector determined by pn
        
      } else {  # if pn(i) == 0 for all i
        mt = T
      }
      
      tt = tt + mt
      
      if (tt >= T) {
        mt = mt - (tt - T)
        stop = 1
      }
      
      # Move the particle a time mt
      X = a * sin(mt) + b * cos(mt)
      V = a * cos(mt) - b * sin(mt)
      
      if (stop) {
        break
      }
      
      # Compute reflected velocity
      qj = F[j, ] %*% V / F2[j]
      V0 = V - 2 * as.vector(qj) * Ft[, j]
      bounce_count = bounce_count + 1
      
    }  # End of inner while loop
    
    # Check that the candidate X satisfies the constraints before accepting it
    if (all(F %*% X + g > 0)) {
      Xs[, i] = X
      last_X = X
      i = i + 1
    } else {
      cat("hmc reject\n")
    }
    
  }  # End of outer while loop
  
  # Transform back to the unwhitened frame
  if (cov) {
    Xs = t(R) %*% Xs + matrix(mu, nrow = d, ncol = L, byrow = TRUE)
  } else {
    Xs = solve(R, Xs) + matrix(mu, nrow = d, ncol = L, byrow = TRUE)
  }
  
  return(list(Xs = Xs, bounce_count = bounce_count))
} # end of function


# cat('\n Running Job with ', dict[1],' data')


for(sim.num in 1:1){
  sim=seed[sim.num]
  cat("\n seed is = ", sim, "for simulation = ", sim.num, '\n')
  set.seed(sim)
  niter=500
  
  ################# process data ##########################
  data[, -25] = as.data.frame(scaledX(data[, -25]))
  # Split the data by class and sample 30 instances per class for training  # For reproducibility
  train_data <- data %>% group_by(Class) %>% slice_sample(n = 50)
  test_data <- anti_join(data, train_data) %>% arrange(Class)  # Assuming "ID" is a unique identifier
  
  # Introduce missing labels in the training set
  train_data <- train_data %>%
    group_by(Class) %>%
    mutate(Class_Missing = ifelse(row_number() <= 10, Class, NA))  # Only 5 observed per class
  
  train_data = as.data.frame(train_data)
  test_data = as.data.frame(test_data)
  
  x = train_data[, -c(25, 26, 1, 2)]
  p = ncol(x)
  print(dim(x))
  train = as.matrix(test_data[, -c(25, 1, 2)])
  # print(dim(X_test))
  Ltrain.true = as.numeric(factor(test_data[, 25])) - 1
  L.true =  as.numeric(factor(train_data[[25]])) - 1 
  L.org = as.numeric(factor(train_data[[26]])) - 1
  L.org[is.na(L.org)] = 2
  # df=data.frame(x, L.org)
  # df.test=data.frame(X_test, L.test)
  # covariates
  zs = as.matrix(train_data[, c(1,2)]) # cbind(1, train_data[, 2])
  train.zs = as.matrix(test_data[, c(1,2)]) # cbind(1, test_data[, 2])
  ncov = ncol(zs)
  
  midpoint=rep(0, 15)
  
  for(J in 8:15){
    knotvec=c(rep(0,3), seq(0, 1, by=1/J), rep(1,3))
    A.bas=splineDesign(knots=knotvec, x = c(1/4, 1/2, 3/4))
    A=rbind(A.bas[2,], A.bas[3,]-A.bas[1,])
    
    J1=dim(A)[2]
    mu.org=qnorm((1:J1-0.375)/(J1+0.25))
    mu.prior=mu.org+t(A)%*%solve(A%*%t(A))%*%(matrix(c(0,1),2,1)-A%*%mu.org)
    Sigma.prior=5*(diag(rep(1, J1))-t(A)%*%solve(A%*%t(A))%*%A)
    
    # solve A*theta=(0 1)^T
    index1=which.max(A[1,])
    index2=which.max(A[2,])
    if(index1==index2){index2=index2+1}
    
    a.index1=A[1,index1]
    a.index2=A[1,index2]
    b.index1=A[2,index1]
    b.index2=A[2,index2]
    
    W=W.org=A[, -c(index1, index2)]
    W[2,]=(W.org[2,]-W.org[1,]*b.index1/a.index1)/-(b.index2-a.index2*b.index1/a.index1)
    W[1,]=-W.org[1,]/a.index1-W[2,]*a.index2/a.index1
    q=c(1/(b.index2-a.index2*b.index1/a.index1)*(-a.index2/a.index1), 1/(b.index2-a.index2*b.index1/a.index1))
    
    # check the calculation
    A[1,c(index1, index2)]%*%W+c(A[1,-c(index1, index2)])
    A[1,c(index1, index2)]%*%q
    A[2,c(index1, index2)]%*%W+c(A[2,-c(index1, index2)])
    A[2,c(index1, index2)]%*%q
    
    mu.prior.r=mu.prior[-c(index1, index2)]
    Sigma.prior.r=Sigma.prior[-c(index1, index2), -c(index1, index2)]
    Sigma.prior.r.inv=chol2inv(chol(Sigma.prior.r))
    
    newindex=function(index){
      count=(index<index1)+(index<index2)
      if(count==2){return(index)}
      if(count==1){return(index-1)}
      if(count==0){return(index-2)}
    }
    
    F.org=cbind(diag(rep(-1, J1-1)), rep(0, J1-1))+cbind(rep(0, J1-1), diag(rep(1, J1-1)))
    F=matrix(NA, 1,J1-2)
    g=rep(0, J1)
    for(i in 2:J1){
      new=rep(0, J1-2)
      if(i==index1){
        new=W[1,]
        g[i]=q[1]
      }else if(i==index2){
        new=W[2,]
        g[i]=q[2]
      }else{
        new[newindex(i)]=1
      }
      if(i-1==index1){
        new=new-W[1,]
        g[i]=g[i]-q[1]
      }else if(i-1==index2){
        new=new-W[2,]
        g[i]=g[i]-q[2]
      }else{
        new[newindex(i-1)]=new[newindex(i-1)]-1
      }
      F=rbind(F, new)
    }
    F=F[-1,]
    g=g[-1]
    
    theta.iter=array(dim=c(niter, J1-2, p))
    B=array(dim=c(n, p, J1-2)) 
    B.star=array(dim=c(n, p, 2))
    for(d in 1:p){
      middle=splineDesign(knots=knotvec, x = x[,d])
      B[,d,] = middle[,-c(index1, index2)]
      B.star[,d,] = middle[,c(index1, index2)]
    }
    
    # initial for theta
    out=gauss.quad(20,kind="hermite",alpha=0,beta=0)
    bspline.matrix.cdf=splineDesign(knots = knotvec, x=pnorm(out$nodes))
    fn.LHS=bspline.matrix.cdf*out$nodes
    fn.LHS.weighted=out$weights*fn.LHS
    LHS.gq=colSums(fn.LHS.weighted)
    
    fn.RHS.weighted=matrix(0, J1, J1)
    RHS.gq=matrix(0, J1, J1)
    for(pts in 1:20){
      for(j in 1:J1){
        for(k in 1:J1){
          fn.RHS.weighted[j,k]=out$weights[pts]*bspline.matrix.cdf[pts, j]*bspline.matrix.cdf[pts, k]
        }
      }
      RHS.gq=RHS.gq+fn.RHS.weighted
    }
    theta.initial=solve.QP(nearPD(RHS.gq%*%RHS.gq)$'mat', as.numeric(LHS.gq%*%RHS.gq), Amat=t(rbind(A, F.org)), bvec=c(0, 1, rep(10^(-4), J1-1)), meq=2)$'solution'
    theta=matrix(rep(theta.initial[-c(index1, index2)], p), ,p)
    
    Y=matrix(NA,n,p)
    for(d in 1:p){
      Y[,d]=(B[,d,]+B.star[,d,]%*%W)%*%theta[,d]+B.star[,d,]%*%matrix(q)
    }
    
    # initial values
    if(sum(L.org==0)==1){
      mu1=Y[L.org==0,]
    }else{
      mu1=colMeans(Y[L.org==0,]) 
    }
    
    if(sum(L.org==1)==1){
      mu2=Y[L.org==1,]
    }else{
      mu2=colMeans(Y[L.org==1,]) 
    }
    
    
    # Assign initial labels according to distance
    L=L.org
    for(i in 1:n){
      if(L.org[i]==2){
        d1=sum((Y[i,]-mu1)^2)
        d2=sum((Y[i,]-mu2)^2)
        if(d1<d2){
          L[i]=0
        }else{
          L[i]=1
        }
      }
    }
    n1=sum(L==0)
    n2=sum(L==1)
    cat("\n No. of missclassifications initially: ", sum(L!=L.true), "for J = ", J)
    
    ################check the original L#################
    
    Sigma1=cov(Y[L==0,])
    Sigma2=cov(Y[L==1,])
    Sigma1inv=chol2inv(chol(Sigma1))
    is.positive.definite(Sigma1inv)
    Sigma2inv=chol2inv(chol(Sigma2))
    is.positive.definite(Sigma2inv)
    Psi=matrix(0, p, ncov)
    
    mu1.iter=matrix(,niter,p)
    Sigma1.iter=array(dim=c(niter, p,p))
    Sigma1inv.iter=array(dim=c(niter, p,p))
    mu2.iter=matrix(,niter,p)
    Sigma2.iter=array(dim=c(niter, p,p))
    Sigma2inv.iter=array(dim=c(niter, p,p))
    lambda.iter=array(niter, 1)
    L.iter=matrix(,niter,n)
    Psi.iter=matrix(,niter,p*ncov)
    
    sigma.psi = 10
    
    mean.prior=mu.prior.r%*%Sigma.prior.r.inv
    mid=array(dim=c(n, p, J1-2)) 
    midsq=array(dim=c(n, p, J1-2, J1-2))
    mid2=matrix(,n, p)
    
    for(d in 1:p){
      mid[,d,]=B[,d,]+B.star[,d,]%*%W
      mid2[,d]=B.star[,d,]%*%matrix(q)
      for(i in 1:n){
        midsq[i,d,,]=(mid[i,d,])%*%t(mid[i,d,])
      }
    }
    
    
    for(iter in 1:niter){
      psi = zs %*% t(Psi) # n x c
      for(d in 1:p){
        Cov=matrix(0, J1-2, J1-2)
        mean1=0
        mean2=0
        var1=matrix(0,J1-2,J1-2)
        var2=matrix(0,J1-2,J1-2)
        for(i in 1:n){
          if(L[i]==0){
            mean1=mean1+mid[i,d,]*((mid2[i,d]-mu1[d]-psi[i,d])*Sigma1inv[d, d]+sum(Sigma1inv[-d, d]*(Y[i,-d]-mu1[-d]-psi[i,-d])))
            var1=var1+midsq[i,d,,]
          }
          if(L[i]==1){
            mean2=mean2+mid[i,d,]*((mid2[i,d]-mu2[d]-psi[i,d])*Sigma2inv[d, d]+sum(Sigma2inv[-d, d]*(Y[i,-d]-mu2[-d]-psi[i,-d])))
            var2=var2+midsq[i,d,,]
          }
        }
        Cov=Sigma1inv[d, d]*var1+Sigma2inv[d, d]*var2+Sigma.prior.r.inv
        if(!is.positive.definite(Cov)){
          Cov=nearPD(Cov)$'mat'
        }
        mean.final=mean.prior-mean1-mean2
        Cov=chol2inv(chol(Cov))
        if(!is.positive.definite(Cov)){
          Cov=nearPD(Cov)$'mat'
        }
        result = HMC_exact(F=F, g=g, M = Cov, mu_r=as.numeric(Cov%*%matrix(mean.final)), cov=TRUE, L=1, initial_X=theta[, d])
        theta[, d] = as.numeric(result$Xs)
        Y[,d]=(B[,d,]+B.star[,d,]%*%W)%*%theta[,d]+B.star[,d,]%*%matrix(q)
      }
      L.iter[iter,]=L
      theta.iter[iter,,]=theta
      # print(dim(as.matrix(kronecker(Sigma1inv, crossprod(zs[L==0,])))))
      # Sample Psi
      Sigma.psi = chol2inv(chol( sigma.psi^(-2)*diag(p*ncov) + 
                                     as.matrix(kronecker(Sigma1inv, crossprod(zs[L==0,]))) + 
                                     as.matrix(kronecker(Sigma2inv, crossprod(zs[L==1,]))) ))
      # print("check")
      # Sigma.psi = as.matrix(nearPD(Sigma.psi)$`mat`)
      mu.psi = Sigma.psi %*% (c(Sigma1inv%*%t(sweep(Y[L==0,], 2, mu1, "-"))%*%zs[L==0,]) 
                              + c(Sigma2inv%*%t(sweep(Y[L==1,], 2, mu2, "-"))%*%zs[L==1,]))
      
      Psi.iter[iter, ] = psi_vec = chol(Sigma.psi) %*% rnorm(p*ncov) + mu.psi
      Psi = matrix(psi_vec, p, ncov)
      
      # print(Psi)
      sigma.psi = sqrt(1/rgamma(1, 1+p*ncov/2, 1+sum(Psi^2)/2))
      # random walk metropolis update for sigma.psi
      # log_lik = function(Y, )
      
      Sigma1inv.iter[iter,,]=Sigma1inv=as.matrix(nearPD(rWishart(1, p+2+n1, chol2inv(chol(cov(Y[L==0,])*(n1-1)+diag(p))))[,,1])$`mat`)
      Sigma1.iter[iter,,]=Sigma1=chol2inv(chol(Sigma1inv))
      mu1.iter[iter,]=  mu1=rmvnorm(1, colMeans(Y[L==0,]), Sigma1/n1)
      Sigma2inv.iter[iter,,]=Sigma2inv=as.matrix(nearPD(rWishart(1, p+2+n2, chol2inv(chol(cov(Y[L==1,])*(n2-1)+diag(p))))[,,1])$`mat`)
      Sigma2.iter[iter,,]=  Sigma2=chol2inv(chol(Sigma2inv))
      mu2.iter[iter,]=  mu2=rmvnorm(1, colMeans(Y[L==1,]), Sigma2/n2)
      lambda.iter[iter] = lambda = rbeta(1, 1+n1, 1+n2) 
      
      for(i in 1:n){
        if(L.org[i]==2){
          log_probs = c(dmvnorm(Y[i,]-zs[i,, drop=FALSE]%*%t(Psi), mu1, Sigma1, log=T), dmvnorm(Y[i,]-zs[i,, drop=FALSE]%*%t(Psi), mu2, Sigma2, log=T))
          log_probs = log_probs + log(c(lambda, 1-lambda))
          prob = exp(log_probs - logSumExp(log_probs))[2]
          L[i] = rbinom(1, 1, prob)
        }
      }
      n1=sum(L==0)
      n2=sum(L==1)
      if(iter%%100==0){
        cat("\n Iterations :===>", iter)
        cat("\n No. of missclassifications: ", sum(L!=L.true))
      }
    }
    theta.final=apply(theta.iter[(niter/5):niter,,], 2:3, mean)
    mu1.final=apply(mu1.iter[(niter/5):niter,], 2, mean)
    Sigma1.final=apply(Sigma1.iter[(niter/5):niter,,], 2:3, mean)
    mu2.final=apply(mu2.iter[(niter/5):niter,], 2, mean)
    Sigma2.final=apply(Sigma2.iter[(niter/5):niter,,], 2:3, mean)
    lambda.final=mean(lambda.iter[(niter/5):niter])
    Psi.final=matrix(apply(Psi.iter[(niter/5):niter,], 2, mean), p, ncov)
    
    for(d in 1:p){
      Y[,d]=(B[,d,]+B.star[,d,]%*%W)%*%theta.final[,d]+B.star[,d,]%*%matrix(q)
    }
    
    p1=dmvnorm(Y-zs%*%t(Psi.final), mu1.final, Sigma1.final)
    p2=dmvnorm(Y-zs%*%t(Psi.final), mu2.final, Sigma2.final)
    r=(lambda.final * p1)/((1-lambda.final) * p2)
    midpoint[J]=sum(r[which(r>1/3)]<3)
  }
  
  ##################### MCMC after selecting J #######################
  niter=10000
  J=7+which.min(midpoint[8:15])
  knotvec=c(rep(0,3), seq(0, 1, by=1/J), rep(1,3))
  A.bas=splineDesign(knots=knotvec, x = c(1/4, 1/2, 3/4))
  A=rbind(A.bas[2,], A.bas[3,]-A.bas[1,])
  
  J1=dim(A)[2]
  mu.org=qnorm((1:J1-0.375)/(J1+0.25))
  mu.prior=mu.org+t(A)%*%solve(A%*%t(A))%*%(matrix(c(0,1),2,1)-A%*%mu.org)
  Sigma.prior=5*(diag(rep(1, J1))-t(A)%*%solve(A%*%t(A))%*%A)
  
  # solve A*theta=(0 1)^T
  index1=which.max(A[1,])
  index2=which.max(A[2,])
  if(index1==index2){index2=index2+1}
  
  a.index1=A[1,index1]
  a.index2=A[1,index2]
  b.index1=A[2,index1]
  b.index2=A[2,index2]
  
  W=W.org=A[, -c(index1, index2)]
  W[2,]=(W.org[2,]-W.org[1,]*b.index1/a.index1)/-(b.index2-a.index2*b.index1/a.index1)
  W[1,]=-W.org[1,]/a.index1-W[2,]*a.index2/a.index1
  q=c(1/(b.index2-a.index2*b.index1/a.index1)*(-a.index2/a.index1), 1/(b.index2-a.index2*b.index1/a.index1))
  
  # check the calculation
  A[1,c(index1, index2)]%*%W+c(A[1,-c(index1, index2)])
  A[1,c(index1, index2)]%*%q
  A[2,c(index1, index2)]%*%W+c(A[2,-c(index1, index2)])
  A[2,c(index1, index2)]%*%q
  
  mu.prior.r=mu.prior[-c(index1, index2)]
  Sigma.prior.r=Sigma.prior[-c(index1, index2), -c(index1, index2)]
  Sigma.prior.r.inv=chol2inv(chol(Sigma.prior.r))
  
  newindex=function(index){
    count=(index<index1)+(index<index2)
    if(count==2){return(index)}
    if(count==1){return(index-1)}
    if(count==0){return(index-2)}
  }
  
  F.org=cbind(diag(rep(-1, J1-1)), rep(0, J1-1))+cbind(rep(0, J1-1), diag(rep(1, J1-1)))
  F=matrix(, 1,J1-2)
  g=rep(0, J1)
  for(i in 2:J1){
    new=rep(0, J1-2)
    if(i==index1){
      new=W[1,]
      g[i]=q[1]
    }else if(i==index2){
      new=W[2,]
      g[i]=q[2]
    }else{
      new[newindex(i)]=1
    }
    if(i-1==index1){
      new=new-W[1,]
      g[i]=g[i]-q[1]
    }else if(i-1==index2){
      new=new-W[2,]
      g[i]=g[i]-q[2]
    }else{
      new[newindex(i-1)]=new[newindex(i-1)]-1
    }
    F=rbind(F, new)
  }
  F=F[-1,]
  g=g[-1]
  
  theta.iter=array(dim=c(niter, J1-2, p))
  B=array(dim=c(n, p, J1-2)) 
  B.star=array(dim=c(n, p, 2))
  for(d in 1:p){
    middle=splineDesign(knots=knotvec, x = x[,d])
    B[,d,] = middle[,-c(index1, index2)]
    B.star[,d,] = middle[,c(index1, index2)]
  }
  
  # initial for theta
  out=gauss.quad(20,kind="hermite",alpha=0,beta=0)
  bspline.matrix.cdf=splineDesign(knots = knotvec, x=pnorm(out$nodes))
  fn.LHS=bspline.matrix.cdf*out$nodes
  fn.LHS.weighted=out$weights*fn.LHS
  LHS.gq=colSums(fn.LHS.weighted)
  
  fn.RHS.weighted=matrix(0, J1, J1)
  RHS.gq=matrix(0, J1, J1)
  for(pts in 1:20){
    for(j in 1:J1){
      for(k in 1:J1){
        fn.RHS.weighted[j,k]=out$weights[pts]*bspline.matrix.cdf[pts, j]*bspline.matrix.cdf[pts, k]
      }
    }
    RHS.gq=RHS.gq+fn.RHS.weighted
  }
  theta.initial=solve.QP(nearPD(RHS.gq%*%RHS.gq)$'mat', as.numeric(LHS.gq%*%RHS.gq), Amat=t(rbind(A, F.org)), bvec=c(0, 1, rep(10^(-4), J1-1)), meq=2)$'solution'
  theta=matrix(rep(theta.initial[-c(index1, index2)], p), ,p)
  
  Y=matrix(,n,p)
  for(d in 1:p){
    Y[,d]=(B[,d,]+B.star[,d,]%*%W)%*%theta[,d]+B.star[,d,]%*%matrix(q)
  }
  
  # initial values
  if(sum(L.org==0)==1){
    mu1=Y[L.org==0,]
  }else{
    mu1=colMeans(Y[L.org==0,]) 
  }
  
  if(sum(L.org==1)==1){
    mu2=Y[L.org==1,]
  }else{
    mu2=colMeans(Y[L.org==1,]) 
  }
  
  
  # Assign initial labels according to distance
  L=L.org
  for(i in 1:n){
    if(L.org[i]==2){
      d1=sum((Y[i,]-mu1)^2)
      d2=sum((Y[i,]-mu2)^2)
      if(d1<d2){
        L[i]=0
      }else{
        L[i]=1
      }
    }
  }
  
  n1=sum(L==0)
  n2=sum(L==1)
  sum(L!=L.true)
  
  ################check the original L#################
  
  Sigma1=cov(Y[L==0,])
  Sigma2=cov(Y[L==1,])
  Sigma1inv=chol2inv(chol(Sigma1))
  is.positive.definite(Sigma1inv)
  Sigma2inv=chol2inv(chol(Sigma2))
  is.positive.definite(Sigma2inv)
  Psi = matrix(0, p, ncov)
  
  mu1.iter=matrix(,niter,p)
  Sigma1.iter=array(dim=c(niter, p,p))
  Sigma1inv.iter=array(dim=c(niter, p,p))
  mu2.iter=matrix(,niter,p)
  Sigma2.iter=array(dim=c(niter, p,p))
  Sigma2inv.iter=array(dim=c(niter, p,p))
  L.iter=matrix(,niter,n)
  Psi.iter=matrix(,niter,p*ncov)
  
  sigma.psi = 10
  
  mean.prior=mu.prior.r%*%Sigma.prior.r.inv
  mid=array(dim=c(n, p, J1-2)) 
  midsq=array(dim=c(n, p, J1-2, J1-2))
  mid2=matrix(,n, p)
  
  for(d in 1:p){
    mid[,d,]=B[,d,]+B.star[,d,]%*%W
    mid2[,d]=B.star[,d,]%*%matrix(q)
    for(i in 1:n){
      midsq[i,d,,]=(mid[i,d,])%*%t(mid[i,d,])
    }
  }
  
  
  for(iter in 1:niter){
    for(d in 1:p){
      Cov=matrix(0, J1-2, J1-2)
      mean1=0
      mean2=0
      var1=matrix(0,J1-2,J1-2)
      var2=matrix(0,J1-2,J1-2)
      for(i in 1:n){
        psi = zs %*% t(Psi)
        if(L[i]==0){
          mean1=mean1+mid[i,d,]*((mid2[i,d]-mu1[d]-psi[i,d])*Sigma1inv[d, d]+sum(Sigma1inv[-d, d]*(Y[i,-d]-mu1[-d]-psi[i,-d])))
          var1=var1+midsq[i,d,,]
        }
        if(L[i]==1){
          mean2=mean2+mid[i,d,]*((mid2[i,d]-mu2[d]-psi[i,d])*Sigma2inv[d, d]+sum(Sigma2inv[-d, d]*(Y[i,-d]-mu2[-d]-psi[i,-d])))
          var2=var2+midsq[i,d,,]
        }
      }
      Cov=Sigma1inv[d, d]*var1+Sigma2inv[d, d]*var2+Sigma.prior.r.inv
      if(!is.positive.definite(Cov)){
        Cov=nearPD(Cov)$'mat'
      }
      mean.final=mean.prior-mean1-mean2
      Cov=chol2inv(chol(Cov))
      if(!is.positive.definite(Cov)){
        Cov=nearPD(Cov)$'mat'
      }
      result = HMC_exact(F=F, g=g, M = Cov, mu_r=as.numeric(Cov%*%matrix(mean.final)), cov=TRUE, L=1, initial_X=theta[, d])
      theta[, d] = as.numeric(result$Xs)
      Y[,d]=(B[,d,]+B.star[,d,]%*%W)%*%theta[,d]+B.star[,d,]%*%matrix(q)
      
    }
    L.iter[iter,]=L
    theta.iter[iter,,]=theta
    
    # print(dim(as.matrix(kronecker(Sigma1inv, crossprod(zs[L==0,])))))
    # Sample Psi
    Sigma.psi = chol2inv(chol( sigma.psi^(-2)*diag(p*ncov) + 
                                 as.matrix(kronecker(Sigma1inv, crossprod(zs[L==0,]))) + 
                                 as.matrix(kronecker(Sigma2inv, crossprod(zs[L==1,]))) ))
    # print("check")
    Sigma.psi = as.matrix(nearPD(Sigma.psi)$`mat`)
    mu.psi = Sigma.psi %*% (c(Sigma1inv%*%t(sweep(Y[L==0,], 2, mu1, "-"))%*%zs[L==0,]) 
                            + c(Sigma2inv%*%t(sweep(Y[L==1,], 2, mu2, "-"))%*%zs[L==1,]))
    
    Psi.iter[iter, ] = psi_vec = chol(Sigma.psi) %*% rnorm(p*ncov) + mu.psi
    Psi = matrix(psi_vec, p, ncov)
    
    # print(Psi)
    sigma.psi = sqrt(1/rgamma(1, 1+p*ncov/2, 1+sum(Psi^2)/2))
    # random walk metropolis update for sigma.psi
    # log_lik = function(Y, )
    
    Sigma1inv.iter[iter,,]=Sigma1inv=rWishart(1, p+2+n1, chol2inv(chol(cov(Y[L==0,])*(n1-1)+diag(p))))[,,1]
    Sigma1.iter[iter,,]=  Sigma1=chol2inv(chol(Sigma1inv))
    mu1.iter[iter,]=  mu1=rmvnorm(1, colMeans(Y[L==0,]), Sigma1/n1)
    Sigma2inv.iter[iter,,]=Sigma2inv=rWishart(1, p+2+n2, chol2inv(chol(cov(Y[L==1,])*(n2-1)+diag(p))))[,,1]
    Sigma2.iter[iter,,]=  Sigma2=chol2inv(chol(Sigma2inv))
    mu2.iter[iter,]=  mu2=rmvnorm(1, colMeans(Y[L==1,]), Sigma2/n2)
    lambda.iter[iter] = lambda = rbeta(1, 1+n1, 1+n2)
    
    
    for(i in 1:n){
      if(L.org[i]==2){
        log_probs = c(dmvnorm(Y[i,]-zs[i,, drop=FALSE]%*%t(Psi), mu1, Sigma1, log=T), dmvnorm(Y[i,]-zs[i,, drop=FALSE]%*%t(Psi), mu2, Sigma2, log=T))
        log_probs = log_probs + log(c(lambda, 1-lambda))
        prob = exp(log_probs - logSumExp(log_probs))[2]
        L[i] = rbinom(1, 1, prob)
      }
    }
    
    n1=sum(L==0)
    n2=sum(L==1)
    if(iter%%1000==0){
      cat("\n Iterations:==========>", iter)
      cat("\n No. of misclassifications := ", sum(L!=L.true))
    }
  }
  theta.final=apply(theta.iter[(niter/5):niter,,], 2:3, mean)
  mu1.final=apply(mu1.iter[(niter/5):niter,], 2, mean)
  Sigma1.final=apply(Sigma1.iter[(niter/5):niter,,], 2:3, mean)
  mu2.final=apply(mu2.iter[(niter/5):niter,], 2, mean)
  Sigma2.final=apply(Sigma2.iter[(niter/5):niter,,], 2:3, mean)
  lambda.final=mean(lambda.iter[(niter/5):niter])
  Psi.final=matrix(apply(Psi.iter[(niter/5):niter,], 2, mean), p, ncov)
  
  ############ testing ###########
  n1.train = sum(Ltrain.true==0)
  n2.train = sum(Ltrain.true==1)
  
  B=array(dim=c(nrow(train), p, J1-2)) 
  B.star=array(dim=c(nrow(train), p, 2))
  Y.train=matrix(,N-n,p)
  for(d in 1:p){
    middle=splineDesign(knots=knotvec, x = train[,d])
    B[,d,] = middle[,-c(index1, index2)]
    B.star[,d,] = middle[,c(index1, index2)]
    Y.train[,d]=(B[,d,]+B.star[,d,]%*%W)%*%theta.final[,d]+B.star[,d,]%*%matrix(q)
  }
  p1=dmvnorm(Y.train-train.zs%*%t(Psi.final), mu1.final, Sigma1.final)
  p2=dmvnorm(Y.train-train.zs%*%t(Psi.final), mu2.final, Sigma2.final)
  r=(lambda.final * p1)<((1-lambda.final)*p2)
  tp = sum(r==1 & Ltrain.true==1)
  tn = sum(r==0 & Ltrain.true==0)
  fp = sum(r==1 & Ltrain.true==0)
  fn = sum(r==0 & Ltrain.true==1)
  mcc[sim.num] = (tn*tp-fn*fp) / sqrt((tp+tn)*(fn+fp)*(tn+fn)*(tp+fp))
  fpr[sim.num] = fp / n1.train
  fnr[sim.num] = fn / n2.train
  
  print(mean(fpr, na.rm=TRUE))
  print(sd(fpr, na.rm=TRUE))
  print(mean(fnr, na.rm=TRUE))
  print(sd(fnr, na.rm=TRUE))
  print(mean(mcc, na.rm=TRUE))
  print(sd(mcc, na.rm=TRUE))
  
}

# save.image(file="ionos.RData")

print(mean(fpr, na.rm=TRUE))
print(sd(fpr, na.rm=TRUE))
print(mean(fnr, na.rm=TRUE))
print(sd(fnr, na.rm=TRUE))
print(mean(mcc, na.rm=TRUE))
print(sd(mcc, na.rm=TRUE))
