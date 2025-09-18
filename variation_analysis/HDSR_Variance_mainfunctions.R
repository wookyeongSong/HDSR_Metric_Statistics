###############################################################################
##### Title: Non-Euclidean Data Analysis With Metric Statistics (HDSR) ########
##### Description: Main Functions Necessary for Variance Analysis #############
##### Date: 09/17/25 ##########################################################
##### Author: Wookyeong Song ##################################################
###############################################################################
library(CovTools)
library(matrixcalc)
library(Matrix)
library(abind)

asym.var = function(dmat){
  
  k = dim(dmat)[1]
  
  Vm = sum(dmat^2) / (2*k*(k-1))
  
  asym.dist=0
  for(i in 1:k){
    asym.dist = asym.dist +(sum(dmat[i,]^2)/(k-1))^2
  }
  
  asym.Vm = asym.dist/k - (2*Vm)^2
  
  asym.sd = sqrt(asym.Vm/k)
  
  return(list(Vm = Vm, asym.Vm = asym.Vm, asym.sd = asym.sd))
  
}


variance.object = function(mPred_aal, lResp_aal){
  
  if(length(lResp_aal) == 0){
    stop("No observation")
  }
  
  case_idx = which(mPred_aal[,'DX_GROUP'] == 1)
  case_mPred_aal = mPred_aal[case_idx,]
  case_lResp_aal = lResp_aal[case_idx]
  
  control_idx = which(mPred_aal[,'DX_GROUP'] == 2)
  control_mPred_aal = mPred_aal[control_idx,]
  control_lResp_aal = lResp_aal[control_idx]
  
  case_k = length(case_lResp_aal)
  control_k = length(control_lResp_aal)
  
  
  case_dmat = matrix(0, nrow = case_k, ncol = case_k)
  for(i in 1:case_k){
    
    for(j in 1:case_k){
      
      if (i == j){next}
      
      case_dmat[i,j] = sqrt(sum((case_lResp_aal[[i]] - case_lResp_aal[[j]])^2))
      
    }
  }
  
  control_dmat = matrix(0, nrow = control_k, ncol = control_k)
  for(i in 1:control_k){
    
    for(j in 1:control_k){
      
      if (i == j){next}
      control_dmat[i,j] = sqrt(sum((control_lResp_aal[[i]] - control_lResp_aal[[j]])^2))
      
    }
  }
  
  
  case_var = asym.var(case_dmat)
  control_var = asym.var(control_dmat)
  
  return(list(case_var = case_var, control_var = control_var))
}

# Main function: Estimate the variances, curvature and test statistics when we have SPD matrices as inputs
cov.curv.est = function(cov.array, metric = NULL) {
  
  if (!metric %in% c("AIRM","Cholesky","Euclidean","LERM","Procrustes.SS","RootEuclidean")){
    stop("metric choice not supported.")
  }
  
  # k: number of observations
  k = dim(cov.array)[3]
  
  # Fmean: Frechet mean
  Fmean = CovMean(cov.array, method=metric)
  
  samples = abind(cov.array,Fmean)
  
  dist.mat = CovDist(samples, method=metric)
  
  # Vf: Frechet variance estimate
  # Vm: Metric variance estimate
  Vf = sum(dist.mat[k+1,1:k]^2)/k
  Vm = sum(dist.mat[1:k,1:k]^2)/(2*k*k)
  
  # sigma.Vf: asymptotic variance of Frechet variance
  # sigma.Vm: asymptotic variance of Metric variance
  # sigma.Vmf: asymptotic variance of cross product of Vf and Vm
  sigma.Vf = mean(dist.mat[k+1,1:k]^4) - Vf^2
  
  asym.dist = 0
  for(i in 1:k){
    asym.dist = asym.dist +(sum(dist.mat[i,1:k]^2)/(k-1))^2
  }
  sigma.Vm = asym.dist/k - (2*Vm)^2
  
  coasym.dist = 0
  for(i in 1:k){
    coasym.dist = coasym.dist + sum(dist.mat[i,1:k]^2)/(k-1)*(dist.mat[k+1,i]^2)
  }
  sigma.Vmf = coasym.dist/k - 2*Vf*Vm
  
  cov = matrix(data = c(sigma.Vm,sigma.Vmf,sigma.Vmf,sigma.Vf),nrow=2)
  cov.normalized = matrix(data = c(sigma.Vm,sigma.Vmf,sigma.Vmf,sigma.Vf)/k,nrow=2)
  
  rho = Vf/Vm - 1
  
  a = c(-Vf/(Vm^2),1/Vm)
  
  sd = sqrt(a %*% cov %*% a)
  
  test.stat = sqrt(k)*rho/sd
  
  pval = 2*min(1-pnorm(test.stat), pnorm(test.stat))
  
  return(list(Fmean = Fmean, Vf = Vf, Vm = Vm, cov = cov, cov.normalized = cov.normalized, test.stat= test.stat, pval = pval, rho= rho, sd = sd))
  
}
