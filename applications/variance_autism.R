###############################################################################
##### Title: Non-Euclidean Data Analysis With Metric Statistics (HDSR) ########
##### Description: Variance Analysis of Autism Brain Imaging Data #############
##### Manuscript reference : Section 5 ########################################
##### Figure reference : Figure 10 ############################################
##### Date: 09/17/25 ##########################################################
##### Author: Wookyeong Song ##################################################
###############################################################################

# load main functions
source("utlis.R")

# -------------------------------------------------------------------------
# Data Preprocessing 
# -------------------------------------------------------------------------
## Predictors
covariates_path = "abide_aal/abide_phenotypic.csv"
aal_pred = read.csv(covariates_path)
sub_ids = aal_pred[,'subject']

## NYU Responses
folder_path <- "abide_aal/aal/"  # Replace with the path to your folder
aal_resp = list()
for (i in 1:length(sub_ids)){
  aal_resp[[i]] = read.csv(paste0(folder_path, sub_ids[i], ".csv"), header = FALSE)
}

## Remove subjects having correlations matrix with NAs: 61, 72, 75, 135, 151
na_cor_ids = c()
for(i in 1:length(aal_resp)){
  if(sum(is.na(cor(aal_resp[[i]]))) != 0){na_cor_ids = c(na_cor_ids, i)}
}

## Pearson correlation to network with trucation alpha = 0.15
mPred_aal = data.frame()
lResp_aal = list()
for (i in 1:length(aal_resp)){
  
  if (i %in% na_cor_ids) {next}
  
  if (25 <= aal_pred[i,'AGE_AT_SCAN'] & aal_pred[i,'AGE_AT_SCAN'] <= 50) {
    
    mPred_aal = rbind(mPred_aal, aal_pred[i,])
    lResp_aal = c(lResp_aal, list(cor(aal_resp[[i]])))
    
  }
  
}

# -------------------------------------------------------------------------
# Figure 10 (Left)
# -------------------------------------------------------------------------
## Case Index
case_idx = which(mPred_aal[,'DX_GROUP'] == 1)
case_mPred_aal = mPred_aal[case_idx,]
case_lResp_aal = lResp_aal[case_idx]
case_k = length(case_lResp_aal)

## Control Index
control_idx = which(mPred_aal[,'DX_GROUP'] == 2)
control_mPred_aal = mPred_aal[control_idx,]
control_lResp_aal = lResp_aal[control_idx]
control_k = length(control_lResp_aal)

## Case distance matrix
case_dmat = matrix(0, nrow = case_k, ncol = case_k)
for(i in 1:case_k){
  
  for(j in 1:case_k){
    
    if (i == j){next}
    
    case_dmat[i,j] = sqrt(sum((case_lResp_aal[[i]] - case_lResp_aal[[j]])^2))
    
  }
}

## Control distance matrix
control_dmat = matrix(0, nrow = control_k, ncol = control_k)
for(i in 1:control_k){
  
  for(j in 1:control_k){
    
    if (i == j){next}
    control_dmat[i,j] = sqrt(sum((control_lResp_aal[[i]] - control_lResp_aal[[j]])^2))
    
  }
}

## Case metric variance and its asymptotic variance
case_Vm = sum(case_dmat^2) / (2*case_k*case_k)
case_var = asym.var(case_dmat)

## Control metric variance
control_Vm = sum(control_dmat^2) / (2*control_k*control_k)
control_var = asym.var(control_dmat)

## Plot Figure 10 (Left)
names.abb = c("Case - Autism", "Control - Normal")
par(mar=c(5,5,3,2))
plotCI(x= 1:2, y= c(case_var$Vm, control_var$Vm), uiw = c(1.96*case_var$asym.sd, 1.96*control_var$asym.sd) , liw = c(1.96*case_var$asym.sd, 1.96*control_var$asym.sd) ,ylab = "Fréchet Variance",axes=FALSE,col = c("cadetblue", "coral3"),scol=c("black", "black"),xlab = "",xlim = c(0.5,2.5),cex=2,pch=19,font.axis=2, cex.axis=2, cex.lab=2, cex.main = 2.5,main="")
axis(side=2,cex.axis=2)         ## add default y-axis (ticks+labels)
axis(side=1,at=c(1,2),  ## add custom x-axis
     label=names.abb,font.axis=1, cex.axis=2, cex.lab=2)
box(bty="l")


# -------------------------------------------------------------------------
# Figure 10 (Right)
# -------------------------------------------------------------------------
## Variance Trends of Case and Control
case_var_overtime = c()
control_var_overtime = c()

range = 5
age = 25
for (j in 1:15){
  
  mPred_subset = data.frame()
  lResp_subset = list()
  for (i in 1:nrow(mPred_aal)){
    
    if ( (age - range) <= mPred_aal[i,'AGE_AT_SCAN'] & mPred_aal[i,'AGE_AT_SCAN'] < (age + range) ) {
      
      mPred_subset = rbind(mPred_subset, mPred_aal[i,])
      lResp_subset = c(lResp_subset, list(lResp_aal[[i]]))
      
    }
  }

  res = variance.object(mPred_subset, lResp_subset)
  case_var_overtime = c(case_var_overtime, res$case_var$Vm)
  control_var_overtime = c(control_var_overtime, res$control_var$Vm)
  
  age = age + 1
}

## Plot Figure 10 (Right)
par(mar=c(5,5,3,2))
plot(case_var_overtime, type = "l", col = 'cadetblue', ylim = c(400, 600),xaxt = 'n', xlab = "Age", ylab = "Fréchet variance", lwd = 2,cex.lab=2,cex.axis=2,cex.main=2.5)
axis(side=1,at=seq(1, 15, length.out = 6),label= c("25","28", "31", "34", "37", "40+"), cex.axis=2, cex.lab=2)
lines(control_var_overtime, type = "l", col = 'coral3', lwd = 2)
legend("topleft", legend = c("Case (Autism)", "Control (Normal)"), col = c("cadetblue", "coral3"),lwd = 2,cex=2)


# -------------------------------------------------------------------------
# Comparison of correlation space with different metrics 1) Affine-invariant, negatively curved, 2) Bures-Wasserstein, positively curved, and 3) Frobenius, flat.
# -------------------------------------------------------------------------
cov.array = simplify2array(lResp_aal)
new.cov.array = array(0, dim = dim(cov.array))
for(i in 1:dim(cov.array)[3]){
  eigenDecom = eigen(cov.array[,,i])
  U = eigenDecom$vectors
  Lambda = pmax(eigenDecom$values, 1e-6)
  qNew = U %*% diag(Lambda) %*% t(U)
  new.cov.array[,,i] <- (qNew + t(qNew))/2
}

## Affine-Invariant Riemannian (AIR) metric
cov.curv.fit.AIRM = cov.curv.est(new.cov.array, metric = "AIRM")
cov.curv.fit.AIRM$Vf
cov.curv.fit.AIRM$Vm

## Procrustes size-and-shape (PSS) metric (a.k.a Bures-Wasserstein (BW) metric)
cov.curv.fit.PSS = cov.curv.est(new.cov.array, metric = "Procrustes.SS")
cov.curv.fit.PSS$Vf
cov.curv.fit.PSS$Vm

## Frobenius metric
cov.curv.fit.Frob = cov.curv.est(new.cov.array, metric = "Euclidean")
cov.curv.fit.Frob$Vf
cov.curv.fit.Frob$Vm

