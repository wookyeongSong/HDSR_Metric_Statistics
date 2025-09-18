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

