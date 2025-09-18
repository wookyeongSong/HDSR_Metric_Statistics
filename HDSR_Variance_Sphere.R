###############################################################################
##### Title: Non-Euclidean Data Analysis With Metric Statistics (HDSR) ########
##### Description: Variance Analysis of Spherical Simulations #################
##### Manuscript reference : Section 3 ########################################
##### Figure reference : Figure 3 ############################################
##### Date: 09/17/25 ##########################################################
##### Author: Wookyeong Song ##################################################
###############################################################################

# load main functions
source("~/Desktop/WK_UCD/Research/Harvard Data Review/HDSR_Variance_mainfunctions.R")


## Set Hyperparameters
rg = pi/4
k = 100
iter = 200

## Frechet and Metric Variances
Vf = c()
Vm = c()

Fmean_df = matrix(0, nrow = iter, ncol = 3)

for(l in 1:iter){
  set.seed( l)
  z <- runif(k,min=-sin(rg),max =sin(rg))          # uniform on [0, 1]
  theta <- 4*rg*runif(k)    # uniform on [0, 2pi]
  x <- cos(theta)*sqrt(1-z^2)  # based on angle
  y <- sin(theta)*sqrt(1-z^2) 
  
  df = as.matrix(data.frame(x,y,z))
  
  df.s = matrix(0,nrow=k,ncol=2)
  for(i in 1:k){
    df.s[i,] = Trans.sph(df[i,])
  }
  
  mu.fr = Trans.Euclid(IntrinsicMean(df.s))
  Fmean_df[l,] = mu.fr
  
  ## Frechet Variance
  Vf = c(Vf, mean(acos(df %*% mu.fr)^2)) 
  
  ## Metric Variance
  dist = 0
  for(i in 1:k){
    dist = dist +sum(acos(df[-(1:i),] %*% df[i,])^2)
  }
  Vm = c(Vm, dist/(k*k))
  
}


# -------------------------------------------------------------------------
# Figure 3
# -------------------------------------------------------------------------
par(mfrow=c(1,1),mar=c(5,6,3,2))
# Create a scatterplot
plot(Vm, Vf, 
     xlab = expression("Metric Variance"~hat(V)[M]), 
     ylab = expression("Fréchet Variance"~hat(V)[F]), 
     main = "", 
     pch = 16,
     col = 'cadetblue',
     xlim = range(c(Vm, Vf)), 
     ylim = range(c(Vm, Vf)),cex.lab=2,cex.axis = 2, cex.main = 2)

# Add the line y = x
abline(a = 0, b = 1, col = "coral3", lty = 2)  # y = x line in red dashed

