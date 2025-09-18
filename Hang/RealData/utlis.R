Get_Inverse <- function(u, y) {
  min_y <- range(y)[1]
  max_y <- range(y)[2]
  min_u <- range(u)[1]
  max_u <- range(u)[2]
  temp_x <- seq(min_y, max_y, length = length(y))
  temp_y <- numeric(length = length(u))
  temp_y <- approx(y, u, temp_x, rule = 2,ties=min)$y
  return(list(x = temp_x, y = temp_y))
}
met_R1<-function(a,b){return(abs(a-b))}
met_sphere<-function(a,b){
  return(acos(min(t(a)%*%b,1) ) )
}
met_FB<-function(A,B){
  return(sqrt(sum(diag(((A-B)%*%t(A-B))) ) ))
}
Dist_profile<-function(w,x,t,Training_set,metric,h){
  Index<-which((Training_set$x>=x-h)&(Training_set$x<=x+h) )
  S <- unlist(lapply(Training_set$y[Index],  function(x) metric(x,w)  )) 
  W<-epanechnikov_kernel((unlist(Training_set$x[Index]) -x)/h )
  Out <- sapply(t, function(x) sum((S <= x)*W)/sum(W))
  return(Out)
}
R_n<-function(y,x,Training_set,metric,F_mean,Q_train,h,workgrid_p){
  Index<-which((Training_set$x>=x-h)&(Training_set$x<=x+h) )
  W<-epanechnikov_kernel((unlist(Training_set$x[Index]) -x)/h )
  Q_train_in<-Q_train[Index,]
  Q_i<-Get_Inverse(workgrid_p,Dist_profile(y,x,workgrid_p,Training_set,metric,h))$y
  R_i<-sum( abs(Q_train_in-rep(1,dim(Q_train_in)[1])%*%t(Q_i))*W )/(length(Q_i))
  return(expit(R_i/sum(W)) )
}
H<-function(z,x,Training_set,Ranks_on_train,h){
  Index<-which((Training_set$x>=x-h)&(Training_set$x<=x+h) )
  W<-epanechnikov_kernel((unlist(Training_set$x[Index]) -x)/h )
  return(sum((Ranks_on_train[Index]<=z)*W)/sum(W) )
}
Get_conf_set<-function(h=0.1,metric=met_sphere,R_depth,Search_grid,Training_set_L,
                       Training_set,workgrid_p,Calibration_L,Calibration_set,alpha=0.9){
  dist_profile <- t(matrix(unlist(lapply(Training_set_L, function(x) Dist_profile(x$y,x$x, workgrid_p, Training_set,metric,h))),
                           length(workgrid_p) ,length(Training_set_L) ) )
  inv_dist_profile <- apply(dist_profile, 1, function(x) Get_Inverse(workgrid_p, x)$y)
  Q_train <- t(inv_dist_profile)
  Ranks_on_train<-unlist(lapply(Training_set_L,function(x) R_depth(x$y,x$x,Training_set,metric,F_mean=NULL,Q_train,h,workgrid_p)))
  Ranks_on_cal<-unlist(lapply(Calibration_L,function(x) R_depth(x$y,x$x,Training_set,metric,F_mean=NULL,Q_train,h,workgrid_p)))
  S_cal_list<-list()
  for (i in 1:length(Calibration_L)) {
    S_cal_list[[i]]<-list(z=Ranks_on_cal[i],x=Calibration_set$x[[i]])
  }
  Scores_on_cal<-unlist(lapply(S_cal_list, function(x) H(x$z,x$x,Training_set,Ranks_on_train,h) ))
  q_H<-quantile(Scores_on_cal,alpha*(1+1/length(Calibration_L)),names = FALSE)
  ##serach grid
  Ranks_on_SG<-unlist(lapply(Search_grid,function(x) R_depth(x$y,x$x,Training_set,metric,F_mean=NULL,Q_train,h,workgrid_p)))
  S_SG_list<-list()
  for (i in 1:length(Search_grid)) {
    S_SG_list[[i]]<-list(z=Ranks_on_SG[i],x=Search_grid[[i]]$x)
  }
  Scores_on_SG<-unlist(lapply(S_SG_list, function(x) H(x$z,x$x,Training_set,Ranks_on_train,h) ))
  C_set<-Search_grid[which(Scores_on_SG<=q_H)] 
  return(C_set)
}
f<-function(x){return((x-1)^2*(x+1) )}
g<-function(x){
  out<-rep(0,length(x))
  id<-which(x>=0)
  out[id]<-sqrt(x[id])
  return(2*out)
}
f_sphere<-function(x){
  out<-matrix(0,length(x),3)
  for (i in which(x<=0)) {
    out[i,]<-c(sin(pi*x[i]/2),cos(pi*x[i]/2),0)
  }
  return(out)
}
g_sphere<-function(x){
  out<-matrix(0,length(x),3)
  for (i in which(x>0)) {
    if(rbinom(1,1,0.5)){
      out[i,]<-c(0,cos(pi*x[i]/2),sin(pi*x[i]/2))
    }else{
      out[i,]<-c(0,cos(pi*x[i]/2),-sin(pi*x[i]/2))
    }
  }
  return(out)
}




belongs_to_set <- function(datapoints, input_set, gap = 1) {
  # Calculate the differences between consecutive elements in the input_set
  diffs <- diff(input_set)
  
  # Find the indices where the difference is greater than the gap
  interval_indices <- c(0, which(diffs > gap+sqrt(.Machine$double.eps)), length(input_set))
  
  # Define a function to check if a single datapoint belongs to any interval
  check_datapoint <- function(datapoint, input_set, interval_indices) {
    start <- input_set[interval_indices[-length(interval_indices)] + 1]
    end <- input_set[interval_indices[-1]]
    any(datapoint >= start-sqrt(.Machine$double.eps) & datapoint <= end+sqrt(.Machine$double.eps))
  }
  
  # Apply the check_datapoint function to each element of datapoints
  belongs <- sapply(datapoints, check_datapoint, input_set, interval_indices)
  
  return(belongs)
}
met_SPD <- function(P,Q) {
  return(norm(logm(solve(sqrtm(P))%*%Q%*%solve(sqrtm(P))) ,"F") )
}
mu_SPD<-function(x){
  if(x<=0){
    if(rbinom(1,1,0.5)){
      return(matrix(c(1,(1+x)/2,(1+x)/2,1),2,2))
    }else{
      return(matrix(c(1,(-x-1)/2,(-x-1)/2,1),2,2))
    }
  }else{
    if(rbinom(1,1,0.5)){
      return(matrix(c(1+x,(1+x)/2,(1+x)/2,1),2,2))
    }else{
      return(matrix(c(1,(-x-1)/2,(-x-1)/2,1+x),2,2))
    }
  }
}
geodesic_distance_SPD <- function(A, B) {
  # Ensure A and B are SPD matrices (symmetric and positive definite)
  if (!all(eigen(A)$values > 0) || !all(eigen(B)$values > 0)) {
    stop("Both matrices must be symmetric positive definite (SPD).")
  }
  
  # Calculate A^(-1/2)
  A_inv_sqrt <- sqrtm(solve(A))
  
  # Compute the matrix logarithm of A^(-1/2) * B * A^(-1/2)
  log_matrix <- logm(A_inv_sqrt %*% B %*% A_inv_sqrt)
  
  # Compute the Frobenius norm of the log matrix
  distance <- norm(log_matrix, type = "F")
  
  return(distance)
}

approximate_matrix <- function(mat) {
  # Calculate eigenvalues and eigenvectors
  eig <- eigen(mat)
  
  # Filter eigenvalues greater than 10^-5
  significant_indices <- which(abs(eig$values) > 1e-5)
  
  # Keep only the significant eigenvalues and corresponding eigenvectors
  significant_values <- eig$values[significant_indices]
  significant_vectors <- eig$vectors[, significant_indices]
  
  # Reconstruct the matrix approximation
  approx_mat <- significant_vectors %*% diag(significant_values) %*% t(significant_vectors)
  
  return(approx_mat)
}



epanechnikov_kernel <- function(u) {
  ifelse(abs(u) <= 1, 3/4 * (1 - u^2), 0)
}





