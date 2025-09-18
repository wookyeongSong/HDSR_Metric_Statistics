library(expm)
Get_Inverse <- function(u, y,xout) {
  temp_y <- numeric(length = length(u))
  temp_y <- approx(y, u, xout, rule = 2,ties=min)$y
  return(list(x = xout, y = temp_y))
}
met_R2<-function(a,b){return(norm(a-b,type = "2"))}
met_R1<-function(a,b){return(abs(a-b))}
met_sphere<-function(a,b){
  return(acos(min(t(a)%*%b,1) ) )
}
met_SPD <- function(P,Q) {
  return(norm(logm(solve(sqrtm(P))%*%Q%*%solve(sqrtm(P))) ,"F") )
}
Dist_profile<-function(w,t,Training_set,metric){
  # S <- apply(Training_set, 1, function(x) metric(x,w)  )
  S <- unlist(lapply(Training_set,  function(x) metric(x,w)  )) 
  n <- length(Training_set)
  Out <- sapply(t, function(x) sum(S <= x)/n)
  return(Out)
}
Depth_dc<-function(x,Training_set,metric,Q_train,workgrid_p,workgrid_q,tau=1){
  w<-exp(-tau*workgrid_q)
  Q_i<-Get_Inverse(workgrid_p,Dist_profile(x,workgrid_p,Training_set,metric),workgrid_q)$y
  R_i<-sum( (Q_train-rep(1,dim(Q_train)[1])%*%t(Q_i) )%*%(w) )/length(Q_i)
  return(expit(R_i/dim(Q_train)[1]))
}
Depth_abs<-function(x,Training_set,metric,Q_train,workgrid_p,workgrid_q,tau=1){
  w<-exp(-tau*workgrid_q)
  Q_i<-Get_Inverse(workgrid_p,Dist_profile(x,workgrid_p,Training_set,metric),workgrid_q)$y
  R_i<-sum( abs(Q_train-rep(1,dim(Q_train)[1])%*%t(Q_i) )%*%(w) )/length(Q_i)
  return(expit(R_i/dim(Q_train)[1]))
}
expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}