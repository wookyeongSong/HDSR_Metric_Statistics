source("utlis.R")
library(ODP)
library(energy)
data_all_chorts <- read.csv("abide_aal/abide_phenotypic.csv")
cohort_names<-unique(data_all_chorts$SITE_ID)
age_min<-25;age_max<-50
metric<-met_FB
Data_control<-list()
Data_disease<-list()
for (k in 1:20) {
  load(paste0(cohort_names[k],".RData"))
  Data_temp<-Data
  for (i in 1:length(Data)) {
    if(Data[[i]]$DSM == 1){
      if(Data[[i]]$age<age_max&Data[[i]]$age>age_min){
        Data_disease[[length(Data_disease)+1]]<-Data[[i]]$mat
      }
    }
    if(Data[[i]]$DSM==0){
      if(Data[[i]]$age<age_max&Data[[i]]$age>age_min){
        Data_control[[length(Data_control)+1]]<-Data[[i]]$mat
      }
    }
  }
}
n1<-length(Data_control)
P_test_L<-c(Data_control,Data_disease)
P_dist_mat <- GetPairDist( data = P_test_L, distfun = metric )
set.seed(2)
prflHomTest( distmat = P_dist_mat, n = n1 )
eqdist.etest(P_dist_mat, sizes=c(n1,dim(P_dist_mat)[1] -n1), distance = T, method="discoF", R=999)





