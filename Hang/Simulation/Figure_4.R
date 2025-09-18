library(ggplot2)
library(MASS)
library(cowplot)
library(KernSmooth)
library(mvtnorm)
source('utils.R')
set.seed(2024)
n<-40000;alpha=0.05
df <- 5  
mu <- c(0, 0)  # Mean vector
sigma <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)  # Scale matrix (covariance structure)
Training_set<-list()
for (i in 1:n) {
  Training_set[[i]]<-rmvt(n = 1, sigma = sigma, df = df, delta = mu)
}
df_train<-data.frame(x=rep(0,n),y=rep(0,n))
for(i in 1:200){
  df_train[i,]<-Training_set[[i]]
}
far_p<-c(4,0)
mid_p<-c(2,0)
cen_p<-c(0,0)
bro_p<-c(2,2)
org_p<-c(-4,4)
### Fig1 1st row
p_u_combined<-ggplot(data = df_train,aes(x=x,y=y))+geom_point(size=0.3,alpha=0.4)+xlim(-5,5)+ylim(-5,5)+
  geom_point(data=data.frame(x=cen_p[1],y=cen_p[2]),aes(x=x,y=y),color="red",size=2)+
  # geom_point(data=data.frame(x=bro_p[1],y=bro_p[2]),aes(x=x,y=y),color="brown",size=2)+
  # geom_point(data=data.frame(x=org_p[1],y=org_p[2]),aes(x=x,y=y),color="orange",size=2)+
  geom_point(data=data.frame(x=mid_p[1],y=mid_p[2]),aes(x=x,y=y),color="purple",size=2)+
  geom_point(data=data.frame(x=far_p[1],y=far_p[2]),aes(x=x,y=y),color="blue",size=2)+theme_minimal()+
  theme(
    axis.title.x = element_text(size = 14),  # Font size for x-axis label
    axis.title.y = element_text(size = 14),  # Font size for y-axis label
    axis.text.x = element_text(size = 10),   # Font size for x-axis tick labels
    axis.text.y = element_text(size = 10)    # Font size for y-axis tick labels
  )

### Fig1 2nd row
workgrid_p<-seq(0,10,length=501)
workgrid_q<-c(0,sort(rbeta(499,0.5,0.5)) ,1)
d_profile_c<-Dist_profile(cen_p, workgrid_p, Training_set,met_R2)
d_profile_m<-Dist_profile(mid_p, workgrid_p, Training_set,met_R2)
d_profile_f<-Dist_profile(far_p , workgrid_p, Training_set,met_R2)

c_cdf_smth<-locpoly(workgrid_p , d_profile_c, bandwidth=0.8, degree=1)
m_cdf_smth<-locpoly(workgrid_p , d_profile_m, bandwidth=0.8, degree=1)
f_cdf_smth<-locpoly(workgrid_p , d_profile_f, bandwidth=0.8, degree=1)


c_cdf_smth$y[c_cdf_smth$y<0]<-0
m_cdf_smth$y[m_cdf_smth$y<0]<-0
f_cdf_smth$y[f_cdf_smth$y<0]<-0


p_m_combined<-ggplot(data=data.frame(x=f_cdf_smth$x,y=f_cdf_smth$y),aes(x=x,y=y))+
  geom_line(color = "blue",size=1.2)+
  theme_minimal()+labs(x="r", y="")+
  geom_line(data=data.frame(x=m_cdf_smth$x,y=m_cdf_smth$y),aes(x=x,y=y),color = "purple",linetype="dashed",size=1.2)+
  geom_line(data=data.frame(x=c_cdf_smth$x,y=c_cdf_smth$y),aes(x=x,y=y),color = "red",size=1.2)+
  # geom_line(data=data.frame(x=b_cdf_smth$x,y=b_cdf_smth$y),aes(x=x,y=y),color = "brown")+
  # geom_line(data=data.frame(x=o_cdf_smth$x,y=o_cdf_smth$y),aes(x=x,y=y),color = "orange")+
  theme(
    axis.title.x = element_text(size = 14),  # Font size for x-axis label
    axis.title.y = element_text(size = 14),  # Font size for y-axis label
    axis.text.x = element_text(size = 10),   # Font size for x-axis tick labels
    axis.text.y = element_text(size = 10)    # Font size for y-axis tick labels
  )
###bottom
d_pdf_c<-diff(d_profile_c);d_pdf_m<-diff(d_profile_m);d_pdf_f<-diff(d_profile_f)
d_pdf_c<-c(0,d_pdf_c);d_pdf_m<-c(0,d_pdf_m);d_pdf_f<-c(0,d_pdf_f)
c_pdf_smth<-locpoly(workgrid_p , d_pdf_c, bandwidth=0.8, degree=1)
m_pdf_smth<-locpoly(workgrid_p , d_pdf_m, bandwidth=0.8, degree=1)
f_pdf_smth<-locpoly(workgrid_p , d_pdf_f, bandwidth=0.8, degree=1)
c_pdf_smth$y<-c_pdf_smth$y/sum(c_pdf_smth$y)*length(c_pdf_smth$x)/10
m_pdf_smth$y<-m_pdf_smth$y/sum(m_pdf_smth$y)*length(m_pdf_smth$x)/10
f_pdf_smth$y<-f_pdf_smth$y/sum(f_pdf_smth$y)*length(f_pdf_smth$x)/10
f_pdf_smth$y[f_pdf_smth$y<0]<-0
m_pdf_smth$y[m_pdf_smth$y<0]<-0


p_b_combined<-ggplot(data=data.frame(x=f_pdf_smth$x,y=f_pdf_smth$y),aes(x=x,y=y))+
  geom_line(color = "blue",size=1.2)+
  labs(x="r", y="")+
  geom_line(data=data.frame(x=m_pdf_smth$x,y=m_pdf_smth$y),aes(x=x,y=y),color = "purple",size=1.2)+
  geom_line(data=data.frame(x=c_pdf_smth$x,y=c_pdf_smth$y),aes(x=x,y=y),color = "red",size=1.2)+
  # geom_line(data=data.frame(x=b_pdf_smth$x,y=b_pdf_smth$y),aes(x=x,y=y),color = "brown")+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 14),  # Font size for x-axis label
    axis.title.y = element_text(size = 14),  # Font size for y-axis label
    axis.text.x = element_text(size = 10),   # Font size for x-axis tick labels
    axis.text.y = element_text(size = 10)    # Font size for y-axis tick labels
  )


plot_grid(p_u_combined,p_m_combined,p_b_combined,nrow = 1)

