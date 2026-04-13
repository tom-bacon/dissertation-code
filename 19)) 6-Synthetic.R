# This script generates a synthetic multi-class curve dataset
# from curvature functions, applies geometric transformations,
# and visualises representative curves and their curvature.

library(signal)
library(ggplot2)
library(patchwork)

set.seed(1)
n<-200
x<-seq(-5,5,length.out=n)
sc<-c(6,8,10)

# Gaussian bump
g<-function(x,m,s) exp(-(x-m)^2/(2*s^2))

# Generate centres with small randomness
gen_centres<-function(K,c0,d){
  base<-seq(-(K-1)/2,(K-1)/2,length.out=K)*d+c0
  jitter<-runif(K,-0.15*d,0.15*d)
  base+jitter
}

# Build curvature as sum of bumps
cmulti<-function(x,a,K){
  s<-runif(1,0.2,0.35)
  c0<-runif(1,-0.5,0.5)
  d<-(2.0+0.6*K)*s
  ms<-gen_centres(K,c0,d)
  if(K==1) ms<-ms+runif(1,-0.5,0.5)
  s_i<-rep(s,K)
  if(K==1){
    signs<-sample(c(-1,1),1)
  } else {
    signs<-rep(c(1,-1),length.out=K)
  }
  w<-runif(K,0.6,1)
  w<-w/sum(w)
  w<-w*signs
  k<-rep(0,length(x))
  for(i in 1:K){
    k<-k+w[i]*g(x,ms[i],s_i[i])
  }
  dx<-x[2]-x[1]
  total_curv<-sum(abs(k))*dx
  target<-runif(1,2.5,4.5)
  k*(target/(total_curv+1e-8))
}

# Reconstruct curve from curvature
rec<-function(k,x){
  ds<-diff(x)[1]
  th<-cumsum(k)*ds
  X<-cumsum(cos(th))*ds
  Y<-cumsum(sin(th))*ds
  X<-X-mean(X)
  Y<-Y-mean(Y)
  list(X=X,Y=Y)
}

# Random transformations
rot<-function(z,a){
  R<-matrix(c(cos(a),-sin(a),sin(a),cos(a)),2,2)
  t(R%*%t(z))
}

tra<-function(z){
  z[,1]<-z[,1]+runif(1,-5,5)
  z[,2]<-z[,2]+runif(1,-5,5)
  z
}

sca<-function(z){
  z*runif(1,0.1,10)
}

# Simple reparameterisation
repf<-function(z){
  dx<-diff(z[,1]);dy<-diff(z[,2])
  s<-c(0,cumsum(sqrt(dx^2+dy^2)))
  s<-s/max(s)
  w<-s^runif(1,0.2,5)
  w<-(w-min(w))/(max(w)-min(w))
  x<-approx(s,z[,1],xout=w)$y
  y<-approx(s,z[,2],xout=w)$y
  cbind(x,y)
}

class_K<-c(1,2,3,4,5,6)
n_class<-6
n_per_class<-100

curves<-vector("list",n_class*n_per_class)
labels<-rep(1:n_class,each=n_per_class)

idx<-1

# Generate base curves
for(cl in 1:n_class){
  for(i in 1:n_per_class){
    a<-sample(sc,1)
    K<-class_K[cl]
    k<-cmulti(x,a,K)
    r<-rec(k,x)
    curves[[idx]]<-cbind(r$X,r$Y)
    idx<-idx+1
  }
}

# Augment dataset with transformations
aug<-vector("list",length(curves))
lab2<-labels

for(i in 1:length(curves)){
  j<-sample(1:length(curves),1)
  z<-curves[[j]]
  z<-rot(z,runif(1,0,2*pi))
  z<-sca(z)
  z<-repf(z)
  z<-tra(z)
  aug[[i]]<-z
}

curves<-c(curves,aug)
labels<-c(labels,lab2)

data<-list(curves=curves,labels=labels)
saveRDS(data,"curve_dataset6class.rds")

# Curvature function for visualisation
curv<-function(x,y){
  n<-length(x)
  t<-seq(0,1,length.out=n)
  w<-max(7,floor(n/10))
  if(w%%2==0) w<-w+1
  xs<-sgolayfilt(x,4,w)
  ys<-sgolayfilt(y,4,w)
  dx<-sgolayfilt(xs,4,w,m=1)
  dy<-sgolayfilt(ys,4,w,m=1)
  ddx<-sgolayfilt(xs,4,w,m=2)
  ddy<-sgolayfilt(ys,4,w,m=2)
  k<-(dx*ddy-ddx*dy)/((dx^2+dy^2)^(3/2)+1e-6)
  list(t=t,k=k)
}

# Larger, cleaner ggplot style
theme_big<-function(){
  theme_minimal(base_size=18)+
    theme(panel.grid=element_blank(),axis.line=element_line(color="black"),
      axis.ticks=element_line(color="black"),
      axis.title=element_text(face="bold"),
      axis.text=element_text(color="black"),
      plot.title=element_text(face="bold",hjust=0.5))
}

plots<-list()

# Build curvature plots for each class
for(cl in 1:6){
  idx<-which(labels==cl)[1]
  z<-curves[[idx]]
  r<-curv(z[,1],z[,2])
  df<-data.frame(t=r$t,k=r$k)
  plots[[cl]]<-ggplot(df,aes(t,k))+geom_line(color="blue",linewidth=1.4)+
    labs(title=paste("Class",cl),x="t",y=expression(kappa))+theme_big()
}

# Combine into grid
final_plot<-(plots[[1]]|plots[[2]]|plots[[3]])/
  (plots[[4]]|plots[[5]]|plots[[6]])

final_plot
