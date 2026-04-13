# This script generates synthetic curves from curvature functions,
# applies random transformations for augmentation, and visualises
# examples of original and transformed curves.

library(ggplot2)
library(patchwork)

n<-100
x<-seq(-2,2,length.out=100)
sc<-c(10,12,14)

# Gaussian basis
g<-function(x,m,s) exp(-(x-m)^2/(2*s^2))

# Curvature type 1
f1<-function(x,a){
  s<-runif(1,0.15,0.2)
  c<-runif(1,-0.5,0.5)
  d<-runif(1,2,3)*s
  m1<-c-d/2
  m2<-c+d/2
  p<-runif(1,0.4,0.6)
  a*(p*g(x,m1,s)-(1-p)*g(x,m2,s))
}

# Curvature type 2
f2<-function(x,a){
  s<-runif(1,0.15,0.2)
  c<-runif(1,-0.5,0.5)
  d<-runif(1,3,5)*s
  m1<-c-d
  m2<-c
  m3<-c+d
  p<-runif(1,0.3,0.45)
  q<-runif(1,0.25,0.35)
  a*(p*g(x,m1,s)-q*g(x,m2,s)+(1-p-q)*g(x,m3,s))
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

# Transformations
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
  s<-runif(1,0.1,10)
  z*s
}

# Arc-length warping
warp<-function(z){
  n<-nrow(z)
  dx<-diff(z[,1])
  dy<-diff(z[,2])
  s<-c(0,cumsum(sqrt(dx^2+dy^2)))
  s<-s/max(s)
  w<-s^runif(1,0.2,5)
  w<-(w-min(w))/(max(w)-min(w))
  x<-approx(s,z[,1],xout=w)$y
  y<-approx(s,z[,2],xout=w)$y
  cbind(x,y)
}

# Generate base curves
curves<-vector("list",200)
lab<-c(rep(1,100),rep(2,100))
k1s<-list()
k2s<-list()
for(i in 1:100){
  a<-sample(sc,1)
  k<-f1(x,a)
  k1s[[i]]<-k
  r<-rec(k,x)
  curves[[i]]<-cbind(r$X,r$Y)
}
for(i in 1:100){
  a<-sample(sc,1)
  k<-f2(x,a)
  k2s[[i]]<-k
  r<-rec(k,x)
  curves[[100+i]]<-cbind(r$X,r$Y)
}

# Data augmentation
aug<-vector("list",200)
lab2<-c(rep(1,100),rep(2,100))
for(i in 1:100){
  j<-sample(1:100,1)
  z<-curves[[j]]
  z<-rot(z,runif(1,0,2*pi))
  z<-sca(z)
  z<-warp(z)
  z<-tra(z)
  aug[[i]]<-z
}
for(i in 1:100){
  j<-sample(101:200,1)
  z<-curves[[j]]
  z<-rot(z,runif(1,0,2*pi))
  z<-sca(z)
  z<-warp(z)
  z<-tra(z)
  aug[[100+i]]<-z
}
curves<-c(curves,aug)
lab<-c(lab,lab2)

# Save dataset
data<-list(curves=curves,labels=lab)
saveRDS(data,"curve_dataset.rds")

# Sample examples for plotting
i<-sample(1:100,1)
j<-sample(1:100,1)
k1<-k1s[[i]]
k2<-k2s[[j]]
r1<-rec(k1,x)
r2<-rec(k2,x)
z1<-cbind(r1$X,r1$Y)
z1<-rot(z1,runif(1,0,2*pi))
z1<-sca(z1)
z1<-warp(z1)
z1<-tra(z1)
z2<-cbind(r2$X,r2$Y)
z2<-rot(z2,runif(1,0,2*pi))
z2<-sca(z2)
z2<-warp(z2)
z2<-tra(z2)

# Data for plotting
d1<-data.frame(s=x,k=k1)
d2<-data.frame(s=x,k=k2)
d3<-data.frame(x=r1$X,y=r1$Y)
d4<-data.frame(x=r2$X,y=r2$Y)
d5<-data.frame(x=z1[,1],y=z1[,2])
d6<-data.frame(x=z2[,1],y=z2[,2])

# Plot style
th<-function(){
  theme_minimal(base_size=20)+theme(panel.grid=element_blank(),
    axis.line=element_line(color="black"),axis.ticks=element_line(color="black"),
    axis.title=element_text(face="bold"),axis.text=element_text(color="black"))
}

# Plots
p1<-ggplot(d1,aes(s,k))+geom_line(linewidth=1.6,color="blue")+
  labs(x="t",y=expression(kappa))+th()
p2<-ggplot(d3,aes(x,y))+geom_path(linewidth=1.6,color="orange")+
  labs(x="x",y="y")+th()
p3<-ggplot(d5,aes(x,y))+geom_path(linewidth=1.6,color="green")+
  labs(x="x",y="y")+th()

p4<-ggplot(d2,aes(s,k))+geom_line(linewidth=1.6,color="blue")+
  labs(x="t",y=expression(kappa))+th()
p5<-ggplot(d4,aes(x,y))+geom_path(linewidth=1.6,color="orange")+
  labs(x="x",y="y")+th()
p6<-ggplot(d6,aes(x,y))+geom_path(linewidth=1.6,color="green")+
  labs(x="x",y="y")+th()

# Combine plots
(p1|p2|p3)/(p4|p5|p6)
