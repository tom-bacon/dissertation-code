# This script extracts features from MNIST-derived curves and
# evaluates classification using LDA across different feature
# representations, with example visualisations.

library(signal)
library(pracma)
library(MASS)

d<-readRDS("mnist_curves.rds")
cs<-d$curves
cl<-factor(d$labels)

# Resample curve to fixed length
resamp<-function(c,n=120){
  if(is.null(c)) return(NULL)
  x<-c[,1];y<-c[,2]
  s<-arc(x,y)
  t<-seq(0,1,length.out=n)
  cbind(approx(s,x,t)$y,approx(s,y,t)$y,t)
}

# Curvature with smoothing
curv<-function(x,y){
  n<-length(x)
  w<-max(7,floor(n/10));if(w%%2==0) w<-w+1
  xs<-sgolayfilt(x,3,w);ys<-sgolayfilt(y,3,w)
  dx<-sgolayfilt(xs,3,w,m=1);dy<-sgolayfilt(ys,3,w,m=1)
  ddx<-sgolayfilt(xs,3,w,m=2);ddy<-sgolayfilt(ys,3,w,m=2)
  d<-(dx^2+dy^2)^(3/2);d[d<1e-6]<-1e-6
  (dx*ddy-dy*ddx)/d
}
arc<-function(x,y){
  s<-c(0,cumsum(sqrt(diff(x)^2+diff(y)^2)))
  if(max(s)==0) return(rep(0,length(s)))
  s/max(s)
}

# Normalise curvature magnitude
normk<-function(k){
  m<-mean(abs(k),na.rm=TRUE)
  if(is.na(m)||m<1e-6) return(k)
  k/m
}

# Extract curvature extrema 
ext<-function(k){
  if(length(k)<10) return(list(p=numeric(0),v=numeric(0)))
  a<-findpeaks(k,minpeakdistance=5,minpeakheight=0.5)
  b<-findpeaks(-k,minpeakdistance=5,minpeakheight=0.5)
  p<-c();v<-c()
  if(!is.null(a)){p<-c(p,a[,2]);v<-c(v,a[,1])}
  if(!is.null(b)){p<-c(p,b[,2]);v<-c(v,-b[,1])}
  o<-order(p)
  list(p=p[o],v=v[o])
}

# Pad feature vectors to equal length
pad<-function(v,n){
  if(length(v)==0) return(rep(0,n))
  while(length(v)<n){
    if(length(v)==1){v<-c(v,v);next}
    g<-abs(diff(v));i<-which.max(g)
    v<-append(v,(v[i]+v[i+1])/2,after=i)
  }
  v[1:n]
}

# Standardise features
sc<-function(a,b){
  m<-colMeans(a,na.rm=TRUE)
  s<-apply(a,2,sd,na.rm=TRUE)
  s[s==0|is.na(s)]<-1
  list(a=scale(a,m,s),b=scale(b,m,s))
}

f1<-list();f2<-list();f3<-list();f4<-list()

# Extract features for each curve
for(i in seq_along(cs)){
  c<-cs[[i]]
  if(is.null(c)||!is.matrix(c)) next
  c<-resamp(c,120)
  x<-c[,1];y<-c[,2]
  k<-curv(x,y);kn<-normk(k)
  e<-ext(kn)
  f1[[i]]<-c(x,y)
  f2[[i]]<-k
  f3[[i]]<-kn
  f4[[i]]<-e$v
}

# Remove invalid entries
keep<-which(!sapply(f1,is.null))
f1<-f1[keep];f2<-f2[keep];f3<-f3[keep];f4<-f4[keep]
cl<-cl[keep]
m1<-do.call(rbind,f1)
m2<-do.call(rbind,f2)
m3<-do.call(rbind,f3)

# Pad extrema features
mx<-max(sapply(f4,length))
m4<-do.call(rbind,lapply(f4,function(v) pad(v,mx)))
res<-matrix(0,10,4)
n<-length(f1)

# Repeated train/test splits
for(i in 1:10){
  id<-sample(1:n,round(0.7*n))
  s1<-sc(m1[id,],m1[-id,])
  s2<-sc(m2[id,],m2[-id,])
  s3<-sc(m3[id,],m3[-id,])
  s4<-sc(m4[id,],m4[-id,])
  y1<-cl[id];y2<-cl[-id]
  p1<-predict(lda(s1$a,grouping=y1),s1$b)$class
  p2<-predict(lda(s2$a,grouping=y1),s2$b)$class
  p3<-predict(lda(s3$a,grouping=y1),s3$b)$class
  p4<-predict(lda(s4$a,grouping=y1),s4$b)$class
  res[i,]<-c(mean(p1==y2),mean(p2==y2),mean(p3==y2),mean(p4==y2))
}

out<-data.frame(method=c("raw","curv","norm","ext"),acc=colMeans(res))
print(out)

#  Example curves and curvature
ix<-c(16,77,238)
par(mfrow=c(3,2),mar=c(3.2,3.2,2.2,1),oma=c(0,0,1,0),cex=1.1,cex.axis=1.1,
    cex.lab=1.2,font.main=2)

for(i in ix){
  c<-cs[[i]]
  if(is.null(c)||!is.matrix(c)) next
  c<-resamp(c,120)
  x<-c[,1];y<-c[,2];t<-c[,3]
  k<-curv(x,y)
  kn<-normk(k);kn[is.na(kn)]<-0
  e<-ext(kn)
  id<-e$p
  plot(x,y,type="l",lwd=2.5,col="black",asp=1,
       xlab="x",ylab="y",main=paste("Curve",i),
       panel.first=grid(col="grey85",lwd=1))
  points(x[id],y[id],pch=19,col="red",cex=1)
  plot(t,kn,type="l",lwd=2.2,col="blue",xlab="t",ylab=expression(kappa),
       main=paste("Curvature",i),panel.first=grid(col="grey85",lwd=1))
  points(t[id],kn[id],pch=19,col="red",cex=1)
}
