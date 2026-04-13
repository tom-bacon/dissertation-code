# This script evaluates classification performance using SVM 
# across different representations.

library(signal)
library(pracma)
library(e1071)

# Load dataset
d<-readRDS("curve_dataset.rds")
curves<-d$curves
cls<-factor(d$labels)

# Curvature with smoothing
curv<-function(x,y,t=NULL){
  n<-length(x)
  if(is.null(t)) t<-seq(0,1,length.out=n)
  w<-max(7,floor(n/10))
  if(w%%2==0) w<-w+1
  xs<-sgolayfilt(x,4,w)
  ys<-sgolayfilt(y,4,w)
  dx_i<-sgolayfilt(xs,4,w,m=1)
  dy_i<-sgolayfilt(ys,4,w,m=1)
  ddx_i<-sgolayfilt(xs,4,w,m=2)
  ddy_i<-sgolayfilt(ys,4,w,m=2)
  dt<-diff(t);dt<-c(dt[1],dt)
  dx<-dx_i/dt;dy<-dy_i/dt
  ddx<-ddx_i/(dt^2);ddy<-ddy_i/(dt^2)
  den<-(dx^2+dy^2)^(3/2)+1e-6
  k<-(dx*ddy-ddx*dy)/den
  list(x=xs,y=ys,k=k,t=t)
}

# Normalised curvature
normk<-function(x,y,k){
  dx<-diff(x);dy<-diff(y)
  l<-sum(sqrt(dx^2+dy^2))
  k*l
}

# Extract extrema
ext<-function(k,t){
  thr<-0.5*sd(k)
  dist<-max(3,floor(length(k)/40))
  m1<-findpeaks(k,minpeakheight=thr,minpeakdistance=dist)
  m2<-findpeaks(-k,minpeakheight=thr,minpeakdistance=dist)
  p1<-if(is.null(m1)) numeric(0) else m1[,2]
  p2<-if(is.null(m2)) numeric(0) else m2[,2]
  pos<-sort(c(p1,p2))
  margin<-max(3,floor(length(k)/50))
  pos<-pos[pos>margin & pos<(length(k)-margin)]
  pos_full<-c(1,pos,length(k))
  list(values=k[pos_full],pos=t[pos_full])
}

# Pad extrema vectors
padinsert<-function(v,maxn){
  if(length(v)==0) v<-c(0,0)
  while(length(v)<maxn){
    gaps<-abs(diff(v))
    i<-which.max(gaps)
    nv<-mean(v[i:(i+1)])
    v<-append(v,nv,after=i)
  }
  list(values=v)
}
n<-length(curves)

# Feature extraction
f1<-list();f2<-list();f3<-list();f4<-list()
for(i in 1:n){
  x<-curves[[i]][,1]
  y<-curves[[i]][,2]
  r<-curv(x,y)
  k<-r$k
  kn<-normk(r$x,r$y,k)
  e<-ext(kn,r$t)
  f1[[i]]<-c(r$x,r$y)
  f2[[i]]<-k
  f3[[i]]<-kn
  f4[[i]]<-e
}

# Convert to matrices
m1<-do.call(rbind,f1)
m2<-do.call(rbind,f2)
m3<-do.call(rbind,f3)

# Pad extrema features
maxn<-max(sapply(f4,function(z) length(z$values)))
f4p<-list()

for(i in 1:n){
  v<-f4[[i]]$values
  tmp<-padinsert(v,maxn)
  f4p[[i]]<-tmp$values
}

m4<-do.call(rbind,f4p)

# SVM evaluation
runs<-50
resall<-matrix(0,nrow=runs,ncol=4)
for(r in 1:runs){
  id<-sample(1:n,round(0.7*n))
  tr1<-m1[id,];te1<-m1[-id,]
  tr2<-m2[id,];te2<-m2[-id,]
  tr3<-m3[id,];te3<-m3[-id,]
  tr4<-m4[id,];te4<-m4[-id,]
  ytr<-cls[id];yte<-cls[-id]
  m1_fit<-svm(tr1,ytr,kernel="radial")
  m2_fit<-svm(tr2,ytr,kernel="radial")
  m3_fit<-svm(tr3,ytr,kernel="radial")
  m4_fit<-svm(tr4,ytr,kernel="radial")
  p1<-predict(m1_fit,te1)
  p2<-predict(m2_fit,te2)
  p3<-predict(m3_fit,te3)
  p4<-predict(m4_fit,te4)
  resall[r,1]<-mean(p1==yte)
  resall[r,2]<-mean(p2==yte)
  resall[r,3]<-mean(p3==yte)
  resall[r,4]<-mean(p4==yte)
}

# Average accuracy
r1<-mean(resall[,1])
r2<-mean(resall[,2])
r3<-mean(resall[,3])
r4<-mean(resall[,4])

res<-data.frame(method=c("curves","curvature","norm_curvature","extrema"),
  accuracy=c(r1,r2,r3,r4))

print(res)
