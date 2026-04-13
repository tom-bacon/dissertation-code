# This script extracts evaluates classification performance using 
# a simple neural network across different representations.

library(signal)
library(pracma)
library(keras3)

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
  v
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
  f4p[[i]]<-padinsert(v,maxn)
}
m4<-do.call(rbind,f4p)

# Standardise data
scale_pair<-function(tr,te){
  mu<-colMeans(tr)
  sdv<-apply(tr,2,sd)
  sdv[sdv==0]<-1
  list(tr=scale(tr,center=mu,scale=sdv),te=scale(te,center=mu,scale=sdv))
}

# Simple neural network
make_model<-function(input_dim,nclass){
  keras_model_sequential()|>
    layer_dense(units=16,activation="relu",input_shape=input_dim)|>
    layer_dense(units=nclass,activation="softmax")|>
    compile(optimizer=optimizer_adam(learning_rate=0.001),
      loss="sparse_categorical_crossentropy",
      metrics="accuracy")
}

runs<-50
resall<-matrix(0,nrow=runs,ncol=4)
y_all<-as.integer(cls)-1
nclass<-length(levels(cls))

# Training loop
for(r in 1:runs){
  id<-sample(1:n,round(0.7*n))
  tr1<-m1[id,];te1<-m1[-id,]
  tr2<-m2[id,];te2<-m2[-id,]
  tr3<-m3[id,];te3<-m3[-id,]
  tr4<-m4[id,];te4<-m4[-id,]
  ytr<-y_all[id];yte<-y_all[-id]
  s1<-scale_pair(tr1,te1)
  s2<-scale_pair(tr2,te2)
  s3<-scale_pair(tr3,te3)
  s4<-scale_pair(tr4,te4)
  m1_fit<-make_model(ncol(s1$tr),nclass)
  m2_fit<-make_model(ncol(s2$tr),nclass)
  m3_fit<-make_model(ncol(s3$tr),nclass)
  m4_fit<-make_model(ncol(s4$tr),nclass)
  m1_fit|>fit(s1$tr,ytr,epochs=20,batch_size=64,verbose=0)
  m2_fit|>fit(s2$tr,ytr,epochs=20,batch_size=64,verbose=0)
  m3_fit|>fit(s3$tr,ytr,epochs=20,batch_size=64,verbose=0)
  m4_fit|>fit(s4$tr,ytr,epochs=20,batch_size=64,verbose=0)
  p1<-max.col(predict(m1_fit,s1$te))-1
  p2<-max.col(predict(m2_fit,s2$te))-1
  p3<-max.col(predict(m3_fit,s3$te))-1
  p4<-max.col(predict(m4_fit,s4$te))-1
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