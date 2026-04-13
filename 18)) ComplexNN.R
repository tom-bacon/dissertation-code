# This script uses a more complex neural network to compare classification
#  performance across different representations and training set sizes.

library(signal)
library(pracma)
library(keras3)

d<-readRDS("curve_dataset.rds")
curves<-d$curves
cls<-factor(d$labels)

n<-length(curves)
y_all<-as.integer(cls)-1
nclass<-length(levels(cls))

# Compute curvature with smoothing
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

# Scale-invariant curvature
normk<-function(x,y,k){
  dx<-diff(x);dy<-diff(y)
  l<-sum(sqrt(dx^2+dy^2))
  k*l
}

# Extract key curvature points (extrema)
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
  list(values=k[pos_full])
}

# Pad vectors so all have same length
padinsert<-function(v,maxn){
  if(length(v)==0) v<-c(0,0)
  while(length(v)<maxn){
    gaps<-abs(diff(v))
    i<-which.max(gaps)
    v<-append(v,mean(v[i:(i+1)]),after=i)
  }
  v
}

# Build feature sets
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
  f4[[i]]<-e$values
}
m1<-do.call(rbind,f1)
m2<-do.call(rbind,f2)
m3<-do.call(rbind,f3)
maxn<-max(sapply(f4,length))
f4p<-lapply(f4,padinsert,maxn=maxn)
m4<-do.call(rbind,f4p)

# Standardise training/test data
scale_pair<-function(tr,te){
  mu<-colMeans(tr)
  sdv<-apply(tr,2,sd)
  sdv[sdv==0]<-1
  list(tr=scale(tr,mu,sdv),te=scale(te,mu,sdv))
}

#  Deeper neural network
make_model_deeper<-function(input_dim,nclass){
  keras_model_sequential()|>
    layer_dense(128,activation="relu",input_shape=input_dim)|>
    layer_batch_normalization()|>
    layer_dropout(0.4)|>
    layer_dense(64,activation="relu")|>
    layer_batch_normalization()|>
    layer_dropout(0.3)|>
    layer_dense(32,activation="relu")|>
    layer_dense(nclass,activation="softmax")|>
    compile(optimizer=optimizer_adam(learning_rate=0.0008),
      loss="sparse_categorical_crossentropy",
      metrics="accuracy")
}

# Run experiment for a given training size
run_experiment<-function(train_frac){
  runs<-50
  res<-matrix(0,runs,4)
  for(r in 1:runs){
    id<-sample(1:n,round(train_frac*n))
    ytr<-y_all[id];yte<-y_all[-id]
    s1<-scale_pair(m1[id,],m1[-id,])
    s2<-scale_pair(m2[id,],m2[-id,])
    s3<-scale_pair(m3[id,],m3[-id,])
    s4<-scale_pair(m4[id,],m4[-id,])
    data_tr<-list(s1$tr,s2$tr,s3$tr,s4$tr)
    data_te<-list(s1$te,s2$te,s3$te,s4$te)
    for(i in 1:4){
      model<-make_model_deeper(ncol(data_tr[[i]]),nclass)
      model|>fit(data_tr[[i]],ytr,epochs=35,batch_size=64,verbose=0)
      pred<-max.col(predict(model,data_te[[i]]))-1
      res[r,i]<-mean(pred==yte)
    }
  }
  colMeans(res)
}

# Compare full vs small training data
res_full<-run_experiment(0.7)
res_small<-run_experiment(0.1)

results<-data.frame(method=c("curves","curvature","norm_curvature","extrema"),
  full_data=res_full,small_data=res_small)

print(results)