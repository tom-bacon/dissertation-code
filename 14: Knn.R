# This script applies k-NN classification, and compares performance
# across features. Example curves are also visualised.

library(signal)
library(pracma)
library(class)
library(ggplot2)
library(patchwork)

# Load dataset
d<-readRDS("curve_dataset.rds")
curves<-d$curves
cls<-factor(d$labels)

# Compute curvature (with smoothing)
curv<-function(x,y,t=NULL){
  n<-length(x)
  if(is.null(t)) t<-seq(0,1,length.out=n)
  w<-max(7,floor(n/10))
  if(w%%2==0) w<-w+1
  xs<-sgolayfilt(x,4,w)
  ys<-sgolayfilt(y,4,w)
  dxi<-sgolayfilt(xs,4,w,m=1)
  dyi<-sgolayfilt(ys,4,w,m=1)
  ddxi<-sgolayfilt(xs,4,w,m=2)
  ddyi<-sgolayfilt(ys,4,w,m=2)
  dt<-diff(t);dt<-c(dt[1],dt)
  dx<-dxi/dt;dy<-dyi/dt
  ddx<-ddxi/(dt^2);ddy<-ddyi/(dt^2)
  den<-(dx^2+dy^2)^(3/2)+1e-6
  k<-(dx*ddy-ddx*dy)/den
  list(x=xs,y=ys,k=k,t=t)
}

# Normalised curvature (scale-invariant)
normk<-function(x,y,k){
  dx<-diff(x);dy<-diff(y)
  l<-sum(sqrt(dx^2+dy^2))
  k*l
}

# Extract extrema of curvature
ext<-function(k,t){
  thr<-0.5*sd(k)
  d<-max(3,floor(length(k)/40))
  m1<-findpeaks(k,minpeakheight=thr,minpeakdistance=d)
  m2<-findpeaks(-k,minpeakheight=thr,minpeakdistance=d)
  p1<-if(is.null(m1)) numeric(0) else m1[,2]
  p2<-if(is.null(m2)) numeric(0) else m2[,2]
  p<-sort(c(p1,p2))
  m<-max(3,floor(length(k)/50))
  p<-p[p>m & p<(length(k)-m)]
  pf<-c(1,p,length(k))
  list(values=k[pf],pos=t[pf])
}

# Pad feature vectors to equal length
pad<-function(v,n){
  if(length(v)==0) v<-c(0,0)
  while(length(v)<n){
    g<-abs(diff(v))
    i<-which.max(g)
    nv<-mean(v[i:(i+1)])
    v<-append(v,nv,after=i)
  }
  list(values=v)
}

padp<-function(v,pos,n){
  if(length(v)==0){
    v<-c(0,0);pos<-c(0,1)
  }
  idx<-c()
  while(length(v)<n){
    g<-abs(diff(v))
    i<-which.max(g)
    nv<-mean(v[i:(i+1)])
    np<-mean(pos[i:(i+1)])
    v<-append(v,nv,after=i)
    pos<-append(pos,np,after=i)
    idx<-c(idx,i+1)
  }
  list(values=v,pos=pos,idx=idx)
}

n<-length(curves)

# Feature containers
f1<-list();f2<-list();f3<-list();f4<-list()

# Extract features
for(i in 1:n){
  z<-curves[[i]]
  r<-curv(z[,1],z[,2])
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
  tmp<-pad(v,maxn)
  f4p[[i]]<-sort(tmp$values)
}

m4<-do.call(rbind,f4p)

# k-NN evaluation
ks<-c(1,3,5,7,9,11)
runs<-50
resall<-array(0,c(runs,length(ks),4))

for(r in 1:runs){
  id<-sample(1:n,round(0.7*n))
  tr1<-m1[id,];te1<-m1[-id,]
  tr2<-m2[id,];te2<-m2[-id,]
  tr3<-m3[id,];te3<-m3[-id,]
  tr4<-m4[id,];te4<-m4[-id,]
  ytr<-cls[id];yte<-cls[-id]
  resall[r,,1]<-sapply(ks,function(k) mean(knn(tr1,te1,ytr,k)==yte))
  resall[r,,2]<-sapply(ks,function(k) mean(knn(tr2,te2,ytr,k)==yte))
  resall[r,,3]<-sapply(ks,function(k) mean(knn(tr3,te3,ytr,k)==yte))
  resall[r,,4]<-sapply(ks,function(k) mean(knn(tr4,te4,ytr,k)==yte))
}

# Average results
r1<-colMeans(resall[,,1])
r2<-colMeans(resall[,,2])
r3<-colMeans(resall[,,3])
r4<-colMeans(resall[,,4])

res<-data.frame(k=ks,curves=r1,curvature=r2,norm_curvature=r3,extrema=r4)
print(res)

# Plot examples
i<-283
j<-320
th<-function(){
  theme_minimal(base_size=16)+theme(panel.grid=element_blank(),
    axis.line=element_line(color="black"),axis.ticks=element_line(color="black"),
    axis.title=element_text(face="bold"),axis.text=element_text(color="black"))
}

plotg<-function(idx){
  z<-curves[[idx]]
  x<-z[,1];y<-z[,2]
  t<-seq(0,1,length.out=length(x))
  r<-curv(x,y,t)
  kn<-normk(r$x,r$y,r$k)
  e<-f4[[idx]]
  padv<-padp(e$values,e$pos,maxn)
  df1<-data.frame(x=x,y=y)
  df2<-data.frame(x=r$x,y=r$y)
  df3<-data.frame(t=t,k=kn)
  df4<-data.frame(t=e$pos,k=e$values)
  df5<-if(length(padv$idx)>0)
    data.frame(t=padv$pos[padv$idx],k=padv$values[padv$idx]) else NULL
  df6<-data.frame(i=1:length(padv$values),v=padv$values)
  
  p1<-ggplot()+
    geom_path(data=df1,aes(x,y),linewidth=1.2,color="grey70")+
    geom_path(data=df2,aes(x,y),linewidth=1.7,color="green")+
    labs(x="x",y="y")+th()
  
  p2<-ggplot(df3,aes(t,k))+
    geom_line(linewidth=1.2,color="black")+
    geom_point(data=df4,aes(t,k),color="red",size=3)+
    {if(!is.null(df5)) geom_point(data=df5,aes(t,k),color="blue",size=3)}+
    labs(x="t",y=expression(kappa))+th()
  
  p3<-ggplot(df6,aes(i,v))+
    geom_line(linewidth=1.2,color="purple")+
    geom_point(size=2,color="purple")+
    labs(x="index",y="value")+th()
  
  p1|p2|p3
}

# Show examples
plotg(i)/plotg(j)
