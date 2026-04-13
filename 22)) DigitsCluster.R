# This script performs unsupervised clustering on MNIST curves
# using different feature representations, and visualises the
# resulting clusters.

library(signal)
library(pracma)

d<-readRDS("mnist_curves.rds")
cs<-d$curves
minsz<-36

# Resample curves to fixed size
resamp<-function(c,n=120){
  x<-c[,1];y<-c[,2]
  s<-arc(x,y)
  t<-seq(0,1,length.out=n)
  cbind(approx(s,x,t)$y,approx(s,y,t)$y)
}

# Ensure consistent orientation
ori<-function(x,y){
  A<-sum(x*c(diff(y),y[1]-y[length(y)])-y*c(diff(x),x[1]-x[length(x)]))
  if(A>0){x<-rev(x);y<-rev(y)}
  list(x=x,y=y)
}

# Curvature with smoothing
curv<-function(x,y){
  n<-length(x)
  w<-max(7,floor(n/10))
  if(w%%2==0) w<-w+1
  xs<-sgolayfilt(x,3,w)
  ys<-sgolayfilt(y,3,w)
  dx<-sgolayfilt(xs,3,w,m=1)
  dy<-sgolayfilt(ys,3,w,m=1)
  ddx<-sgolayfilt(xs,3,w,m=2)
  ddy<-sgolayfilt(ys,3,w,m=2)
  d<-(dx^2+dy^2)^(3/2)
  d[d<1e-6]<-1e-6
  (dx*ddy-dy*ddx)/d
}
arc<-function(x,y){
  s<-c(0,cumsum(sqrt(diff(x)^2+diff(y)^2)))
  if(max(s)==0) return(rep(0,length(s)))
  s/max(s)
}

# Normalise curvature
nk<-function(k){
  m<-mean(abs(k),na.rm=TRUE)
  if(is.na(m)||m<1e-6) return(k)
  k/m
}

# Normalise vector magnitude
nv<-function(v){
  if(length(v)==0) return(v)
  s<-sqrt(sum(v^2))
  if(s<1e-6) return(v)
  v/s
}

# Extract curvature extrema
ext<-function(k){
  if(length(k)<10) return(list(v=numeric(0)))
  thr<-0.3*max(abs(k))
  a<-findpeaks(k,minpeakdistance=5,minpeakheight=thr)
  b<-findpeaks(-k,minpeakdistance=5,minpeakheight=thr)
  p<-c();v<-c()
  if(!is.null(a)){p<-c(p,a[,2]);v<-c(v,a[,1])}
  if(!is.null(b)){p<-c(p,b[,2]);v<-c(v,-b[,1])}
  o<-order(p)
  list(v=v[o])
}

# Pad vectors to equal length
pad<-function(v,n){
  if(length(v)==0) return(rep(0,n))
  while(length(v)<n){
    if(length(v)==1){v<-c(v,v);next}
    g<-abs(diff(v));i<-which.max(g)
    v<-append(v,(v[i]+v[i+1])/2,after=i)
  }
  v[1:n]
}

f1<-list();f2<-list();f3<-list();f4<-list()

# Build feature representations
for(i in seq_along(cs)){
  c<-cs[[i]]
  if(is.null(c)||!is.matrix(c)) next
  c<-resamp(c,120)
  x<-c[,1];y<-c[,2]
  o<-ori(x,y)
  x<-o$x;y<-o$y
  k<-curv(x,y)
  kn<-nk(k)
  e<-ext(kn)
  v<-nv(e$v)
  f1[[i]]<-c(x,y)
  f2[[i]]<-k
  f3[[i]]<-kn
  f4[[i]]<-v
}

keep<-which(!sapply(f1,is.null))
f1<-f1[keep];f2<-f2[keep];f3<-f3[keep];f4<-f4[keep]

# Subsample for efficiency
set.seed(1)
id<-sample(1:length(f1),min(5000,length(f1)))
f1<-f1[id];f2<-f2[id];f3<-f3[id];f4<-f4[id]
m1<-do.call(rbind,f1)
m2<-do.call(rbind,f2)
m3<-do.call(rbind,f3)
mx<-max(sapply(f4,length))
m4<-do.call(rbind,lapply(f4,function(v) pad(v,mx)))

# Normalise feature vectors
normm<-function(x) x/sqrt(rowSums(x^2))
m1<-normm(m1)
m2<-normm(m2)
m3<-normm(m3)
m4<-normm(m4)

# Remove very small clusters
filt<-function(cl,minsz){
  tab<-table(cl)
  small<-names(tab[tab<minsz])
  cl[cl%in%small]<-0
  cl
}

# Hierarchical clustering
run<-function(m,minsz){
  D<-as.dist(1-m%*%t(m))
  hc<-hclust(D,method="average")
  cl<-cutree(hc,h=0.6)
  cl<-filt(cl,minsz)
  list(cl=cl)
}

r1<-run(m1,minsz)
r2<-run(m2,minsz)
r3<-run(m3,minsz)
r4<-run(m4,minsz)

# Print cluster summaries
pr<-function(cl){
  cat("\nclusters:\n")
  print(table(cl))
  cat("\nnum clusters:",length(unique(cl[cl>0])),"\n")
  cat("noise:",sum(cl==0),"\n")
}

pr(r1$cl);pr(r2$cl);pr(r3$cl);pr(r4$cl)

# Plot clusters
plotr<-function(res){
  cl<-res$cl
  u<-sort(unique(cl));u<-u[u>0]
  k<-length(u)
  nr<-floor(sqrt(k))
  nc<-ceiling(k/nr)
  par(mfrow=c(nr,nc),mar=c(0.5,0.5,0.5,0.5))
  for(i in seq_len(k)){
    kk<-u[i]
    idk<-which(cl==kk)
    curves<-lapply(idk,function(i) cs[[keep[id[i]]]])
    curves<-lapply(curves,function(c) cbind(c[,2],-c[,1]))
    x_all<-unlist(lapply(curves,function(c) c[,1]))
    y_all<-unlist(lapply(curves,function(c) c[,2]))
    pad<-0.01
    xlim<-range(x_all);dx<-diff(xlim)
    ylim<-range(y_all);dy<-diff(ylim)
    xlim<-xlim+c(-pad,pad)*dx
    ylim<-ylim+c(-pad,pad)*dy
    plot(0,0,type="n",xlim=xlim,ylim=ylim,asp=1,
         axes=FALSE,xlab="",ylab="")
    col_line<-rgb(0,0,0,0.04)
    for(c in curves){
      lines(c,col=col_line,lwd=1.2,lend=2)
    }
  }
}

plotr(r1)
plotr(r2)
plotr(r3)
plotr(r4)
