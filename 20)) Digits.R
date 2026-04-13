# This script converts MNIST images into curve representations
# by extracting contours, then visualises both the images and
# their corresponding curves.

library(keras)

mnist<-dataset_mnist()
x<-mnist$train$x
y<-mnist$train$y

# Extract main contour from image
imgcurve<-function(img){
  img<-img/max(img)
  b<-img>0.3
  cl<-contourLines(b)
  if(length(cl)==0) return(NULL)
  l<-sapply(cl,function(z) length(z$x))
  c<-cl[[which.max(l)]]
  xx<-c$x-mean(c$x)
  yy<-c$y-mean(c$y)
  cbind(xx,yy)
}

# Resample curve to fixed number of points
resample<-function(curve,n){
  if(is.null(curve)) return(NULL)
  t1<-seq(0,1,length.out=nrow(curve))
  t2<-seq(0,1,length.out=n)
  xx<-approx(t1,curve[,1],t2)$y
  yy<-approx(t1,curve[,2],t2)$y
  cbind(xx,yy)
}

# Simple smoothing using moving average
smooth<-function(curve){
  if(is.null(curve)) return(NULL)
  n<-nrow(curve)
  w<-round(n*0.07)
  if(w%%2==0) w<-w+1
  if(n<w) return(curve)
  p<-floor(w/2)
  xe<-c(curve[(n-p+1):n,1],curve[,1],curve[1:p,1])
  ye<-c(curve[(n-p+1):n,2],curve[,2],curve[1:p,2])
  xs<-stats::filter(xe,rep(1/w,w),sides=2)
  ys<-stats::filter(ye,rep(1/w,w),sides=2)
  xx<-xs[(p+1):(p+n)]
  yy<-ys[(p+1):(p+n)]
  cbind(xx,yy)
}

# Ensure curve is closed
closec<-function(curve){
  if(is.null(curve)) return(NULL)
  if(sum((curve[1,]-curve[nrow(curve),])^2)>1e-6)
    curve<-rbind(curve,curve[1,])
  curve
}

curves<-list()
imgs<-list()
lab<-c()
n<-dim(x)[1]

# Convert images to curves
for(i in 1:n){
  img<-x[i,,]
  c<-imgcurve(img)
  if(is.null(c)) next
  c<-resample(c,199)
  c<-smooth(c)
  c<-closec(c)
  curves[[length(curves)+1]]<-c
  imgs[[length(imgs)+1]]<-img
  lab<-c(lab,y[i])
}

cls<-factor(lab)
saveRDS(list(curves=curves,labels=cls),"mnist_curves.rds")

# Plot one example per digit
idx<-sapply(0:9,function(d) which(cls==d)[1])
par(mfrow=c(5,4),mar=c(1,1,1,1))
for(i in 1:10){
  j<-idx[i]
  img<-imgs[[j]]
  c<-curves[[j]]
  image(1:28,1:28,t(apply(img,2,rev)),col=gray((0:255)/255),
        xlab="",ylab="",axes=FALSE)
  lim<-max(abs(c))
  plot(cbind(c[,2],-c[,1]),type="l",col="blue",lwd=3,asp=1,
       xlim=c(-lim,lim),ylim=c(-lim,lim),xlab="",ylab="",axes=FALSE)
}
