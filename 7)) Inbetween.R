# This script compares two similar curvature functions, finds their
# extrema, and reconstructs the corresponding curves to see
# how small changes in curvature affect shape.

n <- 4000
t <- seq(0,2*pi,length.out=n)
dt <- t[2]-t[1]
eps <- 0.3

# Two curvature functions (one slightly perturbed)
k1 <- sin(t)
k2 <- sin(t)+eps*(cos(t)^2)*sin(2*t)

# Numerical derivatives
d1 <- diff(k1)/dt
d2 <- diff(k2)/dt

# Find extrema via sign changes
ext_idx <- function(d){
  s <- sign(d)
  which(s[-1]!=s[-length(s)])+1
}
e1 <- ext_idx(d1)
e2 <- ext_idx(d2)

# Classify extrema as max or min
class_ext <- function(k,idx){
  out <- data.frame(t=t[idx],k=k[idx],type=NA)
  for(i in seq_along(idx)){
    j <- idx[i]
    if(j>1 && j<length(k)){
      if(k[j]>k[j-1] && k[j]>k[j+1]) out$type[i] <- "max"
      if(k[j]<k[j-1] && k[j]<k[j+1]) out$type[i] <- "min"
    }
  }
  out
}
E1 <- class_ext(k1,e1)
E2 <- class_ext(k2,e2)

# Plot curvature functions
par(mgp=c(2.5,1,0))
plot(t,k1,type="n",xlab=expression(t),ylab=expression(kappa),
     ylim=range(c(k1,k2)),cex.lab=1.7,cex.axis=1.5,
     font.lab=2,font.axis=1)
grid(col="gray88",lty=1)
lines(t,k1,lwd=2.8,col="blue")
lines(t,k2,lwd=2.8,col="red")

# Mark key extrema of k1
points(c(pi/2,3*pi/2),c(1,-1),pch=19,col="black",cex=1.2)
legend("topright",legend=c(expression(kappa[1](t)),expression(kappa[2](t))),
       col=c("blue","red"),lwd=2.8,bty="n",cex=1.5)

# Reconstruct curve from curvature
recon <- function(k){
  th <- cumsum(k)*dt
  x <- cumsum(cos(th))*dt
  y <- cumsum(sin(th))*dt
  list(x=x,y=y,th=th)
}
C1 <- recon(k1)
C2 <- recon(k2)

# Plot reconstructed curves
par(mgp=c(2.5,1,0))
plot(C1$x,C1$y,type="n",asp=1,xlab="x",ylab="y",xlim=range(c(C1$x,C2$x)),
     ylim=range(c(C1$y,C2$y)),cex.lab=1.7,cex.axis=1.5,
     font.lab=2,font.axis=1)
grid(col="gray88",lty=1)
lines(C1$x,C1$y,lwd=2.8,col="blue")
lines(C2$x,C2$y,lwd=2.8,col="red")
legend("topleft",legend=c(expression("Curve from "*kappa[1](t)),
     expression("Curve from "*kappa[2](t))),col=c("blue","red"),
     lwd=2.8,bty="n",cex=1.5)

