# This script shows an example of a curve under different
# transformations: rigid motion, scaling, and reparameterisation.
# The curve and its curvature are plotted side by side.

library(pracma)

t <- seq(0,2*pi,length.out=500)

# Original curve
x <- cos(t)+0.3*cos(3*t)
y <- sin(t)

# Curvaturewith light smoothing
k <- function(x,y,t){
  sx <- splinefun(t,x)
  sy <- splinefun(t,y)
  dx  <- sx(t,deriv=1)
  dy  <- sy(t,deriv=1)
  ddx <- sx(t,deriv=2)
  ddy <- sy(t,deriv=2)
  (dx*ddy-dy*ddx)/(dx^2+dy^2)^(3/2)
}
k0 <- k(x,y,t)

# Rigid motion
th <- pi/4
xr <- x*cos(th) - y*sin(th)+2
yr <- x*sin(th) + y*cos(th)+0.5
kr <- k(xr,yr,t)

# Scaling
a <- 2
xs <- a*x
ys <- a*y
ks <- k(xs,ys,t)

# Reparameterisation
g <- t^1.5
g <- (g-min(g))/(max(g)-min(g))*2*pi
xp <- cos(g)+0.3*cos(3*g)
yp <- sin(g)
kp <- k(xp,yp,g)

# Plot limits
pad <- 0.5
xlim <- range(c(x,xr,xs,xp))
ylim <- range(c(y,yr,ys,yp))
xlim <- xlim+c(-pad,pad)
ylim <- ylim+c(-pad,pad)

# Clean curvature values for plotting
clean <- function(v) v[is.finite(v)]
klim <- range(c(clean(k0), clean(kr), clean(ks), clean(kp)))

layout(matrix(1:8,ncol=2,byrow=TRUE),widths=c(1.2,2))
par(mar=c(4,5,2,1),cex.lab=1.6,cex.axis=1.3,
    cex.main=1.4,font.lab=2,lwd=3,bg="white")

# Original
plot(x,y,type="l",col="blue",xlim=xlim,ylim=ylim,
     main="Original Curve",xlab="x",ylab="y")
plot(t,k0,type="l",col="red",ylim=klim,
     main="Curvature",xlab="t",ylab=expression(kappa))

# Rigid motion
plot(xr,yr,type="l",col="blue",xlim=xlim,ylim=ylim,
     main="Rigid Motion",xlab="x",ylab="y")
plot(t,kr,type="l",col="red",ylim=klim,
     main="Curvature",xlab="t",ylab=expression(kappa))

# Scaling
plot(xs,ys,type="l",col="blue",xlim=xlim,ylim=ylim,
     main="Scaled",xlab="x",ylab="y")
plot(t,ks,type="l",col="red",ylim=klim,
     main="Curvature",xlab="t",ylab=expression(kappa))

# Reparameterisation
plot(xp,yp,type="l",col="blue",xlim=xlim,ylim=ylim,
     main="Reparametrised",xlab="x",ylab="y")
plot(t,kp,type="l",col="red",ylim=klim,
     main="Curvature",xlab="t",ylab=expression(kappa))
