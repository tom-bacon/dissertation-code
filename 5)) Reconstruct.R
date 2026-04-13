# This script defines a circle and its constant curvature,
# then reconstructs the curve from curvature to compare
# the original and reconstructed shapes.

library(ggplot2)
library(fdasrvf)
library(pracma)
library(stats)

R <- 1
n <- 100

# Parameter and original circle
s <- seq(0,2*pi*R,length.out=n)  
x <- R*cos(s)
y <- R*sin(s)

# Constant curvature of a circle
kappa <- rep(1/R,n)
dt <- s[2]-s[1]

# Plot curvature
plot(kappa,type="l",xlab="t",ylab="Curvature")

# Reconstruct curve from curvature
reconstruct <- function(kappa,dt=2*pi/length(kappa),x0=0,y0=0,theta0=0) {
  theta_cum <- theta0+cumsum(kappa*dt)     # Turning angle
  x_recon <- x0+cumsum(cos(theta_cum)*dt)  # x from integration
  y_recon <- y0+cumsum(sin(theta_cum)*dt)  # y from integration
  return(list(x=x_recon,y=y_recon))
}
recon <- reconstruct(kappa,dt,x0=0,y0=-R+1,theta0=0)

# Set plotting limits
ymin <- min(c(y,recon$y))
ymax <- max(c(y,recon$y))

# Plot original vs reconstructed curve
plot(x,y,type="n",col="blue",lwd=4,asp=1,ylim=c(ymin,ymax),xlab="x",ylab="y",
    font.lab=2,cex.lab=1.7,cex.axis=1.5)
grid(col="lightgrey",lty=1) 
lines(x,y,col="blue",lwd=3,lty=1)
lines(recon$x,recon$y,col="red",lwd=4,lty=2)
legend("topright",legend=c("Original","Remade"),col=c("blue","red"),
    lty=c(1,2),lwd=5,bty="n",seg.len=3,cex=1.5)
