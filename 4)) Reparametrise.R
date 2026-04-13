# This script defines a curve and two reparameterisations,
# then compares their curvature to show how it behaves
# under these different parameter changes.

library(ggplot2)

# Reparameterisation functions
g <- function(t) t^2/(2*pi)
p <- function(t) sqrt(t/(2*pi))*2*pi

# Original curve and its derivatives
x1 <- function(t) t
x2 <- function(t) sin(t)
x1p <- function(t) 1
x2p <- function(t) cos(t)
x1pp <- function(t) 0
x2pp <- function(t) -sin(t)
t <- seq(1e-4,2*pi,length.out=400)

# Plot original curve
df1 <- data.frame(t=t,x=x1(t),y=x2(t))
ggplot(df1,aes(x=x,y=y))+geom_path(color="blue",linewidth=1)+
  scale_x_continuous(breaks=c(0,pi/2,pi,3*pi/2,2*pi),
  labels=c("0",expression(pi/2),expression(pi),
  expression(3*pi/2),expression(2*pi)))+coord_equal()+labs(x="x",y="y")+
  theme_minimal(base_size=14)+theme(axis.title=element_text(face="bold"),
  axis.text=element_text(size=16))

# Curvature of original curve
kX <- function(t){
  n <- x1p(t)*x2pp(t)-x2p(t)*x1pp(t)
  d <- (x1p(t)^2+x2p(t)^2)^(3/2)
  n/d
}

# First reparameterisation 
y1 <- function(t) x1(g(t))
y2 <- function(t) x2(g(t))
y1p <- function(t) 2*t*x1p(g(t))
y2p <- function(t) 2*t*x2p(g(t))
y1pp <- function(t) x1pp(g(t))*(2*t)^2+x1p(g(t))*2
y2pp <- function(t) x2pp(g(t))*(2*t)^2+x2p(g(t))*2
kY <- function(t){
  n <- y1p(t)*y2pp(t)-y2p(t)*y1pp(t)
  d <- (y1p(t)^2+y2p(t)^2)^(3/2)
  n/d
}

# Second reparameterisation 
z1 <- function(t) x1(p(t))
z2 <- function(t) x2(p(t))
z1p <- function(t) 0.5/sqrt(t)*x1p(p(t))
z2p <- function(t) 0.5/sqrt(t)*x2p(p(t))
z1pp <- function(t) -0.25/t^(3/2)*x1p(p(t))+(0.5/sqrt(t))^2*x1pp(p(t))
z2pp <- function(t) -0.25/t^(3/2)*x2p(p(t))+(0.5/sqrt(t))^2*x2pp(p(t))
kZ <- function(t){
  n <- z1p(t)*z2pp(t)-z2p(t)*z1pp(t)
  d <- (z1p(t)^2+z2p(t)^2)^(3/2)
  n/d
}

# Combine curvature results
df2 <- data.frame(t=rep(t,3),k=c(kX(t),kY(t),kZ(t)),
  ctype=factor(rep(c("X","XoG","XoP"),each=length(t))))

# Reference curve
df3 <- data.frame(t=t,y=sin(t),ctype="X")

# Plot curvature comparison
ggplot(df2,aes(x=t,y=k,color=ctype))+geom_line(linewidth=2)+
  geom_line(data=df3,aes(x=t,y=y,linetype=ctype),color="black",linewidth=2,
  inherit.aes=FALSE)+scale_x_continuous(breaks=c(0,pi/2,pi,3*pi/2,2*pi),
  labels=c("0",expression(pi/2),expression(pi),expression(3*pi/2),
  expression(2*pi)))+scale_color_manual(values=c("blue","red","green"),
  labels=c(expression(kappa[X](t)),expression(kappa[X*" o "*gamma](t)),
  expression(kappa[X*" o "*phi](t))))+scale_linetype_manual(name="Curve",
  values=c("X"="dashed"))+labs(x="t",y=expression(kappa),color="Curvature")+
  theme_minimal(base_size=14)+theme(axis.title=element_text(size=24,face="bold"),
  axis.text=element_text(size=20),legend.title=element_text(size=22,face="bold"),
  legend.text=element_text(size=20))

