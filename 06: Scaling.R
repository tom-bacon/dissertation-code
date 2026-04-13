# This script defines a curve and a scaled version of it,
# then compares their curvature and shapes to show how scaling
# affects curvature.

library(ggplot2)

# Original curve and derivatives
X1 <- function(t) 2*cos(t)
X2 <- function(t) sin(t)
X1_prime <- function(t) -2*sin(t)
X2_prime <- function(t) cos(t)
X1_double <- function(t) -2*cos(t)
X2_double <- function(t) -sin(t)

# Curvature of original curve
kappa_X <- function(t) {
  num <- X1_prime(t)*X2_double(t)-X2_prime(t)*X1_double(t)
  den <- (X1_prime(t)^2+X2_prime(t)^2)^(3/2)
  num/den
}

# Scale factor
a <- 2

# Scaled curve and derivatives
Y1 <- function(t) a*X1(t)
Y2 <- function(t) a*X2(t)
Y1_prime <- function(t) a*X1_prime(t)
Y2_prime <- function(t) a*X2_prime(t)
Y1_double <- function(t) a*X1_double(t)
Y2_double <- function(t) a*X2_double(t)

# Curvature of scaled curve
kappa_Y <- function(t) {
  num <- Y1_prime(t)*Y2_double(t)-Y2_prime(t)*Y1_double(t)
  den <- (Y1_prime(t)^2+Y2_prime(t)^2)^(3/2)
  num/den
}

t_vals <- seq(0,2*pi,length.out=400)

# Combine curvature data
data <- data.frame(t=rep(t_vals,2),kappa=c(kappa_X(t_vals),kappa_Y(t_vals)),
  curve=factor(rep(c("X(t)","2X(t)"),each=length(t_vals))))
data$curve <- factor(data$curve,levels=c("X(t)","2X(t)"))

# Plot curvature comparison
ggplot(data,aes(x=t,y=kappa,color=curve))+geom_line(size=1.5)+
  scale_color_manual(values=c("blue","red"),
  label=c(expression(kappa[X](t)),expression(kappa[2*X](t))))+
  scale_x_continuous(breaks=c(0,pi/2,pi,3*pi/2,2*pi),labels=c(expression(0),
  expression(pi/2),expression(pi),expression(3*pi/2),expression(2*pi)))+
  theme_minimal(base_size=14)+labs(x="t",y=expression(kappa),color="Curve")+
  theme(axis.title=element_text(size=16,face="bold"),
  axis.text=element_text(size=16),plot.title=element_text(hjust=0.5,face="bold"),
  legend.position="right",legend.title=element_text(size=16,face="bold"),
  legend.text=element_text(size=16))

# Curve shapes
df_X <- data.frame(x=X1(t_vals),y=X2(t_vals),curve="X(t)")
df_Y <- data.frame(x=Y1(t_vals),y=Y2(t_vals),curve="2X(t)")
df_all <- rbind(df_X,df_Y)
df_all$curve <- factor(df_all$curve,levels=c("X(t)","2X(t)"))

# Plot original vs scaled curve
ggplot(df_all, aes(x=x,y=y,color=curve))+geom_path(size = 1.2)+coord_equal()+
  scale_color_manual(values=c("X(t)"="blue","2X(t)"="red"),
  labels=c(expression(X(t)),expression(2*X(t))))+scale_x_continuous(
  breaks=c(-pi,-pi/2,0,pi/2,pi),labels=c(expression(-pi),expression(-pi/2),
  expression(0),expression(pi/2),expression(pi)))+labs(x="x",y="y",
  color="Curve")+theme_minimal(base_size=14)+
  theme(axis.title=element_text(size=16,face="bold"),
  axis.text=element_text(size=16),plot.title=element_text(hjust=0.5,face="bold"),
  legend.position="right",legend.title=element_text(size=16,face="bold"),
  legend.text=element_text(size=16))
