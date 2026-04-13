# This script computes and visualises curvature for two curves:
# a parabola and a cubic, showing both the shape and curvature.

library(ggplot2)
library(viridis)
library(grid)

# Plot style
theme_style<-theme_minimal(base_size=14)+theme(axis.title=
  element_text(size=24,face="bold"),axis.text=element_text(size=20),
  legend.title=element_text(size=24,face="bold"),legend.text=element_text(size=20))

# Curvature colour scale
scale_kappa<-scale_color_viridis_c(limits=c(0,1.3),oob=scales::squish,
  guide=guide_colorbar(direction="vertical",barheight=unit(12,"cm"),
  barwidth=unit(1.2,"cm")))

t<-seq(-2,2,length.out=100)
dt<-t[2]-t[1]

# First curve: parabola
x1<-t
y1<- -t^2
vx1<-c(NA,diff(x1)/dt)
vy1<-c(NA,diff(y1)/dt)
ax1<-c(NA,diff(vx1)/dt)
ay1<-c(NA,diff(vy1)/dt)
curv1<-abs(vx1*ay1-vy1*ax1)/((vx1^2+vy1^2)^(3/2))
df1<-data.frame(t=t,x=x1,y=y1,curvature=curv1)

# Curve + curvature plots
p1<-ggplot(df1,aes(x=x,y=y,color=curvature))+geom_path(linewidth=1.5)+scale_kappa+
  theme_style+labs(x="x",y="y",color=expression(kappa))

p2<-ggplot(df1,aes(x=t,y=curvature,color=curvature))+geom_line(linewidth=1.5)+
  scale_kappa+theme_style+labs(x="t",y=expression(kappa),color=expression(kappa))

# Second curve: cubic
x2<-t
y2<-t^3
vx2<-c(NA,diff(x2)/dt)
vy2<-c(NA,diff(y2)/dt)
ax2<-c(NA,diff(vx2)/dt)
ay2<-c(NA,diff(vy2)/dt)
curv2<-abs(vx2*ay2-vy2*ax2)/((vx2^2+vy2^2)^(3/2))
df2<-data.frame(t=t,x=x2,y=y2,curvature=curv2)

# Curve + curvature plots
p3<-ggplot(df2,aes(x=x,y=y,color=curvature))+geom_path(linewidth=1.5)+scale_kappa+
  theme_style+labs(x="x",y="y",color=expression(kappa))

p4<-ggplot(df2,aes(x=t,y=curvature,color=curvature))+geom_line(linewidth=1.5)+
  scale_kappa+theme_style+labs(x="t",y=expression(kappa),color=expression(kappa))

# Show plots
p1
p2
p3
p4
