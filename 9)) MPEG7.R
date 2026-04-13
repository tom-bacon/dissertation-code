# This script takes example curves, applies smoothing using an
# SG filter and compares curvature before and after smoothing.

library(signal)
library(ggplot2)
library(fdasrvf)
library(pracma)

# Set consistent plot style
theme_set(theme_minimal()+theme(axis.title=element_text(size=24,face="bold"),
      axis.text=element_text(size=20),legend.text=element_text(size=20),
      legend.title=element_text(size=24,face="bold")))

# Load example curves
data("beta")
b<-beta[,,1,1]
c<-beta[,,17,5]

# Plot raw curves
pb0<-ggplot(data.frame(x=b[1,],y=b[2,]),aes(x,y))+
  geom_path(linewidth=1.5,color="blue")+coord_equal()+labs(title="",x="x",y="y")

pc0<-ggplot(data.frame(x=c[1,],y=c[2,]),aes(x,y))+
  geom_path(linewidth=1.5,color="red")+coord_equal()+labs(title="",x="x",y="y")

# Savitzky–Golay smoothing 
sgolay_p<-function(x,p,n,m=0){
  h<-(n-1)/2
  xp<-c(tail(x,h),x,head(x,h))
  yp<-sgolayfilt(xp,p=p,n=n,m=m)
  yp[(h+1):(length(yp)-h)]
}

# Prepare curve b for smoothing
x<-as.numeric(b[1,])
y<-as.numeric(b[2,])
n<-length(x)

# Insert midpoints for better resolution
xm<-(x[-n]+x[-1])/2
ym<-(y[-n]+y[-1])/2
xi<-c(rbind(x[-n],xm),x[n])
yi<-c(rbind(y[-n],ym),y[n])
N<-length(xi)
t<-seq(0,1,length.out=N)
dt<-t[2]-t[1]

# Smoothing parameters
p<-4
w<-21

# Apply smoothing
xs<-sgolay_p(xi,p,w)
ys<-sgolay_p(yi,p,w)

# Enforce closed curve
mx<-(xs[1]+xs[N])/2
my<-(ys[1]+ys[N])/2
xs[c(1,N)]<-mx
ys[c(1,N)]<-my

# Remove drift
dxe<-xs[N]-xs[1]
dye<-ys[N]-ys[1]
xs<-xs-t*dxe
ys<-ys-t*dye

# First and second derivatives
dx<-sgolay_p(xs,p,w,1)/dt
dy<-sgolay_p(ys,p,w,1)/dt
ddx<-sgolay_p(xs,p,w,2)/dt^2
ddy<-sgolay_p(ys,p,w,2)/dt^2

# Curvature from smoothed curve
num<-dx*ddy-dy*ddx
den<-(dx^2+dy^2)^(3/2)
k<-ifelse(den>1e-8,num/den,NA)
s<-sqrt(dx^2+dy^2)

# Plot original vs smoothed curve
p_curve<-ggplot()+geom_path(aes(x=b[1,],y=b[2,],color="Original"),linewidth=1.4)+
  geom_path(aes(x=xs,y=ys,color="Smoothed"),linewidth=1.4,linetype="dashed")+
  scale_color_manual(values=c("Original"="red","Smoothed"="blue"))+
  coord_equal()+labs(title="",x="x",y="y",color="")+
  theme(legend.position=c(0.65, 0.95),legend.justification=c(0,1),
    legend.background=element_rect(fill="white", color=NA))

# Raw curvature (finite differences)
ip1<-c(2:N,1)
im1<-c(N,1:(N-1))
dxr<-(xi[ip1]-xi[im1])/(2*dt)
dyr<-(yi[ip1]-yi[im1])/(2*dt)
ddxr<-(xi[ip1]-2*xi+xi[im1])/(dt^2)
ddyr<-(yi[ip1]-2*yi+yi[im1])/(dt^2)
numr<-dxr*ddyr-dyr*ddxr
denr<-(dxr^2+dyr^2)^(3/2)
kr<-ifelse(denr>1e-8,numr/denr,NA)
dfk<-data.frame(t=t,raw=kr,smooth=k)

# Plot curvature comparison
p_k<-ggplot(dfk,aes(t))+geom_line(aes(y=raw,color="Raw"),linewidth=0.8)+
  geom_line(aes(y=smooth,color="Smoothed"),linewidth=1.2)+
  scale_color_manual(values=c("Raw"="red","Smoothed"="blue"))+
  labs(title="",x="t",y=expression(kappa),color="")+
  theme(legend.position=c(0.05, 0.95),legend.justification=c(0,1),
    legend.background=element_rect(fill="white",color=NA))

# Show plots
pb0
pc0
p_curve
p_k
