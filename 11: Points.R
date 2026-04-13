# This script smooths a curve, computes curvature and identifies
# maxima, minima, and inflection points, which are then visualised.

library(signal)
library(ggplot2)
library(pracma)
library(grid)

data("beta")
b<-beta[,,1,1]

# Savitzky–Golay filter
f<-function(x,p,n,m=0){
  h<-(n-1)/2
  xp<-c(tail(x,h),x,head(x,h))
  yp<-sgolayfilt(xp,p=p,n=n,m=m)
  yp[(h+1):(length(yp)-h)]
}

# Prepare curve
x<-as.numeric(b[1,])
y<-as.numeric(b[2,])
n<-length(x)

# Add midpoints for smoother sampling
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

# Smooth curve
xs<-f(xi,p,w)
ys<-f(yi,p,w)

# Enforce closed curve
mx<-(xs[1]+xs[N])/2
my<-(ys[1]+ys[N])/2
xs[c(1,N)]<-mx
ys[c(1,N)]<-my

# Derivatives
dx<-f(xs,p,w,1)/dt
dy<-f(ys,p,w,1)/dt
ddx<-f(xs,p,w,2)/dt^2
ddy<-f(ys,p,w,2)/dt^2

# Curvature
num<-dx*ddy-dy*ddx
den<-(dx^2+dy^2)^(3/2)
k<-ifelse(den>1e-8,num/den,0)

df<-data.frame(x=xs,y=ys,k=k,t=t)

# Find maxima and minima
rng<-max(abs(k))
thr<-0.15*rng
pk<-findpeaks(k,minpeakheight=thr,minpeakdistance=8)
vk<-findpeaks(-k,minpeakheight=thr,minpeakdistance=8)
imx<-if(!is.null(pk)) pk[,2] else integer(0)
imn<-if(!is.null(vk)) vk[,2] else integer(0)
dfx<-df[imx,]
dfn<-df[imn,]

# Inflection points 
sg<-sign(k)
chg<-which(sg[-length(sg)]!=sg[-1])

# Filter weak changes
keep<-function(id,k){
  if(length(id)==0) return(id)
  g<-abs(diff(k))
  id[g[id]>0.15*max(g)&(abs(k[id])+abs(k[id+1]))>0.1*max(abs(k))]
}
chg<-keep(chg,k)
dfi<-df[chg,]

# Plot limits
L<-quantile(abs(k),0.98)
ylim<-c(-L,1.3*L)

# Plot style
sty<-theme_minimal()+theme(axis.title=element_text(size=24,face="bold"),
  axis.text=element_text(size=20),legend.title=element_text(size=24,face="bold"),
  legend.text=element_text(size=20))

# Curve with extrema
p_ext<-ggplot(df,aes(x,y,color=k))+geom_path(linewidth=2)+
  geom_point(data=dfx,aes(x,y),color="green",size=3)+
  geom_point(data=dfn,aes(x,y),color="red",size=3)+
  scale_color_gradient2(low="blue",mid="grey",high="red",midpoint=0,
  limits=c(-L,L),name=expression(kappa),guide=guide_colorbar(
  barheight=unit(12,"cm"),barwidth=unit(1.2,"cm")))+coord_equal()+sty

# Curve with inflection points
p_inf<-ggplot(df,aes(x,y,color=k))+geom_path(linewidth=2)+
  geom_point(data=dfi,aes(x,y),color="black",size=3)+scale_color_gradient2(
    low="blue",mid="grey",high="red",midpoint=0,limits=c(-L,L),
    name=expression(kappa),guide=guide_colorbar(barheight=unit(12,"cm"),
    barwidth=unit(1.2,"cm")))+coord_equal()+sty

# Curvature data
d2<-data.frame(t=t,k=k)
dfx2<-d2[imx,]
dfn2<-d2[imn,]
dfi2<-data.frame(t=t[chg],k=0)

# Curvature with extrema
p_kext<-ggplot(d2,aes(t,k,color=k))+geom_line(linewidth=1.5)+
  geom_point(data=dfx2,aes(t,k),color="green",size=3)+
  geom_point(data=dfn2,aes(t,k),color="red",size=3)+
  scale_color_gradient2(low="blue",mid="grey",high="red",midpoint=0,
  limits=c(-L,L),name=expression(kappa),guide=guide_colorbar(
  barheight=unit(12,"cm"),barwidth=unit(1.2,"cm")))+
  coord_cartesian(ylim=ylim)+labs(y=expression(kappa))+ sty

# Curvature with inflection points
p_kinf<-ggplot(d2,aes(t,k,color=k))+geom_line(linewidth=1.5)+
  geom_point(data=dfi2,aes(t,k),color="black",size=3)+scale_color_gradient2(
  low="blue",mid="grey",high="red",midpoint=0,limits=c(-L,L),name=expression(kappa),
  guide=guide_colorbar(barheight=unit(12,"cm"),barwidth=unit(1.2,"cm")))+
  coord_cartesian(ylim=ylim)+labs(y=expression(kappa))+sty

# Show plots
p_ext
p_inf
p_kext
p_kinf

