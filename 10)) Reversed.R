# This script smooths a curve, computes its curvature, and
# compares the result with the curve traced in reverse.
# Both the shape and curvature are visualised.

library(signal)
library(ggplot2)
library(fdasrvf)
library(pracma)
library(grid)

# Load curve data
data("beta")
b<-beta[,,1,1]
x<-as.numeric(b[1,])
y<-as.numeric(b[2,])

# Insert midpoints for smoother resolution
n<-length(x)
xm<-(x[-n]+x[-1])/2
ym<-(y[-n]+y[-1])/2
xi<-c(rbind(x[-n],xm),x[n])
yi<-c(rbind(y[-n],ym),y[n])

# Savitzky–Golay smoothing 
sgolay_p<-function(x,p,n,m){
  h<-(n-1)/2
  xp<-c(tail(x,h),x,head(x,h))
  yp<-sgolayfilt(xp,p=p,n=n,m=m)
  yp[(h+1):(length(yp)-h)]
}

# Smooth, differentiate, compute curvature
proc<-function(xi,yi,w=21,p=4){
  N<-length(xi)
  t<-seq(0,1,length.out=N)
  dt<-t[2]-t[1]
  xs<-sgolay_p(xi,p,w,0)
  ys<-sgolay_p(yi,p,w,0)
  # Enforce closed curve
  mx<-(xs[1]+xs[N])/2
  my<-(ys[1]+ys[N])/2
  xs[c(1,N)]<-mx
  ys[c(1,N)]<-my
  # Remove drift
  dx<-xs[N]-xs[1]
  dy<-ys[N]-ys[1]
  xs<-xs-t*dx
  ys<-ys-t*dy
  # Derivatives
  dxs<-sgolay_p(xs,p,w,1)/dt
  dys<-sgolay_p(ys,p,w,1)/dt
  ddxs<-sgolay_p(xs,p,w,2)/dt^2
  ddys<-sgolay_p(ys,p,w,2)/dt^2
  # Curvature
  num<-dxs*ddys-dys*ddxs
  den<-(dxs^2+dys^2)^(3/2)
  k<-ifelse(den>1e-8 & is.finite(den),num/den,NA)
  s<-sqrt(dxs^2+dys^2)
  df<-data.frame(t=t,x=xs,y=ys,k=k)
  list(t=t,x=xs,y=ys,dx=dxs,dy=dys,s=s,k=k,df=df)
}

# Plot curve with curvature colouring and starting direction
p1<-function(o,L,a=70,ttl=""){
  i<-1
  tx<-o$dx[i]/o$s[i]
  ty<-o$dy[i]/o$s[i]
  ggplot(o$df,aes(x=x,y=y,color=k))+geom_path(linewidth=2)+
    geom_point(data=o$df[i,],aes(x=x,y=y),inherit.aes=FALSE,shape=21,size=2.2,
    fill="green",color="green",stroke=0.8)+geom_text(data=o$df[i,],
    aes(x=x,y=y,label="Start"),inherit.aes=FALSE,color="green",size=4,
    nudge_x=8,nudge_y=16,fontface="bold")+geom_segment(data=data.frame(
    x=o$df$x[i],y=o$df$y[i],xend=o$df$x[i]+a*tx,yend=o$df$y[i]+a*ty),
    aes(x=x,y=y,xend=xend,yend=yend),inherit.aes=FALSE,color="green",
    linewidth=1.2,arrow=arrow(type="closed",length=unit(0.25,"cm")))+
    scale_color_gradient2(low="blue",mid="lightgrey",high="red",
    midpoint=0,limits=c(-L,L),oob=scales::squish)+coord_equal()+
    theme_minimal(base_size=16)+labs(x="x",y="y",color=expression(kappa),
    title=ttl)+theme(axis.title=element_text(face="bold"),
    axis.text=element_text(color="black"),legend.key.height=unit(1.5,"cm"),
    legend.title=element_text(face="bold"))
}

# Plot curvature over parameter t
p2<-function(o,ttl=""){
  ggplot(o$df,aes(x=t,y=k))+geom_line(linewidth=1.5)+theme_minimal(base_size=16)+
  labs(x="t",y=expression(kappa),title=ttl)+theme(axis.title=element_text(
  face="bold"),axis.text=element_text(color="black"))
}

# Parameters
w<-21
p<-4

# Original and reversed curve
o<-proc(xi,yi,w,p)
r<-proc(rev(xi),rev(yi),w,p)

# Common colour scale limit
L<-quantile(abs(c(o$df$k,r$df$k)),0.9,na.rm=TRUE)

# Generate plots
p1o<-p1(o,L,70,"")
p2o<-p2(o,"")
p1r<-p1(r,L,70)
p2r<-p2(r)

# Show results
p1o
p2o
p1r
p2r
