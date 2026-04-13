# This script looks at curves in different ways:
# Cartesian vs polar coordinates, velocity, acceleration,
# and curvature. Each is then visualised.

library(ggplot2)

# Parametric definition of a circle
t <- seq(0,2*pi,length.out=1000)
x <- cos(t)
y <- sin(t)

# Approximate velocity and acceleration using finite differences
dt <- t[2]-t[1]
vx <- c(NA,diff(x)/dt)
vy <- c(NA,diff(y)/dt)
ax <- c(NA,diff(vx)/dt)
ay <- c(NA,diff(vy)/dt)

# Remove initial NA values from differentiation
valid_idx <- 3:length(t)
t <- t[valid_idx]
x <- x[valid_idx]
y <- y[valid_idx]
vx <- vx[valid_idx]
vy <- vy[valid_idx]
ax <- ax[valid_idx]
ay <- ay[valid_idx]

# Cartesian coordinates plot
df_cartesian <- data.frame(t=t,x=x,y=y)
plot1 <- ggplot(df_cartesian,aes(x=x,y =y))+geom_path(color="red",size=3)+
  coord_equal()+labs(x="x",y ="y")+theme_minimal(base_size=14)+
  theme(axis.title=element_text(size=34,face="bold"),
        axis.text=element_text(size=30))

# Convert to polar coordinates
r <- sqrt((df_cartesian$x)^2 + (df_cartesian$y)^2)
theta <- atan2(df_cartesian$y,df_cartesian$x) 
df_polar <- data.frame(theta=theta,r=r)

# Polar plot
plot2 <- ggplot(df_polar, aes(theta,r))+geom_path(color="blue",size=3)+
  labs(x="θ",y="r")+scale_x_continuous(breaks=c(-pi,-pi/2,0,pi/2,pi),
                                       labels=c(expression(-pi), expression(-pi/2),"0",
                                                expression(pi/2), expression(pi)))+
  scale_y_continuous(limits=c(0,2),breaks=0:4)+theme_minimal(base_size=14)+
  theme(axis.title=element_text(size=34,face="bold"),
        axis.text=element_text(size=30))

# Velocity vectors along the curve
df_velocity <- data.frame(x=x,y=y,vx=vx,vy=vy)
df_arrows <- df_velocity[seq(1,nrow(df_velocity),by=50), ]
plot3 <- ggplot(df_velocity,aes(x=x,y=y))+geom_path(color="blue",linewidth=1.5)+
  geom_segment(data=df_arrows,aes(xend=x+0.5*vx,yend=y+0.5*vy),
               arrow=arrow(length=unit(0.1,"inches")),color="red",linewidth=1.3)+
  coord_equal()+labs(x="x",y ="y")+theme_minimal(base_size=14)+
  theme(axis.title=element_text(face="bold"),axis.text=element_text(siz =16))

# Acceleration vectors along the curve
df_acceleration <- data.frame(x=x,y=y,ax=ax,ay=ay)
df_arrows <- df_acceleration[seq(1,nrow(df_acceleration),by=50), ]
plot4 <- ggplot(df_acceleration,aes(x=x,y=y))+
  geom_path(color="blue",linewidth=1.5)+geom_segment(data=df_arrows,
                                                     aes(xend=x+0.3*ax,yend=y+0.3*ay),arrow=arrow(length=unit(0.1,"inches")),
                                                     color="red",linewidth=1.3)+coord_equal()+labs(x="x",y ="y")+
  theme_minimal(base_size=14)+theme(axis.title=element_text(face="bold"),
                                    axis.text=element_text(siz =16))

# Curvature calculation from velocity and acceleration
curvature <- abs(vx*ay-vy*ax)/((vx^2+vy^2)^(3/2))
curvature[is.na(curvature)] <- 1.3
df_velocity$curvature <- curvature
range(df_velocity$curvature)

# Plot coloured by curvature
plot5 <- ggplot(df_velocity,aes(x=x,y=y,color=curvature))+
  geom_path(linewidth=2.5)+coord_equal()+
  scale_color_viridis_c(limits=c(0,1.3),oob=scales::squish)+
  labs(x ="x",y ="y",color ="Curvature")+theme_minimal(base_size=14)+
  theme(axis.title=element_text(size=24,face="bold"),
        axis.text=element_text(size=20),legend.position="right",
        legend.key.height=unit(3,"lines"),
        legend.title=element_text(size=24,face="bold"),
        legend.text=element_text(size=20))

