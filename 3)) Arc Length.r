# This script reparameterises a curve by arc length to compare 
# how the motion changes. Both versions are visualised.

library(ggplot2)
library(dplyr)
library(pracma)

# Define curve
t <- seq(-2,2,length.out=1000)
x <- t
y <- tanh(2*t)
dt <- t[2]-t[1]

# Compute velocity and speed
vx <- c(NA, diff(x)/dt)
vy <- c(NA, diff(y)/dt)
speed <- sqrt(vx^2+vy^2)

# Remove NA values from differentiation
valid <- 3:length(t)
t <- t[valid]
x <- x[valid]
y <- y[valid]
vx <- vx[valid]
vy <- vy[valid]
speed <- speed[valid]

df_main <- data.frame(t,x,y,vx,vy,speed)

# Sample points for velocity arrows
df_arrow <- df_main[seq(1,nrow(df_main),by=50), ] |> na.omit()

# Plot original parameterisation with velocity vectors
plot6 <- ggplot(df_main,aes(x,y))+geom_path(color="blue",linewidth=2)+
  geom_segment(data=df_arrow,aes(xend=x+vx,yend=y+vy),
  arrow=arrow(length=unit(0.1,"inches")),color="red",linewidth=1.4)+
  coord_equal()+labs(x="x",y="y")+theme_minimal(base_size=14)+
  theme(axis.title=element_text(size=26,face="bold"),
  axis.text=element_text(size=22))
# Compute arc length using numerical integration
s <- cumtrapz(t,speed)
s <- s-min(s)

# Create evenly spaced arc length values
s_even <- seq(min(s),max(s),length.out=1000)

# Map arc length back to parameter t
t_from_s <- approxfun(s,t)
t_new <- t_from_s(s_even)

# Reconstruct curve using arc length parameter
x_new <- t_new
y_new <- tanh(2*t_new)

ds <- s_even[2]-s_even[1]

# Velocity under arc length parameterisation
vx_new <- c(NA,diff(x_new)/ds)
vy_new <- c(NA,diff(y_new)/ds)
speed_new <- sqrt(vx_new^2+vy_new^2)

df_reparam <- data.frame(s=s_even,x=x_new,y=y_new,
  vx=vx_new,vy=vy_new,speed=speed_new)

# Sample points for arrows
df_sarrow <- df_reparam[seq(1,nrow(df_reparam),by=50), ] |> na.omit()

# Plot arc length parameterisation
plot7 <- ggplot(df_reparam,aes(x,y))+geom_path(color="blue",linewidth=2)+
  geom_segment(data=df_sarrow,aes(xend=x+0.7*vx,yend=y+0.7*vy),
  arrow=arrow(length=unit(0.1,"inches")),color="red",linewidth=1.4)+coord_equal()+
  ylim(-1,2)+labs(x="x",y="y")+theme_minimal(base_size=14)+
  theme(axis.title=element_text(size=26,face="bold"),
  axis.text=element_text(size=22))
