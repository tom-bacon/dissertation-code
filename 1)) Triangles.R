# This script generates random acute and obtuse triangles,
# applies transformations and classifies them using LDA.
# Results are then visualised.

library(MASS)
library(ggplot2)
library(patchwork)

# Calculate the three angles of a triangle
angles=function(p){
  d=as.matrix(dist(p))
  a=d[1,2];b=d[2,3];c=d[1,3]
  # Cosine rule
  cA=(b^2+c^2-a^2)/(2*b*c)
  cB=(a^2+c^2-b^2)/(2*a*c)
  cC=(a^2+b^2-c^2)/(2*a*b)
  acos(pmin(pmax(c(cA,cB,cC),-1),1))
}

# Generate either an acute ("a") or obtuse ("o") triangle
tri=function(type){
  repeat{
    p=matrix(runif(6,-1,1),3,2)
    ang=angles(p)
    d=as.matrix(dist(p))
    edges=d[d>0]
    # Avoid very small or degenerate triangles
    if(min(edges) < 0.25) next
    if(min(ang) < 0.2) next
    if(type=="a"&&all(ang<pi/2)) return(p)
    if(type=="o"&&any(ang>pi/2)) return(p)
  }
}

# Apply random rotation, scaling, translation
move=function(p){
  t=runif(1,0,2*pi)
  R=matrix(c(cos(t),-sin(t),sin(t),cos(t)),2,2)
  p=t(R%*%t(p))
  p=p*runif(1,0.5,3)
  p=sweep(p,2,runif(2,-5,5),"+")
  p[sample(1:3),]  # Shuffle point order
}

# Create dataset of triangles
make=function(n){
  curves=list();y=c()
  for(i in 1:n){curves[[i]]=move(tri("a"));y[i]="acute"}
  for(i in 1:n){curves[[n+i]]=move(tri("o"));y[n+i]="obtuse"}
  list(curves=curves,y=factor(y))
}

# Flatten triangle coordinates into feature vectors
flat=function(curves){
  t(sapply(curves,function(p)as.vector(t(p))))
}

# Build dataframe for plotting triangles
tri_df=function(curves,y,n=50){
  id=sample(1:length(curves),n)
  df=data.frame()
  for(i in id){
    p=curves[[i]]
    p=rbind(p,p[1,])  # Close triangle
    tmp=data.frame(x=p[,1],y=p[,2],group=i,class=y[i])
    df=rbind(df,tmp)
  }
  df
}

# Generate data
d=make(150)
X=flat(d$curves)
y=d$y

# Train/test split
n=nrow(X)
id=sample(1:n,0.7*n)

# Fit LDA model and predict
m=lda(X[id,],grouping=y[id])
pred=predict(m,X[-id,])$class

# Accuracy + confusion matrix
acc=mean(pred==y[-id])
print(acc)
print(table(y[-id],pred))

# Fit on full data for visualisation
m_all=lda(X,grouping=y)
z=predict(m_all)$x[,1]
df1=tri_df(d$curves,y)

# Plot triangles
p1=ggplot(df1,aes(x,y,group=group,color=class))+geom_path(linewidth=1)+
  scale_color_manual(values=c("acute"="blue","obtuse"="red"),
  labels=c("Acute","Obtuse"),name="Class")+coord_equal()+
  theme_classic(base_size=14)+labs(x="x",y="y")+theme(legend.position="none",
  axis.title=element_text(size=16,face="bold"),axis.text=element_text(size=12))

# LDA projection plot
df2=data.frame(ld=z,y=jitter(rep(0,length(z)),0.05),class=y)
p2=ggplot(df2,aes(ld,y,color=class))+geom_point(size=2.5)+
  scale_color_manual(values=c("acute"="blue","obtuse"="red"),
  labels=c("Acute","Obtuse"),name="Class")+theme_classic(base_size=14)+
  labs(x="LDA Score",y="")+theme(legend.position="right",
  axis.title=element_text(size=16,face="bold"),
  axis.text.x=element_text(size=12),axis.text.y=element_blank(),
  axis.ticks.y=element_blank(),legend.title=element_text(face="bold",size=14),
  legend.text=element_text(face="bold",size=12))

# Combine plots
(p1 + p2) + plot_layout(guides="collect")
