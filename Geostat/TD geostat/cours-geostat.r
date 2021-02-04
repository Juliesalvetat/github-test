f.var<-function(x,na.rm=T){
  mean(x^2,na.rm=na.rm)-mean(x,na.rm=na.rm)^2
}
f.cov<-function(x,y,na.rm=T){
  mean(x*y,na.rm=na.rm)-mean(x,na.rm=na.rm)*mean(y,na.rm=na.rm)
}

# creer une grille 2D 
ng<-100
nk1<-10 #fenetre pour moyenne glissante 
nk2<-10
n1<-ng+2*nk1
n2<-ng+2*nk2
set.seed(1234)
z0<-list(x=seq(-nk1,(ng-1)+nk1),
         y=seq(-nk2,(ng-1)+nk2),
               z=matrix(rnorm(n1*n2),ncol=n2,nrow=n1))
image(z0,asp=1) # asp=1 axes x et y de meme taille 

#moyenne glissante 
''
Z<-z0
Z$z<-Z$z*NA # remplace les données par des NA
for (i in (1+nk1):(nk1+ng)){
  for (j in (1+nk2):(nk2+ng)){
    Z$z[i,j]<-mean(z0$z[(i-nk1):(i+nk1),(j-nk2):(j+nk2)])
  }
}
image(Z,asp=1)

#######
nk<-10
n1<- ng + 2*nk
n2<- ng + 2*nk
range<-10
x<-expand.grid(-nk:nk,-nk:nk)
D<-matrix(sqrt(x[,1]^2+x[,2]^2),ncol=1+2*nk,nrow=1+2*nk)
#exponentielle
kern.exp<-exp(-D/(range/3))
#gaussienne
kern.gauss<- exp(-D^2/range)

Z.exp <- z0
Z.exp$z<-Z.exp$z*NA
Z.gauss <- Z.exp

for (i in (1+nk):(nk+ng)){
  for (j in (1+nk):(nk+ng)){
    temp<-z0$z[(i-nk):(i+nk),(j-nk):(j+nk)]
    Z.exp$z[i,j]<-mean(temp*kern.exp)
    Z.gauss$z[i,j]<-mean(temp*kern.gauss)
  }
}

par(mfrow=c(1,2))
image(Z.exp,main='Z.exp')
image(Z.gauss,main='Z.gauss')
hist(Z)
qqnorm(Z.exp$z)
plot(Z.exp$z[(1+1):n1,1:n2],Z.exp$z[1:(n1-1),1:n2])
plot(z0$z[(1+1):n1,1:n2],z0$z[1:(n1-1),1:n2])
plot(Z.exp$z[(1+2):n1,1:n2],Z.exp$z[1:(n1-2),1:n2])

f.grid.cov.spatiale<-function(Z,nlag=10){
  res.x<-rep(NA,nlag+1);res.y<-res.x
  nx<-dim(Z$z)[1];ny<-dim(Z$z)[2]
  for (i in 0:nlag){
    res.x[i+1]<-f.cov(Z$z[(1+i):nx,1:ny],Z$z[1:(nx-i),1:ny])
    res.y[i+1]<-f.cov(Z$z[1:nx,(1+i):ny],Z$z[1:nx,1:(ny-i)])
  }
  list(h=(0:nlag)*(Z$x[2]-Z$x[1]),cov=(res.x+res.y)/2,cov.x=res.x,cov.y=res.y)
}

cov.exp<-f.grid.cov.spatiale(Z=Z.exp,nlag=80)
plot(cov.exp$h,cov.exp$cov.x,type='l',ylim=range(c(cov.exp$cov.x,cov.exp$cov.y)),xaxs='i')
lines(cov.exp$h,cov.exp$cov.y,col="red")
abline(h=0)
#points(x=c(1,4,10),y=c(cov.1,cov.4,cov.10))

cov.gauss<-f.grid.cov.spatiale(Z=Z.gauss,nlag=80)
lines(cov.gauss$h,cov.gauss$cov.x,col="black",lty=2)
lines(cov.gauss$h,cov.gauss$cov.y,col="red",lty=2)

range.gauss<-range/1.73
kern.gauss<-exp(-2*D^2/range.gauss^2)
range.exp<-range/3
kern.exp<-(D/range.exp)^(-1/4*besselK(x=D/range.exp,nu=1/4))
#bidouille
kern.exp[nk+1,nk+1]<-(0.1)^(-1/4)*besselK(0.01,nu=1/4)


Z.exp <- z0
Z.exp$z<-Z.exp$z*NA

for (i in (1+nk):(nk+ng)){
  for (j in (1+nk):(nk+ng)){
    temp<-z0$z[(i-nk):(i+nk),(j-nk):(j+nk)]
    Z.exp$z[i,j]<-mean(temp*kern.exp)
  }
}

par(mfrow=c(1,1))
x11()
image(Z.exp,main='Z.exp')

cov.exp<-f.grid.cov.spatiale(Z=Z.exp,nlag=80)
plot(cov.exp$h,cov.exp$cov.x,type='l',ylim=range(c(cov.exp$cov.x,cov.exp$cov.y)),xaxs='i')
lines(cov.exp$h,cov.exp$cov.y,col="red")
abline(h=0)

Z.exp$z<-Z.exp$z/sqrt(f.var(Z.exp$z))
Z.gauss$z<-Z.gauss$z/sqrt(f.var(Z.gauss$z))

par(mfrow=c(1,2))
image(Z.exp,main='Z.exp')
image(Z.gauss,main='Z.gauss')

cov.exp<-f.grid.cov.spatiale(Z=Z.exp,nlag=80)
x11()
plot(cov.exp$h,cov.exp$cov.x,type='l',ylim=range(c(cov.exp$cov.x,cov.exp$cov.y)),xaxs='i')
lines(cov.exp$h,cov.exp$cov.y,col="red")
abline(h=0)
cov.gauss<-f.grid.cov.spatiale(Z=Z.gauss,nlag=80)
lines(cov.gauss$h,cov.gauss$cov.x,col="black",lty=2)
lines(cov.gauss$h,cov.gauss$cov.y,col="red",lty=2)
############## 2ième jour
par(mfrow=c(1,3))
n.realisation<-50
range <- 30
range.exp <- range/3
nk<-range
ng<-100
n1<- ng + 2*nk
n2<- ng + 2*nk
x<-expand.grid(-nk:nk,-nk:nk)
D<-matrix(sqrt(x[,1]^2+x[,2]^2),ncol=1+2*nk,nrow=1+2*nk)
kern.exp<-(D/range.exp)^(-1/4)*besselK(x=D/range.exp,nu=1/4)
#bidouille
kern.exp[nk+1,nk+1]<-(0.1)^(-1/4)*besselK(0.01,nu=1/4)

res.z<-array(NA,dim=c(n1,n2,n.realisation))
for(k in 1:n.realisation){
  set.seed(k)
  z0<-matrix(rnorm(n1*n2),ncol=n2,nrow=n1)
  for (i in (1+nk):(nk+ng)){
    for (j in (1+nk):(nk+ng)){
      temp<-z0[(i-nk):(i+nk),(j-nk):(j+nk)]
      res.z[i,j,k]<-weighted.mean(temp,kern.exp)    
    }
  }
}

#image(res.z[,,1])
#par(mfrow=c(1,2))
#hist(res.z[37,42,],nclass=20)
#hist(res.z[30,50,],nclass=20)

#par(mfrow=c(1,1))
#plot(res.z[37,42,],res.z[30,50,])
#plot(res.z[37,42,],res.z[37,43,])
#abline(lm(res.z[37,42,]~res.z[37,43,]))
#toto<-lm(res.z[37,42,]~res.z[37,43,])
#cov
#toto$coefficients[2]*f.var(res.z[37,42,])


cov.simulee<-NULL
for(k in 1:n.realisation){
  temp<-f.grid.cov.spatiale(list(x=seq(-nk,(ng-1)+nk),y=seq(-nk,(ng-1)+nk),z=res.z[,,k]),nlag=40)
  cov.simulee<-cbind(cov.simulee,temp$cov)
}

cov.h.1.simulee<-cov.simulee[2,]
mean(cov.h.1.simulee)
plot(c(0,40),range(cov.simulee),type="n",xaxs='i')
abline(h=0)
for(k in 1:n.realisation){
  lines(0:40,cov.simulee[,k],col=1)
}


n.realisation<-50
range <- 10
range.exp <- range/3
nk<-range
ng<-100
n1<- ng + 2*nk
n2<- ng + 2*nk
x<-expand.grid(-nk:nk,-nk:nk)
D<-matrix(sqrt(x[,1]^2+x[,2]^2),ncol=1+2*nk,nrow=1+2*nk)
kern.exp<-(D/range.exp)^(-1/4)*besselK(x=D/range.exp,nu=1/4)
#bidouille
kern.exp[nk+1,nk+1]<-(0.1)^(-1/4)*besselK(0.01,nu=1/4)

res.z<-array(NA,dim=c(n1,n2,n.realisation))
for(k in 1:n.realisation){
  set.seed(k)
  z0<-matrix(rnorm(n1*n2),ncol=n2,nrow=n1)
  for (i in (1+nk):(nk+ng)){
    for (j in (1+nk):(nk+ng)){
      temp<-z0[(i-nk):(i+nk),(j-nk):(j+nk)]
      res.z[i,j,k]<-weighted.mean(temp,kern.exp)    
    }
  }
}

#image(res.z[,,1])
#par(mfrow=c(1,2))
#hist(res.z[37,42,],nclass=20)
#hist(res.z[30,50,],nclass=20)

#par(mfrow=c(1,1))
#plot(res.z[37,42,],res.z[30,50,])
#plot(res.z[37,42,],res.z[37,43,])
#abline(lm(res.z[37,42,]~res.z[37,43,]))
#toto<-lm(res.z[37,42,]~res.z[37,43,])
#cov
#toto$coefficients[2]*f.var(res.z[37,42,])


cov.simulee<-NULL
for(k in 1:n.realisation){
  temp<-f.grid.cov.spatiale(list(x=seq(-nk,(ng-1)+nk),y=seq(-nk,(ng-1)+nk),z=res.z[,,k]),nlag=40)
  cov.simulee<-cbind(cov.simulee,temp$cov)
}

cov.h.1.simulee<-cov.simulee[2,]
mean(cov.h.1.simulee)
plot(c(0,40),range(cov.simulee),type="n",xaxs='i')
abline(h=0)
for(k in 1:n.realisation){
  lines(0:40,cov.simulee[,k],col=1)
}

n.realisation<-50
range <- 5
range.exp <- range/3
nk<-range
ng<-100
n1<- ng + 2*nk
n2<- ng + 2*nk
x<-expand.grid(-nk:nk,-nk:nk)
D<-matrix(sqrt(x[,1]^2+x[,2]^2),ncol=1+2*nk,nrow=1+2*nk)
kern.exp<-(D/range.exp)^(-1/4)*besselK(x=D/range.exp,nu=1/4)
#bidouille
kern.exp[nk+1,nk+1]<-(0.1)^(-1/4)*besselK(0.01,nu=1/4)

res.z<-array(NA,dim=c(n1,n2,n.realisation))
for(k in 1:n.realisation){
  set.seed(k)
  z0<-matrix(rnorm(n1*n2),ncol=n2,nrow=n1)
  for (i in (1+nk):(nk+ng)){
    for (j in (1+nk):(nk+ng)){
      temp<-z0[(i-nk):(i+nk),(j-nk):(j+nk)]
      res.z[i,j,k]<-weighted.mean(temp,kern.exp)    
    }
  }
}

#image(res.z[,,1])
#par(mfrow=c(1,2))
#hist(res.z[37,42,],nclass=20)
#hist(res.z[30,50,],nclass=20)

#par(mfrow=c(1,1))
#
plot(res.z[37,42,],res.z[30,50,])
#plot(res.z[37,42,],res.z[37,43,])
#abline(lm(res.z[37,42,]~res.z[37,43,]))
#toto<-lm(res.z[37,42,]~res.z[37,43,])
#cov
#toto$coefficients[2]*f.var(res.z[37,42,])


cov.simulee<-NULL
for(k in 1:n.realisation){
  temp<-f.grid.cov.spatiale(list(x=seq(-nk,(ng-1)+nk),y=seq(-nk,(ng-1)+nk),z=res.z[,,k]),nlag=40)
  cov.simulee<-cbind(cov.simulee,temp$cov)
}

cov.h.1.simulee<-cov.simulee[2,]
mean(cov.h.1.simulee)
plot(c(0,40),range(cov.simulee),type="n",xaxs='i')
abline(h=0)
for(k in 1:n.realisation){
  lines(0:40,cov.simulee[,k],col=1)
$$}