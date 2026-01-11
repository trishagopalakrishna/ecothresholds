require(mgcv)
test1 <- function(x,z,sx=0.3,sz=0.4) { 
  x <- x*20
  (pi**sx*sz)*(1.2*exp(-(x-0.2)^2/sx^2-(z-0.3)^2/sz^2)+
                 0.8*exp(-(x-0.7)^2/sx^2-(z-0.8)^2/sz^2))
}
n <- 500
#old.par <- par(mfrow=c(2,2))
x <- runif(n)/20;z <- runif(n);
xs <- seq(0,1,length=30)/20;zs <- seq(0,1,length=30)
pr <- data.frame(x=rep(xs,30),z=rep(zs,rep(30,30)))
truth <- matrix(test1(pr$x,pr$z),30,30)
f <- test1(x,z)
y <- f + rnorm(n)*0.2
b1 <- gam(y~s(x,z))
persp(xs,zs,truth);title("truth")
vis.gam(b1);title("t.p.r.s")
b2 <- gam(y~te(x,z))
vis.gam(b2);title("tensor product")


b3 <- gam(y ~ ti(x) + ti(z) + ti(x,z)) ##############################
vis.gam(b3);title("tensor anova")


xgrid <- seq(min(x),max(x),len=200)
zgrid <- seq(min(z),max(z),len=3)

newd <- data.frame( expand.grid(xgrid,zgrid) )

names(newd) <- c("x","z")


pred <- predict(b3,newdata=newd,type="terms")



par(mfrow = c(2, 2))
ind <- 1:200
x1Eff <- rowSums( pred[ind, c(1,3)] )
plot(xgrid,x1Eff, type="l", ylim=c(-0.5,0.5))
ind <- 201:400
x2Eff <- rowSums( pred[ind, c(1,3)] )
lines(xgrid,x2Eff,col="red")
ind <- 401:600
x3Eff <- rowSums( pred[ind, c(1,3)] )
lines(xgrid,x3Eff,col="blue")


# with uncertainty
X <- predict(b3,newdata=newd,type="lpmatrix")


# beta uncertainty
n.sims <- 1000
b_sims <- rmvn(n.sims,coef(b3),b3$Vp)



# pick x1 "effects"
coef(b3)
coef_index <- grep("x",names(coef(b3)),fixed = T)


# compute the spline
x4Eff <- tcrossprod( X[,coef_index] , b_sims[,coef_index] )
dim(x1Eff) #600 values *1000 simulation= 600000



#

plot(xgrid, apply(x4Eff[1:200,],1,mean),type="l", lty = 3, ylim = c(-0.5,0.5))
lines(xgrid, apply(x4Eff[1:200,], 1, quantile,probs = 0.025), type = "l", lty = 3)
lines(xgrid, apply(x4Eff[1:200,], 1, quantile,probs = 0.975), type = "l", lty = 3)


plot(xgrid, apply(x4Eff[201:400,],1,mean),type="l", lty = 3, ylim = c(-0.5,0.5))
lines(xgrid, apply(x4Eff[201:400,], 1, quantile,probs = 0.025), type = "l", lty = 3)
lines(xgrid, apply(x4Eff[201:400,], 1, quantile,probs = 0.975), type = "l", lty = 3)


plot(xgrid, apply(x4Eff[401:600,],1,mean),type="l")
lines(xgrid, apply(x4Eff[401:600,], 1, quantile,probs = 0.025), type = "l", lty = 3)
lines(xgrid, apply(x4Eff[401:600,], 1, quantile,probs = 0.975), type = "l", lty = 3)

 