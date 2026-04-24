#########################
## ST3131 Lecture notes 8
#########################

#######################################################################
### Women Blood pressure example

WBP.dat = read.table("Rsession/WomenBP.txt",header=TRUE)
y=WBP.dat$Pressure; x=WBP.dat$Age.women 


## unweighted fitting
uw.fit = lm(y~x);r = uw.fit$resid;y.h = uw.fit$fitted
plot(y.h,r,xlab="Fitted values",ylab="Residuals")
y.h.uw = uw.fit$fitted

## weighted fitting

wt.fit=lm(log(abs(r))~log(y.h));
s=exp(wt.fit$fitted); w=1/s^2
w.fit=lm(y~x,weight=w); r=w.fit$resid; y.h = w.fit$fitted
r.w=r*sqrt(w); y.w =y.h*sqrt(w)
plot(y.w,r.w,xlab="Weighted fitted values",ylab="Weithed residuals")

summary(uw.fit) 
summary(w.fit)

###########################################################################
### Trial of children with neurological problems
y = c(1.38,1.26,1.51,1.46,1.61,1.59, 1.36,1.28,1.41,1.39,1.51,1.44)
n = c(41,41,33,45,18,23,40,38,35,46,20,23)
s = c(0.22,0.25,0.31,0.28,0.34,0.46,0.28,0.19,0.27,0.33,0.41,0.30)

S=factor(rep(c(1,1,2,2,3,3),2)) #social class
G=factor(rep(c(1,2),6)) # gender
T=factor( c(rep(1,6), rep(2,6))) #treatment

## unweighted fitinng 
fit0=lm(y~T+G+S) #unweighted
r0=fit0$resid; y0=fit0$fitted
plot(y0,r0)

## Fitting with weight n --- constant variance for individuals
w1=n   
fit1=lm(y~T+G+S,weight=w1) 
r1=fit1$resid; y1=fit1$fitted
r1.w=r1*sqrt(w1); y1.w =y1*sqrt(w1)
plot(y1.w,r1.w,xlab="Residual plots with weight =n")


## Fitting with weight n/s^2, non-constant variances for groups
w2=n/s^2 
fit2=lm(y~T+G+S,weight=w2)
r2=fit2$resid; y2=fit2$fitted
r2.w=r2*sqrt(w2); y2.w =y2*sqrt(w2)
plot(y2.w,r2.w,xlab="Residual plots with weight =n/s^2")

## Fitting with functional estimate of weights

plot(y,s)
w.fit = lm(log(s)~log(y)) 
log.s =w.fit$fitted; s3 =exp(log.s)
w3=n/s3^2   #group indivicual variance as function of mean
fit3=lm(y~T+G+S,weight=w3)
r3=fit3$resid; y3=fit3$fitted
r3.w=r3*sqrt(w3); y3.w =y3*sqrt(w3)
plot(y3.w,r3.w,xlab="Residual plots with functional estimated weight")

#######################

v0=round(diag(vcov(fit0)),5)
v1=round(diag(vcov(fit1)),5)
v2=round(diag(vcov(fit2)),5)
v3=round(diag(vcov(fit3)),5)

coefficients =c("Intercept", "T2","G2","S2","S3")
data.frame(coefficients, un.weight =v0, weight1=v1, 
                 weight2=v2, weight3=v3)

rbind(v0,v1,v2,v3)


##############################################################################
## Employee salary example
salary.dat = read.table("Rsession/Salary.txt",header =TRUE)

y=salary.dat$salary
degree =salary.dat$degree
D=factor(degree)
x1=salary.dat$YOE
x2=salary.dat$NSP

wk.dat=data.frame(y,D,x1,x2)

uw.fit=glm(y~D+x1+x2,data=wk.dat)
r = uw.fit$resid
yhat=uw.fit$fitted
#r=abs(r)

par(mfrow=c(2,2))
plot(yhat,r)
plot(degree,r)
plot(x1,r)
plot(x2,r)

hat =	hatvalues(uw.fit,type="diagonal")
  infl = influence(uw.fit, do.coef = FALSE)
rsta = rstandard(uw.fit, infl, type = "pearson")
rstu = rstudent(uw.fit, infl, type = "pearson")
cook.d =cooks.distance(uw.fit, infl,res = infl$pear.res, 
	        dispersion = summary(uw.fit)$dispersion,
	        hat = infl$hat)
par(mfrow=c(2,2))
qqnorm(hat, main="Q-Q plot of hat values")
qqnorm(rstu,main="Q-Q plot of standardized residuals")
qqnorm(cook.d,main="Q-Q plot of cook's distance")

cbind(order(cook.d), cook.d[order(cook.d)], order(hat), hat[order(hat)],order(rstu), rstu[order(rstu)])

##############

wk.dat2 = wk.dat[-43,]

uw.fit=glm(y~D+x1+x2,data=wk.dat2)

hat =	hatvalues(uw.fit,type="diagonal")
	infl = influence(uw.fit, do.coef = FALSE)
rsta = rstandard(uw.fit, infl, type = "pearson")
rstu = rstudent(uw.fit, infl, type = "pearson")
cook.d =cooks.distance(uw.fit, infl,res = infl$pear.res, 
	        dispersion = summary(uw.fit)$dispersion,
	        hat = infl$hat)
par(mfrow=c(2,2))
qqnorm(hat)
qqnorm(rstu)
qqnorm(cook.d)

cbind(order(cook.d), cook.d[order(cook.d)], order(hat), hat[order(hat)],order(rstu), rstu[order(rstu)])


##############

wk.dat3 = wk.dat[-c(3,43),]

uw.fit=glm(y~D+x1+x2,data=wk.dat3)
r = uw.fit$resid
yhat=uw.fit$fitted

weight =lm(log(abs(r))~log(yhat))
s=exp(weight$fitted)
w=1/s^2
w.fit=lm(y~D+x1+x2,weight=w,data=wk.dat3)
r.w=w.fit$resid*sqrt(w)
y.w=w.fit$fitted*sqrt(w)

t.w=lm(abs(r.w)~y.w)
t.uw=lm(abs(r)~yhat)
par(mfrow=c(1,2))
plot(yhat, abs(r))
abline(t.uw)


plot(y.w,abs(r.w))
abline(t.w)

summary(w.fit)


############################################################################################
### Blood pressure example for multicollinearity
bp.dat = read.table("Rsession/bloodpressure.txt",header =TRUE)

X=bp.dat[,-1]
pairs(X)
round(cor(X),3)

fit.full=lm(BP~Age+Weight+BSA+Dur+Pulse+Stress,x=TRUE,data=bp.dat)
summary(fit.full)

fit.p1=lm(BP~Weight+BSA+Dur+Pulse+Stress,data=bp.dat)
summary(fit.p1)


fit.p2=lm(BP~Age+Weight+BSA+Dur+Stress,data=bp.dat)
summary(fit.p2)

## Self defined function vifChen:

vifChen=function(object) {
 ## Input:
  # object -- a fitted lm object.
 ## Output:
  # vif -- variance inflation factors

 X = object$x[,-1]
 V = vcov(object)[-1,-1]
 n = dim(X)[1] 
 sigma=summary(object)$sigma

 v = diag(V)
 S = diag(var(X))*(n-1)
 vif = v*S/sigma^2
 vif
}

### variance infation factor for blood pressure example.

fit.full=lm(BP~Age+Weight+BSA+Dur+Pulse+Stress,x=TRUE,data=bp.dat)

vifChen(fit.full)

regB = lm(BP~Age+BSA+Dur+Pulse+Stress, x=TRUE, data=bp.dat)
vifChen(regB)

### Removing predictor of high VIF as a remedy

fit.full=lm(BP~Age+Weight+BSA+Dur+Pulse+Stress,x=TRUE,data=bp.dat)

vifChen(fit.full)

## Remove Weight

regB = lm(BP~Age+BSA+Dur+Pulse+Stress, x=TRUE, data=bp.dat)
vifChen(regB)

summary(regB)

null.m=glm(BP~1,data=bp.dat)
  forward.fit=stepAIC(null.m,scope = list(upper = ~Age+BSA+Dur+Pulse+Stress, lower = ~1),direction = "forward")
  summary(forward.fit)

## model diagnostics
s.m = glm(BP~Age+BSA+Dur+Pulse+Stress,data=bp.dat)


r = s.m$resid
yhat=s.m$fitted
hat =	hatvalues(s.m,type="diagonal")
	infl = influence(s.m, do.coef = FALSE)
rsta = rstandard(s.m, infl, type = "pearson")
rstu = rstudent(s.m, infl, type = "pearson")
cook.d =cooks.distance(s.m, infl,res = infl$pear.res, 
	        dispersion = summary(s.m)$dispersion,
	        hat = infl$hat)
par(mfrow=c(2,2))
plot(yhat,r,main="Plot of residual vs. fitted value")
qqnorm(hat,main="QQ plot of hat values")
qqnorm(rstu,main="QQ plot of standardized residual")
qqnorm(cook.d,main="QQ plot of cook's distance")




## Ridge regression for blood example.
library(MASS)
regRidge = lm.ridge(BP~Age+Weight+BSA+Dur+Pulse+Stress,data=bp.dat,
                    lambda=seq(0,0.1,0.005))

CV=regRidge$GCV
lambda=regRidge$lambda

par(mfrow = c(1,1))
plot(lambda,CV,type="l")

lambda.best = lambda[order(CV)[1]]
#or
lambda.best=regRidge$lambda[which(regRidge$GCV==min(regRidge$GCV))]

# Final ridge regression fit

ridge.fit = lm.ridge(BP~Age+Weight+BSA+Dur+Pulse+Stress,data=bp.dat,lambda=lambda.best)
round(ridge.fit,4)

ridge.fit$coef
ridge.fit$lambda
ridge.fit$GCV

### PCR #####
library(pls)
regPCR = pcr(BP~Age+Weight+BSA+Dur+Pulse+Stress,
             data=bp.dat, scale=TRUE, validation="CV")
CV = MSEP(regPCR)$val[1,1,-1]
ncomp = 1:length(CV)
plot(ncomp,CV,type="l")
# validationplot(regPCR, val.type="MSEP")
k.best = ncomp[order(CV)[1]]
k.best
pcr.fit = pcr(BP~Age+Weight+BSA+Dur+Pulse+Stress,
              data=bp.dat, scale=TRUE, ncomp=k.best)
coef(pcr.fit, ncomp=k.best)

### End of PCR ####


#############################################################################
# consider the data in R: mtcars with variables 
# mpg,cyl +disp + hp +drat +   wt + qsec+ vs +am +gear+ carb
# Regress  mpg on the ohters

pairs(mtcars)

regLSE = lm(mpg~cyl +disp + hp +drat +   wt + qsec+ vs +am +gear+ carb, x=TRUE, data=mtcars)

# calculate VIF
vifChen(regLSE)

# obviously, there is strong collinearity 

library(MASS)
regRidge = lm.ridge(mpg~cyl +disp + hp +drat +   wt + qsec+ vs +am +gear+ carb, data=mtcars, lambda = seq(1, 50, 0.01))
plot(regRidge$lambda, regRidge$GCV)

# the minimum CV is achieved at 

lambda_best = regRidge$lambda[which(regRidge$GCV==min(regRidge$GCV))]
regRidge_best = lm.ridge(mpg~cyl +disp + hp +drat +   wt + qsec+ vs +am +gear+ carb, data=mtcars, lambda = lambda_best)
regRidge_best



###################################
## senile dimentia patient example
y=c(0.05,0.15,0.35,0.25,0.2,0.05,0.1,0.05,0.3,0.05,0.25,0,0.15,0,0.05,0,0,0.05,0.1)
g=factor(c(rep(1,11),rep(2,8)))

org.fit=lm(y~g)
r.o=org.fit$resid
y.o = org.fit$fitted 

y.t=asin(sqrt(y))
tra.fit=lm(y.t~g) 
r.t=tra.fit$resid
y.t = tra.fit$fitted 

par(mfrow=c(1,2))
plot(y.o,r.o,main="Original fitting")
plot(y.t,r.t,main="Transformed fitting")

summary(org.fit)
summary(tra.fit)

#####################################
## Trial of vaccilation example
y=c(7925,15643,17462,10805,9300,7538,6297,
    3158,3669,5930,5697,8331,11822)
T=factor(c(rep(1,7),rep(2,6)))

org.fit=lm(y~T)
r.o=org.fit$resid
y.o = org.fit$fitted 

y.t=sqrt(y)
tra.fit=lm(y.t~T) 
r.t=tra.fit$resid
y.t = tra.fit$fitted 

par(mfrow=c(1,2))
plot(y.o,r.o,main="Original fitting")
plot(y.t,r.t,main="Transformed fitting")

summary(org.fit)
summary(tra.fit)


###############################
## Miller Lumber Company example

mlc.dat = read.table("D://Rsession/MLC.txt",header =TRUE)

y=mlc.dat$NoC; x1=mlc.dat$HU; x2=mlc.dat$AvI;
x3=mlc.dat$AvA; x4=mlc.dat$CD; x5=mlc.dat$SD

## original fit
orig.fit = glm(y~x1+x2+x3+x4+x5)
r.o = orig.fit$resid
y.o = orig.fit$fitted
infl = influence(orig.fit, do.coef = FALSE)
rsta = rstandard(orig.fit, infl, type = "pearson")

par(mfrow=c(1,1))
plot(y.o,r.o,main="Residual plot of original fitting")


## transformed fit
yy=sqrt(y)
tran.fit = glm(yy~x1+x2+x3+x4+x5)
r.t = tran.fit$resid
y.t = tran.fit$fitted
infl = influence(tran.fit, do.coef = FALSE)
cook.d =cooks.distance(tran.fit, infl,res = infl$pear.res, 
	        dispersion = summary(tran.fit)$dispersion,
	        hat = infl$hat)

par(mfrow=c(1,2))
plot(y.t,r.t,main="Residual plot of transformed fitting")
qqnorm(cook.d,main="Q-Q plot of cook's distance")

cbind(order(cook.d), cook.d[order(cook.d)])

## refitting after deleting outliers
wk.dat1=data.frame(yy,x1,x2,x3,x4,x5)
wk.dat2=wk.dat1[-c(7,37,45),]
tran.fit2 = glm(yy~x1+x2+x3+x4+x5,data=wk.dat2)
r.t = tran.fit2$resid
y.t = tran.fit2$fitted
rsta = rstandard(tran.fit2, infl, type = "pearson")
infl = influence(tran.fit2, do.coef = FALSE)
cook.d =cooks.distance(tran.fit2, infl,res = infl$pear.res, 
	        dispersion = summary(tran.fit2)$dispersion,
	        hat = infl$hat)
par(mfrow=c(2,2))
plot(y.t,r.t,main="Residual plot of transformed fitting")
qqnorm(cook.d,main="Q-Q plot of cook's distance")
qqnorm(rsta,main="Q-Q plot of standardized residual")


summary(tran.fit2)


vifChen=function(object) {
 ## Input:
  # object -- a fitted lm object.
 ## Output:
  # vif -- variance inflation factors

 X = object$x[,-1]
 V = vcov(object)[-1,-1]
 n = dim(X)[1] 
 sigma=summary(object)$sigma

 v = diag(V)
 S = diag(var(X))*(n-1)
 vif = v*S/sigma^2
 vif
}

tran.fit2 = lm(yy~x1+x2+x3+x4+x5,x=TRUE,data=wk.dat2)
vifChen(tran.fit2)

##################################################
## Trial of animal survival example 
 y=c(3.1, 4.5, 4.6, 4.3, 
 8.2, 11.0, 8.8, 7.2, 
 4.3, 4.5, 6.3, 7.6 ,
 4.5, 7.1, 6.6, 6.2, 
 3.6, 2.9, 4.0, 2.3 ,
 9.2, 6.1, 4.9, 12.4, 
 4.4, 3.5, 3.1, 4.0 ,
 5.6, 10.2, 7.1, 3.8, 
 2.2, 2.1, 1.8, 2.3 ,
 3.0, 3.7, 3.8, 2.9,
 2.3, 2.5, 2.4, 2.2 ,
 3.0, 3.6, 3.1, 3.3 )

x=factor(rep(c(1:12),rep(4,12))) 


## Method I for determining alpha
fit.o=lm(y~x)
r.o=fit.o$resid
y.o=fit.o$fitted

par(mfrow=c(1,1))
plot(y.o,r.o,main="Residual plot of the original fit" )

# ungrouped
alpha.fit1=lm(log(abs(r.o))~log(y.o))
summary(alpha.fit1)

# grouped
YY = matrix(y,ncol=4,byrow=T)
ybar = apply(YY,1,mean)
sd = sqrt(apply(YY,1, var))
alpha.fit2 = lm(log(sd)~log(ybar))
summary(alpha.fit2)

## method II 
alpha=c(-1, -0.5,0.5,1,1.5,2,3)
R=NULL
for (i in 1:length(alpha) ) { 
  R[i] = max(sd/ybar^alpha[i])/ min(sd/ybar^alpha[i])
}

A=as.matrix(data.frame(alpha=alpha, R=round(R,2)))
t(A)

### Transformed fit

yy=1/y
fit.inv=lm(yy~x)
r.inv=fit.inv$resid
y.inv=fit.inv$fitted
plot(y.inv,r.inv)
summary(fit.inv)


#####################################################
## peptic ulcer example

YY=matrix(c(0.2, 4.9, 17.6, 0.2, 2.5, 8.8,
0.3, 5.0, 18.9, 0.3, 2.8, 9.1,
0.4, 5.3, 20.7, 0.4, 3.6, 10.3,
1.1, 7.5, 24.0, 0.7, 4.8, 15.6,
2.0, 9.8, 25.4, 1.2, 4.8, 16.1,
2.1, 10.4, 40.0, 1.5, 5.4, 16.5,
3.3, 10.9, 42.2, 1.5, 5.7, 16.7,
3.8, 11.3, 50.0, 1.9, 5.8, 20.0,
4.5, 12.4, 60.0, 2.0, 7.5, 20.7,
4.8, 16.2, 0,2.4, 8.7, 33.0), ncol=6,byrow=T)

y1=c(YY[,1],YY[,2],YY[,3])
y1=y1[-length(y1)]
y2=c(YY[,4],YY[,5],YY[,6])

y1=c(0.2,0.3,0.4,1.1,2.0,2.1,3.3,3.8,4.5,4.8,4.9,5.0,5.3,7.5,9.8,10.4,
     10.9,11.3,12.4,16.2,17.6,18.9,20.7,24.0,25.4,40.0,42.2,50.0,60.0)
y2=c(0.2,0.3,0.4,0.7,1.2,1.5,1.5,1.9,2.0,2.4,2.5,2.8,3.6,4.8,4.8,5.4,5.7,
     5.8,7.5,8.7,8.8,9.1,10.3,15.6,16.1,16.5,16.7,20.0,20.7,33.0)

lambda=c(-1,-0.5,0.001,0.5); s1=NULL; s2=NULL
for (i in 1:4) {
  z1=(y1^lambda[i]-1)/lambda[i]; z2=(y2^lambda[i]-1)/lambda[i]
  s1[i]=sqrt(var(z1)); s2[i]=sqrt(var(z2))
}  
s.ratio = s1/s2

A = as.matrix(data.frame(lambda=lambda,s.ratio =s.ratio))
t(A)

y=c(y1,y2)
y.t=1/y
x=factor(c(rep(1,29),rep(2,30)))

fit.o=lm(y~x); fit.t =lm(y.t~x)

summary(fit.o); summary(fit.t)

##########################
## Worsted yarn experiment example

 Yarn=read.csv(file="D://Rsession/WorstedYarn.csv", header = TRUE)

 y=Yarn[,4]
 x1=factor(Yarn[,1])
 x2=factor(Yarn[,2])
 x3=factor(Yarn[,3])
 yarn.dat = data.frame(y,x1,x2,x3)

## Linear model with y as response
 lfit = lm(y~x1+x2+x3,data=yarn.dat);  summary(lfit)
 y.o=lfit$fitted;  r.o=lfit$resid
 par(mfrow=c(1,1))
   plot(y.o, r.o)


  y.shift = y.o-min(y.o)+1
  summary( lm(abs(r.o)~y.o )

## Linear model with log(y) as response
 lnfit=lm(log(y)~x1+x2+x3,data=yarn.dat)
 summary(lnfit)
 BIC(lnfit)
 lnfitted=lnfit$fitted
 lnres=lnfit$resi
 par(mfrow=c(1,2))
 plot(lnres)
 plot(lnfitted,lnres)

