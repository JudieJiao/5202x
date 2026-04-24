#########################
## ST5202  Lecture notes 7
########################


### Pull strength example
pull.dat=read.csv(file="Rsession/pullstrength.csv",
                  header=TRUE)
glm.fit = glm(y~x1+x2+x3+x4+x5+x6, data=pull.dat)
summary(glm.fit)
fitted = glm.fit$fitted.values
res =	residuals(glm.fit,type="pearson")
hat =	hatvalues(glm.fit,type="diagonal")
infl = influence(glm.fit, do.coef = FALSE)
rsta = rstandard(glm.fit, infl, type = "pearson")
rstu = rstudent(glm.fit, infl, type = "pearson")
cook.d =cooks.distance(glm.fit, infl,res = infl$pear.res, 
                       dispersion = summary(glm.fit)$dispersion,
                       hat = infl$hat)

round(cbind(fitted,res,hat,rsta,rstu,cook.d),4)


## Figure of null pattern
x=rnorm(100,0,1)
ind =1:100

plot(ind,x, xlab="",ylab="")

##Checking for nonlinearity

rider.dat = read.table("Rsession/CH03TA01.txt")
attach(rider.dat)
rider=rider.dat[,1]
m.dist = rider.dat[,2]

fit = lm(rider~m.dist)

fitted = fit$fitted
res = residuals(fit)

par(mfrow=c(2,2), mar = c(3, 3, 0.4, 0.4),
  oma = c(0, 0, 0, 0),
  mgp = c(1.1, 0.3, 0),
  tcl = -0.2)
plot(fitted, res,xlab="Fitted values", ylab="Residuals")
plot(m.dist, res,xlab="Predictor", ylab="Residuals")
plot(m.dist, rider, xlab="Predictor", ylab="Ridership increased")
abline(fit)

##Checking for homogeneity

set.seed(13456)

x=runif(100,0,1)*10
mu=1+2*x
sigma=c(1:100)/10
e =rnorm(100,0,1)*sigma
y=mu+e

ncv=lm(y~x)

fitted = ncv$fitted
res = residuals(ncv)

par(mfrow=c(1,2))
plot(x, res,xlab="Predictor", ylab="Residuals")

plot(fitted, res,xlab="Fitted values", ylab="Residuals")



##Checking for independence

model=list(order=c(1,0,0),ar=c(0.8))
z=arima.sim(model=model,n=100)

y=mu+z


nind=lm(y~x)

fitted = nind$fitted
res = residuals(nind)

index=1:100

par(mfrow=c(1,1))
plot(index, res, xlab="Time", ylab="Residuals")


##Checking for normality 
y=mu+e

nck = glm(y~x)
fitted = nck$fitted.values
res =	residuals(nck,type="pearson")
infl = influence(nck, do.coef = FALSE)
rsta = rstandard(nck, infl, type = "pearson")
par(mfrow=c(2,2))
qqnorm(rnorm(100,0,1), xlab="", ylab="",main="The NULL pattern")
qqnorm(res,xlab="", ylab="",main="Q-Q plot of residuals")
qqnorm(rsta,xlab="", ylab="",main="Q-Q plot of studentized residuals")

z=rnorm(100,0,1)
t=rt(100,1)
g2=rgamma(100,2)
c = rcauchy(100)

par(mfrow=c(2,2))

qqnorm(z,xlab="",ylab="",main="Null pattern: normal")
qqnorm(t,xlab="",ylab="",main="Symetric but heavy tailed")
qqnorm(g2,xlab="",ylab="",main="Skewed to the right")


## Missing predictor

pull.dat=read.csv(file="Rsessions/pullstrength.csv",
                  header=TRUE)
u.fit = glm(y~x4, data=pull.dat)
fitted = u.fit$fitted.values
res =	residuals(u.fit,type="pearson")


x1=pull.dat$x1
x2=pull.dat$x2
x3=pull.dat$x3
x5=pull.dat$x5
x6=pull.dat$x6


par(mfrow=c(2,3))
plot(fitted, res,  xlab="Fitted values",ylab="Residuals") 
plot(x1,res, xlab="x1",ylab="Residuals")
plot(x2,res, xlab="x2",ylab="Residuals")
plot(x3,res, xlab="x3",ylab="Residuals")
plot(x5,res, xlab="x5",ylab="Residuals")
plot(x6,res, xlab="x6",ylab="Residuals")




### Surgical unit example

library(MASS)

surgical.dat=read.table("Rsession/surgical.txt",header=TRUE)

attach(surgical.dat)

par(mfrow=c(2,2))
plot(clotting, survival)
plot(prognostic, survival)
plot(enzyme, survival)
plot(liver, survival)



## scatter plots
par(mfrow=c(2,2))
plot(clotting, log.survival)
plot(prognostic, log.survival)
plot(enzyme, log.survival)
plot(liver, log.survival)

full.fit=lm(log.survival~clotting+prognostic+enzyme+liver,data=surgical.dat)
null.fit= lm(log.survival~1,data=surgical.dat)

forward=stepAIC(null.fit,scope=list(upper=~clotting+prognostic+enzyme+liver,
                                    lower=~1),direction="forward")


backward=stepAIC(full.fit,scope=list(upper=~clotting+prognostic+enzyme+liver,
                                     lower=~1),direction="backward")

summary(forward)
summary(backward)

model1 = lm(log.survival~clotting+prognostic+enzyme,data=surgical.dat)
model2 = lm(log.survival~clotting+prognostic+enzyme+liver,data=surgical.dat)

#### diagnostics for model 1

fitted1=model1$fitted
res1 = residuals(model1)

## Residual  plots for model 1
par(mfrow=c(2,2))


plot(clotting, res1)
plot(prognostic, res1)
plot(enzyme, res1)
plot(liver, res1)

##
par(mfrow=c(1,2))
plot(fitted1, res1, xlab="",ylab="",main="Model 1: Residual vs.fitted values" )
qqnorm(res1,xlab="",ylab="",main="Model 1: Q-Q plot of residual")


#### diagnostics for model 2
fitted2=model2$fitted
res = residuals(model2)
res2 = res

## Residual  plots for model 2
par(mfrow=c(2,2))
plot(clotting, res)
plot(prognostic, res)
plot(enzyme, res)
plot(liver, res)



##
model2 = glm(log.survival~clotting+prognostic+enzyme+liver,data=surgical.dat)

hat =	hatvalues(model1,type="diagonal")
infl = influence(model1, do.coef = FALSE)
rstu = rstudent(model1, infl, type = "pearson")
cook.d =cooks.distance(model1, infl,res = infl$pear.res, 
                       dispersion = summary(model1)$dispersion,
                       hat = infl$hat)

par(mfrow=c(1,2))
plot(fitted2, res2, xlab="",ylab="",main="Pearson residual vs.fitted values" )
qqnorm(rstu,xlab="",ylab="",main="Q-Q plot of studentized residual")

#### Checking on outliers


## Illustration of the rationale of manhananobis distance
A=matrix(c(1,0.5,0.5,1),ncol=2,byrow =T)
set.seed(53212)
X=matrix(rnorm(1000,0,1), ncol=2,byrow=T)

Y=X%*%A 

var(Y)
z1=Y[,1]-mean(Y[,1]) 
z2=Y[,2]-mean(Y[,2]) 



plot(Y[,1],Y[,2])

########################################

model = glm(log.survival~clotting+prognostic+enzyme+liver,data=surgical.dat)

hat =	hatvalues(model,type="diagonal")
infl = influence(model, do.coef = FALSE)
rstu = rstudent(model, infl, type = "pearson")
cook.d =cooks.distance(model, infl,res = infl$pear.res, 
                       dispersion = summary(model)$dispersion,
                       hat = infl$hat)

par(mfrow=c(2,2))
qqnorm(hat,xlab="",ylab="",main="Q-Q  plot of the leverage values")
qqnorm(rstu,xlab="",ylab="",main="Q-Q  plot of the studentized deletion residuals")
qqnorm(cook.d,xlab="",ylab="",main="Q-Q  plot of the Cook's distances")

n=length(rstu)
round(hat[hat>8/54],4)
rstu[order(rstu)][c(1,2,n-1,n-2)]
cook.d[order(cook.d)][(n-3):n]

length(cook.d)

u1=rep(0,n);u2=u3=u4=u5=u6=u1
u1[37]=1; u2[38]=1;u3[9]=1;u4[22]=1;u5[27]=1;u6[28]=1


model.dum = glm(log.survival~clotting+prognostic+enzyme+liver+u1+u2+u3+u4+u5+u6,data=surgical.dat)
summary(model.dum)



