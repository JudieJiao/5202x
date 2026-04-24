### Homework 10 ########
### Remedy: non-constant variance and multicollinearity

##############################################################################
## 1. Employee salary example

salary.dat = read.table("Rsession/Salary.txt",header =TRUE)

y=salary.dat$salary
degree =salary.dat$degree
D=factor(degree)
x1=salary.dat$YOE
x2=salary.dat$NSP

wk.dat=data.frame(y,D,x1,x2)

## unweighted fitting
uw.fit=glm(y~D+x1+x2,data=wk.dat)
r = uw.fit$resid
yhat=uw.fit$fitted

## raw materials for outlier detection
hat = hatvalues(uw.fit,type="diagonal")
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

## remove influential observations and refit
wk.dat3 = wk.dat[-c(3,43),]

uw.fit=glm(y~D+x1+x2,data=wk.dat3)
r = uw.fit$resid
yhat=uw.fit$fitted
par(mfrow=c(1,2))
plot(yhat,r, xlab = "Fitted Value", ylab = "Residuals")
plot(yhat, abs(r), xlab = "Fitted Value", ylab = "Absolute Residuals")

############################################################
## (iii) Estimate weights by log-log regression
## log|r_i| on log(yhat_i)
############################################################

weight.loglog = lm(log(abs(r)) ~ log(yhat))
s.loglog = exp(weight.loglog$fitted)
w.loglog = 1 / s.loglog^2

w.fit.loglog = lm(y ~ D + x1 + x2, weight = w.loglog, data = wk.dat3)

r.loglog = w.fit.loglog$resid
y.loglog = w.fit.loglog$fitted

r.loglog.w = r.loglog * sqrt(w.loglog)
y.loglog.w = y.loglog * sqrt(w.loglog)

par(mfrow = c(1, 2))
plot(y.loglog.w, r.loglog.w,
     xlab = "Weighted fitted values",
     ylab = "Weighted residuals",
     main = "Log-log weights")
plot(y.loglog.w, abs(r.loglog.w),
     xlab = "Weighted fitted values",
     ylab = "|Weighted residuals|",
     main = "Log-log weights")

summary(w.fit.loglog)

############################################################
## (iv) Estimate weights by linear regression
## |r_i| on yhat_i
############################################################

weight.lin = lm(abs(r) ~ yhat)
s.lin = weight.lin$fitted

w.lin = 1 / s.lin^2

w.fit.lin = lm(y ~ D + x1 + x2, weight = w.lin, data = wk.dat3)

r.lin = w.fit.lin$resid
y.lin = w.fit.lin$fitted

r.lin.w = r.lin * sqrt(w.lin)
y.lin.w = y.lin * sqrt(w.lin)

par(mfrow = c(1, 2))
plot(y.lin.w, r.lin.w,
     xlab = "Weighted fitted values",
     ylab = "Weighted residuals",
     main = "Linear weights")
plot(y.lin.w, abs(r.lin.w),
     xlab = "Weighted fitted values",
     ylab = "|Weighted residuals|",
     main = "Linear weights")

summary(w.fit.lin)

############################################################
## (v) Brief comparison
############################################################

se.uw = sqrt(diag(vcov(uw.fit)))
se.loglog = sqrt(diag(vcov(w.fit.loglog)))
se.lin = sqrt(diag(vcov(w.fit.lin)))

data.frame(
  coefficient = names(coef(uw.fit)),
  unweighted = se.uw,
  loglog.weight = se.loglog,
  linear.weight = se.lin
)

############################################################
## Q3. longley data: multicollinearity
############################################################

data(longley)

############################################################
## Full model
############################################################

fit.full = lm(Employed ~ GNP.deflator + GNP + Unemployed +
                Armed.Forces + Population + Year,
              x = TRUE, data = longley)

summary(fit.full)

############################################################
## Variance inflation factor
############################################################

vifChen = function(object) {
  ## Input:
  ##   object -- a fitted lm object. The argument x of lm
  ##             must be set as x = TRUE.
  ## Output:
  ##   vif -- variance inflation factors
  
  X = object$x[, -1]
  V = vcov(object)[-1, -1]
  n = dim(X)[1]
  sigma = summary(object)$sigma
  
  v = diag(V)
  S = diag(var(X)) * (n - 1)
  vif = v * S / sigma^2
  vif
}

vifChen(fit.full)

############################################################
## Remove the predictor with the largest VIF
############################################################

fit.red = lm(Employed ~ GNP.deflator + Unemployed +
               Armed.Forces + Population + Year,
             x = TRUE, data = longley)

summary(fit.red)
vifChen(fit.red)

############################################################
## Remove the predictor with the largest VIF
############################################################

fit.red = lm(Employed ~ GNP.deflator + Unemployed +
               Armed.Forces + Population,
             x = TRUE, data = longley)

summary(fit.red)
vifChen(fit.red)

##############################################################################
## 4. mtcars example: VIF, auxiliary regression, ridge and PCR

## full linear regression model
regLSE = lm(mpg~cyl +disp + hp +drat + wt + qsec + vs + am + gear + carb,
            x=TRUE, data=mtcars)
summary(regLSE)
vifChen(regLSE)

## ridge regression
library(MASS)
regRidge = lm.ridge(mpg~cyl +disp + hp +drat + wt + qsec + vs + am + gear + carb,
                    data=mtcars, lambda = seq(0, 50, 0.01))
plot(regRidge$lambda, regRidge$GCV, type="l")

lambda.best = regRidge$lambda[which(regRidge$GCV==min(regRidge$GCV))]
regRidge.best = lm.ridge(mpg~cyl +disp + hp +drat + wt + qsec + vs + am + gear + carb,
                         data=mtcars, lambda = lambda.best)
regRidge.best
regRidge.best$coef
regRidge.best$lambda
regRidge.best$GCV

## Standardize X in the same way as in regRidge.best
X.new = as.matrix(mtcars[, c("cyl","disp","hp","drat","wt","qsec","vs","am","gear","carb")])
X.new.std = scale(X.new,
                  center = regRidge.best$xm,
                  scale  = regRidge.best$scales)

## Fitted values on original response scale
ridge.pred = as.vector(regRidge.best$ym + X.new.std %*% regRidge.best$coef)
mean((mtcars$mpg - ridge.pred)^2)

## principal component regression
library(pls)
pcr.fit = pcr(mpg~cyl +disp + hp +drat + wt + qsec + vs + am + gear + carb,
              data=mtcars, scale=TRUE, validation="CV")
summary(pcr.fit)
validationplot(pcr.fit, val.type="MSEP")

rmsep.pcr = RMSEP(pcr.fit)$val[1,1,-1]
ncomp.best = which.min(rmsep.pcr)
ncomp.best

pcr.best = pcr(mpg~cyl +disp + hp +drat + wt + qsec + vs + am + gear + carb,
               data=mtcars, scale=TRUE, ncomp=ncomp.best)
summary(pcr.best)

## fitted values from the selected PCR model
pcr.pred = as.numeric(predict(pcr.best, ncomp=ncomp.best))
mean((mtcars$mpg - pcr.pred)^2)

