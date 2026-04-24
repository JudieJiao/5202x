### Leccture 6 ########
### Insurance Example
insur.dat=read.table("insurance.txt",header=TRUE)

#Selection criteria
library(boot)
full.fit = glm(log.Y.~X1+X2+X3+X4+X5+X6+X7+X8,data=insur.dat)

set.seed(123345)
cv.glm(insur.dat,full.fit)$delta
CV.1.f = cv.glm(insur.dat,full.fit)$delta[2]
CV.k.f = cv.glm(insur.dat,full.fit,K=11)$delta[2]
aic.f = AIC(full.fit)
bic.f = BIC(full.fit)

v1=c(CV.1.f,CV.k.f,aic.f,bic.f)

reduce.fit = glm(log.Y.~X1+X2+X3+X8,data=insur.dat)
CV.1.r = cv.glm(insur.dat,reduce.fit)$delta[2]
CV.k.r = cv.glm(insur.dat,reduce.fit,K=11)$delta[2]
aic.r = AIC(reduce.fit)
bic.r = BIC(reduce.fit)

v2=c(CV.1.r,CV.k.r,aic.r,bic.r)

cbind(v1,v2)



#Seclection procecure

# selectin by removing redundant variables.

m1=lm(formula = log.Y. ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8, 
      data = insur.dat)
summary(m1); anova(m1)


# pi= 4*(4*atan(1/5) - atan(1/239)); AIC=54 *log (2.0052/54) + 2*8 + 54 *(log(2*pi) + 1) 
# tmp= lm(formula = log.Y. ~ X1 + X2 + X3  + X5 + X6 + X8, data = insur.dat); AIC(tmp);
# sigma=summary(tmp)$sigma;  AIC=54 *log (sigma^2*47/54) + 2*8 + 54 *(log(2*pi) + 1) 




m2=lm(formula = log.Y. ~ X1 + X2 + X3 + X5 + X6 + X7 + X8,
      data = insur.dat)
summary(m2)

m3=lm(formula = log.Y. ~ X1 + X2 + X3 + X5 + X6 + X8, 
      data = insur.dat)
summary(m3)

m4=lm(formula = log.Y. ~ X1 + X2 + X3 + X6 + X8, data = insur.dat)
summary(m4)

m5=lm(formula = log.Y. ~ X1 + X2 + X3 + X8, data = insur.dat)
summary(m5)

# Sequential procedures.

library(MASS)


# Forward selection
null.m=glm(log.Y.~1,data=insur.dat)
forward.fit=stepAIC(null.m,scope = list(upper = ~X1+X2+X3+X4+X5+X6+X7+X8, lower = ~1),direction = "forward")
summary(forward.fit)

# Backward selection
full.m=lm(log.Y.~X1+X2+X3+X4+X5+X6+X7+X8,data=insur.dat)
backward.fit=stepAIC(full.m,scope = list(upper = ~X1+X2+X3+X4+X5+X6+X7+X8, lower = ~1),direction = "backward")
summary(backward.fit)


# Stepwise selection
null.m=lm(log.Y.~1,data=insur.dat)
stepwise.fit=stepAIC(null.m,scope = list(upper = ~X1+X2+X3+X4+X5+X6+X7+X8, lower = ~1),direction = "both")
summary(stepwise.fit)


full.m=lm(log.Y.~X1+X2+X3+X4+X5+X6,data=insur.dat)

StepDown=stepAIC(full.m,scope = list(upper = ~X1+X2+X3+X4+X5+X6, lower = ~1),direction = "both")
summary(StepDown)

stepwise.fit=stepAIC(null.m,scope = list(upper = ~X1+X2+X3+X4+X5+X6, lower = ~1),direction = "both")


## Pull strength example

pull.dat=read.csv(file="pullstrength.csv",
                  header=TRUE)
fit.1 = lm(y~1, data=pull.dat)
fit.2 = lm(y~x1+x2+x3+x4+x5+x6, data=pull.dat)
forward=stepAIC(fit.1,scope=list(upper=~x1+x2+x3+x4+x5+x6,
                                 lower=~1),direction="forward")
backward=stepAIC(fit.2,scope=list(upper=~x1+x2+x3+x4+x5+x6,
                                  lower=~1),direction="backward")
stepwiseF=stepAIC(fit.1,scope=list(upper=~x1+x2+x3+x4+x5+x6,
                                   lower=~1),direction="both")
stepwiseB=stepAIC(fit.2,scope=list(upper=~x1+x2+x3+x4+x5+x6,
                                   lower=~1),direction="both")


# Penalized likelihood
library(iterators)
library(foreach)
library(Matrix)
library(shape)
library(glmnet)

# insurance example
x = as.matrix(insur.dat[,1:8])
y = as.vector(insur.dat[,10])
fit1 = cv.glmnet(x, y)

coef(fit1, s = fit1$lambda.min)

# Simulated example
set.seed(1010)
n = 1000
p = 100
nzc = trunc(p/10)
x = matrix(rnorm(n * p), n, p); dim(x)
beta = rnorm(nzc)
fx = x[, seq(nzc)] %*% beta
eps = rnorm(n) * 5
y = drop(fx + eps)
cvfit = cv.glmnet(x, y)
beta_lasso = coef(cvfit, s = cvfit$lambda.min)
B <- as.matrix(coef(cvfit, s = cvfit$lambda.min))
beta_init <- B[-1, 1]
names(beta_init) <- rownames(B)[-1]
beta_init[beta_init != 0]

# Adaptive LASSO
# adaptive weights
gamma <- 1
eps <- 1e-6
w <- 1 / (abs(beta_init)^gamma + eps)

cv_adalasso <- cv.glmnet(
  x, y,
  alpha = 1,                 # lasso
  penalty.factor = w,        # adaptive weights
  standardize = TRUE
)

coef(cv_adalasso, s = "lambda.min")
plot(cv_adalasso)

plot(beta_lasso, coef(cv_adalasso, s = "lambda.min"))
abline(a = 0, b = 1, col = 2)
library(ncvreg)

# MCP
fit_mcp <- ncvreg(x, y, penalty = "MCP")
plot(fit_mcp)

cv_mcp <- cv.ncvreg(x, y, penalty = "MCP")
plot(cv_mcp)


coef(cv_mcp, lambda = cv_mcp$lambda.min)
beta.mcp = as.vector(coef(cv_mcp, lambda = cv_mcp$lambda.min))[-1]
coef(cv_mcp, lambda = cv_mcp$lambda.min)[which(beta.mcp != 0)+1]

# SCAD
fit_scad <- ncvreg(x, y, penalty = "SCAD")
plot(fit_scad)

cv_scad <- cv.ncvreg(x, y, penalty = "SCAD")
plot(cv_scad)

coef(cv_scad, lambda = cv_scad$lambda.min)
beta.scad = as.vector(coef(cv_scad, lambda = cv_scad$lambda.min))[-1]
coef(cv_scad, lambda = cv_scad$lambda.min)[which(beta.scad != 0)+1]
beta
