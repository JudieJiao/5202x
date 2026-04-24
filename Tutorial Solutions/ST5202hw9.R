###########################
## ST5202 Tutorial 9
###########################


#######################################################
## Question 1: riboflavin data: log-transformed riboflavin production rate as response

library(hdi)
library(glmnet)
library(ncvreg)

# -----------------------------
# 1. Load data
# -----------------------------
data(riboflavin)

X <- unclass(riboflavin$x)
y <- riboflavin$y

cat("Data dimension of X:", dim(X)[1], "x", dim(X)[2], "\n")
cat("Length of y:", length(y), "\n\n")

# -----------------------------
# Fixed train/test split
# -----------------------------
set.seed(123)
n <- length(y)
train_id <- sample(seq_len(n), round(0.7 * n))
test_id <- setdiff(seq_len(n), train_id)

X_train <- as.matrix(X[train_id, , drop = FALSE])
X_test  <- as.matrix(X[test_id, , drop = FALSE])
y_train <- y[train_id]
y_test  <- y[test_id]

cat("Training size:", nrow(X_train), "\n")
cat("Test size:", nrow(X_test), "\n\n")

# -----------------------------
# Fixed 10-fold CV assignment
# -----------------------------
set.seed(123)
foldid <- sample(rep(1:10, length.out = length(y_train)))

# -----------------------------
# Lasso
# -----------------------------
cv_lasso <- cv.glmnet(
  X_train, y_train,
  foldid = foldid
)

pdf("lasso_cv_plot.pdf", width = 6, height = 5)
plot(cv_lasso)
dev.off()

lasso_lambda <- cv_lasso$lambda.min
lasso_lambda
lasso_fit <- glmnet(X_train, y_train, lambda = lasso_lambda)
lasso_pred <- as.numeric(predict(lasso_fit, newx = X_test))
lasso_mse <- mean((y_test - lasso_pred)^2)
lasso_mse
summary(lasso_fit)
lasso_coef <- as.matrix(coef(lasso_fit))
lasso_nz <- sum(lasso_coef[-1, 1] != 0)
lasso_nz

# -----------------------------
# MCP
# -----------------------------

cv_mcp <- cv.ncvreg(
  X_train, y_train,
  penalty = "MCP",
  fold = foldid
)

pdf("mcp_cv_plot.pdf", width = 6, height = 5)
plot(cv_mcp)
dev.off()

mcp_lambda <- cv_mcp$lambda.min
mcp_lambda
mcp_pred <- as.numeric(predict(cv_mcp$fit, X_test, lambda = mcp_lambda))
mcp_mse <- mean((y_test - mcp_pred)^2)
mcp_mse
mcp_coef <- coef(cv_mcp$fit, lambda = mcp_lambda)
mcp_nz <- sum(mcp_coef[-1] != 0)
mcp_nz
# -----------------------------
# SCAD
# -----------------------------
cv_scad <- cv.ncvreg(
  X_train, y_train,
  penalty = "SCAD",
  fold = foldid
)

pdf("scad_cv_plot.pdf", width = 6, height = 5)
plot(cv_scad)
dev.off()

scad_lambda <- cv_scad$lambda.min
scad_lambda
scad_pred <- as.numeric(predict(cv_scad$fit, X_test, lambda = scad_lambda))
scad_mse <- mean((y_test - scad_pred)^2)
scad_mse
scad_coef <- coef(cv_scad$fit, lambda = scad_lambda)
scad_nz <- sum(scad_coef[-1] != 0)
scad_nz

# -----------------------------
# Results table
# -----------------------------
result <- data.frame(
  Method = c("Lasso", "MCP", "SCAD"),
  Lambda = c(lasso_lambda, mcp_lambda, scad_lambda),
  Nonzero = c(lasso_nz, mcp_nz, scad_nz),
  Test_MSE = c(lasso_mse, mcp_mse, scad_mse)
)

print(result)

best_mse_method <- result$Method[which.min(result$Test_MSE)]
most_sparse_method <- result$Method[which.min(result$Nonzero)]

cat("\nMethod with smallest test MSE:", best_mse_method, "\n")
cat("Most sparse method:", most_sparse_method, "\n")

#######################################################
## Question 2: SMSA data: number of physicians as response

  SMSA.dat = read.table("D://Teaching/22-23IST3131/Rsession/SMSA.txt",header=TRUE)

  y=SMSA.dat$doctor; x1=SMSA.dat$Land; x2=SMSA.dat$T.p; x3=SMSA.dat$P.city; x4=SMSA.dat$p.65; 
  x5=SMSA.dat$beds; x6=SMSA.dat$h.sch;x7=SMSA.dat$labor; x8=SMSA.dat$income; x9=SMSA.dat$crimes;
  x10=factor(SMSA.dat$region);
  p.dat=data.frame(y,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)


## (i) selecting models
  library(MASS)

# Starting models for the selection procedure:
  null1 = lm(y~1, data=p.dat); null2 = lm(y~x10, data=p.dat)

# Selection with null model as starting model
  f1=stepAIC(null1,scope=list(upper=~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,lower=~1),direction="forward")
  summary(f1)
 
# Selection with starting model containg x10
  f2=stepAIC(null2,scope=list(upper=~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,
                lower=~x10),direction="forward")
  summary(f2)

  AIC(f1); AIC(f2)
## f1, which has the smaller AIC, is the selected model. The select model: y.p~x1 + x2 + x5 + x8 + x9 + x10.


## (ii) Diagnostics
# compute raw materials
  fit1=glm(y~x1 + x2 + x5 + x8 + x9 + x10,data=p.dat)
  yhat = fit1$fitted.values
  r =	residuals(fit1,type="pearson")
  h =	hatvalues(fit1,type="diagonal")
	infl = influence(fit1, do.coef = FALSE)
  rsta = rstandard(fit1, infl, type = "pearson")
  rstu = rstudent(fit1, infl, type = "pearson")
  cook =cooks.distance(fit1, infl,res = infl$pear.res, 
	        dispersion = summary(fit1)$dispersion,hat = infl$hat)
#(iia) Check for linearity
  par(mfrow=c(2,3))
  plot(x1, r,xlab="Residual vs. x1")
  plot(x2, r,xlab="Residual vs. x2")
  plot(x5, r,xlab="Residual vs. x5")
  plot(x8, r,xlab="Residual vs. x8")
  plot(x9, r,xlab="Residual vs. x9")
  plot(yhat,r,xlab="Residual vs. fitted value")


# (iib) Check for homogeneity and normality
  par(mfrow=c(1,2))
  plot(yhat, r,xlab="Residual vs. fitted values")
  qqnorm(rsta,xlab="Q-Q plot of standardized residuals")

## (iii) Check for outliers
# (iiia)
  par(mfrow=c(2,2))
  qqnorm(h,xlab="",ylab="",main="Q-Q  plot of the hat values") # checking leverage
  qqnorm(rstu,xlab="",ylab="",main="Q-Q  plot of the studentized deletion residuals") # checking consistency
  qqnorm(cook,xlab="",ylab="",main="Q-Q  plot of the Cook's distances") # checking influence

 # (iiib)
  p=9 
  n=length(h)
  round(h[h>=2*p/n],3)

# (iiic) 
  d = sqrt(var(rstu))
  round(rstu[abs(rstu)>2*d],3)

# (iiid) 
  round(cook[cook>0.1],3)

## (iv) 
  u1=u2=u3=u4=u5=u6=u7=u8=rep(0,n)
  u1[1]=1; u2[3]=1; u3[4]=1;u4[5]=1; u5[7]=1; u6[8]=1; u7[9]=1; u8[27]=1

  fit2=glm(y~x1 + x2 + x5 + x8 + x9 + x10 +u1+u2+u3+u4+u5+u6+u7+u8)
  summary(fit2)

  qt(0.05/(2*8),124,lower.tail=FALSE)

 

#############################################################################
## Question 3. Kidney 


kidney.dat = read.table("D://Teaching/22-23IST3131/Rsession/kidney.txt",header=TRUE)

pairs(kidney.dat)

attach(kidney.dat)
 y=C.clear; x1=serum; x2=Age.y; x3=weight


 fit1=glm(C.clear~serum+Age.y+weight,data=kidney.dat)
  yhat = fit1$fitted.values
  r =	residuals(fit1,type="pearson")
  h =	hatvalues(fit1,type="diagonal")
	infl = influence(fit1, do.coef = FALSE)
  rsta = rstandard(fit1, infl, type = "pearson")
  rstu = rstudent(fit1, infl, type = "pearson")
  cook =cooks.distance(fit1, infl,res = infl$pear.res, 
	        dispersion = summary(fit1)$dispersion,hat = infl$hat)



# (iia) Check for homogeneity and normality
  par(mfrow=c(1,2))
  plot(yhat, r,xlab="Residual vs. fitted values")
  qqnorm(rsta,xlab="Q-Q plot of standardized residuals")

# (iib)
  par(mfrow=c(2,2))
  qqnorm(h,xlab="",ylab="",main="Q-Q  plot of the hat values") # checking leverage
  qqnorm(rstu,xlab="",ylab="",main="Q-Q  plot of the studentized deletion residuals") # checking consistency
  qqnorm(cook,xlab="",ylab="",main="Q-Q  plot of the Cook's distances") # checking influence


##################
n = length(r)
rstu[order(rstu)][c(1,2,3,n)] 
cook[order(cook)][(n-3):n]

# (iic)
  n=length(r)
  u1=rep(0,n);u2=u1;u3=u1;u4=u5=u1;
  u1[16]=1;u2[20]=1;u3[21]=1;u4[26]=1; u5[29]=1
  outlier=glm(C.clear~serum+Age.y+weight+u1+u2+u3+u4+u5,data=kidney.dat)
  summary(outlier)
		

