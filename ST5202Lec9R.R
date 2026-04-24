
##############################################
# R for Lecture 9: Generalized Linear Model #
##############################################

### Challenger Data ###
#######################

### Create the data object within R
  ct=matrix(c(
     53,	1,	70,	1,
     56,	1,	70,	1,
     57,	1,	72,	0,
     63,	0,	73,	0,
     66,	0,	75,	0,
     67,	0,	75,	1,
     67,	0,	76,	0,
     67,	0,	76,	0,
     68,	0,	78,	0,
     69,	0,	79,	0,
     70,	0,	80,	0,
     70,	1,	81,	0), ncol=2,byrow=T)

Challenger=data.frame(temp = ct[,1],oring.f=ct[,2])

### Fit the ligistic model to data
  EX1.fit=glm(oring.f~temp,family=binomial(logit),data=Challenger)
  summary(EX1.fit)

### Construct confidence intervals
  b=EX1.fit$coef   # extract the estimated coefficients
  V=vcov(EX1.fit)  # extract the covariance matrix of estimated coefficients
  sd = sqrt(diag(V))

  alpha = 0.05
  q.alpha = qnorm((1-alpha/2))

  b.ci=cbind(b-q.alpha*sd, b+q.alpha*sd) # C.I. for beta 

  oddsR=exp(b.ci)             # C.I. for the odds ratio

  oddsR
###################################################################################



### Aneurysm Data ###
#####################
  Aneurysm.dat=read.csv(file="Rsession/Aneurysm.csv", header = TRUE)  # Read the data
  
  glm(comp30~control+male+early+diameter,family=binomial(logit),data=Aneurysm.dat)

  tt=na.omit(Aneurysm.dat)      # Eliminate NA records from the data frame
  EX2.fit=glm(comp30~control+male+early+diameter,family=binomial(logit),data=tt)  

  summary(EX2.fit)
  b=EX2.fit$coef[-1]               # estimates of coefficients (excluding the intercept)
  V=vcov(EX2.fit)[-1,-1]           # covariance matrix of estimates (excluding intercept)
  sd = sqrt(diag(V))
  resi = resid(EX2.fit,type="pearson")   #  Extract the Pearson residuals

  attach(tt)     # make the variables in tt accessable in the current session
  n=length(comp30)
  p = length(b)+1
  disp=sum(resi^2)/(n-p)

### Construct 100(1-alpha)% C.I. for the coefficients (excluding the intercept)
  alpha=0.05
  q=qnorm(1-alpha/2)
  b.ci = cbind(b - q*sd, b+q*sd)     # C.I for coefficients

  exp(b.ci)                          # C.I for odds ratio
##########################################################################################

### Mine Fracture data ###
###########################

  mine.dat=read.csv(file="Rsession/MineFracture.csv", header = TRUE)  # Read the data

### Stepwise selection
  library(MASS)
  null.fit=glm(y~1,family=poisson(log),data=mine.dat)
  full.fit=glm(y~.,family=poisson(log),data=mine.dat)

  #Forward selection
  forward=stepAIC(null.fit, scope=list(upper=~x1+x2+x3+x4, lower=~1), direction="forward")
  summary(forward)
  
  #Backward selection
  backward=stepAIC(full.fit, scope=list(upper=~x1+x2+x3+x4, lower=~1),	direction="backward")
  summary(backward)
  
## Fitting the selected model
  mine.fit=glm(y~x1+x2+x4,family=poisson(log),data=mine.dat)
  summary(mine.fit)
  
##############################################################################################