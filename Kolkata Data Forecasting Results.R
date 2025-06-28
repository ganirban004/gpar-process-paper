##### IMPORT NECESSARY LIBRARIES #####
library(schoolmath)
library(MASS)
library(NlcOptim)
library(matrixcalc)
library(tidyverse)
library(randomForest)
library(e1071)





##### IMPORT DATASET #####
d=read.csv(file.choose(), header = T)





##### MTD(1), PAR(1) AND GPAR(1) MODELS #####
# Mode of a conditional distribution

mode.pred = function(pred.dist){
  a=c(0:(length(pred.dist)-1))
  return(min(a[which(pred.dist==max(pred.dist))]))
} 

# PTP(h) measures

ptp = function(test.sample, pred.sample){
  c = test.sample-pred.sample
  return(length(c[c==0])/length(c)*100)
}

y=d$Category
y
### the generalized kernel

tr=max(y)
k=tr+1

# the generalized kernels

g=function(a,b,c) #GPAR
{
  return(((c^(abs(a-b)))*(1-c))/(1+c-c^(b+1)-c^(tr-b+1)))
}

g1=function(a,b) #PAR
{
  if(a==b)
    return(1)
  else 
    return(0)
}

n=length(y) # sample size
  
  # general kernel matrix
  
  ker_mat=function(th)
  {
    ss1=max(y);ss2=max(y)+1
    M=matrix(0,ss2,ss2)
    for(jj in 0:ss1)
    {
      for(ii in 0:ss1)
      {
        M[jj+1,ii+1]=g(ii,jj,th)
      }
    }
    return(M)
  }
  
  # Yt's marginal
  
  id=diag(1,nrow=max(y)+1,ncol=max(y)+1)
  
  l=length(y);l
  ss=366
  
  y.training=y[1:(l-ss)];length(y.training)
  y.test=y[(l-(ss-1)):l];length(y.test)
  
  n_s=length(y.training)
  
  kk=max(y.training)+1;kk
  
  ### preparation of likelihood functions
  
  J=matrix(0, kk,1)
  for(i in 1:kk){
    for(t in 1:n_s){
      if(y.training[t]==(i-1))
        J[i]=J[i]+1
    }
  }
  
  ### preparation for MTD forecasting
  
  # calculating Q
  
  N=matrix(0,ncol=(tr+1),nrow=(tr+1))
  Q=matrix(0,ncol=(tr+1),nrow=(tr+1))
  le=length(y.training)
  p=y.training[1:(le-1)]
  q=y.training[2:le]
  for(j in 0:tr)
  {
    for(k in 0:tr)
    {
      cou=0
      for(i in 1:(le-1))
      {
        if(p[i]==j&&q[i]==k)
          cou=cou+1
        else
          cou=cou
      }
      N[j+1,k+1]=cou
    }
  }
  
  for(j in 0:tr)
  {
    for(k in 0:tr)
    {
      if(sum(N[j+1,])==0)
        Q[j+1,k+1]=N[j+1,k+1]
      else
        Q[j+1,k+1]=N[j+1,k+1]/sum(N[j+1,])  
    }
  } 
  
  # h-step ahead forecasting distribution for mtd(1)
  
  # powA_q=function(n1)
  # {
  #   if(n1==1) return(Q)
  #   if(n1==2) return(Q%*%Q)
  #   if(n1>2) return(Q%*%powA_q(n1-1))
  # }
  
  predMTD=function(uu,h)
  {
    ma=matrix.power(Q,h)
    pp=NULL
    for(i in 0:tr)
    {
      pp[i+1]=ma[uu+1,i+1]    
    }
    return(pp)
  }
  
  # predicted error for a single data set
  
  predictedErrorMTD = function(u, training.set, test.set){ 
    n = length(training.set)
    m = length(test.set)
    
    predictedErrorMTD_h = function(u, h, test.set){         
      pred.dist = sapply(1:(m-h+1), function(j) predMTD(u[n+j-1], h))
      mode = sapply(1:(m-h+1), function(j) mode.pred(pred.dist[,j]))
      return(c(ptp(test.set[h:m], mode)))  
    }      
    
    return( t(sapply(1:3, function(h) 
      predictedErrorMTD_h(u, h, test.set))) )
  }
  
  
  # initial value for par process
  
  ph.start=abs(cor(y.training[1:length(y.training)-1],
                   y.training[2:length(y.training)]))
  
  ### preparation for PAR forecasting
  
  # h-step ahead conditional distribution p_h(j|yn)
  
  conditionalPAR = function(ynh, yn, h, phi, p_v){
    q = phi^h * g1(ynh,yn) + (1-phi^h) * p_v[ynh+1]
    return(q)   
  }
  
  # mle function
  
  mlePAR=function(yy)
  {
    n=length(yy)
    likelihood_par=function(xx)
    {
      phi_hat=xx[1]
      
      f=log(J[yy[1]+1]/n)
      
      for(i in 1:(n-1))
      {
        f=f+log(phi_hat*g1(yy[i+1],yy[i])+(1-phi_hat)*(J[yy[i+1]+1]/n))
      }  
      return(-f)
    }
    opt2=optim(fn=likelihood_par, par=c(ph.start), 
               lower=c(0.01), upper=c(0.99), method="L-BFGS-B")
    return(opt2$par)
  }
  
  # distribution of Y_(n+h) given Y_n = (p_h(0), p_h(1), ..., p_h(l))
  
  predPAR = function(yn, h, pp){
    phi_hat=pp[1]
    J_eps=J/n_s
    l=tr
    return(sapply(0:l, function(j) conditionalPAR(j, yn, h, phi_hat, J_eps)))
  }  
  
  # predicted error for a single data set
  
  predictedErrorPAR = function(u, training.set, test.set){ 
    n = length(training.set)
    m = length(test.set)
    
    predictedErrorPAR_h = function(u, h, beta, test.set){         
      pred.dist = sapply(1:(m-h+1), function(j) predPAR(u[n+j-1], h, beta))
      mode = sapply(1:(m-h+1), function(j) mode.pred(pred.dist[,j]))
      return(c(ptp(test.set[h:m], mode)))  
    }      
    
    mle = mlePAR(training.set)          
    return( t(sapply(1:3, function(h) 
      predictedErrorPAR_h(u, h, mle, test.set))) )
  }
  
  ### preparation for GPAR forecasting
  
  sum_fn=function(a,v1,v2,h,yv)
  {
    s=0
    for(i in 1:(h-1))
    {
      s=s+a^(h-i)*v1%*%cbind(matrix.power(v2,h-i)[,yv+1])
    }
    return(s)
  }
  
  # ker1=function(y,th)
  # {
  #   tr1=max(y);k1=tr1+1
  #   
  #   M1=matrix(0,k1,k1)
  #   for(j in 0:tr1)
  #   {
  #     for(i in 0:tr1)
  #     {
  #       M1[j+1,i+1]=g(i,j,th)
  #     }
  #   }
  #   return(M1)
  # }
  
  # h-step ahead conditional distribution p_h(j|yn)
  
  conditionalGPAR = function(ynh, yn, h, phi, ker, p_v){
    if(h==1){return(phi*(matrix.power(ker,1)[yn+1,ynh+1])+(1-phi)*p_v[ynh+1])}
    
    else {return(phi^h*matrix.power(ker,h)[yn+1,ynh+1] +
                   (1-phi)*sum_fn(phi,t(p_v),ker,h,ynh) + 
                   (1-phi)*p_v[ynh+1])}
  }
  
  # mle function
  
  mleGPAR=function(yy)
  {
    n=length(yy)
    
    likelihood_gpar=function(x)
    {
      phi_hat=x[1]
      theta_hat=x[2]
      
      J1=((id-phi_hat*t(ker_mat(theta_hat)))*(1/(1-phi_hat)))%*%(J/n)
      
      f=log(J[yy[1]+1]/n)
      
      for(i in 1:(n-1))
      {
        if(is.finite(log(phi_hat*g(yy[i+1],yy[i],theta_hat)+(1-phi_hat)*(J1[yy[i+1]+1]))))
        {
          f=f+log(phi_hat*g(yy[i+1],yy[i],theta_hat)+(1-phi_hat)*(J1[yy[i+1]+1]))
        }
        else{}  
      }
      return(-f)
    }
    
    # opt1=optim(fn=likelihood_gpar, par=c(phi-0.1,theta+0.1), lower=c(0.01,0.01), 
    #            upper=c(0.99,0.99), method="L-BFGS-B")
    # return(opt1$par) 
    
    confun=function(x){
      f=NULL
      f=(id-x[1]*t(ker_mat(x[2])))%*%(J/n) 
      return(list(ceq=NULL,c=-f))       # constraint fn: (I-\phi*K')p_y>=0
    }
    
    sol=solnl(c(0.5,0.5),objfun=likelihood_gpar,
              confun=confun,lb=c(0.01,0.01),ub=c(0.99,0.99))
    
    return(sol$par) # estimates
  }
  mleGPAR(y.training)
  # distribution of Y_(n+h) given Y_n = (p_h(0), p_h(1), ..., p_h(l))
  
  predGPAR = function(yn, h, par){
    phi_hat=par[1]
    th_hat=par[2]
    J_eps=((id-phi_hat*t(ker_mat(th_hat)))*(1/(1-phi_hat)))%*%(J/n_s)
    l=tr
    return(sapply(0:l, function(j) conditionalGPAR(j, yn, h, phi_hat, ker_mat(th_hat),
                                                   J_eps)))
  }  
  
  # predicted error for a single data set
  
  predictedErrorGPAR = function(u, training.set, test.set){ 
    n = length(training.set)
    m = length(test.set)
    
    predictedErrorGPAR_h = function(u, h, phi, test.set){         
      pred.dist = sapply(1:(m-h+1), function(j) predGPAR(u[n+j-1], h, phi))
      mode = sapply(1:(m-h+1), function(j) mode.pred(pred.dist[,j]))
      return(c(ptp(test.set[h:m], mode)))  
    }      
    
    mle = mleGPAR(training.set)          
    return( t(sapply(1:3, function(h) 
      predictedErrorGPAR_h(u, h, mle, test.set))) )
  }
  
# storing the PTPs
round(predictedErrorMTD(y,y.training,y.test),2) # for MTD
round(predictedErrorPAR(y,y.training,y.test),2) # for PAR
round(predictedErrorGPAR(y,y.training,y.test),2) # for GPAR





##### ML MODELS #####
dml <- d
dml$Category <- factor(dml$Category, levels = 0:5, ordered = T)

### Ordinal Logistic Regression ###
model_data <- data.frame(Predictor = dml[-nrow((dml)), ]$Category,
                         Response = dml[-1, ]$Category)

# splitting train-test set
n <- nrow(model_data)
train_data <- model_data[1:(n - ss), ]
test_data <- model_data[(n - (ss - 1)):n, ]

# function for computing TPM
compute_TPM <- function(model, levels = 0:5) {
  n_levels <- length(levels)
  Q <- matrix(0, nrow = n_levels, ncol = n_levels)
  colnames(Q) <- rownames(Q) <- levels
  for (i in levels) {
    new_data <- data.frame(Predictor = factor(i, levels = levels, ordered = TRUE))
    probs <- predict(model, newdata = new_data, type = "probs")
    Q[i + 1, ] <- probs
  }
  return(Q)
}
  
# function for getting PTP of h-step ahead forecasting
predictedErrorLogistic_h <- function(Q, test_data, h) {
  Q_h <- matrix.power(Q, h)
  test_predictor_indices <- as.integer(as.character(test_data$Predictor)) + 1
  predicted <- as.numeric(apply(Q_h[test_predictor_indices, ], 1, which.max) - 1)
  actual <- as.numeric(as.character(test_data$Response))
  return(ptp(actual[(h + 1):length(actual)], predicted[1:(length(predicted) - h)]))
}

# function for getting PTP for a single dataset
predictedErrorLogistic <- function(train_data, test_data) {
  
  # Fit the ordinal logistic regression model
  log_model <- polr(Response ~ Predictor, data = train_data, Hess = TRUE)
  
  # Compute Q
  Q <- compute_TPM(log_model, levels = 0:5)
  
  # Return accuracies for h = 1, 2, 3
  return( t(sapply(1:3, function(h) 
    predictedErrorLogistic_h(Q, test_data, h))) )
}

# storing PTP's
round(predictedErrorLogistic(train_data, test_data), 2)



### RANDOM FOREST MODEL ###
model_data <- data.frame(Predictor1 = lag(dml$Category, n = 1), 
                         Predictor2 = lag(dml$Category, n = 2), 
                         Predictor3 = lag(dml$Category, n = 3), 
                         Predictor4 = lag(dml$Category, n = 4), 
                         Predictor5 = lag(dml$Category, n = 5), 
                         Response = dml$Category)
model_data <- na.omit(model_data)

# splitting train-test set
n <- nrow(model_data)
train_data <- model_data[1:(n - ss), ]
test_data <- model_data[(n - (ss - 1)):n, ]

# fitting rf
#set.seed(123)
rf_model_1 <- randomForest(Response ~ ., data = train_data)
rf_model_2 <- randomForest(Response ~ Predictor2 + Predictor3 + Predictor4 + Predictor5, data = train_data)
rf_model_3 <- randomForest(Response ~ Predictor3 + Predictor4 + Predictor5, data = train_data)

# storing PTP's
pred_1 <- as.numeric(predict(rf_model_1, test_data)) - 1
pred_2 <- as.numeric(predict(rf_model_2, test_data)) - 1
pred_3 <- as.numeric(predict(rf_model_3, test_data)) - 1

round(ptp(as.numeric(test_data$Response) - 1, pred_1), 2)
round(ptp(as.numeric(test_data$Response) - 1, pred_2), 2)
round(ptp(as.numeric(test_data$Response) - 1, pred_3), 2)



### SUPPORT VECTOR MACHINE ###
model_data <- data.frame(Predictor1 = lag(dml$Category, n = 1), 
                         Predictor2 = lag(dml$Category, n = 2), 
                         Predictor3 = lag(dml$Category, n = 3), 
                         Predictor4 = lag(dml$Category, n = 4), 
                         Predictor5 = lag(dml$Category, n = 5), 
                         Response = dml$Category)
model_data <- na.omit(model_data)

# splitting train-test set
n <- nrow(model_data)
train_data <- model_data[1:(n - ss), ]
test_data <- model_data[(n - (ss - 1)):n, ]

# fitting svm
svm_model_1 <- svm(Response ~ ., data = train_data, type = 'C')
svm_model_2 <- svm(Response ~ Predictor2 + Predictor3 + Predictor4 + Predictor5, data = train_data, type = 'C')
svm_model_3 <- svm(Response ~ Predictor3 + Predictor4 + Predictor5, data = train_data, type = 'C')

# storing PTP's
pred_1 <- as.numeric(predict(svm_model_1, test_data)) - 1
pred_2 <- as.numeric(predict(svm_model_2, test_data)) - 1
pred_3 <- as.numeric(predict(svm_model_3, test_data)) - 1

round(ptp(as.numeric(test_data$Response) - 1, pred_1), 2)
round(ptp(as.numeric(test_data$Response) - 1, pred_2), 2)
round(ptp(as.numeric(test_data$Response) - 1, pred_3), 2)
