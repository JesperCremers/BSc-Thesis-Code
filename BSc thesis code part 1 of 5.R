# Code Comparing the Partial Measures of Causlity of Hosoya and Geweke in Simulated VAR models 


# Importing the libraries
library(vars)

# This sets various graphical parameters.
setpar<-function(...) {
  par(mar=c(2.5,2.5,0.5,0.5),mgp=c(1.4,0.2,0), # margins and distance
      cex=0.7,bty="l",las=1,lwd=0.5,tck=0.01,...)  # text size 
}

# create a function that generates samples from a multivariate normal distribution.
# n = number of samples, m = mean vector, S = covariance matrix
rmultinorm<-function(n,m,S) {
  d<-ifelse(is.null(nrow(m)),length(m),nrow(m))
  m+t(chol(S))%*%matrix(rnorm(d*n),d,n)
}

# r: decay rate, the magnitude of the autocorrelation decreases by 1.6 for each lag unit
# theta: phase shift, strength and timing of the sinusoidal components at each lag
ar.coeff<-function(r,theta) {    
  p<-sum(ifelse(theta==0,1,0)+ifelse(theta==pi,1,0)) # returns 1 if theta = 0 or returns 1 if theta = pi. 
  p<-p+2*(length(theta)-p)                           # number of basis functions needed for the expansion of the radial function.
  phi<-rep(0,p+1)                                    # create vector of length p+1
  phi[1]<-1                                          # set first element equal to 1
  for (i in (1:length(theta))) {                     # loop through theta
    if ((theta[i]==0)||(theta[i]==pi)) {             #
      q<-c(1,-1/r[i],0)
    } else {
      q<-c(1,-2/r[i]*cos(theta[i]),1/r[i]^2)	
    }
    phiold<-phi
    phi[2:(p+1)]<-phi[2:(p+1)]+q[2]*phiold[1:(p)]
    phi[3:(p+1)]<-phi[3:(p+1)]+q[3]*phiold[1:(p-1)]
  }
  phi<--phi[2:(p+1)]
  return(phi)
}

#### Calculate spectral density matrix ####
var.spec=function(f,A,S) {                   # takes the arguments f, A, S
  D=nrow(A)                                  # D is number of rows in A  (D=3)
  p=ncol(A)/D                                # p is number of columns of A divided by D  (p=2) 
  Q=array(0,c(D,D,length(f)))                # Create array of dimensions D by D by length(f) (hence 100 3x3 matrices)
  F=array(0,c(D,D,length(f)))                # Create array of dimensions D by D by length(f)
  I=array(0,c(D,D,length(f)))                # Create array of dimensions D by D by length(f)
  for (j in 1:D) {
    Q[j,j,]=1  # for all matrices, the entry on the first row and first column are 1, because the AR(0) is just the identity matrix
    for (k in 1:D) {
      for (i in (1:p)) {
        Q[j,k,]=Q[j,k,]-A[j,(i-1)*D+k]*exp(complex(imaginary=f*i)) #here i is i from the loop. we have e^(imaginary*f)^p
      }   # note that in the above we have Q[j,k,]- A[...]... because we are here dealing with the AR case, it should be negative to have it mean reverting and hence stationary. the charasteristic equation is 1-thetaX(t-1)-theta^2X(t-2) and it it still on the left side, when we then inverse this matrix, it becomes positive and then it is on the right side.
    }
  }
  for (k in 1:length(f)) {
    I[,,k]=solve(Q[,,k])           # inverse because we have an VAR(2) model, hence the matrix is on the left side and we take it to the right side
    F[,,k]=I[,,k]%*%S%*%Conj(t(I[,,k]))/(2*pi)            # WHY MULTIPLY WITH S, we have to include the covariance structure of the time series. Because we specified non-zero covariance, if the covariance wasa 0 this might not have been neccesary.
  }
  return(list(spec=F,factor=I))      # we return the spectral density, and the factor function.
}

# Matrix exponentiation with matrix x and exponent n
"%^%" <- function(x, n) 
  with(eigen(x), vectors %*% (values^n * t(vectors)))

#### Hosoya code in a method ####
hosoyaTriv=function(f,A,S,i,j) {
  Lambda1=var.spec(f,A,S)$factor  #calculate the factor matrix of the spectral decomposition. Recall lambda here is the big lambda from the factorization.
  I=c(i,j)  # vector containing the specified indices
  J=(1:d)[-I]  # integer containing the not specified index
  P1=diag(d)  # P1 is the identity  matrix of dimensions 3x3
  P1[J,I]=-S[J,I]%*%solve(S[I,I])  # calculates -Sigma_3. x Sigma_..^-1, hence we get that P1 is the right matrix
  P2=matrix(0,d,d)       # specify a 3x3 matrix 
  P2[I,I]=S[I,I]%^%(-0.5)     # this is the same as asked in the paper, it is Sigma_..^-1/2
  P2[J,J]=(S[J,J]-S[J,I]%*%solve(S[I,I])%*%S[J,I])%^%(-0.5)  #this is simply (Sigma_33:-) ^-1/2
  P=P2%*%P1      # we multiply both to obtain the full matrix. Here P is PI
  Lambda2=array(0,dim(Lambda1))   # 100 times (frequencies) a 3x3 matrix (3 series)
  for (k in 1:length(f)) {   # obtain lambda~ = lambda * P1^-1 for all frequencies
    Lambda2[,,k]=Lambda1[,,k]%*%solve(P1)          #  note that we can keep Lambda1(e^(if) as we are interested in this quantity at the end too)
  }
  Lambda3=Lambda2[I,I,] # lambda3 is lambda~_..
  Xi=diag(2)  # create a 2x2 identity matrix
  Xi[2,1]=-S[j,i]/S[i,i]  # this is for the XI -> -Sigma_21 x Sigma_11^-1. 
  Gamma=array(0,dim(Lambda3)) # create 2x2 matrices f times
  for (k in 1:length(f)) { # we do it for all frequencies
    Gamma[,,k]=Lambda3[,,k]%*%solve(Xi)   # is the last part on page 3 of my notes
  }
  PM=log(1+abs(Gamma[1,2,])^2/abs(Gamma[1,1,])^2) #this is the partial measure
  return(PM)
}

#### Geweke code in a method ####
geweke = function(X,i,j,k,f) { #effect of j on i, given k
  dd_1 = cbind(X[,i], X[,k])   
  dd_2 = cbind(X[,i], X[,j], X[,k])
  
  model1 = VAR(dd_1, p = 2)
  coef_model1 = array(0, dim=c(model1$K, model1$K, model1$p))
  add = array(0, dim=c(model1$K,model1$K,model1$p,length(f)))
  ADD = array(0,dim=c(model1$K,model1$K,length(f)))
  ADD_x1 = array(0,dim=c(model1$K,model1$K,length(f)))
  H1_1=array(0,dim=c(model1$K,model1$K,length(f)))
  
  model2 = VAR(dd_2, p=2)
  coef_model2=array(0,dim=c(model2$K,model2$K,model2$p))
  add_2<-array(0,dim=c(model2$K,model2$K,model2$p,length(f)))
  ADD_2=array(0,dim=c(model2$K,model2$K,length(f)))
  ADD_x2=array(0,dim=c(model2$K,model2$K,length(f)))
  H1_2=array(0,dim=c(model2$K,model2$K,length(f)))
  
  G_2<-array(0,dim=c(model2$K,model2$K,length(f)))
  Q1<-array(0,dim=c(model2$K,model2$K,length(f)))
  
  for (p in 1:model1$K){
    for (k in 1:model1$p){
      coef_model1[1,p,k]=coef(model1)$y1[p+(k-1)*model1$K,1]}}
  for (p in 1:model1$K){
    for (k in 1:model1$p){
      coef_model1[2,p,k]=coef(model1)$y2[p+(k-1)*model1$K,1]}}
  for (p in 1:model2$K){
    for (k in 1:model2$p){
      coef_model2[1,p,k]=coef(model2)$y1[p+(k-1)*model2$K,1]}}
  for (p in 1:model2$K){
    for (k in 1:model2$p){
      coef_model2[2,p,k]=coef(model2)$y2[p+(k-1)*model2$K,1]}}
  for (p in 1:model2$K){
    for (k in 1:model2$p){
      coef_model2[3,p,k]=coef(model2)$y3[p+(k-1)*model2$K,1]}}
  
  for (l in 1:length(f)) {
    for (k in 1:model1$p) {
      add[,,k,l]<-coef_model1[,,k]*exp(-k*(1i)*f[l])
    }
  }
  
  for (l in 1:length(f)) {
    for (k in 1:model1$p) {
      ADD[,,l]=ADD[,,l]+add[,,k,l]
    }
  }
  
  #Sigma_model1<-summary(model1)$cov
  Sigma_model1<-ar.yw(dd_1,aic=FALSE,order.max=2)$var.pred # for storage issues
  P_model1<-matrix(c(1,-Sigma_model1[1,2]/Sigma_model1[1,1],0,1),2,2)
  
  for (l in 1:length(f)){
    ADD_x1[,,l]=(P_model1)%*%(-ADD[,,l]+diag(x=1,model1$K,model1$K)*1)
  }
  for (l in 1:length(f)){
    H1_1[,,l]=solve(ADD_x1[,,l])
  }
  
  #Sigma_model2<-summary(model2)$cov
  Sigma_model2 <- ar.yw(dd_2, aic=FALSE, order.max=2)$var.pred # for storage issues
  
  P1_model2<-matrix(c(1,0,0,-Sigma_model2[2,1]*solve(Sigma_model2[1,1]),1,0,-Sigma_model2[3,1]*solve(Sigma_model2[1,1]),0,1),3,3,byrow=T)
  P2_model2<-matrix(c(1,0,0,0,1,0,0,-(Sigma_model2[3,2]-Sigma_model2[3,1]*solve(Sigma_model2[1,1])*Sigma_model2[1,2])*solve(Sigma_model2[2,2]-Sigma_model2[2,1]*solve(Sigma_model2[1,1])*Sigma_model2[1,2]),1),3,3,byrow=T)
  
  P_model2<-P1_model2%*%P2_model2
  
  for (l in 1:length(f)){
    for (k in 1:model2$p)
    {add_2[,,k,l]<-coef_model2[,,k]*exp(-k*(1i)*f[l])}
  }
  
  
  for (l in 1:length(f)){
    for (k in 1:model2$p){
      ADD_2[,,l]=ADD_2[,,l]+add_2[,,k,l]}
  }
  
  for (l in 1:length(f)){
    ADD_x2[,,l]=(P_model2)%*%(-ADD_2[,,l]+diag(x=1,model2$K,model2$K)*1)
  }
  
  for (l in 1:length(f)){
    H1_2[,,l]=solve(ADD_x2[,,l])}
  
  for (l in 1:length(f)){
    G_2[1,1,l]<-H1_1[1,1,l]
    G_2[1,3,l]<-H1_1[1,2,l]
    G_2[3,1,l]<-H1_1[2,1,l]
    G_2[3,3,l]<-H1_1[2,2,l]
    G_2[2,2,l]<-1
  }
  
  for (l in 1:length(f)){
    Q1[,,l]<-solve(G_2[,,l])%*%H1_2[,,l]
  }
  
  astr_too1<-vector("numeric",length(f))
  astr_too2<-vector("numeric",length(f))
  astr_too3<-vector("numeric",length(f))
  astr_too<-vector("numeric",length(f))
  boh_too<-vector("numeric",length(f))
  for (l in 1:length(f)){
    astr_too1[l]<-Q1[1,1,l]*Sigma_model2[1,1]*Conj(Q1[1,1,l])
    astr_too2[l]<-Q1[1,2,l]*Sigma_model2[2,2]*Conj(Q1[1,2,l])
    astr_too3[l]<-Q1[1,3,l]*Sigma_model2[3,3]*Conj(Q1[1,3,l])
    astr_too[l]<-astr_too1[l]+astr_too2[l]+astr_too3[l]
    boh_too[l]<-log(abs(astr_too[l])/abs(astr_too1[l]))
  }
  return(boh_too)
}


#### DGP Simulation ####
T=5000 # real sample
B=100  # burn in part
A=array(0,c(2,3,3)) # 2 lag matrices, both matrices have dimension 3x3 (X,Y,Z), actually, we should have an extra 0'th matrix for the coefficients from the not lags, however, this isan identity matrix and is captured in the lines for the spectal matrix (lines 24-26)
A[,1,1]= c(0.7, -0.4)      #ar.coeff(1.6,0.9)
A[,2,2]= c(0.8, -0.5)      #ar.coeff(1.5,1.2)
A[,3,3]= c(0.6, -0.4)     #ar.coeff(1.3,1.9)
A[1,3,2]=0.7
A[1,1,3]=0.6
S=diag(3)
S[1,1]=0.8
S[2,2]=0.7
S[3,3]=0.6

p=2 # we model a VAR(2) process, that is we model for 2 lags.
d=3 # we specify 3  time series
X=matrix(rmultinorm((T+B),c(0,0,0),S),T+B,3,byrow=T) # create a matrix X with dimensions (T+B) x 3 using the rmultinorm with mean (0,0,0) and covariance matrix S.
for (t in 3:(T+B)) {           # start at 3 because we have 2 lags.
  X[t,]=X[t,]+A[1,,]%*%X[t-1,]+A[2,,]%*%X[t-2,] # creates the AR(2) series with A[1,,] being the first lag matrix and A[2,,] being the second lag matrix. Here the first X[t,] aftere the = is that X[t,] would represent the prediction error
}
X=X[(B+1):(T+B),]     # removing burn in samples

# Specify frequency range
f = seq(0, pi, length=100)

# Calculate Hosoya parameters
res=ar.yw(X,aic=FALSE,order.max=2)
Sh = res$var.pred   # the covaariance matrix of the the residuals using this predicted VAR(2) model. THIS IS BIG SIGMA !
St=S 

Ah=matrix(0,d,d*p) #p=2, d=3
At=matrix(0,d,d*p)
for (i in 1:p) {  # extracts coefficients for lag "i" from res$ar and place them i the correct columns of Ah
  Ah[,(i-1)*d+(1:3)]=res$ar[i,,] # Ah is just a 3x6 matrix with the first 3x3 being res$ar[1,,] and the last 33 matrix being res$ar[2,,]
  At[,(i-1)*d+(1:3)]=A[i,,] # At is just a 3x6 matrix with the first 3x3 being A[1,,] and the last 33 matrix being A[2,,]
}

#### partial measures of causality for estimated model (Hosoya) ####
PM.hos=array(0,c(d,d,length(f)))
for (i in 1:d) {
  for (j in 1:d) {
    if (i==j) PM.hos[i,j,]=0
    else PM.hos[i,j,]=hosoyaTriv(f,Ah,Sh,i,j)
  }
}

#### partial measures of causality for estimated model (Geweke) ####
PM.gew=array(0,c(d,d,length(f)))
PM.gew[1,2,]=geweke(X,1,2,3,f)
PM.gew[1,3,]=geweke(X,1,3,2,f)
PM.gew[2,1,]=geweke(X,2,1,3,f)
PM.gew[2,3,]=geweke(X,2,3,1,f)
PM.gew[3,1,]=geweke(X,3,1,2,f)
PM.gew[3,2,]=geweke(X,3,2,1,f)

# plot all the measures 
setpar(mfrow=c(d,d-1))
for (i in 1:d) {
  for (j in 1:d) {
    if (j!=i) 
      plot(f,PM.hos[i,j,],type="l",col=4, main= paste(j, " to ", i),
           xlab = "frequency", ylab = "PM Hosoya")   # blue
    abline(h = qchisq(0.05, df=1), col = "red", lty = 2 )
  }
}

setpar(mfrow=c(d,d-1))
for (i in 1:d) {
  for (j in 1:d) {
    if (j!=i) 
      plot(f,PM.gew[i,j,],type="l",col=3, main= paste(j, " to ", i),
           xlab = "frequency", ylab = "PM Geweke")   # blue
  }
}



# upper and lower bounds
n = 500

X.rep = list()
res.rep = list()
A.rep = list()
A.rep.new = list()
S.rep = list()
PM1on2 = list()
PM2on1 = list()
PM1on3 = list()
PM3on1 = list()
PM2on3 = list()
PM3on2 = list()

PM1on2.G = list()
PM2on1.G = list()
PM1on3.G = list()
PM3on1.G = list()
PM2on3.G = list()
PM3on2.G = list()
for (i in 1:n){
  X.rep[[i]]=matrix(rmultinorm((T+B),c(0,0,0),S),T+B,3,byrow=T) # create a matrix X with dimensions (T+B) x 3 using the rmultinorm with mean (0,0,0) and covariance matrix S.
  for (t in 3:(T+B)) {           # start at 3 because we have 2 lags.
    X.rep[[i]][t,]=as.matrix(X.rep[[i]][t,]+A[1,,]%*%X.rep[[i]][t-1,]+A[2,,]%*%X.rep[[i]][t-2,]) # creates the AR(2) series with A[1,,] being the first lag matrix and A[2,,] being the second lag matrix. Here the first X[t,] aftere the = is that X[t,] would represent the prediction error
  }
  X.rep[[i]]=X.rep[[i]][(B+1):(T+B),]
  res.rep[[i]] = ar.yw(X.rep[[i]],aic=FALSE,order.max=2)
  A.rep[[i]]= res.rep[[i]]$ar
  S.rep[[i]]= res.rep[[i]]$var.pred
  A.rep.new[[i]]=matrix(0,d,d*p)
  for (c in 1:p) {
    A.rep.new[[i]][,(c-1)*d+(1:3)]= A.rep[[i]][c,,]
  }
  PM1on2[[i]]= hosoyaTriv(f, A.rep.new[[i]], S.rep[[i]],2,1)
  PM2on1[[i]]= hosoyaTriv(f, A.rep.new[[i]], S.rep[[i]],1,2)
  PM1on3[[i]]= hosoyaTriv(f, A.rep.new[[i]], S.rep[[i]],3,1)
  PM3on1[[i]]= hosoyaTriv(f, A.rep.new[[i]], S.rep[[i]],1,3)
  PM2on3[[i]]= hosoyaTriv(f, A.rep.new[[i]], S.rep[[i]],3,2)
  PM3on2[[i]]= hosoyaTriv(f, A.rep.new[[i]], S.rep[[i]],2,3)
  
  PM1on2.G[[i]] = geweke(X.rep[[i]],2,1,3,f)
  PM2on1.G[[i]] = geweke(X.rep[[i]],1,2,3,f)
  PM1on3.G[[i]] = geweke(X.rep[[i]],3,1,2,f)
  PM3on1.G[[i]] = geweke(X.rep[[i]],1,3,2,f)
  PM2on3.G[[i]] = geweke(X.rep[[i]],3,2,1,f)
  PM3on2.G[[i]] = geweke(X.rep[[i]],2,3,1,f)
}

# loop for matrix of means
PM.rep.mean=array(0,c(d,d,length(f)))
for (i in 1:d){
  PM.rep.mean[i,i,]=0
}
for (i in 1:length(f)){
  PM.rep.mean[1,2,i] = mean((sapply(PM2on1, "[[", i)), probs = 0.975)
  PM.rep.mean[1,3,i] = mean((sapply(PM3on1, "[[", i)), probs = 0.975)
  PM.rep.mean[2,1,i] = mean((sapply(PM1on2, "[[", i)), probs = 0.975)
  PM.rep.mean[2,3,i] = mean((sapply(PM3on2, "[[", i)), probs = 0.975)
  PM.rep.mean[3,1,i] = mean((sapply(PM1on3, "[[", i)), probs = 0.975)
  PM.rep.mean[3,2,i] = mean((sapply(PM2on3, "[[", i)), probs = 0.975)
}

PM.rep.mean.G=array(0,c(d,d,length(f)))
for (i in 1:d){
  PM.rep.mean.G[i,i,]=0
}
for (i in 1:length(f)){
  PM.rep.mean.G[1,2,i] = mean((sapply(PM2on1.G, "[[", i)), probs = 0.975)
  PM.rep.mean.G[1,3,i] = mean((sapply(PM3on1.G, "[[", i)), probs = 0.975)
  PM.rep.mean.G[2,1,i] = mean((sapply(PM1on2.G, "[[", i)), probs = 0.975)
  PM.rep.mean.G[2,3,i] = mean((sapply(PM3on2.G, "[[", i)), probs = 0.975)
  PM.rep.mean.G[3,1,i] = mean((sapply(PM1on3.G, "[[", i)), probs = 0.975)
  PM.rep.mean.G[3,2,i] = mean((sapply(PM2on3.G, "[[", i)), probs = 0.975)
}


# loop for matrix of upper bounds
PM.rep.upper=array(0,c(d,d,length(f)))
for (i in 1:d){
  PM.rep.upper[i,i,]=0
}
for (i in 1:length(f)){
  PM.rep.upper[1,2,i] = quantile((sapply(PM2on1, "[[", i)), probs = 0.975)
  PM.rep.upper[1,3,i] = quantile((sapply(PM3on1, "[[", i)), probs = 0.975)
  PM.rep.upper[2,1,i] = quantile((sapply(PM1on2, "[[", i)), probs = 0.975)
  PM.rep.upper[2,3,i] = quantile((sapply(PM3on2, "[[", i)), probs = 0.975)
  PM.rep.upper[3,1,i] = quantile((sapply(PM1on3, "[[", i)), probs = 0.975)
  PM.rep.upper[3,2,i] = quantile((sapply(PM2on3, "[[", i)), probs = 0.975)
}

PM.rep.upper.G=array(0,c(d,d,length(f)))
for (i in 1:d){
  PM.rep.upper.G[i,i,]=0
}
for (i in 1:length(f)){
  PM.rep.upper.G[1,2,i] = quantile((sapply(PM2on1.G, "[[", i)), probs = 0.975)
  PM.rep.upper.G[1,3,i] = quantile((sapply(PM3on1.G, "[[", i)), probs = 0.975)
  PM.rep.upper.G[2,1,i] = quantile((sapply(PM1on2.G, "[[", i)), probs = 0.975)
  PM.rep.upper.G[2,3,i] = quantile((sapply(PM3on2.G, "[[", i)), probs = 0.975)
  PM.rep.upper.G[3,1,i] = quantile((sapply(PM1on3.G, "[[", i)), probs = 0.975)
  PM.rep.upper.G[3,2,i] = quantile((sapply(PM2on3.G, "[[", i)), probs = 0.975)
}

# loop for matrix of lower bounds
PM.rep.lower=array(0,c(d,d,length(f)))
for (i in 1:d){
  PM.rep.lower[i,i,]=0
}
for (i in 1:length(f)){
  PM.rep.lower[1,2,i] = quantile((sapply(PM2on1, "[[", i)), probs = 0.025)
  PM.rep.lower[1,3,i] = quantile((sapply(PM3on1, "[[", i)), probs = 0.025)
  PM.rep.lower[2,1,i] = quantile((sapply(PM1on2, "[[", i)), probs = 0.025)
  PM.rep.lower[2,3,i] = quantile((sapply(PM3on2, "[[", i)), probs = 0.025)
  PM.rep.lower[3,1,i] = quantile((sapply(PM1on3, "[[", i)), probs = 0.025)
  PM.rep.lower[3,2,i] = quantile((sapply(PM2on3, "[[", i)), probs = 0.025)
}

PM.rep.lower.G=array(0,c(d,d,length(f)))
for (i in 1:d){
  PM.rep.lower.G[i,i,]=0
}
for (i in 1:length(f)){
  PM.rep.lower.G[1,2,i] = quantile((sapply(PM2on1.G, "[[", i)), probs = 0.025)
  PM.rep.lower.G[1,3,i] = quantile((sapply(PM3on1.G, "[[", i)), probs = 0.025)
  PM.rep.lower.G[2,1,i] = quantile((sapply(PM1on2.G, "[[", i)), probs = 0.025)
  PM.rep.lower.G[2,3,i] = quantile((sapply(PM3on2.G, "[[", i)), probs = 0.025)
  PM.rep.lower.G[3,1,i] = quantile((sapply(PM1on3.G, "[[", i)), probs = 0.025)
  PM.rep.lower.G[3,2,i] = quantile((sapply(PM2on3.G, "[[", i)), probs = 0.025)
}


# plot all the measures 
setpar(mfrow=c(d,d-1))
par(mar = c(3.7, 3.7, 1.7, 1.7))
labels <- c("x", "y", "z")
for (i in 1:d) {
  for (j in 1:d) {
    if (j != i) {
      y_max <- max(c(PM.rep.upper[i, j, ], PM.rep.mean[i, j, ], PM.rep.lower[i, j, ]))
      y_min <- min(c(PM.rep.upper[i, j, ], PM.rep.mean[i, j, ], PM.rep.lower[i, j, ]))
      
      plot(f, PM.rep.upper[i, j, ], type = "l", main = "",
           xlab = "frequency", ylab = "", lty = 3, ylim = c(y_min, y_max))
      text(x = max(f), y = y_max, labels = paste(labels[j], " to ", labels[i]), pos = 2, cex = 1.2)
    }
    
    lines(f, PM.rep.mean[i, j, ])
    lines(f, PM.rep.lower[i, j, ], lty = 3)
  }
}

setpar(mfrow = c(d, d-1))
par(mar = c(3.7, 3.7, 1.7, 1.7))
labels <- c("x", "y", "z")
for (i in 1:d) {
  for (j in 1:d) {
    if (j != i) {
      y_max <- max(c(PM.rep.upper.G[i, j, ], PM.rep.mean.G[i, j, ], PM.rep.lower.G[i, j, ]))
      y_min <- min(c(PM.rep.upper.G[i, j, ], PM.rep.mean.G[i, j, ], PM.rep.lower.G[i, j, ]))
      
      plot(f, PM.rep.upper.G[i, j, ], type = "l", lty = 3, main = "",
           xlab = "frequency", ylab = "", ylim = c(y_min, y_max))
      text(x = max(f), y = y_max, labels = paste(labels[j], " to ", labels[i]), pos = 2, cex = 1.2)
    }
    
    lines(f, PM.rep.mean.G[i, j, ])
    lines(f, PM.rep.lower.G[i, j, ], lty = 3)
  }
}


