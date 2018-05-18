#setwd("~/Desktop/DTCC_Liquidity")
setwd("C:/Users/rstrepparava/Documents/Projects_2017/DTCC_liquidity_model_development/DTCC_liquidity_results/SARIMAX_Model")

library(forecast)
library(CombMSC)
library(tseries)
library(expm)
library(matrixcalc)

################################################
######      read the data     ######
################################################
data = read.csv("cover1.csv",header=TRUE)
data = data[,c(2,3,5,6,7)]

data[,1] = as.Date(data[,1])
total_raw = data[,3] + data[,4] + data[,5]
new = data.frame(data[,1],total_raw)
colnames(new) = c("date","total_raw")

####remove "Sandy"
sandy = which(new[,2]==0)
new = new[-206,]

plot(new,type="l",main="true data",ylab="RLN")

#new = ts(new)
##### center
new[,2] = as.vector(scale(new[,2],center=TRUE,scale=FALSE))/(10^10)

plot(new,type="l",main="true data",ylab="centered RLN")






################################################
######      seperate seasonal component     ######
################################################
first = 500
total = nrow(new)-first+1
sigmavec = rep(0,total)

#for (i in 1:length(sigmavec)){
#  print(i)
#  nstart = first+(i-1)
#  newcut = ts(new[1:nstart,2],frequency = 21)
  
  for (i in 1:total){
    print(i)
    nstart = first+(i-1)
    newcut = ts(new[i:nstart,2],frequency = 21)
  ################################################
  ######      seperate seasonal component     ######
  ################################################
  dnew = stl(newcut, "periodic")
  seasonal   <- dnew$time.series[,1]
  trend     <- dnew$time.series[,2]
  random  <- dnew$time.series[,3]
  
  dseason = newcut-seasonal
  #plot(dseason,type="l")
  
  dseason2 = ts(dseason,frequency=1)
  #plot(dseason2,type="l")
  
  ######################################################
  ######      represent into AR infinity model     ######
  ######################################################
  #armafit <- arma(dseason2, order = c(1, 1))
  armafit = arima(dseason2, order=c(1,1,3))
  
  #########Method1: directly use the estimated arma coefficients
  
  arorder = armafit$coef[1]
  maorder = armafit$coef[2:4]
  
  #########Method2: use the values from acf and pacf (not correct)
  
  #Switch the phi and the theta in the ARMAtoMA module to get the PIs of the infinite AR process
  
  n = 10
  coeff = ARMAtoMA(ar = c(maorder), ma = c(arorder), lag.max=n)
  
  
  ######################################################
  ######      estimate epsilon t     ######
  ######################################################
  
  x = new[i:nstart,2]
  
  ct = seasonal[length(seasonal)]
  pastn = as.matrix(x[(length(seasonal)-n):(length(seasonal)-1)])
  et = x[length(x)]-t(as.matrix(coeff))%*%pastn-ct
  
  
  ######################################################
  ######      construct vector and matrix     ######
  ######################################################
  
  cvec = as.matrix(c(ct,rep(0,(n-1))))
  C = rbind(coeff,cbind(diag((n-1)),0))
  evec = as.matrix(c(et,rep(0,(n-1))))
  
  
  ######################################################
  ######      estimation     ######
  ######################################################
  #H = -logm(C,method="Higham08")
  H = -logm(C,method="Eigen")
  mu = solve(diag(n)-C) %*% cvec
  StS = C %*% cov(evec %*% t(evec)) %*% t(C)
  
  ###vec function returns a column vector that is a stack of the columns of x, an m by n matrix.
  
  middle = kronecker(H,diag(n))+kronecker(diag(n),H)
  vecS = solve(middle) %*% vec(StS)
  
  S = matrix(vecS,nrow=n,byrow=FALSE)
  eigens = eigen(S)$values
  
  lambda = sqrt(max(eigens[eigens>0]))
  
  
  
  ######################################################
  ######      final    ######
  ######################################################
  #mufinal = mu[1,]
  ##sigmafinal = vecS[1,]
  sigmavec[i] = lambda
  
}


################################################################################################
################################################################################################
################################################################################################
######      divide by 10^10 and center it   ######
################################################################################################
################################################################################################
################################################################################################

################################################
######      read the data     ######
################################################
################################################
################################################

data = read.csv("cover1.csv",header=TRUE)
data = data[,c(2,3,5,6,7)]

data[,1] = as.Date(data[,1])
total_raw = data[,3] + data[,4] + data[,5]
new = data.frame(data[,1],total_raw)
colnames(new) = c("date","total_raw")

####remove "Sandy"
sandy = which(new[,2]==0)
new = new[-206,]

#### divide by 10^10
#new[,2] = new[,2]/(10^10)

#plot(new,type="l",main="true data",ylab="RLN in 10 billion")

##### center
new[,2] = as.vector(scale(new[,2],center=TRUE,scale=FALSE))/(10^10)

plot(new,type="l",main="true data",ylab="centered RLN")

################################################
######      seperate seasonal component     ######
################################################
newts = ts(new[,2], frequency = 21)
#dsimu = decompose(simu, "additive")
dnew = stl(newts, "periodic")
seasonal   <- dnew$time.series[,1]
trend     <- dnew$time.series[,2]
random  <- dnew$time.series[,3]

#plot(as.ts(dsimu$seasonal))
#plot(as.ts(dsimu$trend))
#plot(as.ts(dsimu$random))
#plot(dsimu)


plot(seasonal)
plot(trend)
plot(random)

#seasonal = as.ts(dsimu$seasonal)
#trend = as.ts(dsimu$trend)
#random = as.ts(dsimu$random)

dseason = newts-seasonal
plot(dseason,type="l")

dseason2 = ts(dseason,frequency=1)
plot(dseason2,type="l")

######################################################
######      represent into AR infinity model     ######
######################################################
#armafit <- arma(dseason2, order = c(0,3))
armafit = arima(dseason2, order=c(1,1,3))

#########Method1: directly use the estimated arma coefficients

arorder = armafit$coef[1]
maorder = armafit$coef[2:4]

#########Method2: use the values from acf and pacf (not correct)

#Switch the phi and the theta in the ARMAtoMA module to get the PIs of the infinite AR process

n = 10
coeff = ARMAtoMA(ar = c(maorder), ma = c(arorder), lag.max=n)


######################################################
######      estimate epsilon t     ######
######################################################

x = new[,2]

ct = seasonal[length(seasonal)]
pastn = as.matrix(x[(length(seasonal)-n):(length(seasonal)-1)])
et = x[length(x)]-t(as.matrix(coeff))%*%pastn-ct



######################################################
######      construct vector and matrix     ######
######################################################

cvec = as.matrix(c(ct,rep(0,(n-1))))
C = rbind(coeff,cbind(diag((n-1)),0))
evec = as.matrix(c(et,rep(0,(n-1))))


######################################################
######      estimation     ######
######################################################
H = -logm(C,method="Eigen")
#H = -logm(C,method="Higham08")
mu = solve(diag(n)-C) %*% cvec
StS = C %*% cov(evec %*% t(evec)) %*% t(C)

###vec function returns a column vector that is a stack of the columns of x, an m by n matrix.

middle = kronecker(H,diag(n))+kronecker(diag(n),H)
vecS = solve(middle) %*% vec(StS)

S = matrix(vecS,nrow=n,byrow=FALSE)
eigens = eigen(S)$values

lambda = sqrt(max(eigens[eigens>0]))



######################################################
######      final    ######
######################################################
mufinal = mu[1,]
##sigmafinal = vecS[1,]
sigma = lambda

muvec = rep(mufinal,length(seasonal))

#####plot:
plot(new[,2],type="l",main="True data with 95%",ylab="RLN in 10 billion",ylim=c(-3,3))
lines(muvec,col="red",type="l")
lines(muvec+qnorm(0.975)*sigma,col="blue",type="l")
lines(muvec-qnorm(0.975)*sigma,col="blue",type="l")

###95% confidence interval
mufinal = mu[1,]
muvec = rep(mufinal,length(seasonal))

plot(new[,1:2],type="l",main="True data with 95% Confidence Interval (max error)",ylab="RLN in 10 billion",ylim=c(-3,3))
lines(cbind(new[,1],muvec),col="red",type="l")
lines(cbind(new[,1],muvec+qnorm(0.975)*max(sigmavec)),col="blue",type="l")
lines(cbind(new[,1],muvec-qnorm(0.975)*max(sigmavec)),col="blue",type="l")


###95% confidence interval

mufinal = mu[1,]
muvec = rep(mufinal,length(seasonal))

plot(new[,1:2],type="l",main="True data with 95% Confidence Interval (95% Percentile)",ylab="RLN in 10 billion",ylim=c(-3,3))
lines(cbind(new[,1],muvec),col="red",type="l")
lines(cbind(new[,1],muvec+qnorm(0.975)*as.numeric(quantile(sigmavec, 0.95))),col="blue",type="l")
lines(cbind(new[,1],muvec-qnorm(0.975)*as.numeric(quantile(sigmavec, 0.95))),col="blue",type="l")


try1 = c(rep(0,499),sigmavec)

plot(new[,2],type="l",main="True data with 95%",ylab="RLN in 10 billion",ylim=c(-3,3))
lines(muvec,col="red",type="l")
lines(muvec+qnorm(0.975)*try1,col="blue",type="l")
lines(muvec-qnorm(0.975)*try1,col="blue",type="l")

############################################################################################################
######      backtest     ############################################################
############################################################################################################
new = data.frame(new)
new$year = as.numeric(format(new$date,"%y"))

#table(new$year)

sigvecmax95 = rep(0,5)
sigvecmax70 = rep(0,5)
sigvecmax = rep(0,5)

subindex = c("12","13","14","15","16")

for (k in 1:(length(subindex))){
if (k==length(subindex)){
  sub1 = subset(new,(year=="16")|(year=="17"))
} else {
  
yearn = subindex[k]  
sub1 = subset(new,year==yearn)
}
  
first = 220
total = nrow(sub1)-first+1
sigmavec = rep(0,total)

#for (i in 1:length(sigmavec)){
#  print(i)
#  nstart = first+(i-1)
#  newcut = ts(new[1:nstart,2],frequency = 21)

for (i in 1:total){
 # print(i)
  nstart = first+(i-1)
  subcut = ts(sub1[i:nstart,2],frequency = 21)
  ################################################
  ######      seperate seasonal component     ######
  ################################################
  dnew = stl(subcut, "periodic")
  seasonal   <- dnew$time.series[,1]
  trend     <- dnew$time.series[,2]
  random  <- dnew$time.series[,3]
  
  dseason = subcut-seasonal
  #plot(dseason,type="l")
  
  dseason2 = ts(dseason,frequency=1)
  #plot(dseason2,type="l")
  
  ######################################################
  ######      represent into AR infinity model     ######
  ######################################################
  #armafit <- arma(dseason2, order = c(1, 1))
  armafit = arima(dseason2, order=c(1,1,3))
  
  #########Method1: directly use the estimated arma coefficients
  
  arorder = armafit$coef[1]
  maorder = armafit$coef[2:4]
  
  #########Method2: use the values from acf and pacf (not correct)
  
  #Switch the phi and the theta in the ARMAtoMA module to get the PIs of the infinite AR process
  
  n = 10
  coeff = ARMAtoMA(ar = c(maorder), ma = c(arorder), lag.max=n)
  
  
  ######################################################
  ######      estimate epsilon t     ######
  ######################################################
  
  x = sub1[i:nstart,2]
  
  ct = seasonal[length(seasonal)]
  pastn = as.matrix(x[(length(seasonal)-n):(length(seasonal)-1)])
  et = x[length(x)]-t(as.matrix(coeff))%*%pastn-ct
  
  
  ######################################################
  ######      construct vector and matrix     ######
  ######################################################
  
  cvec = as.matrix(c(ct,rep(0,(n-1))))
  C = rbind(coeff,cbind(diag((n-1)),0))
  evec = as.matrix(c(et,rep(0,(n-1))))
  
  
  ######################################################
  ######      estimation     ######
  ######################################################
  #H = -logm(C,method="Higham08")
  H = -logm(C,method="Eigen")
  mu = solve(diag(n)-C) %*% cvec
  StS = C %*% cov(evec %*% t(evec)) %*% t(C)
  
  ###vec function returns a column vector that is a stack of the columns of x, an m by n matrix.
  
  middle = kronecker(H,diag(n))+kronecker(diag(n),H)
  vecS = solve(middle) %*% vec(StS)
  
  S = matrix(vecS,nrow=n,byrow=FALSE)
  eigens = eigen(S)$values
  
  lambda = sqrt(max(eigens[eigens>0]))
  
  
  
  ######################################################
  ######      final    ######
  ######################################################
  #mufinal = mu[1,]
  ##sigmafinal = vecS[1,]
  sigmavec[i] = lambda
  
}

sigvecmax95[k] = as.numeric(quantile(sigmavec, 0.95))
sigvecmax70[k] = as.numeric(quantile(sigmavec, 0.7))
sigvecmax[k] = max(sigmavec)
#print(sigvecmax)
}

#c(as.numeric(table(new$year)))

try1 = rep(sigvecmax95,c(247,250,250,250,336))
try2 = rep(sigvecmax70,c(247,250,250,250,336))

plot(new[,1:2],type="l",main="True data with 95% (In-sample)",ylab="RLN in 10 billion",ylim=c(-3,3))
lines(cbind(new[,1],muvec),col="red",type="l")
lines(cbind(new[,1],muvec+qnorm(0.975)*try1),col="blue",type="l")
lines(cbind(new[,1],muvec-qnorm(0.975)*try1),col="blue",type="l")
lines(cbind(new[,1],muvec+qnorm(0.975)*try2),col="green",type="l")
lines(cbind(new[,1],muvec-qnorm(0.975)*try2),col="green",type="l")

legend( x="bottomleft", 
        legend=c("95 perc","70 perc"),
        col=c("blue","green"), lwd=1, lty=c(1,1),  merge=FALSE)


############################################################################################################
######      table     ############################################################
############################################################################################################
data = read.csv("cover1.csv",header=TRUE)
data = data[,c(2,3,5,6,7)]

data[,1] = as.Date(data[,1])
total_raw = data[,3] + data[,4] + data[,5]
new = data.frame(data[,1],total_raw)
colnames(new) = c("date","total_raw")

####remove "Sandy"
sandy = which(new[,2]==0)
new = new[-206,]

(mufinal+qnorm(0.975)*sigvecmax95)*(10^10)+mean(new[,2])


# ############################################################################################################
# ######      volatility     ############################################################
# ############################################################################################################
# 
# new = data.frame(new)
# new$year = as.numeric(format(new$date,"%y"))
# 
# volatility = rep(0,nrow(new))
# 
# 
# subindex = c("12","13","14","15","16")
# 
# for (k in 1:(length(subindex))){
#   print(k)
#   if (k==length(subindex)){
#     sub1 = subset(new,(year=="16")|(year=="17"))
#     index = which((new$year=="16")|(new$year=="17"))
#   } else {
#     yearn = subindex[k]  
#     sub1 = subset(new,year==yearn)
#     index = which(new$year==yearn)
#   }
#   
#   for (i in 31:length(index)){  
# volatility[index[i]] = sd(sub1[1:i,2])
# }
# }
# 
# volatility[which(volatility==0)]="NA"
# 
# index = which(volatility !="NA")
# 
# plot(volatility,col="red",type="l",main="Volatility in each year")
# 
# plot(new[index,1],volatility[index],type="l",main="True data with 95%",ylab="RLN in 10 billion")
# 
# 
# plot(new[,1:2],type="l",main="True data with 95% (Out-of-Sample)",ylab="RLN in 10 billion",ylim=c(-3,3))
# lines(cbind(new[,1],muvec),col="red",type="l")
# lines(cbind(new[248:nrow(new),1],muvec[248:nrow(new)]+qnorm(0.975)*try1),col="blue",type="l")


############################################################################################################
######      out-of-sample     ############################################################
############################################################################################################


new = data.frame(new)
new$year = as.numeric(format(new$date,"%y"))

#table(new$year)

sigvecmax95 = rep(0,5)
sigvecmax70 = rep(0,5)

subindex = c("12","13","14","15","16")

for (k in 1:(length(subindex))){
  if (k==length(subindex)){
    sub1 = subset(new,(year=="16")|(year=="17"))
  } else {
    
    yearn = subindex[k]  
    sub1 = subset(new,year==yearn)
  }
  
  first = 220
  total = nrow(sub1)-first+1
  sigmavec = rep(0,total)
  
  #for (i in 1:length(sigmavec)){
  #  print(i)
  #  nstart = first+(i-1)
  #  newcut = ts(new[1:nstart,2],frequency = 21)
  
  for (i in 1:total){
    # print(i)
    nstart = first+(i-1)
    subcut = ts(sub1[i:nstart,2],frequency = 21)
    ################################################
    ######      seperate seasonal component     ######
    ################################################
    dnew = stl(subcut, "periodic")
    seasonal   <- dnew$time.series[,1]
    trend     <- dnew$time.series[,2]
    random  <- dnew$time.series[,3]
    
    dseason = subcut-seasonal
    #plot(dseason,type="l")
    
    dseason2 = ts(dseason,frequency=1)
    #plot(dseason2,type="l")
    
    ######################################################
    ######      represent into AR infinity model     ######
    ######################################################
    #armafit <- arma(dseason2, order = c(1, 1))
    armafit = arima(dseason2, order=c(1,1,3))
    
    #########Method1: directly use the estimated arma coefficients
    
    arorder = armafit$coef[1]
    maorder = armafit$coef[2:4]
    
    #########Method2: use the values from acf and pacf (not correct)
    
    #Switch the phi and the theta in the ARMAtoMA module to get the PIs of the infinite AR process
    
    n = 10
    coeff = ARMAtoMA(ar = c(maorder), ma = c(arorder), lag.max=n)
    
    
    ######################################################
    ######      estimate epsilon t     ######
    ######################################################
    
    x = sub1[i:nstart,2]
    
    ct = seasonal[length(seasonal)]
    pastn = as.matrix(x[(length(seasonal)-n):(length(seasonal)-1)])
    et = x[length(x)]-t(as.matrix(coeff))%*%pastn-ct
    
    
    ######################################################
    ######      construct vector and matrix     ######
    ######################################################
    
    cvec = as.matrix(c(ct,rep(0,(n-1))))
    C = rbind(coeff,cbind(diag((n-1)),0))
    evec = as.matrix(c(et,rep(0,(n-1))))
    
    
    ######################################################
    ######      estimation     ######
    ######################################################
    #H = -logm(C,method="Higham08")
    H = -logm(C,method="Eigen")
    mu = solve(diag(n)-C) %*% cvec
    StS = C %*% cov(evec %*% t(evec)) %*% t(C)
    
    ###vec function returns a column vector that is a stack of the columns of x, an m by n matrix.
    
    middle = kronecker(H,diag(n))+kronecker(diag(n),H)
    vecS = solve(middle) %*% vec(StS)
    
    S = matrix(vecS,nrow=n,byrow=FALSE)
    eigens = eigen(S)$values
    
    lambda = sqrt(max(eigens[eigens>0]))
    
    
    
    ######################################################
    ######      final    ######
    ######################################################
    #mufinal = mu[1,]
    ##sigmafinal = vecS[1,]
    sigmavec[i] = lambda
    
  }
  
  sigvecmax95[k] = as.numeric(quantile(sigmavec, 0.95))
  sigvecmax70[k] = as.numeric(quantile(sigmavec, 0.7))
  #print(sigvecmax)
}

#c(as.numeric(table(new$year)))

try1 = rep(sigvecmax95[1:4],c(250,250,250,336))
#try2 = rep(sigvecmax70[1:4],c(250,250,250,336))

plot(new[,1:2],type="l",main="True data with 95% (Out-of-Sample)",ylab="RLN in 10 billion",ylim=c(-3,3))
lines(cbind(new[,1],muvec),col="red",type="l")
lines(cbind(new[248:nrow(new),1],muvec[248:nrow(new)]+qnorm(0.975)*try1),col="blue",type="l")
lines(cbind(new[248:nrow(new),1],muvec[248:nrow(new)]-qnorm(0.975)*try1),col="blue",type="l")
#lines(cbind(new[248:nrow(new),1],muvec[248:nrow(new)]+qnorm(0.975)*try2),col="green",type="l")
#lines(cbind(new[248:nrow(new),1],muvec[248:nrow(new)]-qnorm(0.975)*try2),col="green",type="l")

legend( x="bottomleft", 
        legend=c("95 perc","70 perc"),
        col=c("blue","green"), lwd=1, lty=c(1,1),  merge=FALSE)



############################################################################################################
######      in-sample and out-of-sample for only 2012 to 2015     ############################################################
############################################################################################################

try1 = rep(sigvecmax95[1:4],c(247,250,250,250))
try2 = rep(sigvecmax[1:4],c(247,250,250,250))

index = which(new$year<=15)

plot(new[index,1:2],type="l",main="True data with 95% (In-sample)",ylab="RLN in 10 billion",ylim=c(-3,3))
lines(cbind(new[index,1],muvec[index]),col="red",type="l")
lines(cbind(new[index,1],muvec[index]+qnorm(0.975)*try1),col="blue",type="l")
lines(cbind(new[index,1],muvec[index]-qnorm(0.975)*try1),col="blue",type="l")
lines(cbind(new[index,1],muvec[index]+qnorm(0.975)*try2),col="green",type="l")
lines(cbind(new[index,1],muvec[index]-qnorm(0.975)*try2),col="green",type="l")

legend( x="bottomleft", 
        legend=c("95 perc","max"),
        col=c("blue","green"), lwd=1, lty=c(1,1),  merge=FALSE)



try1 = rep(sigvecmax95[1:3],c(250,250,250))
try2 = rep(sigvecmax[1:3],c(250,250,250))

index1 = which((new$year<=15))
index = which((new$year<=15)&(new$year>=13))

plot(new[index1,1:2],type="l",main="True data with 95% (Out-of-sample)",ylab="RLN in 10 billion",ylim=c(-3,3))
lines(cbind(new[index1,1],muvec[index1]),col="red",type="l")
lines(cbind(new[index,1],muvec[index]+qnorm(0.975)*try1),col="blue",type="l")
lines(cbind(new[index,1],muvec[index]-qnorm(0.975)*try1),col="blue",type="l")
lines(cbind(new[index,1],muvec[index]+qnorm(0.975)*try2),col="green",type="l")
lines(cbind(new[index,1],muvec[index]-qnorm(0.975)*try2),col="green",type="l")

legend( x="bottomleft", 
        legend=c("95 perc","max"),
        col=c("blue","green"), lwd=1, lty=c(1,1),  merge=FALSE)


############################################################################################################
######      get the values     roberto - start here to the end #############################################
############################################################################################################
###### improve the model ###################################################################################
## https://math.stackexchange.com/questions/89030/expectation-of-the-maximum-of-gaussian-random-variables ##
## https://www.gwern.net/docs/conscientiousness/2008-nadarajah.pdf #########################################
######## in particular the function to compute the pdf of the max ##########################################
############################################################################################################

data = read.csv("cover1.csv",header=TRUE)
data = data[,c(2,3,5,6,7)]

data[,1] = as.Date(data[,1])
total_raw = data[,3] + data[,4] + data[,5]
new = data.frame(data[,1],total_raw)
colnames(new) = c("date","total_raw")

####remove "Sandy"
sandy = which(new[,2]==0)
new = new[-206,]
new[,2] = as.vector(scale(new[,2],center=TRUE,scale=FALSE))/(10^10)
meanrln = mean(new[,2])

##in-sample 95%:
(mufinal + qnorm(0.975)*sigvecmax95)*(10^10)+meanrln

##out-of-sample 95%:
(mufinal + qnorm(0.975)*sigvecmax95[1:4])*(10^10)+meanrln

####only 2012 to 2015
##in-sample 95%:
(mufinal + qnorm(0.975)*sigvecmax95[1:4])*(10^10)+meanrln
##in-sample max:
(mufinal + qnorm(0.975)*sigvecmax[1:4])*(10^10)+meanrln



############################################################################################################
######      filtered EWMA     ############################################################
############################################################################################################


new = data.frame(new)
new$year = as.numeric(format(new$date,"%y"))

  first = 500
  total = nrow(new)-first+1
  ###using first 150 as average
  sigmavecori = rep(0,total)
  sigmavecewma = rep(0,total)
  lambdaewma = 0.97
  avernum = 250
  
  for (i in 1:total){
    # print(i)
    nstart = first+(i-1)
    subcut = ts(new[i:nstart,2],frequency = 21)
    ################################################
    ######      seperate seasonal component     ######
    ################################################
    dnew = stl(subcut, "periodic")
    seasonal   <- dnew$time.series[,1]
    trend     <- dnew$time.series[,2]
    random  <- dnew$time.series[,3]
    
    dseason = subcut-seasonal
    #plot(dseason,type="l")
    
    dseason2 = ts(dseason,frequency=1)
    #plot(dseason2,type="l")
    
    ######################################################
    ######      represent into AR infinity model     ######
    ######################################################
    #armafit <- arma(dseason2, order = c(1, 1))
    armafit = arima(dseason2, order=c(1,1,3))
    
    #########Method1: directly use the estimated arma coefficients
    
    arorder = armafit$coef[1]
    maorder = armafit$coef[2:4]
    
    #########Method2: use the values from acf and pacf (not correct)
    
    #Switch the phi and the theta in the ARMAtoMA module to get the PIs of the infinite AR process
    
    n = 10
    coeff = ARMAtoMA(ar = c(maorder), ma = c(arorder), lag.max=n)
    
    
    ######################################################
    ######      estimate epsilon t     ######
    ######################################################
    
    x = new[i:nstart,2]
    
    ct = seasonal[length(seasonal)]
    pastn = as.matrix(x[(length(seasonal)-n):(length(seasonal)-1)])
    et = x[length(x)]-t(as.matrix(coeff))%*%pastn-ct
    
    
    ######################################################
    ######      construct vector and matrix     ######
    ######################################################
    
    cvec = as.matrix(c(ct,rep(0,(n-1))))
    C = rbind(coeff,cbind(diag((n-1)),0))
    evec = as.matrix(c(et,rep(0,(n-1))))
    
    
    ######################################################
    ######      estimation     ######
    ######################################################
    #H = -logm(C,method="Higham08")
    H = -logm(C,method="Eigen")
    mu = solve(diag(n)-C) %*% cvec
    StS = C %*% cov(evec %*% t(evec)) %*% t(C)
    
    ###vec function returns a column vector that is a stack of the columns of x, an m by n matrix.
    
    middle = kronecker(H,diag(n))+kronecker(diag(n),H)
    vecS = solve(middle) %*% vec(StS)
    
    S = matrix(vecS,nrow=n,byrow=FALSE)
    eigens = eigen(S)$values
    
    lambda = sqrt(max(eigens[eigens>0]))
    
    
    
    ######################################################
    ######      final    ######
    ######################################################
    #mufinal = mu[1,]
    ##sigmafinal = vecS[1,]
    sigmavecori[i] = lambda
    
  }
  
  sigmavecewma[avernum] = quantile(sigmavecori[1:avernum],0.95)
  for (i in (avernum+1):length(sigmavecewma)){
    sigmavecewma[i] = max(sigmavecewma[avernum],lambdaewma*sigmavecori[i]+(1-lambdaewma)*sigmavecewma[i-1])
  }
  
  sigmavecewma2 = rep(0,total)
  sigmavecewma2[avernum] = quantile(sigmavecori[1:avernum],0.95)
  for (i in (avernum+1):length(sigmavecewma2)){
    sigmavecewma2[i] = max(sigmavecewma2[avernum],0.93*sigmavecori[i]+(1-0.93)*sigmavecewma2[i-1])
  }
  
  
plot(new[(first+avernum-1):nrow(new),1:2],type="l",main="True data with 95% Confidence Interval (EWMA-0.97)",ylab="RLN in 10 billion",ylim=c(-3,3))
lines(cbind(new[(first+avernum-1):nrow(new),1],muvec[(first+avernum-1):nrow(new)]),col="red",type="l")
lines(cbind(new[(first+avernum-1):nrow(new),1],muvec[(first+avernum-1):nrow(new)]+qnorm(0.975)*sigmavecewma[avernum:length(sigmavecewma)]),col="blue",type="l")
lines(cbind(new[(first+avernum-1):nrow(new),1],muvec[(first+avernum-1):nrow(new)]-qnorm(0.975)*sigmavecewma[avernum:length(sigmavecewma)]),col="blue",type="l")
lines(cbind(new[(first+avernum-1):nrow(new),1],muvec[(first+avernum-1):nrow(new)]+qnorm(0.975)*as.numeric(quantile(sigmavecori, 0.95))),col="orange",type="l")
lines(cbind(new[(first+avernum-1):nrow(new),1],muvec[(first+avernum-1):nrow(new)]-qnorm(0.975)*as.numeric(quantile(sigmavecori, 0.95))),col="orange",type="l")

legend( x="bottomleft", 
        legend=c("EWMA","Asymptotic Mean","Rolling"),
        col=c("blue","red","orange"), lwd=1, lty=c(1,1,1),  merge=FALSE,cex=0.75)

plot(new[(first+avernum-1):nrow(new),1:2],type="l",main="True data with 95% Confidence Interval (EWMA-0.93)",ylab="RLN in 10 billion",ylim=c(-3,3))
lines(cbind(new[(first+avernum-1):nrow(new),1],muvec[(first+avernum-1):nrow(new)]),col="red",type="l")
lines(cbind(new[(first+avernum-1):nrow(new),1],muvec[(first+avernum-1):nrow(new)]+qnorm(0.975)*sigmavecewma2[avernum:length(sigmavecewma2)]),col="blue",type="l")
lines(cbind(new[(first+avernum-1):nrow(new),1],muvec[(first+avernum-1):nrow(new)]-qnorm(0.975)*sigmavecewma2[avernum:length(sigmavecewma2)]),col="blue",type="l")
lines(cbind(new[(first+avernum-1):nrow(new),1],muvec[(first+avernum-1):nrow(new)]+qnorm(0.975)*as.numeric(quantile(sigmavecori, 0.95))),col="orange",type="l")
lines(cbind(new[(first+avernum-1):nrow(new),1],muvec[(first+avernum-1):nrow(new)]-qnorm(0.975)*as.numeric(quantile(sigmavecori, 0.95))),col="orange",type="l")

legend( x="bottomleft", 
        legend=c("EWMA","Asymptotic Mean","Rolling"),
        col=c("blue","red","orange"), lwd=1, lty=c(1,1,1),  merge=FALSE,cex=0.75)

# ###############no 2016 considered
# 
# newcut = subset(new,year<16)
# 
# first = 200
# total = nrow(newcut)-first+1
# ###using first 150 as average
# sigmavecori = rep(0,total)
# sigmavecewma = rep(0,total)
# lambdaewma = 0.97
# avernum = 250
# 
# for (i in 1:total){
#   # print(i)
#   nstart = first+(i-1)
#   subcut = ts(newcut[i:nstart,2],frequency = 21)
#   ################################################
#   ######      seperate seasonal component     ######
#   ################################################
#   dnew = stl(subcut, "periodic")
#   seasonal   <- dnew$time.series[,1]
#   trend     <- dnew$time.series[,2]
#   random  <- dnew$time.series[,3]
#   
#   dseason = subcut-seasonal
#   #plot(dseason,type="l")
#   
#   dseason2 = ts(dseason,frequency=1)
#   #plot(dseason2,type="l")
#   
#   ######################################################
#   ######      represent into AR infinity model     ######
#   ######################################################
#   #armafit <- arma(dseason2, order = c(1, 1))
#   armafit = arima(dseason2, order=c(1,1,3))
#   
#   #########Method1: directly use the estimated arma coefficients
#   
#   arorder = armafit$coef[1]
#   maorder = armafit$coef[2:4]
#   
#   #########Method2: use the values from acf and pacf (not correct)
#   
#   #Switch the phi and the theta in the ARMAtoMA module to get the PIs of the infinite AR process
#   
#   n = 10
#   coeff = ARMAtoMA(ar = c(maorder), ma = c(arorder), lag.max=n)
#   
#   
#   ######################################################
#   ######      estimate epsilon t     ######
#   ######################################################
#   
#   x = newcut[i:nstart,2]
#   
#   ct = seasonal[length(seasonal)]
#   pastn = as.matrix(x[(length(seasonal)-n):(length(seasonal)-1)])
#   et = x[length(x)]-t(as.matrix(coeff))%*%pastn-ct
#   
#   
#   ######################################################
#   ######      construct vector and matrix     ######
#   ######################################################
#   
#   cvec = as.matrix(c(ct,rep(0,(n-1))))
#   C = rbind(coeff,cbind(diag((n-1)),0))
#   evec = as.matrix(c(et,rep(0,(n-1))))
#   
#   
#   ######################################################
#   ######      estimation     ######
#   ######################################################
#   #H = -logm(C,method="Higham08")
#   H = -logm(C,method="Eigen")
#   mu = solve(diag(n)-C) %*% cvec
#   StS = C %*% cov(evec %*% t(evec)) %*% t(C)
#   
#   ###vec function returns a column vector that is a stack of the columns of x, an m by n matrix.
#   
#   middle = kronecker(H,diag(n))+kronecker(diag(n),H)
#   vecS = solve(middle) %*% vec(StS)
#   
#   S = matrix(vecS,nrow=n,byrow=FALSE)
#   eigens = eigen(S)$values
#   
#   lambda = sqrt(max(eigens[eigens>0]))
#   
#   
#   
#   ######################################################
#   ######      final    ######
#   ######################################################
#   #mufinal = mu[1,]
#   ##sigmafinal = vecS[1,]
#   sigmavecori[i] = lambda
#   
# }
# 
# sigmavecewma[avernum] = quantile(sigmavecori[1:avernum],0.95)
# for (i in (avernum+1):length(sigmavecewma)){
#   sigmavecewma[i] = max(sigmavecewma[avernum],lambdaewma*sigmavecori[i]+(1-lambdaewma)*sigmavecewma[i-1])
# }
# 
# 
# plot(newcut[(first+avernum-1):nrow(newcut),1:2],type="l",main="True data with 95% Confidence Interval (EWMA)",ylab="RLN in 10 billion",ylim=c(-3,3))
# lines(cbind(newcut[(first+avernum-1):nrow(newcut),1],muvec[(first+avernum-1):nrow(newcut)]),col="red",type="l")
# lines(cbind(newcut[(first+avernum-1):nrow(newcut),1],muvec[(first+avernum-1):nrow(newcut)]+qnorm(0.975)*sigmavecewma[avernum:length(sigmavecewma)]),col="blue",type="l")
# lines(cbind(newcut[(first+avernum-1):nrow(newcut),1],muvec[(first+avernum-1):nrow(newcut)]-qnorm(0.975)*sigmavecewma[avernum:length(sigmavecewma)]),col="blue",type="l")
# lines(cbind(newcut[(first+avernum-1):nrow(newcut),1],sigmavecori[avernum:length(sigmavecori)]),col="green",type="l")
# 
# 
# 
# plot(new[first:nrow(new),1:2],type="l",main="True data with 95% (In-sample)",ylab="RLN in 10 billion",ylim=c(-10,10))
# lines(cbind(new[first:nrow(new),1],sigmavec),col="green",type="l")


