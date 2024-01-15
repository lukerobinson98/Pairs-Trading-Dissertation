install.packages("xts")
install.packages("MASS")
install.packages("MTS")
install.packages("KFAS")
install.packages("TTR")

library(quantmod)
library(xts)
library(MASS)
library(MTS)
library(KFAS)
library(TTR)

###TRADING ALGORITHM###

begin="2020-01-01"
end="2021-01-01"

stocks=corp=c('TSLA','GM','HMC','TM','AMZN','AAPL','MSFT','IBM','ORCL','INTC','HPQ','ADBE','FB','TWTR','GOOGL','TGT','WMT','KO','PEP','BP','C','XOM','CVX','NG', 'AXP','GZPFY','JPM','GS','MS','WFC')

stockPrices=xts()

for (stock_index in 1:length(stocks))
  stockPrices=cbind(stockPrices, Ad(getSymbols(stocks[stock_index],from = begin, to = end, auto.assign = FALSE)))

colnames(stockPrices)=stocks
indexClass(stockPrices)="Date"
Y=log(stockPrices)

set.seed(252)
synth_t=1000

trueMu=xts(c(rep(0.6, synth_t/4), rep(0.4, synth_t/2), rep(0.6, synth_t/4)), order.by = as.Date("2020-01-01") + 0:(synth_t-1))
colnames(trueMu)= "True Mu"

trueGamma=xts(c(rep(0.2, synth_t/2), rep(0.8, synth_t/2)), index(trueMu))
colnames(trueGamma)="True Gamma"

trueMu[]=filter(trueMu, rep(1, 50)/50, sides = 1)
trueMu[]=filter(trueMu, rep(1, 50)/50, sides = 1)
trueMu=na.locf(trueMu, fromLast = TRUE)
trueGamma[]=filter(trueGamma, rep(1, 50)/50, sides = 1)
trueGamma[]=filter(trueGamma, rep(1, 50)/50, sides = 1)
trueGamma=na.locf(trueGamma, fromLast = TRUE)
plot(cbind(trueMu, trueGamma), col = c("light blue", "orange"),legend.loc = "bottomleft", main = "True values for Mu and Gamma")

rwalkDailyVol=0.5/sqrt(252)
rwalk_sd=0.2*rwalkDailyVol
rwalk_trend=cumsum(MASS::mvrnorm(synth_t, 0, rwalk_sd^2))

residDailyVol=0.4/sqrt(252)
resid_sd=0.1*residDailyVol
resid=MTS::VARMAsim(synth_t, arlags = 1, phi = 0.9*diag(2), sigma = resid_sd^2*diag(2))$series  

synth_Y=cbind(trueMu, 0) + cbind(trueGamma, 1)*rwalk_trend + resid
colnames(synth_Y)=c("Y1 = Mu + Gamma*X + W1", "Y2 = X + W2") 
plot(synth_Y, legend.loc = "left", main = "Synthetic data")

#Least Squares Regression
LS_coef=coef(lm(synth_Y[1:(0.3*synth_t), 1] ~ synth_Y[1:(0.3*synth_t), 2]))
Mu_LS <- trueMu
Mu_LS[] <- LS_coef[1]
colnames(Mu_LS) <- "Mu from LS Regression"
Gamma_LS <- trueGamma
Gamma_LS[] <- LS_coef[2]
colnames(Gamma_LS) <- "Gamma from LS Regression"

par(mfrow = c(2, 1))
plot(cbind(trueMu, Mu_LS), legend.loc = "bottom", main = "Least Squares estimation of Mu")
plot(cbind(trueGamma, Gamma_LS), legend.loc = "topleft", main = "Least Squares estimation of Gamma")

KalFiltGamma=KalFiltMu=xts(rep(NA, synth_t), index(trueMu))
colnames(KalFiltMu)="Mu: Kalman-filtering"
colnames(KalFiltGamma)="Gamma: Kalman-filtering"
KalSmoothGamma=KalSmoothMu=xts(rep(NA, synth_t), index(trueMu))
colnames(KalSmoothMu)="Mu: Kalman-smoothing"
colnames(KalSmoothGamma)="Gamma: Kalman-smoothing"

Tt=diag(2)
Rt=diag(2)
Qt=diag(c(1e-3, 1e-2))  
Zt=array(as.vector(t(cbind(1, as.matrix(synth_Y[, 2])))), dim = c(1, 2, synth_t))
Ht=matrix(1.3e-7)

a1=matrix(c(0.6, 0.2), 2, 1)
P1=1e-12*diag(2) 
P1inf=0*diag(2)

kalmanModel=SSModel(as.matrix(synth_Y[, 1]) ~ 0 + SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, a1=a1, P1=P1, P1inf=P1inf), H=Ht)

KFSmodel=KFS(kalmanModel)

KalFiltMu[]=KFSmodel$a[-1, 1] 
KalFiltGamma[]=KFSmodel$a[-1, 2]
KalSmoothMu[]=KFSmodel$alphahat[, 1]  
KalSmoothGamma[]=KFSmodel$alphahat[, 2]

KalFiltMu[]=filter(KalFiltMu, rep(1, 5)/5, sides = 1)
KalFiltMu=na.locf(KalFiltMu, fromLast = TRUE)
KalFiltGamma[]=filter(KalFiltGamma, rep(1, 5)/5, sides = 1)
KalFiltGamma=na.locf(KalFiltGamma, fromLast = TRUE)

compute_spread <- function(Y, Gamma, Mu, name = NULL) {
  port_spread <- cbind(1, -Gamma)/cbind(1+Gamma, 1+Gamma)
  spread <- rowSums(Y * port_spread) - Mu/(1+Gamma)
  colnames(spread) <- name
  return(spread)
}

trueSpread <- compute_spread(synth_Y, trueGamma, trueMu, "True Spread")
spread_LS <- compute_spread(synth_Y, Gamma_LS, Mu_LS, "Least Squares Regression Spread")
spread_Kalman <- compute_spread(synth_Y, KalFiltGamma, KalFiltMu, "Kalman Spread")

plot(cbind(trueSpread, spread_LS, spread_Kalman), legend.loc = "bottomright", main = "Spreads")

LS_estimates=function(Y,training_prop=0.3){
  T=nrow(Y)
  T_train=round(training_prop*T)
  
  coeffs_LS=coef(lm(Y[1:T_train, 1] ~ Y[1:T_train, 2]))
  Mu=xts(rep(coeffs_LS[1], T), index(Y))
  colnames(Mu)="Mu: Least Squares"
  Gamma=xts(rep(coeffs_LS[2], T), index(Y))
  colnames(Gamma)="Gamma: Least Squares"
  return(list(Mu = Mu, Gamma = Gamma))
}

Kal_estimates=function(Y){
  T=nrow(Y)
  
  KalFiltGamma=KalFiltMu=xts(rep(NA, T), index(Y))
  colnames(KalFiltMu)="Mu: KFS"
  colnames(KalFiltGamma)="Gamma: KFS"
  
  Tt=diag(2)
  Rt=diag(2)
  Qt=1e-5*diag(2)  
  Zt=array(as.vector(t(cbind(1, as.matrix(Y[, 2])))), dim = c(1, 2, T))
  Ht=matrix(1e-3)
  
  init=LS_estimates(Y)
  
  a1 =matrix(c(init$mu[1], init$Gamma[1]), 2, 1)
  P1=1e-5*diag(2)
  P1inf=0*diag(2)
  
  kalmanModel=SSModel(as.matrix(Y[, 1]) ~ 0 + SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, a1=a1, P1=P1, P1inf=P1inf), H=Ht)
  
  KFSmodel=KFS(kalmanModel)
  KalFiltMu[]=KFSmodel$a[-1, 1]
  KalFiltGamma[]=KFSmodel$a[-1, 2]
  
  L=30
  KalFiltMu[]=filter(KalFiltMu, rep(1, L)/L, sides = 1)
  KalFiltMu=na.locf(KalFiltMu, fromLast = TRUE)
  KalFiltGamma[]=filter(KalFiltGamma, rep(1, L)/L, sides = 1)
  KalFiltGamma=na.locf(KalFiltGamma, fromLast = TRUE)
  return(list(Mu = KalFiltMu, Gamma = KalFiltGamma))
}  

Y_=Y["2020-01-01::2021-01-01", c("GZPFY", "BP")]
if(anyNA(Y_)) 
  Y_=na.approx(Y_)
plot(Y_, legend.loc = "bottomleft", main = "Log-prices")

Kal=Kal_estimates(Y_)
LS=LS_estimates(Y_)

spread_LS <- compute_spread(Y_, LS$Gamma, LS$Mu, "Least Squares Spread")
spread_Kalman <- compute_spread(Y_, Kal$Gamma, Kal$Mu, "Kalman Spread")

plot(cbind(spread_LS, spread_Kalman), legend.loc = "topleft", main = "Spreads")

generate_Z_score_EMA <- function(spread, n = 120) {
  
  spreadMean <- EMA(spread, n)
  spreadMean <- na.locf(spreadMean, fromLast = TRUE)
  spreadDemeaned <- spread - spreadMean
  
  spreadVar <- EMA(spreadDemeaned^2, n)
  spreadVar <- na.locf(spreadVar, fromLast = TRUE)
  
  zScore <- spreadDemeaned/sqrt(spreadVar)
  return(zScore)
}

generateSignal <- function(Z_score, threshold_long, threshold_short) {
  signal <- Z_score
  colnames(signal) <- "signal"
  signal[] <- NA
  
  signal[1] <- 0
  if (Z_score[1] <= threshold_long[1]) {
    signal[1] <- 1
  } else if (Z_score[1] >= threshold_short[1])
    signal[1] <- -1
  
  for (t in 2:nrow(Z_score)) {
    if (signal[t-1] == 0) {  #no position
      if (Z_score[t] <= threshold_long[t]) {
        signal[t] <- 1
      } else if(Z_score[t] >= threshold_short[t]) {
        signal[t] <- -1
      } else signal[t] <- 0
    } else if (signal[t-1] == 1) {  #long position
      if (Z_score[t] >= 0) signal[t] <- 0
      else signal[t] <- signal[t-1]
    } else {  #short position
      if (Z_score[t] <= 0) signal[t] <- 0
      else signal[t] <- signal[t-1]
    }
  }
  return(signal)
}

trade <- function(Y, Gamma, Mu, name = NULL, threshold = 0.7, plot = FALSE) {
  
  port_spread <- cbind(1, -Gamma)/cbind(1+Gamma, 1+Gamma)
  spread <- rowSums(Y * port_spread) - Mu/(1+Gamma)
  
  Z_score <- generate_Z_score_EMA(spread)
  threshold_long <- threshold_short <- Z_score
  threshold_short[] <- threshold
  threshold_long[] <- -threshold
  
  signal <- generateSignal(Z_score, threshold_long, threshold_short)
  
  port <- port_spread * lag(cbind(signal, signal))
  
  X <- diff(Y)
  port_ret <- xts(rowSums(X * port), index(X))
  port_ret[is.na(port_ret)] <- 0
  colnames(port_ret) <- name
  
  if (plot) {
    tmp <- cbind(Z_score, signal)
    colnames(tmp) <- c("Z-score", "Signal")
    par(mfrow = c(2, 1))
    { plot(tmp, legend.loc = "bottomright",main = paste("Z-score and trading on spread based on", name))
      lines(threshold_short, lty = 2)
      print(lines(threshold_long, lty = 2)) }
    print(plot(cumprod(1 + port_ret), main = paste("Cumulative Profit & Loss when trading with GZPFY and BP for spread based on", name)))
  }
  
  return(port_ret)
}

ret_Kal=trade(Y_["::2021-01-01"], Kal$Gamma["::2021-01-01"], Kal$Mu["::2021-01-01"], 
              "Kalman", plot = TRUE)

