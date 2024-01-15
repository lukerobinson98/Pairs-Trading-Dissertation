library(tidyverse)
library(tseries)
library(quantmod)

###PAIRS SELECTION###

corp=c('TSLA','GM','HMC','TM','AMZN','AAPL','MSFT','IBM','ORCL','INTC','HPQ','ADBE','FB','TWTR','GOOGL','TGT','WMT','KO','PEP','BP','C','XOM','CVX','NG', 'AXP','GZPFY','JPM','GS','MS','WFC')

corpstock=lapply(corp,function(x){getSymbols(x,from = "2017-01-01",to="2021-01-01", periodicity="daily",auto.assign=FALSE)})

names(corp)=corpstock

AdjClosePrice=lapply(corpstock,Cl)
AdjClosePrice=do.call(merge,AdjClosePrice)

names(AdjClosePrice)=sub("\\.Close","",names(AdjClosePrice))
head(AdjClosePrice)


trainset=log(AdjClosePrice[1:504])

testset=log(AdjClosePrice[505:1007])

left=NULL
right=NULL
correlation=NULL
beta=NULL
pValue=NULL

for(i in 1:length(corp)){
  for(j in 1:length(corp)){
    if(i>j){
      left=c(left,corp[i])
      right=c(right,corp[j])
      correlation=c(correlation,cor(trainset[,corp[i]],trainset[,corp[j]]))
      
      reg=lm(trainset[,corp[i]]~trainset[,corp[j]]-1)
      beta=c(beta,as.numeric(coef(reg)[1]))
      
      spread=residuals(reg)
      
      pValue=c(pValue,adf.test(spread,alternative="stationary",k=0)$p.value)
    }
  }
}

results=data.frame(left,right,correlation,beta,pValue)
cointegratedPairs=results%>%filter(pValue<0.05)%>%arrange(pValue)
correlatedPairs=results%>%filter(correlation>0.95)%>%arrange(correlation)
bestPairs=results%>%filter(pValue<0.05,correlation>0.875)%>%arrange(pValue)
cointegratedPairs
correlatedPairs
bestPairs

