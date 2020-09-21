
##
##  smaniplot.r                                                                  R. Schlittgen    07.01.2020
##
##  functions for simple time series manipulation and plots

 
#' \code{BoxCox} determines the power of a Box-Cox transformation to stabilize the variance of a time series
#   additionaly, a diagram can be plotted
#'  
#' @param  y    the series, a vector or a time series  
#' @param  seg  scalar, number of segments 
#' @param  Plot logical, should a plot be produced?
#' @return l    scalar, the power of the Box-Cox transformation
#'
#' @examples 
#' data(INORDER)
#' lambda <-BoxCox(INORDER,6,Plot=FALSE)
#'
#' @export 


BoxCox<-function(y,seg,Plot=FALSE){
 y1 <- y
 n<-length(y1)
 breit <- ceiling(n/seg)
 if ( min(y1)<=0 ){ y1 <- y1-min(y1)+1 }
 y1 <- c(y1,rep(NA,breit*seg-n)) 
 ta <- matrix(y1,breit,seg)
 m <- log(apply(ta,2,mean,na.rm=TRUE))
 s <-log(apply(ta,2,sd,na.rm=TRUE))
 reg<-lm(s~m)
 beta<-reg$coefficients  
 lambda <- 1-as.numeric(beta[2])
 if(Plot==TRUE){
   plot(m,s,type="p",col="blue",cex=1.5,pch=16,xlab=expression(paste("ln(",bar(x)[i],")")),ylab=expression(paste("ln(",s[i],")")))
   abline(beta, col="red")
   lambda<-round(lambda*1000)/1000
   text(median(m),median(s),labels=expression(paste(hat(lambda)," = ") )) 
   text(median(m)+0.07*(max(m)-min(m)),median(s)-0.01*(max(s)-min(s)),labels=lambda)
 }
 return(lambda)
} 


#' \code{statcheck} determines the means, standard deviations and acf's of segmets of a time series 
#'                        and plots the acf's for the segments.
#'  
#' @param  y  the series, a vector or a time series  
#' @param  d  scalar, number of segments 
#' @return out list with components:
#' \item{ms}{matrix with means and standard deviations of the segments  }
#' \item{ac}{matrix with acf's, the first column: acf of the series, the others: acf's  of the segments  }
#' @examples 
#' \donttest{
#' data(COFFEE)  
#' out <- statcheck(COFFEE,4)   }
#' @export 
  
statcheck<-function(y,d){  
  n <- length(y) 
  k <- floor(n/d) 
  x <- y[c(1:k)] 
  a <- acf(y[c(1:k)],type = "correlation",plot = FALSE)
  m <- length(a$acf)
  ac <- a$acf[1:m]
  plot(c(0:(m-1)),a$acf[1:m],ylim=c(-1,1),type="l",xlab="Lag",ylab="ACF") 
  lines(c(0,m-1),c(0,0))
  for (i in 2:d){   
  x<-   y[(i-1)*k+c(1:k)]  
  a <- acf(y[(i-1)*k+c(1:k)], type = "correlation",plot = FALSE)
  lines(c(0:(m-1)),a$acf[1:m]) 
  ac <- cbind(ac,a$acf[1:m])
  }  
  m<-apply(matrix(y[1:(k*d)],k,d),2,mean)
  s<-apply(matrix(y[1:(k*d)],k,d),2,sd)
  ms <- cbind(m,s)  
  colnames(ms)<-c("segment-mean","segment-sd")
  out <- list(ms=ms, ac=ac) 
  return(out)
  } 


#' \code{vartable} determines table of variate differences 
#'  
#' @param  y      the series, a vector or a time series  
#' @param  season scalar, period of seasonal component  
#' @return d      matrix with ratios of variances for differend numbers of simple and seasonal differencing 
#'
#' @examples 
#' data(GDP)
#' out <- vartable(GDP,4) 
#'
#' @export  

vartable<-function(y,season){ 
 y1<-y
 ta<-matrix(0,4,4)
 
 ta[1,1]<-var(y1)  
 for (j in 2:4){ ta[1,j]<-var(diff(y1,lag=season,differences = j-1)) }
 for (i in 2:4)
   {
    y1<-diff(y,lag=1,differences = i-1)
    ta[i,1]<-var(y1)
    for (j in 2:4){ ta[i,j]<-var(diff(y1,lag=season,differences = j-1)) }
   }    
  ta<-ta/var(y) 
  Saisdiff0<-ta[,1] 
  Saisdiff1<-ta[,2] 
  Saisdiff2<-ta[,3] 
  Saisdiff3<-ta[,4]
  Einfachdiff<-c(0:3) 
 d<-cbind(Einfachdiff,Saisdiff0,Saisdiff1,Saisdiff2,Saisdiff3)
 return(d)
}  
 

#' \code{acfpacf} produces a plot of the acf and the pacf of a time series 
#'  
#' @param x   the series, a vector or a time series  
#' @param lag scalar, maximal lag to be plotted  
#' @param HV character, controls division of graphic window: "H" horizontal, "V" vertical, default is "H" 
#'
#' @examples 
#' \donttest{
#'   data(LYNX)
#'   acfpacf(log(LYNX),15,HV="H") }
#' @export  
 
acfpacf<-function(x,lag,HV = "H"){
 x <- as.vector(x)
 m <- mean(x,na.rm = TRUE)
 n <- length(na.omit(x))
 y <- x - m
 o <- c(1:length(y))
 y[o[is.na(y)]] <- 0
 a <- acf(y,lag,plot=FALSE) 
 a1 <- a$acf[-1]
 bart <- cumsum(c(1,2*a1^2))  
 oldpar <- par(no.readonly = TRUE) 
 on.exit(par(oldpar)) 
 if(HV == "V"){  par(mfrow=c(1,2),mex=0.8) }
 else{  par(mfrow=c(2,1),mex=0.8) }
 
 plot(c(0:lag),a$acf,type="h",lwd=2,ylim=c(-1,1),xlab="Lag",ylab="")   # ylim=c(min(c(-2.1*sqrt(bart/n),a1)),1),xlab="Lag",ylab="")
 title("ACF with Bartlett-95%-bounds")
 lines(c(0,lag),rep(0,2))
 lines(c(1:lag),2*sqrt(bart[-length(bart)]/n))
 lines(c(1:lag),-2*sqrt(bart[-length(bart)]/n))
 p <- pacf(y,lag,plot=FALSE) 
 p <- p$acf
 plot(c(0:lag),c(0,p),type="h",lwd=2,ylim=c(-1,1) ,xlab="Lag",ylab="") #   ylim=c(min(c(p,-2.1*sqrt(1/n))),max(p,2.1*sqrt(1/n))),xlab="Lag",ylab="")
 title("PACF with 95%-bounds")
 lines(c(0,lag),rep(0,2)) 
 lines(c(1,lag),rep(2*sqrt(1/n),2))
 lines(c(1,lag),rep(-2*sqrt(1/n),2)) 
 par(mfrow=c(1,1))
 }


#' \code{LjungBoxPierceTest} determines the test statistic and p values for several lags  for a residual series
#'  
#' @param  y      the series of residuals, a vector or a time series  
#' @param  n.par  number of parameters which had been estimated 
#' @param  maxlag maximal lag up to which the test statistic is computed, default is maxlag = 48
#' @return BT     matrix with columns:  
#'                  lags,   degrees of freedom,    test statistic,   p-value  
#' @examples 
#' data(COFFEE)
#' out <- arima(COFFEE,order=c(1,0,0))
#' BT <- LjungBoxPierceTest(out$residuals,1,20)
#'
#' @export  

LjungBoxPierceTest <- function(y,n.par=0,maxlag=48){
la <- seq(6,maxlag,6) 
BT<-matrix(NA,length(la),4) 
for (i in c(1:length(la))){
if(la[i]>n.par){
   bt <- Box.test(y,lag=la[i], type ="Ljung-Box",fitdf = n.par)  
   BT[i,]<-round(c(la[i],bt$parameter,bt$statistic,bt$p.value),3)
   }
 }
colnames(BT)<-c("lags","df","statistic","p-value")
return(BT)  
}     

#' \code{symplot} produces a symmetry plot 
#'  
#' @param  y   the series, a vector or a time series    
#'
#' @examples 
#' \donttest{
#' data(LYNX)
#' symplot(LYNX)   }
#' @export 
 
symplot<-function(y){
y1<-sort(y)
n<-length(y)
m<-median(y)
unten<-rev(m-y1[1:floor((n+1)/2)])
oben<-y1[ceiling((n+1)/2):n]-m
plot(unten,oben,type="p",pch=16,col="blue")
abline(0,1,col="red",lwd=2)
}
   
 
#' \code{bandfilt} does a bandpass filtering of a time series
#'  
#' @param  y  the series, a vector or a time series  
#' @param  q  scalar, half of length of symmetric weights
#' @param  pl scalar, lower periodicity ( >= 2 )
#' @param  pu scalar, upper periodicity ( > pl )
#' @return yf (n,1) vector, the centered filtered time series with NA's at beginning and ending
#'
#' @examples 
#' data(GDP)
#' yf <- bandfilt(GDP,5,2,6)
#' \donttest{plot(GDP); lines(yf+mean(GDP),col="red")   }
#' @export   
 
bandfilt <- function(y,q,pl,pu){
 lc <- (1/pl+1/pu)/2     # the frequency at which the filter is centered
 l0 <- (1/pl-1/pu)/2     # half the desired broadness of the filter (in frequency terms)
 s <- c(1:q) 
 f <- sin(2*pi*l0*s)/(pi*s) 
 f <- f*sin(2*pi*s/(2*q+1))/(2*pi*s/(2*q+1))
 f <- c(rev(f),(2*l0),f) 
 f <- f*2 # /sum(abs(f))
 f <- f*cos(2*pi*lc*c(-q:q))
 f <- filter(y, f, method = "convolution", sides = 2, circular = FALSE)
 f <- f-mean(f,na.rm = TRUE)
 return(f) 
}
   