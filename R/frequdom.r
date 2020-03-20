##
##  frequdm.r                                                                  R. Schlittgen    07.01.2020
##
##     functions  for analysis in the frequency domain 

 
#' \code{armathspec} determines the theoretical spectrum of an arma process 
#'  
#' @param  a   ar-coefficients
#' @param  b   ma-coefficients
#' @param  nf   scalar, the number of equally spaced frequencies
#' @param  s    variance of error process
#' @param  pl   logical,  if TRUE, the spectrum is plotted, FALSE for no plot 
#' @return  out   (nf+1,2) matrix, the frequencies and the spectrum  
#'
#' @examples 
#' out <-armathspec(c(0.3,-0.5),c(-0.8,0.7),50,s=1,pl=FALSE) 
#'
#' @export 
 
armathspec<-function(a,b,nf,s=1,pl = FALSE){  
p <- length(a)
q <- length(b)
f <- seq(0,0.5,length.out = nf)
spec <- rep(0,nf)  
for(j in 1:nf){ 
               spec[j] <- s*(abs(1-sum(b*exp(1i*2*pi*f[j]*c(1:q))))^2)/
                            (abs(1-sum(a*exp(1i*2*pi*f[j]*c(1:p))))^2)
              }
if(pl == TRUE){ plot(f,spec,type="l") }
return( cbind (f,spec) )
}

 
#' \code{periodogram} determines the periodogram of a time series
#'  
#' @param  y   (n,1) vector, the time series or  an  acf at lags 0,1,...,n-1
#' @param  ACF  logical, FALSE, if y is ts, TRUE, if y is acf
#' @param  nf   scalar, the number of equally spaced frequencies
#' @param  type  c("cov","cor"), area under spectrum, can be variance or normed to 1. 
#' @return  out  (floor(nf/2)+1,2) matrix, the frequencies and the periodogram
#'
#' @examples 
#' data(WHORMONE)
#' out <-periodogram(WHORMONE,50,ACF=FALSE,type="cov")
#'
#' @export 
 

periodogram<-function(y,nf,ACF=FALSE,type="cov"){ 
   n <- length(y) 
if (ACF==FALSE){  
   y1 <- y - mean(y)
   if(type == "cor"){ y1 <- y1/sd(y1) }    
   n <- length(y1)
   f <- c(0:nf)/(2*nf)
   per <- rep(0,nf+1) 
   for (j in c(1:(nf+1))){
       per[j] <- sum(exp((0+1i)*2*pi*f[j]*c(1:n))*y1 ) 
       }
   per<- (abs(per)^2)/n   
   }  
if (ACF==TRUE){    
   if(type == "cor"){ y <- y/y[1] } 
   f <- c(0:nf)/(2*nf)
   per <- rep(0,nf+1) 
   for (j in c(1:(nf+1))){
        per[j] <- y[1]+2*sum(cos(2*pi*f[j]*c(1:(n-1)))*y[2:n]) 
       }
    }   
cbind(f,per)  
}
 

#' \code{perwinpa} Parzen's window for direct spectral estimation
#'  
#' @param  e   equal bandwidth  (at most n frequencies are used for averaging)
#' @param  n   length of time series 
#' @return w   weights (symmetric)
#'
#' @examples 
#' data(WHORMONE)
#' w <- perwinpa(0.1,length(WHORMONE)) 
#'
#' @export  

perwinpa<-function(e,n){ 
  if((e >= 0.5)|(e < 0)){ stop("e must be 0<=e<0.5")  }
  q <- 1 
  w <- 1  
  while ((1/(n*sum(w^2)) < e) & (q < n/2)){ 
    q <- q+1 
    w <- c(1:q)/q 
    w <- c((1-6*(w[1:q/2]/q)^2 + 6*(w[1:q/2]/q)^3) ,(2*(1-w[(1+q/2):q-1])^3)) 
    w <- c(rev(w),1,w)
    w <- w/sum(w) 
  }
return(w) 
}

#' \code{perwinba} Bartlett-Priestley  window for direct spectral estimation
#'  
#' @param  e   equal bandwidth  (at most n frequencies are used for averaging)
#' @param  n   length of time series 
#' @return  w   weights (symmetric)
#'
#' @examples 
#' data(WHORMONE)
#' w <- perwinba(0.1,length(WHORMONE)) 
#'
#' @export 
  
perwinba<-function(e,n){ 
  if((e >= 0.5)|(e < 0)){ stop("e must be 0 <=e <0.5") }
  q <- 1 
  w <- 1 
  while ((1/(n*sum(w^2)) < e)&(q < n/2)){
    q <- q+1;
    w <- (3*q/(4*q^2-1))*(1-(c(1:(q-1))/q)^2) 
    w <- c(rev(w),(3*q/(4*q^2-1)),w)
    w <- w/sum(w) 
  }
return(w)
} 

#' \code{specest}  direct spectral estimation of series y
#'                              using periodogram window win
#'  
#' @param  y   (n,1) vector, the ts
#' @param  nf  number of equally spaced frequencies
#' @param  e   equal bandwidth
#' @param  win string, name of periodogram window
#'                             (possible: "perwinba", "perwinpa")
#' @param conf  scalar, the level for confidence intervals
#' @param type  c("cov","cor"), area under spectrum is variance or is normed to 1.
#'
#' @return  est  (nf+1,2)- or (nf+1,4)-matrix:   
#' \item{column 1:}{frequencies 0, 1/n, 2/n, ..., m/n }
#' \item{column 2:}{the estimated spectrum  }
#' \item{column 3+4:}{the confidence bounds   }
#' @examples 
#' data(WHORMONE)
#' est <- specest(WHORMONE,50,0.05,win = c("perwinba","perwinpa"),conf=0,type="cov") 
#'
#' @export  

specest <- function(y,nf,e,win = c("perwinba","perwinpa"), conf=0,type="cov"){ 
  y <-y-mean(y) 
  n <- length(y)  
  if(win[1]=="perwinpa"){  w <- perwinpa(e,n)  }
  else{ w <- perwinba(e,n)  }
  m <- (length(w)-1)/2 
  if(type=="cov"){ fp <- periodogram(y,nf,type="cov") }
  if(type=="cor"){ fp <- periodogram(y,nf,type="cor") }
  p <- fp[,2] 
  p1 <- c(rev(p[2:(m+1)]),p,rev(p[(nf+1-m):(nf)])) 
  for(j in 0:nf){ 
   p[j+1]<- sum(w*p1[j+c(1:(2*m+1))])   
   }   
  if ((conf > 0) & (conf < 1)){
   nu <- 2/sum(w^2) 
   u <- p*nu/qchisq(1-(1-conf)/2,nu) 
   o <- p*nu/qchisq((1-conf)/2,nu) 
   p <- cbind(p,u,o)
  }
 return(  cbind(fp[,1],p) )
}  


#' \code{specplot}    plot of spectral estimate 
#'  
#' @param   s  (n,2) or (n,4) matrix, output of specest 
#' @param   Log, logical, if TRUE, the logs of the spectral estimates are shown
#'
#' @examples 
#' \donttest{
#' data(WHORMONE)
#' est <- specest(WHORMONE,50,0.05,win = c("perwinba","perwinpa"),conf=0,type="cov") 
#' specplot(est,Log=FALSE) }
#' @export 
 
specplot <- function(s,Log=FALSE){
   f <- s[,1]
   sp <-s[,2]
   u <- sp
   o <- sp
   if(ncol(s)==4){
      u <- s[,3]
      o <- s[,4]
      }
   if(Log==FALSE){ 
   plot(f,sp,type="l",ylim=c(0,max(sp,u,o)),xlab="Frequenz",ylab="Spektrum")
   if(ncol(s)==4){
      lines(f,u,lty=2)
      lines(f,o,lty=2)
      } 
   } 
   if(Log==TRUE){ 
   plot(f,log(sp),type="l",ylim=c(0,max(log(sp),log(u),log(o))),xlab="Frequenz",ylab="Logspektrum")
   if(ncol(s)==4){
      lines(f,log(u),lty=2)
      lines(f,log(o),lty=2)
      } 
   }     
   }      

#' \code{wntest}   graphical test for white noise for a time series  or a series of regression residuals 
#'  
#' @param  e  vector, the time series (k = 0) or residuals (k > 0)
#' @param  a  scalar, level of significance
#' @param  k  scalar >= 0, number of regressors used to compute e as residuals

#' @return tp vector,  value of test statistic and  p-value
#'
#' @examples 
#' \donttest{
#' data(WHORMONE)
#' out <- wntest(WHORMONE,0.05,0) }
#' @export 
 
wntest <- function(e,a,k=0){ 
    if (a<0 | a>1){ stop("error probability must be inside (0,1)")  }
    n<-length(e);   m<-trunc(n/2);    s<-rep(0,m); p<-NA
    cf <- acf(e,lag.max=n-1,type="covariance",plot=F)
    c0 <- cf$acf[1]; cf <- cf$acf[2:n];
    se<-c(1:(n-1))/n 
    s[1] <- c0 + 2*sum(cf*cos(2*pi*se)) 
    i <- 1 
    while (i < m){
       i <- i+1;
       s[i] = s[i-1] + c0 + 2*sum(cf*cos(2*pi*i*se));
    }
    s <- c(0,s/s[m])
    if (k == 0){
       nen <- sqrt(m-1) + 0.2 + 0.68/sqrt(m-1) 
       c <- sqrt(-log(a/2)/2)/nen - 0.4/(m-1) 
       x <- seq(0,1,1/m)
       t <- max(abs(s-x)) 
       p <- min(c(1,(2*exp(-2*(t+0.4/(m-1))*(t+0.4/(m-1))*nen*nen)))) 
       k1 <- x-c;#  k1=miss(k1.*(k1.>=0),0);
       k2 <- x+c; # k2=miss(k2.*(k2.<=1),0);
       plot(x,x,type="l",xlab=paste("test statistic: ",round(t,3),",  p-value : ",round(p,3)),ylab="S")    
   #    title("bartletts kolmogorov test for white noise" ); 
       lines(x[-1],s[-1],lwd=2 )
       lines(x,k1,lwd=2,lty=2 )
       lines(x,k2,lwd=2,lty=2 )  
       } 
     if (k > 0){
       c <- sqrt(-log(a/2)/2)/(sqrt(m-1) + 0.2 + 0.68/sqrt(m-1)) - 0.4/(m-1) 
       x=seq(0,1,1/m)
       t <- max(abs(s-x)) 
       mm <- (n-k)/2 
       cc <- c-(k-1)/(2*mm)  
       plot(x,x,type="l",xlab=paste("test statistic: ",round(t,3)),ylab="S")    
   #    title(paste("bartletts kolmogorov test for white noise with ",k," regressors") ) 
       lines(x[-1],s[-1],lwd=2 )
       lines(c(0,mm/m*(1-c)),  c(c,1) ,lty=2,lwd=2 )
       lines(c(0,mm/m*(1-cc) ),c(cc,1) ,lty=2,lwd=2)
       lines(c(c,1) , c(0,mm/m*(1-c))  ,lty=2 ,lwd=2 )
       lines(c(cc,1), c(0,mm/m*(1-cc) ) ,lty=2,lwd=2 )
    }
 return( c(t,p) )
 }


#' \code{taper} taper modification of a time series 
#'  
#' @param  y    the time series
#' @param  part scalar, 0 <= part <= 0.5, part of modification (at each end of y)   
#' @return tp   tapered time series
#'
#' @examples 
#' data(WHORMONE)
#' out <-taper(WHORMONE,0.3)
#' \donttest{ 
#' plot(WHORMON) 
#' lines(out,col="red") }
#' @export 
 
taper <- function(y,part){  
    if((part <  0)|(part > 0.5)){ stop("part must be 0 <= part <= 0.5") } 
    n <- length(y) 
    if(part==0){ tap <- rep(1,n) }
    if(part>0){ 
    tap <-  0.5*(1-cos(c(1:(part*n))*pi/(part*n))) 
    tap <- c(tap,rep(1,n-2*length(tap)),rev(tap))
    }
    return(mean(y)+(y-mean(y))*tap) 
    } 


#' \code{lagwinba} Bartlett's Lag-window for indirect spectrum estimation
#'  
#' @param  NL   number of lags used for estimation 
#' @return win vector, one-sided weights
#'
#' @examples 
#' win <-lagwinba(5) 
#' @export 
 
lagwinba <- function(NL){  
  win <- (1-c(0:NL)/(NL+1)) 
  return(win)
 } 

#' \code{lagwinpa} Parzen's Lag-window for indirect spectrum estimation
#'  
#' @param  NL   number of lags used for estimation 
#' @return win vector, one-sided weights
#'
#' @examples 
#' win <- lagwinpa(5)  
#' @export 
 
lagwinpa <- function(NL){  
  k <- c(1:T) 
  win  <- c(1,(1 - 6*(k[1:(NL/2)]/NL)^2 + 6*(k[1:(NL/2)]/NL)^3),(2*(1-k[floor(1+NL/2):(NL-1)]/NL)^3))  
  return(win)
 }  

#' \code{lagwintu} Tukey's Lag-window for indirect spectrum estimation
#'  
#' @param  NL   number of lags used for estimation 
#' @return win vector, one-sided weights
#'
#' @examples 
#' win <- lagwintu(5)   
#' @export 

lagwintu <- function(NL){  win <- 0.54 + 0.46*cos(pi*c(0:NL)/(NL+1)) ;   return(win)  }
 

#' \code{dynspecest}  performs a dynamic  spectrum estimation  
#'
#' @param y     time series or vector
#' @param nseg  number of segments for which the spectrum is estimated
#' @param nf    number of equally spaced frequencies
#' @param e     equal bandwidth
#' @param theta azimuthal viewing direction, see R function persp 
#' @param phi   colatitude viewing direction, see R function persp
#' @param d     a value to vary the strength of the perspective transformation, see R function persp
#' @param Plot  logical, schould a plot be generated?    
#'
#' @return out list with components
#' \item{f}{frequencies, vector of length nf }
#' \item{t}{time, vector of length nseg }
#' \item{spec}{the spectral estimates, (nf,nt)-matrix    }
#' @examples
#' data(IBM) 
#' y <- diff(log(IBM))
#' out <- dynspecest(y,60,50,0.2,theta=0,phi=15,d=1,Plot=FALSE)
#' @export  

dynspecest <- function( y,nseg,nf,e,theta = 0, phi = 15,d,Plot=FALSE ){ 
   n <- length(y)  
   spec <- NULL
   t <- NULL
   lseg <- max(n/nseg  ,25) 
   segstart  <- seq(from=1,to=n-lseg,by=trunc(n/nseg))         
   I <- length(segstart)
   i<-1
   while( (segstart[i] + lseg) < n ){
     out <- specest(y[segstart[i]:(segstart[i]+lseg)],50,0.2)
     f <- out[,1]
     t <- cbind(t,segstart[i])
     spec <- cbind(spec,out[,2]) 
     i <- i+1
     }
     out<- list(f=f,t=t,spec=spec)
     if(Plot == TRUE){
         persp(f,t,spec,ticktype="detailed",xlab="Frequency",ylab="Time",
               zlab="Spectrum",theta=theta,d=d) 
       }
      return(out)
      } 
#' \code{periodotest} computes the p-value of the test for a hidden periodicity
#'
#' @param y vector, the time series
#'
#' @return pval the p-value of the test
#' @examples   
#' data(PIGPRICE)
#' y <- PIGPRICE
#' out <- stl(y,s.window=6)  
#' e <- out$time.series[,3]
#' out <- periodotest(e)

periodotest <- function(y){
          n <- length(y)    
          if(n < 100){ warning("pvalue is appropriate only for n >= 100.") }
          p <- periodogram(y - mean(y), trunc(n/2))  # periodogram at Fourier frequencies 
          m <- mean(p[,2])   
          pval <- 1-(1-exp(-max(p[,2])/m))^trunc((n-1)/2)
          return(pval)
          }
