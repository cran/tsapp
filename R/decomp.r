
##
##  decomp.r                                                                  R. Schlittgen    07.01.2020
##
##  functions for decomposition of time series 


#' \code{moveav} smoothes a time series by moving averages
#'  
#' @param  y the series, a vector or a time series  
#' @param  q scalar,  span of moving average
#' @return g vector,  smooth component
#'
#' @examples 
#' data(GDP)
#' g <- moveav(GDP,12) 
#' \donttest{ plot(GDP) ; lines(g,col="red") } 
#' @export 
 
 
moveav <- function(y,q){
  if ( q<1 ) stop("Length of moving average must be >=1")
  if (!(q%%2)){ fi<- c(0.5,rep(1,q-1),0.5)/q }
  if (q%%2) { fi <- rep(1,q)/q }
  out <- filter(y,fi, method="convolution",sides=2)
  return(out)
  }


#' \code{movemed} smoothes a time series by moving medians
#'  
#' @param  y  the series, a vector or a time series  
#' @param  q  scalar,  span of moving median
#' @return g  vector,  smooth component
#'
#' @examples 
#' data(BIP)
#' g <- movemed(GDP,12) 
#' \donttest{ plot(GDP) ; t <- seq(from = 1970, to = 2009.5,by=0.25) ; lines(t,g,col="red") } 
#' @export 
  

movemed <- function(y,q){ 
  n <- length(y) 
  g <- rep(NA,n)  
  if (1-q%%2){ 
      g1  <- apply(tsmat(y[-1],q),1,median)
      g2  <- apply(tsmat(y[-n],q),1,median)
      g <- (g1+g2)/2
   }else{
      g <- apply(tsmat(y,q),1,median)
   }
   g <- c(rep(NA,q/2),g,rep(NA,q/2))
   g   
} 


#' \code{smoothls} smoothes a time series by Whittaker graduation.  
#' The function depends on  the package Matrix. 
#'  
#' @param  y     the series, a vector or a time series  
#' @param  beta  smoothing parameter  >=0 (the larger beta is, the smoother will g be)
#' @return g     vector,  smooth component
#'
#' @examples  
#' data(GDP)
#' g <- smoothls(GDP,12)
#' \donttest{
#'  plot(GDP)   
#'  t <- seq(from = tsp(GDP)[1], to = tsp(GDP)[2],by=1/tsp(GDP)[3]) ; lines(t,g,col="red") } 
#' @export   

smoothls <- function(y,beta=0){  
   x <- as.vector(y) 
   n <- length(y)  
   if(beta==0){
     g <- movemed(x,12)
     beta <- var(na.omit(x-g))
     }
   IP <- Matrix::bandSparse(2*n-2, n, 0, diagonals=list(rep(1,n)))  
   P <- Matrix::bandSparse(n-2, n, c(0,1,2), 
                                        diagonals=list(rep(beta,n-2),rep(-2*beta,n-2),rep(beta,n-2)))  
   IP[((n+1):(2*n-2)),1:n] <- P   
   IP[((n+1):(2*n-2)),1:n] <- P   
   sx  <- Matrix::sparseMatrix(c(1:n),rep(1,n),x=x,dims=c(2*n-2,1)) 
   g <- Matrix::solve(Matrix::crossprod(IP),Matrix::crossprod(IP,sx))   
  return(g)
}


#' \code{smoothrb} smoothes a time series robustly by using Huber's psi-function. 
#' The initialisation uses a moving median. 
#'  
#' @param  y       the series, a vector or a time series  
#' @param  beta smoothing parameter (The larger beta is, the smoother will the smooth component g be.)
#' @param  q      length of running median which is used to get initial values
#' @return g       vector,  the smooth component
#'
#' @examples 
#' data(GDP)
#' g  <- smoothrb(GDP,8,q=8)
#' \donttest{ 
#'  plot(GDP) ; t <- seq(from = 1970, to = 2009.5,by=0.25) ; lines(t,g,col="red") } 
#' @export  
 
smoothrb <- function(y,beta=0,q=NA){
   x <- y
   n <- length(x)  
   g <- movemed(x,q)
   g[c(1:(1+q/2))] <- median(x[1:q])
   g[c(ceiling(n-q/2):n)] <- median(x[(n+1-q):n])  
   if(beta == 0){
     beta <- (1.4*mad(x-g))^2
     } 
   sigm <- sd(x-g)   
   IP <- Matrix::bandSparse(2*n-2, n, 0, diagonals=list(rep(1,n)))  
   P <- Matrix::bandSparse(n-2, n, c(0,1,2), diagonals=list(rep(1,n-2),rep(-2,n-2),rep(1,n-2)))  
   IP[((n+1):(2*n-2)),1:n] <- beta*P  
   g0 <- rep(0,n)    
   while(sqrt(mean(g-g0)^2)>(0.0001*sigm)){ 
     g0 <- g
     b <- c(sigm*psihuber((x-g0)/sigm),rep(0,n-2)) - c(as.vector(beta^2*Matrix::crossprod(P)%*%matrix(g,n,1)),rep(0,n-2)) 
     g <- Matrix::solve(Matrix::crossprod(IP),Matrix::crossprod(IP,b)) 
     g <- as.vector(g) 
     g <- g0 + g
     }
   return(as.vector(g))
 }


#' \code{splinedecomp} decomposes a time series into trend, season and irregular component
#' by spline approach.
#'  
#' @param  x  the series, a vector or a time series  
#' @param  d  seasonal period
#' @param  alpha smoothing parameter for trend component   (The larger alpha is, the smoother will the smooth component g be.)
#' @param  beta  smoothing parameter for  seasonal component
#' @param   Plot  logical, should a plot be produced?
#' @return    out (n,3) matrix:
#' \item{1. column}{smooth component }
#' \item{2. column}{seasonal component }
#' \item{3. column}{irregular component    }
#' @examples 
#' data(GDP)
#' out  <- splinedecomp(GDP,4,2,4,Plot=FALSE) 
#' @export  
 

splinedecomp <- function(x,d,alpha,beta,Plot=FALSE){  
   n <- length(x)  
   IPR <- Matrix::bandSparse(3*n-2-d+1, 2*n, n, diagonals=list(rep(1,n)))  
   IPR[1:n,1:n] <- Matrix::bandSparse(n,n,0,diagonals=list(rep(1,n)))     
   P <- Matrix::bandSparse(n-2, n, c(0,1,2), diagonals=list(rep(1,n-2),rep(-2,n-2),rep(1,n-2)))  
   IPR[((n+1):(2*n-2)),1:n]<-alpha*P 
   R <- Matrix::bandSparse(n-d+1, n, 0, diagonals=list(rep(1,n-d+1)))   
   for(i in 1:(n-d+1)){ R[i,i:(i+d-1)] <- 1 }
   IPR[(2*n-2+1):(3*n-2-d+1),(n+1):(2*n)]<-beta*R 
   gs <- Matrix::solve(Matrix::crossprod(IPR),Matrix::crossprod(IPR,c(x,rep(0,2*n-2-d+1)))) 
   trend <- gs[1:n] 
   season <- gs[(n+1):(2*n)] 
   residual <- x - (trend+season) 
   gsr <- cbind(trend,season,residual)  
   if(Plot==TRUE){
     if( is.ts(x) == FALSE){ x <- ts(x,start=1,frequency=d) }
     gsr <- ts(gsr,start = tsp(x)[1] , frequency = d )
     oldpar <- par(no.readonly = TRUE) 
     on.exit(par(oldpar)) 
     par(mfrow=c(3,1),mar=c(2,4,1,1)) 
     season <- ts(season, start <- tsp(x)[1] , frequency = d )
     residual <- ts(residual,start <- tsp(x)[1] , frequency = d )
     plot(x,lwd=1.5,ylab="y, Trend",xlab="")
     lines(trend,lwd=1.5)
     plot(season,ylab="Season",xlab="",lwd=1.5)
     lines(c(tsp(x)[1],tsp(x)[1]+length(x)/d),c(0,0),lwd=1.5)
     par(mar=c(4,4,1,1))
     plot(residual,ylab="Residual",xlab="t",lwd=1.5)  
    } 
   return( as.matrix(gsr) )
}  


#' \code{robsplinedecomp} decomposes a vector into trend, season and irregular  component
#' by robustified spline approach; a time series attribute is lost
#'  
#' @param  y  the series, a vector or a time series  
#' @param  d  seasonal period
#' @param  alpha smoothing parameter for trend component   (the larger alpha is, the smoother will the smooth component g be) 
#' @param  beta  smoothing parameter for  seasonal  component
#' @param  Plot  logical, should a plot be produced?
#' @return out   list with the elements trend, season, residual 
#'
#' @examples 
#' data(GDP) 
#' out  <- robsplinedecomp(GDP,4,2,10,Plot=FALSE) 
#' @export  
 

robsplinedecomp <- function(y,d,alpha,beta,Plot=FALSE){
n <- length(y)

if( is.ts(y) == TRUE){ y <- as.vector(y) }  #  not-yet-implemented method for -(<ts>, <dgeMatrix>)

# determination of start values gm, sm

y1 <- y
n1 <- d*floor(n/d) 
if( n1 < n ){ 
  r <- n%%d
  y1 <- c(y[1:n],y[(n1-(d-r)+1):n1]) 
  }
 
teily <-  matrix(y1,d,length(y1)/d)  

for(s in c(1:d)){
  y1 <- teily[s,]
  teily[s,] <- smoothrb(y1,beta,q=min(5,floor(length(y1)/2)))
  }

y1 <- as.vector(teily)[1:n]  

out  <- splinedecomp(y1,d,alpha,beta)  
gm1 <- out[,1]
sm1 <- out[,2]
sigmau <- sd(out[,3])  

P <- Matrix::bandSparse(n-2,n,k=c(0:2),diagonals=list(rep(1,n),rep(-2,n),rep(1,n))) 
R <- Matrix::bandSparse(n-d+1,n,k=c(0:(d-1)),diagonals=matrix(1,n-d+1,d))
I <- Matrix::bandSparse(n,n,k=0,diagonals=list(rep(1,n))) 
PP <- Matrix::crossprod(P) 
A <- I + (alpha^2)*PP
RR <- Matrix::crossprod(R)
B <- I + (beta^2)*RR
S <- A%*%B - I

gm <- 0
sm <- 0  
it <- 0
 
while( (it < 50) & (max(abs(c(as.vector(gm - gm1),as.vector(sm-sm1)))) > (0.001*sigmau)) ){ 
it <- it+1
gm <- gm1
sm <- sm1

psi1m <- sigmau*psihuber((y - gm - sm)/sigmau) - (alpha^2)*PP%*%matrix(gm,n,1)
psi2m <- sigmau*psihuber((y - gm - sm)/sigmau) - (beta^2)*RR%*%matrix(sm,n,1) 
b  <- A%*%matrix(psi2m,length(psi2m),1) - psi1m 
dsm <- Matrix::solve(S,b) 
b <- psi1m - dsm
dgm <- Matrix::solve(A,b)

gm1 <- gm + dgm
sm1 <- sm + dsm
} 

out <- list(trend=gm1,season=sm1, residual = y - gm1 - sm1) 
   if(Plot==TRUE){
     oldpar <- par(no.readonly = TRUE)   
     on.exit(par(oldpar))  
     par(mfrow=c(3,1),mar=c(2,4,1,1))  
     plot(c(1:n),y,type="l",lwd=1.5,ylab="y, Trend",xlab="")
     lines(c(1:n),out$trend,lwd=1.5)
     plot(c(1:n),out$season,type="l",ylab="Season",xlab="",lwd=1.5)
     lines(c(1,n),c(0,0),lwd=1.5)
     par(mar=c(4,4,1,1))
     plot(c(1:n),out$residual,type="l",ylab="Residual",xlab="t",lwd=1.5) 
     lines(c(1,n),c(0,0),lwd=1.5) 
    }  
return(out)
}
 

#' \code{simpledecomp} decomposes a vector into trend, season and irregular component
#' by linear regression approach
#'  
#' @param  y      the series, a vector or a time series  
#' @param  trend  order of trend polynomial 
#' @param  season period of seasonal component  
#' @param  Plot    logical, should a plot be produced?
#' @return  out:  (n,3) matrix
#' \item{1. column}{smooth component }
#' \item{2. column}{seasonal component }
#' \item{3. column}{irregular component  }
#' @examples 
#' data(GDP)
#' out  <- simpledecomp(GDP,trend=3,season=4,Plot=FALSE) 
#' @export  
  

simpledecomp <- function(y,trend=0,season=0,Plot=FALSE){ 
  p <- trend
  s <- season
  n <- length(y) 
  if ((p<0)|(!(p%%1)==0)){  stop("Wrong input for trend (>=0, integer)")  }
  if ((s<0)|(!(s%%1)==0)){  stop("Wrong input for season (>=0,  integer)")  }
  if (n<s){ stop("Length of season must not be greater than the length of the series.") }
  trend <- matrix(mean(y),n,1)  
  season <- matrix(0,n,1)
  if (p>0) {
            xx<-matrix(1,n,p)  
            for (i in 1:p){xx[,i]<-seq(1,n)^i}
            out <- lm(y ~ xx )
            trend <- out$coefficients[1]+xx%*%out$coefficients[2:(p+1)]
           }        
  if (s>1){  
            ss<-matrix(0,n,s-1)  
            for (t in 1:n){ 
                            if ((t%%s)>0) { ss[t,(t%%s)] <- 1 }
                            if ((t%%s)==0){ ss[t,1:s-1] <- -1 }
                          }  
            y.ohne.t <- y-trend
            out <- lm(y.ohne.t ~ ss-1) 
            season <- ss%*%out$coefficients            
            }  
  if ((p>0) & (s<=1)) { residual <- y - (trend+season) }   
  if ((p==0) & (s>1)) { residual <- y - (trend+season) }     
  if ((p>0) & (s>1))  { 
                        regressor <- matrix(0,n,p+s-1)
                        regressor[,1:p] <- xx
                        regressor[,(p+1):(p+s-1)] <- ss
  out <- lm(y~regressor) 
  trend <- out$coefficients[1]+xx%*%out$coefficients[2:(p+1)]
  season<-ss%*%out$coefficients[(p+2):(p+s)]  
  } 
  residual <- y - (trend+season)
  if(is.ts(y) == FALSE){y <- ts(y,start=1,frequency=1)}
  trend <- ts(trend, start=tsp(y)[1], frequency=tsp(y)[3]) 
  season <- ts(season, start=tsp(y)[1], frequency=tsp(y)[3])
  residual <- ts(residual, start=tsp(y)[1], frequency=tsp(y)[3])
  if(Plot == TRUE){
    oldpar <- par(no.readonly = TRUE)   
    on.exit(par(oldpar))  
    par(mfrow=c(3,1),mar=c(2,4,1,1))
    plot(y,lwd=1.5,ylab="y, Trend")
    lines(seq(tsp(y)[1],tsp(y)[2],1/tsp(y)[3]),trend,lwd=1.5)
    plot(season,type="h",ylab="Season",xlab="",lwd=1.5)
    lines(c(tsp(y)[1],tsp(y)[1]+length(y)/tsp(y)[3]),c(0,0),lwd=1.5)
    par(mar=c(4,4,2,1))
    plot(residual,xlab="t",ylab="Residual",lwd=1.5) 
  }
  out <- cbind(trend,season,residual)
  return(out)
  }
