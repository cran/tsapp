##
##  specialfun.r                                                                  R. Schlittgen    07.01.2020
##
##    special functions for  time series analysis  


#' \code{ldrec} does Levinson-Durbin recursion for determing all coefficients a(i,j)
#'  
#' @param a (p+1,1)-vector of acf of a time series: acov(0),...,acov(p)
#'                              or  1,acor(1),..,acor(p)
#' @return  mat (p,p+2)-matrix, coefficients in lower triangular,
#'  pacf in colum p+2 and Q(p) in colum p+1
#'
#' @examples 
#' data(HEARTBEAT)
#' a <- acf(HEARTBEAT,5,plot=FALSE)
#' mat <- ldrec(a$acf)
#' @export  


ldrec <- function(a){ 
  p <- length(a)-1 
  mat <- matrix(0,p,p+2) 
  acor <- a/a[1]
  acor <- acor[-1] 
  mat[1,1] <- acor[1];  mat[1,p+1] <- 1-acor[1]^2;   mat[1,p+2] <- acor[1] 
  i <- 1 
  while(i < p){
    mat[i+1,i+1] <- (acor[i+1] - sum(mat[i,1:i]*acor[i:1]))/mat[i,p+1] 
    mat[i+1,1:i] <- mat[i,1:i]-mat[i+1,i+1]*mat[i,i:1]   
    mat[i+1,p+2] <- mat[i+1,i+1] 
    mat[i+1,p+1] <- mat[i,p+1]*(1-mat[i+1,p+2]^2) 
    i <- i+1;
 }
return(mat)
}


#' \code{outidentify} performs one iteration of Wei's iterative procedure to identify impact, locations and type 
#'  of outliers in arma processes
#'  
#' @param  x            vector, the time series 
#' @param  object    output of a model fit with the function  arima  (from stats)
#' @param  alpha      the level of the tests for deciding which value is to be considered an outlier 
#' @param robust     logical, should the standard error be computed robustly?
#' @return  out list with   elements
#' \item{outlier}{matrix with time index (ind), type of outlier  (1 = AO, 2 = IO) and value of test statistic (lambda)} 
#' \item{arima.out}{output of final arima model where the outliers are incorporated as fixed regressors    }
#' @examples 
#' data(SPRUCE)
#' out <- arima(SPRUCE,order=c(2,0,0))
#' out2 <- outidentify(SPRUCE,out,alpha=0.05, robust = FALSE)
#' @export  
 

outidentify <- function(x,object,alpha=0.05, robust = FALSE){
 
 e <- residuals(object)
 n <- length(e)
 outlier <- matrix(0,1,3)  
 colnames(outlier) <-c("ind","type","lambda")
 I <- NULL 
 colnamI <- NULL
 
#     step 1      
 
  piwt <- ARMAtoMA(ar = -object$model$theta, ma = -object$model$phi, lag.max = length(e) - 1)   
  piwt <- c(1, piwt)
  if (robust) { sigma = sqrt(pi/2) * mean(abs(e), na.rm = TRUE)
     }else{ sigma = object$sigma2^0.5  }

#     step 2       
 
  lengaus  <- 0
  lengausalt <- -1 
  while(lengausalt < lengaus){  
    lengausalt <- lengaus
    omega = filter(c(0 * e[-1], rev(e)), filter = piwt, sides = 1, method = "convolution")
    omega = omega[!is.na(omega)]
    rho2 = 1/cumsum(piwt^2)
    omega = omega * rho2
    lambda1T = omega/(sigma*sqrt(rho2)) 
    lambda1T = rev(lambda1T)                        # ao
  
    lambda2T = e/sigma                              # io
    
    if(any(c(max(abs(lambda1T)),max(abs(lambda2T))) >= qnorm(1 - alpha/2))){
      ltyp <- which.max(c(max(abs(lambda1T)),max(abs(lambda2T))))  
      if(ltyp==1){ 
          TT <- which.max(abs(lambda1T))  
          if(!any(outlier[,1] == TT))
          { 
            outlier <- as.matrix(rbind(outlier,c(TT,ltyp,lambda1T[TT]))) 
            e  <- e - omega[TT]*c(rep(0,TT-1),piwt)[1:n] 
            I <- cbind(I,c(rep(0,TT-1),piwt)[1:n])
            colnamI <-c(colnamI,"AO")
            }
         }else{ 
          TT <- which.max(abs(lambda2T)) 
        if(!any(outlier[,1] == TT))
        { 
          outlier <- as.matrix(rbind(outlier,c(TT,ltyp,lambda2T[ TT ])))   
          e[TT] <- 0 
          i <- rep(0,n)
          i[TT] <- 1
          I <- cbind(I,i)
          colnamI <- c(colnamI,"IO")
         } 
        }
      }
   if (robust) { sigma = sqrt(pi/2) * mean(abs(e), na.rm = TRUE)
      }else{ sigma = sqrt(sum(e^2)/n)   }  
   lengaus <-  length(outlier[,1])   
 }   
 outlier <- outlier[-1,,drop=F] 
 o <- order(outlier[,1] )
 outlier <- outlier[o,,drop=F] 
 ind <- outlier[,1]     
 colnames(I) <- as.character(ind)
 arima.out <- arima(x,order=c(length(object$model$phi),0,length(object$model$theta)),xreg=as.matrix(I)) 
 out <- list(outlier = outlier, arima.out = arima.out) 
 return(out)
}


#' \code{missar} Substitution of missing values in a time series by 
#'    conditional  exspectations of AR(p) models
#'  
#' @param x vector, the time series 
#' @param p integer, the maximal order of ar polynom   0 < p < 18,  
#' @param   iterout  if = 1, iteration history is printed
#' @return  out list with   elements
#' \item{a}{(p,p)-matrix, estimated ar coefficients for ar-models   }
#' \item{y}{(n,1)-vector, completed time series }
#' \item{iterhist}{matrix, NULL or the iteration history}
#' @source Miller R.B., Ferreiro O. (1984) <doi.org/10.1007/978-1-4684-9403-7_12> "A Strategy to Complete a Time Series with Missing Observations" 
#' @examples 
#' data(HEARTBEAT)
#' x <- HEARTBEAT
#' x[c(20,21)] <- NA
#' out <- missar(x,2)
#' @export  
  

missar <- function(x,p,iterout=0){  
 iterhist <- NULL
 tol <- 0.001 
 n <- length(x)  
 a <- 0 
 voll <- x  
 if(length(x[!is.na(x)]) == n){ stop("no missing observation")  } 
 if(is.na(x[1])|is.na(x[n])){ stop("First and last observation must not be missing")  }   
 if( (p <=  0)|(p >= 18) ){ stop("order p of ar process must be 0<p<18")  } 

# ... indexf : indexes of missing values
 indexf <- c(1:n)
 indexf <- indexf[is.na(x)] 
 m <- length(indexf) 
 
# ... centering of existing values
  mu <- mean(x,na.rm = TRUE) 
  xcent <- x-mu 
  xcent[indexf] <- 0

# ... start values:  yule-walker-estimates for a(1),...,a(p)  


   a <- acf(xcent,p,plot=FALSE) 
   a <- ldrec(a$acf) 
   a <- a[,1:p,drop=FALSE]  
 
 
# ... indexes of the beginnings and endings of the gaps

if(m == 1)                                                # ... only one missing 
   { luecke <- matrix(indexf,1,2) }
else if( (indexf[m]-indexf[1]) == (m-1) )   # ... one sequence of missings of length m 
   { luecke <- matrix(c(indexf[1],indexf[m]),1,2) } 
else                                       # ... else
   { TT <- as.matrix(cbind(indexf[2:m],indexf[1:(m-1)]))
     TT <- TT[(indexf[2:m]-indexf[1:m-1]) != 1,] 
     luecke <- as.matrix(cbind(c(indexf[1],TT[,1]),c(TT[,2],indexf[m])))  
   } 
lz <- nrow(luecke)  
 
# .... iteration  
 
  xcent <- x-mu  
  abbruch <- 0
  wieder <- 0 
  while(abbruch == 0){ 
    xneu <- xcent 
   
    # .... PE-step 
   
    f <- matrix(0,p,p) 
    i <- 0 
    while(i < p){
      i <- i+1  
      acov <-  as.vector(ARMAacf(ma = -a[i,], lag.max = p))   
      f[i,] <- acov[-1]                      # IACF  
      }  
      amiss <- indexf[1]  
      j <- 1 
      while(j < lz){ 
         p1 <- luecke[j+1,1]-1-luecke[j,2] 
         p1 <- min(c(p,p1,(amiss-1)))  
         if(p1 > 1){
           p1 <- min(c(p,p1)) 
           bmiss <- luecke[j,2] 
           xneu[amiss:bmiss] <- pestep(matrix(f[p1,1:p1],p1,1),xneu[(amiss-p1):(bmiss+p1)]) #  !!!! some problems
          }
         else                     # ... linear interpolation  
          { bmiss <- luecke[j,2] 
            linint <- (xneu[bmiss+1]-xneu[amiss-1])/(bmiss-amiss+2) 
            xneu[amiss:bmiss] <- xneu[amiss-1] + c(1:(bmiss-amiss+1))*linint                         
          } 
         amiss <- luecke[j+1,1] 
         j <- j+1 
      } 

      # ......... last gap 
      
      bmiss <- indexf[m] 
      if((n-bmiss) > 1){ 
        p1 <- min(c(p,(n-bmiss))) 
 
        xneu[amiss:bmiss] <- pestep(matrix(f[p1,1:p1],p1,1),xneu[(amiss-p1):(bmiss+p1)]) 
        }
      else  
        { linint <- (xcent[bmiss+1]-xcent[amiss-1])/(bmiss-amiss+2) 
         xneu[amiss:bmiss] <- xneu[amiss-1] + c(1:(bmiss-amiss+1))*linint  
        }  

   # ......... completed series 

   voll <- xneu 

      # .... PM-step 
      aalt <- a  
      a <- acf(voll,p,plot=FALSE) 
      a <- a$acf 
      a <- ldrec(a)[,1:p]   
      # .... check for interruption  
 
      if( max(abs(aalt-a)) <= tol ){ abbruch <- 1  }

      wieder <- wieder+1  
      if( wieder > 20) {  abbruch <- 1  } 

if(iterout  == 1){  
 iterhist <- rbind(iterhist,c( c("Iteration step",wieder),c("PACF:",round(diag(a),5)),c("Max. param-diff. from last iteration:", round(max(max(abs(aalt-a))),6 )))) 
 } 
} 
out <-  list(a=a,y=voll+mu,iterhist = iterhist) 
return(out) 
}
 
#' \code{pestep} help function for missar 
#'  
#' @param  f  IACF, inverse ACF  
#' @param xt  segment of the time series
#' @return xt new version of xt
#' 
#' @export  
  
pestep <- function(f,xt){  
 
   p <- nrow(f)     
   lrand <- xt[p:1]   
   lrand <- matrix(lrand,length(lrand),1)
   rrand <- xt[(length(xt)-p+1):length(xt)] 
   rrand <- matrix(rrand,length(rrand),1)
   xt <- xt[(p+1):(length(xt)-p)]         
   m <- length(xt)  
   indfehl <- c(1:m)              ## Indexes of missing values 
   indfehl <- indfehl[is.na(xt)]  

   if(m > 1){ 
      if( m <= (p+1) ){ arg <- f[1:(m-1)] }
      else{ arg <- c(f,rep(0,m-p-1)) }  

      mat1 <- toeplitz(c(1,arg))  
      if( p > 1){ 
         Tu <- toeplitz(rev(-f)) 
         tt <- row(Tu)<=col(Tu)    
         tt <- t(tt*Tu)[p:1,]  
         b1 <- tt%*%lrand 
         b2 <- (Tu*(row(Tu)>=col(Tu)))%*%rrand
         }
      else{  
         b1 <- -f%*%lrand 
         b2 <- -f%*%rrand 
      } 

      if(m <= p){ 
         b <- b1[1:m] + b2[(p-m+1):p] 
        }else{  
         b <- c(b1,rep(0,m-p)) + c(rep(0,m-p),b2) 
        } 

      ## possibly existing values ########################/
       

      if( length(indfehl) < m ){
         ind <- c(1:m)
         ind <- ind[complete.cases(xt)]  
         b <- b - mat1[,ind]%*%matrix(xt[ind],length(xt[ind]),1) 
         neu <- lm.fit(mat1[indfehl,indfehl],b[indfehl])
         neu <- neu$coefficients 
         } 
         else{ 
         neu <- lm.fit(mat1,b)
         neu <- neu$coefficients  
         } 

   }else{      ##  only one missing value ##################*/
         neu <- sum(-f*(lrand+rrand)) 
        }  
   xt[indfehl] <- neu  
   return(xt) 
} 

#' \code{missls} substitutes missing values in a time series using the LS approach with ARMA models  
#'  
#' @param  x     vector, the time series 
#' @param  p     integer, the  order of polynom   alpha(B)/beta(B)    
#' @param tol    tolerance that can be set; it enters via tol*sd(x,na.rm=TRUE)
#' @param theo   (k,1)-vector, prespecified Inverse ACF, IACF (starting at lag 1) 
#' @return  y    completed time series
#' @source  S. R. Brubacher and G. Tunnicliffe Wilson (1976)   <https://www.jstor.org/stable/2346678>  "Interpolating Time Series with Application to the 
#'    Estimation of Holiday Effects on Electricity Demand Journal of the Royal Statistical Society"
#' @examples 
#' data(HEARTBEAT)
#' x <- HEARTBEAT
#' x[c(20,21)] <- NA
#' out <-  missls(x,p=2,tol=0.001,theo=0)
#' @export  


 missls <- function(x,p=0,tol=0.001,theo=0){ 
 
  n <- length(x)   
  if(length(x[!is.na(x)]) == n){ stop("no missing observation")  } 
  if(is.na(x[1])|is.na(x[n])){ stop("First and last observation must not be missing")  }    
     
  mu <- mean(x, na.rm = TRUE)   
  xcent <- x-mu 
  tol <- tol*sd(x, na.rm = TRUE) 
 
  if(theo==0){                            # fitting of an AR[p] model    
   if(p==0) { p <- trunc(n/10) } 
                                                # estimation of ACF
   indexf <- c(1:n)
   indexf <- indexf[is.na(x)]  
   y <- xcent
   y[indexf] <- 0           
   g <- 1*(!is.na(xcent))     
   ycov <- acf(y,p,type="covariance",demean=FALSE,plot=FALSE) 
   ycov <- ycov$acf  
   gcov <- acf(g,p,type="covariance",demean=FALSE,plot=FALSE)  
   gcov <- gcov$acf  
   xcov <- ycov/gcov 
   xcor <- xcov/xcov[1]  
   mat <- ldrec(xcor)                   # Compute Levinson-Durbin recursion   
   a <- mat[p,1:p]                        # select AR coefficients    
   rho <- as.vector(ARMAacf(ma = -a, lag.max = p))     #  iacf   
   rho <- rho[-1] 
 }else{ 
    rho <- theo  
    p <- length(rho) 
 }  
  wieder <- 0 
  abbruch <- 0  
  while(abbruch == 0){  
       z <- interpol(rho,xcent)   
       if(theo == 0){      
         aneu <- ar(z, aic = FALSE, order.max = p, method= "yule-walker")  
         aneu <- aneu$ar
         if (max(abs(a-aneu)) < tol) { abbruch <- 1  }          
          else{
           a <- aneu  
           rho <- as.vector(ARMAacf(ma = -a, lag.max = p))     # new iacf 
           rho <- rho[-1]   
          }
        }else{ abbruch <- 1 } 
       wieder <- wieder+1 
       if(wieder > 20)  {abbruch <- 1}
       }  
    out <- z+mu 
    return(out) 
 } 

##################################################################*/

#' \code{interpol} help function for missls 
#'  
#' @param rho   autocorrelation function  
#' @param xcent centered time series
#' @return z  new version of xcent
#' 
#' @export  

interpol <- function(rho,xcent){ 
 
   n <- length(xcent)                     
   p <- length(rho)  
   fehl <- c(1:n)
   fehl <- fehl[is.na(xcent)]      
   m <- length(fehl)  
   z <- xcent
   z[fehl] <- 0  
   zt <- rep(0,m)                     # \tilde{z}  
   s <- fehl[1]  
   k <- 1 
   while( k <= m ){
      i <- fehl[k]-s 
      bis1 <- min(c(n-i-s,p))  
      bis2 <- min(c(s+i-1,p))  
      zt[k] <- -( sum(rho[1:bis1]*z[(s+i+1):(bis1+s+i)])
                 + sum(rho[1:bis2]*z[(s+i-1):(s+i-bis2)]) )                
      k <- k+1 
   } 
   mat <- diag(rep(1,m)) 
   k <- 1 
   while( k < m ){
      dist <- fehl[(k+1):m]-fehl[k] 
      if( min(dist) <= p ){
         lp <- c(1:length(dist))
         lp <- lp[(dist > 0)&(dist <= p)] 
         mat[k,k+lp] <- t(rho[dist[lp]])  
         mat[k+lp,k] <- t(mat[k,k+lp])  
      } 
      k <- k+1 
   }    
   neu.lm <- lm.fit(mat,zt)
   z[fehl] <- neu.lm$coefficients 
   return(z)
} 


#' \code{RS} rescaled adjusted range statistic
#'  
#' @param   x  univariate time series
#' @param   k length of the segments for which the statistic is computed. Starting with t=1, the segments do not overlap. 
#' @return  (l,3)-matrix, 1. column: k, second column: starting time of segment, third column: value of RS statistic.
#'
#' @examples  
#'  data(TREMOR)
#'  R <- RS(TREMOR,10)  
#' @export 

RS <- function(x,k){
   y <- cumsum(x)
   n <- length(x)
   nk <- length(k) 
   R <- c(NA,NA,NA) 
   for(ik in c(1:nk)){  
     Ik <- c(0:k[ik])/k[ik]
     t <- 1- k[ik]  
     while(t < (n-2*k[ik])){
      t <- t+k[ik]  
       zw <- y[t:(t+k[ik])] - y[t] - Ik*(y[t+k[ik]] - y[t])
       zw <- max(zw) - min(zw)
       zw <- zw/(sqrt((k[ik]-1)/k[ik])*sd(x[(t+1):(t+k[ik])])) 
       R  <-  rbind(R,c(k[ik],t,zw))
       }
   }    
   return(R[-1,])
   }