

##
##  auxilary.r                                                                  R. Schlittgen    07.01.2020
##
##  auxilary functions

#' \code{tsmat} constructs a (n-p+1,p) matrix from a time series 
#' where the first column is the shortened series y[p],...,y[n], the second is y[p-1],...,y[n-1],  etc.
#'  
#' @param  y the series, a vector or a time series of length n
#' @param  p desired number of columns 
#' @return mat (n-p+1,p) matrix  
#'
#' @examples 
#' out <- tsmat(c(1:20),4)
#'
#' @export 
  
tsmat <- function(y,p){
         n <- length(y) 
         mat <- matrix(0,n-p+1,p) 
         mat[,p] <- y[1:(n-p+1)]
         y1 <- y 
         for(j in c((p-1):1)){
         y1 <- y1[-1]   
         mat[,j] <- y1[1:(n-p+1)]  
         }     
       return(mat) 
       }  
         
         
 
#' \code{subsets} determines all subsets of  a set of n elements (labelled by  1,2,...,n ).   
#'  
#' @param  n   scalar, integer >= 1 
#' @return mat (2^n,n)-matrix, each row gives the membership indicators of the elements 1,2,...,n   
#'
#' @examples 
#' out <- subsets(4)
#'
#' @export 
 

subsets <- function(n){  
  mat <- matrix(c(0,1),2,1)  
  i <- 1
  while (i < n){
    i<-i+1  
    co <- ncol(mat)
    ro <- nrow(mat)
    mat2 <- matrix(0,2*ro,co+1)
    mat2[1:ro,1:co] <- mat 
    mat2[(ro+1):(2*ro),1:co] <- mat 
    mat2[(ro+1):(2*ro),co+1] <- 1 
    mat<-mat2
    }
 mat
}

#' \code{psihuber} is a psi-function for robust estimation  
#'  
#' @param  u   vector  
#' @return out transformed vector 
#'
#' @examples 
#' out <- psihuber(c(3.3,-0.7,2.1,1.8)) 
#'
#' @export 
 

psihuber <- function(u){
  u[u > 1.5] <- 1.5
  u[u < -1.5] <- -1.5
  return(u)
}

#' \code{psifair} is a psi-function for robust estimation  
#'  
#' @param  u    vector  
#' @return out  transformed vector 
#'
#' @examples  
#' out <- psifair(c(3.3,-0.7,2.1,1.8))
#'
#' @export 

psifair <- function(u){ 
  return( u/(1+abs(u/1.4))  )
}  
  

#' \code{polymake} generates the coefficients of an AR process given the zeros of the
#' characteristic polynomial. The norm of the roots must be greater than  one  for stationary processes.       
#'  
#' @param   r      vector, the zeros of the characteristic polynomial  
#' @return   C      coefficients (a[1],a[2],...,a[p]) of the polynomial  1 - a[1]z -a[2]z^2 -...- a[p]z^p  
#'  
#' @examples 
#' C <- polymake(c(2,-1.5,3))
#'
#' @export 
 

polymake<-function(r){
 n <- length(r)
 C <- rep(0,n+1)
 j <- 1
 C[1] <- 1
 while(j <= n){
	     C <- C - r[j]*c(0,C[1:n])  
	     j <- j+1
	     }
	    C <- rev(C)/C[n+1]
	    C <- -C[c(2:(n+1))]
	 return(C)
 } 

#' \code{aspectratio}  determines the aspect ratio to plot a time series       
#'  
#' @param   y      time series
#' @return   a      scalar, the aspect ratio 
#'  
#' @examples 
#' data(GDP)
#' a <- aspectratio(GDP) 
#' @export 

aspectratio<-function(y){
   n <- length(y)
   a <- 2
   s <- y[-1]-y[-n] 
   ry<- max(y)-min(y)
   rx <- n-1 
   s <- (s/ry)*rx 
   w <- sqrt(s^2+1)
   plhelp1 <- function(a){ a1 <- (weighted.mean(abs(atan(s[s>0]/a)),w[s>0])-1)^2 ; a1 }
   plhelp2 <- function(a){ a1 <- (weighted.mean(abs(atan(s[s<0]/a)),w[s<0])-1)^2 ; a1 }  
     a1 <- nlm(plhelp1,1) 
     a2 <- nlm(plhelp2,1) 
   if (sd(atan(s[s>0]/a1$estimate)) <= sd(atan(s[s<0]/a2$estimate)) ){ a <- a1$estimate } else { a <- a2$estimate }
   return(a)
  }


         