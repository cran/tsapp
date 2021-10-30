##
##  multivarts.r                                                                  R. Schlittgen    07.01.2020
##
##     functions  for  analysis of multivariate series


#' \code{bispeces} performs indirect bivariate spectral estimation of two series y1, y2 using lagwindows
#'  
#' @param  y1   vector, the first time series  
#' @param  y2   vector, the second time series  
#' @param  q    number of covariances used for indirect spectral estimation
#' @param  win  lagwindow  (possible: "bartlett", "parzen", "tukey")
#' @return out  data frame with columns:
#' \item{f}{frequencies 0, 1/n, 2/n, ...  (<= 1/2 )  }
#' \item{coh}{estimated coherency at Fourier frequencies 0,1/n, ...}
#' \item{ph}{estimated phase at Fourier frequencies 0,1/n, ...} 
#' @examples 
#' data(ICECREAM)
#' y <- ICECREAM
#' out <- bispeces(y[,1],y[,2],8,win="bartlett") 
#' @export 

bispeces <- function(y1,y2,q,win="bartlett"){           
n <- length(y1) 
if(win=="bartlett"){  w <- lagwinba(q) }            
if(win=="parzen"){  w <- lagwinpa(q) }               
if(win=="tukey"){  w <- lagwintu(q) } 
a1 <- acf(y1,q,plot=FALSE)                            
a1 <- a1$acf*w 
a2 <- acf(y2,q,plot=FALSE)                            
a2 <- a2$acf*w 
out1 <- periodogram(a1,n/2,ACF=TRUE)                     
s1 <- out1[,2] 
out2 <- periodogram(a2,n/2,ACF=TRUE)                  
s2 <- out2[,2] 
cc <-  ccf(y1,y2,lag.max=q,plot=FALSE)                
cc <- cc$acf*c(rev(w),w[-1])                          
lags <- c(-q:q) 
j <- -1                                               
f <- NULL 
coh <- NULL                                           
ph <- NULL 
while(j < (n/2-1)){                                   
  j = j+1; 
  f <- append(f,(j/n))                                 
  s <- sum(cc*exp((1i)*2*pi*lags*j/n))                 
  coh <- append(coh,(abs(s)/sqrt(s1[j+1]*s2[j+1])))    
  ph <- append(ph,atan(-Im(s)/Re(s))) } 
return(cbind(f,coh,ph)) 
} 
             
             
#' \code{acfmat} computes a sequence of autocorrelation matrices for a multivariate time series 
#'  
#' @param  y        multivariate time series  
#' @param  lag.max  maximum number of lag  
#' @return  out   list with components: 
#' \item{M}{ array with autocovariance matrices }
#' \item{M1}{array with indicators if autocovariances are significantly greater (+),    
#'                  lower (-) than the critical value or insignificant (.) at 95 percent level        } 
#' @examples 
#' data(ICECREAM)
#' out <- acfmat(ICECREAM,7) 
#' @export              
  

acfmat <-function(y,lag.max){
P <- lag.max
y <- as.matrix(unname(y))
komp <- ncol(y)
n <- nrow(y)
M <- array(0,c(komp,komp,P)) 
M1 <- M 
  for (i in c(1:komp)){ 
    for (j in c(1:komp)){ 
      if (i != j){ 
        cc <- ccf(y[,i],y[,j],P,plot=FALSE) 
        M[j,i,] <- rev(cc$acf[1:P])
        } 
       else{ 
        a <- acf(y[,i],P,plot=FALSE) 
        M[j,i,] <- a$acf[-1]
        } 
       for (ord in c(1:P)){
         if( M[j,i,ord] > 2/sqrt(n-ord) ){ M1[j,i,ord] <- "+"  }
         else if( M[j,i,ord] < -2/sqrt(n-ord) ){ M1[j,i,ord] <- "-"  } 
         else{M1[j,i,ord] <- "."  }  
      }
   } }   
  return(list(M=M,M1=M1))
} 

             
#' \code{pacfmat} sequence of partial autocorrelation matrices and related statistics for a multivariate time series
#'  
#' @param  y    multivariate time series  
#' @param  lag.max  maximum number of lag  
#' @return  out  list with components: 
#'  \item{M}{array with matrices of partial autocovariances divided by their standard error}
#'  \item{M1}{array with indicators if partial autocovariances are significantly greater (+), lower (-) than the   
#'          critical value or insignificant (.)  }
#'  \item{R}{array with matrices of partial autocovariances  }
#'  \item{S}{matrix of diagonals of residual covariances (row-wise) }
#'  \item{Test}{test statistic  }
#'  \item{pval}{p value of test} 
#' @examples 
#' data(ICECREAM)
#' out <- pacfmat(ICECREAM,7) 
#' @export              
  

pacfmat <- function(y,lag.max){
P <- lag.max
CN <- colnames(y)
if( is.null(colnames(CN) ) ){ 
    CN <- "Y1"
    for(i in c(2:ncol(y))){ CN <- c(CN,paste("Y",i,sep="")) } 
    }
y <- as.matrix(unname(y))
komp <- ncol(y)
n <- nrow(y)
M <- array(0,c(komp,komp,P))         # becomes quotient pacf/stderror 
M1 <- M                                          # becomes indicator-matrix for part. correl.
R <- M                                             # becomes pacf-matrices
Test  <- rep(0,P+1)                         # becomes test statistic 
S<- matrix(0,P,ncol(y))                     # becomes diagonal  of residual covmat (rowwise)
Test[1] <- det(t(scale(y,center=TRUE, scale = FALSE))%*%scale(y,center=TRUE, scale = FALSE))
colnames(y) <- CN
for (ord in c(1:P)){  
   out <- vars::VAR(y,p=ord)   
   S[ord,] <- diag(var(resid(out)) ) 
   Test[ord+1] <- det(t(as.matrix(resid(out)))%*%as.matrix(resid(out)))
   for (i in c(1:komp)){  
      m <- matrix(unlist(coef(out)[i]),komp*ord+1,4)  
      R[i,,ord] <- m[(komp*(ord-1)+1):(komp*ord),1] 
      M[i,,ord] <- m[(komp*(ord-1)+1):(komp*ord),1]/m[(komp*(ord-1)+1):(komp*ord),2] 
      } 
    }      
  M1[M > 2] <- "+"  
  M1[M < -2] <- "-"
  M1[(M >= -2)&(M <= 2)] <- "."
  Test <- Test[-1]/Test[-(P+1)]
  Test <- -(n-0.5-c(1:P)*komp)*log(Test) 
  pval <- 1-pchisq(Test,komp^2)
  out <- list(M=M,M1=M1,R=R,S=S,Test=Test,pval=pval)
  return(out)
}  

#' \code{Grangercaus}  determines three values of BIC from a twodimensional VAR process 
#'  
#' @param  x  first time series  
#' @param  y  second time series  
#' @param  p  maximal order of VAR process
#' @return  out  list with components
#'  \item{BIC}{vector of length 3: }
#'   \tabular{ll}{
#'    BIC1 \tab minimum aic value for all possible lag structures \cr
#'    BIC2 \tab minimum aic value when Y is not included as regressor in the equation for X \cr
#'    BIC3 \tab minimum aic value when X is not included as regressor in the equation for Y  }
#'  \item{out1}{output of function lm for regression equation for x-series }
#'  \item{out2}{output of function lm for regression equation for y-series  }
#' 
#' @examples 
#' \donttest{
#' data(ICECREAM)
#' out <- Grangercaus(ICECREAM[,1],ICECREAM[,2],3) 
#' }
#' @export    

Grangercaus <- function(x,y,p){ 

X <- x
Y <- y
x<- as.vector(x)
y<- as.vector(y)
mat<-subsets(p)
p2<-2^p
ind <- c(1:p) 
 
tsmat1 <- tsmat(x,p+1)[,-1]
tsmat2 <- tsmat(y,p+1)[,-1] 
x <- x[-c(1:p)]
y <- y[-c(1:p)]

namx <- NULL
for (j in (1:p)){namx <- append(namx, paste("X[t-",j,"]",sep = "")) }
colnames(tsmat1) <- namx
namy <- NULL
for (j in (1:p)){namy <- append(namy, paste("Y[t-",j,"]",sep = "")) }
colnames(tsmat2) <- namy 
 
n <- length(x)

BIC1 <- NA
BIC2 <- NA
BIC3 <- NA

for (i in c(1:p2)){
for (j in c(1:p2)){
for (k in c(1:p2)){
for (l in c(1:p2)){
 
dat1<-as.data.frame(cbind(x,tsmat1[,ind*mat[i,]],tsmat2[,ind*mat[j,]]))
dat2<-as.data.frame(cbind(y,tsmat1[,ind*mat[k,]],tsmat2[,ind*mat[l,]]))
out1 <- lm(x ~ .,data=dat1)
out2 <- lm(y ~ .,data=dat2)
D <- det(var(cbind(out1$residuals,out2$residuals)))
m <- log(D) + log(n-p)*(ncol(dat1)+ncol(dat2))/(n-p)

if (is.na(BIC1)) { BIC1 <- m 
                   ort1 <- c(i,j,k,l) }
if (m < BIC1){ BIC1 <- m
               ort1 <- c(i,j,k,l)
             } 
if (is.na(BIC2)&(j==1)) { BIC2 <- m 
                           ort2 <- c(i,j,k,l) }            
if ((m < BIC2)&(j==1)) { BIC2 <- m
                                ort2 <- c(i,j,k,l)
                               }
if (is.na(BIC3)&(k==1)) { BIC3 <- m
                         ort3 <- c(i,j,k,l) }    
if ((m < BIC3)&(k==1)) { BIC3 <- m
                                ort3 <- c(i,j,k,l)
                               }                       
}}}} 
 
dat1<-as.data.frame(cbind(x,tsmat1[,ind*mat[ort1[1],]],tsmat2[,ind*mat[ort1[2],]]))
colnames(dat1)<-c("x",namx[ind*mat[ort1[1],]],namy[ind*mat[ort1[2],]])
dat2<-as.data.frame(cbind(y,tsmat1[,ind*mat[ort1[3],]],tsmat2[,ind*mat[ort1[4],]]))
colnames(dat2)<-c("y",namx[ind*mat[ort1[3],]],namy[ind*mat[ort1[4],]])
out1 <- lm(x ~ .,dat1)
out2 <- lm(y ~ .,dat2)
out <- list(BIC=c(BIC1,BIC2,BIC3),out1=out1,out2=out2)
return(out)
}   
