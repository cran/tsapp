## HAC.r                                                                  C. Sattarhoff   07.01.2020
#'  HAC Covariance Matrix Estimation

#' \code{HAC}  computes the central quantity (the meat) in the HAC covariance matrix estimator, also called
#' sandwich estimator. HAC is the abbreviation for "heteroskedasticity and autocorrelation consistent".
#'  
#' @param  mcond    a q-dimensional multivariate time series. In the case of OLS regression with q regressors mcond 
#'                  contains the series of the form regressor*residual (see example below).
#' @param  method   kernel function, choose between "Truncated", "Bartlett", "Parzen", "Tukey-Hanning", "Quadratic Spectral".
#' @param  bw       bandwidth parameter, controls the number of lags considered in the estimation. 
#' @return mat      a (q,q)-matrix
#'
#' @source  Heberle, J. and Sattarhoff, C. (2017) <doi:10.3390/econometrics5010009>  "A Fast Algorithm for the Computation of HAC
#'    Covariance Matrix Estimators" 
#' @examples  
#' data(MUSKRAT)
#' y <- ts(log10(MUSKRAT))
#' n <- length(y)
#' t <- c(1:n)
#' t2 <- t^2
#' out2 <- lm(y ~ t +t2)
#' mat_xu <- matrix(c(out2$residuals,t*out2$residuals, t2*out2$residuals),nrow=62,ncol=3)
#' hac <- HAC(mat_xu, method="Bartlett", 4)
#' 
#' mat_regr<- matrix(c(rep(1,62),t,t2),nrow=62,ncol=3)
#' mat_q <- t(mat_regr)%*%mat_regr/62
#' vcov_HAC <- solve(mat_q)%*%hac%*%solve(mat_q)/62
#' # vcov_HAC is the HAC covariance matrix estimation for the OLS coefficients. 
#' @export  


HAC <- function(mcond, method="Bartlett", bw){
  
  # dimensions
  dimmcond    <- dim(mcond)
  Nlen        <- dimmcond[1]
  qlen        <- dimmcond[2]
  
  ww <- kweightsHAC(kernel = method, Nlen, bw)
  
  # step 1: eigenvalues of w* using FFT
  ww <- c(1, ww[1:(Nlen-1)], 0, ww[(Nlen-1):1])
  ww <- Re(fftwtools::fftw(ww))
  
  # step 2: build F*
  FF <- rbind(mcond, matrix(0, Nlen, qlen)) # F*
  
  # step 3 (a) - (c)
  if(dim(FF)[2] == 1){
    FF <- matrix(fftwtools::fftw(FF), ncol = 1)
  }else{
    FF <- fftwtools::mvfftw(FF)
  }    
  FF <- FF * matrix(rep(ww, qlen), ncol=qlen) # LambdaV*F*
  FF <- Re(fftwtools::mvfftw(FF, inverse = TRUE)) / (2*Nlen) # CF*
  
  # step 4: build T(w)F
  FF <- FF[1:Nlen,] # TF
  
  # step 5: S^hat
  return((t(mcond) %*% FF) / Nlen) # S^hat
  
  
  # if(dim(FF)[2] == 1){
  #     FF <- matrix(fftw(FF), ncol = 1)
  # }else{
  #     FF <- mvfftw(FF)
  # }
  
}

######################################################################################################################################

#' \code{kweightsHAC}  help function for HAC
#'  
#' @param  kernel kernel function, choose between "Truncated", "Bartlett", "Parzen", "Tukey-Hanning", "Quadratic Spectral".
#' @param  dimN   number of observations
#' @param  bw bandwidth parameter
#' @return  ww weights
#' @export  

kweightsHAC <- function(kernel = c("Truncated", "Bartlett", "Parzen", "Tukey-Hanning", "Quadratic Spectral"), dimN, bw){
  ww <- numeric(dimN)
  switch(kernel,
         Truncated = {
           ww[1:bw] <- 1
         },
         Bartlett = {
           ww[1:bw] <- 1 - (seq(1,bw) / (bw+1))
         },
         Parzen = {
           seq1 <- (seq(1,floor(bw/2))) / bw
           seq2 <- (seq(floor(bw/2)+1,bw)) / bw
           ww[1:length(seq1)]      <- 1 - 6*seq1^2 + 6*seq1^3
           ww[(length(seq1)+1):bw] <- 2*(1-seq2)^3
         },
         `Tukey-Hanning` = {
           ww[1:bw] <- (1 + cos(pi*((seq(1,bw))/bw))) / 2
         },
         `Quadratic Spectral` = {
           aa <- pi*((seq(1,dimN))/bw)/5
           ww <- 1/(12*aa^2) * (sin(6*aa) / (6*aa) - cos(6*aa))
         }
  )
  return(ww)
}
