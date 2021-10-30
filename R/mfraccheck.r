## mfraccheck.r                                                                  C. Sattarhoff   07.10.2020
#'  multifractal check

#' \code{mfraccheck}  computes the absolute empirical moments of the differenced series for various lags
#'                    and moment orders. E.g. for lag = 3 and moment order = 1 the average absolute value of 
#'                    the differences with lag 3 will be computed. By default, the maximum lag is determined 
#'                    so that the differenced series contains at lest 50 observations. 
#'  
#' @param  p        the series 
#' @param  q_max    maximum moment order 
#' @return out list with components:
#' \item{moments}{matrix with lagmax raws and q_max columns containing the values of the absolute empirical moments }
#' \item{lagmax}{the maximum lag for differencing}
#' @examples 
#' data(NIKKEI)
#' p <- NIKKEI
#' out <- mfraccheck(log(p),5)
#' mom <- ts(out$moments,start=1)
#' ts.plot(mom, log ="xy",xlab="lag",ylab="abs. empirical moments", lty=c(1:5)) 
#' @export   

mfraccheck <- function(p,q_max){
  
  T <- length(p)
  for (d in c(1:(T-1))) {
    n <- (T-1-((T-1) %% d))/d
    if (n < 50) {break}
  }  
  d_max <- d-1
  
  nulls <- rep(0,5*d_max)
  S <- matrix(nulls,d_max,q_max)
  M <- matrix(nulls,d_max,q_max)
  
  for (d in c(1:d_max)){  
    n <- (T-1-((T-1) %% d))/d
    k <- c(1:(n+1))
    i <- 1+(k-1)*d
    r <- diff(p[i])
    
    for (q in c(1:q_max)){
      S[d,q] <- sum((abs(r))^q) 
      M[d,q] <- S[d,q]/n; #abs. empirical moments
    }   
  }  
  return(list(moments = M, lagmax = d_max))  
}