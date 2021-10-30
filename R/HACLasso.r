## HACLasso.r     

#' \code{rlassoHAC} performs Lasso estimation under heteroscedastic and autocorrelated non-Gaussian disturbances.
#'
#' @param  x Regressors (vector, matrix or object can be coerced to matrix).
#' @param  y Dependent variable (vector, matrix or object can be coerced to matrix).
#' @param  kernel Kernel function, choose between "Truncated", "Bartlett" (by default), "Parzen", 
#'         "Tukey-Hanning", "Quadratic Spectral".
#' @param  bands Bandwidth parameter with default bands=10.
#' @param  bns Block length with default bns=10.
#' @param  lns Number of blocks with default lns = floor(T/bns).
#' @param  nboot Number of bootstrap iterations with default nboot=5000.
#' @param  post Logical. If TRUE (default), post-Lasso estimation is conducted, i.e. a refit of the model with the selected variables.
#' @param  intercept Logical. If TRUE, intercept is included which is not penalized.
#' @param  model Logical. If TRUE (default), model matrix is returned. 
#' @param X.dependent.lambda Logical, TRUE, if the penalization parameter depends on the design
#'                                     of the matrix x. FALSE (default), if independent of the design matrix.
#' @param c Constant for the penalty, default value is 2.
#' @param gamma Constant for the penalty, default gamma=0.1/log(T) with T=data length.   
#' @param  numIter Number of iterations for the algorithm for the estimation of the variance and  
#'                 data-driven penalty, ie. loadings.  
#' @param  tol  Constant tolerance for improvement of the estimated variances. 
#' @param  threshold Constant applied to the final estimated lasso coefficients. Absolute values 
#'                   below the threshold are set to zero.
#' @param ...  further parameters 
#' @return  
#' rlassoHAC returns an object of class "rlasso". An object of class "rlasso" is a list containing at least the
#' following components:
#' \item{coefficients}{Parameter estimates.}
#' \item{beta}{Parameter estimates (named vector of coefficients without intercept). }
#' \item{intercept}{Value of the intercept. }
#' \item{index}{Index of selected variables (logical vector). }
#' \item{lambda}{Data-driven penalty term for each variable, product of lambda0 (the penalization parameter) and the loadings. }
#' \item{lambda0}{Penalty term. }
#' \item{loadings}{Penalty loadings, vector of lenght p (no. of regressors). }
#' \item{residuals}{Residuals, response minus fitted values. }
#' \item{sigma}{Root of the variance of the residuals. }
#' \item{iter}{Number of iterations. }
#' \item{call}{Function call. }
#' \item{options}{Options. }
#' \item{model}{Model matrix (if model = TRUE in function call). }
#'
#' @source  Victor Chernozhukov, Chris Hansen, Martin Spindler (2016). hdm: High-Dimensional Metrics,
#' R Journal, 8(2), 185-199. URL https://journal.r-project.org/archive/2016/RJ-2016-040/index.html. 
#' @examples  
#' \donttest{
#' set.seed(1)
#' T = 100 #sample size
#' p = 20 # number of variables
#' b = 5 # number of variables with non-zero coefficients
#' beta0 = c(rep(10,b), rep(0,p-b))
#' rho = 0.1 #AR parameter
#' Cov = matrix(0,p,p)
#' for(i in 1:p){
#'   for(j in 1:p){
#'      Cov[i,j] = 0.5^(abs(i-j))
#'   }
#' } 
#' C <- chol(Cov)
#' X <- matrix(rnorm(T*p),T,p)%*%C
#' eps <- arima.sim(list(ar=rho), n = T+100)
#' eps <- eps[101:(T+100)] 
#' Y = X%*%beta0 + eps
#' reg.lasso.hac1 <- rlassoHAC(X, Y,"Bartlett") #lambda is chosen independent of regressor 
#'                                              #matrix X by default.
#'
#' bn = 10 # block length
#' bwNeweyWest = 0.75*(T^(1/3))
#' reg.lasso.hac2 <- rlassoHAC(X, Y,"Bartlett", bands=bwNeweyWest, bns=bn, nboot=5000,
#'                             X.dependent.lambda = TRUE, c=2.7) 
#' }
#' @export  

rlassoHAC <- function(x, y, kernel="Bartlett", bands=10, bns=10, lns=NULL, nboot=5000, post = TRUE, intercept = TRUE, 
                      model = TRUE, X.dependent.lambda = FALSE, c = 2, gamma=NULL, numIter = 15, tol = 10^-5, threshold = NULL, ...) {
  
  homoscedastic = FALSE
  lambda.start = NULL
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  if(is.null(lns)) { lns=floor(n/bns) }
  if(is.null(gamma)) { gamma=0.1/log(n) }
  
  if (is.null(colnames(x)))
    colnames(x) <- paste("V", 1:p, sep = "")
  ind.names <- 1:p
  
  penalty <- list(homoscedastic = homoscedastic, X.dependent.lambda = X.dependent.lambda, 
      lambda.start = lambda.start, c = c, gamma = gamma)
  
  # set options to default values if missing
  #if (!exists("homoscedastic", where = penalty))  penalty$homoscedastic = "FALSE"
  if (!exists("X.dependent.lambda", where = penalty))  penalty$X.dependent.lambda = "FALSE"
  if (!exists("gamma", where = penalty))  penalty$gamma = 0.1/log(n)
  #if (penalty$homoscedastic=="TRUE") stop("The argument penalty$homoscedastic must be set to FALSE!")
  
  control <- list(numIter = numIter, tol = tol, threshold = threshold)
  
  # checking input numIter, tol
  if (!exists("numIter", where = control)) {
    control$numIter = 15
  }
  
  if (!exists("tol", where = control)) {
    control$tol = 10^-5
  }
  
  #if (post==FALSE & (!exists("c", where = penalty) | is.na(match("penalty", names(as.list(match.call)))))) {
  if (post==FALSE & (!exists("c", where = penalty))) {  
    penalty$c = 0.5
  }
  
  default_pen <-  list(homoscedastic = FALSE, X.dependent.lambda = FALSE, lambda.start = NULL,c = penalty$c, gamma = .1/log(n))
  if (post==FALSE &  isTRUE(all.equal(penalty, default_pen))) {  
    penalty$c = 0.5
  }
  
  # Intercept handling and scaling
  if (intercept) {
    meanx <- colMeans(x)
    x <- scale(x, meanx, FALSE)
    mu <- mean(y)
    y <- y - mu
  } else {
    meanx <- rep(0, p)
    mu <- 0
  }
  
  normx <- sqrt(apply(x, 2, var))
  # Psi <- apply(x, 2, function(x) mean(x^2)) 
  Psi <- sqrt(apply(x,2, function(x) mean(x^2)))
  ind <- rep(FALSE, p) #
  
  XX <- crossprod(x)
  Xy <- crossprod(x, y)
  startingval <- init_values(x,y)$residuals
  
  #==============================================
  # The HAC section starts here 
  #==============================================
  pen <- lambdaCalculationHAC(X.dependent.lambda, c, gamma, kernel, bands, bns, lns,
                              nboot,y = startingval, x = x)
  lambda <- pen$lambda
  Ups0 <- Ups1 <- pen$Ups0
  lambda0 <- pen$lambda0
  #==============================================
  # The HAC section ends here  
  #==============================================

  mm <- 1
  s0 <- sqrt(var(y))
  
  while (mm <= control$numIter) {
    # calculation parameters
    # coefTemp <- hdm::LassoShooting.fit(x, y, lambda, XX = XX, Xy = Xy)$coefficients
    #xn <- t(t(x)/as.vector(Ups1))
    if (mm==1 && post) {
      coefTemp <- hdm::LassoShooting.fit(x, y, lambda/2, XX = XX, Xy = Xy)$coefficients
      # lasso.reg <- glmnet::glmnet(xn, y, family = c("gaussian"), alpha = 1,
      #                            lambda = lambda0/(2*n)/2, standardize = FALSE, intercept = FALSE)
      # lasso.reg <- glmnet::glmnet(x, y, family = c("gaussian"), alpha = 1,
      #                           lambda = lambda0/(2*n)/2, standardize = FALSE, intercept = FALSE,  penalty.factor = Ups1)
    } else {
      coefTemp <- hdm::LassoShooting.fit(x, y, lambda, XX = XX, Xy = Xy)$coefficients
      #lasso.reg <- glmnet::glmnet(xn, y, family = c("gaussian"), alpha = 1,
      #                           lambda = lambda0/(2*n), standardize = FALSE, intercept = FALSE)
      #lasso.reg <- glmnet::glmnet(x, y, family = c("gaussian"), alpha = 1,
      #                           lambda = lambda0/(2*n), standardize = FALSE, intercept = FALSE,  penalty.factor = Ups1)
    }
    #coefTemp <- as.vector(lasso.reg$beta)
    #names(coefTemp) <- colnames(x)
    coefTemp[is.na(coefTemp)] <- 0
    ind1 <- (abs(coefTemp) > 0)
    x1 <- as.matrix(x[, ind1, drop = FALSE])
    if (dim(x1)[2] == 0) {
      if (intercept) {
        intercept.value <- mean(y + mu)
        coef <- rep(0,p+1)
        names(coef) <- c("intercept", colnames(x)) #c("intercept", names(coefTemp))
      } else {
        intercept.value <- mean(y)
        coef <- rep(0,p)
        names(coef) <- colnames(x) #names(coefTemp)
      }
      
      est <- list(coefficients = coef, beta=rep(0,p), intercept=intercept.value, index = rep(FALSE, p),
                  lambda = lambda, lambda0 = lambda0, loadings = Ups0, residuals = y -
                    mean(y), sigma = var(y), iter = mm, call = match.call(),
                  options = list(post = post, intercept = intercept, ind.scale=ind, 
                                 control = control, mu = mu, meanx = meanx))
      if (model) {
        est$model <- x
      } else {
        est$model <- NULL
      }
      est$tss <- est$rss <- sum((y - mean(y))^2)
      est$dev <- y - mean(y)
      class(est) <- "rlasso"
      return(est)
    }
    
    # refinement variance estimation
    if (post) {
      reg <- lm(y ~ -1 + x1)
      coefT <- coef(reg)
      coefT[is.na(coefT)] <- 0
      e1 <- y - x1 %*% coefT
      coefTemp[ind1] <- coefT
    }
    if (!post) {
      e1 <- y - x1 %*% coefTemp[ind1]
    }
    s1 <- sqrt(var(e1))
    
    # residuals e1
    #==============================================
    # The HAC section starts here 
    #==============================================
      #Ups1 <- 1/sqrt(n) * sqrt(t(t(e1^2) %*% (x^2)))
      #lambda <- pen$lambda0 * Ups1
    lc <- lambdaCalculationHAC(X.dependent.lambda, c, gamma,kernel, bands, bns, lns,
                               nboot,y = e1, x = x) 
    Ups1 <- lc$Ups0
    lambda <- lc$lambda
    lambda0 <- lc$lambda0
    #==============================================
    # The HAC section ends here 
    #==============================================
    
    mm <- mm + 1
    if (abs(s0 - s1) < control$tol) {
      break
    }
    s0 <- s1
  }
  
  
  if (dim(x1)[2] == 0) {
    coefTemp = NULL
    ind1 <- rep(0, p)
  }
  coefTemp <- as.vector(coefTemp)
  coefTemp[abs(coefTemp) < control$threshold] <- 0
  ind1 <- as.vector(ind1)
  coefTemp <- as.vector(as.vector(coefTemp))
  names(coefTemp) <- names(ind1) <- colnames(x)
  if (intercept) {
    if (is.null(mu)) mu <-0
    if (is.null(meanx))  meanx <-  rep(0, length(coefTemp))  #<- 0
    if (sum(ind)==0) {
      intercept.value <- mu - sum(meanx*coefTemp)
    } else {
      intercept.value <- mu - sum(meanx*coefTemp) #sum(meanx[-ind]*coefTemp)
    }
  } else {
    intercept.value <- NA
  }
  
  #if (intercept) {
  #  e1 <- y - x1 %*% coefTemp[ind1] - intercept.value 
  #} else {
  #  e1 <- y - x1 %*% coefTemp[ind1]
  #}
  if (intercept) {
    beta <- c(intercept.value, coefTemp)
    names(beta)[1] <- "(Intercept)"
  } else {
    beta <- coefTemp
  }
  
  s1 <- sqrt(var(e1))
  est <- list(coefficients = beta, beta=coefTemp, intercept=intercept.value, index = ind1, lambda = lambda,
              lambda0 = lambda0, loadings = Ups1, residuals = as.vector(e1), sigma = s1,
              iter = mm, call = match.call(), options = list(post = post, intercept = intercept,
                                                             control = control, penalty = penalty, ind.scale=ind,
                                                             mu = mu, meanx = meanx, kernel=kernel, bands=bands,
                                                             bns=bns, lns=lns, nboot=nboot), model=model)
  if (model) {
    x <- scale(x, -meanx, FALSE)
    est$model <- x
  } else {
    est$model <- NULL
  }
  est$tss <- sum((y - mean(y))^2)
  est$rss <- sum(est$residuals^2)
  est$dev <- y - mean(y)
  class(est) <- "rlasso"
  return(est)
}

###############################################################################################################
#' \code{lambdaCalculationHAC} is an auxiliary function for rlassoHAC; it calculates the penalty parameters.
#'  
#' @param X.dependent.lambda Logical, TRUE, if the penalization parameter depends on the design 
#'        of the matrix x. FALSE, if independent of the design matrix (default).  
#' @param c Constant for the penalty with default c = 2 . 
#' @param gamma Constant for the penalty with default gamma=0.1. 
#' @param  kernel String  kernel function, choose between "Truncated", "Bartlett", "Parzen", 
#'        "Tukey-Hanning", "Quadratic Spectral".
#' @param  bands Constant bandwidth parameter.
#' @param  bns Block length.
#' @param  lns Number of blocks.
#' @param  nboot Number of bootstrap iterations.
#' @param  y Residual which is used for calculation of the variance or the data-dependent loadings.  
#' @param  x Regressors (vector, matrix or object can be coerced to matrix).
#' @return  
#' \item{lambda0}{Penalty term}
#' \item{Ups0}{Penalty loadings, vector of length p (no. of regressors)}
#' \item{lambda}{This is lambda0 * Ups0}
#' \item{penalty}{Summary of the used penalty function.}
#' 
#' @source  Victor Chernozhukov, Chris Hansen, Martin Spindler (2016). hdm: High-Dimensional Metrics,
#' R Journal, 8(2), 185-199. URL https://journal.r-project.org/archive/2016/RJ-2016-040/index.html. 
#' @export  

lambdaCalculationHAC <- function(X.dependent.lambda = FALSE,c = 2, gamma = 0.1,
                                 kernel, bands, bns, lns, nboot, y = NULL, x = NULL) {
  homoscedastic = FALSE
  lambda.start = NULL
  penalty = list(homoscedastic = homoscedastic, X.dependent.lambda = X.dependent.lambda, lambda.start = lambda.start, 
                 c = c, gamma = gamma)
  # checkmate::checkChoice(penalty$X.dependent.lambda, c(TRUE, FALSE, NULL))
  # # checkmate::checkChoice(penalty$homoscedastic, c(TRUE, FALSE, "none"))
  # checkmate::checkChoice(penalty$homoscedastic, c(FALSE))
  
  #if (!exists("homoscedastic", where = penalty))  penalty$homoscedastic = "FALSE"
  if (!exists("X.dependent.lambda", where = penalty))  penalty$X.dependent.lambda = "FALSE"
  # if (!exists("c", where = penalty) & penalty$homoscedastic!="none") {
  #   penalty$c = 1.1
  # }
  # if (!exists("gamma", where = penalty) & penalty$homoscedastic!="none") {
  #   penalty$gamma = 0.1
  # }
  
  # heteroscedastic and X-independent
  if (penalty$homoscedastic==FALSE && penalty$X.dependent.lambda == FALSE) {
    p <- dim(x)[2]
    n <- dim(x)[1]
    #lambda0 <- 2*penalty$c*sqrt(n)*sqrt(2*log(2*p*log(n)/penalty$gamma))
    lambda0 <- 2 * penalty$c * sqrt(n) * qnorm(1 - penalty$gamma/(2 * p * 1))  # 1=num endogenous variables
    
    #==============================================================================================================
    # The HAC section starts here 
    #==============================================================================================================
    
    #Ups0 <- 1/sqrt(n) * sqrt(t(t(y^2) %*% (x^2)))
    load = rep(0,p)
    for(i in 1:p){
      load[i] = sqrt(HAC(matrix(x[,i]*y-mean(x[,i]*y), ncol=1), kernel, bands))
    }
    Ups0 <- load
    #==============================================================================================================
    # The HAC section ends here  
    #==============================================================================================================
    
    lambda <- lambda0 * Ups0
  }
  
  # heteroscedastic and X-dependent
  if (penalty$homoscedastic==FALSE && penalty$X.dependent.lambda == TRUE) {
    
    #==============================================================================================================
    # The HAC section starts here  # compute lambda
    #==============================================================================================================
    p <- dim(x)[2]
    n <- dim(x)[1]
    Smatrix.boot  = array(0,dim=c(p,nboot))
    cv <- rep(0,nboot)
    load = rep(0,p)
    
    for (j in 1:nboot){
      res.boot = rnorm(lns) 
      for (i in 1:p){
        Smatrix.boot[i,j] = sum(c(rep(res.boot, each=bns),rep(rnorm(1),n-lns*bns))*x[,i]*y)/sqrt(n)
      }
    }
    
    #==============================================================================================================
    # The HAC section starts here  # compute the penalty loadings
    #==============================================================================================================
    #Ups0 <- 1/sqrt(n) * sqrt(t(t(y^2) %*% (x^2)))
    
    for(i in 1:p){
      load[i] = sqrt(HAC(matrix(x[,i]*y-mean(x[,i]*y), ncol=1),kernel, bands))
    }
    
    for(j in 1:nboot){
      cv[j]=max(abs(Smatrix.boot[,j]/load))
    }
    lambda.boot = 2*quantile(cv, 0.9)*sqrt(n)*penalty$c
    
    
    #==============================================================================================================
    # The HAC section ends here 
    #==============================================================================================================
    Ups0 <- load  
    lambda0 <- lambda.boot
    lambda <- lambda0 * Ups0
    
  }
  
  return(list(lambda0 = lambda0, lambda = lambda, Ups0 = Ups0, method = penalty))
}

####################################################################################
#' \code{init_values} is an auxiliary function for rlassoHAC, for fitting linear models with
#' the method of least squares where only the variables in X with highest correlations
#' are considered; taken from package hdm.
#'
#' @param  X Regressors (matrix or object can be coerced to matrix).
#' @param  y Dependent variable(s).
#' @param  number How many regressors in X should be considered.
#' @param  intercept Logical. If TRUE, intercept is included which is not penalized. 
#' @return  
#' init_values returns a list containing the
#' following components:
#' \item{residuals}{Residuals.}
#' \item{coefficients}{Estimated coefficients.}
#' 
#' @source  Victor Chernozhukov, Chris Hansen, Martin Spindler (2016). hdm: High-Dimensional Metrics,
#' R Journal, 8(2), 185-199. URL https://journal.r-project.org/archive/2016/RJ-2016-040/index.html. 
#' @export  

init_values <- function(X, y, number = 5, intercept = TRUE) {
  suppressWarnings(corr <- abs(cor(y, X)))
  kx <- dim(X)[2]
  index <- order(corr, decreasing = T)[1:min(number, kx)]
  coefficients <- rep(0, kx)
  if (intercept == TRUE) {
    reg <- lm(y ~ X[, index, drop = FALSE])
    coefficients[index] <- coef(reg)[-1]
  } else {
    reg <- lm(y ~ -1 + X[, index, drop = FALSE])
    coefficients[index] <- coef(reg)
  }
  coefficients[is.na( coefficients)] <- 0
  res <- list(residuals = reg$residuals, coefficients = coefficients)
  return(res)
}

#########################################################################################
#' \code{rlassoLoad} performs Lasso estimation under heteroscedastic and autocorrelated non-Gaussian disturbances
#'                   with predefined penalty loadings.
#'
#' @param  x Regressors (vector, matrix or object can be coerced to matrix).
#' @param  y Dependent variable (vector, matrix or object can be coerced to matrix).
#' @param  load Penalty loadings, vector of length p (no. of regressors).
#' @param  bns Block length with default bns=10.
#' @param  lns Number of blocks with default lns = floor(T/bns).
#' @param  nboot Number of bootstrap iterations with default nboot=5000.
#' @param  post Logical. If TRUE (default), post-Lasso estimation is conducted, i.e. a refit of the model with the selected variables.
#' @param  intercept Logical. If TRUE, intercept is included which is not penalized.
#' @param  model Logical. If TRUE (default), model matrix is returned. 
#' @param X.dependent.lambda Logical, TRUE, if the penalization parameter depends on the design
#'                                     of the matrix x. FALSE (default), if independent of the design matrix.
#' @param c Constant for the penalty default is 2.
#' @param gamma Constant for the penalty default gamma=0.1/log(T) with T=data length.   
#' @param  numIter Number of iterations for the algorithm for the estimation of the variance and data-driven penalty.  
#' @param  tol  Constant tolerance for improvement of the estimated variances. 
#' @param  threshold Constant applied to the final estimated lasso coefficients. Absolute values 
#'                   below the threshold are set to zero.
#' @param ...  further parameters 
#' @return  
#' rlassoLoad returns an object of class "rlasso". An object of class "rlasso" is a list containing at least the
#' following components:
#' \item{coefficients}{Parameter estimates.}
#' \item{beta}{Parameter estimates (named vector of coefficients without intercept). }
#' \item{intercept}{Value of the intercept. }
#' \item{index}{Index of selected variables (logical vector). }
#' \item{lambda}{Data-driven penalty term for each variable, product of lambda0 (the penalization parameter) and the loadings. }
#' \item{lambda0}{Penalty term. }
#' \item{loadings}{Penalty loadings, vector of lenght p (no. of regressors). }
#' \item{residuals}{Residuals, response minus fitted values. }
#' \item{sigma}{Root of the variance of the residuals. }
#' \item{iter}{Number of iterations. }
#' \item{call}{Function call. }
#' \item{options}{Options. }
#' \item{model}{Model matrix (if model = TRUE in function call). }
#'
#' @source  Victor Chernozhukov, Chris Hansen, Martin Spindler (2016). hdm: High-Dimensional Metrics,
#' R Journal, 8(2), 185-199. URL https://journal.r-project.org/archive/2016/RJ-2016-040/index.html. 
#' @examples 
#' \donttest{  
#' set.seed(1)
#' T = 100 #sample size
#' p = 20 # number of variables
#' b = 5 # number of variables with non-zero coefficients
#' beta0 = c(rep(10,b), rep(0,p-b))
#' rho = 0.1 #AR parameter
#' Cov = matrix(0,p,p)
#' for(i in 1:p){
#'   for(j in 1:p){
#'      Cov[i,j] = 0.5^(abs(i-j))
#'   }
#' } 
#' C <- chol(Cov)
#' X <- matrix(rnorm(T*p),T,p)%*%C
#' eps <- arima.sim(list(ar=rho), n = T+100)
#' eps <- eps[101:(T+100)] 
#' Y = X%*%beta0 + eps
#' 
#' fit1 =  rlasso(X, Y, penalty = list(homoscedastic = "none",
#'               lambda.start = 2*0.5*sqrt(T)*qnorm(1-0.1/(2*p))), post=FALSE)
#' beta = fit1$beta
#' intercept = fit1$intercept
#' res = Y - X %*% beta - intercept * rep(1, length(Y))
#'
#' load = rep(0,p)
#' for(i in 1:p){
#'   load[i] = sqrt(lrvar(X[,i]*res)*T)
#'   }
#' reg.lasso.load1 <- rlassoLoad(X,Y,load) #lambda is chosen independent of regressor 
#'                                              #matrix X by default.
#'
#' bn = 10 # block length
#' reg.lasso.load2 <- rlassoLoad(X, Y,load, bns=bn, nboot=5000,
#'                             X.dependent.lambda = TRUE, c=2.7)
#' } 
#' @export  

rlassoLoad <- function(x, y, load, bns=10, lns=NULL, nboot=5000, post = TRUE, intercept = TRUE, 
                       model = TRUE, X.dependent.lambda = FALSE, c = 2, gamma=NULL, numIter = 15, tol = 10^-5, threshold = NULL, ...) {
  
  homoscedastic = FALSE
  lambda.start = NULL
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  if(is.null(lns)) { lns=floor(n/bns) }
  if(is.null(gamma)) { gamma=0.1/log(n) }
  
  if (is.null(colnames(x)))
    colnames(x) <- paste("V", 1:p, sep = "")
  ind.names <- 1:p
  
  penalty <- list(homoscedastic = homoscedastic, X.dependent.lambda = X.dependent.lambda, 
                  lambda.start = lambda.start, c = c, gamma = gamma)
  
  # set options to default values if missing
  #if (!exists("homoscedastic", where = penalty))  penalty$homoscedastic = "FALSE"
  if (!exists("X.dependent.lambda", where = penalty))  penalty$X.dependent.lambda = "FALSE"
  if (!exists("gamma", where = penalty))  penalty$gamma = 0.1/log(n)
  #if (penalty$homoscedastic=="TRUE") stop("The argument penalty$homoscedastic must be set to FALSE!")
  
  control <- list(numIter = numIter, tol = tol, threshold = threshold)
  
  # checking input numIter, tol
  if (!exists("numIter", where = control)) {
    control$numIter = 15
  }
  
  if (!exists("tol", where = control)) {
    control$tol = 10^-5
  }
  
  #if (post==FALSE & (!exists("c", where = penalty) | is.na(match("penalty", names(as.list(match.call)))))) {
  if (post==FALSE & (!exists("c", where = penalty))) {  
    penalty$c = 0.5
  }
  
  default_pen <-  list(homoscedastic = FALSE, X.dependent.lambda = FALSE, lambda.start = NULL,c = penalty$c, gamma = .1/log(n))
  if (post==FALSE &  isTRUE(all.equal(penalty, default_pen))) {  
    penalty$c = 0.5
  }
  
  # Intercept handling and scaling
  if (intercept) {
    meanx <- colMeans(x)
    x <- scale(x, meanx, FALSE)
    mu <- mean(y)
    y <- y - mu
  } else {
    meanx <- rep(0, p)
    mu <- 0
  }
  
  normx <- sqrt(apply(x, 2, var))
  # Psi <- apply(x, 2, function(x) mean(x^2)) 
  Psi <- sqrt(apply(x,2, function(x) mean(x^2)))
  ind <- rep(FALSE, p) #
  
  XX <- crossprod(x)
  Xy <- crossprod(x, y)
  startingval <- init_values(x,y)$residuals
  
  #==============================================
  # The loadings section starts here 
  #==============================================
  pen <- lambdaCalculationLoad(X.dependent.lambda, c, gamma, load, bns, lns,
                               nboot,y = startingval, x = x)
  lambda <- pen$lambda
  Ups0 <- Ups1 <- pen$Ups0
  lambda0 <- pen$lambda0
  #==============================================
  # The loadings section ends here  
  #==============================================
  
  mm <- 1
  s0 <- sqrt(var(y))
  
  while (mm <= control$numIter) {
    # calculation parameters
    # coefTemp <- hdm::LassoShooting.fit(x, y, lambda, XX = XX, Xy = Xy)$coefficients
    #xn <- t(t(x)/as.vector(Ups1))
    if (mm==1 && post) {
      coefTemp <- hdm::LassoShooting.fit(x, y, lambda/2, XX = XX, Xy = Xy)$coefficients
      # lasso.reg <- glmnet::glmnet(xn, y, family = c("gaussian"), alpha = 1,
      #                            lambda = lambda0/(2*n)/2, standardize = FALSE, intercept = FALSE)
      # lasso.reg <- glmnet::glmnet(x, y, family = c("gaussian"), alpha = 1,
      #                           lambda = lambda0/(2*n)/2, standardize = FALSE, intercept = FALSE,  penalty.factor = Ups1)
    } else {
      coefTemp <- hdm::LassoShooting.fit(x, y, lambda, XX = XX, Xy = Xy)$coefficients
      #lasso.reg <- glmnet::glmnet(xn, y, family = c("gaussian"), alpha = 1,
      #                           lambda = lambda0/(2*n), standardize = FALSE, intercept = FALSE)
      #lasso.reg <- glmnet::glmnet(x, y, family = c("gaussian"), alpha = 1,
      #                           lambda = lambda0/(2*n), standardize = FALSE, intercept = FALSE,  penalty.factor = Ups1)
    }
    #coefTemp <- as.vector(lasso.reg$beta)
    #names(coefTemp) <- colnames(x)
    coefTemp[is.na(coefTemp)] <- 0
    ind1 <- (abs(coefTemp) > 0)
    x1 <- as.matrix(x[, ind1, drop = FALSE])
    if (dim(x1)[2] == 0) {
      if (intercept) {
        intercept.value <- mean(y + mu)
        coef <- rep(0,p+1)
        names(coef) <- c("intercept", colnames(x)) #c("intercept", names(coefTemp))
      } else {
        intercept.value <- mean(y)
        coef <- rep(0,p)
        names(coef) <- colnames(x) #names(coefTemp)
      }
      
      est <- list(coefficients = coef, beta=rep(0,p), intercept=intercept.value, index = rep(FALSE, p),
                  lambda = lambda, lambda0 = lambda0, loadings = Ups0, residuals = y -
                    mean(y), sigma = var(y), iter = mm, call = match.call(),
                  options = list(post = post, intercept = intercept, ind.scale=ind, 
                                 control = control, mu = mu, meanx = meanx))
      if (model) {
        est$model <- x
      } else {
        est$model <- NULL
      }
      est$tss <- est$rss <- sum((y - mean(y))^2)
      est$dev <- y - mean(y)
      class(est) <- "rlasso"
      return(est)
    }
    
    # refinement variance estimation
    if (post) {
      reg <- lm(y ~ -1 + x1)
      coefT <- coef(reg)
      coefT[is.na(coefT)] <- 0
      e1 <- y - x1 %*% coefT
      coefTemp[ind1] <- coefT
    }
    if (!post) {
      e1 <- y - x1 %*% coefTemp[ind1]
    }
    s1 <- sqrt(var(e1))
    
    # residuals e1
    #==============================================
    # The loadings section starts here 
    #==============================================
    #Ups1 <- 1/sqrt(n) * sqrt(t(t(e1^2) %*% (x^2)))
    #lambda <- pen$lambda0 * Ups1
    lc <- lambdaCalculationLoad(X.dependent.lambda, c, gamma,load, bns, lns,
                                nboot,y = e1, x = x) 
    Ups1 <- lc$Ups0
    lambda <- lc$lambda
    lambda0 <- lc$lambda0
    #==============================================
    # The loadings section ends here 
    #==============================================
    
    mm <- mm + 1
    if (abs(s0 - s1) < control$tol) {
      break
    }
    s0 <- s1
  }
  
  
  if (dim(x1)[2] == 0) {
    coefTemp = NULL
    ind1 <- rep(0, p)
  }
  coefTemp <- as.vector(coefTemp)
  coefTemp[abs(coefTemp) < control$threshold] <- 0
  ind1 <- as.vector(ind1)
  coefTemp <- as.vector(as.vector(coefTemp))
  names(coefTemp) <- names(ind1) <- colnames(x)
  if (intercept) {
    if (is.null(mu)) mu <-0
    if (is.null(meanx))  meanx <-  rep(0, length(coefTemp))  #<- 0
    if (sum(ind)==0) {
      intercept.value <- mu - sum(meanx*coefTemp)
    } else {
      intercept.value <- mu - sum(meanx*coefTemp) #sum(meanx[-ind]*coefTemp)
    }
  } else {
    intercept.value <- NA
  }
  
  #if (intercept) {
  #  e1 <- y - x1 %*% coefTemp[ind1] - intercept.value 
  #} else {
  #  e1 <- y - x1 %*% coefTemp[ind1]
  #}
  if (intercept) {
    beta <- c(intercept.value, coefTemp)
    names(beta)[1] <- "(Intercept)"
  } else {
    beta <- coefTemp
  }
  
  s1 <- sqrt(var(e1))
  est <- list(coefficients = beta, beta=coefTemp, intercept=intercept.value, index = ind1, lambda = lambda,
              lambda0 = lambda0, loadings = Ups1, residuals = as.vector(e1), sigma = s1,
              iter = mm, call = match.call(), options = list(post = post, intercept = intercept,
                                                             control = control, penalty = penalty, ind.scale=ind,
                                                             mu = mu, meanx = meanx, load = load,
                                                             bns=bns, lns=lns, nboot=nboot), model=model)
  if (model) {
    x <- scale(x, -meanx, FALSE)
    est$model <- x
  } else {
    est$model <- NULL
  }
  est$tss <- sum((y - mean(y))^2)
  est$rss <- sum(est$residuals^2)
  est$dev <- y - mean(y)
  class(est) <- "rlasso"
  return(est)
}

###############################################################################################################
#' \code{lambdaCalculationLoad} is an auxiliary function for rlassoLoad; it calculates the penalty parameters
#'                               with predefined loadings.
#'  
#' @param X.dependent.lambda Logical, TRUE, if the penalization parameter depends on the design 
#'        of the matrix x. FALSE, if independent of the design matrix (default).  
#' @param c Constant for the penalty with default c = 2 . 
#' @param gamma Constant for the penalty with default gamma=0.1. 
#' @param  load Penalty loadings, vector of length p (no. of regressors).
#' @param  bns Block length.
#' @param  lns Number of blocks.
#' @param  nboot Number of bootstrap iterations.
#' @param  y Residual which is used for calculation of the variance or the data-dependent penalty.  
#' @param  x Regressors (vector, matrix or object can be coerced to matrix).
#' @return  
#' \item{lambda0}{Penalty term}
#' \item{Ups0}{Penalty loadings, vector of length p (no. of regressors)}
#' \item{lambda}{This is lambda0 * Ups0}
#' \item{penalty}{Summary of the used penalty function}
#' 
#' @source  Victor Chernozhukov, Chris Hansen, Martin Spindler (2016). hdm: High-Dimensional Metrics,
#' R Journal, 8(2), 185-199. URL https://journal.r-project.org/archive/2016/RJ-2016-040/index.html. 
#' @export  

lambdaCalculationLoad <- function(X.dependent.lambda = FALSE,c = 2, gamma = 0.1,
                                  load, bns, lns, nboot, y = NULL, x = NULL) {
  homoscedastic = FALSE
  lambda.start = NULL
  penalty = list(homoscedastic = homoscedastic, X.dependent.lambda = X.dependent.lambda, lambda.start = lambda.start, 
                 c = c, gamma = gamma)
  # checkmate::checkChoice(penalty$X.dependent.lambda, c(TRUE, FALSE, NULL))
  # # checkmate::checkChoice(penalty$homoscedastic, c(TRUE, FALSE, "none"))
  # checkmate::checkChoice(penalty$homoscedastic, c(FALSE))
  
  #if (!exists("homoscedastic", where = penalty))  penalty$homoscedastic = "FALSE"
  if (!exists("X.dependent.lambda", where = penalty))  penalty$X.dependent.lambda = "FALSE"
  # if (!exists("c", where = penalty) & penalty$homoscedastic!="none") {
  #   penalty$c = 1.1
  # }
  # if (!exists("gamma", where = penalty) & penalty$homoscedastic!="none") {
  #   penalty$gamma = 0.1
  # }
  
  # heteroscedastic and X-independent
  if (penalty$homoscedastic==FALSE && penalty$X.dependent.lambda == FALSE) {
    p <- dim(x)[2]
    n <- dim(x)[1]
    #lambda0 <- 2*penalty$c*sqrt(n)*sqrt(2*log(2*p*log(n)/penalty$gamma))
    lambda0 <- 2 * penalty$c * sqrt(n) * qnorm(1 - penalty$gamma/(2 * p * 1))  # 1=num endogenous variables
    
    #==============================================================================================================
    # The loadings section starts here 
    #==============================================================================================================
    
    # #Ups0 <- 1/sqrt(n) * sqrt(t(t(y^2) %*% (x^2)))
    # load = rep(0,p)
    # for(i in 1:p){
    #   load[i] = sqrt(HAC(matrix(x[,i]*y-mean(x[,i]*y), ncol=1), kernel, bands))
    # }
    Ups0 <- load
    #==============================================================================================================
    # The loadings section ends here  
    #==============================================================================================================
    
    lambda <- lambda0 * Ups0
  }
  
  # heteroscedastic and X-dependent
  if (penalty$homoscedastic==FALSE && penalty$X.dependent.lambda == TRUE) {
    
    #==============================================================================================================
    # The loadings section starts here  # compute lambda
    #==============================================================================================================
    p <- dim(x)[2]
    n <- dim(x)[1]
    Smatrix.boot  = array(0,dim=c(p,nboot))
    cv <- rep(0,nboot)
    # load = rep(0,p)
    
    for (j in 1:nboot){
      res.boot = rnorm(lns) 
      for (i in 1:p){
        Smatrix.boot[i,j] = sum(c(rep(res.boot, each=bns),rep(rnorm(1),n-lns*bns))*x[,i]*y)/sqrt(n)
      }
    }
    
    # # ==============================================================================================================
    # # The HAC section starts here  # compute the penalty loadings
    # # ==============================================================================================================
    # #Ups0 <- 1/sqrt(n) * sqrt(t(t(y^2) %*% (x^2)))
    # 
    # for(i in 1:p){
    #   load[i] = sqrt(HAC(matrix(x[,i]*y-mean(x[,i]*y), ncol=1),kernel, bands))
    # }
    
    for(j in 1:nboot){
      cv[j]=max(abs(Smatrix.boot[,j]/load))
    }
    lambda.boot = 2*quantile(cv, 0.9)*sqrt(n)*penalty$c
    
    
    #==============================================================================================================
    # The loadings section ends here 
    #==============================================================================================================
    Ups0 <- load  
    lambda0 <- lambda.boot
    lambda <- lambda0 * Ups0
    
  }
  
  return(list(lambda0 = lambda0, lambda = lambda, Ups0 = Ups0, method = penalty))
}
