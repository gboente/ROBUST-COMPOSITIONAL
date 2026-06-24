library(compositions)
library(Compositional)
library(mvtnorm)
library(robustbase)


nucleo.chacon <- function(x,xi, H){

	d <- length(x)
	xstar <- ilr(x)
	xistar <- ilr(xi)

	y <- as.vector(xstar-xistar)
	nucleo <- dmvnorm(y,   sigma = H)

	return(nucleo)

	}


nucleo.marzio<- function(x,xi, H){
	determinante <- det(H)
	d <- length(x)
	xstar <- ilr(x)
	xistar <- ilr(xi)
	v <- as.vector(xstar-xistar)

	u <- solve(H)%*%v

	norma <- t(u)%*% u

	nucleo <- exp(- norma/2)/ determinante 

	return(nucleo)

	}

#######################################################
# x vector of dimension d in the simplex
# X  n*d  matrix where X[i,] are in simplex
# y response, length  n
# H smoothing matrix
#######################################################

local.constant <-function(x, X, y, H){

	n<-length(y)
	pesos <-  rep(NA, length=n)

	
	for ( i in 1:n){
		xi <- X[i,]
		pesos[i] <- nucleo.marzio(x,xi, H)
		
	}

	pesos <- pesos/sum(pesos)

	emehat <- sum(pesos*y)

	return(emehat)
}

local.poli <-function(x, X, y, H){

	n<-length(y)
	d <- dim(X)[2]
	pesos <- rep(NA, length=n)
	
	z <- matrix(NA, nrow=n, ncol=(d-1))

	for ( i in 1:n){
		xi <- X[i,]
		pesos[i] <- nucleo.marzio(x,xi, H)
		xstar <- ilr(x)
		xistar <- ilr(xi)
		z[i,] <- as.vector(xstar-xistar)
	}

	matriz.design <- cbind(1,z)
	
	pesos <- pesos/sum(pesos)

	Ka <- diag(pesos)
 
	inversa<- solve(t(matriz.design) %*% Ka %*% matriz.design)
	
	resulta <- inversa%*% t(matriz.design) %*% Ka %*%y

	vectori <- rep(0,length=d)
	vectori[1] <- 1

	emehat<- vectori %*% resulta

	return(emehat)
}



################################################
# Tukey's Psi function
################################################

psi.tukey <- function(r, k=4.685){
  u <- abs(r/k)
  w <- r*((1-u)*(1+u))^2
  w[u>1] <- 0
  return(w)
}

################################################
# Tukey's weight function "Psi(r)/r"
################################################
w.tukey  <- function(r, k= 4.685){
  u <- abs(r/k)
  w <- ((1 + u) * (1 - u))^2
  w[u > 1] <- 0
  return(w)
}

################################################
# Huber's Psi function
################################################
psi.huber <- function(r, k=1.345)
  pmin(k, pmax(-k, r))

################################################
# Huber's weight function "Psi(r)/r"
################################################

w.huber  <- function(r, k=1.345)
  pmin(1, k/abs(r))


 
##################################################
# To compute the local median
##################################################

weighted.fractile<-function (y, w, p){
    w <- w/sum(w)
    a <- 1 - p
    b <- p
    ox <- order(y)
    y <- y[ox]
    w <- w[ox]
    k <- 1
    low <- cumsum(c(0, w))
    up <- sum(w) - low
    df <- a * low - b * up
    repeat {
        if (df[k] < 0) 
            k <- k + 1
        else if (df[k] == 0) 
            return((w[k] * y[k] + w[k - 1] * y[k - 1])/(w[k] + 
                w[k - 1]))
        else return(y[k - 1])
    }
}

########################################
# LOCAL MEDIAN
#########################################
mediana.local <- function(x, X, y, H){

	n<-length(y)
	pesos <- rep(NA, length=n)

	for ( i in 1:n){
		xi <- X[i,]
		pesos[i] <- nucleo.marzio(x,xi, H)
		
	}
 
	emehat <- unname(weighted.fractile(y,pesos,0.5) )

	return(emehat)
}

########################################
# GLOBAL SCALE ESTIMATOR
#########################################

Rho <- function(r, tuning.rho=1.54764) Mchi(x=r, cc=tuning.rho, psi="bisquare", deriv=0)
Rhop <- function(r, tuning.rho=1.54764) Mchi(x=r, cc=tuning.rho,  psi="bisquare", deriv=1)

s.scale <- function(r, tuning.rho=1.54764, b=.5, max.it=1000, ep=1e-4){
    s1 <- mad(r)
    if(abs(s1)<1e-10) return(s1)
    s0 <- s1 + 1
    it <- 0
    while( ( abs(s0-s1) > ep ) && (it < max.it) ) {
    it <- it + 1
    s0 <- s1
    s1 <- s0*mean(Rho(r/s0,tuning.rho=tuning.rho))/b
}
return(s1)
}

########################################
# LOCAL SCALE ESTIMATOR
# LOCAL MAD
#########################################

escala.rob <- function(x, X, y, H){

	n<-length(y)
	 
	emehat <- mediana.local(x, X, y, H)

	residuo <- abs(y-emehat)

	escala <-  mediana.local(x, X, residuo, H)
	return(escala)
}


########################################
# ESTIMADOR DE ESCALA GLOBAL
# MODELO HOMOCEDASTICO
# USO UN S-estimador de los residuos
#########################################

escala.rob.global <- function(X, y, H){

	n<-length(y)
	
	residuo <- rep(NA, length=n)
	for (i in 1:n){
		x<- X[i,]
		emehat <- mediana.local(x, X, y, H)
		residuo[i] <- y[i]-emehat
	}

	escala <-  s.scale(residuo)
	return(escala)
}


########################################
# LOCAL M-SMOOTHER 
#########################################
 
################################################
# The function M.local.constant.homo
# computes the local M-smoother
# assuming an homoscedastic model, that is,
# the scale is constant.
# The scale needs to be previously computed 
# and given as sigmahat
################################################

M.local.constant.homo <-function(x, X, y, H, max.it=20, eta=0.001, type="Tukey", sigmahat=NULL){

	if(is.null(sigmahat)){sigmahat <- escala.rob.global(X, y, H)}
	
	n<-length(y)
	pesos <- rep(NA, length=n)

	for ( i in 1:n){
		xi <- X[i,]
		pesos[i] <- nucleo.marzio(x,xi, H)
		
	}

	pesos <- pesos/sum(pesos)

	muhat<- emehat.ini <- mediana.local(x, X, y, H)
	
	corte = 10 * eta
	it<- 1
	while( (corte > eta) && (it<  max.it)){

		res <- aux1 <- aux2 <- rep(NA,length=n)

		for(i in 1:n){
			res[i] = ( y[i] - muhat) / sigmahat 
			aux1[i] = pesos[i] * w.tukey(res[i], k= 4.685) * y[i] 
			aux2[i] = pesos[i] * w.tukey(res[i], k= 4.685) 
		} 
		muold = muhat 
		muhat = sum(aux1) / sum(aux2)
		corte =  abs(muold-muhat) / (abs(muold) + eta ) 
		it = it + 1 
	}

	estimador <- muhat

	return(list(Mest=estimador, mediana=emehat.ini, sigmahat=sigmahat))

}

 
################################################
# The function M.local.poli.homo
# computes the local linear M-smoother
# assuming an homoscedastic model, that is,
# the scale is constant.
# The scale needs to be previously computed 
# and given as sigmahat
################################################

M.local.poli.homo <-function(x, X, y, H, max.it=20, eta=0.001, type="Tukey", sigmahat=NULL){

	if(is.null(sigmahat)){sigmahat <- escala.rob.global(X, y, H)}

	n<-length(y)
	d <- dim(X)[2]

	vectori <- rep(0,length=d)
	vectori[1] <- 1

	pesos <- rep(NA, length=n)
	
	z <- matrix(NA, nrow=n, ncol=(d-1))

	for ( i in 1:n){
		xi <- X[i,]
		pesos[i] <- nucleo.marzio(x,xi, H)
		xstar <- ilr(x)
		xistar <- ilr(xi)
		z[i,] <- as.vector(xstar-xistar)
	}

	matriz.design <- cbind(1,z)
	
	pesos <- pesos/sum(pesos)


	emehat.ini <- mediana.local(x, X, y, H)
	b0hat<- emehat.ini
	b1hat <- rep(0, length=(d-1))
	bhat<- as.vector(c(b0hat, b1hat))
	
	corte = 10 * eta
	it<- 1
	while( (corte > eta) && (it<  max.it)){

		res <- auxiliar <- rep(NA,length=n)

		for(i in 1:n){
			res[i] = ( y[i] - b0hat - b1hat%*%  z[i,]) / sigmahat 
			auxiliar[i] = pesos[i] * w.tukey(res[i], k= 4.685)  
		} 

		Ka <- diag(auxiliar)
		
 
		inversa<- solve(t(matriz.design) %*% Ka %*% matriz.design)
	
		bhat <- inversa%*% t(matriz.design) %*% Ka %*%y

		 

		b.old <- bhat
		b0hat <- bhat[1]
		b1hat <- bhat[-1]  
		corte =  sqrt(sum((b.old -bhat)^2))/ (sqrt(sum(b.old^2)) + eta ) 
		it = it + 1 
	}

	emehat<- vectori %*% bhat

	
	return(list(Mest=emehat, inicial=emehat.ini, sigmahat=sigmahat))

}
