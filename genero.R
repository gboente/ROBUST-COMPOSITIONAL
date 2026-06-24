##########################################################
# Computes the inner product between 
# two compositional vectors x and y
##########################################################

prod.interno.atch <- function(x,y){
	d <- length(x)
	producto=NA
	if(length(y)!=d){
		print("ERROR")
	}
	if(length(y)==d){
		if(min(x)!=0 & min(y)!=0){
			matriz.dif<- matrix(0, nrow=d,ncol=d)
			d1 <- d-1
			for(i in 1:d1){
				i1 <- i+1
				for(j in i1:d){
				matriz.dif[i,j] <- log(x[i]/x[j])*log(y[i]/y[j])  

				}
			}
		producto<- sum(matriz.dif)/d 
		}
	}
		return(producto)
}

##########################################################
# b is a d-dimensional vector
# Regression function  
##########################################################
funciong <- function(X,b){
	bstar <- log(b)- mean(log(b))
	V <- log(X)-apply(log(X),1,mean)

	n=dim(V)[1]

	regre1 <- rep(NA, length=n)
	for (i in 1:n){
		vi <- V[i,]
		regre1[i] <- sum(vi*bstar)	
		}
	funciong <- sin(regre1)
	return(funciong)
}


###################################################
# IN this function the inner product is used
# This is the funcion selected in the simulation
###################################################

funciong.bis <- function(X,b){
	
	n=dim(X)[1]

	regre1 <- rep(NA, length=n)
	for (i in 1:n){
		xi <- X[i,]
		regre1[i] <- prod.interno.atch(xi,b)	
		}
	funciong <- sin(regre1)
	return(funciong)
}

####################################################################
## semilla: seed to be used
## n	 : sample size
## d	 : dimension of the covariates
## a	 : parameters of the Dirichlet distribution
## sigma : standard deviation of the errors
## b	 : vector on the simplex to define the regression function
## conta  : list of contamination parameters
## 		type  : 'C0','C1' 
## 		ratio : contamination ratio
##		ee_mu : mean
##		ee_sd : sd
## 		type  : 'CAU'
##	
####################################################################


genero <- function(semilla, n, d, a=c(5, 7, 1, 3, 10, 2, 4), sigma, b,conta, tipo.regre='bis'){

	set.seed(semilla)
	if(length(a)!=d){print('ERROR')}

	if(length(a)==d){
		X <- rdiri( n, a)
 
	}

	epsilon <- sigma*rnorm(n,0,1)

	indices.out <- NULL

 	############################
    	## OUTLIERS in the response    
    	############################

	if(conta$type== 'C1') {
        	out <- rbinom(n, 1, conta$ratio)
        	epsilon[out == 1] <- rnorm(sum(out), conta$ee_mu, conta$ee_sd)
		indices.out <- (1:n)[out == 1]
     	}

	if(conta$type== 'CAU') {
        	epsilon  <- rt(n, df=1)
     	}

	if( tipo.regre!='bis'){
		y <- funciong(X,b) + epsilon
	}
	if( tipo.regre=='bis'){
		y <- funciong.bis(X,b) + epsilon
	}

	return(list(y=y, X=X, iout=indices.out))

}
	
