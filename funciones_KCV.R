############################################
# Function to predict m(point) given the 
# sample (y,X) with X in the simplex
# It assumes an homoscedastic model, that is,
# the scale is constant.
# The scale needs to be previously computed 
# and given as sigmahat
#############################################
prediccion.rob.local.poli.homo <-function(point, X, y, H, max.it=20, eta=0.001, type="Tukey", sigmahat=NULL){

	np <- dim(point)[1]
	prediccion <- rep(0,np)
  	 

	if(is.null(sigmahat)){sigmahat <- escala.rob.global(X, y, H)}

	for(k in 1:np){
		punto <- point[k,]
		pepe<- M.local.poli.homo(punto, X, y, H, max.it=max.it, eta=eta, type=type, sigmahat=sigmahat)
		prediccion[k] <- pepe$Mest
	}

	return(list(predic.vector=prediccion))
}


prediccion.cl.local.poli <-function(point, X, y, H){

	np <- dim(point)[1]
	prediccion <- rep(0,np) 

	
	for(k in 1:np){
		punto <- point[k,]
		pepe<- local.poli(punto, X, y, H)
		prediccion[k] <- pepe 
	}

	return(list(predic.vector=prediccion))
}


#################################
# ROBUST KCV
#################################
#################################
# KCV function
#################################

h.rob.KCV <- function(kfold=5, X, y, ventanas,  max.it=20, eta=0.001, type="Tukey", sigmahat=NULL, semilla=1234) {
  # does k-fold CV and returns "robust root mean-squared prediction error"
  	
	set.seed(semilla)

	n <- length(y)
  	k1 <- floor(n/kfold)
  	ids <- rep(1:kfold, each=k1)
  	if( length(ids) < n ) ids <- c(ids, 1:(n%%kfold))

  	ids <- sample(ids)

	pe <- dim(X)[2]
	identidad <- diag(pe-1)

	lv <- length(ventanas)
 	RCV <- rep(0, lv)

	predichos <- matrix(NA,   ncol=n,nrow=lv)
   	for(jv in 1:lv){
		ventana <- ventanas[jv]
		preds <- rep(NA, n)
		
		for(j in 1:kfold) {
    			XX <- X[ids!=j,]
    			yy <- y[ids!=j]
      	    	tmp <- try(prediccion.rob.local.poli.homo(point=X[ids==j,], X=XX, y=yy,  H=ventana*identidad,  max.it=max.it, eta=eta,  type=type, sigmahat=sigmahat) )
    			if( class(tmp) != 'try-error') {
      			preds[ids==j] <-  tmp$predic.vector 
    			}
  		}
  		mun <- median( (preds-y), na.rm=TRUE )
		sn <- mad( (preds-y), na.rm=TRUE )
		RCV[jv] <-sn^2 + mun^2  
		predichos[jv,]=preds
		print(jv)
	}
	 

	h.opt<-ventanas[which.min(RCV)]

	return(list(h.opt=h.opt,RCV=RCV, predichos=predichos))
}
   

#################################
# CLASSICAL KCV
#################################

#################################
# KCV function
#################################


h.cl.KCV<- function(kfold=5, X, y, ventanas, semilla=1234)  {
  # does k-fold CV and returns "robust root mean-squared prediction error"
  	set.seed(semilla)

	n <- length(y)
  	k1 <- floor(n/kfold)
  	ids <- rep(1:kfold, each=k1)
  	if( length(ids) < n ) ids <- c(ids, 1:(n%%kfold))
  	ids <- sample(ids)
  
	pe <- dim(X)[2]
	identidad <- diag(pe-1)

	lv <- length(ventanas)
 	CV <- rep(0, lv)
	predichos <- matrix(NA,   ncol=n,nrow=lv)

   	for(jv in 1:lv){
		ventana <- ventanas[jv]
		preds <- rep(NA, n)

  		for(j in 1:kfold) {
    			XX <- X[ids!=j,]
    			yy <- y[ids!=j]
      	    	tmp <- try(prediccion.cl.local.poli(point=X[ids==j,], X=XX, y=yy,  H=ventana*identidad) )
    			if( class(tmp) != 'try-error') {
      			preds[ids==j] <-  tmp$predic.vector 
    			}
  		}
  		CV[jv] <-mean((preds-y)^2, na.rm=TRUE ) 
		predichos[jv,]=preds
		print(jv)
	}
 
	h.opt<-ventanas[which.min(CV)]

	return(list(h.opt=h.opt,CV=CV, predichos=predichos))
}

 	 

 