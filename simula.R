
simula.homo <- function(nrep, n,ntest=500, H,  d=3, a=c(5, 7, 1), sigma, b,conta, tipo.regre='bis'){

	####################################
	# CAMBIO A LA CARPETA DE GUARDADO
	####################################
	carpeta <-paste("SALIDAS_HOMO_a-",a[1],"-",a[2],"-",a[3],"_n",n,"_funciong_",tipo.regre, sep="")    

	if (!dir.exists(carpeta)) dir.create(carpeta)

	setwd(carpeta)

	set.seed(1234)
	x.predict <- rdiri(ntest, a)

	mhat.cl <- matrix(NA, nrow=nrep, ncol=ntest)
	mhat.lineal.cl <- matrix(NA, nrow=nrep, ncol=ntest)
	mhat.rob <- matrix(NA, nrow=nrep, ncol=ntest)
	mhat.lineal.rob <- matrix(NA, nrow=nrep, ncol=ntest)

	media.cl <- media.lineal.cl<- media.rob <- media.lineal.rob <- sigmahat.rob<- rep(NA, length=nrep)

	################################################################
	# PRIMERO CALCULO LA VERDAD
	################################################################
	m.real <- funciong.bis(x.predict,b)
	 

        ARCHIVO.mreal<- paste("mreal_ntest",ntest,"_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  conta$ee_mu,".txt",sep="") 
	
	ARCHIVO.xpred<- paste("xpredict_n",ntest,"_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  conta$ee_mu,".txt",sep="")
	

	if(conta$type== 'C0'){
		ARCHIVO.mreal<- paste("mreal_ntest",ntest,"_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,".txt",sep="") 
		ARCHIVO.xpred<- paste("xpredict_n",ntest,"_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,".txt",sep="")
		}

	if(conta$type== 'CAU'){
		ARCHIVO.mreal<- paste("mreal_ntest",ntest,"_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,".txt",sep="") 
		ARCHIVO.xpred<- paste("xpredict_n",ntest,"_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,".txt",sep="")
		}
	
	write(t(m.real),file=ARCHIVO.mreal,ncolumns=length(m.real),append=T)
	
	for (i in 1:ntest){
 		write(t(x.predict[i,]),file=ARCHIVO.xpred,ncolumns=length(x.predict[i,]),append=T)
	}

	
	################################
	#GUARDO FUNCION ESTIMADA  
	################################

	ARCHIVO_ROB <- paste("M_ROB_HOMO_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  conta$ee_mu,".txt",sep="")   
	ARCHIVO_LINEAL_ROB <-  paste("M_LINEAL_ROB_HOMO_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,"_delta_",100*conta$ratio,
			"_shift_",  conta$ee_mu,".txt",sep="")
	ARCHIVO_CL <- paste("M_CL_HOMO_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  conta$ee_mu,".txt",sep="")   
	ARCHIVO_LINEAL_CL <-  paste("M_LINEAL_CL_HOMO_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  conta$ee_mu,".txt",sep="")
   
	if(conta$type== 'C0'){
		ARCHIVO_ROB <-paste("M_ROB_HOMO_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,".txt",sep="")
		ARCHIVO_LINEAL_ROB <-paste("M_LINEAL_ROB_HOMO_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,".txt",sep="")   
		ARCHIVO_CL  <-paste("M_CL_HOMO_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,".txt",sep="")
		ARCHIVO_LINEAL_CL <-  paste("M_LINEAL_CL_HOMO_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,".txt",sep="")   
  	}

	if(conta$type== 'CAU'){
		ARCHIVO_ROB <-paste("M_ROB_HOMO_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,".txt",sep="")
		ARCHIVO_LINEAL_ROB <-paste("M_LINEAL_ROB_HOMO_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,".txt",sep="")   
		ARCHIVO_CL  <-paste("M_CL_HOMO_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,".txt",sep="")
		ARCHIVO_LINEAL_CL <-  paste("M_LINEAL_CL_HOMO_n",n,"_funciong_",tipo.regre,"_cont_",conta$type,".txt",sep="")   
  	}
 
	t1=Sys.time()
	for(irep in 1:nrep){
		print(c('replicacion= ', irep))

		semilla <- irep+2345

		datos<-genero(semilla, n, d, a, sigma, b,conta, tipo.regre=tipo.regre)

		y <- datos$y
		X <- datos$X
		outliers <- datos$iout

		sigmahat.rob[irep] <- escala.rob.global(X, y, H)

		for (ipred in 1:ntest){

			mhat.cl[irep,ipred ] <- local.constant(x.predict[ipred,], X, y, H)

			mhat.lineal.cl[irep,ipred ] <- local.poli(x.predict[ipred,], X, y, H)

			mhat.rob[irep,ipred ] <- M.local.constant.homo(x.predict[ipred,], X, y, H, sigmahat=sigmahat.rob[irep])$Mest

			mhat.lineal.rob[irep,ipred ] <- M.local.poli.homo(x.predict[ipred,], X, y, H, sigmahat=sigmahat.rob[irep])$Mest
		}
		
 

		media.cl[irep] <- mean((mhat.cl[irep,]-m.real)^2)
		media.lineal.cl[irep] <- mean((mhat.lineal.cl[irep,]-m.real)^2)
		media.rob[irep] <- mean((mhat.rob[irep,]-m.real)^2)
		media.lineal.rob[irep] <- mean((mhat.lineal.rob[irep,]-m.real)^2)

		vector_ROB<- c(irep,media.rob[irep], sigmahat.rob[irep],mhat.rob[irep,]) # OJO GUARDO TAMBIEN EL SIGMAHAT
 		lvec_ROB=length(vector_ROB)
 		
 		vector_CL<- c(irep,media.cl[irep], mhat.cl[irep,])
 		lvec_CL=length(vector_CL)

		vector_CL_lineal<- c(irep,media.lineal.cl[irep], mhat.lineal.cl[irep,])
 		lvec_CL_lineal=length(vector_CL_lineal)

		vector_ROB_lineal<- c(irep,media.lineal.rob[irep], sigmahat.rob[irep], mhat.lineal.rob[irep,]) # OJO GUARDO TAMBIEN EL SIGMAHAT
 		lvec_ROB_lineal=length(vector_ROB_lineal)

 		write(t(vector_ROB),file=ARCHIVO_ROB,ncolumns=lvec_ROB,append=T)
		write(t(vector_ROB_lineal),file=ARCHIVO_LINEAL_ROB,ncolumns=lvec_ROB,append=T)
 
 		write(t(vector_CL),file=ARCHIVO_CL,ncolumns=lvec_CL,append=T)
		write(t(vector_CL_lineal),file=ARCHIVO_LINEAL_CL,ncolumns=lvec_CL_lineal,append=T)
	}

	ecm.cl <- mean(media.cl)

	ecm.lineal.cl <- mean(media.lineal.cl)

	ecm.rob <- mean(media.rob)

	ecm.lineal.rob <- mean(media.lineal.rob)


 	t2=Sys.time()

	tiempo=t2-t1

	print(c('tardo= ', tiempo))

	archivo <- paste("simu-CL-ROB-HOMO-n-",n,"_funciong_",tipo.regre,"_cont_",conta$type,"_delta_",100*conta$ratio,"_shift_",  conta$ee_mu,".RData", sep="")
   
   
	if(conta$type== 'C0'){
		archivo <- paste("simu-CL-ROB-HOMO-n-",n,"_funciong_",tipo.regre,"_cont_",conta$type,".RData", sep="")
	}

        if(conta$type== 'CAU'){
		archivo <- paste("simu-CL-ROB-HOMO-n-",n,"_funciong_",tipo.regre,"_cont_",conta$type,".RData", sep="")
	}

#	save.image(archivo)

	save(x.predict,m.real,mhat.cl,mhat.lineal.cl, mhat.rob, mhat.lineal.rob,
		media.cl, media.lineal.cl, media.rob, media.lineal.rob, sigmahat.rob, 
		ecm.cl, ecm.lineal.cl, ecm.rob, ecm.lineal.rob,
		n, d, a, sigma, b,conta, tipo.regre,file=archivo) 

	return(list(ecm.cl=ecm.cl, ecm.lineal.cl=ecm.lineal.cl, ecm.rob=ecm.rob, ecm.lineal.rob=ecm.lineal.rob))
}