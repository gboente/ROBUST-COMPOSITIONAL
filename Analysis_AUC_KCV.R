rm(list=ls())

library(UsingR)
library(robustbase)   

#####################################
# LIBRARIES FOR COMPOSITIONAL DATA
# AND TO PLOT THE TERNARY DIAGRAM
#####################################
           
library("Ternary")
library(ggtern)

library(compositions)
library(Compositional)

#####################################
# LIBRARIES TO PLOT THE SURFACES
##################################### 
library(lattice)
library(animation)

library(plot3D)
require("graphics")


library('rgl')
library('magick')


##############################################################
# CODES TO COMPUTE THE Kernel estimators
##############################################################
source('funciones.R')
source('funciones_KCV.R')



datos.gude<- read.table("data_example.txt",header=T)
names(datos.gude)
 
 
XCOMP<- cbind(datos.gude$PROTEINS,datos.gude$LIPIDS,datos.gude$CARBOH)[datos.gude$dm==0,]
y <- datos.gude$AUC[datos.gude$dm==0]
x <- ilr(XCOMP)

#############################################
# PLOTS OF THE COMPOSITIONAL COVARIATES
# in the Ternary Diagram
#############################################
 
pdf('ternary-plot-covariates.pdf',bg='transparent')
 

par(mar = rep(0.2, 4))

TernaryPlot(alab = "Proteins", blab = "Lipids",
 clab ="Carbohydrates" ,  
atip.pos=4, atip.rotate=0, 
cex=1.8, lab.col=c("blue4","red4","green4"),grid.minor.lines=0)

AddToTernary(graphics::points, XCOMP, pch = 16, cex = 1.2,,col='black')

dev.off()
 

minimo1=min(x[,1])
maximo1=max(x[,1])


minimo2=min(x[,2])
maximo2=max(x[,2])


#############################################
# PLOTS OF THE COMPOSITIONAL COVARIATES
# in the ilr space
#############################################

pdf('ilr-plot-covariates.pdf',bg='transparent')

par(mar=c(5,6,3,3))
 
plot(x[,1], x[,2], 
	xlab=expression(x[1]^"*" ), ylab=expression(x[2]^"*"),
	xlim=c(minimo1,maximo1), ylim=c(minimo2,maximo2),
	col = "blue4", pch=16, cex=1.2,cex.axis=1.5,
	cex.lab=1.8, mgp=c(3.5,1,0), las=1.2) 

dev.off()

  

##############################################################
# Kernel estimators
##############################################################
 


######################################
# WE CHOOSE THE BANDWIDTH USING KCV
######################################

#####################################
# CLASSICAL ESTIMATOR
# THE GRID VARIES BETWEEN 0.1 AND 1
#####################################

ventanas <- seq(0.1,1,by=0.1)

 
t1=Sys.time()

opt.cl <- h.cl.KCV(kfold=5, X=XCOMP, y=y, ventanas=ventanas,semilla=1234)

t2=Sys.time()
t2-t1

h.opt.cl <- opt.cl$h.opt
CV.cl <- opt.cl$CV

h.opt.cl
# 0.4

resultado.1<- cbind(ventanas, CV.cl)

resultado.1
 

nombre='CV.pdf'
pdf(nombre,bg='transparent')

plot(ventanas, CV.cl)

dev.off()


##############################################
# WE REFINE THE GRID AROUND THE MINIMUM
##############################################

ventanas <- seq(0.2,0.6,by=0.02)

 
 t1=Sys.time()

opt.cl <- h.cl.KCV(kfold=5, X=XCOMP, y=y, ventanas=ventanas,semilla=1234)

t2=Sys.time()
t2-t1

h.opt.cl <- opt.cl$h.opt
CV.cl <- opt.cl$CV

h.opt.cl
#0.44
 
 

#############################################
# ROBUST CROSS VALIDATION
#############################################
# WE COMPUTE THE SCALE ESTIMATOR
# USING A GIVEN BANDWIDTH H=0.2
#############################################

H= 0.2*diag(2)

X=XCOMP
sigmahat.rob  <- escala.rob.global(X, y, H)

sigmahat.rob 

#############################################
# WE SELECT THE BANDWIDTH FOR 
# THE REGRESSION ESTIMATOR USING 
# ROBUST CROSS VALIDATION
# THE GRID VARIES BETWEEN 0.1 AND 1
#############################################

ventanas <- seq(0.1,1,by=0.1)


set.seed(1234)
 t1=Sys.time()


opt.rob <- h.rob.KCV(kfold=5, X=XCOMP, y=y, ventanas=ventanas,  max.it=20, eta=0.001, type="Tukey", sigmahat=sigmahat.rob,semilla=1234 )


t2=Sys.time()
t2-t1

h.opt.rob <- opt.rob$h.opt

CV.rob <-  opt.rob$RCV
h.opt.rob
#0.2

 
##############################################
# WE REFINE THE GRID AROUND THE MINIMUM
# WE GO BEYOND 1 TO SEE IF THERE IS A LOCAL
# MINIMUM LARGER THAN 1
##############################################


ventanas <- seq(0.18,1.4,by=0.02)

 
 t1=Sys.time()

opt.rob <- h.rob.KCV(kfold=5, X=XCOMP, y=y, ventanas=ventanas,  max.it=20, eta=0.001, type="Tukey", sigmahat=sigmahat.rob ,semilla=1234)


t2=Sys.time()
t2-t1
 

h.opt.rob <- opt.rob$h.opt

CV.rob <-  opt.rob$RCV
h.opt.rob
#0.22


########################################
# DEFINE THE BANDWIDTH MATRIX
########################################

H.rob= h.opt.rob*diag(2)

H.cl= h.opt.cl*diag(2)

#######################################
#  COMPUTE THE ESTIMATORS AT EACH X[i,]
#######################################

X=XCOMP
x.predict <- X
 
ntest <- dim(X)[1]

mhat.cl <- rep(NA, length=ntest)
mhat.lineal.cl <- rep(NA, length=ntest)

mhat.rob <- mhat.lineal.rob <- rep(NA, length=ntest)

sigmahat.rob  <- escala.rob.global(X, y, H.rob)

sigmahat.rob 

for (ipred in 1:ntest){
	mhat.lineal.cl[ipred ] <- local.poli(x.predict[ipred,], X, y, H.cl)
	mhat.lineal.rob[ipred ] <- M.local.poli.homo(x.predict[ipred,], X, y, H.rob, sigmahat=sigmahat.rob)$Mest
		}


#######################################
# OUTLIER IDENTIFICATION
# THROUGH LARGE RESIDUALS VIA BOXPLOT
#######################################
names(y)=1:length(y)
residuos.lineal.local.rob =y-mhat.lineal.rob
atipicos.lineal.local.rob =as.numeric(names(boxplot(residuos.lineal.local.rob)$out))
atipicos.lineal.local.rob
#  13  23  33  35  39  71  72  98 107 112 135 164 175 245 357 413 478

#######################################
# LABELS IN THE ORIGINAL DATA BASE
#######################################

id.new <- datos.gude$id[datos.gude$dm==0 ]

identifico.out<-  id.new[atipicos.lineal.local.rob ]
identifico.out
# 23   41   58   65   71  133  134  190  208  227  277  358  391  611 1008 1178 1382  


###############################################
# WE REMOVE THE DETECTED OUTLIERS 
# TO COMPUTE THE CLASSICAL ESTIMATORS
# WITHOUT OUTLIERS
##############################################
 
XCOMP.so  <- XCOMP[-atipicos.lineal.local.rob,]

y.so <- y[-atipicos.lineal.local.rob]

######################################
# WE CHOOSE THE BADWIDTH USING 
# CLASSICAL KCV
######################################

ventanas <- seq(0.1,1,by=0.1)

  
 t1=Sys.time()

opt.cl.so <- h.cl.KCV(kfold=5,X=XCOMP.so, y=y.so, ventanas=ventanas ,semilla=1234)


t2=Sys.time()
t2-t1

h.opt.cl.so <- opt.cl.so$h.opt
CV.cl.so <- opt.cl.so$CV

h.opt.cl.so
# 0.3   
  

######################################
# WE REFINE THE GRID
######################################

hmin=max(0.1, h.opt.cl.so - 0.1)
hmax=min(1, h.opt.cl.so + 0.1)

ventanas <- seq(hmin,hmax,by=0.02)
 

 t1=Sys.time()

opt.cl.so <- h.cl.KCV(kfold=5,X=XCOMP.so, y=y.so, ventanas=ventanas ,semilla=1234)


t2=Sys.time()
t2-t1

h.opt.cl.so <- opt.cl.so$h.opt
CV.cl.so <- opt.cl.so$CV

h.opt.cl.so
#0.26
 
########################################
# DEFINE THE BANDWIDTH MATRIX
# AND COMPUTE THE LOCAL LINEAR 
# ESTIMATORS AT X.so[i,] 
########################################

H.cl.so= h.opt.cl.so*diag(2)


X.so=XCOMP.so
x.predict.so <- X.so


ntest.so <- dim(X.so)[1]


mhat.lineal.cl.so <- rep(NA, length=ntest.so)
 


for (ipred in 1:ntest.so){
	mhat.lineal.cl.so[ipred ] <- local.poli(x.predict.so[ipred,], X.so, y.so, H.cl.so)
	}

 

#######################################
# SOME PLOTS
#######################################

nombre='boxplot-resid.pdf'
 
pdf(nombre,bg='transparent')

boxplot(residuos.lineal.local.rob, col='blue')

dev.off()


nombre='rob-resid-vs-predicted.pdf'
pdf(nombre,bg='transparent')

par(mar=c(5,5,3,3))
plot(mhat.lineal.rob, residuos.lineal.local.rob, 
   xlab=expression(hat(y)), ylab=expression(hat(r)),pch=16, lwd=2,cex=1.3, cex.lab=1.3)

dev.off()
 

 
################################################################
# COMPUTE:
# THE CLASSICAL ESTIMATORS  
# THE ROBUST ESTIMATORS  
# THE CLASSICAL ESTIMATORS WITH THE DATA SET WITHOUT OUTLIERS
# All AT THE DEFINED GRID
# WE ALSO SAVE THE GRID OF COMPOSITIONAL POINTS CREATED 
#################################################################

 
largo=50

ejex=seq(0,1,length=largo)
ejey=seq(-0.3,1,length=largo)

 
mhat.lineal.cl.ilr <- matrix(NA, ncol=largo,nrow=largo)

mhat.lineal.rob.ilr <- matrix(NA, ncol=largo,nrow=largo)

mhat.lineal.cl.ilr.so  <- matrix(NA, ncol=largo,nrow=largo)

data.grilla <- matrix(0,nrow=length(ejex)*length(ejey),ncol=3)

sigmahat.rob.ilr  <- sigmahat.rob
 

nx <- 0;
 
for(i in 1:length(ejex)){
	for(j in 1:length(ejey)){
		punto<- c(ejex[i],ejey[j])
		X.pred <- ilrInv(punto)
		data.grilla [j+nx,1:3] <-  X.pred 

 		mhat.lineal.cl.ilr[i,j] <- local.poli(X.pred, X, y, H.cl)

 		mhat.lineal.rob.ilr[i,j] <- M.local.poli.homo(X.pred, X, y, H.rob, sigmahat=sigmahat.rob)$Mest
					
		mhat.lineal.cl.ilr.so[i,j] <- local.poli(X.pred, X.so, y.so, H.cl.so)
	}
	nx <- nx+length(ejey);
}

####################################
# SOME PLOTS
####################################
#######################################################
# PLOTTING THE GRID AND THE DATA COVARIATES
# IN THE TERNARY
########################################################

nombre= 'ternary_data_grid.pdf' 
 

pdf(nombre,bg='transparent')

par(mar = rep(0.2, 4))
TernaryPlot(alab = "Proteins", blab = "Lipids",
 clab ="Carbohydrats" ,  cex=1.8, lab.col=c("blue4","red4","green4"))

AddToTernary(graphics::points, XCOMP, pch = 16, cex = 1.8)

AddToTernary(graphics::points,  data.grilla, col="gold", pch = 1, cex = 0.2)

dev.off()


#######################################################
# PLOTTING THE GRID AND THE DATA COVARIATES
# IN THE ILR SPACE
########################################################

matriz<- c()
for(i in 1:length(ejex)){
  for(j in 1:length(ejey)){
    matriz<- rbind(matriz,c(ejex[i],ejey[j]))
    }
  }
 


pdf('ilr_data_grid.pdf',bg='transparent')

par(mar=c(5,6,3,3))
plot(x[,1], x[,2], 
	xlab=expression(x[1]^"*" ), ylab=expression(x[2]^"*"),
	xlim=c(minimo1,maximo1), ylim=c(minimo2,maximo2),
	col = "blue4", pch=16, cex=1.8,cex.axis=1.5,
	cex.lab=1.8, mgp=c(3.5,1,0), las=1.2) 

points(matriz[,1] ,matriz[,2], col="gold",pch=1, cex=0.8)

dev.off()


##############################################
# PLOTS OF THE SURFACES IN THE TERNARY
##############################################
 

######################################################################
# THIS FUNCTION IS AS movie3d BUT HELPS IN LATEX 
# TO SAVE EACH FRAME AS A .PNG FILE
######################################################################

 peli<-   function (f, duration, dev = cur3d(), ..., fps = 10, movie = "movie", 
    frames = movie, dir =  getwd(), convert = NULL, clean = TRUE, 
    verbose = TRUE, top = !rgl.useNULL(), type = "gif", startTime = 0, 
    webshot = TRUE) 
{
    olddir <- setwd(dir)
    on.exit(setwd(olddir))
    for (i in round(startTime * fps):(duration * fps)) {
        time <- i/fps
        if (cur3d() != dev) 
            set3d(dev)
        stopifnot(cur3d() != 0)
        args <- f(time, ...)
        subs <- args$subscene
        if (is.null(subs)) 
            subs <- currentSubscene3d(dev)
        else args$subscene <- NULL
        for (s in subs) par3d(args, subscene = s)
        filename <- paste(frames,i, '.png', sep='')
        if (verbose) {
            cat(gettextf("Writing '%s'\r", filename))
            flush.console()
        }
        if (top) 
            rgl.bringtotop()
        snapshot3d(filename = filename, webshot = webshot)
    }
    cat("\n")
    if (.Platform$OS.type == "windows") 
             system <- shell
    if (is.null(convert) && requireNamespace("magick", quietly = TRUE)) {
        m <- NULL
        for (i in round(startTime * fps):(duration * fps)) {
            filename <- paste(frames,i, '.png', sep='')
            frame <- magick::image_read(filename)
            if (is.null(m)) 
                m <- frame
            else m <- c(m, frame)
            if (clean) 
                unlink(filename)
        }
        m <- magick::image_animate(m, fps = fps, loop = 1, dispose = "previous")
        magick::image_write(m, paste0(movie, ".", type))
        return(invisible(m))
    }
    else if (is.null(convert)) {
        warning("R package 'magick' is not installed; trying external package.")
        convert <- TRUE
    }
    if (is.logical(convert) && convert) {
        progname <- "magick"
        version <- try(system2(progname, "--version", stdout = TRUE, 
            stderr = TRUE), silent = TRUE)
        if (inherits(version, "try-error") || !length(grep("ImageMagick",
version))) {
            progname <- "convert"
            version <- try(system2(progname, "--version", stdout = TRUE, 
                stderr = TRUE), silent = TRUE)
        }
        if (inherits(version, "try-error") || !length(grep("ImageMagick", 
            version))) 
            stop("'ImageMagick' not found")
        filename <- paste0(movie, ".", type)
        if (verbose) 
            cat(gettextf("Will create: %s\n", file.path(dir, 
                filename)))
        convert <- paste(progname, "-delay 1x%d %s*.png %s.%s")
    }
    if (is.character(convert)) {
        convert <- sprintf(convert, fps, frames, movie, type, 
            duration, dir)
        if (verbose) {
            cat(gettextf("Executing: '%s'\n", convert))
            flush.console()
        }
        system(convert)
        if (clean) {
            if (verbose) 
                cat(gettext("Deleting frames\n"))
            for (i in 0:(duration * fps)) {
                filename <- sprintf("%s%03d.png", frames, i)
                unlink(filename)
            }
        }
    }
    invisible(convert)
} 


#################################################
# DEFINE SOME COORDINATES FOR THE PLOTS
#################################################

coord.tern.grillax <- matrix(0,nrow=length(ejex) ,ncol=length(ejey)) 
coord.tern.grillay <- matrix(0,nrow=length(ejex) ,ncol=length(ejey)) 

coordx <- coordy<- rep(NA, length=length(ejex)*length(ejey))

coord.compo.grillax <- coord.compo.grillay<- coord.compo.grillaz<- matrix(NA,nrow=length(ejex),ncol=length(ejey))
 

nx <- 0;
j <- 1;
while(j<=length(ejex)){
	for (k in 1:length(ejey)){
		punto <- c(ejex[j], ejey[k]);
		X.pred <- as.vector( unname(ilrInv(punto)))
		pepe <- as.vector(CoordinatesToXY(X.pred))
		coord.tern.grillax[j,k] <- pepe[1]
		coord.tern.grillay[j,k] <- pepe[2]

		coordx[k+nx]<- pepe[1]
		coordy[k+nx]<-pepe[2]

		coord.compo.grillax[j,k] <- X.pred[1] # PROTEIN 
		coord.compo.grillay[j,k] <- X.pred[2] # LIPIDS
		coord.compo.grillaz[j,k] <- X.pred[3] #CARBOHYDRATS
	}
	nx <- nx+length(ejey);
   	j <- j+1;
}



ejex_triang<- seq(-4,4,length=50)
ejey_triang<- seq(-4,4,length=50)

coordx_triang <- coordy_triang<- rep(NA, length=length(ejex)*length(ejey))

nx <- 0;
j <- 1;

while(j<=length(ejex_triang)){
	for (k in 1:length(ejey_triang)){
 		punto <- c(ejex_triang[j], ejey_triang[k]);
		X.pred <- as.vector( unname(ilrInv(punto)))
		pepe <- as.vector(CoordinatesToXY(X.pred))
 
		coordx_triang[k+nx]<- pepe[1]
		coordy_triang[k+nx]<-pepe[2]
	}
	nx <- nx+length(ejey_triang);
	j <- j+1;
}
 
 
  
###################################
# USES library('rgl')
##################################

nrz <- nrow(mhat.lineal.cl.ilr)
ncz <- ncol(mhat.lineal.cl.ilr)

ang.theta <- 20
ang.phi <- 50

# Create a function interpolating colors in the range of specified colors

###################################
# COLORS FOR CLASSICAL ESTIMATOR
###################################

jet.colors <- colorRampPalette( c("red", "red4") )

# Generate the desired number of colors from this palette

nbcol <- 100
color <- jet.colors(nbcol)

# Compute the z-value at the facet centres

zfacet <- (mhat.lineal.cl.ilr[-1, -1] +mhat.lineal.cl.ilr[-1, -ncz] 
	+ mhat.lineal.cl.ilr[-nrz, -1] + mhat.lineal.cl.ilr[-nrz, -ncz])/4

# Recode facet z-values into color indices

facetcol <- cut(zfacet, nbcol)
 


###################################
# COLORS FOR ROBUST ESTIMATOR
###################################

jet.colors.rob <- colorRampPalette( c("blue", "green") )

color.rob <- jet.colors.rob(nbcol)

# Compute the z-value at the facet centres

zfacet.rob <- (mhat.lineal.rob.ilr[-1, -1] +mhat.lineal.rob.ilr[-1, -ncz] 
	+ mhat.lineal.rob.ilr[-nrz, -1] + mhat.lineal.rob.ilr[-nrz, -ncz])/4

# Recode facet z-values into color indices

facetcol.rob <- cut(zfacet.rob, nbcol)
  
#############################################
# PLOTTING CLASSICAL AND ROBUST SURFACES
# IN THE TERNARY
############################################
 
par3d(windowRect = c(15, 15, 500, 500))
    
  
persp3d(coord.tern.grillax,coord.tern.grillay,mhat.lineal.cl.ilr,
	theta = ang.theta, phi = ang.phi,  
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2,  col = color[facetcol],
	#screen = list(z = 20, x = -70),
	alpha=0.6,
	#xlim=c(-0.6,0.6),ylim=c(0,1.2),
	zlim=c(80,105),xlim=c(-0.5,0.5), ylim=c(0,1),
	axes=FALSE,  
	breaks = seq(80,105, length=11),
	#colkey = list(side = 2, length = 0.5,dist=-0.1,shift=0.17),
	facets=TRUE,   plot=T) 

persp3d(coord.tern.grillax,coord.tern.grillay,mhat.lineal.cl.ilr,
	theta = ang.theta, phi = ang.phi,  col="red4",
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2, front = "lines", back = "lines",
	alpha=0.4,
 	lit = FALSE, add = TRUE)



persp3d(coord.tern.grillax,coord.tern.grillay,mhat.lineal.rob.ilr,
	theta = ang.theta, phi = ang.phi,  
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2,  col = color.rob[facetcol.rob],
	#screen = list(z = 20, x = -70),
	alpha=0.6,
	#xlim=c(-0.6,0.6),ylim=c(0,1.2),
	zlim=c(80,105),xlim=c(-0.5,0.5), ylim=c(0,1),
	axes=FALSE,  
	breaks = seq(80,105, length=11),
	#colkey = list(side = 2, length = 0.5,dist=-0.1,shift=0.17),
	facets=TRUE,   plot=T) 

persp3d(coord.tern.grillax,coord.tern.grillay,mhat.lineal.rob.ilr,
	theta = ang.theta, phi = ang.phi,  col="blue4",
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2, front = "lines", back = "lines",
	alpha=0.4,
 	lit = FALSE, add = TRUE)


  
zeta= rep(80, length=length(coordy))
points3d(coordx,coordy, zeta,  col="gray60", pch=".",cex=0.5, add=T)
points3d(coordx_triang,coordy_triang, zeta,   col="gray", pch=".",cex=0.5, add=T)

   

 axes3d(c( 'z'))

arrow3d( c(0.5,-0.05, 80),c(-0.5,-0.05, 80),n=2.1,width=0.2, s=0.2,  type="rotation",col="black",add=T)

arrow3d( c(-0.55,0, 80), c(-0.05,0.9, 80),  n=2.1,width=0.2,s=0.2,  type="rotation",col="black",add=T)

arrow3d( c(0.05,0.9, 80), c(0.55,0, 80), n=2.1,width=0.2, s=0.2, type="rotation",col="black",add=T)

text3d(0,0,78, "Carbohydrates",add=T)
text3d(-0.3,0.5,80, "Proteins",add=T,adj=c(1,0))
text3d(0.6,0.4,80, "Lipids",add=T,adj=c(1,0) )


par3d( zoom=0.67 )
 
rgl.snapshot('AUC_cl_rob_KCV_3dplot-bis.png', fmt = 'png')
   


#############################################
# PLOTTING: ROBUST SURFACES
# AND CLASSICAL WITHOUT OUTLIERS
# IN THE TERNARY
############################################

###################################
# COLORS FOR CLASSICAL ESTIMATOR
# COMPUTED WITHOUT THE OUTLIERS
###################################

jet.colors.ls.so <- colorRampPalette( c("yellow","gold") )

#col2rgb(paste0("gold", 1:3)) #colorRampPalette( c("gold", "yellow") )

# Compute the z-value at the facet centres

zfacet.ls.so <- (mhat.lineal.cl.ilr.so[-1, -1] +mhat.lineal.cl.ilr.so[-1, -ncz] 
	+ mhat.lineal.cl.ilr.so[-nrz, -1] + mhat.lineal.cl.ilr.so[-nrz, -ncz])/4

# Recode facet z-values into color indices

facetcol.ls.so <- cut(zfacet.ls.so, nbcol)
  

######################################
# PLOTS OF THE SURFACES:
# ONLY THE ROBUST ONE AND
# THE CLASSICAL WITHOUT OUTLIERS
##########################################

par3d(windowRect = c(15, 15, 500, 500))


persp3d(coord.tern.grillax,coord.tern.grillay,mhat.lineal.rob.ilr,
	theta = ang.theta, phi = ang.phi,  
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2,  col = color.rob[facetcol.rob],
	#screen = list(z = 20, x = -70),
	alpha=0.6,
	#xlim=c(-0.6,0.6),ylim=c(0,1.2),
	zlim=c(80,105),xlim=c(-0.5,0.5), ylim=c(0,1),
	axes=FALSE,  
	breaks = seq(80,105, length=11),
	#colkey = list(side = 2, length = 0.5,dist=-0.1,shift=0.17),
	facets=TRUE,   plot=T) 

persp3d(coord.tern.grillax,coord.tern.grillay,mhat.lineal.rob.ilr,
	theta = ang.theta, phi = ang.phi,  col="blue4",
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2, front = "lines", back = "lines",
	alpha=0.4,
 	lit = FALSE, add = TRUE)

persp3d(coord.tern.grillax,coord.tern.grillay,mhat.lineal.cl.ilr.so,
	theta = ang.theta, phi = ang.phi,  
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2,  col = color[facetcol.ls.so],
	#screen = list(z = 20, x = -70),
	alpha=0.6,
	#xlim=c(-0.6,0.6),ylim=c(0,1.2),
	zlim=c(80,105),xlim=c(-0.5,0.5), ylim=c(0,1),
	axes=FALSE,  
	breaks = seq(80,105, length=11),
	#colkey = list(side = 2, length = 0.5,dist=-0.1,shift=0.17),
	facets=TRUE,add = TRUE,   plot=T) 

persp3d(coord.tern.grillax,coord.tern.grillay,mhat.lineal.cl.ilr.so,
	theta = ang.theta, phi = ang.phi,  col="gold3",
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2, front = "lines", back = "lines",
	alpha=0.4,
 	lit = FALSE, add = TRUE)


zeta= rep(80, length=length(coordy))
points3d(coordx,coordy, zeta,  col="gray60", pch=".",cex=0.5, add=T)
points3d(coordx_triang,coordy_triang, zeta,   col="gray", pch=".",cex=0.5, add=T)


axes3d(c( 'z'))

arrow3d( c(0.5,-0.05, 80),c(-0.5,-0.05, 80),n=2.1,width=0.2, s=0.2,  type="rotation",col="black",add=T)

arrow3d( c(-0.55,0, 80), c(-0.05,0.9, 80),  n=2.1,width=0.2,s=0.2,  type="rotation",col="black",add=T)

arrow3d( c(0.05,0.9, 80), c(0.55,0, 80), n=2.1,width=0.2, s=0.2, type="rotation",col="black",add=T)


text3d(0,0,78, "Carbohydrates",add=T)
text3d(-0.3,0.5,80, "Proteins",add=T,adj=c(1,0))
text3d(0.6,0.4,80, "Lipids",add=T,adj=c(1,0) )

par3d( zoom=0.67 )
 
rgl.snapshot('AUC_rob_clSO_KCV_3dplot-bis.png', fmt = 'png')
 

####################################################
# PLOT OF THE ESTIMATORS VERSUS TWO OF THE
# COMPOSITIONAL COMPONENTS
####################################################

nrz <- nrow(mhat.lineal.cl.ilr)
ncz <- ncol(mhat.lineal.cl.ilr)

# Create a function interpolating colors in the range of specified colors

################################################################
# COLORS FOR THE CLASSICAL ESTIMATOR
################################################################

jet.colors <- colorRampPalette( c("red", "red4") )

# Generate the desired number of colors from this palette

nbcol <- 100
color <- jet.colors(nbcol)

# Compute the z-value at the facet centres

zfacet <- (mhat.lineal.cl.ilr[-1, -1] +mhat.lineal.cl.ilr[-1, -ncz] 
	+ mhat.lineal.cl.ilr[-nrz, -1] + mhat.lineal.cl.ilr[-nrz, -ncz])/4

# Recode facet z-values into color indices

facetcol <- cut(zfacet, nbcol)
 

################################################################
# COLORS FOR THE ROBUST ESTIMATOR
################################################################

jet.colors.rob <- colorRampPalette( c("blue", "green") )

color.rob <- jet.colors.rob(nbcol)

# Compute the z-value at the facet centres

zfacet.rob <- (mhat.lineal.rob.ilr[-1, -1] +mhat.lineal.rob.ilr[-1, -ncz] 
	+ mhat.lineal.rob.ilr[-nrz, -1] + mhat.lineal.rob.ilr[-nrz, -ncz])/4

# Recode facet z-values into color indices

facetcol.rob <- cut(zfacet.rob, nbcol)


################################################################
# CLASSICAL AND ROBUST ESTIMATORS   
# VERSUS PROTEINS AND LIPIDS
# IN THE PAPER WE CHOOSE THE FRAME 82 FOR BETTER VISUALIZATION
################################################################  

par3d(windowRect = c(15, 15, 500, 500))

persp3d(coord.compo.grillax,coord.compo.grillay,mhat.lineal.cl.ilr,
	theta = ang.theta, phi = ang.phi,  
	ltheta = 120,   ticktype = "detailed", 
	xlab=" P ",ylab="L", zlab = "", 
	expand = 2,  
      col= "red", 
	  aspect=c(1,1), # y-size/x-size and z-size/x-size
      alpha=0.6,  
	 zlim=c(80,105), 
	 breaks = seq(80,105, length=11),
	#colkey = list(side = 2, length = 0.5,dist=-0.1,shift=0.17),
	facets=TRUE,   plot=T) 

persp3d(coord.compo.grillax,coord.compo.grillay,mhat.lineal.cl.ilr,
	theta = ang.theta, phi = ang.phi,  col="red4",
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2, front = "lines", back = "lines",
	alpha=0.4,
 	lit = FALSE, add = TRUE)

 

persp3d(coord.compo.grillax,coord.compo.grillay,mhat.lineal.rob.ilr,
	theta = ang.theta, phi = ang.phi,  
	ltheta = 120,   ticktype = "detailed", 
	xlab=" P ",ylab="L", zlab = "", 
	expand = 2,  
      col= "cyan", 
	  aspect=c(1,1), # y-size/x-size and z-size/x-size
      alpha=0.6,  
	 zlim=c(80,105), 
	 breaks = seq(80,105, length=11),
	#colkey = list(side = 2, length = 0.5,dist=-0.1,shift=0.17),
	facets=TRUE,   plot=T, add = TRUE) 

persp3d(coord.compo.grillax,coord.compo.grillay,mhat.lineal.rob.ilr,
	theta = ang.theta, phi = ang.phi,  col="blue4",
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2, front = "lines", back = "lines",
	alpha=0.4,
 	lit = FALSE, add = TRUE)


par3d( zoom=0.8 )

rgl.snapshot('AUC_cl_rob_Prot_LIP.png', fmt = 'png')


peli( spin3d(rpm=3), movie = "AUC_CL_ROB_Prot_LIP", duration=20, clean=FALSE , convert=TRUE)


################################################################
# CLASSICAL AND ROBUST ESTIMATORS   
# VERSUS PROTEINS AND CARBOHYDRATS
# IN THE PAPER WE CHOOSE THE FRAME 17 FOR BETTER VISUALIZATION
################################################################ 


par3d(windowRect = c(15, 15, 500, 500))

persp3d(coord.compo.grillax,coord.compo.grillaz,mhat.lineal.cl.ilr,
	theta = ang.theta, phi = ang.phi,  
	ltheta = 120,   ticktype = "detailed", 
	xlab=" P ",ylab="C", zlab = "", 
	expand = 2,  
      col= "red", 
	  aspect=c(1,1), # y-size/x-size and z-size/x-size
      alpha=0.6,  
	 zlim=c(80,105), 
	 breaks = seq(80,105, length=11),
	#colkey = list(side = 2, length = 0.5,dist=-0.1,shift=0.17),
	facets=TRUE,   plot=T) 

persp3d(coord.compo.grillax,coord.compo.grillaz,mhat.lineal.cl.ilr,
	theta = ang.theta, phi = ang.phi,  col="red4",
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2, front = "lines", back = "lines",
	alpha=0.4,
 	lit = FALSE, add = TRUE)


persp3d(coord.compo.grillax,coord.compo.grillaz,mhat.lineal.rob.ilr,
	theta = ang.theta, phi = ang.phi,  
	ltheta = 120,   ticktype = "detailed", 
	xlab=" P ",ylab="C", zlab = "", 
	expand = 2,  
      col= "cyan", 
	  aspect=c(1,1), # y-size/x-size and z-size/x-size
      alpha=0.6,  
	 zlim=c(80,105), 
	 breaks = seq(80,105, length=11),
	#colkey = list(side = 2, length = 0.5,dist=-0.1,shift=0.17),
	facets=TRUE,   plot=T,add=T  ) 


persp3d(coord.compo.grillax,coord.compo.grillaz,mhat.lineal.rob.ilr,
	theta = ang.theta, phi = ang.phi,  col="blue4",
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2, front = "lines", back = "lines",
	alpha=0.4,
 	lit = FALSE, add = TRUE)

par3d( zoom=0.8 )
 
rgl.snapshot('AUC_CL_ROB_Prot_CARB.png', fmt = 'png')


peli( spin3d(rpm=3), movie = "AUC_CL_ROB_Prot_CARB", duration=20, clean=FALSE , convert=TRUE)
 


################################################################
# CLASSICAL AND ROBUST ESTIMATORS   
# VERSUS LIPIDS AND CARBOHYDRATS
# IN THE PAPER WE CHOOSE THE FRAME 93 FOR BETTER VISUALIZATION
################################################################

par3d(windowRect = c(15, 15, 500, 500))

persp3d(coord.compo.grillay,coord.compo.grillaz,mhat.lineal.cl.ilr,
	theta = ang.theta, phi = ang.phi,  
	ltheta = 120,   ticktype = "detailed", 
	xlab=" L ",ylab="C", zlab = "", 
	expand = 2,  
      col= "red", 
	  aspect=c(1,1), # y-size/x-size and z-size/x-size
      alpha=0.6,  
	 zlim=c(80,105), 
	 breaks = seq(80,105, length=11),
	#colkey = list(side = 2, length = 0.5,dist=-0.1,shift=0.17),
	facets=TRUE,   plot=T) 

persp3d(coord.compo.grillay,coord.compo.grillaz,mhat.lineal.cl.ilr,
	theta = ang.theta, phi = ang.phi,  col="red4",
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2, front = "lines", back = "lines",
	alpha=0.4,
 	lit = FALSE, add = TRUE)

 

persp3d(coord.compo.grillay,coord.compo.grillaz,mhat.lineal.rob.ilr,
	theta = ang.theta, phi = ang.phi,  
	ltheta = 120,   ticktype = "detailed", 
	xlab=" L ",ylab="C", zlab = "", 
	expand = 2,  
      col= "cyan", 
	  aspect=c(1,1), # y-size/x-size and z-size/x-size
      alpha=0.6,  
	 zlim=c(80,105), 
	 breaks = seq(80,105, length=11),
	#colkey = list(side = 2, length = 0.5,dist=-0.1,shift=0.17),
	facets=TRUE,   plot=T, add = TRUE) 

persp3d(coord.compo.grillay,coord.compo.grillaz,mhat.lineal.rob.ilr,
	theta = ang.theta, phi = ang.phi,  col="blue4",
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2, front = "lines", back = "lines",
	alpha=0.4,
 	lit = FALSE, add = TRUE)


par3d( zoom=0.8 )

  rgl.snapshot('AUC_cl_rob_LIP_CARB.png', fmt = 'png')


peli( spin3d(rpm=3), movie = "AUC_CL_ROB_LIP_CARB", duration=20, clean=FALSE , convert=TRUE)


####################################
# ONLY ROBUST ESTIMATOR  
####################################
 
#################################################################
# ROBUST ESTIMATOR   VERSUS PROTEINS AND LIPIDS
# IN THE PAPER WE CHOOSE THE FRAME 82 FOR BETTER VISUALIZATION
#################################################################

par3d(windowRect = c(15, 15, 500, 500))

persp3d(coord.compo.grillax,coord.compo.grillay,mhat.lineal.rob.ilr,
	theta = ang.theta, phi = ang.phi,  
	ltheta = 120,   ticktype = "detailed", 
	xlab=" P ",ylab="L", zlab = "", 
	expand = 2,  
      col= "cyan", 
	  aspect=c(1,1), # y-size/x-size and z-size/x-size
      alpha=0.6,  
	 zlim=c(80,105), 
	 breaks = seq(80,105, length=11),
	#colkey = list(side = 2, length = 0.5,dist=-0.1,shift=0.17),
	facets=TRUE,   plot=T  ) 

persp3d(coord.compo.grillax,coord.compo.grillay,mhat.lineal.rob.ilr,
	theta = ang.theta, phi = ang.phi,  col="blue4",
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2, front = "lines", back = "lines",
	alpha=0.4,
 	lit = FALSE, add = TRUE)

  
par3d( zoom=0.8 )
 
rgl.snapshot('AUC_ROB_Prot_LIP.png', fmt = 'png')

 
peli( spin3d(rpm=3), movie = "AUC_ROB_Prot_LIP", duration=20, clean=FALSE , convert=TRUE)
 
###################################################################
# ROBUST ESTIMATOR  VERSUS PROTEINS AND CARBOHYDRATS
# IN THE PAPER WE CHOOSE THE FRAME 17 FOR BETTER VISUALIZATION
###################################################################


par3d(windowRect = c(15, 15, 500, 500))

persp3d(coord.compo.grillax,coord.compo.grillaz,mhat.lineal.rob.ilr,
	theta = ang.theta, phi = ang.phi,  
	ltheta = 120,   ticktype = "detailed", 
	xlab=" P ",ylab="C", zlab = "", 
	expand = 2,  
      col= "cyan", 
	  aspect=c(1,1), # y-size/x-size and z-size/x-size
      alpha=0.6,  
	 zlim=c(80,105), 
	 breaks = seq(80,105, length=11),
	#colkey = list(side = 2, length = 0.5,dist=-0.1,shift=0.17),
	facets=TRUE,   plot=T) 


persp3d(coord.compo.grillax,coord.compo.grillaz,mhat.lineal.rob.ilr,
	theta = ang.theta, phi = ang.phi,  col="blue4",
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2, front = "lines", back = "lines",
	alpha=0.4,
 	lit = FALSE, add = TRUE)


par3d( zoom=0.8 )
 
rgl.snapshot('AUC_ROB_Prot_CARB.png', fmt = 'png')


peli( spin3d(rpm=3), movie = "AUC_ROB_Prot_CARB", duration=20, clean=FALSE , convert=TRUE)
 

####################################################
# ROBUST ESTIMATOR VERSUS LIPIDS AND CARBOHYDRATS
# IN THE PAPER WE CHOOSE THE FRAME 93 FOR BETTER VISUALIZATION
#########################################################


par3d(windowRect = c(15, 15, 500, 500))

persp3d(coord.compo.grillay,coord.compo.grillaz,mhat.lineal.rob.ilr,
	theta = ang.theta, phi = ang.phi,  
	ltheta = 120,   ticktype = "detailed", 
	xlab=" L ",ylab="C", zlab = "", 
	expand = 2,  
      col= "cyan", 
	  aspect=c(1,1), # y-size/x-size and z-size/x-size
      alpha=0.6,  
	 zlim=c(80,105), 
	 breaks = seq(80,105, length=11),
	#colkey = list(side = 2, length = 0.5,dist=-0.1,shift=0.17),
	facets=TRUE,   plot=T  ) 


persp3d(coord.compo.grillay,coord.compo.grillaz,mhat.lineal.rob.ilr,
	theta = ang.theta, phi = ang.phi,  col="blue4",
	ltheta = 120,   ticktype = "detailed", 
	xlab=" ",ylab=" ", zlab = "", 
	expand = 2, front = "lines", back = "lines",
	alpha=0.4,
 	lit = FALSE, add = TRUE)


par3d( zoom=0.8 )
 
rgl.snapshot('AUC_ROB_LIP_CARB.png', fmt = 'png')


peli(spin3d(rpm=3), movie = "AUC_ROB_LIP_CARB", duration=20, clean=FALSE , convert=TRUE)
 

###################################################
# PLOTS OF THE ESTIMATES IN THE ILR SPACE
# FOR FIXED VALUES OF THE ILR-COORDINATES
###################################################

###################################
# x_2* =ejey[j] FIXED
##########################################
  


nombre= 'estimadores_distintos_x2star_rayas.pdf')
  
pdf(nombre,bg='transparent')

par(mfrow=c(3,3))

indices.ejey= seq(from=1, to=50, by=5)
for(k in indices.ejey[1:9]){
	par(mar=c(4, 3, 2, 2) )
	plot(ejex, mhat.lineal.cl.ilr[ ,k],type="l", col="red",lwd=2, 
		ylim=c(85,102),xlab=expression(x[1]^"*"), ylab="" )
	lines(ejey, mhat.lineal.rob.ilr[ ,k],type="l",col="blue",lwd=3)
	lines(ejey, mhat.lineal.cl.ilr.so[ ,k],type="l",col="gold",lwd=3, lty=2)
	valor=round(ejey[k],2)
	title(main=substitute(paste( x[2]^"*", " =  ", nn), list(nn=valor)),
	cex=0.6, col.main="green4")

}
 
dev.off()

###################################
# x_1* =ejex[j]  FIXED
##########################################

nombre= 'estimadores_distintos_x1star_rayas.pdf')

  
pdf(nombre,bg='transparent')

par(mfrow=c(3,3))

indices.ejex= seq(from=1, to=50, by=5)
for(j in indices.ejex[1:9]){
	par(mar=c(4, 3, 2, 2) )
	plot(ejey, mhat.lineal.cl.ilr[j,],type="l", col="red",lwd=2, 
		ylim=c(88,102),xlab=expression(x[2]^"*"), ylab="" )
	lines(ejey, mhat.lineal.rob.ilr[j,],type="l",col="blue",lwd=3)
	lines(ejey, mhat.lineal.cl.ilr.so[j,],type="l",col="gold",lty=2, lwd=3)
	valor=round(ejex[j],2)
	title(main=substitute(paste( x[1]^"*", " =  ", nn), list(nn=valor)),
	cex=0.6, col.main="green4")
}


dev.off()
