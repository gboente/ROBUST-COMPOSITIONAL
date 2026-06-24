##############################################
# FUNCTION TO RUN THE SIMULATION
# UNDER C_{1,0.10, 10}
# 10% CONTAMINATION WITH ERRORS WITH MEAN 10
###############################################

rm(list=ls())
 
source('funciones.R')
source('genero.R')
source('simula.R')

d=3
a=c(5, 7, 1)
sigma=1
H=2*diag(2)
b = c( 0.05920067, 0.7193872, 0.2214121)
 
tipo.regre='bis'

nrep=1000 

n=100
ntest=100


ee_mu = 10
ee_sd  = 0.1
conta=  list(type = 'C1',ratio=0.1,  ee_mu =ee_mu, ee_sd =ee_sd )

simula.homo(nrep=nrep, n=n,ntest=ntest, H,  d=d, a=a, sigma=sigma, b=b,conta=conta, tipo.regre=tipo.regre)


