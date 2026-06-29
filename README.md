This repository contains the code that implements the robust estimators for nonparametric regression introduced in

Ana M. Bianco, Graciela Boente, Wenceslao González--Manteiga, Francisco Gude Sampedro, Ana Pérez--González (2025). Robust Nonparametric Regression for Compositional Data: the Simplicial-Real case. Available at (https://arxiv.org/pdf/2405.12924)

The data set AEGIS used in the paper is given in the file data_example.txt.

The file Analysis_AUC_KCV.R, used to analize the data, calls the two source files, funciones.R and funciones_KCV.R, containing the functions needed to compute the local linear estimators and to choose the bandwidth through K-fold cross-validation.

For the classical procedure least squares cross-validation is used, while for the robust one the robust $K-$fold cross-validation is considered, as defined in Section 2.4 of the above mentioned paper.
