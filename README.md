The data set data_example.txt corresponds to the subset of A Estrada Glycation and Inflammation Study (AEGIS), trial NCT01796184, analysed in 

Ana M. Bianco, Graciela Boente, Wenceslao González--Manteiga, Francisco Gude Sampedro, Ana Pérez--González (2025). Robust Nonparametric Regression for Compositional Data: the Simplicial-Real case.  Available at  (https://arxiv.org/pdf/2405.12924)

The code to analyse the data corresponds to the file Analysis_AUC_KCV.R which calls the two source file containing the functions needed to compute the local linear estimators and to choose the bandwidth through K-fold cross-validation.

For the classical procedure least squares cross-validation is used, while for the robust one the robust $K-$fold cross-validation  defined in Section 2.4 is considered.
