# flm.fr

Code intended to implement a goodness-of-fit test for regression on Principal Components (PCs) for the Functional Linear Model (FLM) with Functional Response.

- flm.fourier.R: example of regression on Fourier basis expansions for the FLM with functional response.
- fregre.lr.pc.ex.R: example of regression on FPCs for the FLM with functional response.
- wb_fourier.R: implements the wild bootstrap and the simulation study on the projected residuals in a Fourier basis.
- wb_fpca.R: implements the wild bootstrap and the simulation study on the projected residuals in a FPCs basis.

The R directory contains exclusively .R functions, such as:

- flm_test.R: performs the goodness-of-fit test.
- fourier_expansion.R: computes the projection of a given functional variable onto a Fourier basis.
- fpc.R: Computes the FPCs for a given dataset.
- fregre.sr.pc.R: regression on FPCs for the FLM with scalar response.
- integrateSimp1D.R: implements the Simpson's rule in one dimension for equispaced data. Extracted from fda.usc
- integrateSimp2D.R: implements the Simpson's rule in two dimensions for equispaced data. Extracted from fda.usc
- linear_model.R: genertes a linear/non-linear model from provided X, the surface, the noise and the deviation are given.
- PCvM_statistic.R: implementation of the test statistic.
- pseudoinverse.R: computes a pseudo-inverse by means of a singular value decompsoition (SVD), needed for the estimatin of the model when N >> n.
- trap1D_unequal.R: implements the trapezoidal rule in one dimensions for non-equispaced data. Useful for extensions of the work
