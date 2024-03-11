This instruction aims to reproduce the results in the paper “*Stacking designs: designing multi-fidelity computer experiments with target predictive accuracy*”
by [Sung, Ji, Mak, Wang, and Tang
(2024)](https://epubs.siam.org/doi/full/10.1137/22M1532007), SIAM/ASA Journal on Uncertainty Quantification, 12(1), 157-181.

The code folder reproduces the result in Section 5 of the manuscript.

-   `currin_example.R` reproduces the results of Section 5.1.
-   `poisson_example.R` reproduces the results of Section 5.2.
-   `blade_example.R` reproduces the results of Section 5.3.

The following R packages are required to run the code:

-   `matlabr (>=1.5.2)`
-   `randtoolbox (>=2.0.2)`
-   `RColorBrewer (>=1.1-3)`
-   `ggplot2 (>=3.3.6)`
-   `gridExtra (>=2.3)`
-   `plgp (>=1.1-11)`

`MATLAB` needs to be installed in order to run the finite element
simulations. The MATLAB code can be accessed from the folder
`Rmatlab_files`.
