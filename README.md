This folder consists of the R code for the manuscript “*Stacking designs: designing multi-fidelity computer experiments with confidence*” by Sung, Ji, Tang, and Mak. 

The code folder reproduces the result in Section 5 of the manuscript. 

- `currin_example.R` reproduces the results of Section 5.1. 
- `poisson_example.R` reproduces the results of Section 5.2. 
- `blade_example.R` reproduces the results of Section 5.3. 
	- Note: this code could take a long time because it runs very high fidelity simulations.
	

The following R packages are required to run the code:
  - `matlabr (>=1.5.2)`
  - `randtoolbox (>=2.0.2)`
  - `RColorBrewer (>=1.1-3)`
  - `ggplot2 (>=3.3.6)`
  - `gridExtra (>=2.3)`
  - `plgp (>=1.1-11)`

The reproducibility is demonstrated as follows.
``` r
source("currin_example.R")
source("poisson_example.R")
source("blade_example.R")
```
