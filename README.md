# DIANMF: Combined MS and MS/MS deconvolution of SWATH DIA data with the DIA-NMF workflow for comprehensive annotation in metabolomics


<!-- badges: start -->

<!-- [![Codecov test coverage](https://codecov.io/gh/odisce/mineMS2/graph/badge.svg)](https://app.codecov.io/gh/odisce/mineMS2) -->

<!-- badges: end -->

## Description

**DIANMF** is an open-source R package for deconvolving complex SWATH-DIA metabolomics data using sparse non-negative matrix factorization (NMF). It is the first method to jointly unmix MS1 and MS2 spectra in a fully untargeted manner, without relying predefined peak models or on spectral libraries. **DIA-NMF** enables precursor-level interpretation by recovering pure MS1 spectra and enriched, unmixed MS2 fragmentation patterns from all relavent isolation windows. This improves compound identification, especially for co-eluting and low-intensity metabolites.

The workflow detects MS1 peaks, extracts mixed MS1 and MS2 signals within minimally overlapping retention time windows, aligns them, and jointly unimixes these to recover pure precursor and fragment patterns (MS1 and MS2 spectra, respectively).
        
![DIA-NMF Workflow](vignettes/figures/dianmf_workflow.png)

        

## Installation

The package can be installed from GitHub with:

``` r
#install.packages("devtools")
devtools::install_github("odisce/DIANMF")
```

## Tutorials (vignettes)

One vignettes detail how to **deconvolute MS1 and SWATH-DIA data** at vignettes/Process_SWATH_DIA_Data.Rmd.


## Dataset

This package includes a sample dataset used in the accompanying examples and vignettes. It consists of selected regions from replicates 1 and 2 of a 10 ng/mL spiked human plasma sample, focusing on the 280–320 m/z range and 280–320 seconds retention time window (Barbier et al. 2020).

## Citation

Diana Karaki, Annelaure Damon, Antoine Souloumiac, François Fenaille, Etienne A. Thévenot and Sylvain Dechaumet. Combined MS and MS/MS deconvolution of SWATH DIA data with the DIA-NMF workflow for comprehensive annotation in metabolomics. DOI:soon. 

## Contacts

[diana.karaki\@cea.fr](mailto:diana.karaki@cea.fr), [sylvain.dechaumet\@cea.fr](mailto:sylvain.dechaumet@cea.fr), and [etienne.thevenot\@cea.fr](mailto:etienne.thevenot@cea.fr)

## Licence

[CeCILL licenses](http://www.cecill.info/index.en.html)

## References
  
Karaki Diana, Dechaumet Sylvain *et al.* (2024) Non-Negative Matrix Factorization of SWATH DIA Data Improves Global Metabolite Identification. *EUSIPCO*, **38**, 1967–1993. DOI:10.23919/EUSIPCO63174.2024.10715181.

Wang R. *et al.* (2019) Advancing untargeted metabolomics using data-independent acquisition mass spectrometry technology. *Analytical and bioanalytical chemistry*, 4349--4357 DOI:10.1007/s00216-019-01709-1. 

Barbier *et al.* (2020) Comparative Evaluation of Data Dependent and Data Independent Acquisition Workflows Implemented on an Orbitrap Fusion for Untargeted Metabolomics. Metabolites, 2218-1989. DOI:10.3390/metabo10040158.