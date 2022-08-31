# Temporal variability declines with increasing trophic levels and spatial scales in freshwater ecosystems

Code and data to reproduce the results in Siqueira et al. (submitted) published as a Preprint (https://ecoevorxiv.org/mpf5x)

# Analysis

The full set of results, including those made available as supplementary material, can be reproduced by running five scripts in the **R_codes** folder following this sequence:

- 00_Renv_help_reproduc.R
- 01_Dataprep_stability_metrics.R 
- 02_SEM_analyses.R
- 03_Stab_figs.R
- 04_Stab_supp_m.R
- 05_Sensit_analysis.R

and using the data available in the **Input_data** folder.

This is a collaborative effort and not all authors are allowed to share their raw data. For example, one data set (LEPAS) was not made available due to data sharing policies of The Ohio Division of Wildlife (ODOW). The orginal raw data include the abundance (individual counts, biomass, coverage area) of a given taxon, at a given site, in a given year. See details here https://ecoevorxiv.org/mpf5x

Thus, here, instead of starting the analyses with raw data, we start with data that has been generated with the code: 01_Dataprep_stability_metrics.R (available under request, together with the full set of data here https://zenodo.org/record/6591419)

These data include variability and synchrony components estimated using the methods described in Wang et al. (2019 Ecography; doi/10.1111/ecog.04290), diversity metrics (alpha and gamma diversity), and some variables describing the data. 

# Reproducibility

To improve reproducibility, I’m using the R package renv to register the versions of the R packages I use and to manage a local library that doesn’t affect the rest of my system. 

renv::init() is used to prepare a lockfile that records the exact versions of R packages I used in this project.

Anyone who wants to reproduce the results described in the preprint can just download the whole R project (that incluses code and data) and run codes from 00 to 05. 

The renv folder is compressed here. One needs extract the renv folder within the same folder as the R project.

Alternatevely, I am making the whole R project folder (with everything needed to reproduce the results) available as a compressed file:
Stab_R_project_full.zip
