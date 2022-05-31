# Ecological stability propagates across spatial scales and trophic levels in freshwater ecosystems

Code and data to reproduce the results in Siqueira et al. (submitted) published as a Preprint (https://ecoevorxiv.org/mpf5x)

# Analysis

The full set of results, including those made available as supplementary material, can be reproduced by running four scripts in the **R-codes** folder following this sequence:

- 02_Siqueira_etal_SEM_analyses.R
- 03_Siqueira_etal_stab_figs.R 
- 04_Siqueira_etal_stab_supp_m.R

and using the data available in the **Input_data** folder.

This is a collaborative effort and not all authors are allowed to share their raw data. For example, one data set (LEPAS) was not made available due to data sharing policies of The Ohio Division of Wildlife (ODOW). The orginal raw data include the abundance (individual counts, biomass, coverage area) of a given taxon, at a given site, in a given year. See details here https://ecoevorxiv.org/mpf5x

Thus, here, instead of starting the analyses with raw data, we start with data that has been generated with the code: 01_Siqueira_etal_dataprep_stability_metrics.R (available under request, together with the full set of data here https://zenodo.org/record/6591419)

These data include variability and synchrony components estimated using the methods described in Wang et al. (2019 Ecography; doi/10.1111/ecog.04290), diversity metrics (alpha and gamma diversity), and some variables describing the data. 
