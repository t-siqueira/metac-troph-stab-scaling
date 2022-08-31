## Code for reproducing results in Siqueira et al. (2022):
# "Ecological stability propagates across spatial scales and 
# trophic levels in freshwater ecosystems"
# Tadeu Siqueira, Rio Claro, SP, Brazil, 2022-08-23

# To improve reproducibility, I’m using renv to register the versions of the 
# R packages I use and to manage a local library that doesn’t affect the rest 
# of my system. 

# renv::init() is used to prepare a lockfile that records the exact versions 
# of R packages I used in this project

# Anyone who wants to reproduce the results described in the preprint 
# could download the whole R project (that incluses code and data), run: 

renv::restore() 

# and have an R environment very similar to the one I use.

# For more information on the renv package see, for example:
# https://mran.microsoft.com/snapshot/2019-11-29/web/packages/renv/vignettes/renv.html
# https://www.r-bloggers.com/2022/01/how-renv-restores-packages-from-r-universe-for-reproducibility-or-production/
#====



