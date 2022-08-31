## Code for reproducing results in Siqueira et al. (2022):
# "Ecological stability propagates across spatial scales and trophic 
# levels in freshwater ecosystems"
# Tadeu Siqueira, Rio Claro, 2022-05-01

# Here we run structural equation models (SEM) using data prepared with code
# in Siqueira_etal_dataprep_stability_metrics.R and data on environmental 
# predictors. 

# Packages
library(tidyverse)
library(piecewiseSEM)
library(lme4)
library(nlme)
library(igraph)

#===============================
## Data preparation

# Import regional and local stability metrics
# These were created with code: 
# "01_Siqueira_etal_dataprep_stability_metrics.R"

meta.metr <- read_csv("Input_data/meta_stab_metrics.csv", T)
site.metr <- read_csv("Input_data/site_stab_metrics.csv", T)

# Import env predictors
# See methods and supplementary information in Siqueira et al. (2022) for
# details about how environmental and spatial predictors were obtained.

local.preds <- read_csv("Input_data/site_env_preds.csv", T) 
meta.preds <- read_csv("Input_data/meta_env_spa_preds.csv", T) 

# Merge these tibbles and reorder trophic gr levels

local.metr.full <- site.metr %>% 
  right_join(local.preds, by = "Site_troph") %>% 
  mutate(New_tr_g = factor(New_tr_g, levels = c("Producers","Primary",
                                                "Secondary", "Tertiary"))) 

meta.metr.full <- meta.metr %>% 
  right_join(meta.preds, by = "Meta_troph") %>% 
  mutate(New_tr_g = factor(New_tr_g, levels = c("Producers","Primary",
                                                "Secondary", "Tertiary")))

# Remove pond and isolate lakes data
# This was a request made by one of the reviewers - and it makes sense

local.metr.full1 <- local.metr.full %>% 
  filter (Metacom != "ELA_1", Metacom != "ELA_2", Metacom != "ELA_3", 
          Metacom != "NTL", Metacom != "M_Hill")

meta.metr.full1 <- meta.metr.full %>% 
  filter (Metacom != "ELA_1", Metacom != "ELA_2", Metacom != "ELA_3", 
          Metacom != "NTL", Metacom != "M_Hill")

# Get some descriptive stats of the data
with(local.metr.full1, unique(Metacom))
with(local.metr.full1, unique(Site))
with(meta.metr.full1, unique(Metacom))

local.metr.full1 %>% 
  select (Ecosys, Metacom) %>% 
  distinct () %>% 
  table () %>% 
  rowSums ()

meta.metr.full1 %>% 
  select (nSites, Metacom) %>%
  distinct () %>% 
  summary (nSites)

meta.metr.full1 %>% 
  select (Time_step, Metacom) %>%
  distinct () %>% 
  summary (Time_step)

meta.metr.full1 %>%
  group_by (New_tr_g) %>% 
  summarise (meta_tr = length(unique(Meta_troph)))

meta.metr.full1 %>%
  select(Meta_troph) %>% 
  distinct ()

#===============================
# For the local scale SEM, I must check with there are NAs as site variability and
# synchrony metrics might not have been calculated for sites that had one 
# species in a given time step 

which(is.na(local.metr.full1), arr.ind=TRUE)
names(local.metr.full1)

# I will also remove sites with less than 2 (median) species
# and make sure that this does not lead to metas with less than 4 sites in 
# the data set

local.metr.full2 <- local.metr.full1 %>% 
  left_join(distinct(select(meta.metr.full1, Metacom, Time_step)), 
            by = "Metacom") %>% 
  filter(!is.na(cv_comm_site)) %>% 
  filter(S >= 2) 

new.nSites <- local.metr.full2 %>% 
  group_by(Meta_troph) %>% 
  summarise(nSites = length(Site_troph)) %>% 
  ungroup()

# Merge and exclude metas with less than 4 sites

local.metr.full3 <- local.metr.full2 %>% 
  left_join(new.nSites, by = "Meta_troph") %>% 
  filter(nSites >= 4)

meta.metr.full2 <- new.nSites %>% 
  right_join(select(meta.metr.full1, -nSites), by = "Meta_troph") %>% 
  ungroup() %>% 
  filter(nSites >= 4)
#====

## Prepare data for the regional scale SEM
# including calculating spatial predictors

which(is.na(meta.metr.full2), arr.ind=TRUE)
# NAs might appear as we changed the sites we kept within each metacom

meta.metr.full3 <- meta.metr.full2 %>% 
  filter(S_gamma >= 2) 

# Spatial predictors from igraph

list.1 = split(local.metr.full3, f = local.metr.full3$Meta_troph)

resul.igraph <-
  lapply (list.1, function (x) {
    g = make_full_graph(nrow(x))
    E(g)$weight = as.vector(dist(select(x, Lat, Long)))
    resul = tibble(Site_troph = x$Site_troph, Clos = closeness(g, normalized = T))
    return(resul)
  })

meta.metr.full4 <- right_join(local.metr.full3, do.call(rbind, resul.igraph), 
                        by = "Site_troph") %>% 
  group_by(Meta_troph) %>% 
  summarise(Meta_clos = mean(Clos)) %>% 
  right_join(meta.metr.full2, by = "Meta_troph")

# Get more descriptive stats of the data
with(local.metr.full3, unique(Site_troph))
with(local.metr.full3, unique(Meta_troph))

# These will be exported to be used in figure prep.

write_csv(meta.metr.full4, "Input_data/meta_metr_figs.csv")
write_csv(local.metr.full3, "Input_data/site_metr_figs.csv")

### Local scale SEM
## Run SEM and model selection
# log transform the before SEM

local.metr.full4 <- local.metr.full3 %>%
  mutate(cv_comm_site = log(cv_comm_site), 
         synchrony_comm_site = log(synchrony_comm_site),
         mean_cv_species_site = log(mean_cv_species_site), 
         Simp = log(Simp), S = log(S)) 

# Basic model with env predictors + spatial predictors
# Other models will be variations of this one
# Although previous tests indicated that a glmmPQL (family = Gamma(link=log))
# would result in better looking residuals, psem does not produces standardized
# coefficients for this family. As we are interested in comparing coeffs. among
# trophic groups, I decided to use lme even if this resulted in some
# structure in residuals.

# Here meta_troph are within meta, so we need to specify this in the random
# model structure

names(local.metr.full4)

# The first model is a representation of our concept model
# This might need to change if d-separation tests indicate a path should 
# be included to improve fit

site.cv0 <- psem(
  lme(cv_comm_site ~ synchrony_comm_site + mean_cv_species_site + bio15 + Time_step + S, 
      random = ~ 1|Metacom/Meta_troph, data = local.metr.full4),
  lme(synchrony_comm_site ~ S, 
      random = ~ 1|Metacom/Meta_troph, data = local.metr.full4),
  lme(mean_cv_species_site ~ S + bio15 + Time_step, 
      random = ~ 1|Metacom/Meta_troph, data = local.metr.full4),
  data = as.data.frame(local.metr.full4))

summary(site.cv0) # good fit, but only after I include bio15 + S as predictors
# of cv_comm_site, as suggested by the d-separation test

(mt.cv0 <- multigroup(site.cv0, group = "New_tr_g", standardize = "range"))

# Alternative model with Simpson diversity
site.cv1 <- psem(
  lme(cv_comm_site ~ synchrony_comm_site + mean_cv_species_site + bio15 + Simp + Time_step, 
      random = ~ 1|Metacom/Meta_troph, data = local.metr.full4),
  lme(synchrony_comm_site ~ Simp, 
      random = ~ 1|Metacom/Meta_troph, data = local.metr.full4),
  lme(mean_cv_species_site ~ Simp + bio15 + Time_step, 
      random = ~ 1|Metacom/Meta_troph, data = local.metr.full4),
  data = as.data.frame(local.metr.full4))

summary(site.cv1) # not a good fit

# I wont include Simpson diversity in any of the subsequent alternative models

# Replace bio15 by bio4
site.cv2 <- psem(
  lme(cv_comm_site ~ synchrony_comm_site + mean_cv_species_site + bio4 + S + Time_step, 
      random = ~ 1|Metacom/Meta_troph, data = local.metr.full4),
  lme(synchrony_comm_site ~ S, 
      random = ~ 1|Metacom/Meta_troph, data = local.metr.full4),
  lme(mean_cv_species_site ~ S + bio4 + Time_step, 
      random = ~ 1|Metacom/Meta_troph, data = local.metr.full4),
  data = as.data.frame(local.metr.full4))

summary(site.cv2) # Not a good fit

# Use both

site.cv3 <- psem(
  lme(cv_comm_site ~ synchrony_comm_site + mean_cv_species_site + S + Time_step + bio15 + bio4, 
      random = ~ 1|Metacom/Meta_troph, data = local.metr.full4),
  lme(synchrony_comm_site ~ S, 
      random = ~ 1|Metacom/Meta_troph, data = local.metr.full4),
  lme(mean_cv_species_site ~ S + bio4 + bio15 + Time_step, 
      random = ~ 1|Metacom/Meta_troph, data = local.metr.full4),
  bio4 %~~% bio15,
  data = as.data.frame(local.metr.full4))

# Check fit
summary(site.cv3) # marginal good fit
(mt.cv3 <- multigroup(site.cv3, group = "New_tr_g", standardize = "range")) 
# not a good fit

# Check residuals
hist(residuals(site.cv0)[,1])
hist(residuals(site.cv0)[,2])
hist(residuals(site.cv0)[,3]) # not horrible

# Multigroup analyses
mt.cv0$Cstat
mt.cv0$group.coefs # I need to use the coefs of this model where 
# Model-wide Interactions are associated with p < 0.05
summary(site.cv0)$coeff # and of this model when they are constrained by 
# the global model

#====================
## Regional scale SEM
# I will use the community route (sensu Wang et al. 2019); it is also possible
# to use the metapopulation route.

# Log transform response variables
meta.metr.full4 <- meta.metr.full4 %>%
  mutate(CV_C_R = log(CV_C_R), phi_S2C_R = log(phi_S2C_R), S_gamma = log(S_gamma),
         Simp_gamma = log(Simp_gamma))  %>% 
  select(-S, -Simp)

# Full model
# Because I have already modeled CV_C_L above, I will only include it as a 
# predictor here, not as a response. These are different models though. 
# Here I only need to inform Metas as random 
# I will compare models with different predictor variables representing
# env synchrony and spatial isolation and diversity metrics 

# Start with only one env sync predictor

names(meta.metr.full4)

metac.sem0 <- psem(
  lme(CV_C_R ~ phi_C_L2R + CV_C_L, random = ~ 1|Metacom,
     data = meta.metr.full4),
  lme(phi_C_L2R ~ Meta_clos + Tmin_sync, random = ~ 1|Metacom,
     data = meta.metr.full4),
  data = as.data.frame(meta.metr.full4))

summary(metac.sem0) # good fit even if d-separation tests indicate a path that
# was not part of our conceptual model

(mt.cv.r0 <- multigroup(metac.sem0, group = "New_tr_g", standardize = "range"))
# good fit

# Alternative models

metac.sem2 <- psem(
  lme(CV_C_R ~ phi_C_L2R + CV_C_L, random = ~ 1|Metacom,
      data = meta.metr.full4),
  lme(phi_C_L2R ~ Meta_clos + Tmax_sync, random = ~ 1|Metacom,
      data = meta.metr.full4),
  data = as.data.frame(meta.metr.full4))

summary(metac.sem2) # not a good fit
(mt.cv.r2 <- multigroup(metac.sem2, group = "New_tr_g", standardize = "range"))

#
metac.sem4 <- psem(
  lme(CV_C_R ~ phi_C_L2R + CV_C_L, random = ~ 1|Metacom,
      data = meta.metr.full4),
  lme(phi_C_L2R ~ Meta_clos + Precip_sync, random = ~ 1|Metacom,
      data = meta.metr.full4),
  data = as.data.frame(meta.metr.full4))

summary(metac.sem4) # good fit even if d-separation tests indicate a path that
# was not part of our conceptual model

(mt.cv.r4 <- multigroup(metac.sem4, group = "New_tr_g", standardize = "range"))
# moderate good fit

# Model selection including models with good fit only (here, all of them)

all.metac.sems <- unlist(c(summary(metac.sem0)$IC[2], summary(metac.sem2)$IC[2],
                          summary(metac.sem4)$IC[2]))

names(all.metac.sems) = c("mod0", "mod2", "mod4")

wiqid::AICtable(all.metac.sems)

# All are equally plausible.
# Let's check R2 values
summary(metac.sem0)$R2
summary(metac.sem4)$R2 # this is a bit better
summary(metac.sem2)$R2 


# Check residuals
hist(residuals(metac.sem0)[,1]) # not very good
hist(residuals(metac.sem0)[,2]) # not horrible
hist(residuals(metac.sem4)[,1])# not very good
hist(residuals(metac.sem4)[,2])
hist(residuals(metac.sem2)[,1])# not very good
hist(residuals(metac.sem2)[,2])

# I will interpret the one with the highest R2
mt.cv.r4$Cstat
mt.cv.r4$group.coefs # I need to use the coefs from this model where 
# Model-wide Interactions are associated with p < 0.05
summary(metac.sem4)$coeff # and of this model when they are constrained by 
# the global model

# End
# ==============================================================================

