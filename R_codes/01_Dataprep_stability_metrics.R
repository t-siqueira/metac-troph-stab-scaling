## Code for reproducing results in Siqueira et al. (2022):
# "Ecological stability propagates across spatial scales and 
# trophic levels in freshwater ecosystems"
# Tadeu Siqueira, Rio Claro, SP, Brazil, 2022-05-01

# Here we do some data preparation, estimate stability and diversity metrics,
# and save these data for subsequent analyses with other code:
# "02_Siqueira_etal_SEM_analyses.R"

# All data sets, except one, used in the paper are entirely available.
# The one data set (LEPAS) that was not made available - due to data sharing 
# policies of The Ohio Division of Wildlife (ODOW) - can still be used in this
# code. We made them available not in the format of raw numbers, but as 
# pre-processed variables that we used in our models (e.g., stability and 
# diversity metrics, trophic groups, etc.). 

## Packages
library(tidyverse)
library(vegan)

## Code to partition variability and synchrony (Wang et al. 2019; Ecography)
source ("R_codes/var_part_wang.R") 

#===============================
## Data preparation 

## Import all selected data
# Each row here is the abundance of a given species in a given site and time
# Various other columns describe these row (lat, long, sample id, etc.) 

# This data table has been pre-processed to replace NAs in species abundances
# by the median of the abundance of that species in other surveys.

full.bio.data <- as_tibble(
  read.csv("Input_data/data_stab_analysis_except_lepas.csv", T))
# Not sure why read_csv read the data weirdly

# Load the pre-processed LEPAS data sets that will be merged with the rest of
# the data sets later. The raw version of these data is not available.
# See details in supplementary information (Siqueira et al. 2022)

lepas.site.var <- read_csv("Input_data/lepas_site_stab_metrics.csv")
lepas.meta.var <- read_csv("Input_data/lepas_meta_stab_metrics.csv")  

## Prepare data according to Wang et al. (2019) 
# Wang's code requires data to be organized as arrays

# "The input array "metacomm_tsdata" is an N*T*M array. 
# The first dimension represents N species, the second represents 
# time-series observations of length T, and the third represents 
# M local communities".

# Data should be arranged by site/time/species

full.bio.data %>% 
  filter(is.na(New_tr_g)) %>% 
  select(Metacom) %>% 
  table()

# Remove data with no trophic info
full.data1 <- full.bio.data %>%
  mutate(across(where(is.factor), as.character)) %>%
  filter (!is.na(New_tr_g))

full.data1 %>% 
  filter(is.na(New_tr_g)) %>% 
  select(Metacom) %>% 
  table()

# Each of these below is a independent metacommunity - a set of sites sampled
# by data providers that are thought to be linked by dispersal

with(full.data1, unique(Metacom))

# Remove pond and isolate lakes data
# These were included in the first version of the preprint, but we decided to
# remove them after a reviewer suggested we kept only metacommunities
# that were made of communities potentially phisically connected.

full.data.no.ponds <- full.data1 %>% 
  filter (Metacom != "ELA_1", Metacom != "ELA_2", Metacom != "ELA_3", 
          Metacom != "NTL", Metacom != "M_Hill")

with(full.data.no.ponds, unique(Metacom))


## Create data tables that represent trophic levels
# Because I have more than one trophic group per metacommunity I need a 
# new variable.
# As I will likely need a site_troph variable, I will create this too.

full.data2 <- full.data.no.ponds %>%
  mutate (Meta_troph = paste0(Metacom, "_", New_tr_g),
          Site_troph = paste0(Site, "_", New_tr_g))

# Re-organize
full.data2 <- full.data2 %>%
  select(-Data_set_ID) %>% 
  relocate(Metacom,Meta_troph,Site,Site_troph, Sample_id,Year,Frequency,
           Time_step,Lat,Long,Datum,Coord_system,Bio_gr,Ecosys,Trophic_gr,
           New_tr_g, Species, Abundance)

## Export this as it will be used in a sensibility analysis described in
# code "05_Siqueira_etal_sensit_analysis.R"

write_csv(full.data2, "Input_data/Df_sens_analysis.csv")

# I will use a list from now on as I think it is easier to use lapply than 
# various for loops, and I am still incompetent to purr functions.

long.median2 <- split(full.data2, f = full.data2$Meta_troph, drop = T)

length(long.median2)
names(long.median2) # here, each element is a unique combination of 
                    # a metacommunity and trophic group

## Arrange it the way Wang's function requires

long.median2 <- lapply(long.median2, function (x){
  arrange(x, Site_troph, Time_step, Species)})

## Create the arrays

all.array <- lapply(long.median2, function (x){ 
  array(data = x$Abundance,
        dim = c(length(unique(x$Species)), 
                length(unique(x$Time_step)), 
                length(unique(x$Site_troph))),
        dimnames = list(unique(x$Species),
                        unique(x$Time_step),
                        unique(x$Site_troph)))
})

## Check if data is well organized as array
all.array[[1]][,,1] # this is site 1 along time

# Test
wide.test = as_tibble(
  reshape2::dcast(select(long.median2[[1]], -Trophic_gr, -New_tr_g), 
                    formula = ... ~ Species, 
                    value.var = "Abundance", fun = sum))

matrix.real =
  wide.test %>%
  arrange (Site_troph, Time_step) %>%
  filter (Site_troph == "1701_015-06-IS_Primary") %>%
  select (-c(1:14)) %>% #pay attention to the number of columns
  t()

# If things went well the line above should result in a matrix of TRUEs
matrix.real == all.array[[1]][,,1]

#===============================
## Run partitioning analysis following Wang et al. 2019  

part.W <- lapply(all.array, function (x) var.partition(x))
part.W

# The name of each metric follows Wang et al. (2019)

## Organize results
resul.all <- do.call(rbind, part.W) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Meta_troph") %>% 
  as_tibble() 
resul.all

# Organize per metacom
df.info.all <- lapply(long.median2, function (x) {
  data.frame(Metacom = unlist(unlist(as.character(unique(na.omit(x$Metacom))))),
             Meta_troph = unlist(unlist(as.character(unique(na.omit(x$Meta_troph))))),
             Freq = unlist(unlist(as.character(unique(na.omit(x$Frequency))))),
             Time_step = unlist(unlist(as.numeric(unique(max(na.omit(x$Time_step)))))),
             Bio_gr = unlist(unlist(as.character(unique(na.omit(x$Bio_gr))))),
             Trophic_gr = unlist(unlist(as.character(unique(na.omit(x$Trophic_gr))))),
             New_tr_g = unlist(unlist(as.character(unique(na.omit(x$New_tr_g))))),
             Ecosys = unlist(unlist(as.character(unique(na.omit(x$Ecosys))))),
             nSites = unlist(unlist(as.numeric(length(unique(x$Site))))))})

df.info.all

resul.all2.metac <- Reduce(rbind, df.info.all) %>% 
  left_join(resul.all, by = "Meta_troph") %>% 
  as_tibble()

resul.all2.metac

#================================
## Extract site level stability metrics

# This function was provided by Lise Comte
site_level_vars <- function (metacomm_tsdata) {
  
  ###calculate community synchrony (species synchrony in Wang & Loreau 2019)
  #
  ts_patch <- apply(metacomm_tsdata,c(2,3),sum) #total biomass of local comms per time
  sd_patch_k <- apply(ts_patch,2,sd) #temporal var. in comm biomass 
  sd_species_patch_ik <- apply(metacomm_tsdata,c(1,3),sd) #temporal var. in spp. biomass per comm.
  summed_sd_sp <- apply(sd_species_patch_ik,2,sum)  #summed SD across species
  synchrony_comm_site <- sd_patch_k/summed_sd_sp
  
  ###calculate mean species variability
  mean_species_patch_ik <-apply(metacomm_tsdata,c(1,3),mean)
  cv_species_site <- sd_species_patch_ik/mean_species_patch_ik
  mean_cv_species_site <- apply(cv_species_site, 2, mean, na.rm=T)
  
  ###calculate community variability
  #
  mean_patch_k <- apply(ts_patch,2,mean)
  cv_comm_site <- sd_patch_k/mean_patch_k
  
  return (tibble(Site_troph = names(cv_comm_site), synchrony_comm_site, 
                 mean_cv_species_site, cv_comm_site))
}

# Run the function
resul.site.metrics <- Reduce(rbind,
  lapply(all.array, site_level_vars))

# Organize results
df.site.info <- Reduce(rbind,
                      lapply(long.median2, function (x) 
                        na.omit(distinct(x, Metacom, Meta_troph, Site, Site_troph,
                                         Bio_gr, Ecosys, New_tr_g))))

resul.site.metrics2 <- left_join(resul.site.metrics, df.site.info, 
                                 by = "Site_troph")

# This data table should have the same number of rows as the number of sites in
# your study

resul.site.metrics2

#================================
## Extract biodiversity predictors

# I will use Simpson diversity index and observed species richness
# Create a list of lists where each element is a list of sites 

list.of.lists1 <- lapply(long.median2, function (x) split(x, f = x$Site_troph))

# e.g.,
with(list.of.lists1$`1701_Primary`$`1701_015-06-IS_Primary`, 
  table(Site, Time_step)) 
# 28 is the number of species names, including those with zero abundance
# and that's why all time_steps have same number

# Function to estimate Simpson diversity
# This will result in a measure of div per site per time step
# I am using the median of values across time steps as alpha diversity
# This is not computationally efficient, but I want to keep the same data 
# structure I used to estimate variability and synchrony

func_simpson <- function (x) {
  unlist(lapply(
    lapply(x, function (xx) {
      xx %>%
        pivot_wider(names_from = Species, values_from = Abundance) %>%
        select(-Metacom:-New_tr_g) %>%
        vegan::diversity(index = "simpson")}), 
    median))}

# Run it
alpha.resul1 <- lapply(list.of.lists1, func_simpson)

# Test if the function worked was expected
alpha.resul1$`1701_Primary` 

median(apply(X = list.of.lists1$`1701_Primary`$`1701_015-06-IS_Primary` %>% 
        pivot_wider(names_from = Species, values_from = Abundance) %>% 
        select(-Metacom:-New_tr_g),
      MARGIN = 1,
      FUN = function (x) vegan::diversity(x, index = "simpson")))
#

## Now a function to calculate observed raw species richness
func_s <- function (x) {
  unlist(lapply(
    lapply(x, function (xx) {
      xx %>%
        pivot_wider(names_from = Species, values_from = Abundance) %>%
        select(-Metacom:-New_tr_g) %>%
        vegan::specnumber()}), 
    median))}

alpha.resul1.2 <- lapply(list.of.lists1, func_s)

# Test if the function worked was expected
alpha.resul1.2$`1701_Primary` 

median(apply(X = list.of.lists1$`1701_Primary`$`1701_015-06-IS_Primary` %>% 
               pivot_wider(names_from = Species, values_from = Abundance) %>% 
               select(-Metacom:-New_tr_g),
             MARGIN = 1,
             FUN = function (x) vegan::specnumber(x)))
#

# Organize this so that it matches the format of the stability metric tibble
alpha.resul2 <- lapply(alpha.resul1, function (x) {
  tibble (Site_troph = names(x), Simp = x)
})

alpha.resul2.1 <- lapply(alpha.resul1.2, function (x) {
  tibble (Site_troph = names(x), S = x)
})

# Join this with site metrics data

resul.final.site <- do.call(rbind, alpha.resul2) %>%
  right_join(resul.site.metrics2, by = "Site_troph") %>% 
  right_join(do.call(rbind, alpha.resul2.1), by = "Site_troph") 
 
## Gamma diversity
# Again I will use Simpson index and species richness

# Prepare data
list1.gamma <- lapply(lapply(long.median2, function (x) { 
  x %>%
    group_by(Time_step, Species) %>%
    summarise(Abundance_meta = sum(Abundance)) %>%
    pivot_wider(id_cols = Species, names_from = Time_step, 
                values_from = Abundance_meta) %>%
    select(-Species) %>%
    vegan::decostand(method = "pa")}), function (z) z[rowSums(z)>0,])

# Estimates
resul.gamma <- as_tibble(do.call(rbind, lapply(list1.gamma, function(x) {
  median(vegan::diversity(t(x), index = "simp"))}))) %>%
  mutate (Meta_troph = resul.all2.metac$Meta_troph, Simp_gamma = V1) %>%
  select(-V1)

resul.gamma2 <- as_tibble(do.call(rbind, lapply(list1.gamma, function(x) {
  median(vegan::specnumber(t(x)))}))) %>%
  mutate (Meta_troph = resul.all2.metac$Meta_troph, S_gamma = V1) %>%
  select(-V1)

# Create a Metacom data for this
resul.final.metac <- resul.final.site %>%
  select(Simp, S, Meta_troph) %>%
  group_by(Meta_troph) %>%
  summarise(across(.fns = mean)) %>%
  ungroup() %>%
  right_join(resul.all2.metac, by = "Meta_troph") %>%
  right_join(resul.gamma, by = "Meta_troph") %>% 
  right_join(resul.gamma2, by = "Meta_troph")

#=============================================================================
## Final organization and file export
# Bind the rows of the Lepas data sets

resul.final.site1 <- full.data2 %>% 
  select(Site_troph, Lat, Long) %>% 
  distinct() %>% 
  left_join(resul.final.site, by = "Site_troph") %>% 
  bind_rows(lepas.site.var)

write_csv(resul.final.site1, "Input_data/site_stab_metrics.csv")

resul.final.metac %>% 
  bind_rows(lepas.meta.var) %>% 
  write_csv("Input_data/meta_stab_metrics.csv")

## End
#=============================================================================








