## Code for reproducing results in Siqueira et al. (2022):
# "Ecological stability propagates across spatial scales and trophic 
# levels in freshwater ecosystems"
# Tadeu Siqueira, Rio Claro, 2022-08-02

# Here we prepare sensitivity analysis we show in supp. information. 

## Packages
library(tidyverse)
library(patchwork)
source("R_codes/Split_viol_plot.R")

## Code to partition variability and synchrony (Wang et al. 2019; Ecography)
source ("R_codes/var_part_wang.R")


## Input data
# Original full data
full.data <- read.csv("Input_data/Df_sens_analysis.csv", T)
with(full.data, unique(Metacom))

# Data to be compared with rarefied data
meta.metr <- read_csv("Input_data/meta_metr_figs.csv", T)
with(meta.metr, unique(Metacom)) # here we have Lepas also
# but this is not a problem as it is not included in the sensitivity analysis

# Let's find out the number of sites to be used in the sensitivity analysis

full.data %>% 
  group_by(New_tr_g, Metacom) %>% 
  summarise(N_sites = length(unique(Site))) %>% 
  group_by(New_tr_g) %>% 
  summarise(Max_sites = max(N_sites))

full.data %>% 
  group_by(New_tr_g, Metacom) %>% 
  summarise(N_sites = length(unique(Site))) %>% 
  group_by(New_tr_g) %>% 
  summarise(Median_sites = median(N_sites), 
            Mean_sites = mean(N_sites))

# For secondary and tertiary consumers, I will get all data sets with more than
# 7 sites and rarefy to 8 sites 1000 times. At each run I will estimate
# variability and synchrony metrics. Finally I will average them and compare
# with the full data.

# Sample a reduced number of sites
# Get only secondary and tertiary consumers with more than 13 sites

sec.ter <- full.data %>% 
  filter (New_tr_g == "Secondary" | New_tr_g == "Tertiary") %>% 
  group_by (New_tr_g, Metacom) %>% 
  summarise (N_sites = length(unique(Site))) %>% 
  filter (N_sites > 7) %>%
  ungroup () %>% 
  select (-N_sites) %>% 
  inner_join (full.data, by = c("Metacom", "New_tr_g"))

unique(sec.ter$New_tr_g)


# I will use a list from now on as I think it is easier to use lapply than 
# various for loops, and I am still incompetent to purr functions.
# But we need to create a variable that combines site with troph group first

sec.ter.list <- split(sec.ter, f = sec.ter$Meta_troph, drop = T)

names(sec.ter.list) # here, each element is a unique combination of 
# a metacommunity and trophic group

## Arrange it the way Wang's function requires

sec.ter.list <- lapply(sec.ter.list, function (x){
  arrange(x, Site_troph, Time_step, Species)})

# Rarefy the n of sites here

func_raref <- function (a.list, raref_n, runs) {
  
  resul.list <- list()
  
  for (x in 1:runs) {
    resul.list[[x]] <- a.list[a.list$Site %in% 
                                sample(unique(a.list$Site), raref_n),]
  }
  
  list.array <- lapply(resul.list, function (x){ 
    array(data = x$Abundance,
          dim = c(length(unique(x$Species)), 
                  length(unique(x$Time_step)), 
                  length(unique(x$Site_troph))),
          dimnames = list(unique(x$Species),
                          unique(x$Time_step),
                          unique(x$Site_troph)))
  })
  
  resul.part <- lapply(list.array, function (x) var.partition(x))
  resul.part2 <- apply(do.call(rbind, resul.part), 2, mean)  
}


# Apply the rarefy function to the sec.ter.list

resul.raref <- lapply(sec.ter.list, function (x) {
  func_raref(x, 8, 1000)})

resul.raref

resul.raref2 <- do.call(rbind, resul.raref) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Meta_troph") 

resul.raref3 <- resul.raref2 %>% 
  mutate(New_tr_g =
           str_split(string = resul.raref2$Meta_troph, 
                     pattern = "_", simplify = T)[,2]) %>% 
  pivot_longer(cols = CV_S_L:phi_S2C_R, names_to = "Partition", 
               values_to = "Value_part") %>% 
    mutate(Sens_ana = "Rarefied") 

obs.raref <- meta.metr %>% 
  select(-Meta_clos:-Trophic_gr, -Ecosys, -Simp_gamma:-Precip_sync) %>% 
  pivot_longer(cols = c(CV_S_L:phi_S2C_R),
               names_to = "Partition", values_to = "Value_part") %>% 
  mutate(Sens_ana = "Full") %>% 
  right_join(select(resul.raref3, Meta_troph, Partition), 
             by = c("Meta_troph", "Partition")) %>% 
  rbind(resul.raref3) 
#

# Correlation between full and rarefied
obs.raref %>% 
  filter(!is.na(Value_part)) %>% 
  filter(Partition == "CV_S_L" |  Partition == "CV_C_L" | 
           Partition == "CV_C_R") %>% 
  pivot_wider(names_from = Sens_ana,
              values_from = Value_part) %>% 
  group_by(New_tr_g, Partition) %>% 
  summarise(corr = cor(Full, Rarefied)) %>% 
  select(corr) %>% 
  summary()

obs.raref %>% 
  filter(!is.na(Value_part)) %>% 
  filter(Partition == "phi_S2C_L" |  Partition == "phi_C_L2R") %>% 
  pivot_wider(names_from = Sens_ana,
              values_from = Value_part) %>% 
  group_by(New_tr_g, Partition) %>% 
  summarise(corr = cor(Full, Rarefied)) %>% 
  select(corr) %>% 
  summary()
  

#==================================
## Figures

col_meta_cv <- "#888888"
col_com_cv <- "#798234"
col_pop_cv <- "#834ba0"

col_syn_pop <- "#834ba0"
col_syn_com <- "#798234"

s7.1 <- obs.raref %>% 
  filter (!is.na(Sens_ana)) %>%
  filter (New_tr_g == "Tertiary") %>% 
  filter (Partition == "CV_S_L" | Partition == "CV_C_L" | Partition == "CV_C_R") %>%
  mutate (Partition = factor(Partition, levels = c("CV_S_L", "CV_C_L", "CV_C_R"))) %>%
  ggplot (aes(y=Value_part, x=Sens_ana)) +
  geom_dotplot (aes(fill = Partition, color = Partition), 
                alpha = 0.7, binaxis = "y", stackdir = "center", 
               position = "dodge", dotsize = 1.5, binwidth = 1/17) +
  geom_violin (aes(fill = Partition),
               trim = F, draw_quantiles = c(0.5), alpha = 0.6, scale = "count",
              bw = 0.2) +
  scale_fill_manual(name = "", 
                    labels = c("Population", "Community", "Metacommunity"),
                    values = c(col_pop_cv, col_com_cv, col_meta_cv)) +
  scale_color_manual(name = "", 
                     labels = c("Population", "Community", "Metacommunity"),
                     values = c(col_pop_cv, col_com_cv, col_meta_cv)) +
  xlab ("Sensitivity analysis data") +
  ggtitle("Tertiary consumers") +
  theme_classic() +
  theme(plot.title = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12), 
        legend.position = c(x=0.8, y=0.95),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.75,"line")) 

s7.2 <- obs.raref %>% 
  filter (!is.na(Sens_ana)) %>%
  filter (New_tr_g == "Secondary") %>% 
  filter (Partition == "CV_S_L" | Partition == "CV_C_L" | Partition == "CV_C_R") %>%
  mutate (Partition = factor(Partition, levels = c("CV_S_L", "CV_C_L", "CV_C_R"))) %>%
  ggplot (aes(y=Value_part, x=Sens_ana)) +
  geom_dotplot (aes(fill = Partition, color = Partition), 
                alpha = 0.7, binaxis = "y", stackdir = "center", 
                position = "dodge", dotsize = 1.5, binwidth = 1/17) +
  geom_violin (aes(fill = Partition),
               trim = F, draw_quantiles = c(0.5), alpha = 0.6, scale = "count",
               bw = 0.2) +
  scale_fill_manual(name = "", 
                    labels = c("Population", "Community", "Metacommunity"),
                    values = c(col_pop_cv, col_com_cv, col_meta_cv)) +
  scale_color_manual(name = "", 
                     labels = c("Population", "Community", "Metacommunity"),
                     values = c(col_pop_cv, col_com_cv, col_meta_cv)) +
  xlab ("Sensitivity analysis data") + ylab ("Temporal variability (CV)") +
  ggtitle("Secondary consumers") +
  theme_classic() +
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12), 
        legend.position = "none") 

s7.3 <- obs.raref %>%   
  filter (!is.na(Sens_ana)) %>%
  filter (New_tr_g == "Tertiary") %>%
  filter (Partition == "phi_S2C_L" | Partition == "phi_C_L2R") %>%
  mutate (Partition = factor(Partition, levels = c("phi_S2C_L", "phi_C_L2R"))) %>%
  ggplot (aes(y=Value_part, x=Sens_ana)) +
  geom_dotplot (aes(fill = Partition, color = Partition), 
                alpha = 0.7, binaxis = "y", stackdir = "center", 
                position = "dodge", dotsize = 1.5, binwidth = 1/17) +
  geom_violin (aes(fill = Partition),
               trim = F, draw_quantiles = c(0.5), alpha = 0.6, scale = "count",
               bw = 0.2) +
  scale_fill_manual(values = c(col_syn_pop, col_syn_com)) +
  scale_color_manual(values = c(col_syn_pop, col_syn_com)) +
  xlab ("Sensitivity analysis data") +
  ggtitle("Tertiary consumers") +
  theme_classic() +
  theme(plot.title = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        legend.position = "none")  

s7.4 <- obs.raref %>%   
  filter (!is.na(Sens_ana)) %>%
  filter (New_tr_g == "Secondary") %>%
  filter (Partition == "phi_S2C_L" | Partition == "phi_C_L2R") %>%
  mutate (Partition = factor(Partition, levels = c("phi_S2C_L", "phi_C_L2R"))) %>%
  ggplot (aes(y=Value_part, x=Sens_ana)) +
  geom_dotplot (aes(fill = Partition, color = Partition), 
                alpha = 0.7, binaxis = "y", stackdir = "center", 
                position = "dodge", dotsize = 1.5, binwidth = 1/17) +
  geom_violin (aes(fill = Partition),
               trim = F, draw_quantiles = c(0.5), alpha = 0.6, scale = "count",
               bw = 0.2) +
  scale_fill_manual(values = c(col_syn_pop, col_syn_com)) +
  scale_color_manual(values = c(col_syn_pop, col_syn_com)) +
  ylab("Synchrony") + xlab ("Sensitivity analysis data") +
  ggtitle("Secondary consumers") +
  theme_classic() +
  theme(plot.title = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12), 
        legend.position = "none") 

s7.2 + s7.1
ggsave("Output/Fig_rarefy1.pdf", width=24,height=12, units = "cm")

s7.4 + s7.3
ggsave("Output/Fig_rarefy2.pdf", width=24,height=12, units = "cm")
#===============================================

## Let's do something similar with number of years
# Because we need to the keep the time series continuous, there is no way
# we can rarefy them. We will just have to make some of them shorter.
# and the problem here is the opposite of what we have above - i.e.,
# producers and primary consumers have the longest time series.

# Let's explore the number of years in the data sets

full.data %>% 
  group_by(New_tr_g, Metacom) %>% 
  summarise(N_years = length(unique(Time_step))) %>% 
  group_by(New_tr_g) %>% 
  summarise(Max_years = max(N_years))

# Let's be conservative and try with only 11 years

full.data.red.series <- full.data %>% 
  filter (Time_step <= 11)
  

# I will use a list from now on as I think it is easier to use lapply than 
# various for loops, and I am still incompetent to purr functions.

long.median3 <- split(full.data.red.series, 
                      f = full.data.red.series$Meta_troph, drop = T)

length(long.median3)
names(long.median3) # here, each element is a unique combination of 
# a metacommunity and trophic group

## Arrange it the way Wang's function requires

long.median3 <- lapply(long.median3, function (x){
  arrange(x, Site_troph, Time_step, Species)})

## Create the arrays

all.array2 <- lapply(long.median3, function (x){ 
  array(data = x$Abundance,
        dim = c(length(unique(x$Species)), 
                length(unique(x$Time_step)), 
                length(unique(x$Site_troph))),
        dimnames = list(unique(x$Species),
                        unique(x$Time_step),
                        unique(x$Site_troph)))
})

#===============================
## Run partitioning analysis following Wang et al. 2019  

part.W2 <- lapply(all.array2, function (x) var.partition(x))
part.W2

# The name of each metric follows Wang et al. (2019)

## Organize results
resul.all2 <- Reduce(rbind, part.W2)
resul.all2

# Organize per metacom
df.info.all2 <- lapply(long.median3, function (x) {
  data.frame(Metacom = unlist(unlist(as.character(unique(na.omit(x$Metacom))))),
             Meta_troph = unlist(unlist(as.character(unique(na.omit(x$Meta_troph))))),
             Freq = unlist(unlist(as.character(unique(na.omit(x$Frequency))))),
             Time_step = unlist(unlist(as.numeric(unique(max(na.omit(x$Time_step)))))),
             Bio_gr = unlist(unlist(as.character(unique(na.omit(x$Bio_gr))))),
             Trophic_gr = unlist(unlist(as.character(unique(na.omit(x$Trophic_gr))))),
             New_tr_g = unlist(unlist(as.character(unique(na.omit(x$New_tr_g))))),
             Ecosys = unlist(unlist(as.character(unique(na.omit(x$Ecosys))))),
             nSites = unlist(unlist(as.numeric(length(unique(x$Site))))))})

df.info.all2

meta.metr.red <- as_tibble(data.frame(Reduce(rbind, df.info.all2), 
                                         resul.all2))
#

## Estimate correlations between full and reduced times series 
meta.metr.full <- read_csv("Input_data/meta_metr_figs.csv", T) %>% 
  select(Meta_troph, CV_S_L:phi_S2C_R) %>% 
  mutate(Full_red = "full")

meta.metr.red2 <- meta.metr.red %>% 
  select(Meta_troph, CV_S_L:phi_S2C_R) %>% 
  mutate(Full_red = "reduced")

meta.metr.full %>% 
  bind_rows(meta.metr.red2) %>% 
  pivot_longer(cols = CV_S_L:phi_S2C_R, names_to = "Partition", 
               values_to = "Value_part") %>% 
  filter(Partition == "CV_S_L" |  Partition == "CV_C_L" | 
                                                  Partition == "CV_C_R") %>% 
  pivot_wider(names_from = Full_red,
              values_from = Value_part) %>% 
  filter(!is.na(full), !is.na(reduced)) %>% 
  group_by(Partition) %>% 
  summarise(corr = cor(full, reduced)) 

meta.metr.full %>% 
  bind_rows(meta.metr.red2) %>% 
  pivot_longer(cols = CV_S_L:phi_S2C_R, names_to = "Partition", 
               values_to = "Value_part") %>% 
  filter(Partition == "phi_S2C_L" |  Partition == "phi_C_L2R") %>% 
  pivot_wider(names_from = Full_red,
              values_from = Value_part) %>% 
  filter(!is.na(full), !is.na(reduced)) %>% 
  group_by(Partition) %>% 
  summarise(corr = cor(full, reduced)) 


#===============================
# Prep figure with metacom partitions

## Order levels
meta.metr.red$New_tr_g <- factor(meta.metr.red$New_tr_g, 
                             levels = c("Producers","Primary",
                                        "Secondary", "Tertiary"))

# Move to long format
meta.metr.long <- 
  pivot_longer(meta.metr.red,
               cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R, phi_S_L2R, phi_C_L2R,
                        phi_S2C_L, phi_S2C_R),
               names_to = "Partition", values_to = "Value_part")
#
meta.metr.long$New_tr_g <- 
  factor(meta.metr.long$New_tr_g, 
         levels = c("Producers","Primary",
                    "Secondary", "Tertiary"))

with(meta.metr.long, unique(Partition))

# Choose colors

col_meta_cv <- "#888888"
col_com_cv <- "#798234"
col_pop_cv <- "#834ba0"

col_syn_pop <- "#834ba0"
col_syn_com <- "#798234"

#======
# Violin plots

s9 <- meta.metr.long %>%
  filter (Partition == "CV_S_L" | Partition == "CV_C_L" | Partition == "CV_C_R") %>%
  mutate (Partition = factor(Partition, levels = c("CV_S_L", "CV_C_L", "CV_C_R"))) %>%
  ggplot(aes(y=Value_part, x=New_tr_g)) +
  geom_dotplot(aes(fill = Partition, color = Partition),
               alpha = 0.7, binaxis = "y", stackdir = "center", 
               position = "dodge", dotsize = 3, binwidth = 1/17) +
  geom_violin(aes(fill = Partition), 
              trim = F, draw_quantiles = c(0.5), alpha = 0.6, scale = "count",
              bw = 0.2) +
  scale_fill_manual(name = "", 
                    labels = c("Population", "Community", "Metacommunity"),
                    values = c(col_pop_cv, col_com_cv, col_meta_cv)) +
  scale_color_manual(name = "", 
                     labels = c("Population", "Community", "Metacommunity"),
                     values = c(col_pop_cv, col_com_cv, col_meta_cv)) +
  ylab("Variability (CV)") + xlab ("") +
  ggtitle("A (reduced time-series)") +
  theme_classic() +
  theme(plot.title = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12), 
        legend.position = c(x=0.8, y=0.95),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.75,"line")) 

s10 <- meta.metr.long %>%   
  filter (Partition == "phi_S2C_L" | Partition == "phi_C_L2R") %>%
  mutate (Partition = factor(Partition, levels = c("phi_S2C_L", "phi_C_L2R"))) %>%
  ggplot(aes(y=Value_part, x=New_tr_g, fill = Partition, color = Partition)) +
  geom_dotplot(alpha = 0.7, binaxis = "y", stackdir = "center", 
               position = "nudge", dotsize = 3, binwidth = 1/50) +
  geom_split_violin(aes(y=Value_part, x=New_tr_g, fill = Partition), 
                    trim = F, draw_quantiles = c(0.5), alpha = 0.6, 
                    scale = "count", inherit.aes = F) +
  scale_fill_manual(values = c(col_syn_pop, col_syn_com)) +
  scale_color_manual(values = c(col_syn_pop, col_syn_com)) +
  ylab("Synchrony") + xlab ("Trophic levels") +
  scale_x_discrete(labels=c("Producers", "Primary", 
                            "Secondary", "Tertiary")) +
  ylim(c(0, 1.25)) +
  ggtitle("B (reduced time-series)") +
  theme_classic() +
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12), 
        legend.position = "none") 

s9/s10

ggsave("Output/Fig_red_time_series2.pdf", width=12,height=12, units = "cm")
