## Code for reproducing results in Siqueira et al. 2022
# "Ecological stability propagates across spatial scales and trophic 
# levels in freshwater ecosystems"
# Tadeu Siqueira, Rio Claro, 2022-05-01

# Here we prepare the figures and table we show in the main text. 

# Packages
library(tidyverse)
library(patchwork)
library(rcartocolor)
library(emmeans)
source("R_codes/Split_viol_plot.R")

## Code to partition variability and synchrony (Wang et al. 2019; Ecography)
source ("R_codes/var_part_wang.R") 

#===============================
## Data preparation
# Import data frame with results (stability partitions and diversity metrics)
# These were created with code: 
# "02_Siqueira_etal_SEM_analyses.R"

meta.metr <- read_csv("Input_data/meta_metr_figs.csv", T)
site.metr <- read_csv("Input_data/site_metr_figs.csv", T)

## Order levels
meta.metr$New_tr_g <- factor(meta.metr$New_tr_g, 
                               levels = c("Producers","Primary",
                                          "Secondary", "Tertiary"))

site.metr$New_tr_g <- factor(site.metr$New_tr_g, 
                            levels = c("Producers","Primary",
                                       "Secondary", "Tertiary"))

#===============================

## Figures
## Prep figure with metacom partitions

# Move to long format
meta.metr.long <- 
  pivot_longer(meta.metr,
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
display_carto_all(colorblind_friendly = T)

# e.g.:
carto_pal(12, "Safe")

col_meta_cv <- "#888888"
col_com_cv <- "#798234"
col_pop_cv <- "#834ba0"

col_syn_pop <- "#834ba0"
col_syn_com <- "#798234"

#======
# Violin plots

a1 <- meta.metr.long %>%
  filter (Partition == "CV_S_L" | Partition == "CV_C_L" | Partition == "CV_C_R") %>%
  mutate (Partition = factor(Partition, levels = c("CV_S_L", "CV_C_L", "CV_C_R"))) %>%
  ggplot(aes(y=Value_part, x=New_tr_g, fill = Partition)) +
  geom_dotplot(alpha = 0.7, binaxis = "y", stackdir = "center", 
               position = "dodge", dotsize = 3, binwidth = 1/17, aes(color = Partition)) +
  geom_violin(trim = F, draw_quantiles = c(0.5), alpha = 0.7, scale = "count",
              bw = 0.2) +
  scale_fill_manual(name = "", 
                    labels = c("Population", "Community", "Metacommunity"),
                    values = c(col_pop_cv, col_com_cv, col_meta_cv)) +
  scale_color_manual(name = "", 
                    labels = c("Population", "Community", "Metacommunity"),
                    values = c(col_pop_cv, col_com_cv, col_meta_cv)) +
  ylab("Variability (CV)") + xlab ("") +
  ggtitle("(a)") +
  theme_classic() +
  theme(plot.title = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12), 
        legend.position = c(x=0.8, y=0.95),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.75,"line")) +
  annotate("text", x=1, y=2.75, label= "a1",
           col="#798234", size=3.5, parse=TRUE) +
  annotate("text", x=2, y=3.4, label= "a1",
           col="#798234", size=3.5, parse=TRUE) +
  annotate("text", x=1.3, y=1.95, label= "a2",
           col="#888888", size=3.5, parse=TRUE) +
  annotate("text", x=2.3, y=3.12, label= "a2",
           col="#888888", size=3.5, parse=TRUE) +
  annotate("text", x=3.3, y=1.76, label= "a3",
           col="#888888", size=3.5, parse=TRUE) +
  annotate("text", x=4.3, y=1.4, label= "a3",
           col="#888888", size=3, parse=TRUE)
  
  

a2 <- meta.metr.long %>%   
  filter (Partition == "phi_S2C_L" | Partition == "phi_C_L2R") %>%
  mutate (Partition = factor(Partition, levels = c("phi_S2C_L", "phi_C_L2R"))) %>%
  ggplot(aes(y=Value_part, x=New_tr_g, fill = Partition)) +
  geom_dotplot(alpha = 0.7, binaxis = "y", stackdir = "center", 
               position = "nudge", dotsize = 3, binwidth = 1/50, aes(color = Partition)) +
  geom_split_violin(trim = F, draw_quantiles = c(0.5), alpha = 0.7, 
                    scale = "count") +
  scale_fill_manual(values = c(col_syn_pop, col_syn_com)) +
  scale_color_manual(values = c(col_syn_pop, col_syn_com)) +
  ylab("Synchrony") + xlab ("Trophic levels") +
  scale_x_discrete(labels=c("Producers", "Primary", 
                            "Secondary", "Tertiary")) +
  ylim(c(0, 1.25)) +
  ggtitle("(b)") +
  theme_classic() +
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12), 
        legend.position = "none") +
  annotate("text", x=1.9, y=1.1, label= "b1",
           col="#834ba0", size=3.5, parse=TRUE) +
  annotate("text", x=2.9, y=1.1, label= "b1",
           col="#834ba0", size=3.5, parse=TRUE) +
  annotate("text", x=1.25, y=0.9, label= "b2 b3 b4",
           col="#798234", size=3.5, parse=F) +
  annotate("text", x=2.1, y=1.1, label= "b2",
           col="#798234", size=3.5, parse=TRUE) +
  annotate("text", x=3.2, y=0.9, label= "b3 b5",
           col="#798234", size=3.5, parse=F) +
  annotate("text", x=4.2, y=0.75, label= "b4 b5",
           col="#798234", size=3.5, parse=F)
 
  
  

a1/a2

ggsave("Output/Fig1.pdf", width=18,height=18, units = "cm")

# Test the differences suggested in the plot

df.cvs <- meta.metr.long %>%
  filter (Partition == "CV_S_L" | Partition == "CV_C_L" | Partition == "CV_C_R") %>%
  mutate (Partition = factor(Partition, levels = c("CV_S_L", "CV_C_L", "CV_C_R")))

mod.cvs <- lm(log(Value_part) ~ Partition*New_tr_g, data = df.cvs)
plot(mod.cvs) # They look OK
anova(mod.cvs) # No interaction
summary(mod.cvs)

# Pair-wise specific comparisons
# I am not interested in comparing Mv vs Cv, or Cv vs Pv
# I am interested in comparing each partition among trophic levels

t1 <- df.cvs %>%
  filter(Partition == "CV_C_R") %>%
  rstatix::emmeans_test(log(Value_part) ~ New_tr_g, p.adjust.method = "holm",
                        detailed = T) #%>% 
  #get_emmeans()

t2 <- df.cvs %>%
  filter(Partition == "CV_C_L") %>%
  rstatix::emmeans_test(log(Value_part) ~ New_tr_g, p.adjust.method = "holm",
                        detailed = T)

t3 <- df.cvs %>%
  filter(Partition == "CV_S_L") %>%
  rstatix::emmeans_test(log(Value_part) ~ New_tr_g, p.adjust.method = "holm",
                        detailed = T)

write.table(x=tibble(bind_rows(t1,t2,t3), 
                     part = rep(c("CV_C_R", "CV_C_L", "CV_S_L"), each = 6)), 
            file = "Output/cvs_pair_comps.csv", sep = ",", 
            quote = FALSE, row.names = F)

### Do the same with synchrony
df.sync <- meta.metr.long %>%   
  filter (Partition == "phi_S2C_L" | Partition == "phi_C_L2R") %>%
  mutate (Partition = factor(Partition, levels = c("phi_S2C_L", "phi_C_L2R")))

mod.syn <- lm(log(Value_part) ~ Partition*New_tr_g, data = df.sync)
plot(mod.syn) # Not great
anova(mod.syn) # There is an interaction
summary(mod.syn)

# Pair-wise specific comparisons

t4 <- df.sync %>%
  filter(Partition == "phi_S2C_L") %>%
  rstatix::emmeans_test(log(Value_part) ~ New_tr_g, p.adjust.method = "holm",
                        detailed = T) 

t5 <- df.sync %>%
  filter(Partition == "phi_C_L2R") %>%
  rstatix::emmeans_test(log(Value_part) ~ New_tr_g, p.adjust.method = "holm",
                        detailed = T)


write.table(x=tibble(bind_rows(t4,t5), 
                     part = rep(c("phi_S2C_L", "phi_C_L2R"), each = 6)), 
            file = "Output/syns_pair_comps.csv", sep = ",", 
            quote = FALSE, row.names = F)

#statistic: Test statistic (t.ratio) used to compute the p-value
#estimate: estimate of the effect size, that is the difference between 
#the two emmeans (estimated marginal means).

meta.metr %>% 
  mutate(Pcv_Ccv = CV_C_L/CV_S_L, Ccv_Mcv = CV_C_R/CV_C_L) %>% 
  select(New_tr_g, Pcv_Ccv, Ccv_Mcv) %>% 
  group_by(New_tr_g) %>% 
  summarise(m1 = mean(Pcv_Ccv), m2 = mean(Ccv_Mcv))

#===============================

# Figure 2 is based on a subset of the whole data
# I will show how to create the subset first and then prepare the figure

## Import the complete data set 
# This is the same we used in script 01
# Each row here is the abundance of a given species in a given site and time
# Various other columns describe these row (lat, long, sample id, etc.) 

# This data table has been pre-processed to replace NAs in species abundances
# by the median of the abundance of that species in other surveys.

full.bio.data <- as_tibble(
  read.csv("Input_data/data_stab_analysis_except_lepas.csv", T))
# Not sure why read_csv read the data weirdly

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

full.data.no.ponds <- full.data1 %>% 
  filter (Metacom != "ELA_1", Metacom != "ELA_2", Metacom != "ELA_3", 
          Metacom != "NTL", Metacom != "M_Hill")

with(full.data.no.ponds, unique(Metacom))
with(full.data.no.ponds, unique(Site))

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

# Find out which data sets include info on multiple trophic levels
# I know beforehand that none include info on producers to tertiary consumers 
# I found out that only two data sets include info on primary to tertiary cons. 

# So, I decided to gather data sets that have primary and secondary, 
# secondary and tertiary, primary, secondary, tertiary 
# and analyze them separately

# This first one is special as these are data sets sampled in the same 
# region, but with different methods for phyto and zoops.

prim.sec <- full.data2 %>% 
  group_by (Metacom) %>% 
  filter(all(c('Primary', 'Secondary') %in% New_tr_g) & 
           all(New_tr_g %in% c('Primary', 'Secondary'))) %>% 
  ungroup () %>% 
  mutate(Mult_tr_g = "prim_sec") 

sec.ter <- full.data2 %>% 
  group_by (Metacom) %>% 
  filter(all(c('Secondary', 'Tertiary') %in% New_tr_g) & 
           all(New_tr_g %in% c('Secondary', 'Tertiary'))) %>% 
  ungroup () %>% 
  mutate(Mult_tr_g = "sec_ter") 

prim.sec.ter <- full.data2 %>% 
  group_by (Metacom) %>% 
  filter(all(c('Primary', 'Secondary', 'Tertiary') %in% New_tr_g) & 
           all(New_tr_g %in% c('Primary', 'Secondary', 'Tertiary'))) %>% 
  ungroup () %>% 
  mutate(Mult_tr_g = "prim_sec_ter") 

mult.tr <- bind_rows(prim.sec, sec.ter, prim.sec.ter)

#=============

# I will use a list for subsequent analyses

long.median2 <- split(mult.tr, f = mult.tr$Meta_troph, drop = T)

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

#===============================
## Run partitioning analysis following Wang et al. 2019  

part.W <- lapply(all.array, function (x) var.partition(x))
part.W

# The name of each metric follows Wang et al. (2019)

## Organize results
resul.all <- Reduce(rbind, part.W)
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
             Mult_tr_g = unlist(unlist(as.character(unique(na.omit(x$Mult_tr_g))))),
             nSites = unlist(unlist(as.numeric(length(unique(x$Site))))))})

df.info.all

meta.metr <- as_tibble(data.frame(Reduce(rbind, df.info.all), 
                                  resul.all))

#===========================================
## Figure 2
## Order levels
meta.metr$New_tr_g <- factor(meta.metr$New_tr_g, 
                             levels = c("Primary",
                                        "Secondary", "Tertiary"))

# Move to long format
meta.metr.long <- 
  pivot_longer(meta.metr,
               cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R, phi_S_L2R, phi_C_L2R,
                        phi_S2C_L, phi_S2C_R),
               names_to = "Partition", values_to = "Value_part")
#
meta.metr.long$New_tr_g <- 
  factor(meta.metr.long$New_tr_g, 
         levels = c("Primary",
                    "Secondary", "Tertiary"))
#===========================================

### Variability

p.all <- meta.metr.long %>%
  filter (Partition == "CV_S_L") %>%
  ggplot(aes(y=Value_part, x=New_tr_g)) +
  geom_point(alpha = 0.7, size = 3, color = col_pop_cv) +
  geom_line(aes(group = Metacom, linetype = Mult_tr_g), color = "darkgrey",
            alpha = .6) +
  scale_linetype_manual(values = c(2,3,1)) +
  ylab("Variability (CV)") + xlab ("Trophic level") +
  ggtitle("(a)", subtitle = "Population") +
  ylim (0.10,3) + # change this accordingly
  theme_classic() +
  theme(plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 12, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12), 
        axis.title.x = element_blank(),
        legend.position = "none") +
  annotate("text", x=1, y=2.9, label= "a1",
           col="#834ba0", size=3.5, parse=TRUE) +
  annotate("text", x=2, y=1.95, label= "a1",
           col="#834ba0", size=3.5, parse=TRUE) 

# Paired t-tests
t1 <- t.test (Value_part ~ New_tr_g, 
              data = meta.metr.long %>%
                filter (Partition == "CV_S_L", Mult_tr_g == "sec_ter"), 
              paired = TRUE, 
              alternative = "greater")
t1$stderr

t1.1 <- t.test (Value_part ~ New_tr_g, 
              data = meta.metr.long %>%
                filter (Partition == "CV_S_L", Mult_tr_g == "prim_sec"), 
              paired = TRUE, 
              alternative = "greater")

t1.1$stderr

# Just for curiosity, I tried a mixed-effects model and got the same t value
# Here is the reason:
# https://medium.com/@marco_laube/the-paired-t-test-and-linear-mixed-models-185a084d7813

#t1.2 <- lme4::lmer(Value_part ~ New_tr_g + (1|Metacom), 
             #data = meta.metr.long %>%
               #filter (Partition == "CV_S_L", Mult_tr_g == "sec_ter"))
#summary(t1.1)
#

c.all <- meta.metr.long %>%
  filter (Partition == "CV_C_L") %>%
  ggplot(aes(y=Value_part, x=New_tr_g, fill = New_tr_g, color = New_tr_g)) +
  geom_point(alpha = 0.7, size = 3, color = col_com_cv) +
  geom_line(aes(group = Metacom, linetype = Mult_tr_g), color = "darkgrey",
            alpha = .6) +
  scale_linetype_manual(values = c(2,3,1)) +
  xlab ("Trophic level") +
  ggtitle(label = "", subtitle = "Community") +
  ylim (0.10,3) +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 12, hjust=0.5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  annotate("text", x=1, y=2.8, label= "a2",
           col="#798234", size=3.5, parse=TRUE) +
  annotate("text", x=2, y=1.6, label= "a2",
           col="#798234", size=3.5, parse=TRUE)

#
t2 <- t.test (Value_part ~ New_tr_g, 
              data = meta.metr.long %>%
                filter (Partition == "CV_C_L", Mult_tr_g == "sec_ter"), 
              paired = TRUE, 
              alternative = "greater")

t2$stderr

t2.1 <- t.test (Value_part ~ New_tr_g, 
              data = meta.metr.long %>%
                filter (Partition == "CV_C_L", Mult_tr_g == "prim_sec"), 
              paired = TRUE, 
              alternative = "greater")

t2.1$stderr

#

mc.all <- meta.metr.long %>%
  filter (Partition == "CV_C_R") %>%
  ggplot(aes(y=Value_part, x=New_tr_g, fill = New_tr_g, color = New_tr_g)) +
  geom_point(alpha = 0.7, size = 3, color = col_meta_cv) +
  geom_line(aes(group = Metacom, linetype = Mult_tr_g), color = "darkgrey",
            alpha = .6) +
  scale_linetype_manual(values = c(2,3,1)) +
  xlab ("Trophic level") +
  ggtitle(label = "", subtitle = "Metacommunity") +
  ylim (0.10,3) +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 12, hjust=0.5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  annotate("text", x=1, y=2.5, label= "a3",
           col="#888888", size=3.5, parse=TRUE) +
  annotate("text", x=2, y=1.3, label= "a3 a4",
           col="#888888", size=3.5, parse=F) +
  annotate("text", x=3, y=0.8, label= "a4",
           col="#888888", size=3.5, parse=T)
  
  
#
t3 <- t.test (Value_part ~ New_tr_g, 
              data = meta.metr.long %>%
                filter (Partition == "CV_C_R", Mult_tr_g == "sec_ter"), 
              paired = TRUE, 
              alternative = "greater")

t3$stderr

t3.1 <- t.test (Value_part ~ New_tr_g, 
              data = meta.metr.long %>%
                filter (Partition == "CV_C_R", Mult_tr_g == "prim_sec"), 
              paired = TRUE, 
              alternative = "greater")

t3.1$stderr
#

p.all + c.all + mc.all


# Synchrony components

psy.prim <- meta.metr.long %>%
  filter (New_tr_g == "Primary") %>%
  filter (Partition == "phi_S2C_L" | Partition == "phi_C_L2R") %>%
  mutate (Partition = factor(Partition, levels = c("phi_S2C_L", "phi_C_L2R"))) %>%
  ggplot(aes(y=Value_part, x=Partition, color = Partition)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_line(aes(group = Metacom, linetype = Mult_tr_g), color = "darkgrey",
            alpha = .6) +
  ylab("Synchrony") + xlab ("Organizational level") +
  ggtitle("(b)", subtitle = "Primary") +
  ylim(c(0, 1.25)) +
  scale_x_discrete(labels=c("Population", "Community")) +
  scale_color_manual(values = c(col_syn_pop, col_syn_com)) +
  scale_linetype_manual(values = c(2,3)) +
  theme_classic() +
  theme(plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 12, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  annotate("text", x=1, y=1.1, label= "b1",
           col="#834ba0", size=3.5, parse=TRUE) +
  annotate("text", x=2, y=1.1, label= "b1",
         col="#798234", size=3.5, parse=TRUE)
  
#
# Paired t-test
t4 <- t.test (Value_part ~ Partition, 
              data = meta.metr.long %>%
                filter (New_tr_g == "Primary") %>%
                filter (Partition == "phi_S2C_L" | Partition == "phi_C_L2R"), 
              paired = TRUE, 
              alternative = "less")

t4$stderr
#

psy.sec <- meta.metr.long %>%
  filter (New_tr_g == "Secondary") %>%
  filter (Partition == "phi_S2C_L" | Partition == "phi_C_L2R") %>%
  mutate (Partition = factor(Partition, levels = c("phi_S2C_L", "phi_C_L2R"))) %>%
  ggplot(aes(y=Value_part, x=Partition, color = Partition)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_line(aes(group = Metacom, linetype = Mult_tr_g), color = "darkgrey",
            alpha = .6) +
  ylab("Synchrony") + xlab ("Organizational level") +
  ggtitle(label = "", subtitle = "Secondary") +
  ylim(c(0, 1.25)) +
  scale_x_discrete(labels=c("Population", "Community")) +
  scale_color_manual(values = c(col_syn_pop, col_syn_com)) +
  scale_linetype_manual(values = c(2,3,1)) +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 12, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

#
t5 <- t.test (Value_part ~ Partition, 
              data = meta.metr.long %>%
                filter (New_tr_g == "Secondary") %>%
                filter (Partition == "phi_S2C_L" | Partition == "phi_C_L2R"), 
              paired = TRUE, 
              alternative = "less")

t5$stderr
#

psy.ter <- meta.metr.long %>%
  filter (New_tr_g == "Tertiary") %>%
  filter (Partition == "phi_S2C_L" | Partition == "phi_C_L2R") %>%
  mutate (Partition = factor(Partition, levels = c("phi_S2C_L", "phi_C_L2R"))) %>%
  ggplot(aes(y=Value_part, x=Partition, color = Partition)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_line(aes(group = Metacom, linetype = Mult_tr_g), color = "darkgrey",
            alpha = .6) +
  xlab ("Organizational level") +
  ggtitle(label = "", subtitle = "Tertiary") +
  ylim(c(0, 1.25)) +
  scale_x_discrete(labels=c("Population", "Community")) +
  scale_color_manual(values = c(col_syn_pop, col_syn_com)) +
  scale_linetype_manual(values = c(3,1)) +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 12, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

#
t6 <- t.test (Value_part ~ Partition, 
              data = meta.metr.long %>%
                filter (New_tr_g == "Tertiary") %>%
                filter (Partition == "phi_S2C_L" | Partition == "phi_C_L2R"), 
              paired = TRUE, 
              alternative = "less")

t6$stderr
#

psy.prim + psy.sec + psy.ter

(p.all + c.all + mc.all) /
  (psy.prim + psy.sec + psy.ter)

ggsave("Output/NewFig_2.pdf", width=21,height=18, units = "cm")

# Table with Paired t-tests results

write.table(x=tibble(bind_rows(
  c(t1$statistic, p = t1.1$p.value, t1.1$estimate, t1.1$parameter),
  c(t1$statistic, p = t1$p.value, t1$estimate, t1$parameter),
  c(t1$statistic, p = t2.1$p.value, t2.1$estimate, t2.1$parameter),
  c(t2$statistic, p = t2$p.value, t2$estimate, t2$parameter),
  c(t1$statistic, p = t3.1$p.value, t3.1$estimate, t3.1$parameter),
  c(t3$statistic, p = t3$p.value, t3$estimate, t3$parameter),
  c(t4$statistic, p = t4$p.value, t4$estimate, t4$parameter),
  c(t5$statistic, p = t5$p.value, t5$estimate, t5$parameter),
  c(t6$statistic, p = t6$p.value, t6$estimate, t6$parameter)),
  Comparison = c(rep("Variability (CV)", 6), rep("Synchrony", 3)),
  "Paired comparison" = c(rep("Primary vs. Secondary", 3),
                          rep("Secondary vs. Tertiary", 3), 
                          rep("Population vs. Community", 3)),
  Condition = c(rep(c("Population", "Community", "Metacommunity"), 2), 
                c("Primary", "Secondary", "Tertiary"))), 
  file = "Output/paited_t_tests.csv", sep = ",", 
  quote = FALSE, row.names = F)
