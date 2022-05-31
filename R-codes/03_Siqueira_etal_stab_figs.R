## Code for reproducing results in Siqueira et al. (2022):
# "Ecological stability propagates across spatial scales and trophic 
# levels in freshwater ecosystems"
# https://ecoevorxiv.org/mpf5x
# Tadeu Siqueira, Rio Claro, 2022-05-31

# Here we prepare the figures and table we show in the main text. 

# Packages
library(tidyverse)
library(patchwork)
library(rcartocolor)
library(emmeans)
library(here)
source(here("R-codes", "Split_viol_plot.R"))

#===============================
## Data preparation
# Import data frame with results (stability partitions and diversity metrics)
# These were created with code: 
# "02_Siqueira_etal_SEM_analyses.R"

meta.metr <- read_csv(here("Input_data", "meta_metr_figs.csv"), T)
site.metr <- read_csv(here("Input_data", "site_metr_figs.csv"), T)

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
  ggplot(aes(y=Value_part, x=New_tr_g, fill = Partition, color = Partition)) +
  geom_dotplot(alpha = 0.7, binaxis = "y", stackdir = "center", 
               position = "dodge", dotsize = 3, binwidth = 1/17) +
  geom_violin(trim = F, draw_quantiles = c(0.5), alpha = 0.7, scale = "count",
              bw = 0.2) +
  scale_fill_manual(name = "", 
                    labels = c("Population", "Community", "Metacommunity"),
                    values = c(col_pop_cv, col_com_cv, col_meta_cv)) +
  scale_color_manual(name = "", 
                    labels = c("Population", "Community", "Metacommunity"),
                    values = c(col_pop_cv, col_com_cv, col_meta_cv)) +
  ylab("Variability (CV)") + xlab ("") +
  ggtitle("A") +
  theme_classic() +
  theme(plot.title = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12), 
        legend.position = c(x=0.8, y=0.95),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.75,"line")) 

a2 <- meta.metr.long %>%   
  filter (Partition == "phi_S2C_L" | Partition == "phi_C_L2R") %>%
  mutate (Partition = factor(Partition, levels = c("phi_S2C_L", "phi_C_L2R"))) %>%
  ggplot(aes(y=Value_part, x=New_tr_g, fill = Partition, color = Partition)) +
  geom_dotplot(alpha = 0.7, binaxis = "y", stackdir = "center", 
               position = "nudge", dotsize = 3, binwidth = 1/50) +
  geom_split_violin(trim = F, draw_quantiles = c(0.5), alpha = 0.7, 
                    scale = "count") +
  scale_fill_manual(values = c(col_syn_pop, col_syn_com)) +
  scale_color_manual(values = c(col_syn_pop, col_syn_com)) +
  ylab("Synchrony") + xlab ("Trophic levels") +
  scale_x_discrete(labels=c("Producers", "Primary", 
                            "Secondary", "Tertiary")) +
  ylim(c(0, 1.25)) +
  ggtitle("B") +
  theme_classic() +
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12), 
        legend.position = "none") 

a1/a2

ggsave(here("Figures_tables", "Fig2.pdf"), width=12,height=12, units = "cm")

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
  rstatix::emmeans_test(log(Value_part) ~ New_tr_g, p.adjust.method = "holm") #%>% 
  #get_emmeans()

t2 <- df.cvs %>%
  filter(Partition == "CV_C_L") %>%
  rstatix::emmeans_test(log(Value_part) ~ New_tr_g, p.adjust.method = "holm")

t3 <- df.cvs %>%
  filter(Partition == "CV_S_L") %>%
  rstatix::emmeans_test(log(Value_part) ~ New_tr_g, p.adjust.method = "holm")

write.table(x=tibble(bind_rows(t1,t2,t3), 
                     part = rep(c("CV_C_R", "CV_C_L", "CV_S_L"), each = 6)), 
            file = here("Figures_tables", "cvs_pair_comps.csv"), sep = ",", 
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
  rstatix::emmeans_test(log(Value_part) ~ New_tr_g, p.adjust.method = "holm") 

t5 <- df.sync %>%
  filter(Partition == "phi_C_L2R") %>%
  rstatix::emmeans_test(log(Value_part) ~ New_tr_g, p.adjust.method = "holm")


write.table(x=tibble(bind_rows(t4,t5), 
                     part = rep(c("phi_S2C_L", "phi_C_L2R"), each = 6)), 
            file = here("Figures_tables", "syns_pair_comps.csv"), sep = ",", 
            quote = FALSE, row.names = F)

#statistic: Test statistic (t.ratio) used to compute the p-value
#estimate: estimate of the effect size, that is the difference between 
#the two emmeans (estimated marginal means).

## Now site and metacom metrics
# CV local
# global plot

# Re-organize trophic levels for the sake of figure legends
site.metr$New_tr_g <- factor(
  ifelse(site.metr$New_tr_g == "Primary", "Primary consumer",
      ifelse(site.metr$New_tr_g == "Secondary", "Secondary consumer",
             ifelse(site.metr$New_tr_g == "Tertiary", "Tertiary consumer",
                    "Producer"))))

# Reverse the trophic level order as follow
site.metr$New_tr_g <- factor(site.metr$New_tr_g, 
                              levels = c("Tertiary consumer",
                                         "Secondary consumer",
                                         "Primary consumer",
                                         "Producer"))
# Set colors
col_prod <- "#661100"
col_prim <- "#f7945d" 
col_sec <- "#6699CC" 
col_ter <- "#117733"

## Metacom scale
# Reverse the trophic level order as follow
meta.metr$New_tr_g <- factor(
  ifelse(meta.metr$New_tr_g == "Primary", "Primary consumer",
         ifelse(meta.metr$New_tr_g == "Secondary", "Secondary consumer",
                ifelse(meta.metr$New_tr_g == "Tertiary", "Tertiary consumer",
                       "Producer"))))

meta.metr$New_tr_g <- factor(meta.metr$New_tr_g, 
                              levels = c("Tertiary consumer",
                                         "Secondary consumer",
                                         "Primary consumer",
                                         "Producer"))
## Figures

p7 <- ggplot(meta.metr, 
             aes(x = phi_S2C_L, y = phi_C_L2R, 
                 color = New_tr_g) ) +
  geom_point(size = 3, alpha = 0.7, aes(shape = New_tr_g, fill = New_tr_g)) + 
  ylab ("Community spatial synchrony") +
  ylim (0, 1) + xlim (0,1) + 
  geom_abline(intercept = 0, color = "darkgrey") +
  xlab ("Population synchrony") +
  scale_shape_manual(values=c(15:17, 23), 
                     name = "Trophic levels") +
  scale_color_manual(values=c(col_prod, col_prim, col_sec, col_ter), 
                     name = "Trophic levels") +
  scale_fill_manual(values=c(col_prod, col_prim, col_sec, col_ter), 
                     name = "Trophic levels") +
  ggtitle("B")+
  theme_classic() +
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = "none") 

p7  

p8 <- ggplot(meta.metr, 
             aes(x = CV_S_L, y = CV_C_L, 
                 color = New_tr_g) ) +
  geom_point(size = 3, alpha = 0.7, aes(shape = New_tr_g, fill = New_tr_g)) + 
  ylab ("Community variability") +
  ylim (0, 3) + xlim (0, 3) + 
  geom_abline(intercept = 0, color = "darkgrey") +
  xlab ("Population variability") +
  scale_shape_manual(values=c(15:17, 23), 
                     name = "Trophic levels") +
  scale_color_manual(values=c(col_prod, col_prim, col_sec, col_ter), 
                     name = "Trophic levels") +
  scale_fill_manual(values=c(col_prod, col_prim, col_sec, col_ter), 
                    name = "Trophic levels") +
  ggtitle("A") +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = c(x=0.30, y=0.85),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10, margin = margin(t = 1)),
        legend.spacing.y = unit(0.01, 'cm'),
        legend.spacing.x = unit(0.01, 'cm'),
        plot.title = element_text(size = 12)) 
  
p8
#

p9 <- ggplot(meta.metr, 
             aes(x = CV_C_L, y = CV_C_R, 
                 color = New_tr_g) ) +
  geom_point(size = 3, alpha = 0.7, aes(shape = New_tr_g, fill = New_tr_g)) + 
  ylab ("Metacommunity variability") +
  ylim (0, 3) + xlim (0, 3) + 
  geom_abline(intercept = 0, color = "darkgrey") +
  xlab ("Community variability") +
  scale_shape_manual(values=c(15:17, 23), 
                     name = "Trophic levels") +
  scale_color_manual(values=c(col_prod, col_prim, col_sec, col_ter), 
                     name = "Trophic levels") +
  scale_fill_manual(values=c(col_prod, col_prim, col_sec, col_ter), 
                    name = "Trophic levels") +
  ggtitle ("C") +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(size = 12))

p8 / p7 / p9
ggsave(here("Figures_tables", "Fig3.pdf"), width=12,height=22, units = "cm")

#===============================
## Table 1

meta.metr %>% 
  mutate(Pcv_Ccv = CV_C_L/CV_S_L, Ccv_Mcv = CV_C_R/CV_C_L,
         Pcv_Mcv = CV_C_R/CV_S_L) %>% 
  select(New_tr_g, Pcv_Ccv, Ccv_Mcv, Pcv_Mcv) %>% 
  group_by(New_tr_g) %>% 
  summarise(m1 = mean(Pcv_Ccv), m2 = mean(Ccv_Mcv),
            m3 = mean(Pcv_Mcv))

### End
#===============================

