## Code for reproducing results in Siqueira et al. (2022):
# "Ecological stability propagates across spatial scales and trophic 
# levels in freshwater ecosystems"
# Tadeu Siqueira, Rio Claro, 2022-05-01

# Here we prepare the figures we show as supp. information. 

# Packages
library(tidyverse)
library(patchwork)
library(ggdist)
library(sf)
library(spData)        # load geographic data
library(spDataLarge)   # load larger geographic data
library(tmap)

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
#===============================
## Figures as supp. material

# Map

# options
tmap_options(check.and.fix = TRUE)

# vector
occ_vect <- site.metr %>% 
  select (Metacom, Lat, Long) %>% 
  group_by(Metacom) %>% 
  summarise(Lat = mean(Lat), Long = mean(Long)) %>% 
  sf::st_as_sf(coords = c("Long", "Lat"), crs = 4326) 

occ_vect

occ_vect_moll <- sf::st_transform(occ_vect, crs = "+proj=moll")
occ_vect_moll

# world
world_moll <- sf::st_transform(world, crs = "+proj=moll")
world_moll

# map option 1

world_mollweide_gr <- st_graticule(lat = c(-89.9, seq(-80, 40, 20), 89.9)) %>%
  lwgeom::st_transform_proj(crs = "+proj=moll")

map1 <- tm_shape(world_mollweide_gr) +
  tm_lines(col = "grey") +
  tm_shape(world_moll) +
  tm_polygons() +
  tm_shape(occ_vect_moll) +
  tm_dots(col = "black", size = .2, alpha = .5) +
  tm_layout(frame = FALSE)
map1

# save
tmap_save(tm = map1, filename = "Output/map1.png", wi = 30, he = 20, 
          un = "cm", dpi = 300)

# map option 2

map2 <- tm_shape(world, 
                 bbox = st_bbox(c(xmin = -130, xmax = 40, ymax = 90, ymin = -60), 
                                crs = st_crs(4326))) +
  tm_graticules() +
  tm_polygons() +
  tm_shape(occ_vect) +
  tm_dots(size = .4, alpha = .45) +
  tm_compass(position = c("left", "bottom")) +
  tm_scale_bar(position = c("left", "bottom"), breaks = c(0, 1000, 2000, 3000))
map2

# save
tmap_save(tm = map2, filename = "Output/map2.png", wi = 9, he = 8, 
          un = "in", dpi = 300)


# End map prep. --------------------------------------------------------------
#
## Figures to show that the propagation of variability and synchrony across
# trophic levels did not depend on ecosystem type (lotic vs lentic)

meta.metr1.supp <- meta.metr %>%
  mutate(New_ecos_tr = paste0(Ecosys, "_", New_tr_g)) %>%
  mutate(New_ecos_tr = case_when(
    New_ecos_tr == "Stream_Primary consumer" ~ "Lot_PriCon",
    New_ecos_tr == "Stream_Secondary consumer" ~ "Lot_SecCon",
    New_ecos_tr == "Stream_Tertiary consumer" ~ "Lot_TerCon",
    New_ecos_tr == "Stream_Producer" ~ "Lot_Prod",
    New_ecos_tr == "Lake_Producer" ~ "Len_Prod",
    New_ecos_tr == "Lake_Primary consumer" ~ "Len_PriCon")) %>% 
  mutate(New_ecos_tr = factor(New_ecos_tr,
                              levels = rev(c("Lot_TerCon", "Lot_SecCon",
                                             "Lot_PriCon", "Lot_Prod",
                                             "Len_PriCon", "Len_Prod"))))

ps1 <- meta.metr1.supp %>% 
ggplot(aes(y = New_ecos_tr, x = log(CV_C_R), fill = New_tr_g)) +
  stat_slab(aes(thickness = stat(pdf*n)), scale = 0.7) +
  stat_dotsinterval(side = "bottom", scale = 0.7, slab_size = NA,
                    dotsize = 2, alpha = 0.75) +
  ylab ("Trophic level per ecosystem type") +
  xlab ("Metacommunity variability (Mcv)") +
  ggtitle ("A") +
  theme_classic() +
  theme(axis.text = element_text(size = 17),
        axis.title = element_text(size = 20),
        legend.position = c(x=0.75, y=0.9),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18)) 

#

ps2 <- meta.metr1.supp %>% 
  ggplot(aes(y = New_ecos_tr, x = log(CV_C_L), fill = New_tr_g)) +
  stat_slab(aes(thickness = stat(pdf*n)), scale = 0.6) +
  stat_dotsinterval(side = "bottom", scale = 0.7, slab_size = NA,
                    dotsize = 2, alpha = 0.75) +
  ylab ("") +
  xlab ("Community variability (Ccv)") +
  theme_classic() +
  ggtitle ("B") +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size = 17),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_blank(),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18))
#

ps3 <- meta.metr1.supp %>% 
  ggplot(aes(y = New_ecos_tr, x = log(CV_S_L), fill = New_tr_g)) +
  stat_slab(aes(thickness = stat(pdf*n)), scale = 0.7) +
  stat_dotsinterval(side = "bottom", scale = 0.7, slab_size = NA,
                    dotsize = 2, alpha = 0.75) +
  ylab ("") +
  xlab ("Population variability (Ccv)") +
  theme_classic() +
  ggtitle ("C") +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size = 17),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_blank(),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18))
#
(ps1 + ps2 + ps3)  
#

ps4 <- meta.metr1.supp %>% 
  ggplot(aes(y = New_ecos_tr, x = log(phi_C_L2R), fill = New_tr_g)) +
  stat_slab(aes(thickness = stat(pdf*n)), scale = 0.7) +
  stat_dotsinterval(side = "bottom", scale = 0.7, slab_size = NA,
                    dotsize = 3.5, alpha = 0.75) +
  ylab ("Trophic level per ecosystem type") +
  xlab ("Community spatial synchrony (Csy)") +
  theme_classic() +
  ggtitle ("D") +
  theme(legend.position="none") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18))
#

ps5 <- meta.metr1.supp %>% 
  ggplot(aes(y = New_ecos_tr, x = log(phi_S2C_L), fill = New_tr_g)) +
  stat_slab(aes(thickness = stat(pdf*n)), scale = 0.7) +
  stat_dotsinterval(side = "bottom", scale = 0.7, slab_size = NA,
                    dotsize = 1.5, alpha = 0.75) +
  ylab ("") +
  xlab ("Local population synchrony (Psy)") +
  theme_classic() +
  ggtitle ("E") +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_blank(),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18))
#
(ps4 + ps5)
#

# Site metrics
site.metr1.supp <- site.metr %>%
  mutate(New_ecos_tr = paste0(Ecosys, "_", New_tr_g)) %>%
  mutate(New_ecos_tr = case_when(
    New_ecos_tr == "Stream_Primary consumer" ~ "Lot_PriCon",
    New_ecos_tr == "Stream_Secondary consumer" ~ "Lot_SecCon",
    New_ecos_tr == "Stream_Tertiary consumer" ~ "Lot_TerCon",
    New_ecos_tr == "Stream_Producer" ~ "Lot_Prod",
    New_ecos_tr == "Lake_Producer" ~ "Len_Prod",
    New_ecos_tr == "Lake_Primary consumer" ~ "Len_PriCon")) %>% 
  mutate(New_ecos_tr = factor(New_ecos_tr,
                              levels = rev(c("Lot_TerCon", "Lot_SecCon",
                                             "Lot_PriCon", "Lot_Prod",
                                             "Len_PriCon", "Len_Prod"))))

ps6 <-  site.metr1.supp %>% 
  ggplot(aes(y = New_ecos_tr, x = log(cv_comm_site), fill = New_tr_g)) +
  stat_slab(aes(thickness = stat(pdf*n)), scale = 0.3) +
  stat_dotsinterval(side = "bottom", scale = 0.7, slab_size = NA,
                    dotsize = 5, alpha = 0.75) +
  ylab ("Trophic level per ecosystem type") +
  xlab ("Community variability (Ccv)") +
  ggtitle ("A") +
  theme_classic() +
  theme(axis.text = element_text(size = 17),
        axis.title = element_text(size = 20),
        legend.position = "none",
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18))
#

ps7 <- site.metr1.supp %>% 
  ggplot(aes(y = New_ecos_tr, x = log(mean_cv_species_site), fill = New_tr_g)) +
  stat_slab(aes(thickness = stat(pdf*n)), scale = 0.5) +
  stat_dotsinterval(side = "bottom", scale = 0.7, slab_size = NA,
                    dotsize = 5, alpha = 0.75) +
  ylab ("") +
  xlab ("Population variability (Ccv)") +
  theme_classic() +
  ggtitle ("B") +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size = 17),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_blank(),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18))
#

ps8 <- site.metr1.supp %>% 
  ggplot(aes(y = New_ecos_tr, x = log(synchrony_comm_site), fill = New_tr_g)) +
  stat_slab(aes(thickness = stat(pdf*n)), scale = 0.5) +
  stat_dotsinterval(side = "bottom", scale = 0.7, slab_size = NA,
                    dotsize = 5, alpha = 0.75) +
  ylab ("") +
  xlab ("Local population synchrony (Psy)") +
  theme_classic() +
  ggtitle ("C") +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size = 17),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_blank(),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18)) 
#
ps6+ps7+ps8
#
#==============================================================================

## Figures to explore relationships between variability 
# (and synchrony) metrics with potential confounding variables 
# (i.e., # of sites; # of time steps) 

anova(lm(log(meta.metr1.supp$CV_C_R) ~ meta.metr1.supp$Time_step *
           meta.metr1.supp$New_tr_g))

ps9 <- meta.metr1.supp %>%
  ggplot(aes(y = log(CV_C_R), x = Time_step)) +
  geom_point(aes(color = New_tr_g), size = 5, alpha = 0.75) +
  geom_line(stat="smooth",method = "lm", formula = y ~ x, size = 1.5,
            color = "darkgrey") +
  ylim (-2.5, 2.5) +
  xlab ("Number of years") +
  ylab ("Metacommunity variability (Mcv)") +
  annotate("text", x = 25, y = -2, label = "p-value (interaction) = 0.89",
           col = "black", size = 5) +
  ggtitle ("A") +
  theme_classic() +
  theme(axis.text = element_text(size = 17),
        axis.title = element_text(size = 20),
        legend.position = c(x=0.75, y=0.9),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18)) 
#
ps9
#

anova(lm(log(meta.metr1.supp$CV_C_L) ~ meta.metr1.supp$Time_step *
           meta.metr1.supp$New_tr_g))

ps10 <- meta.metr1.supp %>%
  ggplot(aes(y = log(CV_C_L), x = Time_step)) +
  geom_point(aes(color = New_tr_g), size = 5, alpha = 0.75) +
  geom_line(stat="smooth",method = "lm", formula = y ~ x, size = 1.5,
            color = "darkgrey") +
  ylim (-2.5, 2.5) +
  xlab ("Number of years") +
  ylab ("Community variability (Ccv)") +
  annotate("text", x = 25, y = -2, label = "p-value (interaction) = 0.80",
           col = "black", size = 5) +
  ggtitle ("B") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text = element_text(size = 17),
        axis.title = element_text(size = 20),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18))
#
ps10
#

anova(lm(log(meta.metr1.supp$CV_S_L) ~ meta.metr1.supp$Time_step *
             meta.metr1.supp$New_tr_g))

ps11 <- meta.metr1.supp %>%
  ggplot(aes(y = log(CV_S_L), x = Time_step)) +
  geom_point(aes(color = New_tr_g), size = 5, alpha = 0.75) +
  geom_line(stat="smooth", method = "lm", formula = y ~ x, size = 1.5,
            color = "darkgrey") +
  ylim (-2.5, 2.5) +
  xlab ("Number of years") +
  ylab ("Population variability (Pcv)") +
  annotate("text", x = 25, y = -2, label = "p-value (interaction) = 0.78",
           col = "black", size = 5) +
  ggtitle ("C") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text = element_text(size = 17),
        axis.title = element_text(size = 20),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18))
#
ps11
#
(ps9 + ps10 + ps11) 

#

anova(lm(log(meta.metr1.supp$phi_C_L2R) ~ meta.metr1.supp$Time_step *
           meta.metr1.supp$New_tr_g))

ps12 <- meta.metr1.supp %>%
  ggplot(aes(y = log(phi_C_L2R), x = Time_step,  color = New_tr_g)) +
  geom_point(size = 5, alpha = 0.75) +
  ylim (-1.5, 0.15) +
  xlab ("Number of years") +
  ylab ("Community spatial synchrony (Csy)") +
  annotate("text", x = 25, y = -1.25, label = "p-value (interaction) = 0.50",
           col = "black", size = 5) +
  ggtitle ("D") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text = element_text(size = 17),
        axis.title = element_text(size = 20),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18))
#
ps12
#

anova(lm(log(meta.metr1.supp$phi_S2C_L) ~ meta.metr1.supp$Time_step *
           meta.metr1.supp$New_tr_g))

ps13 <- meta.metr1.supp %>%
  ggplot(aes(y = log(phi_S2C_L), x = Time_step,  color = New_tr_g)) +
  geom_point(size = 5, alpha = 0.75) +
  ylim (-1.5, 0.15) +
  xlab ("Number of years") +
  ylab ("Local population synchrony (Psy)") +
  annotate("text", x = 25, y = -1.25, label = "p-value (interaction) = 0.10",
           col = "black", size = 5) +
  ggtitle ("E") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text = element_text(size = 17),
        axis.title = element_text(size = 20),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18))
#
ps13
#
ps12+ps13

###

# Now with number of sites on the x-axis

anova(lm(log(meta.metr1.supp$CV_C_R) ~ meta.metr1.supp$nSites *
           meta.metr1.supp$New_tr_g))

ps14 <- meta.metr1.supp %>%
  ggplot(aes(y = log(CV_C_R), x = nSites)) +
  geom_point(aes(color = New_tr_g), size = 5, alpha = 0.75) +
  geom_line(stat="smooth",method = "lm", formula = y ~ x, size = 1.5,
            color = "darkgrey") +
  ylim (-2.5, 2.5) +
  xlab ("Number of sites") +
  ylab ("Metacommunity variability (Mcv)") +
  annotate("text", x = 25, y = -2, label = "p-value (interaction) = 0.24",
           col = "black", size = 5) +
  ggtitle ("A") +
  theme_classic() +
  theme(axis.text = element_text(size = 17),
        axis.title = element_text(size = 20),
        legend.position = c(x=0.75, y=0.9),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18)) 
#
ps14
#

anova(lm(log(meta.metr1.supp$CV_C_L) ~ meta.metr1.supp$nSites *
           meta.metr1.supp$New_tr_g))

ps15 <- meta.metr1.supp %>%
  ggplot(aes(y = log(CV_C_L), x = nSites)) +
  geom_point(aes(color = New_tr_g), size = 5, alpha = 0.75) +
  geom_line(stat="smooth",method = "lm", formula = y ~ x, size = 1.5,
            color = "darkgrey") +
  ylim (-2.5, 2.5) +
  xlab ("Number of sites") +
  ylab ("Community variability (Ccv)") +
  annotate("text", x = 25, y = -2, label = "p-value (interaction) = 0.13",
           col = "black", size = 5) +
  ggtitle ("B") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text = element_text(size = 17),
        axis.title = element_text(size = 20),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18))
#
ps15
#

anova(lm(log(meta.metr1.supp$CV_S_L) ~ meta.metr1.supp$nSites *
           meta.metr1.supp$New_tr_g))

ps16 <- meta.metr1.supp %>%
  ggplot(aes(y = log(CV_S_L), x = nSites)) +
  geom_point(aes(color = New_tr_g), size = 5, alpha = 0.75) +
  geom_line(stat="smooth", method = "lm", formula = y ~ x, size = 1.5, 
            color = "darkgrey") +
  ylim (-2.5, 2.5) +
  xlab ("Number of sites") +
  ylab ("Population variability (Pcv)") +
  annotate("text", x = 25, y = -2, label = "p-value (interaction) = 0.33",
           col = "black", size = 5) +
  ggtitle ("C") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text = element_text(size = 17),
        axis.title = element_text(size = 20),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18))
#
ps16
#
(ps14 + ps15 + ps16) 

#

anova(lm(log(meta.metr1.supp$phi_C_L2R) ~ meta.metr1.supp$nSites *
           meta.metr1.supp$New_tr_g))

ps17 <- meta.metr1.supp %>%
  ggplot(aes(y = log(phi_C_L2R), x = nSites)) +
  geom_point(aes(color = New_tr_g), size = 5, alpha = 0.75) +
  geom_line(stat="smooth",method = "lm", formula = y ~ x, size = 1.5,
            color = "darkgrey") +
  ylim (-1.5, 0.15) +
  xlab ("Number of sites") +
  ylab ("Community spatial synchrony (Csy)") +
  annotate("text", x = 25, y = -1.25, label = "p-value (interaction) = 0.27",
           col = "black", size = 5) +
  ggtitle ("D") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text = element_text(size = 17),
        axis.title = element_text(size = 20),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18))
#
ps17
#

anova(lm(log(meta.metr1.supp$phi_S2C_L) ~ meta.metr1.supp$nSites *
           meta.metr1.supp$New_tr_g))

ps18 <- meta.metr1.supp %>%
  ggplot(aes(y = log(phi_S2C_L), x = nSites,  color = New_tr_g)) +
  geom_point(size = 5, alpha = 0.75) +
  ylim (-1.5, 0.15) +
  xlab ("Number of sites") +
  ylab ("Local population synchrony (Psy)") +
  annotate("text", x = 25, y = -1.25, label = "p-value (interaction) = 0.18",
           col = "black", size = 5) +
  ggtitle ("E") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text = element_text(size = 17),
        axis.title = element_text(size = 20),
        plot.title = element_text(vjust = -8, hjust = 0.025, size = 18))
#
ps18
#
ps17+ps18

##=========================================================================

