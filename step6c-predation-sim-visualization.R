###########################################################################
# Script Information ------------------------------------------------------
###########################################################################

### project: MCF Special Issue - Snook in Americas - Energetic Needs
### author(s): MW, WRJ, and JL
### goal: vizualizing predation simulations

# Housekeeping ------------------------------------------------------------

### load necessary packages -----
# install.packages("librarian")
librarian::shelf(tidyverse, readxl, readr, performance, ggpubr, 
                 mgcv, glmmTMB, MuMIn, corrplot, ggeffects)

# read in necessary data and wrangle/summarize ----------------------------
ind <- read_csv('data/prey-numeric-needs-2025.csv') |> 
  mutate(diet = factor(diet,
                       levels = c("0% fish/100% invert", "25% fish/75% invert",
                                  "50% fish/50% invert", "75% fish/25% invert",
                                  "100% fish/0% invert"))) |> 
  mutate(weight_factor = as.character(weight_kg),
         weight_factor = case_when(
           weight_factor == '0.5' ~ '0.5 kg',
           weight_factor == '3.3' ~ '3.3 kg',
           weight_factor == '6.2' ~ '6.2 kg'
         ))

bm <- read_csv('data/prey-bm-needs-2025.csv') |> 
  mutate(diet = factor(diet,
                       levels = c("0% fish/100% invert", "25% fish/75% invert",
                                  "50% fish/50% invert", "75% fish/25% invert",
                                  "100% fish/0% invert"))) |> 
  mutate(weight_factor = as.character(weight_kg),
         weight_factor = case_when(
           weight_factor == '0.5' ~ '0.5 kg',
           weight_factor == '3.3' ~ '3.3 kg',
           weight_factor == '6.2' ~ '6.2 kg'
         ))

### prey biomass needed ---
a <- ggplot(bm, aes(x = as.factor(temp), y = mean_bm_needs, color = diet)) +
  geom_point(size = 4, position = position_dodge(width = 1)) +                             
  scale_color_viridis_d(option = 'turbo', end = 1) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),  
                width = 0,
                linewidth = 1.5,
                position = position_dodge(width = 1)) +   
  facet_wrap(~weight_factor) + 
  labs(x = "Water Temperature (°C)", y = "Prey Demand (g/day)", color = 'Diet Type', fill = 'Diet Type') +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        # legend.position = "right",
        # legend.text = element_text(size = 12, color = "black", face = 'bold'),
        # legend.title = element_text(size = 12, color = "black", face = 'bold'),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5))

### individual prey items needed ---
b <- ggplot(ind, aes(x = as.factor(temp), y = mean_num_needs, color = diet)) +
  geom_point(size = 4, position = position_dodge(width = 1)) +                             
  scale_color_viridis_d(option = 'turbo', end = 1) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),  
                width = 0,
                linewidth = 1.5,
                position = position_dodge(width = 1)) +   
  facet_wrap(~weight_factor) + 
  labs(x = "Water Temperature (°C)", y = "Prey Demand (#/day)", color = 'Diet Type', fill = 'Diet Type') +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        # legend.position = "right",
        # legend.text = element_text(size = 12, color = "black", face = 'bold'),
        # legend.title = element_text(size = 12, color = "black", face = 'bold'),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5))

arranged <- ggarrange(a, b,
          ncol = 1, align = "h")

arranged |> 
  annotate_figure(
    bottom = text_grob(label = "Water Temperature (°C)",
                       just = 'centre', 
                       color = "black", 
                       face = "bold",
                       size = 15))

# ggsave('plots/2025/fig4.png',
#        dpi = 600, units= 'in', height = 10, width = 12)

### legend ---

ggplot(ind, aes(x = as.factor(temp), y = mean_num_needs, color = diet)) +
  geom_point(size = 4, position = position_dodge(width = 1)) +                             
  scale_color_viridis_d(option = 'turbo', end = 1) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),  
                width = 0,
                linewidth = 1.5,
                position = position_dodge(width = 1)) +   
  facet_wrap(~weight_factor) + 
  labs(x = "Water Temperature (°C)", y = "Prey Demand (#/day)", color = 'Diet Type', fill = 'Diet Type') +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # legend.position = "none",
        legend.position = "right",
        legend.text = element_text(size = 12, color = "black", face = 'bold'),
        legend.title = element_text(size = 12, color = "black", face = 'bold'),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5))

# ggsave('plots/2025/fig4-legend.png',
#        dpi = 600, units= 'in', height = 10, width = 12)

