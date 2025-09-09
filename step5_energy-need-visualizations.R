###########################################################################
# Script Information ------------------------------------------------------
###########################################################################

### project: MCF Special Issue - Snook in Americas - Energetic Needs
### author(s): MW, WRJ, and JL
### goal: visualizing snook caloric need estimates (fig 3)

# Housekeeping ------------------------------------------------------------

### load necessary packages -----
# install.packages("librarian")
librarian::shelf(tidyverse, readxl, readr, akima, plotly, ggpubr, cowplot, patchwork)

### read in data and set sequence of temp and trophic data ---
# for summary information ----
pop_df <- read_csv("data/manuscript/simulated-snook-population-short.csv")
# for figure 2a generation ----
dt_test <- read_csv("data/manuscript/metabolic-rates-temp15thru14by1.csv")
# for figure 2b generation ----
met_df <- read_csv("data/manuscript/metabolic-rates-temp15thru14by5.csv") |> filter(act == 2.5, troph == 4.5)

# simulated population summary  -------------------------------------------
test <- pop_df |> 
  mutate(sl_mm = -14.7684 + (0.9338*fl_mm_vbf),
         sl_cm = sl_mm/10)

summ <- test |> 
  summarize(length_m = mean(sl_cm),
            length_sd = sd(sl_cm),
            weight_m = mean(weight_kg),
            weight_sd = sd(weight_kg),
            age_m = mean(age_in_years),
            age_sd = sd(age_in_years))


# figure 2a ---------------------------------------------------------------
met_constrained <- dt_test |> 
  filter(act == 2.5 & troph == 4.5) |> 
  mutate(rel_cog = Cost_growth_j_g/Total_metabolic_rate_j_d,
         sl_mm = -14.7684 + (0.9338*fl_mm),
         sl_cm = sl_mm/10,
         Total_metabolic_rate_kj_d = Total_metabolic_rate_j_d/1000,
         Total_kilocalories_d = Total_metabolic_rate_kj_d/4.184)
summary(met_constrained)

interp_kcal <- with(met_constrained, interp(
  x = weight_g / 1000,
  y = temp,
  z = Total_kilocalories_d,
  duplicate = "mean"
))

# Convert to long format
interp_df <- interp2xyz(interp_kcal, data.frame = TRUE)

# Filter out edges where NA might linger (buffer of 0.1 around extremes)
interp_df <- interp_df |>
  filter(!is.na(z),
         x > min(x) + 0.1, x < max(x) - 0.1,
         y > min(y) + 0.5, y < max(y) - 0.5)

a <- interp_df |> 
  ggplot(aes(x = x, y = y, fill = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(name = expression(bold("Total kcal"~"·"~day^-1))) +
  labs(
    x = "Weight (kg)", 
    y = expression(bold("Water Temperature ("*degree*C*")"))
  ) +
  theme(axis.text = element_text(size = 12, face = "bold", colour = "black"), 
        axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 14, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        legend.text = element_text(size = 9, face = "bold", colour = "black"))
a

# figure 2b ---------------------------------------------------------------
met_df_constrained <- met_df |> 
  filter(act == 2.5 & troph == 4.5) |> 
  mutate(rel_cog = Cost_growth_j_g/Total_metabolic_rate_j_d,
         sl_mm = -14.7684 + (0.9338*fl_mm),
         sl_cm = sl_mm/10)

# Smoothing the data
smooth_met_df <- met_df_constrained |> 
  mutate(mmr_kj_d = Total_metabolic_rate_j_d/1000,
         smr_kj_d = Resting_metabolic_rate_j_d/1000,
         mmr_gO2_d = mmr_kj_d/14.07,
         smr_gO2_d = smr_kj_d/14.07,
         mmr_mgO2_hr = (mmr_gO2_d*1000)/24,
         smr_mgO2_hr = (smr_gO2_d*1000)/24,
         mmr_cal_d = Total_metabolic_rate_j_d/4.184,
         smr_cal_d = Resting_metabolic_rate_j_d/4.184)

summ_met_df <- smooth_met_df |> 
  summarize(mmr_cal_d_m = mean(mmr_cal_d),
            mmr_cal_d_sd = sd(mmr_cal_d))

summary <- smooth_met_df |> 
  mutate(weight_kg = weight_g/1000,
         weight_kg = round(weight_kg, digits = 3),
         mmr_cal_d = mmr_cal_d/1000) |> 
  filter(weight_kg %in% c(0.509, 3.282, 6.174)) |> 
  group_by(weight_kg) |> 
  summarize(kcal_d_m = mean(mmr_cal_d),
            kcal_d_sd = sd(mmr_cal_d))

summary1 <- smooth_met_df |>   
  mutate(weight_kg = weight_g/1000,
         weight_kg = round(weight_kg, digits = 3),
         mmr_cal_d = mmr_cal_d/1000) |> 
  summarize(kcal_d_m = mean(mmr_cal_d),
            kcal_d_sd = sd(mmr_cal_d))

summary2 <- smooth_met_df |> 
  mutate(weight_kg = weight_g/1000,
         weight_kg = round(weight_kg, digits = 3),
         mmr_cal_d = mmr_cal_d/1000) |> 
  filter(weight_kg %in% c(0.509, 3.282, 6.174),
         temp %in% c(15,25,35)) |> 
  group_by(weight_kg, temp) |> 
  summarize(kcal_d_m = mean(mmr_cal_d),
            kcal_d_sd = sd(mmr_cal_d))

summary3 <- summary2 |> 
  group_by(weight_kg) |> 
  mutate(inc = max(kcal_d_m)/min(kcal_d_m))

summary4 <- smooth_met_df |> 
  mutate(weight_kg = weight_g/1000,
         weight_kg = round(weight_kg, digits = 1),
         mmr_cal_d = mmr_cal_d/1000) |> 
  filter(weight_kg %in% c(0.5, 3.3, 6.2)) 

b <- summary4 |> 
  ggplot(aes(x=temp, y=mmr_cal_d, group = as.factor(weight_kg), color = as.factor(weight_kg))) +
  geom_smooth(se = FALSE, size = 2) +
  scale_color_viridis_d(option = 'turbo', end = 1) +
  labs(y = expression(bold("Total kcals"~"·"~day^-1)),
       color = 'Weight (kg)',
       x = expression(bold("Water Temperature ("*degree*C*")"))) +
  theme(axis.text = element_text(size = 12, face = "bold", colour = "black"), 
        axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 14, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        legend.text = element_text(size = 9, face = "bold", colour = "black"))
b

### patch it together ----
a <- a  + theme(strip.text = element_blank(), plot.tag.position = c(0.03, 1.03), plot.tag = element_text(size = 14, face = "bold"))
b <- b  + theme(strip.text = element_blank(), plot.tag.position = c(0.03, 1.03), plot.tag = element_text(size = 14, face = "bold"))

final <- (a / b) + 
  plot_layout(heights = c(1, 0.75)) +
  plot_annotation(tag_levels = 'a')
final

# ggsave('plots/2025/fig2-twopanel-kcalneeds.png',
#        dpi = 600, units= 'in', height = 6, width = 6)
