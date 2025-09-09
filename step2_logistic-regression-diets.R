###########################################################################
# Script Information ------------------------------------------------------
###########################################################################

### project: MCF Special Issue - Snook in Americas - Energetic Needs
### author(s): MW & WRJ
### goal: conduct log regression for diet analyses (fig 2)

# Housekeeping ------------------------------------------------------------

### load necessary packages -----
# install.packages("librarian")
librarian::shelf(tidyverse, readxl, readr, performance, ggpubr, 
                 mgcv, glmmTMB, MuMIn, corrplot, ggeffects, dataRetrieval,
                 slider, grid, patchwork, writexl)

palette <- c("Fish"="#44AA99",
             "Invertebrate"="#D55E00",
             "Empty"="#8b6b93")

### read in datasets -----
stage <- read_csv("../MAP/data/hydrology/ssr_waterlevels.csv") |> 
  mutate(date = mdy(date),
         year = year(date),
         month = month(date),
         usace_wyear = if_else(month >=5, year+1, year))

diet <- read_csv("data/diet-data/primary_snookdiet_06062024.csv") |> 
  select(-date) |> 
  mutate(date = date_good,
         year = year(date),
         month = month(date),
         usace_wyear = if_else(month >=5, year+1, year)) |> 
  select(-date_good) |> 
  rename(id = ID) |> 
  mutate(fish_diet_binary = if_else(coarse_3 == "vert", 1, 0),
         invert_diet_binary = if_else(coarse_3 == "invert", 1, 0),
         empty_diet_binary = if_else(coarse_3 == "empty", 1, 0)) |> 
  select(usace_wyear, date, site, id, sl, fish_diet_binary, invert_diet_binary, empty_diet_binary)
glimpse(diet) 

diet_dates <- diet |> select(date) |> distinct() |> mutate(sample = 'diet')
glimpse(diet_dates)

### Bottle Creek Temp near Rookery Branch, Everglades NPS 
siteNumber <- "022908295" # found here https://waterdata.usgs.gov/monitoring-location/022908295/#parameterCode=00065&period=P7D
parameterCd <- "00010" # parameter codes listed here -> http://water.nv.gov/hearings/past/Spring%20Valley%202006/exhibits/SNWA/5__Hydrochemistry_and_Geochemistry/Data/USGS/USGS_NWIS/ParameterCodes.htm

bc_temperature <- readNWISdv(
  siteNumber, parameterCd,
  "2000-01-01", "2024-12-31"
)

botcreektemp <- bc_temperature |> 
  mutate(temp_c = X_00010_00003) |> 
  mutate(date = Date) |> 
  select(date, temp_c) |> 
  arrange(date) |> 
  complete(date = seq.Date(min(date), max(date), by = "day")) |> 
  arrange(date) |> 
  mutate(temp_c = zoo::na.approx(temp_c, x = date, na.rm = FALSE)) |> 
  filter(date >= "2010-12-01")
summary(botcreektemp)

env <- botcreektemp |> 
  left_join(stage) |> 
  left_join(diet_dates)
summary(env)
glimpse(env)

df <- env |>
  arrange(env) |> 
  mutate(
    stage3d = slide_dbl(stage_cm, mean, .before = 3, .complete = TRUE),
    stage7d = slide_dbl(stage_cm, mean, .before = 7, .complete = TRUE),
    stage14d = slide_dbl(stage_cm, mean, .before = 14, .complete = TRUE),
    stage30d = slide_dbl(stage_cm, mean, .before = 30, .complete = TRUE),
    temp3d = slide_dbl(temp_c, mean, .before = 3, .complete = TRUE),
    temp7d = slide_dbl(temp_c, mean, .before = 7, .complete = TRUE),
    temp14d = slide_dbl(stage_cm, mean, .before = 14, .complete = TRUE),
    temp30d = slide_dbl(temp_c, mean, .before = 30, .complete = TRUE)
  ) |> 
  group_by(year, month) |> 
  mutate(stagem = mean(stage_cm, na.rm = TRUE),
         tempm = mean(temp_c, na.rm = TRUE)) |> 
  filter(sample == 'diet')
print(df)
glimpse(df)
summary(df)

mod_dta <- diet |> 
  left_join(df, by = c("usace_wyear", "date")) |> 
  rename(stage0d = stage_cm,
         temp0d = temp_c) |> 
  select(-c('stage_feet', 'year', 'month', 'sample'))
print(mod_dta)
glimpse(mod_dta)
summary(mod_dta)

# lag env factors 0 days and replot -----------------------------------

mod_dt <- mod_dta |> 
  mutate(stage_cm = stage0d,
         temp_c = temp0d) |> 
  select(-temp0d:-tempm)
glimpse(mod_dt)

global_fish_poly <- glm(
  fish_diet_binary ~ 
    poly(stage_cm, 2) + 
    stage_cm + 
    poly(temp_c, 2) +
    temp_c +
    sl + 
    usace_wyear,
  data = mod_dt,
  na.action = "na.fail",
  family = binomial(link = "logit")
)

model_set_poly_fish <- dredge(global_fish_poly,
                              subset = !(`temp_c` && `poly(temp_c, 2)`) & 
                                !(`stage_cm` && `poly(stage_cm, 2)`),
                              eval = TRUE) |> filter(delta < 4)

global_invert_poly <- glm(
  invert_diet_binary ~ 
    poly(stage_cm, 2) + 
    stage_cm + 
    poly(temp_c, 2) +
    temp_c +
    sl + 
    usace_wyear,
  data = mod_dt,
  na.action = "na.fail",
  family = binomial(link = "logit"),
)

model_set_poly_invert <- dredge(global_invert_poly,
                                subset = !(`temp_c` && `poly(temp_c, 2)`) & 
                                  !(`stage_cm` && `poly(stage_cm, 2)`),
                                eval = TRUE
) |> filter(delta < 4)

global_empty_poly <- glm(
  empty_diet_binary ~ 
    poly(stage_cm, 2) + 
    stage_cm + 
    poly(temp_c, 2) +
    temp_c +
    sl + 
    usace_wyear,
  data = mod_dt,
  na.action = "na.fail",
  family = binomial(link = "logit"),
)

model_set_poly_empty <- dredge(global_empty_poly,
                               subset = !(`temp_c` && `poly(temp_c, 2)`) & 
                                 !(`stage_cm` && `poly(stage_cm, 2)`),
                               eval = TRUE
) |> filter(delta < 4)

### fit selected models ----
m1_poly_fish <- glm(
  fish_diet_binary ~ stage_cm + temp_c,
  data = mod_dt,
  na.action = "na.fail",
  family = binomial(link = "logit"),
)

performance(m1_poly_fish) # PCP = 0.62, Tjur's R2 = 0.163
summary(m1_poly_fish) #all significant
#                     Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       -1.787394   0.589802  -3.030  0.00244 ** 
#   stage_cm          -0.018152   0.008798  -2.063  0.03910 *  
#   poly(temp_c, 2)1 -12.325786   2.782164  -4.430 9.41e-06 ***
#   poly(temp_c, 2)2  -6.647713   2.558556  -2.598  0.00937 ** 
#   sl                 0.032880   0.014504   2.267  0.02339 *  

m1_poly_invert <- glm(
  invert_diet_binary ~ stage_cm + sl,
  data = mod_dt,
  na.action = "na.fail",
  family = binomial(link = "logit"),
)

performance(m1_poly_invert) # PCP = 0.62, Tjur's R2 = 0.10
summary(m1_poly_invert) #all significant
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.921472   0.642348   1.435 0.151419    
# stage_cm     0.027129   0.007004   3.874 0.000107 ***
#   sl          -0.058186   0.017302  -3.363 0.000771 ***

m1_poly_empty <- glm(
  empty_diet_binary ~ temp_c,
  data = mod_dt,
  na.action = "na.fail",
  family = binomial(link = "logit"),
)

performance(m1_poly_empty) # PCP = 0.589, Tjur's R2 = 0.09
summary(m1_poly_empty) #all significant
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -7.03182    1.37044  -5.131 2.88e-07 ***
#   temp_c       0.25875    0.05462   4.738 2.16e-06 ***

# marginal effects plots --------------------------------------------------

### identify bounds for visualizations ----
summary(mod_dt)

### fish in diet ---
stage_effect_fish <- ggemmeans(m1_poly_fish, 'stage_cm[-15:58 by=0.1]')
plot(stage_effect_fish)
temp_effect_fish <- ggemmeans(m1_poly_fish, 'temp_c[15:32 by=0.1]')
plot(temp_effect_fish)

### inverts in diet ---
stage_effect_invert <- ggemmeans(m1_poly_invert, 'stage_cm[-15:58 by=0.1]')
plot(stage_effect_invert)
sl_effect_invert <- ggemmeans(m1_poly_invert, 'sl[23:75 by=0.1]')
plot(sl_effect_invert)

### empty diet
temp_effect_empty <- ggemmeans(m1_poly_empty, 'temp_c[15:32 by=0.1]')
plot(temp_effect_empty)

### revised diet figure -----

### temp effects ---
temp_effect_fish$diet <- "Fish"
temp_effect_empty$diet <- "Empty"
temp <- rbind(temp_effect_fish, temp_effect_empty)

a <- temp |> 
  ggplot(aes(x = x, y = predicted, color = diet)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = diet), alpha = 0.5)+
  geom_line(size = 1.5) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by =0.25))+
  scale_x_continuous(limits = c(15,32), breaks = seq(16,32, by = 4))+
  labs(x = expression(bold("Water Temperature ("*degree*C*")")), y = 'Probability of Occurence') +
  scale_color_manual(values = palette)+
  scale_fill_manual(values = palette) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title.x = element_text(size = 15, face = "bold", colour = "black"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
a

### stage effects ---
stage_effect_fish$diet <- "Fish"
stage_effect_invert$diet <- "Invertebrate"
stage <- rbind(stage_effect_fish, stage_effect_invert)

b <- stage |> 
  ggplot(aes(x = x, y = predicted, color = diet)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = diet), alpha = 0.5)+
  geom_line(size = 1.5) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by =0.25))+
  scale_x_continuous(limits = c(-15,60), breaks = seq(-15,60, by = 10))+
  labs(x = 'Marsh Stage (cm)', y = 'Probability of Occurence') +
  scale_color_manual(values = palette)+
  scale_fill_manual(values = palette) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title.x = element_text(size = 15, face = "bold", colour = "black"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
b

### size effects ---
sl_effect_invert$diet <- "Invertebrate"
size <- sl_effect_invert

c <- size |> 
  ggplot(aes(x = x, y = predicted, color = diet)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = diet), alpha = 0.5)+
  geom_line(size = 1.5) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by =0.25)) +
  scale_x_continuous(limits = c(23,75), breaks = seq(25,75, by = 10)) +
  labs(x = 'Standard Length (cm)', y = 'Probability of Occurence') +
  scale_color_manual(values = palette)+
  scale_fill_manual(values = palette) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title.x = element_text(size = 15, face = "bold", colour = "black"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
c

rm(list = setdiff(ls(), c('mod_dt', 'palette', 'a', 'b', 'c')))

# visualize predictions ---------------------------------------------------

### create dummy legend ----
dummy_data <- data.frame(
  x = c(1, 2, 3),
  y = c(0.1, 0.5, 0.9),
  diet = factor(c("Fish", "Invertebrate", "Empty"), levels = c("Fish", "Invertebrate", "Empty"))
)

legend_plot <- ggplot(dummy_data, aes(x = x, y = y, color = diet, fill = diet)) +
  geom_line() +
  geom_ribbon(aes(ymin = y - 0.05, ymax = y + 0.05), alpha = 0.5) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(fill = 'Diet', color = 'Diet') +
  theme_minimal() +
  theme(legend.position = "right",
        legend.text = element_text(size = 14, face = "bold", colour = "black"), 
        legend.title = element_text(size = 15, face = "bold", colour = "black"))
legend_shared <- cowplot::get_legend(legend_plot)

### clean up individual figures and prepare for patch ----
a <- a + theme(plot.tag.position = c(0.02, 1.02), plot.tag = element_text(size = 14, face = "bold"))
b <- b + theme(plot.tag.position = c(0.00, 1.02), plot.tag = element_text(size = 14, face = "bold"), axis.text.y = element_blank())
c <- c + theme(plot.tag.position = c(0.00, 1.02), plot.tag = element_text(size = 14, face = "bold"), axis.text.y = element_blank())

final_plot <- (a + b + c) +
  plot_layout(widths = c(1, 1, 1)) +
  plot_annotation(tag_levels = 'a') &
  ggplot2::theme(
    legend.position = 'none',
    plot.margin = ggplot2::margin(t = 10, r = 10, b = 5, l = 5)
  )
final_plot

### patch it together ----

grid::grid.newpage()
grid::grid.draw(
  patchwork::wrap_elements(full = final_plot) +
    
    # Add shared y-axis label
    patchwork::inset_element(
      grid::textGrob('Diet Probability', rot = 90,
                     gp = grid::gpar(fontsize = 15, fontface = 'bold')),
      left = -0.01, bottom = 0.25, right = 0.01, top = 0.75, align_to = 'full'
    ) +
    
    # Add shared legend
    patchwork::inset_element(
      legend_shared,
      left = 0.85, bottom = 0.83, right = 0.95, top = 0.95, align_to = 'full'
    )
)

# ggsave('plots/2025/fig2-revised062025.png',
#        dpi = 600, units= 'in', height = 6, width = 12)

summary <- mod_dt |> 
  group_by(usace_wyear) |> 
  summarize(
    n = n_distinct(id),
    meanlength = mean(sl),
    sdlength = sd(sl),
    meantemp = mean(temp_c),
    sdtemp = sd(temp_c),
    meanstage = mean(stage_cm),
    sdstage = sd(stage_cm),
    .groups = 'drop'
  ) |> 
  mutate(
    across(c(meanlength, sdlength, meantemp, sdtemp, meanstage, sdstage), ~round(.x, 1)),
    sdlength = ifelse(is.na(sdlength), 0, sdlength),
    sdtemp = ifelse(is.na(sdtemp), 0, sdtemp),
    sdstage = ifelse(is.na(sdstage), 0, sdstage),
    length_pm = paste0(meanlength, " ± ", sdlength),
    temp_pm = paste0(meantemp, " ± ", sdtemp),
    stage_pm = paste0(meanstage, " ± ", sdstage)
  ) |> 
  select(usace_wyear, n, length_pm, temp_pm, stage_pm)

view(summary)

# write_xlsx(summary, "data/daily-env-summary.xlsx")

summary1 <- mod_dt |> 
  summarize(
    n = n_distinct(id),
    meanlength = mean(sl),
    sdlength = sd(sl),
    meantemp = mean(temp_c),
    sdtemp = sd(temp_c),
    meanstage = mean(stage_cm),
    sdstage = sd(stage_cm),
    .groups = 'drop'
  ) |> 
  mutate(
    across(c(meanlength, sdlength, meantemp, sdtemp, meanstage, sdstage), ~round(.x, 1)),
    sdlength = ifelse(is.na(sdlength), 0, sdlength),
    sdtemp = ifelse(is.na(sdtemp), 0, sdtemp),
    sdstage = ifelse(is.na(sdstage), 0, sdstage),
    length_pm = paste0(meanlength, " ± ", sdlength),
    temp_pm = paste0(meantemp, " ± ", sdtemp),
    stage_pm = paste0(meanstage, " ± ", sdstage)
  ) |> 
  select(n, length_pm, temp_pm, stage_pm)

view(summary1)

### model selection tables -----

### fish ---
global_fish_poly <- glm(
  fish_diet_binary ~ 
    poly(stage_cm, 2) + 
    stage_cm + 
    poly(temp_c, 2) +
    temp_c +
    sl + 
    usace_wyear,
  data = mod_dt,
  na.action = "na.fail",
  family = binomial(link = "logit")
)

model_set_poly_fish <- dredge(global_fish_poly,
                              subset = !(`temp_c` && `poly(temp_c, 2)`) & 
                                !(`stage_cm` && `poly(stage_cm, 2)`),
                              eval = TRUE)

model_set_poly_fish = model_set_poly_fish |> 
  as_tibble() |> mutate(response = 'fish')

### inverts ---
global_invert_poly <- glm(
  invert_diet_binary ~ 
    poly(stage_cm, 2) + 
    stage_cm + 
    poly(temp_c, 2) +
    temp_c +
    sl + 
    usace_wyear,
  data = mod_dt,
  na.action = "na.fail",
  family = binomial(link = "logit"),
)

model_set_poly_invert <- dredge(global_invert_poly,
                                subset = !(`temp_c` && `poly(temp_c, 2)`) & 
                                  !(`stage_cm` && `poly(stage_cm, 2)`),
                                eval = TRUE
)

model_set_poly_invert = model_set_poly_invert |> 
  as_tibble() |> mutate(response = 'invertebrate')

### empty ---
global_empty_poly <- glm(
  empty_diet_binary ~ 
    poly(stage_cm, 2) + 
    stage_cm + 
    poly(temp_c, 2) +
    temp_c +
    sl + 
    usace_wyear,
  data = mod_dt,
  na.action = "na.fail",
  family = binomial(link = "logit"),
)

model_set_poly_empty <- dredge(global_empty_poly,
                               subset = !(`temp_c` && `poly(temp_c, 2)`) & 
                                 !(`stage_cm` && `poly(stage_cm, 2)`),
                               eval = TRUE
)

model_set_poly_empty = model_set_poly_empty |> 
  as_tibble() |> mutate(response = 'empty')

all <- rbind(model_set_poly_fish, model_set_poly_invert, model_set_poly_empty)
glimpse(all)
# writexl::write_xlsx(all, 'tables/full-model-selection-all-response.xlsx')

### model summaries ---
m1_poly_fish <- glm(
  fish_diet_binary ~ stage_cm + temp_c,
  data = mod_dt,
  na.action = "na.fail",
  family = binomial(link = "logit"),
)

m1_poly_invert <- glm(
  invert_diet_binary ~ stage_cm + sl,
  data = mod_dt,
  na.action = "na.fail",
  family = binomial(link = "logit"),
)

m1_poly_empty <- glm(
  empty_diet_binary ~ temp_c,
  data = mod_dt,
  na.action = "na.fail",
  family = binomial(link = "logit"),
)

performance(m1_poly_fish)
summary(m1_poly_fish) 

performance(m1_poly_invert) 
summary(m1_poly_invert)   

performance(m1_poly_empty)
summary(m1_poly_empty)
