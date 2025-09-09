###########################################################################
# Script Information ------------------------------------------------------
###########################################################################

### project: MCF Special Issue - Snook in Americas - Energetic Needs
### author(s): MW, WRJ, and JL
### goal: generate estimates of snook caloric needs using fishflux + senitivity analyses

# Housekeeping ------------------------------------------------------------

### load necessary packages -----
# install.packages("librarian")
# devtools::install_github("nschiett/fishflux", dependencies=TRUE)
# library(fishflux)
librarian::shelf(tidyverse, readxl, readr, fishflux, plotly, sensitivity, patchwork, cowplot)

# set.seed(50)

### read in data and set sequence of temp and trophic data ---
sim <- read_csv("data/manuscript/simulated-snook-population-short.csv")
temp = seq(from = 12, to = 33, length.out = 70)
troph = seq(from = 2.5, to = 5, by = 0.5)
act = seq(from = 1, to = 4, by = 0.5)
a = c(0.5964, 0.6732, 0.75)

### skip to line 82 of script unless metabolic rate needs recalculated

perciforms <- data.frame(
  family = c("Apogonidae", "Centrarchidae", "Cichlidae", "Gobiidae",
             "Pomacentridae", "Centropomidae")
)

ecolog <- crossing(temp, troph, act, a)

all_step1 <- perciforms |> cross_join(ecolog) |>  as_tibble()

met_coeff <- all_step1 |>
  mutate(i = row_number()) |>
  group_by(i) |>
  nest() |>
  mutate(met_coeffs = map(data, \(data)
                              metabolism(family = data$family, temp = data$temp, troph_m = data$troph)
  ))

met_coeff2 <- met_coeff |>
  unnest(data) |>
  unnest_wider(met_coeffs)

met_coeff3 <- met_coeff2 |>
  group_by(temp, troph , act) |>
  summarize(f0 = mean(f0_m),
            a = mean(alpha_m),
            B0 = mean(b0_m))

met_coeff4 <- met_coeff3 |>
  select(-a, -B0)

sim_slim <- sim |>
  select(fl_mm_rounded, weight_g, daily_growth_g) |>
  rename(fl_mm = fl_mm_rounded)

all_step2 <- sim_slim |> cross_join(met_coeff4) |>  as_tibble()
all_step3 <- all_step2 |> cross_join(as.data.frame(a)) |> as_tibble()

met_rates <- all_step3 |>
  mutate(i = row_number()) |>
  group_by(i) |>
  nest() |>
  mutate(metabolic_rate = map(data, \(data)
                              ### only use f0 for B0 if accounted for temp and taxa as done above
                              metabolic_rate(temp = data$temp, troph = data$troph, asp = 1.90693, B0 = data$f0,
                                             m_max = 7868.990, m = data$weight_g, growth_g_day = data$daily_growth_g,
                                             a = data$a, f = data$act)
  ))

met_rates2 <- met_rates |>
  unnest(data) |>
  unnest_wider(metabolic_rate)

# read out data -----------------------------------------------
# write_csv(met_rates2, "data/manuscript/final-MR-model-dataset-june2025.csv")

# read data back in -------------------------------------------------------
met_rates2 <- read_csv('data/manuscript/final-MR-model-dataset-june2025.csv')

df <- met_rates2
rm(list = setdiff(ls(), c('df', 'sim')))
glimpse(df)

# sensitivity analyses ----------------------------------------------------

snook_kcal <- function(X) {
  apply(X, 1, function(x) {
    names(x) <- c("temp", "troph", "m", "growth_g_day", "f0", "a", "f")
    result <- metabolic_rate(
      temp = x["temp"],
      troph = x["troph"],
      asp = 1.91,
      B0 = x["f0"],
      m_max = 7869,
      m = x["m"],
      growth_g_day = x["growth_g_day"],
      a = x["a"],
      f = x["f"]
    )
    result$Total_metabolic_rate_j_d / 4.184 / 1000  # Convert to kcal/day
  })
}

factors <- c("temp", "troph", "m", "growth_g_day", "f0", "a", "f")
lower_bounds <- c(15, 1.5, 500, 0.1, 0.00069, 0.5, 1.5)
upper_bounds <- c(35, 5.0, 6600, 2, 0.007, 0.75, 4)

### sobol global sensitivity approach ----
set.seed(20)
sobol_design <- sensitivity::sobol2002(
  model = NULL,
  X1 = as.data.frame(lhs::randomLHS(1000, length(factors))),
  X2 = as.data.frame(lhs::randomLHS(1000, length(factors))),
  order = 1
)

for (i in seq_along(factors)) {
  sobol_design$X1[[i]] <- lower_bounds[i] + (upper_bounds[i] - lower_bounds[i]) * sobol_design$X1[[i]]
  sobol_design$X2[[i]] <- lower_bounds[i] + (upper_bounds[i] - lower_bounds[i]) * sobol_design$X2[[i]]
}

sobol_result <- sensitivity::sobol2002(
  model = snook_kcal,
  X1 = sobol_design$X1,
  X2 = sobol_design$X2,
  # nboot = 100
)

sobol_df <- data.frame(
  parameter = factors,
  first_order = sobol_result$S,
  total_order = sobol_result$T
) |> 
  rename(first_order = original,
         total_order = original.1)

sobol_df_long <- sobol_df |>
  tidyr::pivot_longer(cols = c(first_order, total_order),
                      names_to = "order", values_to = "sensitivity") |> 
  mutate(parameter_long = case_when(
    parameter == 'a' ~ 'Mass-Scaling Exponent (MSE)',
    parameter == 'growth_g_day' ~ 'Cost of Daily Growth (G)',
    parameter == 'f0' ~ 'Metabolic Normalization Constant (MNC)',
    parameter == 'f' ~ 'Activity Scope (AS)',
    parameter == 'troph' ~ 'Trophic Level (TL)',
    parameter == 'm' ~ 'Wet Mass (WM)',
    parameter == 'temp' ~ 'Environmental Temperature (T)')
  ) |> 
  mutate(parameter_short = case_when(
    parameter == 'a' ~ 'MSE',
    parameter == 'growth_g_day' ~ 'G',
    parameter == 'f0' ~ 'MNC',
    parameter == 'f' ~ 'AS',
    parameter == 'troph' ~ 'TL',
    parameter == 'm' ~ 'WM',
    parameter == 'temp' ~ 'T')
  ) |> 
  mutate(parameter_short = factor(parameter_short, 
                                  levels = c('MSE', 'WM', 'MNC', 'AS',
                                             'T', 'TL', 'G'))) |> 
  mutate(indices = case_when(
    order == 'first_order' ~ 'First-Order',
    order == 'total_order' ~ 'Total-Order'
  )) |> 
  mutate(indices = factor(indices, levels = c('Total-Order', 'First-Order')))

m <- lm(total_order ~ first_order, data = sobol_df)
summary(m)

total <- sobol_df_long |> filter(order == 'total_order')

a <- sobol_df_long |>
  filter(parameter_short %in% c('MSE', 'MNC', 'AS', 'WM')) |> 
  ggplot(aes(x = parameter_short, y = sensitivity, fill = indices)) +
  geom_col(position = "dodge") +
  labs(title = "Sobol Global Sensitivity Analysis",
       x = 'Model Parameter',
       y = "Propotion of Variance Explained",
       fill = 'Indices') +
  theme_classic() +
  scale_y_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title.x = element_text(size = 15, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 15, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.85,0.9),
        legend.title = element_text(size = 15, face = "bold", colour = "black"),
        legend.text = element_text(size = 14, face = "bold", colour = "black"))
a
# ggsave('plots/2025/supplemental-sensitivity-analysis.png',
#        dpi = 600, units= 'in', height = 6, width = 6)

### morris approach ----
set.seed(20)
morris_result <- sensitivity::morris(
  model = snook_kcal,
  factors = factors,
  r = 1000,
  design = list(type = "oat", levels = 4, grid.jump = 2),
  binf = lower_bounds,
  bsup = upper_bounds,
  scale = TRUE
)

summary_df <- data.frame(
  parameter = factors,
  mu = apply(morris_result$ee, 2, mean),
  sigma = apply(morris_result$ee, 2, sd)
)

### creating the second part of the plot ----
glimpse(df)
df1 <- df |> mutate(kcal_day = Total_metabolic_rate_j_d / 4.184 / 1000)
glimpse(df1)
summary(df1)
cb_palette <- c("#E69F00", "#56B4E9", "#009E73")

b <- df1 |> 
  mutate(temp = case_when(
    temp >= 22.5 & temp <= 22.7 ~ 22.5,
    TRUE ~ temp
  )) |> 
  filter(temp %in% c(12, 22.5, 33)) |>
  mutate(
    temp_label = factor(
      temp,
      levels = c(12, 22.5, 33),
      labels = c("12°C", "22.5°C", "33°C")
    )
  ) |> 
  filter(a %in% c(0.5964, 0.6732, 0.7500)) |> 
  filter(act == 2.5) |> 
  filter(troph == 4.5) |> 
  mutate(temp = as.factor(temp),
         a = as.factor(a)) |> 
  ggplot(aes(x = weight_g/1000, y = kcal_day, group = a)) + 
  geom_line(aes(color = a), size = 1.5) + 
  scale_color_manual(values = cb_palette) +
  labs(x = 'Weight (kg)', y = expression(bold("Total kcals"~"·"~day^-1)), 
       color = 'Alpha', title = 'Mass-Scaling Exponent (MSE)') + 
  facet_wrap(~temp_label) + 
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title.x = element_text(size = 15, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 15, face = "bold", colour = "black"),
        plot.title = element_text(size = 15, face = "bold", colour = "black", hjust = 0.5), 
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.105,0.72),
        # legend.title = element_text(size = 15, face = "bold", colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold", colour = "black"))
b

c <- df1 |> 
  mutate(temp = case_when(
    temp >= 22.5 & temp <= 22.7 ~ 22.5,
    TRUE ~ temp
  )) |> 
  filter(temp %in% c(12, 22.5, 33)) |> 
  mutate(
    temp_label = factor(
      temp,
      levels = c(12, 22.5, 33),
      labels = c("12°C", "22.5°C", "33°C")
    )
  ) |> 
  filter(a == 0.6732) |> 
  filter(act == 2.5) |> 
  filter(troph %in% c(2.5, 3.5, 4.5)) |> 
  mutate(temp = as.factor(temp),
         troph = as.factor(troph)) |> 
  ggplot(aes(x = weight_g/1000, y = kcal_day, group = troph)) + 
  geom_line(aes(color = troph), size = 1.5) + 
  scale_color_manual(values = cb_palette) +
  labs(x = 'Weight (kg)', y = expression(bold("Total kcals"~"·"~day^-1)), 
       color = 'Trophic Level', title = 'Trophic Level (TL)') + 
  facet_wrap(~temp_label) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0,100)) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title.x = element_text(size = 15, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 15, face = "bold", colour = "black"),
        plot.title = element_text(size = 15, face = "bold", colour = "black", hjust = 0.5), 
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.08,0.72),
        # legend.title = element_text(size = 15, face = "bold", colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold", colour = "black"))
c

d <- df1 |> 
  mutate(temp = case_when(
    temp >= 22.5 & temp <= 22.7 ~ 22.5,
    TRUE ~ temp
  )) |> 
  filter(temp %in% c(12, 22.5, 33)) |> 
  mutate(
    temp_label = factor(
      temp,
      levels = c(12, 22.5, 33),
      labels = c("12°C", "22.5°C", "33°C")
    )
  ) |> 
  filter(a == 0.6732) |> 
  filter(act %in% c(1.0, 2.5, 4.0)) |> 
  filter(troph == 4.5) |> 
  mutate(temp = as.factor(temp),
         act = as.factor(act)) |> 
  ggplot(aes(x = weight_g/1000, y = kcal_day, group = act)) + 
  geom_line(aes(color = act), size = 1.5) +
  scale_color_manual(values = cb_palette) +
  labs(x = 'Weight (kg)', y = expression(bold("Total kcals"~"·"~day^-1)), 
       color = 'Activity', title = 'Activity Scope (AS)') + 
  facet_wrap(~temp_label) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0,100)) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title.x = element_text(size = 15, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 15, face = "bold", colour = "black"),
        plot.title = element_text(size = 15, face = "bold", colour = "black", hjust = 0.5), 
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.08,0.72),
        # legend.title = element_text(size = 15, face = "bold", colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold", colour = "black"))
d

### patch it together ----
a <- a + theme(plot.tag.position = c(0.05, 1), plot.tag = element_text(size = 14, face = "bold"))
left <- a
left

b <- b + labs(x = NULL, y = NULL) + theme(plot.tag.position = c(0.05, 1), plot.tag = element_text(size = 14, face = "bold"))
c <- c + labs(x = NULL) + theme(strip.text = element_blank(), plot.tag.position = c(0.05, 1), plot.tag = element_text(size = 14, face = "bold"))
d <- d + labs(y = NULL) + theme(strip.text = element_blank(), plot.tag.position = c(0.05, 1), plot.tag = element_text(size = 14, face = "bold"))

final <- (left | (b / c / d)) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = 'a')
final

# ggsave('plots/2025/sensitivity-four-panel.png',
#        dpi = 600, units= 'in', height = 6, width = 12)