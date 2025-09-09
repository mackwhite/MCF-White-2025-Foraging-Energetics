###########################################################################
# Script Information ------------------------------------------------------
###########################################################################

### project: MCF Special Issue - Snook in Americas - Energetic Needs
### author(s): MW, WRJ, and JL
### goal: simulate snook population for caloric need estimate integration

# Housekeeping ------------------------------------------------------------

### load necessary packages -----
# install.packages("librarian")
# devtools::install_github("nschiett/fishflux", dependencies=TRUE)
# library(fishflux)
librarian::shelf(tidyverse, readxl, readr, fishflux)

### TL range seen in POR USACE-MAP data for Common Snook
size <- data.frame(seq(from = 25, to = 90, length.out = 50)) |> 
  rename(TL = seq.from...25..to...90..length.out...50.) |> 
  mutate(fl_cm = -14.8606 + 0.9512*TL,
         fl_mm = fl_cm*10)

fork_length_range <- data.frame(
  seq(from = 10, to = 900, length.out = 90)
) |> 
  rename(fl_mm = seq.from...10..to...900..length.out...90.)

# a = 10^-5.0821 = 8.277515e-06
# b = 3.1117

dt <- fork_length_range |> 
  mutate(log10_weight_g = -5.0821 + 3.1117*log10(fl_mm),
         weight_g = 10^log10_weight_g) |> 
  select(-log10_weight_g)

### set parameters from Taylor et al 2000
### asymptotic length estimates
L_asymp_high <- 947.3 + 32.15
L_asymp <- 947.3
L_asymp_low <- 947.3 - 32.15
K_high <- 0.175 + 0.0155
K <- 0.175
K_low <- 0.175 - 0.0155
t0_high <- -1.352 + 0.1714
t0 <- -1.352
t0_low <- -1.352 - 0.1714

dt_1 <- dt |> 
  mutate(age = t0 - 1/K*log(1-fl_mm/L_asymp),
         age_low = t0_high - 1/K_high*log(1-fl_mm/L_asymp_high),
         age_high = t0_low - 1/K_low*log(1-fl_mm/L_asymp_low)) |> 
  mutate(age_days = age*365.25) |> 
  filter(fl_mm >= 210) |> 
  select(-age_low, -age_high) |> 
  mutate(id = row_number()) |> 
  select(id, fl_mm, age)

date_sequence <- seq(as.Date("2021-01-01"), as.Date("2021-12-31"), by = "day")

expanded_data <- dt_1 |> 
  tidyr::crossing(date = date_sequence)  |> 
  mutate(
    days_since_start = as.numeric(date - min(date_sequence)),
    additional_years = days_since_start / 365.25,
    age_in_years = age + additional_years,
    age_in_days = age_in_years*365.25
  )

calculate_length <- function(L_inf, K, t0, t) {
  L_inf * (1 - exp(-K * (t - t0)))
}

expanded_data <- expanded_data |> 
  mutate(fl_mm_vbf = calculate_length(L_asymp, K, t0, age_in_years)) |> 
  mutate(log10_weight_g = -5.3564 + 3.1117*log10(fl_mm_vbf),
         weight_g = 10^log10_weight_g) |> 
  select(id,date,age_in_years,age_in_days,fl_mm_vbf,fl_mm,weight_g) |> 
  rename(starting_size = fl_mm)


df <- expanded_data |> 
  arrange(id, date) |> 
  group_by(id) |> 
  mutate(daily_growth_g = weight_g - lag(weight_g, default = first(weight_g))) |> 
  filter(daily_growth_g > 0) |> 
  mutate(weight_kg = weight_g/1000) |> 
  ungroup()

short_pop <- df |> 
  rename(fl_mm_rounded = starting_size) |> 
  group_by(id) |> 
  slice(1) |> 
  ungroup()

# read out data -----------------------------------------------
# write_csv(short_pop, "data/manuscript/simulated-snook-population-short.csv")