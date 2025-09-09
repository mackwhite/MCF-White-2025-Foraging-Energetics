###########################################################################
# Script Information ------------------------------------------------------
###########################################################################

### project: MCF Special Issue - Snook in Americas - Energetic Needs
### author(s): MW, WRJ, and JL
### goal: simulating snook prey (individual) needs (for fig 4)

# Housekeeping ------------------------------------------------------------

### load necessary packages -----
# install.packages("librarian")

librarian::shelf(tidyverse, readxl, readr, performance, ggpubr, 
                 mgcv, glmmTMB, MuMIn, corrplot, ggeffects)

# read in necessary data and wrangle/summarize ----------------------------
rates <- read_csv("data/manuscript/metabolic-rates-temp15thru14by5.csv") |> 
  filter(act == 2.5,
         troph == 4.5)

prey <- read_csv("data/manuscript/prey-kj-data.csv") |> 
  rename(Jpi_m = Jpi)

diet <- c('100% fish/0% invert', '75% fish/25% invert', 
          '50% fish/50% invert', '25% fish/75% invert', 
          '0% fish/100% invert')

df <- tibble(rates) |> 
  slice(rep(1:n(), each=length(diet))) |> 
  mutate(diet = rep(diet, times=nrow(rates)))

# Assuming df is your first dataset and prey is your second dataset
library(dplyr)

# Function to generate 100 new rows based on diet
generate_rows <- function(row, prey) {
  # Convert the row (list) to a data frame
  row_df <- as.data.frame(as.list(row))
  
  # Initialize a list to store the new rows
  new_rows <- list()
  
  for (i in 1:100) {
    if (row_df$diet == '100% fish/0% invert') {
      # Sample 100 times from Vertebrate
      sampled_Jpi <- sum(sample(prey$Jpi_m[prey$v_class == "Vertebrate"], 100, replace = TRUE))
    } else if (row_df$diet == '75% fish/25% invert') {
      # Sample 75 times from Vertebrate and 25 times from Invertebrate
      sampled_Jpi_vert <- sample(prey$Jpi_m[prey$v_class == "Vertebrate"], 75, replace = TRUE)
      sampled_Jpi_invert <- sample(prey$Jpi_m[prey$v_class == "Invertebrate"], 25, replace = TRUE)
      sampled_Jpi <- sum(c(sampled_Jpi_vert, sampled_Jpi_invert))
    } else if (row_df$diet == '50% fish/50% invert') {
      # Sample 50 times from Vertebrate and 50 times from Invertebrate
      sampled_Jpi_vert <- sample(prey$Jpi_m[prey$v_class == "Vertebrate"], 50, replace = TRUE)
      sampled_Jpi_invert <- sample(prey$Jpi_m[prey$v_class == "Invertebrate"], 50, replace = TRUE)
      sampled_Jpi <- sum(c(sampled_Jpi_vert, sampled_Jpi_invert))
    } else if (row_df$diet == '25% fish/75% invert') {
      # Sample 25 times from Vertebrate and 75 times from Invertebrate
      sampled_Jpi_vert <- sample(prey$Jpi_m[prey$v_class == "Vertebrate"], 25, replace = TRUE)
      sampled_Jpi_invert <- sample(prey$Jpi_m[prey$v_class == "Invertebrate"], 75, replace = TRUE)
      sampled_Jpi <- sum(c(sampled_Jpi_vert, sampled_Jpi_invert))
    } else if (row_df$diet == '0% fish/100% invert') {
      # Sample 100 times from Invertebrate
      sampled_Jpi <- sum(sample(prey$Jpi_m[prey$v_class == "Invertebrate"], 100, replace = TRUE))
    }
    
    # Create a new row with the original data and the new Jpi_prey column
    new_row <- row_df %>%
      mutate(Jpi_prey = sampled_Jpi)
    
    # Add the new row to the list of new rows
    new_rows[[i]] <- new_row
  }
  
  # Combine the list of new rows into a data frame
  bind_rows(new_rows)
}

# Apply the function to each row of df and combine the results
new_df <- df |> 
  rowwise() |> 
  do(generate_rows(., prey))

# View the resulting dataset
head(new_df)

df_long <- new_df |> 
  mutate(Jpi_prey1 = Jpi_prey/100,
         num_needs = Total_metabolic_rate_j_d/Jpi_prey1,
         diet = factor(diet,
                       levels = c("0% fish/100% invert", "25% fish/75% invert",
                                  "50% fish/50% invert", "75% fish/25% invert",
                                  "100% fish/0% invert")))

df_summ <- df_long |> 
  mutate(weight_kg = weight_g/1000) |> 
  group_by(temp, weight_kg, diet) |> 
  summarize(mean_num_needs = mean(num_needs),
            sd_num_needs = sd(num_needs),
            sd_low = mean_num_needs - sd_num_needs,
            sd_high = mean_num_needs + sd_num_needs,
            ci_low = quantile(num_needs, probs = 0.025),
            ci_high = quantile(num_needs, probs = 0.975))

df_summ |> 
  ggplot(aes(weight_kg, mean_num_needs)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = diet),
              alpha = 0.5)+
  geom_line(aes(color = diet)) +
  scale_color_viridis_d(option = 'turbo', end = 1) +
  scale_fill_viridis_d(option = 'turbo', end = 1) +
  facet_wrap(~temp) +
  labs(x = "Weight (kg)", y = "Individual Prey Needs (#)", color = 'Diet Type', fill = 'Diet Type') +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.text = element_text(face = 'bold'),
        legend.title = element_text(face = 'bold'))

###saving for publication
# ggsave("plots/ind-prey-needs-randomsampling.tiff", units = "in", width = 16,
#        height = 6, dpi =  600, compression = "lzw")

# testing out new version of figure six -----------------------------------

df_summ_new <- df_summ |> 
  filter(temp %in% c(15, 25, 35)) |> 
  mutate(weight_kg = round(weight_kg, digits = 3)) |> 
  filter(weight_kg %in% c(0.509, 3.282, 6.174)) |> 
  mutate(weight_kg = round(weight_kg, digits = 1)) |> 
  filter(weight_kg %in% c(0.5, 3.3, 6.2))

write_csv(df_summ_new,"data/prey-numeric-needs-2025.csv")

b <- ggplot(df_summ_new, aes(x = as.factor(temp), y = mean_num_needs, color = diet)) +
  geom_point(size = 3, position = position_dodge(width = 1)) +                             
  scale_color_viridis_d(option = 'turbo', end = 1) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),  
                width = 0,
                linewidth = 1,
                position = position_dodge(width = 1)) +   
  facet_wrap(~weight_kg) + 
  labs(x = "Water Temperature (°C)", y = "Prey Demand (#/day)", color = 'Diet Type', fill = 'Diet Type') +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
        panel.grid.major = element_blank(),
        # axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        # panel.border = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 12, color = "black", face = 'bold'),
        legend.title = element_text(size = 12, color = "black", face = 'bold')) 

figure <- ggarrange(a, b,
          labels = c('a)','b)'),
          ncol = 1, hjust = -0.5, vjust = 1.5, align = "h", legend = 'right',
          common.legend = TRUE)

figure2 <- annotate_figure(figure, 
                           bottom = text_grob(label = "Water Temperature (°C)",
                           just = 'centre', color = "black", face = "bold",
                           size = 15))

#saving for publication
# ggsave("plots/fig6-two-panel-prey-needs.tiff", units = "in", width = 12,
#        height = 10, dpi =  600, compression = "lzw")
