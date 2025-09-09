###########################################################################
# Script Information ------------------------------------------------------
###########################################################################

### project: MCF Special Issue - Snook in Americas - Energetic Needs
### author(s): MW
### goal: generate map for manuscript (fig 1)

# Housekeeping ------------------------------------------------------------

### load necessary packages -----
# install.packages("librarian")
librarian::shelf(tidyverse, readr, janitor, zoo, 
                 lubridate, openintro, maps, ggmap,
                 ggthemes, shapefiles, broom, sf, ggspatial, maptools, 
                 ggsn, GISTools)
### set theme ---
theme_set(theme_minimal())

### set up google map key ---
api_secret <- 'TOPSECRET :P'
register_google(key = api_secret)
has_google_key()

### read in necessary data ---
dat <- read_csv("mapping/efishing_monitoringstation_coords.csv") |> 
  janitor::clean_names() |> 
  filter(type %in% c('Marsh Stage', 'Temperature', 'Electrofishing')) |> 
  filter(site_name %in% c("RB8B1", "RB8B2", "RB8B3", 
                          "RB9B1", "RB9B2", "RB9B3",
                          "RB10B1", "RB10B2", "RB10B3",
                          "RB11B1", "RB11B2", "RB11B3",
                          "RB12B1", "RB12B2", "RB12B3",
                          "RB13B1", "RB13B2", "RB13B3",
                          'M0215', 'Bottle Creek')) |> 
  mutate(type = case_when(
    type == 'Electrofishing' ~ 'Diet Sampling',
    type == 'Temperature' ~ 'Water Temperature',
    TRUE ~ type
  ))

layer_three <- get_map(
  ### determine bounding box: https://www.openstreetmap.org/#map=5/25.304/-69.412
  c(left = -80.96, bottom = 25.39, right = -80.82, top = 25.51),
  maptype = 'satellite',
  source = 'google',
  api_key = api_secret
)

ggmap(layer_three) +
  geom_point(data = dat,
             aes(x = longitude, y = latitude, color = 'white'),
             size = 5.5,
             alpha = 1.0) +
  geom_point(data = dat,
             aes(x = longitude, y = latitude, color = type),
             size = 3.5) +
  scale_color_manual(values = c('Diet Sampling' = '#ffcc00',
                                'Water Temperature' = '#e41a1c',
                                'Marsh Stage' = '#377eb8')) +
  scale_y_continuous(limits =c(25.420,25.475)) +
  annotation_north_arrow(location = 'tl',
                         style = north_arrow_fancy_orienteering(text_col = 'white',
                                                                fill = 'white',
                                                                line_col = 'white',
                                                                text_face = "bold",
                                                                text_size = 18)) +
  annotation_scale(location = 'br', width_hint = 0.5, text_cex = 1.25,
                   text_face = "bold", text_col = 'white') +
  coord_sf(crs = 4326) +
  theme_map() +
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face = 'bold',
                                   color = 'white',
                                   size = 16))

### save for publication
ggsave("mapping/map-layer-one.tiff", units = "in", width = 12,
       height = 6, dpi =  600, compression = "lzw")

layer_two <- get_map(
  ### determine bounding box: https://www.openstreetmap.org/#map=5/25.304/-69.412
  c(left = -81.5804, bottom = 24.9163, right = -80.4268, top = 25.8147),
  maptype = 'satellite',
  source = 'google',
  api_key = api_secret
)

ggmap(layer_two) +
  theme_map() +
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face = 'bold',
                                   color = 'white',
                                   size = 16))

ggsave("mapping/test_map_layer_two.tiff", units = "in", width = 12,
       height = 6, dpi =  600, compression = "lzw")
