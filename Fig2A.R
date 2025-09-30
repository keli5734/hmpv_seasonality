library(dplyr)
library(usmap)
library(ggplot2)

# approximate label positions for each HHS region
# approximate label positions for each HHS region
hhs_labels <- tibble::tribble(
  ~region, ~lon,  ~lat,
  1, -71, 44,   # New England
  2, -75, 43,   # NY/NJ
  3, -78, 40,   # Mid-Atlantic
  4, -84, 34,   # Southeast
  5, -89, 43,   # Great Lakes
  6, -97, 32,   # South Central
  7, -96, 40,   # Central Plains
  8, -105, 45,  # Rockies
  9, -118, 36,  # California/NV/AZ
  10, -120, 46   # Pacific NW
)

# transform coords to match plot_usmap projection
hhs_labels <- usmap::usmap_transform(hhs_labels)

# extract x/y from sf geometry and drop geometry column
hhs_labels <- hhs_labels %>%
  mutate(
    x = sf::st_coordinates(.)[, 1],
    y = sf::st_coordinates(.)[, 2]
  ) %>%
  sf::st_drop_geometry()


#hhs_polys <- usmap::usmap_transform(hhs_polys)  # convert to usmap projection if needed

 
library(dplyr)
library(sf)
library(tigris)
library(units)
library(usmap)

# --- mapping table for states → HHS region ---
hhs_map <- tibble::tribble(
  ~state, ~region,
  "CT",1, "ME",1, "MA",1, "NH",1, "RI",1, "VT",1,
  "NJ",2, "NY",2, 
  "DE",3, "DC",3, "MD",3, "PA",3, "VA",3, "WV",3,
  "AL",4, "FL",4, "GA",4, "KY",4, "MS",4, "NC",4, "SC",4, "TN",4,
  "IL",5, "IN",5, "MI",5, "MN",5, "OH",5, "WI",5,
  "AR",6, "LA",6, "NM",6, "OK",6, "TX",6,
  "IA",7, "KS",7, "MO",7, "NE",7,
  "CO",8, "MT",8, "ND",8, "SD",8, "UT",8, "WY",8,
  "AZ",9, "CA",9, "HI",9, "NV",9,
  "AK",10, "ID",10, "OR",10, "WA",10
)

# --- build region polygons and drop tiny islands ---
hhs_polys <- states(cb = TRUE, year = 2022) %>%
  select(STUSPS, geometry) %>%
  left_join(hhs_map, by = c("STUSPS" = "state")) %>%
  dplyr:: mutate(area = set_units(st_area(.), km^2)) %>%
  dplyr:: filter(area > set_units(17000, km^2)) %>%   
  group_by(region) %>% 
  summarise() %>%         # union states per region
  usmap_transform() %>%       
  dplyr:: filter(!is.na(region))

# # --- draw only the borders ---
# Fig2A <- Fig2A +
#   geom_sf(data = st_boundary(hhs_polys),
#           color = "black", linewidth = 1.2, inherit.aes = F)
