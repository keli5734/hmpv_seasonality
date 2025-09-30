
setwd(dir = "/Users/ke/Library/CloudStorage/OneDrive-YaleUniversity/Postdoc projects/hmpv seasonality/hmpv/two_strain_model/")

 
 

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)  
library(deSolve)
library(patchwork)   # optional, much cleaner for arranging plots
library(ggsci)

source("hmpv_transmission_model_viral_interference.R")
clean_data <-  readRDS("hmpv_data_combined.rds")
RSV_data <- readRDS("nrevss_raw_scaled_RSV.rds")

plot_list <- list()  # empty list to store plots
CG_list <- list()


#regions = c(1,4,8,9)
#regions = c(2,3,5,6,7,10)
regions = 1:10
 

beta_vec = c()
alpha_vec = c()
phi_vec = c()
interaction = c()


region_states <- c(
  "1"  = "CT, ME, MA, NH, RI, VT",
  "2"  = "NJ, NY, PR, VI",
  "3"  = "DE, DC, MD, PA, VA, WV",
  "4"  = "AL, FL, GA, KY, MS, NC, SC, TN",
  "5"  = "IL, IN, MI, MN, OH, WI",
  "6"  = "AR, LA, NM, OK, TX",
  "7"  = "IA, KS, MO, NE",
  "8"  = "CO, MT, ND, SD, UT, WY",
  "9"  = "AZ, CA, HI, NV",
  "10" = "AK, ID, OR, WA"
)

scaled_epiweek <- function(week) {
  len.week <- length(week)
  Epiweek1 <- c()
  for(i in 1:len.week){
    if(week[i] > 26){
      Epiweek1[i] <- week[i] - 26
    }
    else{
      Epiweek1[i] <- week[i] - 26 + 52
    }
    
  }
  return(Epiweek1)
}


 

for(region in regions){
  
  climate_parameter = readRDS(paste0("parm_",8,"_climate.rds"))
  
  source("parameter_setting.R")
  my_region = region
  
  data_pet <- readRDS("pet.rds") %>% dplyr::filter(region == my_region) %>% pull(pet)
  data_ppt <- readRDS("ppt.rds") %>% dplyr::filter(region == my_region) %>% pull(ppt)
  data_tmin <- readRDS("tmin.rds") %>% dplyr::filter(region == my_region) %>% pull(tmin)
  data_vap <- readRDS("vap.rds") %>% dplyr::filter(region == my_region) %>% pull(vap)
   
  parm = readRDS(paste0("parm_viral_interference_3_updated/parm_",region,"_viral_interference.rds"))
   
  RSV_incidence <- RSV_data %>%  dplyr::filter(regions == region) %>% pull(scaled_cases)

  RSV_average <- RSV_data %>%
    dplyr::filter(regions == region) %>%
    group_by(epi_week_cdc) %>%
    summarise(rsv_ave = mean(scaled_cases)) %>%
    pull(rsv_ave)

  beta_vec[region] = parm$baseline.txn.rate
  alpha_vec[region] = parm$Amp
  phi_vec[region] = parm$phi
  interaction[region] = parm$interaction

   

  # # Cumulative sum every 4 elements
  # rs <- tapply(RSV_incidence,
  #                  (seq_along(RSV_incidence) - 1) %/% 20,
  #                  cumsum)
  #
  # # Flatten into a single vector if needed
  # rs <- unlist(rs, use.names = FALSE)


  #parm_for_fit$data_rsv = c(rep(RSV_incidence, length.out = length(B_part2)), RSV_incidence)
  parm_for_fit$data_rsv = c(rep(RSV_average, length.out = length(B_part2)),RSV_incidence)
  parm_for_fit$data_rsv = (parm_for_fit$data_rsv - (min(parm_for_fit$data_rsv))) / (max(parm_for_fit$data_rsv) - min(parm_for_fit$data_rsv))
  #parm_for_fit$reporting_fraction = parm$reporting_fraction
  
  parm_for_fit$data_pet = data_pet
  parm_for_fit$data_ppt = data_ppt
  parm_for_fit$data_tmin = data_tmin
  parm_for_fit$data_vap = data_vap
  
  
  amp_pet =  0
  amp_ppt =  0  
  amp_tmin = 0
  amp_vap =  climate_parameter$amp_vap.b2 
  
  results <-
    ode(y = yinit.vector,
        t = my_times,
        func = hmpv_transmission_model_viral_interference,
        parms = c(parm_for_fit,
                  Amp = climate_parameter$Amp.b1 ,  
                  phi =  climate_parameter$phi.trans ,
                  baseline.txn.rate = parm$baseline.txn.rate,
                  f_int = parm$interaction.f * 0,
                  amp_pet = amp_pet,
                  amp_ppt = amp_ppt,
                  amp_tmin = amp_tmin,
                  amp_vap = amp_vap),
        atol = 1e-6,
        rtol = 1e-6)

  burnN <- t_burn_in
  results.burned <- results[-c(1:burnN),]


  

  delta1=c(rep(.4,12), # < 6 mos
           rep(.2,4), # 1-5 years
           rep(.05,2), # 5-18 yr
           rep(.05,2), # 18-64 yr
           rep(.2,1)) # > 65 yr



  #proportion of second infections that are LRI
  delta2=.5*delta1
  #proportion of third infections that are LRI
  delta3=.5*delta2

  q <- 1
  beta0 <- parm$baseline.txn.rate/(parm_for_fit$dur.days1/7)
  beta <- (beta0)/(sum(yinit.matrix)^(1-q))*parm_for_fit$contact
  Amp <-  parm$Amp
  phi <-  parm$phi

  t0 = nrow(results.burned)

  I1 <- results.burned[,grep('I1', colnames(results.burned))]
  I2 <- results.burned[,grep('I2', colnames(results.burned))]
  I3 <- results.burned[,grep('I3', colnames(results.burned))]
  I4 <- results.burned[,grep('I4', colnames(results.burned))]
  S0 <- results.burned[,grep('S0', colnames(results.burned))]
  S1 <- results.burned[,grep('S1', colnames(results.burned))]
  S2 <- results.burned[,grep('S2', colnames(results.burned))]
  S3 <- results.burned[,grep('S3', colnames(results.burned))]


  lambda1=matrix(0,nrow=t0,ncol=N_ages)
  for (t in 1:t0) {lambda1[t,]<-as.vector((1+Amp*cos(2*pi*(t-phi*52.1775)/52.1775) + 
                                             amp_pet * data_pet[t] + 
                                             amp_ppt * data_ppt[t] + 
                                             amp_tmin * data_tmin[t] + 
                                             amp_vap * data_vap[t])*
                                            ((I1[t,]+rho1*I2[t,]+rho2*I3[t,]+rho2*I4[t,])
                                             %*%beta)/sum(results.burned[t,], na.rm = T))}

  H1=matrix(0,nrow=t0,ncol=N_ages)
  reporting_fraction = c(rep(parm$reporting_fraction1, 12),
                        rep(parm$reporting_fraction2, 4),
                        rep(parm$reporting_fraction3, 2),
                        rep(parm$reporting_fraction4,2),
                        rep(parm$reporting_fraction5)
                        )
  
  for (i in 1:N_ages){
    H1[,i]=
      reporting_fraction[i]*(delta1[i]*S0[,i]*lambda1[,i]+
      delta2[i]*sigma1*S1[,i]*lambda1[,i]+
      delta3[i]*sigma2*S2[,i]*lambda1[,i]+
      delta3[i]*sigma3*S3[,i]*lambda1[,i])}
  
  

  H_true <- rowSums(H1,na.rm = T)


  agedist = c(sum(colSums(H1, na.rm = T)[1:12]),
              sum(colSums(H1, na.rm = T)[13:16]),
              sum(colSums(H1, na.rm = T)[17:18]),
              sum(colSums(H1, na.rm = T)[19:20]),
              sum(colSums(H1, na.rm = T)[21])) / sum(H1, na.rm = T)

  agedist = c(agedist[1], agedist[2], agedist[3],
              agedist[4] + 1*agedist[5]/8,
              7*agedist[5]/8)

  hmpv <- clean_data$ts_data %>%
    dplyr::filter(regions == region)  %>%
    pull(raw_cases)

  my_date <- clean_data$ts_data %>%
    dplyr::filter(regions == region)  %>%
    pull(date)

  hmpv2 <- clean_data$ts_data %>%
    dplyr::filter(regions == region) %>%
    pull(scaled_cases)

  age_dist <-  clean_data$age_dist %>%
    dplyr::filter(Regions == region) %>%
    pull(prop)

  my_data = list(hmpv = hmpv,
                 prop = round(age_dist * sum(hmpv)))

  agedist_vec <- agedist
  #agedist_vec <- c(agedist[1], agedist[2], agedist[3], agedist[4] + agedist[5])

  normalized_vec <- my_data$prop / sum(my_data$prop)
  #normalized_vec <- c(normalized_vec[1], normalized_vec[2], normalized_vec[3], normalized_vec[4] + normalized_vec[5])

  line_df <- data.frame(
    Date = my_date,
    Raw = hmpv,
    Rescaled = hmpv2,
    Fit = H_true #* parm$reporting_fraction,
    #RSV = RSV_incidence/6
  ) %>%
    pivot_longer(-Date, names_to = "Type", values_to = "Count") %>%
    mutate(Type = factor(Type, levels = c("Raw", "Rescaled", "Fit")))

  
  
  line_df$Type <- factor(line_df$Type, levels = c(levels(line_df$Type), "Model fit"))
  
  # Replace "Fit" with "Model fit"
  line_df$Type[line_df$Type == "Fit"] <- "Model fit"  
  
   
  
  line_df  <-  
    line_df %>% 
    mutate(epi_week_cdc = lubridate::epiweek(Date)) %>%
    mutate(epi_week_cdc = scaled_epiweek(epi_week_cdc)) %>%
    mutate(year = lubridate::year(Date))  
  
  
  
  
  line_df <- line_df %>%  mutate(season = case_when(
    (epi_week_cdc %in% 1:26 & year == 2005) | (epi_week_cdc %in% 27:52 & year == 2006)  ~ "2005-2006",
    (epi_week_cdc %in% 1:26 & year == 2006) | (epi_week_cdc %in% 27:52 & year == 2007)  ~ "2006-2007",
    (epi_week_cdc %in% 1:26 & year == 2007) | (epi_week_cdc %in% 27:52 & year == 2008)  ~ "2007-2008",
    (epi_week_cdc %in% 1:26 & year == 2008) | (epi_week_cdc %in% 27:52 & year == 2009)  ~ "2008-2009",
    (epi_week_cdc %in% 1:26 & year == 2009) | (epi_week_cdc %in% 27:52 & year == 2010)  ~ "2009-2010",
    (epi_week_cdc %in% 1:26 & year == 2010) | (epi_week_cdc %in% 27:52 & year == 2011)  ~ "2010-2011",
    (epi_week_cdc %in% 1:26 & year == 2011) | (epi_week_cdc %in% 27:52 & year == 2012)  ~ "2011-2012",
    (epi_week_cdc %in% 1:26 & year == 2012) | (epi_week_cdc %in% 27:52 & year == 2013)  ~ "2012-2013",
    (epi_week_cdc %in% 1:26 & year == 2013) | (epi_week_cdc %in% 27:52 & year == 2014)  ~ "2013-2014",
    (epi_week_cdc %in% 1:26 & year == 2014) | (epi_week_cdc %in% 27:52 & year == 2015)  ~ "2014-2015",
    (epi_week_cdc %in% 1:26 & year == 2015) | (epi_week_cdc %in% 27:52 & year == 2016)  ~ "2015-2016",
    (epi_week_cdc %in% 1:26 & year == 2016) | (epi_week_cdc %in% 27:52 & year == 2017)  ~ "2016-2017",
    (epi_week_cdc %in% 1:26 & year == 2017) | (epi_week_cdc %in% 27:52 & year == 2018)  ~ "2017-2018",
    (epi_week_cdc %in% 1:26 & year == 2018) | (epi_week_cdc %in% 27:52 & year == 2019)  ~ "2018-2019",
    (epi_week_cdc %in% 1:26 & year == 2019) | (epi_week_cdc %in% 27:52 & year == 2020)  ~ "2019-2020",
    (epi_week_cdc %in% 1:26 & year == 2020) | (epi_week_cdc %in% 27:52 & year == 2021)  ~ "2020-2021",
    (epi_week_cdc %in% 1:26 & year == 2021) | (epi_week_cdc %in% 27:52 & year == 2022)  ~ "2021-2022",
    (epi_week_cdc %in% 1:26 & year == 2022) | (epi_week_cdc %in% 27:52 & year == 2023)  ~ "2022-2023",
    (epi_week_cdc %in% 1:26 & year == 2023) | (epi_week_cdc %in% 27:52 & year == 2024)  ~ "2023-2024",
    (epi_week_cdc %in% 1:26 & year == 2024) | (epi_week_cdc %in% 27:52 & year == 2025)  ~ "2024-2025",
  )) %>% 
    dplyr:: filter(!is.na(season))
  
  
 regional_CG = line_df %>% 
   dplyr:: filter(Type != "Rescaled") %>% 
   group_by(season, Type) %>%
   summarise(CG = sum(epi_week_cdc * Count) / (sum(Count)) ) %>% 
   mutate(regions = region)
   
  
  
  p1 <- ggplot() +
    geom_line(data = line_df, aes(x = as.Date(Date), y = Count, color = Type), size = 1) +
    scale_color_manual(values = c("Raw" = "grey", "Rescaled" = "black", "Model fit" = "red", "RSV" = "transparent" )) +
    labs(x = "Date", y = "The number of positive tests", color = NULL) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top")

  age_labels <- c("<1", "1-5", "5-18", "18-65", ">65")

  bar_df <- data.frame(
    Age = factor(age_labels,levels = age_labels),
    Prediction = agedist_vec,
    Data = normalized_vec
  ) %>%
    pivot_longer(cols = c(Prediction, Data), names_to = "Source", values_to = "Proportion")
  
  bar_df$Source[bar_df$Source == "Prediction"] <- "Model fit"  

  mypal <- pal_npg("nrc", alpha = .7)(9)
  p2 <- ggplot(bar_df, aes(x = Age, y = Proportion, fill = Source)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("Model fit" = mypal[1], "Data" = mypal[4])) +
    labs(y = "Proportion", x = NULL, fill = NULL) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top")
  p2
   
  #title <- textGrob(paste0("Region ", region), gp = gpar(fontsize = 16, fontface = "bold"))
  title <- textGrob(paste0("Region ", region, ": ", region_states[as.character(region)]), gp = gpar(fontsize = 13, fontface = "bold"))

  plot_combined <- arrangeGrob(
    title,
    arrangeGrob(p1, p2, ncol = 2),
    heights = c(0.1, 1)
  )

  plot_list[[region]] <- plot_combined
  
  #CG_list[[region]] <- regional_CG
}

do.call(grid.arrange, c(compact(plot_list), nrow = 5, ncol = 2))  # 2 columns, 5 rows



#par(mfrow = c(1, 1))  # Back to single plot layout

#CG_list[[1]] %>% group_by(Type) %>% summarise(mean(CG))

# CG_combined_df <- do.call(rbind, CG_list)
# 
# CG_combined_df <- CG_combined_df %>% group_by(Type, regions) %>% summarise(x = mean(CG))
# 
# CG_wide <- CG_combined_df %>%
#   pivot_wider(
#     names_from = Type,
#     values_from = x
#   )
# 
# cor_test <- cor.test(CG_wide$Raw, CG_wide$`Model fit`, use = "complete.obs")
# 
# r_val <- cor_test$estimate
# p_val <- cor_test$p.value
# 
#  CG_wide %>% mutate(regions = factor(regions, levels = 1:10)) %>% 
#   ggplot() +
#   geom_point(aes(x = Raw, y = `Model fit`, color = regions), size = 3) +
#   theme_bw() +
#   xlab("The center of gravity of HMPV (observations)") +
#   ylab("The center of gravity of HMPV (model fit)") +
#   theme(
#     axis.text.x = element_text(color = "black", size = 10, angle = 0, vjust = 0.5),
#     axis.text.y = element_text(color = "black", size = 10, angle = 0),
#     text = element_text(size = 13)
#   ) +
#   scale_color_d3() +
#   annotate(
#     "text",
#     x = min(CG_wide$Raw -1.5, na.rm = TRUE),
#     y = max(CG_wide$`Model fit`, na.rm = TRUE),
#     label = paste0("r = ", round(r_val, 2), 
#                    ", p = ", signif(p_val, 3)),
#     hjust = 0, vjust = 1, size = 5
#   )  + 
#   ylim(32,36) + 
#   xlim(32,36)











library(png)
library(grid)
library(ggpubr)

img <- png::readPNG("Picture2.png")
Fig3A <- rasterGrob(img, interpolate=TRUE)

Fig3AB <- ggarrange(Fig3A, Fig3B,
                       nrow = 1,
                       ncol = 2,
                       widths = c(0.3, 0.7),
                       labels = c("A", "B"),
                  font.label = list(size = 20, color = "black"))

Fig3AB

  
# b <- beta_vec #/ (parm_for_fit$dur.days1/7)
# # q depends on transmission type q = 1 here
# q = parm_for_fit$q
# # c2 is the contact matrix transmission probability per unit time in each age group
# contact = parm_for_fit$contact
# # transmission rate
# beta <- numeric(10)  # preallocate vector of length 10
# for (i in c(1:10)) {
#   beta[i] <- max(eigen(b[i] * contact / (parm_for_fit$dur.days1))$values)
# }

beta = readRDS("R0.rds")
alpha_vec = readRDS("alpha_vec.rds")



clean_data <-  readRDS("hmpv_data_combined.rds")

hmpv_TS <- clean_data$ts_data %>% 
  mutate(even_odd = ifelse(season %in% c("2008-2009",
                                         "2010-2011", 
                                         "2012-2013", 
                                         "2014-2015",
                                         "2016-2017", 
                                         "2018-2019"), "even", "odd") )

 
hmpv_TS_mean_region <- hmpv_TS %>% 
  group_by(regions, epi_week_cdc) %>% 
  summarise(mean_cases = mean(scaled_cases))


 


hmpv_CG <- hmpv_TS %>% 
  #dplyr::filter(even_odd == "even") %>% 
  group_by(regions, season) %>% 
  summarise(gravity = sum(epi_week_cdc * scaled_cases) / (sum(scaled_cases)) ) %>% 
  ungroup()

 


# 1. Mean gravity per region
mean_by_region_hmpv <- hmpv_CG %>% 
  #dplyr::filter(season == "2009-2010") %>% 
  group_by(regions) %>% 
  summarise(mean_gravity = mean(gravity, na.rm = TRUE)) %>% 
  mutate(regions = as.character(regions))  %>% 
  ungroup() 

mean_by_region_hmpv$mean_gravity = beta

# 2. Map HHS regions to states
region_states <- list(
  "1" = c("CT","ME","MA","NH","RI","VT"),
  "2" = c("NJ","NY"),
  "3" = c("DE","DC","MD","PA","VA","WV"),
  "4" = c("AL","FL","GA","KY","MS","NC","SC","TN"),
  "5" = c("IL","IN","MI","MN","OH","WI"),
  "6" = c("AR","LA","NM","OK","TX"),
  "7" = c("IA","KS","MO","NE"),
  "8" = c("CO","MT","ND","SD","UT","WY"),
  "9" = c("AZ","CA","HI","NV"),
  "10"= c("AK","ID","OR","WA")
)

state_region <- tibble(
  regions = rep(names(region_states), times = sapply(region_states, length)),
  abbr    = unlist(region_states)
)

# 3. Join means and rename to required column
plot_df <- state_region %>% 
  left_join(mean_by_region_hmpv, by = "regions") %>% 
  rename(state = abbr)

plot_long <- plot_df %>%
  pivot_longer(
    cols = c(mean_gravity),
    names_to = "virus",
    values_to = "mean_gravity"
  ) %>%
  mutate(
    virus = recode(virus,
                   mean_gravity = "HMPV")
                    
  )  
source("Fig2A.R")

 Fig3C  = plot_usmap(data = plot_long, values = "mean_gravity", color = "NA") +
   scale_fill_gradient(
    name = "The basic \n reproduction number",
    low = "#F1F1F1FF",          # start
    #mid = "#B695BCFF",
    high = "#172869FF", 
    limits = c(6.4, 7.6),
    breaks = seq(6.4, 7.6, .2),
    labels = seq(6.4, 7.6, .2),
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold")
  )+ 
   geom_sf(data = st_boundary(hhs_polys),
           color = "black", linewidth = 1.2, inherit.aes = F)
  

 library(tidyverse)
library(usmap)
library(viridis)
library(signal)    # for periodogram if needed

# Set frequency in weeks
annual_freq   <- 1 / 52
biennial_freq <- 1 / 104

# Helper function to get amplitude ratio from Fourier spectrum
get_fourier_ratio <- function(x, freq_per_year = 52.18) {
  n <- length(x)
  fft_out <- fft(x - mean(x))  # center the signal
  amp <- Mod(fft_out)[1:(n/2)]  # only positive frequencies
  
  freqs <- (0:(n/2 - 1)) / n * freq_per_year  # in cycles/year
  
  # Find nearest bins to 1/year and 1/2 year
  idx_annual   <- which.min(abs(freqs - 1))
  idx_biennial <- which.min(abs(freqs - 0.5))
  
  ratio <- amp[idx_biennial] / amp[idx_annual]
  return(ratio)
}

# Apply per region
fourier_ratios.hmpv <- hmpv_TS %>%
  group_by(regions) %>%
  summarise(ratio = get_fourier_ratio(raw_cases), .groups = "drop")

fourier_ratios.hmpv$ratio = alpha_vec

plot_df2 <- state_region %>%
  left_join(fourier_ratios.hmpv, by = "regions") %>%
  rename(state = abbr)


plot_long2 <- plot_df2 %>%
  pivot_longer(
    cols = c(ratio),
    names_to = "virus",
    values_to = "ratio"
  )  

Fig3D =  plot_usmap(data = plot_long2, values = "ratio", color = "NA") +
  scale_fill_gradient(
    name = "Seasonal amplitude",
    #colours = c("#4575b4", "#74add1", "#dfc27d", "#8c510a"),
    low = "#F1F1F1FF",          # start
    high = "#01665EFF",   
    limits = c(0.1, 0.24),
    breaks = seq(0.1, 0.24, 0.04),
    na.value = "grey90"  # optional, for missing data
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  theme(strip.text = element_blank()) + 
  geom_sf(data = st_boundary(hhs_polys),
          color = "black", linewidth = 1.2, inherit.aes = F)

Fig3D
Fig3CD  = ggarrange(Fig3C, Fig3D,
                              nrow = 1,
                              ncol = 2,
                              widths = c(0.5, 0.5),
                              labels = c("C", "D"),
                              font.label = list(size = 20, color = "black"))

Fig3CD

Fig3  = ggarrange(Fig3AB, Fig3CD,
                 nrow = 2,
                 ncol = 1,
                 heights = c(0.7, 0.3),
                 labels = c("", ""),
                 font.label = list(size = 20, color = "black"))

Fig3
#ggsave("./Fig3.tiff", width = 16, height = 9, units = "in", bg = "white", dpi = 200)

# for(region in regions){
# 
#   source("parameter_setting.R")
#   my_region = region
# 
#   parm = readRDS(paste0("parm_viral_interference_3/parm_",region,"_viral_interference.rds"))
# 
#   RSV_incidence <- RSV_data %>%
#     dplyr::filter(regions == region) %>% pull(scaled_cases)
# 
#   RSV_average <- RSV_data %>%
#     dplyr::filter(regions == region) %>%
#     group_by(epi_week_cdc) %>%
#     summarise(rsv_ave = mean(scaled_cases)) %>%
#     pull(rsv_ave)
# 
#   # Normalize RSV input
#   parm_for_fit$data_rsv = c(rep(RSV_average, length.out = length(B_part2)), RSV_incidence)
#   parm_for_fit$data_rsv = (parm_for_fit$data_rsv - min(parm_for_fit$data_rsv)) /
#     (max(parm_for_fit$data_rsv) - min(parm_for_fit$data_rsv))
#   parm_for_fit$reporting_fraction = parm$reporting_fraction
# 
#   ## --- Function to run model ---
#   run_model <- function(f_int_value, label){
#     results <- ode(
#       y = yinit.vector,
#       t = my_times,
#       func = hmpv_transmission_model_viral_interference,
#       parms = c(parm_for_fit,
#                 Amp = parm$Amp,
#                 phi = parm$phi,
#                 baseline.txn.rate = parm$baseline.txn.rate,
#                 f_int = f_int_value),
#       atol = 1e-6, rtol = 1e-6
#     )
# 
#     results.burned <- results[-c(1:t_burn_in),]
# 
#     I1 <- results.burned[,grep('I1', colnames(results.burned))]
#     I2 <- results.burned[,grep('I2', colnames(results.burned))]
#     I3 <- results.burned[,grep('I3', colnames(results.burned))]
#     I4 <- results.burned[,grep('I4', colnames(results.burned))]
#     S0 <- results.burned[,grep('S0', colnames(results.burned))]
#     S1 <- results.burned[,grep('S1', colnames(results.burned))]
#     S2 <- results.burned[,grep('S2', colnames(results.burned))]
#     S3 <- results.burned[,grep('S3', colnames(results.burned))]
# 
#     t0 = nrow(results.burned)
#     lambda1 = matrix(0,nrow=t0,ncol=N_ages)
#     for (t in 1:t0) {
#       lambda1[t,] <- as.vector((1+parm$Amp*cos(2*pi*(t - parm$phi*52.1775)/52.1775)) *
#                                  ((I1[t,]+rho1*I2[t,]+rho2*I3[t,]+rho2*I4[t,]) %*% parm_for_fit$contact) /
#                                  sum(results.burned[t,], na.rm = T))
#     }
# 
#     H1 = matrix(0,nrow=t0,ncol=N_ages)
#     for (i in 1:N_ages){
#       H1[,i] = delta1[i]*S0[,i]*lambda1[,i] +
#         delta2[i]*sigma1*S1[,i]*lambda1[,i] +
#         delta3[i]*sigma2*S2[,i]*lambda1[,i] +
#         delta3[i]*sigma3*S3[,i]*lambda1[,i]
#     }
#     H_true <- rowSums(H1, na.rm = T)
# 
#     data.frame(Date = my_date,
#                Simulation = H_true * parm$reporting_fraction,
#                Scenario = label)
#   }
# 
#   ## --- Run both scenarios ---
#   df_fit   <- run_model(parm$interaction, "Fit")
#   df_null  <- run_model(0, "NoInteraction")
# 
#   ## --- Combine with observed ---
#   line_df <- data.frame(
#     Date = my_date,
#     Raw = clean_data$ts_data %>% dplyr::filter(regions==region) %>% pull(raw_cases),
#     Rescaled = clean_data$ts_data %>% dplyr::filter(regions==region) %>% pull(scaled_cases),
#     RSV = RSV_incidence/6
#   ) %>%
#     pivot_longer(-Date, names_to = "Type", values_to = "Count") %>%
#     mutate(Type = factor(Type, levels = c("Raw","Rescaled","RSV")))
# 
#   sim_df <- bind_rows(
#     df_fit %>% rename(Count = Simulation) %>% mutate(Type = Scenario),
#     df_null %>% rename(Count = Simulation) %>% mutate(Type = Scenario)
#   )
# 
#   line_df <- bind_rows(line_df, sim_df)
# 
#   ## --- Plot ---
#   p1 <- ggplot() +
#     geom_line(data = line_df, aes(x = as.Date(Date), y = Count, color = Type), size = 1) +
#     scale_color_manual(values = c("Raw"="grey", "Rescaled"="black",
#                                   "RSV"="transparent", "Fit"="red",
#                                   "NoInteraction"="blue")) +
#     labs(x = "Time", y = "Count", color = NULL) +
#     theme_minimal(base_size = 13) +
#     theme(legend.position = "top")
# 
#   print(p1)  ## Or save into plot_list
# }
