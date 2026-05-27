
setwd(dir = "/Users/ke/Library/CloudStorage/OneDrive-YaleUniversity/Postdoc projects/hmpv seasonality/hmpv/epic_data_cleaning/")

library(factoextra)
library(ggdendro)
library(ISOweek)
library(readxl)
library(tidyverse)
library(viridis)  
library(zoo)
library(openxlsx)
library(lubridate)
library(MASS)
library(wesanderson)
library(RColorBrewer)
library(proxy)
library(dendextend)
library(ggdendro)
library(plotly)
library(gridExtra)
library(grid)
library(lattice)
library(cowplot)
library(ggpubr)


source('epic_fraction_import.R')
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

#EPIC ED all cause

file1 <- "hmpv_test_data/2025_08_26_hMPVtests_RespEDvisits_State_Week_Age_2020.csv"
file2 <- "hmpv_test_data/2025_08_26_hMPVtests_RespEDvisits_State_Week_Age_2021.csv"
file3 <- "hmpv_test_data/2025_08_26_hMPVtests_RespEDvisits_State_Week_Age_2022.csv"
file4 <- "hmpv_test_data/2025_08_26_hMPVtests_RespEDvisits_State_Week_Age_2023.csv"
file5 <- "hmpv_test_data/2025_08_26_hMPVtests_RespEDvisits_State_Week_Age_2024.csv"
file6 <- "hmpv_test_data/2025_08_26_hMPVtests_RespEDvisits_State_Week_Age_2025.csv"

file7 <- "hmpv_test_data/2025_08_26_hMPVtests_RespEDvisits_State_Week_Age_2016.csv"
file8 <- "hmpv_test_data/2025_08_26_hMPVtests_RespEDvisits_State_Week_Age_2017.csv"
file9 <- "hmpv_test_data/2025_08_26_hMPVtests_RespEDvisits_State_Week_Age_2018.csv"
file10 <- "hmpv_test_data/2025_08_26_hMPVtests_RespEDvisits_State_Week_Age_2019.csv"

epic_ed_1 <- epic_fraction_import(ds_name = file1) 
epic_ed_2 <- epic_fraction_import(ds_name = file2) 
epic_ed_3 <- epic_fraction_import(ds_name = file3) 
epic_ed_4 <- epic_fraction_import(ds_name = file4) 
epic_ed_5 <- epic_fraction_import(ds_name = file5) 
epic_ed_6 <- epic_fraction_import(ds_name = file6) 

epic_ed_all <- rbind(epic_ed_1, epic_ed_2, epic_ed_3, epic_ed_4,epic_ed_5,epic_ed_6) %>% 
  mutate(age = factor(age, 
                      levels = c("Less than 1", 
                                 "≥ 1 and < 5", 
                                 "≥ 5 and < 18", 
                                 "≥ 18 and < 50",
                                 "≥ 50 and < 65",
                                 "65 or more")))
#   mutate(epi_week_cdc = lubridate::epiweek(date)) %>%
#   mutate(epi_week_cdc = scaled_epiweek(epi_week_cdc)) %>%
#   mutate(year = year(date))  %>%
#   mutate(season = case_when(
#     (epi_week_cdc %in% 1:26 & year == 2008) | (epi_week_cdc %in% 27:52 & year == 2009)  ~ "2008-2009",
#     (epi_week_cdc %in% 1:26 & year == 2009) | (epi_week_cdc %in% 27:52 & year == 2010)  ~ "2009-2010",
#     (epi_week_cdc %in% 1:26 & year == 2010) | (epi_week_cdc %in% 27:52 & year == 2011)  ~ "2010-2011",
#     (epi_week_cdc %in% 1:26 & year == 2011) | (epi_week_cdc %in% 27:52 & year == 2012)  ~ "2011-2012",
#     (epi_week_cdc %in% 1:26 & year == 2012) | (epi_week_cdc %in% 27:52 & year == 2013)  ~ "2012-2013",
#     (epi_week_cdc %in% 1:26 & year == 2013) | (epi_week_cdc %in% 27:52 & year == 2014)  ~ "2013-2014",
#     (epi_week_cdc %in% 1:26 & year == 2014) | (epi_week_cdc %in% 27:52 & year == 2015)  ~ "2014-2015",
#     (epi_week_cdc %in% 1:26 & year == 2015) | (epi_week_cdc %in% 27:52 & year == 2016)  ~ "2015-2016",
#     (epi_week_cdc %in% 1:26 & year == 2016) | (epi_week_cdc %in% 27:52 & year == 2017)  ~ "2016-2017",
#     (epi_week_cdc %in% 1:26 & year == 2017) | (epi_week_cdc %in% 27:52 & year == 2018)  ~ "2017-2018",
#     (epi_week_cdc %in% 1:26 & year == 2018) | (epi_week_cdc %in% 27:52 & year == 2019)  ~ "2018-2019",     
#     (epi_week_cdc %in% 1:26 & year == 2019) | (epi_week_cdc %in% 27:52 & year == 2020)  ~ "2019-2020",  
#     (epi_week_cdc %in% 1:26 & year == 2020) | (epi_week_cdc %in% 27:52 & year == 2021)  ~ "2020-2021",  
#     (epi_week_cdc %in% 1:26 & year == 2021) | (epi_week_cdc %in% 27:52 & year == 2022)  ~ "2021-2022",   
#     (epi_week_cdc %in% 1:26 & year == 2022) | (epi_week_cdc %in% 27:52 & year == 2023)  ~ "2022-2023",  
#     (epi_week_cdc %in% 1:26 & year == 2023) | (epi_week_cdc %in% 27:52 & year == 2024)  ~ "2023-2024",
#     (epi_week_cdc %in% 1:26 & year == 2024) | (epi_week_cdc %in% 27:52 & year == 2025)  ~ "2024-2025" )) %>% 
#    filter(!is.na(season))   %>% 
#   ungroup()



age_colors <- c(
  "Less than 1"     = "#02818a",
  "≥ 1 and < 5"     = "#3690c0",
  "≥ 5 and < 18"    = "#67a9cf",
  "≥ 18 and < 50"   = "#a6bddb",
  "≥ 50 and < 65"   = "#d0d1e6",
  "65 or more"      = "#ece2f0"
)
# 
# epic_ed_all %>% 
#   mutate(date = my(paste(month, year))) %>% 
#   filter(geography %in% c("Illinois")) %>% 
#   ggplot(aes(x = date, y = N_cases, color = age)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   facet_wrap(~ age, scales = "free_y") +
#   labs(
#     title = "Age Distribution of Cases Over Time by State",
#     x = "Date",
#     y = "Number of Cases",
#     fill = "Age Group"
#   ) +
#   scale_color_manual(values = age_colors)+
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "none"   
#   )


region1 <- c("Connecticut", "Maine", "Massachusetts", "New Hampshire", "Rhode Island", "Vermont")
region2 <- c("New Jersey", "New York", "Puerto Rico", "U.S. Virgin Islands")
region3 <- c("Delaware", "District of Columbia", "Maryland", "Pennsylvania", "Virginia", "West Virginia")
region4 <- c("Alabama", "Florida", "Georgia", "Kentucky", "Mississippi", "North Carolina", "South Carolina", "Tennessee")
region5 <- c("Illinois", "Indiana", "Michigan", "Minnesota", "Ohio", "Wisconsin")
region6 <- c("Arkansas", "Louisiana", "New Mexico", "Oklahoma", "Texas")
region7 <- c("Iowa", "Kansas", "Missouri", "Nebraska")
region8 <- c("Colorado", "Montana", "North Dakota", "South Dakota", "Utah", "Wyoming")
region9 <- c("Arizona", "California", "Hawaii", "Nevada", "American Samoa", "Guam", "Northern Mariana Islands", 
             "Federated States of Micronesia", "Marshall Islands", "Palau")
region10 <- c("Alaska", "Idaho", "Oregon", "Washington")


age_distribution <- epic_ed_all %>% 
  mutate(Regions = ifelse(geography %in% region1, "1", 
                          ifelse(geography %in% region2, "2", 
                                 ifelse(geography %in% region3, "3", 
                                        ifelse(geography %in% region4, "4",
                                               ifelse(geography %in% region5, "5", 
                                                      ifelse(geography %in% region6, "6", 
                                                             ifelse(geography %in% region7, "7", 
                                                                    ifelse(geography %in% region8, "8",
                                                                           ifelse(geography %in% region9, "9", "10")))))))))) %>% 
  dplyr::filter(!is.na(Regions))  %>% 
  mutate(age = ifelse(age == "Less than 1", "< 1", 
                      ifelse(age == "≥ 1 and < 5", "[1,5)", 
                             ifelse(age == "≥ 5 and < 18", "[5,18)", 
                                    ifelse(age == "65 or more", ">= 65", "[18,65)")))))  %>% 
  group_by(Regions, age) %>% 
  summarise(hmpv_tests_average = mean(as.numeric(gsub("%", "", hmpv_tests)), na.rm = T),
            positive_tests_per = mean(as.numeric(gsub("%", "", Positive_tests_per)), na.rm = T)) %>% 
  ungroup() %>% 
  # group_by(Regions) %>% 
  # mutate(N_total = sum(N_group)) %>% 
  # ungroup() %>% 
  # mutate(prop = N_group / N_total) %>% 
  mutate(Regions = factor(Regions, 
                          levels = paste0(1:10)),
         age = factor(age, 
                      levels = c("< 1"  ,
                                 "[1,5)" , 
                                 "[5,18)",
                                 "[18,65)",
                                 ">= 65" ))) %>% 
  arrange(Regions, age)



age_colors <- c(
  "< 1"     = "#02818a",
  "[1,5)"     = "#3690c0",
  "[5,18)"    = "#67a9cf",
  "[18,65)"   = "#a6bddb",
  ">= 65"    = "#d0d1e6" 
  
)

library(ggsci)


age_distribution$age <- factor(age_distribution$age, 
                               levels = rev(levels(factor(age_distribution$age))))

age_distribution %>% 
  group_by(age, Regions) %>% 
  summarise(x = mean(hmpv_tests_average)) %>% 
  dplyr::filter(Regions == 2)

#age_distribution %>% group_by(age) %>% summarise(x = mean(positive_tests_per, na.rm = T))


Fig1B <- ggplot(age_distribution, aes(x = Regions, y = hmpv_tests_average/100, fill = age)) +
  geom_bar(stat = "identity", position = "fill", alpha = 1) +  # stacked bar normalized to 1
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = age_colors)+
  #scale_fill_jco()+
  labs(
    title = "ED Visits with HMPV Tests Performed (2016-2025)",
    x = "HHS Region",
    y = "Proportion",
    fill = "Age group"
  ) +
  guides(fill = guide_legend(reverse = T)) +  # reverse legend order
  theme_minimal() +
theme(
  legend.position = "top",
  plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
  axis.title = element_text(size = 20),
  axis.text.y = element_text(size = 18),
  axis.text.x = element_text(size = 20)
)

Fig1B


 
 
age_distribution$age <- factor(age_distribution$age, 
                               levels = rev(levels(factor(age_distribution$age))))

Fig1C =  ggplot(data = age_distribution, aes(x = age, y = hmpv_tests_average/100, fill = age)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = scales::percent(hmpv_tests_average/100, accuracy = 0.1)),
            vjust = -0.1, size = 5)+
  facet_wrap(~Regions, 
               nrow =  5,
               labeller = labeller(Regions =
                                     c("1" = "Region 1: CT, ME, MA, NH, RI, VT",
                                       "2" = "Region 2: NJ, NY, PR, VI",
                                       "3" = "Region 3: DE, DC, MD, PA, VA, WV",
                                       "4" = "Region 4: AL, FL, GA, KY, MS, NC, SC, TN",
                                       "5" = "Region 5: IL, IN, MI, MN, OH, WI",
                                       "6" = "Region 6: AR, LA, NM, OK, TX",
                                       "7" = "Region 7: IA, KS, MO, NE",
                                       "8" = "Region 8: CO, MT, ND, SD, UT, WY ",
                                       "9" = "Region 9: AZ, CA, HI, NV",
                                       "10" = "Region 10: AK, ID, OR, WA")))+
  scale_y_continuous(labels = scales::percent, 
                     breaks = c(0, 0.2, 0.4),
                     limits = c(0, 0.5)) +
  scale_fill_manual(values = age_colors)+
  labs(
    title = "ED Visits with HMPV Tests Performed (2016-2025)",
    x = "Age group",
    y = "Proportion",
    fill = "Age group"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 20)
  )   + 
  theme(strip.text = element_text(size = 12, margin = margin(), face = "bold"),
          strip.background = element_rect(fill = "white", color = "transparent"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8))  

Fig1C
 


age_distribution2 <- epic_ed_all %>% 
  mutate(Regions = ifelse(geography %in% region1, "1", 
                          ifelse(geography %in% region2, "2", 
                                 ifelse(geography %in% region3, "3", 
                                        ifelse(geography %in% region4, "4",
                                               ifelse(geography %in% region5, "5", 
                                                      ifelse(geography %in% region6, "6", 
                                                             ifelse(geography %in% region7, "7", 
                                                                    ifelse(geography %in% region8, "8",
                                                                           ifelse(geography %in% region9, "9", "10")))))))))) %>% 
  dplyr::filter(!is.na(Regions)) %>% 
  mutate(age = ifelse(age == "Less than 1", "< 1", 
                      ifelse(age == "≥ 1 and < 5", "[1,5)", 
                             ifelse(age == "≥ 5 and < 18", "[5,18)", 
                                    ifelse(age == "65 or more", ">= 65", "[18,65)")))))  %>% 
  group_by(Regions, age, week_start) %>% 
  summarise(hmpv_tests_average = mean(as.numeric(gsub("%", "", hmpv_tests)), na.rm = T),
            positive_tests_per = mean(as.numeric(gsub("%", "", Positive_tests_per)), na.rm = T)) %>% 
  ungroup() %>% 
  # group_by(Regions) %>% 
  # mutate(N_total = sum(N_group)) %>% 
  # ungroup() %>% 
  # mutate(prop = N_group / N_total) %>% 
  mutate(Regions = factor(Regions, 
                          levels = paste0(1:10)),
         age = factor(age, 
                      levels = c("< 1"  ,
                                 "[1,5)" , 
                                 "[5,18)",
                                 "[18,65)",
                                 ">= 65" ))) %>% 
  arrange(Regions, age)

age_distribution2$age <- factor(age_distribution2$age, 
                                levels = rev(levels(factor(age_distribution2$age))))


 ggplot(age_distribution2) +
  geom_line(aes(x = week_start, y = hmpv_tests_average/100, color = age)) +  # stacked bar normalized to 1
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values = age_colors)+
  facet_wrap(~Regions)+
  labs(
    title = "HMPV Testing (2020-2025)",
    x = "Date",
    y = "% of ED encounters with a test performed for hMPV", 
    color = "Age group"
  ) +
  guides(color = guide_legend(reverse = TRUE)) +  # reverse legend order
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    strip.text = element_text(size = 15, face = "bold")
  )   

 

# Fig1D <- ggplot(age_distribution, aes(x = Regions, y = positive_tests_per, fill = age)) +
#   geom_bar(stat = "identity", position = "fill", alpha = 1) +  # stacked bar normalized to 1
#   scale_y_continuous(labels = scales::percent) +
#   scale_fill_manual(values = age_colors)+
#   #scale_fill_jco()+
#   labs(
#     title = "Age Distribution of HMPV Testing (2020-2025)",
#     x = "HHS Region",
#     y = "Proportion",
#     fill = "Age group"
#   ) +
#   guides(fill = guide_legend(reverse = TRUE)) +  # reverse legend order
#   theme_minimal() +
#   theme(
#     legend.position = "top",
#     plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
#     axis.title = element_text(size = 15),
#     axis.text = element_text(size = 15),
#     legend.title = element_text(size = 15),
#     legend.text = element_text(size = 15)
#   )   
# 
# Fig1D
 
  