
library(tidyverse)
library(ggpubr)
library(cowplot)
setwd("/Users/ke/Desktop/OneDrive - Yale University/Postdoc projects/hmpv seasonality/hmpv/")

hmpv = readRDS("hmpv_model2/nrevss_raw_scaled_hmpv2.rds")
RSV = readRDS("nrevss_raw_scaled_RSV.rds")



rsv_hmpv = data.frame(date = hmpv$date, 
                      regions = hmpv$regions,
                      hmpv = hmpv$scaled_cases, 
                      rsv = RSV$scaled_cases)
S4A= rsv_hmpv %>% 
  ggplot() + 
  geom_line(aes(x = date, y = hmpv * 3, color = "hmpv")) + 
  geom_line(aes(x = date, y = rsv, color = "rsv")) + 
  scale_color_manual(values = c("rsv" = "#1b9e77", "hmpv" = "#7570b3"), name = "Virus") +
  ylab("# of cases") + 
  facet_wrap(~ regions, scales = "free_y", nrow = 4,
             labeller = labeller(regions =
                                   c("1" = "Region 1: CT, ME, MA, NH, RI, VT",
                                     "2" = "Region 2: NJ, NY, PR, VI",
                                     "3" = "Region 3: DE, DC, MD, PA, VA, WV",
                                     "4" = "Region 4: AL, FL, GA, KY, MS, NC, SC, TN",
                                     "5" = "Region 5: IL, IN, MI, MN, OH, WI",
                                     "6" = "Region 6: AR, LA, NM, OK, TX",
                                     "7" = "Region 7: IA, KS, MO, NE",
                                     "8" = "Region 8: CO, MT, ND, SD, UT, WY ",
                                     "9" = "Region 9: AZ, CA, HI, NV",
                                     "10" = "Region 10: AK, ID, OR, WA"))) + 
  theme_bw() 
 

hmpv <- readRDS("nrevss_raw_scaled_hmpv_full.rds")

 

Fig1A  <- hmpv %>%
  filter(regions %in% c(1,4,8,9)) %>% 
  dplyr:: filter(date >= "2008-07-01") %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = raw_cases, group = regions, color = "Raw Positive Tests"), linewidth = 1) +
  geom_line(aes(y = scaled_cases, group = regions, color = "Scaled Positive Tests"), linewidth = 1) +
  geom_line(aes(y = MA_52w.hmpv / 20, group = regions, color = "52-Week Moving Average of Total Tests"), linewidth = 1, alpha = .7) +
  facet_wrap(~regions, scales = "free",
             nrow =  2,
             labeller = labeller(regions =
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
  scale_color_manual(
    name = "",
    values = c(
      "Raw Positive Tests" = "#bdbdbd",
      "Scaled Positive Tests" = "#252525",
      "52-Week Moving Average of Total Tests" = "orange"
    )
  ) +
    theme_half_open() + 
    theme(strip.text = element_text(size = 18, margin = margin(), face = "bold"),
          strip.background = element_rect(fill = "white", color = "transparent"))  +
    theme(legend.key = element_rect(fill = "white", colour = "white"), legend.position = "top",legend.justification = "center") +
    
    theme(panel.grid.minor  = element_blank()) +
    labs(title = "",
         x = "Date",
         y = "Number of positive tests") +
    theme(axis.text.x = element_text(color="black",
                                     size=20, angle=0,vjust = 0.5),
          axis.text.y = element_text(color="black",
                                     size=20, angle=0),
          text = element_text(size = 20))



Fig1BC <- ggarrange(Fig1B, Fig1C,
                    nrow = 1, 
                    ncol = 2,
                    labels = c("B", "C"),  
                    widths = c(0.4, 0.6),
                    #heights = c(0.6, 0.4),
                    font.label = list(size = 20, color = "black"))
Fig1BC


 

Fig1 <- ggarrange(Fig1A, Fig1BC,
                    nrow = 2, 
                    ncol = 1,
                    labels = c("A", " "),  
                    #widths = c(0.7, 0.3),
                    heights = c(0.4, 0.6),
                    font.label = list(size = 20, color = "black"))


Fig1

ggsave("Figures/main_figures/Fig1.tiff", plot = Fig1, device = "tiff",
       width = 20, height = 12, units = "in", dpi = 300,
       compression = "lzw")  # lossless compression
 

