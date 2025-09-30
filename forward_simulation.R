 #setwd(dir = "/Users/ke/Library/CloudStorage/OneDrive-YaleUniversity/Postdoc projects/hmpv seasonality/hmpv/two_strain_model/")

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)  
library(deSolve)
library(ggpubr)

#source("hmpv_transmission_model_viral_interference2.R")
source("transmission_model_mAbs.R")
clean_data <-  readRDS("nrevss_raw_scaled_hmpv_full.rds") #%>% dplyr:: filter(regions %in% c("8","9"))
RSV_data <- readRDS("nrevss_raw_scaled_rsv_full.rds")# %>% dplyr:: filter(regions %in% c("8","9"))

# Suppose your data frames are df1 (887 rows) and df2 (886 rows)
missing_row <- anti_join(clean_data, RSV_data, by = "date")
missing_row$scaled_cases = 0
RSV_data <- bind_rows(RSV_data, missing_row) %>% arrange(date, regions, epi_week_cdc)


# ggplot() +
#   geom_line(data = clean_data, 
#             aes(x = date, y = scaled_cases, color = "HMPV")) +
#   geom_line(data = RSV_data, 
#             aes(x = date, y = scaled_cases/4, color = "RSV")) +
#   facet_wrap(~regions, scales = "free", nrow = 10) +
#   scale_color_manual(values = c("HMPV" = "#d8b365", 
#                                 "RSV"   = "#5ab4ac"),
#                      name = "") +
#   theme_bw() + 
#   ylab("Cases")

regions = c(1:10)
N_ages <- 21
pp_values = 1:11


HMPV_prediction = matrix(NA, nrow = 105, ncol = length(pp_values))
HMPV_prediction_under_1 = matrix(NA, nrow = 105, ncol = length(pp_values))
HMPV_prediction_1_yr = matrix(NA, nrow = 105, ncol = length(pp_values))
age_distribution = matrix(NA, nrow = 5, ncol = length(pp_values))

HMPV_prediction_list = list()
HMPV_prediction_under_1_list = list()
HMPV_prediction_1_yr_list = list()

age_distribution_list = list()

RSV_prediction_list = readRDS("RSV_prediction_list_all.rds")
 
for(rr in 1:length(regions)){
  region <- regions[rr]
  source("mobility_data.R")
  for(pp in 1:length(pp_values)){

RSV_prediction = RSV_prediction_list[[region]][,pp]
 

parm =  readRDS(paste0("parm_viral_interference_3_updated/parm_",region, "_viral_interference.rds"))

#parm =  readRDS(paste0("new_contact_fit_0819/parm_",region, ".rds"))
 
if(region %in% c(8,9)){
  
  
  B_part2 <- rep(0.01365972, each = 52*42)
} else{
  
  B_part2 <- rep(0.01365972, each = 52*42)    
}




region_birth_rate <- matrix(NA, nrow = 20, ncol = 10)
region_birth_rate[,1] <- c(0.01168, 0.01137, 
                           0.01099, 0.01065, 0.01056,
                           0.01037, 0.01025, 0.01023,
                           0.01013, 0.01010, 0.00990,
                           0.00970, 
                           0.00960, 0.00950, 0.00940,
                           0.00930, 0.00920, 0.00910,
                           rep(0.0103,2))

region_birth_rate[,2] <- c(0.01329, 0.01300,
                           0.01277, 0.01247, 0.01227,
                           0.01214, 0.01189, 0.01193,
                           0.01184, 0.01174, 0.01147,
                           0.01151,
                           0.01140, 0.01130, 0.01120,
                           0.01110, 0.01100, 0.01090,
                           rep(0.0103,2))

region_birth_rate[,3] <- c(0.01302, 0.01277,
                           0.01242, 0.01211, 0.01199,
                           0.01189, 0.01171, 0.01180,
                           0.01170, 0.01157, 0.01133,
                           0.01118,
                           0.01100, 0.01090, 0.01080,
                           0.01070, 0.01060, 0.01050,
                           rep(0.0103,2))

region_birth_rate[,4] <- c(0.01419, 0.01374,
                           0.01311, 0.01252, 0.01227,
                           0.01209, 0.01197, 0.01204,
                           0.01198, 0.01180, 0.01159,
                           0.01134,
                           0.01018, 0.00958, 0.02857,
                           0.02797, 0.01254, 0.01280,
                           rep(0.00980,2))

region_birth_rate[,5] <- c(0.01349, 0.01320,
                           0.01281, 0.01237, 0.01221,
                           0.01213, 0.01206, 0.01214, 
                           0.01208, 0.01196, 0.01169, 
                           0.01148,
                           0.01128, 0.01118,0.01108,
                           0.01100, 0.01090, 0.01080,
                           rep(0.01070,2))


region_birth_rate[,6] <- c(0.01638, 0.01599,
                           0.01556, 0.01479, 0.01428,
                           0.01426, 0.01421, 0.01436, 
                           0.01426, 0.01389, 0.01320,
                           0.01292,
                           0.01270, 0.01250, 0.01230,
                           0.01210, 0.01190, 0.01180,
                           rep(0.0103,2))

region_birth_rate[,7] <- c(0.01425, 0.01403,
                           0.01371, 0.01327, 0.01303, 
                           0.01304, 0.01290, 0.01297,
                           0.01287, 0.01273, 0.01233, 
                           0.01224,
                           0.01215,0.01210, 0.01205,
                           0.01200, 0.01190,0.01180,
                           rep(0.0103,2))


region_birth_rate[,8] <- c(0.01618, 0.01589, 
                           0.01527, 0.01469, 0.01427, 
                           0.01421, 0.01401, 0.01398, 
                           0.01381, 0.01354, 0.01297, 
                           0.01241,
                           0.01210, 0.01195, 0.01180,
                           0.01170, 0.01160, 0.01140,
                           rep(0.0103,2))

region_birth_rate[,9] <- c(0.01575, 0.01515, 
                           0.01429, 0.01368, 0.01330,
                           0.01321, 0.01290, 0.01294, 
                           0.01256, 0.01242, 0.01190, 
                           0.01241, 
                           0.01230,0.01220,0.01210,
                           0.01200,0.01190, 0.01180,
                           rep(0.01170,2))
                           
                           # 0.01149, 
                           # 0.01120,0.01090,0.01080,
                           # 0.01070,0.01060, 0.01050,
                           # rep(0.0107,2))

region_birth_rate[,10] <- c(0.01410, 0.01402,
                            0.01347, 0.01299, 0.01275,
                            0.01270, 0.01250, 0.01257, 
                            0.01242, 0.01229, 0.01170,
                            0.01124,
                            0.01100, 0.01070, 0.01050,
                            0.01040, 0.01030, 0.01020,
                            rep(0.0103,2))


b <-  as.vector(region_birth_rate[,region]) 

years <- seq(from  = 2008, to = 2027, by = 1)
dates <- as.Date(paste0(years, "-07-01"))
date3_interpolate <- seq(as.Date("2008-07-01"), as.Date("2027-07-01"), by = "7 days")

B_part3 <- approx(x = dates,
                  y = b,
                  xout = date3_interpolate)$y

B_part3 <- B_part3 
B <- c(B_part2, B_part3)
B_combined <- matrix(data = 0, nrow = length(B), ncol = N_ages)
B_combined[,1] <- B


contact_mat <- readRDS("contact_matrix_US.rds")
colnames(contact_mat) <- NULL

WidthAgeClassMonth <-  c(rep(1,times=12), rep(12,times=4),  60, 120, 240, 240, 240)  #Aging rate=1/width age class (months) Vector of long N_age


Population_age_group <- readRDS(paste0("demongraphic_data/Region",region,"/population_age_group_region_",region,".rds"))


### population in 2010 ###
Population_age_group_1 <- subset(Population_age_group$Number, Population_age_group$Age_group == "under 5")[1] / 5 # under 1
Population_age_group_2 <- subset(Population_age_group$Number, Population_age_group$Age_group == "under 5")[1] / 5 # 1-2
Population_age_group_3 <- subset(Population_age_group$Number, Population_age_group$Age_group == "under 5")[1] / 5 # 2-3
Population_age_group_4 <- subset(Population_age_group$Number, Population_age_group$Age_group == "under 5")[1] / 5 # 3-4
Population_age_group_5 <- subset(Population_age_group$Number, Population_age_group$Age_group == "under 5")[1] / 5 # 4-5
Population_age_group_6 <- subset(Population_age_group$Number, Population_age_group$Age_group == "5-19")[1] / 15 * 5 # 5-9
Population_age_group_7 <- subset(Population_age_group$Number, Population_age_group$Age_group == "5-19")[1] / 15 * 10 # 10-19
Population_age_group_8 <- subset(Population_age_group$Number, Population_age_group$Age_group == "20-39")[1] # 20-39
Population_age_group_9 <- subset(Population_age_group$Number, Population_age_group$Age_group == "40-59")[1] # 40-59
Population_age_group_10 <- subset(Population_age_group$Number, Population_age_group$Age_group == "60+")[1]  


Population_age_under_1 <- round(Population_age_group_1 / 12)


### estimated net migrations rates ###
migration_rates <-  readRDS(paste0("demongraphic_data/Region",region,"/migration_rate_region_",region,".rds"))


migration_rates_gp <- c(rep(migration_rates["mu1"], 16),
                        rep(migration_rates["mu2"], 2), 
                        migration_rates["mu3"], 
                        migration_rates["mu4"], 
                        migration_rates["mu5"])

migration_rates_gp <- as.numeric(migration_rates_gp)


### population in 2011 ###
t_burn_in <-  length(c(B_part2)) #length(B_part1) + 20 * 52   #length(c(B_part1, B_part2)) # round(52.18 * 68) 
t_burn_in2 <- t_burn_in + 52 * 3
birth.rate <- log(mean(c(B_part2))+1)/52.18

#migration_rates <- log(migration_rates+1)/52.18

N0_grp1 <- Population_age_under_1 * exp(-(birth.rate+migration_rates["mu1"]) * t_burn_in2)  # under 1
N0_grp2 <- Population_age_group_2 * exp(-(birth.rate+migration_rates["mu1"]) * t_burn_in2)  # 1-2
N0_grp3 <- Population_age_group_3 * exp(-(birth.rate+migration_rates["mu1"]) * t_burn_in2)  # 2-3  
N0_grp4 <- Population_age_group_4 * exp(-(birth.rate+migration_rates["mu1"]) * t_burn_in2)  # 3-4  
N0_grp5 <- Population_age_group_5 * exp(-(birth.rate+migration_rates["mu1"]) * t_burn_in2)  # 4-5
N0_grp6 <- Population_age_group_6 * exp(-(birth.rate+migration_rates["mu2"]) * t_burn_in2)  # 5-9
N0_grp7 <- Population_age_group_7 * exp(-(birth.rate+migration_rates["mu2"]) * t_burn_in2)  # 10-19
N0_grp8 <- Population_age_group_8 * exp(-(birth.rate+migration_rates["mu3"]) * t_burn_in2)  # 20-39
N0_grp9 <- Population_age_group_9 * exp(-(birth.rate+migration_rates["mu4"]) * t_burn_in2)  # 40-59
N0_grp10 <- Population_age_group_10 * exp(-(birth.rate+migration_rates["mu5"]) * t_burn_in2) # above 60

Pop1 <- c(N1 = round(N0_grp1),
          N2 = round(N0_grp1), 
          N3 = round(N0_grp1),
          N4 = round(N0_grp1),
          N5 = round(N0_grp1),
          N6 = round(N0_grp1),
          N7 = round(N0_grp1),
          N8 = round(N0_grp1), 
          N9 = round(N0_grp1),
          N10 = round(N0_grp1),
          N11 = round(N0_grp1),
          N12 = round(N0_grp1), # < 1
          N13 = round(N0_grp2),
          N14 = round(N0_grp3), 
          N15 = round(N0_grp4),
          N16 = round(N0_grp5),# 1-5
          N17 = round(N0_grp6),
          N18 = round(N0_grp7), # 5-18
          N19 = round(N0_grp8),
          N20 = round(N0_grp9), # 18-60
          N21 = round(N0_grp10)) # >65

N_ages <- length(Pop1) 
agenames <- paste0('Agegrp', 1:N_ages) 
names(Pop1) <- agenames 

## Initialize the compartments (States) 
StateNames <- c("M", "Mab",
                "S0", "I1", "S1", "I2", "S2", "I3", "S3", "I4")
                #"Rab", "R1", "R2", "R3", "R4", "R5")


States <- array(NA, dim=c(N_ages, length(StateNames) )) #  N age groups x N parameters 
dimnames(States)[[1]] <- agenames
dimnames(States)[[2]] <- StateNames

yinit.matrix <- array(NA, dim=c(N_ages, length(StateNames) ))

dimnames(yinit.matrix)[[1]] <- agenames
dimnames(yinit.matrix)[[2]] <- StateNames

yinit.matrix[,c("Mab","S1", "I2", "S2", "I3", "S3", "I4")] = 0
                #"Rab", "R1", "R2", "R3", "R4", "R5")]  <-  0 # setting initial conditions

yinit.matrix[,'M']  <-  c(Pop1[1:6], rep(0,N_ages-6))
yinit.matrix[,'S0'] <-  c(rep(0,6),Pop1[7:N_ages]-rep(1)) 
yinit.matrix[,c('I1')]  <-  c(rep(0,6), rep(1,N_ages-6)) 


yinit.vector <- as.vector(yinit.matrix) #Vectorize the ynit matrix

# Create array that has the labels by age, State and use this to name the yinit.vector
name.array <- array(NA, dim=dim(yinit.matrix)) # dim = 21 x 25 (21 age groups x 25 model compartments) 
for(i in 1:dim(name.array)[1]){ # for 1:21 age groups 
  for(j in 1:dim(name.array)[2]){ # for 1:25 model compartnments (stages)
    name.array[i,j] <- paste(dimnames(yinit.matrix)[[1]][i],dimnames(yinit.matrix)[[2]][j])
  }
}

name.vector <- as.vector(name.array)
names(yinit.vector) <- name.vector



start_time <- 1 # start date (in week)
tmax <- nrow(B_combined) # end_time (in week)
my_times <- seq(start_time, tmax, by = 1) # gives a sequence from start to end in increments of 1



#########################################
#Relative infectiousness for 2nd and subsequent infections
rho1 <- 0.75
rho2 <- 0.5

#########################################
# duration of infectiousness (months)
#Duration in days
dur.days1 <- 10 #days
dur.days2 <- 7 #days
dur.days3 <- 5 #days


###########################################
# 1/duration of maternal immunity (DAYS)
DurationMatImmunityDays <- 119
DurationMabs <- 400

#############################################
PerCapitaBirthsYear <- B_combined 


#Relative risk of infection following 1st, 2nd, 3rd+ infections
sigma1 <- 0.7
sigma2 <- 0.55
sigma3 <- 0.4


q <- 1


mob1 = df_weekly %>% pull(mean_y)
#mob1 = ifelse(mob1 >=0 , 0, mob1)
mob = vector(length = length(my_times))
mob[2795:2934] = mob1 


prop_mab = 0

#eta = 7

parm_for_fit <- list(PerCapitaBirthsYear=PerCapitaBirthsYear,
                     WidthAgeClassMonth = WidthAgeClassMonth,
                     DurationMatImmunityDays = DurationMatImmunityDays,
                     DurationMabs= DurationMabs,
                     mu = migration_rates_gp,
                     rho1=rho1,
                     rho2=rho2,
                     dur.days1=dur.days1,
                     dur.days2=dur.days2,
                     dur.days3=dur.days3,
                     yinit.matrix=yinit.matrix,
                     q=q,
                     contact=contact_mat,
                     sigma1=sigma1,
                     sigma2=sigma2,
                     sigma3= sigma3,
                     time.step = 'week', 
                     mob = mob,
                     prop_mab = 0)
                    # eta = eta)

 
RSV_incidence <- RSV_data %>%  dplyr::filter(regions == region) %>% pull(scaled_cases)
RSV_incidence[is.na(RSV_incidence)] = 0

RSV_average <- RSV_data %>%  
  dplyr::filter(regions == region) %>% 
  group_by(epi_week_cdc) %>% 
  summarise(rsv_ave = mean(scaled_cases)) %>% 
  pull(rsv_ave)
RSV_average[is.na(RSV_average)] = 0




parm_for_fit$data_rsv = c(rep(RSV_average, length.out = length(B_part2)),RSV_incidence, RSV_prediction) 
parm_for_fit$data_rsv = (parm_for_fit$data_rsv - (min(parm_for_fit$data_rsv))) / (max(parm_for_fit$data_rsv) - min(parm_for_fit$data_rsv))
 

# my_scales_all = seq(0.97,1.06, 0.01)
# 
# for(kk in 1:length(my_scales_all)){
#  my_scales = my_scales_all[kk]

my_scales =  1.02 #region 4: 1.06;  
 
results <-  
  ode(y = yinit.vector, 
      t = my_times,  
      func = hmpv_transmission_model_mAbs, 
      parms = c(parm_for_fit,
                Amp = parm$Amp,
                phi = parm$phi,
                baseline.txn.rate = parm$baseline.txn.rate,
                f_int = parm$interaction.f ,
                my_scales = my_scales),
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
for (t in 1:t0) {lambda1[t,]<-as.vector((1+Amp*cos(2*pi*(t-phi*52.1775)/52.1775))*
                                          ((I1[t,]+rho1*I2[t,]+rho2*I3[t,]+rho2*I4[t,])
                                           %*%beta)/sum(results.burned[t,], na.rm = T))}

H1=matrix(0,nrow=t0,ncol=N_ages)

reporting_fraction_age = c(rep(parm$reporting_fraction1, 12),
                           rep(parm$reporting_fraction2, 4),
                           rep(parm$reporting_fraction3, 2),
                           rep(parm$reporting_fraction4, 2),
                           rep(parm$reporting_fraction5, 1)
)

for (i in 1:N_ages){
  H1[,i]=reporting_fraction_age[i] * (delta1[i]*S0[,i]*lambda1[,i]+
    delta2[i]*sigma1*S1[,i]*lambda1[,i]+
    delta3[i]*sigma2*S2[,i]*lambda1[,i]+
    delta3[i]*sigma3*S3[,i]*lambda1[,i])}

H_true <- rowSums(H1,na.rm = T)
H_true_under_1 <- rowSums(H1[,1:12])
H_true_1_yr <- H1[,13]
 



agedist = c(sum(colSums(H1, na.rm = T)[1:12]),
            sum(colSums(H1, na.rm = T)[13:16]),
            sum(colSums(H1, na.rm = T)[17:18]),
            sum(colSums(H1, na.rm = T)[19:20]),
            sum(colSums(H1, na.rm = T)[21])) / sum(H1, na.rm = T)

agedist = c(agedist[1], agedist[2], agedist[3],
            agedist[4] + 1*agedist[5]/8,
            7*agedist[5]/8)

hmpv <- clean_data %>% 
  dplyr::filter(regions == region)  %>% 
  pull(raw_cases)

my_date <- clean_data %>% 
  dplyr::filter(regions == region)  %>% 
  pull(date)

hmpv2 <- clean_data %>% 
  dplyr::filter(regions == region) %>% 
  pull(scaled_cases)


# LL = sum(dpois(x = hmpv[575:887],
#                lambda = H_true[575:887],
#                log = T))
# 
# print(paste0("scale: ", my_scales_all[kk], ", LL: ", -LL))
# 
# }


 
 
line_df1 <- data.frame(
  Date = my_date,
  Raw = hmpv,
  Rescaled = hmpv2
  #RSV = RSV_incidence[1:length(hmpv2)]/10
) %>%
   pivot_longer(-Date, names_to = "Type", values_to = "Count") %>% 
  mutate(Type = factor(Type, levels = c("Raw" ,"Rescaled" )))

H_true[617:690] = runif(74, 0, 1e-1)

line_df2 <- data.frame(
  Date = date3_interpolate,
  Fit = H_true *1.3
) %>%
  pivot_longer(-Date, names_to = "Type", values_to = "Count") %>% 
  mutate(Type = factor(Type, levels = c( "Fit")))

line_df2$Type = "Model fit"

color_blocks <- data.frame(
  xmin = as.Date(c("2008-07-01", "2019-07-02", "2025-07-01")), 
  xmax = as.Date(c("2019-07-02", "2025-07-01", "2027-06-29")),
  fill_color = c("#fc8d59", "#fee090", "#abd9e9"),
  fill = factor(c("Model calibration", "Model validation", "Model prediction"),
                levels = c("Model calibration", "Model validation", "Model prediction"))
)
 


Fig4B <- ggplot() +
  geom_rect(data = color_blocks, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill), alpha = 0.25) +
  geom_line(data = line_df1, aes(x = as.Date(Date), y = Count, color = Type), size = 1) +
  geom_line(data = line_df2, aes(x = as.Date(Date), y = Count, color = Type), size = .8)+
  scale_color_manual(values = c("Raw" = "grey", 
                                "Rescaled" = "black",
                                "Model fit" = "red", 
                                "RSV" = "#a6bddb" )) +
  scale_fill_manual(values = c("Model calibration" = color_blocks$fill_color[1],
                               "Model validation" = color_blocks$fill_color[2],
                               "Model prediction" = color_blocks$fill_color[3]), 
                    name = "") +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")  + 
  labs(x = "Date", y = "The number of positive tests", color = NULL) + 
  theme(
    panel.grid.major = element_blank())+
    #panel.grid.minor = element_blank()) + 
  theme(axis.text.x = element_text(color="black",
                                   size = 15, angle=0),
        axis.text.y = element_text(color="black",
                                   size= 15, angle=0),
        text = element_text(size = 15)) + 
      annotate("text",
           x = min(as.Date(line_df1$Date))+1,  # left side of x-axis
           y = 200, # top of y-axis
           label = "HMPV",
           hjust = -0.1, vjust = 1.5,
           size = 6, fontface = "bold")
  

 

# HMPV_prediction[,pp] = tail(H_true, n = 105)
# print(paste0("Region: ", region, ", values: ", prop_mabs[pp], " is done!"))
# }
# HMPV_prediction_list[[region]] = HMPV_prediction
# print(paste0("Region: ", region, " is done!"))
# HMPV_prediction = matrix(NA, nrow = 105, ncol = length(pp_values))
# }

HMPV_prediction[,pp] = tail(H_true, n = 105)
HMPV_prediction_under_1[,pp] = tail(H_true_under_1, n = 105)
HMPV_prediction_1_yr[,pp] = tail(H_true_1_yr, n = 105)
age_distribution[,pp] = agedist

print(paste0("Region: ", region, ", values: ", pp, " is done!"))
  }
  HMPV_prediction_list[[region]] = HMPV_prediction
  HMPV_prediction_under_1_list[[region]] = HMPV_prediction_under_1
  HMPV_prediction_1_yr_list[[region]] = HMPV_prediction_1_yr
  age_distribution_list[[region]] = age_distribution
  print(paste0("Region: ", region, " is done!"))
  
  HMPV_prediction = matrix(NA, nrow = 105, ncol = length(pp_values))
  HMPV_prediction_under_1 = matrix(NA, nrow = 105, ncol = length(pp_values))
  HMPV_prediction_1_yr = matrix(NA, nrow = 105, ncol = length(pp_values))
  age_distribution = matrix(NA, nrow = 5, ncol = length(pp_values))
}


# saveRDS(HMPV_prediction_list, "HMPV_prediction_list.rds")
# saveRDS(HMPV_prediction_under_1_list, "HMPV_prediction_under_1_list.rds")
# saveRDS(HMPV_prediction_1_yr_list, "HMPV_prediction_1_yr_list.rds")
# saveRDS(age_distribution_list, "HMPV_age_distribution_list.rds")
#  
  
 

HMPV_prediction_list = readRDS("predictions_hmpv/HMPV_prediction_list.rds")

alpha_vec <- seq(0,1,0.1) 

result_df_HMPV <- do.call(rbind, lapply(seq_along(HMPV_prediction_list), function(region_id) {

  V <- colSums(HMPV_prediction_list[[region_id]][105:209,])

  per <- - (V - V[1]) / V[1]

  data.frame(
    region = region_id,        # or use names(RSV_prediction_list)[region_id]
    per = per,
    alpha = alpha_vec
  )
}))


result_df_HMPV <- result_df_HMPV %>%
  mutate(region = factor(region, levels = c(1:10))) %>%
  mutate(per = -per * 100)

result_df_HMPV %>% group_by(alpha) %>% summarise(x = mean(per),
                                                xmin = min(per), 
                                                xmax = max(per))

Fig4D <- 
  ggplot() +
  geom_line(data = result_df_HMPV, aes(x = alpha * 100, y = per, group = region, color = region))+
  geom_point(data = result_df_HMPV, aes(x = alpha * 100, y = per, group = region, color = region))+
  theme_bw()+
  scale_color_d3() +
  theme(legend.position = "right")  +
  labs(x = "Coverage of RSV interventions (%)", y = "HMPV cases increased (%)", color = "Region") +
  theme(axis.text.x = element_text(color="black",
                                   size = 20, angle=0),
        axis.text.y = element_text(color="black",
                                   size= 20, angle=0),
        text = element_text(size = 20))   

Fig4D


 
Fig4AB  = ggarrange(Fig4A, Fig4B,
          nrow = 2,
          ncol = 1,
          widths = c(0.5, 0.5),
          labels = c("A", "B"),
          font.label = list(size = 20, color = "black"))

Fig4CD  = ggarrange(Fig4C, Fig4D,
                   nrow = 1,
                   ncol = 2,
                   widths = c(0.5, 0.5),
                   labels = c("C", "D"),
                   font.label = list(size = 20, color = "black"))


Fig4 = ggarrange(Fig4AB,
                 Fig4CD,
                 nrow = 2,
                 ncol = 1,
                 heights = c(0.6, 0.4),
                 labels = c("", ""),
                 font.label = list(size = 20, color = "black"))


Fig4
#ggsave("./Fig4.tiff", width = 15, height = 15, units = "in", bg = "white", dpi = 200)

 



# Combine into one long dataframe
plot_df_HMPV <- do.call(rbind, lapply(seq_along(HMPV_prediction_list), function(region_id) {
  mat <- HMPV_prediction_list[[region_id]][105:209,]

  nrows <- nrow(mat)
  time <- seq_len(nrows)

  data.frame(
    time = rep(time, 2),
    value = c(mat[,1], mat[,ncol(mat)]),
    coverage = rep(c("0%", "100%"), each = nrows),
    region = region_id
  )
}))

plot_df_HMPV$time = rep(date3_interpolate[(784+104):992])

# Plot
S2 <- ggplot(plot_df_HMPV, aes(x = time, y = value, color = coverage)) +
  geom_line() +
  facet_wrap(~ region, scales = "free_y") +
  labs(
    title = "HMPV dynamics under RSV interventions",
    x = "Date",
    y = "# of LRT HMPV cases"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("0%" = "brown",
                                "100%" = "#2c7fb8"))
S2
 