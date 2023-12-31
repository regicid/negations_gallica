###############################################
### PACKAGES, INITIALISATIONS, DATA IMPORTS ###
###############################################

# Packages
library(httr)
library(jsonlite)
library(dplyr)
library(ggplot2)
library(zoo)
library(lubridate)
library(gridExtra)
library(DHARMa)
library(tidyr)
library(zoo)
library(plotly)

# Custom path
path_projet = "~/code/negations_gallica/"
# Analysis parameters
time_frame = c("1850", "1920")
liste_corpus = c("huma",
                 "figaro",
                 "moniteur",
                 "temps",
                 "petit_journal",
                 "journal_des_debats",
                 "la_presse",
                 "petit_parisien")
scrap_data = T # If the data are already downloaded on your machine set to FALSE

# Custom functions
# API retreiving function call
data_call <- function(corpus = "presse", mot, debut, fin, resolution) {
  call = paste0("https://shiny.ens-paris-saclay.fr/guni/query?corpus=", corpus,
                "&mot=", mot,
                "&from=", debut,
                "&to=", fin,
                "&resolution=", resolution)
  call = GET(url = call, timeout(100000))
  df = content(call, "text", encoding = "UTF-8")
  df <- read.table(text = df, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  return(df)
}
# Rolling weighted mean without outliers
compute_mean_without_outliers <- function(temp, probs, weights = 1) {
  if (weights == 1) {
    weights = rep(1, length(temp))
  } else {
    if (weights == 'time_lin') { # For a weighted rolling mean. Carefull, it only works with centered rolling windows
      window_size = length(temp)
      middle = as.integer(window_size/2)
      if (middle == window_size/2) {pad = NULL} else{pad = 0}
      weights = c(c(1:middle), 
                  c(pad), 
                  c(middle:1))
      weights = weights+1
    }
  }
  q <- quantile(na.omit(temp), probs = probs/2) 
  BOOL1 = temp >= q
  q <- quantile(na.omit(temp), probs = 1-(probs/2))  
  BOOL2 = temp <= q
  temp <- temp[BOOL1 & BOOL2]
  weights <- weights[BOOL1 & BOOL2]
  mean_without_outliers <- weighted.mean(x = temp, 
                                         w = weights, 
                                         na.rm = TRUE)
  return(mean_without_outliers)
}
# Rolling weighted std without outliers
compute_var_without_outliers <- function(temp, probs, weights = 1) {
  if (weights == 1) {
    weights = rep(1, length(temp))
  } else {
    if (weights == 'time_lin') { # For a weighted rolling mean. Carefull, it only works with centered rolling windows
      window_size = length(temp)
      middle = as.integer(window_size/2)
      if (middle == window_size/2) {pad = NULL} else{pad = 0}
      weights = c(c(1:middle), 
                  c(pad), 
                  c(middle:1))
      weights = weights+1
    }
  }
  q <- quantile(na.omit(temp), probs = probs/2) 
  BOOL1 = temp >= q
  q <- quantile(na.omit(temp), probs = 1-(probs/2))  
  BOOL2 = temp <= q
  temp <- temp[BOOL1 & BOOL2]
  weights <- weights[BOOL1 & BOOL2]
  weights_square <- weights**2
  sd_without_outliers <- weighted.mean(x = temp, 
                                       w = weights_square, 
                                       na.rm = TRUE) * sum(weights_square)
  sd_without_outliers = sd_without_outliers / sum(weights)**2
  return(sd_without_outliers)
}

# Gallicagram data import
if (scrap_data) {
  gallicagram_data = 0
  for (corpus in liste_corpus) {
    print(corpus)
    df = data_call(mot = 'pas', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df$journal = corpus
    if (length(gallicagram_data) == 1) {
      gallicagram_data = df
    } else {
      gallicagram_data = gallicagram_data %>% rbind(df)
    }
  }
  gallicagram_data$date <- as.Date(paste(gallicagram_data$annee, 
                                         gallicagram_data$mois, 
                                         gallicagram_data$jour, sep = "-"), 
                                   format = "%Y-%m-%d")
  write.csv(gallicagram_data, 
            paste0(path_projet, "\\gallicagram_data.csv"), 
            fileEncoding = 'UTF-8', row.names = FALSE)
} else {
  gallicagram_data = read.csv(paste0(path_projet, "\\gallicagram_data.csv"), encoding = 'UTF-8')
}

gallicagram_data$ratio = gallicagram_data$n / gallicagram_data$total

# Specific important subsets
BOOL = (gallicagram_data$journal == 'moniteur') & (year(gallicagram_data$date) > 1869)
gallicagram_data = gallicagram_data[!BOOL,]

##########################
### VARIANCE COMPUTING ###
##########################

# Our data points are not comparable to one another
# This is caused by the fact that our measure n/total is a random variable with varying precision
# The bigger the total, the less the variable will be subject to noise and the more efficiently we will able to measure its correlation to specific events
# To understand our method, see paper for details

# We measure empirical variance by creating bins of the total variable and measuring the variance in each of these
reg_data = gallicagram_data %>% 
  mutate(year = year(date)) %>% 
  select(journal, ratio, n, total, year) %>% 
  na.omit() %>% 
  mutate(total_bin_100 = exp(round(log(total*100)/100)))
reg_data_total_bin = reg_data %>%
  group_by(total_bin_100, journal, year) %>%
  filter(ratio <= quantile(ratio, 0.99),
         ratio >= quantile(ratio, 0.01)) %>%
  # slice_sample(n = 30) %>% 
  summarise(var_ratio = var(ratio, na.rm = TRUE),
            ratio = mean(ratio, na.rm = TRUE),
            total = mean(total, na.rm = TRUE),
            group_size = n()) %>%
  ungroup() %>% 
  filter(group_size >= 30,
         journal != 'moniteur', # These newspapers have too much annual fluctuations
         journal != 'presse')

# We validate empirically our model using a train / test split
train_ind <- sample(seq_len(nrow(reg_data_total_bin)), size = as.integer(nrow(reg_data_total_bin)*0.75))
# reg_data_total_bin_train = reg_data_total_bin[train_ind,]
# reg_data_total_bin_test = reg_data_total_bin[-train_ind,]
reg_data_total_bin_train = reg_data_total_bin %>% filter(journal != 'figaro')
reg_data_total_bin_test = reg_data_total_bin%>% filter(journal == 'figaro')

custom_var_computer <- function(ratio, total, b_hat_T, b_hat_R, ess, conf_int = 0.5) {
  T = total
  m_90 = ratio + qnorm(conf_int, mean = 0, sd = 1) * sqrt(ratio/(ess*T))
  V_X = m_90*(1-m_90)
  V_R = (ess * V_X / T) + (ess**2) * (b_hat_T*(1/log(T)) + b_hat_R*V_X)
  return(V_R)
}
cov_ij_computer <- function(var, ess, ratio, total) {
  cov_ij = (1/(ess**2)) * (var - ((ess**2)*ratio/total))
  return(cov_ij)
}

ess = quantile(reg_data$ratio, 0.99)
reg_data_total_bin_train = reg_data_total_bin_train %>% 
  mutate(empirical_cov_ij = cov_ij_computer(var = reg_data_total_bin_train$var_ratio,
                                            ess = ess,
                                            ratio = reg_data_total_bin_train$ratio,
                                            total = reg_data_total_bin_train$total))
plot(reg_data_total_bin_train$total, reg_data_total_bin_train$empirical_cov_ij)
plot(log(reg_data_total_bin_train$total), reg_data_total_bin_train$empirical_cov_ij)
plot(1/reg_data_total_bin_train$total, reg_data_total_bin_train$empirical_cov_ij)
plot(1/log(reg_data_total_bin_train$total), reg_data_total_bin_train$empirical_cov_ij)
plot(reg_data_total_bin_train$ratio, reg_data_total_bin_train$empirical_cov_ij)
reg <- lm(empirical_cov_ij ~ 0 + I(1/log(total)) + ratio,
          data = reg_data_total_bin_train)
summary(reg)
b_hat_T = coef(reg)['I(1/log(total))']
b_hat_R = coef(reg)['ratio']
conf_int = 0.5
# b_hat_T = confint(reg, level = 0.90)['I(1/log(total))', '95 %']
# b_hat_R = confint(reg, level = 0.90)['ratio', '5 %']
# conf_int = 0.9
reg_data_total_bin_test = reg_data_total_bin_test %>% 
  mutate(var_math = custom_var_computer(ratio = ratio,
                                        total = total,
                                        b_hat_T = b_hat_T,
                                        b_hat_R = b_hat_R,
                                        ess = ess,
                                        conf_int = conf_int)) %>% 
  mutate(erreur_pred = (var_ratio - var_math) / var_ratio)
plot(reg_data_total_bin_test$var_ratio, reg_data_total_bin_test$var_math)
abline(a = 0, b = 1, col = 'red')
plot(density(reg_data_total_bin_test$erreur_pred, na.rm = TRUE))
abline(v = mean(reg_data_total_bin_test$erreur_pred, na.rm = TRUE), col = 'red')
reg <- lm(var_ratio ~ 0 + var_math, data = reg_data_total_bin_test)
summary(reg)
# We compare with a classical linear regression model estimation
reg <- lm(log(var_ratio) ~ log(total) + log(ratio),
          data = reg_data_total_bin_train)
summary(reg)
reg_data_total_bin_test = reg_data_total_bin_test %>% 
  mutate(var_math2 = exp(predict(reg, reg_data_total_bin_test))) %>% 
  mutate(erreur_pred2 = (var_ratio - var_math2) / var_ratio)
plot(reg_data_total_bin_test$var_ratio, reg_data_total_bin_test$var_math2)
abline(a = 0, b = 1, col = 'red')
plot(density(reg_data_total_bin_test$erreur_pred2, na.rm = TRUE))
abline(v = mean(reg_data_total_bin_test$erreur_pred2, na.rm = TRUE), col = 'red')
reg <- lm(var_ratio ~ 0 + var_math2, data = reg_data_total_bin_test)
summary(reg)

# We compute estimated stds and we compute associated weights for each of our journal-date observations
# Var math
reg_data_total_bin = reg_data_total_bin %>% 
  mutate(empirical_cov_ij = cov_ij_computer(var = var_ratio,
                                            ess = ess,
                                            ratio = ratio,
                                            total = total))
reg <- lm(empirical_cov_ij ~ 0 + I(1/log(total)) + ratio,
          data = reg_data_total_bin)
summary(reg)
b_hat_T = coef(reg)['I(1/log(total))']
b_hat_R = coef(reg)['ratio']
temp = gallicagram_data %>% 
  select(total, journal, date, ratio) %>% 
  mutate(var_math = custom_var_computer(ratio = ratio,
                                        total = total,
                                        b_hat_T = b_hat_T,
                                        b_hat_R = b_hat_R,
                                        ess = ess,
                                        conf_int = 0.50)) %>% 
# Var ESS
  mutate(var_ess = ratio/(total*max(ratio, na.rm = TRUE)))
plot(temp$var_math, temp$var_ess)
abline(a = 0, b = 1, col = 'red')
# Var log log
reg <- lm(log(var_ratio) ~ log(total) + log(ratio),
          data = reg_data_total_bin)
summary(reg)
temp = temp %>% 
  mutate(var_loglog = exp(predict(reg, temp)),
# Var benoit
         var_benoit = 3.5*ratio/total)
  
# To prevent absurd extreme values
min = quantile(temp$var_math, 0.05)
temp$var_math[temp$var_math < min] = min
gallicagram_data$var_math = NULL
min = quantile(temp$var_ess, 0.05)
temp$var_ess[temp$var_ess < min] = min
gallicagram_data$var_ess = NULL
min = quantile(temp$var_loglog, 0.05)
temp$var_loglog[temp$var_loglog < min] = min
gallicagram_data$var_loglog = NULL
min = quantile(temp$var_benoit, 0.05)
temp$var_benoit[temp$var_benoit < min] = min
gallicagram_data$var_benoit = NULL
gallicagram_data = gallicagram_data %>% 
  left_join(temp %>% select(journal, date, 
                            var_math, var_ess, var_loglog, var_benoit)) %>% 
  mutate(weights_math = 1/sqrt(var_math),
         weights_ess = 1/sqrt(var_ess),
         weights_loglog = 1/sqrt(var_loglog),
         weights_benoit = 1/sqrt(var_benoit),
         weights_avg = (weights_math + weights_ess + weights_loglog + weights_benoit)/4,
         weights = 1/sqrt(var_math))
plot(density(gallicagram_data$weights_math, na.rm = TRUE))
plot(density(gallicagram_data$weights_ess, na.rm = TRUE))
plot(density(gallicagram_data$weights_loglog, na.rm = TRUE))
plot(density(gallicagram_data$weights_benoit, na.rm = TRUE))
plot(density(gallicagram_data$weights_avg, na.rm = TRUE))

ggplot(gallicagram_data) +
  geom_density(aes(x = weights_math), fill = "blue", alpha = 0.5) +
  geom_density(aes(x = weights_ess), fill = "red", alpha = 0.5) +
  geom_density(aes(x = weights_loglog), fill = "green", alpha = 0.5) +
  geom_density(aes(x = weights_benoit), fill = "black", alpha = 0.5)

#############################
### SEASONALITY COMPUTING ###
#############################
#### Trend Computing ####
# Before doing further computations, we need to remove the mean
# Trend level
# The variable used to control for the trend level is a moving average centered on date with a width of several years
# NB : For each trend / seasonality, outliers are removed as the can be considered as unidentified potential events non-related to seasonality 

trend_data = data.frame()
# !! 1 minute * nb of journal to run this loop
for (journal in liste_corpus) {
  
  print(journal) # To keep track
  
  # For loop initialization
  data = gallicagram_data[gallicagram_data$journal == journal,]
  data = data[order(data$date),]
  
  window = as.integer(365.25*3)
  outliers = 0.1
  moving_average <- rollapply(data$ratio, 
                              width = window,
                              align = "center",
                              FUN = compute_mean_without_outliers, 
                              weights = "time_lin",
                              probs = outliers, 
                              fill = NA)
  data$y5_trend = moving_average
  data$ratio_detrend = data$ratio - data$y5_trend
  
  # End of the loop, data storing
  trend_data = trend_data %>% rbind(data)
}

#### Seasonality computing ####

clean_data = data.frame()
# !! 1 minute * nb of journal to run this loop
for (journal in liste_corpus) {
  
  print(journal) # To keep track
  
  # For loop initialisation
  data = trend_data[trend_data$journal == journal,]
  data = data[order(data$date),]
  data = data %>% mutate(date = as.Date(date, format = "%Y-%m-%d"))
  data$ratio_desais = data$ratio_detrend
  
  # Month seasonality  
  # The variable used to control for the month seasonality is the average values of the same month as the considered date on a 5 year window
  window_size = 30
  outliers = 0.20
  moving_average <- rollapply(data$ratio_detrend, # We first compute the moving window of the month average
                              width = window_size, 
                              align = "center", 
                              FUN = compute_mean_without_outliers, 
                              probs = outliers, 
                              weights = 'time_lin',
                              fill = NA)
  temp = data
  temp$s_mois = moving_average
  temp1 = temp %>% mutate(date = date + years(1))
  temp2 = temp %>% mutate(date = date + years(2))
  temp3 = temp %>% mutate(date = date - years(1))
  temp4 = temp %>% mutate(date = date - years(2))
  temp$s_mois = (temp1$s_mois + temp2$s_mois*(1/2) + temp3$s_mois + temp4$s_mois*(1/2)) / (1 + (1/2) + 1 + (1/2)) # We then compute the month seasonality as the average of the month neighbouring years
  data$s_mois = NULL
  data = data %>% left_join(temp %>% select(s_mois, date))
  data = data %>% 
    mutate(ratio_desais = ratio_desais - s_mois) # We now remove the seasonality from the data
  
  # Week seasonality
  # The variable computes the average of the day effects within 3 years
  data$jour_s = weekdays(data$date)
  temp1 = data.frame()
  for (day in unique(data$jour_s)) {
    temp2 = data[data$jour_s == day,]
    window_size = 3*365.25
    outliers = 0.20
    moving_average <- rollapply(temp2$ratio_desais, 
                                width = as.integer(window_size/7), 
                                align = "center", 
                                FUN = compute_mean_without_outliers, 
                                probs = outliers, 
                                weights = 'time_lin',
                                fill = NA)
    temp2$jour_s = moving_average
    temp2 = temp2[c("date", "jour_s")]
    temp1 = temp1 %>% rbind(temp2)
  }
  data = data %>% select(-jour_s) %>% left_join(temp1, by = "date")
  data = data %>% 
    mutate(ratio_desais = ratio_desais - jour_s) # We now remove the seasonality from the data
  
  # End of the loop, data storing
  clean_data = clean_data %>% rbind(data)
}

#### Computing rollmean, sd ####

# !! Takes a few minutes
clean_data = clean_data %>% 
  # Ratio with trend deaseasonalised
  mutate(ratio_desais_trend = ratio_desais + y5_trend) %>% 
  mutate(ratio_desais_trend_90 = qnorm(0.9, 
                                      mean = y5_trend, 
                                      sd = sqrt(var_math)),
         ratio_desais_trend_10 = qnorm(0.1, 
                                      mean = y5_trend, 
                                      sd = sqrt(var_math))) %>% 
  # Ratio with trend smoothed on 30 days
  mutate(ratio_rolled_30 = rollapply(ratio, 
                                     width = 30,
                                     compute_mean_without_outliers,
                                     probs = 0.01, weights = 'time_lin',
                                     align = 'center', fill = NA, na.pad = TRUE),
         var_pred_rolled_30 = rollapply(var_math, # Intermediate stage to compute rolling variance
                                        width = 30,
                                        compute_var_without_outliers,
                                        probs = 0.01, weights = 'time_lin',
                                        align = 'center', fill = NA, na.pad = TRUE),
         ratio_rolled_30_90 = qnorm(0.9, 
                                    mean = y5_trend, 
                                    sd = sqrt(var_pred_rolled_30)),
         ratio_rolled_30_10 = qnorm(0.1, 
                                  mean = y5_trend, 
                                  sd = sqrt(var_pred_rolled_30))) %>% 
  # Ratio with trend deaseasonalised smoothed on 30 days
  mutate(ratio_desais_trend_rolled_30 = rollapply(ratio_desais_trend, 
                                                  width = 30,
                                                  compute_mean_without_outliers,
                                                  probs = 0.01, weights = 'time_lin',
                                                  align = 'center', fill = NA, na.pad = TRUE),
         ratio_desais_trend_rolled_30_90 = qnorm(0.9, 
                                                 mean = y5_trend, 
                                                 sd = sqrt(var_pred_rolled_30)),
         ratio_desais_trend_rolled_30_10 = qnorm(0.1, 
                                                 mean = y5_trend, 
                                                 sd = sqrt(var_pred_rolled_30))) %>% 
  # Ratio with trend deaseasonalised smoothed on 90 days
  mutate(ratio_desais_trend_rolled_90 = rollapply(ratio_desais_trend, 
                                                  width = 90,
                                                  compute_mean_without_outliers,
                                                  probs = 0.01, weights = 'time_lin',
                                                  align = 'center', fill = NA, na.pad = TRUE),
         var_pred_rolled_90 = rollapply(var_math, # Intermediate stage to compute rolling variance
                                        width = 90,
                                        compute_var_without_outliers,
                                        probs = 0.01, weights = 'time_lin',
                                        align = 'center', fill = NA, na.pad = TRUE),
         ratio_desais_trend_rolled_90_90 = qnorm(0.9, 
                                                 mean = y5_trend, 
                                                 sd = sqrt(var_pred_rolled_90)),
         ratio_desais_trend_rolled_90_10 = qnorm(0.1, 
                                                 mean = y5_trend, 
                                                 sd = sqrt(var_pred_rolled_90)))

# Saving data since all this was computationally intense 🥵
write.csv(clean_data, 
          paste0(path_projet, "\\clean_data.csv"), 
          fileEncoding = 'UTF-8', row.names = FALSE)

#### Visualization ####

clean_data = read.csv(paste0(path_projet, "\\clean_data.csv"), encoding = 'UTF-8')
clean_data = clean_data %>% mutate(date = as.Date(date, format = "%Y-%m-%d"))

# Let's plot this new data
plot_list = list()
for (journal in "figaro") {
  BOOL = (clean_data$journal == journal) & (year(clean_data$date) > 1880) & (year(clean_data$date) < 1883)
  data = clean_data[BOOL,]
  new_plot = ggplot(data, aes(x = date)) +
    geom_line(aes(y = ratio_desais_trend_rolled_30), col = 'blue') +
    geom_line(aes(y = ratio_desais_trend_rolled_30_90), col = 'red') +
    geom_line(aes(y = ratio_desais_trend_rolled_30_10), col = 'red')
  new_plot = new_plot +
    labs(title = journal) +
    labs(x = NULL, y = NULL) +
    scale_x_date(
      date_breaks = "1 year",  
      date_labels = "%Y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot_list = c(plot_list, list(new_plot))
}
grid.arrange(grobs = plot_list, ncol = 1)

plot_ly(clean_data,
        x=~date,
        y=~ratio_desais_trend,
        mode="lines")


#############################
### EVENT REGRESSION TEST ###
#############################

# Analysing wether events have an influence on time series
# We will be more precise than in the the above graphic section as we will not resort on smoothing out the data
for (liste_corpus in c("temps","figaro","journal_des_debats","petit_journal","petit_parisien")){
reg_data = filter(clean_data,journal==liste_corpus)

#### Event Creation ####
# We create the events we want to test
# Regime
reg_data$Regime = 'Other'
BOOL = (as.Date("1851-12-14", format = "%Y-%m-%d") <= reg_data$date) & (reg_data$date <= as.Date("1870-10-27", format = "%Y-%m-%d"))
reg_data$Regime[BOOL] = 'EmpII'
BOOL = (as.Date("1870-09-04", format = "%Y-%m-%d") <= reg_data$date) & (reg_data$date <= as.Date("1940-07-10", format = "%Y-%m-%d"))
reg_data$Regime[BOOL] = 'RepIII'
# Franco Prussian War
# Commune and political instability
reg_data$guerre70 = 0 
BOOL = (as.Date("1870-07-19", format = "%Y-%m-%d") <= reg_data$date) & (reg_data$date <= as.Date("1871-05-28", format = "%Y-%m-%d"))
reg_data$guerre70[BOOL] = 1
# WWI
reg_data$WWI = 0 
BOOL = (as.Date("1914-07-28", format = "%Y-%m-%d") <= reg_data$date) & (reg_data$date <= as.Date("1918-11-11", format = "%Y-%m-%d"))
reg_data$WWI[BOOL] = 1
# 14 July
reg_data$July14 = 0 
BOOL = (month(reg_data$date) == 7) & (day(reg_data$date) == 14)
reg_data$July14[BOOL] = 1

#### Event proxies ####
# Some of the events do not have strict timelines and occur with varying intensity
# To take this into account, we create variables that measure this intensity
# The ratio_ variables define the mean occurence of certain words inside the journal_list corpus
# They are then smoothed and set to 0 on irrelevant periods

### Dreyfus ###
if (scrap_data) {
  df_dreyfus = 0
  for (corpus in liste_corpus) {
    print(corpus)
    df = data_call(mot = 'dreyfus', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df$journal = corpus
    if (length(df_dreyfus) == 1) {
      df_dreyfus = df
    } else {
      df_dreyfus = df_dreyfus %>% rbind(df)
    }
  }
  df_dreyfus = df_dreyfus %>% 
    mutate(ratio_dreyfus = n / total,
           date =  as.Date(paste(annee, mois, jour, sep = "-"), format = "%Y-%m-%d")) %>% 
    group_by(date) %>% 
    summarise(ratio_dreyfus = mean(ratio_dreyfus, na.rm = TRUE)) %>% ungroup()
  df_dreyfus$ratio_dreyfus = rollmean(df_dreyfus$ratio_dreyfus, 30, align = "center", na.rm = TRUE, na.pad = TRUE)
  BOOL = (as.Date("1894-09-01", format = "%Y-%m-%d") <= df_dreyfus$date) & (df_dreyfus$date <= as.Date("1902-00-00", format = "%Y-%m-%d"))
  df_dreyfus$ratio_dreyfus[!BOOL] = 0
  
  write.csv(df_dreyfus, 
            paste0(path_projet, "\\df_dreyfus.csv"), 
            fileEncoding = 'UTF-8', row.names = FALSE)
  } else {
    df_dreyfus = read.csv(paste0(path_projet, "\\df_dreyfus.csv"), encoding = 'UTF-8')
    df_dreyfus = df_dreyfus %>% mutate(date = as.Date(date, format = "%Y-%m-%d"))
}
reg_data = reg_data %>% left_join(df_dreyfus %>% select(date, ratio_dreyfus))

### Panama scandal ###
if (scrap_data) {
  df_panama = 0
  for (corpus in liste_corpus) {
    print(corpus)
    df = data_call(mot = 'procès', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df2 = data_call(mot = 'avocat', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df3 = data_call(mot = "l'avocat", debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df4 = data_call(mot = 'audience', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df5 = data_call(mot = "l'audience", debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df$n = df$n + df2$n + df3$n + df4$n + df5$n
    df$journal = corpus
    if (length(df_panama) == 1) {
      df_panama = df
    } else {
      df_panama = df_panama %>% rbind(df)
    }
  }
  df_panama = df_panama %>% 
    mutate(ratio_panama = n / total,
           date =  as.Date(paste(annee, mois, jour, sep = "-"), format = "%Y-%m-%d")) %>% 
    group_by(date) %>% 
    summarise(ratio_panama = mean(ratio_panama, na.rm = TRUE)) %>% ungroup()
  df_panama$ratio_panama = rollmean(df_panama$ratio_panama, 30, align = "center", na.rm = TRUE, na.pad = TRUE)
  #BOOL = (as.Date("1889-02-01", format = "%Y-%m-%d") <= df_panama$date) & (df_panama$date <= as.Date("1898-05-30", format = "%Y-%m-%d"))
  #df_panama$ratio_panama[!BOOL] = 0
  
  write.csv(df_panama, 
            paste0(path_projet, "\\df_panama.csv"), 
            fileEncoding = 'UTF-8', row.names = FALSE)
  } else {
    df_panama = read.csv(paste0(path_projet, "\\df_panama.csv"), encoding = 'UTF-8')
    df_panama = df_panama %>% mutate(date = as.Date(date, format = "%Y-%m-%d"))
}
reg_data = reg_data %>% left_join(df_panama %>% select(date, ratio_panama))

### Strikes ###
if (scrap_data) {
  df_greve = 0
  for (corpus in liste_corpus) {
    print(corpus)
    df = data_call(mot = 'grève', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df2 = data_call(mot = 'grèves', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df$n = df$n + df2$n
    df$journal = corpus
    if (length(df_greve) == 1) {
      df_greve = df
    } else {
      df_greve = df_greve %>% rbind(df)
    }
  }
  df_greve = df_greve %>% 
    mutate(ratio_greve = n / total,
           date =  as.Date(paste(annee, mois, jour, sep = "-"), format = "%Y-%m-%d")) %>% 
    group_by(date) %>% 
    summarise(ratio_greve = mean(ratio_greve, na.rm = TRUE)) %>% ungroup()
  df_greve$ratio_greve = rollmean(df_greve$ratio_greve, 30, align = "center", na.rm = TRUE, na.pad = TRUE)
  
  write.csv(df_greve, 
            paste0(path_projet, "\\df_greve.csv"), 
            fileEncoding = 'UTF-8', row.names = FALSE)
} else {
  df_greve = read.csv(paste0(path_projet, "\\df_greve.csv"), encoding = 'UTF-8')
  df_greve = df_greve %>% mutate(date = as.Date(date, format = "%Y-%m-%d"))
}
reg_data = reg_data %>% left_join(df_greve %>% select(date, ratio_greve))

### Cholera epidemics ###
if (scrap_data) {
  df_cholera = 0
  for (corpus in liste_corpus) {
    print(corpus)
    df = data_call(mot = 'choléra', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df$journal = corpus
    if (length(df_cholera) == 1) {
      df_cholera = df
    } else {
      df_cholera = df_cholera %>% rbind(df)
    }
  }
  df_cholera = df_cholera %>% 
    mutate(ratio_cholera = n / total,
           date =  as.Date(paste(annee, mois, jour, sep = "-"), format = "%Y-%m-%d")) %>% 
    group_by(date) %>% 
    summarise(ratio_cholera = mean(ratio_cholera, na.rm = TRUE)) %>% ungroup()
  df_cholera$ratio_cholera = rollmean(df_cholera$ratio_cholera, 30, align = "center", na.rm = TRUE, na.pad = TRUE)
  
  write.csv(df_cholera, 
            paste0(path_projet, "\\df_cholera.csv"), 
            fileEncoding = 'UTF-8', row.names = FALSE)
} else {
  df_cholera = read.csv(paste0(path_projet, "\\df_cholera.csv"), encoding = 'UTF-8')
  df_cholera = df_cholera %>% mutate(date = as.Date(date, format = "%Y-%m-%d"))
}
reg_data = reg_data %>% left_join(df_cholera %>% select(date, ratio_cholera))

### War ###
if (scrap_data) {
  df_guerre = 0
  for (corpus in liste_corpus) {
    print(corpus)
    df = data_call(mot = 'guerre', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df$journal = corpus
    if (length(df_guerre) == 1) {
      df_guerre = df
    } else {
      df_guerre = df_guerre %>% rbind(df)
    }
  }
  df_guerre = df_guerre %>% 
    mutate(ratio_guerre = n / total,
           date =  as.Date(paste(annee, mois, jour, sep = "-"), format = "%Y-%m-%d")) %>% 
    group_by(date) %>% 
    summarise(ratio_guerre = mean(ratio_guerre, na.rm = TRUE)) %>% ungroup()
  df_guerre$ratio_guerre = rollmean(df_guerre$ratio_guerre, 30, align = "center", na.rm = TRUE, na.pad = TRUE)
  
  write.csv(df_guerre, 
            paste0(path_projet, "\\df_guerre.csv"), 
            fileEncoding = 'UTF-8', row.names = FALSE)
} else {
  df_guerre = read.csv(paste0(path_projet, "\\df_guerre.csv"), encoding = 'UTF-8')
  df_guerre = df_guerre %>% mutate(date = as.Date(date, format = "%Y-%m-%d"))
}
reg_data = reg_data %>% left_join(df_guerre %>% select(date, ratio_guerre))

### 1900 Exposition universelle ###
if (scrap_data) {
  df_expo = 0
  for (corpus in liste_corpus) {
    print(corpus)
    df = data_call(mot = 'exposition', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df$journal = corpus
    if (length(df_expo) == 1) {
      df_expo = df
    } else {
      df_expo = df_expo %>% rbind(df)
    }
  }
  df_expo = df_expo %>% 
    mutate(ratio_exposition = n / total,
           date =  as.Date(paste(annee, mois, jour, sep = "-"), format = "%Y-%m-%d")) %>% 
    group_by(date) %>% 
    summarise(ratio_exposition = mean(ratio_exposition, na.rm = TRUE)) %>% ungroup()
  df_expo$ratio_exposition = rollmean(df_expo$ratio_exposition, 30, align = "center", na.rm = TRUE, na.pad = TRUE)
  BOOL = (as.Date("1900-01-01", format = "%Y-%m-%d") <= df_expo$date) & (df_expo$date <= as.Date("1900-12-31", format = "%Y-%m-%d"))
  df_expo$ratio_exposition[!BOOL] = 0
  
  write.csv(df_expo, 
            paste0(path_projet, "\\df_expo.csv"), 
            fileEncoding = 'UTF-8', row.names = FALSE)
} else {
  df_expo = read.csv(paste0(path_projet, "\\df_expo.csv"), encoding = 'UTF-8')
  df_expo = df_expo %>% mutate(date = as.Date(date, format = "%Y-%m-%d"))
}

reg_data = reg_data %>% left_join(df_expo %>% select(date, ratio_exposition))

### Elections ###
if (scrap_data) {
  df_election = 0
  for (corpus in liste_corpus) {
    print(corpus)
    df = data_call(mot = 'élection', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df2 = data_call(mot = 'élections', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df$n = df$n + df2$n
    df$journal = corpus
    if (length(df_election) == 1) {
      df_election = df
    } else {
      df_election = df_election %>% rbind(df)
    }
  }
  df_election = df_election %>% 
    mutate(ratio_election = n / total,
           date =  as.Date(paste(annee, mois, jour, sep = "-"), format = "%Y-%m-%d")) %>% 
    group_by(date) %>% 
    summarise(ratio_election = mean(ratio_election, na.rm = TRUE)) %>% ungroup()
  df_election$ratio_election = rollmean(df_election$ratio_election, 30, align = "center", na.rm = TRUE, na.pad = TRUE)
  
  write.csv(df_election, 
            paste0(path_projet, "\\df_election.csv"), 
            fileEncoding = 'UTF-8', row.names = FALSE)
} else {
  df_election = read.csv(paste0(path_projet, "\\df_election.csv"), encoding = 'UTF-8')
  df_election = df_election %>% mutate(date = as.Date(date, format = "%Y-%m-%d"))
}
reg_data = reg_data %>% left_join(df_election %>% select(date, ratio_election))

### Parlimentary representation ###
if (scrap_data) {
  df_assemblee = 0
  for (corpus in liste_corpus) {
    print(corpus)
    df = data_call(mot = 'assemblée', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df2 = data_call(mot = 'chambre', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df3 = data_call(mot = 'législatif', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df4 = data_call(mot = 'sénat', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
    df$n = df$n + df2$n + df3$n + df4$n
    df$journal = corpus
    if (length(df_assemblee) == 1) {
      df_assemblee = df
    } else {
      df_assemblee = df_assemblee %>% rbind(df)
    }
  }
  df_assemblee = df_assemblee %>% 
    mutate(ratio_assemblee = n / total,
           date =  as.Date(paste(annee, mois, jour, sep = "-"), format = "%Y-%m-%d")) %>% 
    group_by(date) %>% 
    summarise(ratio_assemblee = mean(ratio_assemblee, na.rm = TRUE)) %>% ungroup()
  df_assemblee$ratio_assemblee = rollmean(df_assemblee$ratio_assemblee, 30, align = "center", na.rm = TRUE, na.pad = TRUE)
  
  write.csv(df_assemblee, 
            paste0(path_projet, "\\df_assemblee.csv"), 
            fileEncoding = 'UTF-8', row.names = FALSE)
} else {
  df_assemblee = read.csv(paste0(path_projet, "\\df_assemblee.csv"), encoding = 'UTF-8')
  df_assemblee = df_assemblee %>% mutate(date = as.Date(date, format = "%Y-%m-%d"))
}
reg_data = reg_data %>% left_join(df_assemblee %>% select(date, ratio_assemblee))

# Other data sources are retrieved

### French economic growth ###
# Retrieved on https://mail.google.com/mail/u/0/?tab=rm&ogbl#search/eurostar/FMfcgzGtwWHnmzLdXtwbxbGJnWlkKkMQ
data_eco = read.csv(paste0(path_projet, "data croissance.csv"))
data_eco = data_eco %>% 
  mutate(gdppc = as.numeric(gdppc),
         pop= as.numeric(pop)) %>% 
  mutate(gdp = gdppc*pop,
         year = as.numeric(year)) %>% 
  arrange(year) %>%
  mutate(gdp_growth = (gdp - lag(gdp))*100 / lag(gdp)) 
reg_data$year = year(reg_data$date)
reg_data = reg_data %>% left_join(data_eco %>% select(year, gdp_growth))

### President of council resignation ###
data_dem = read.csv(paste0(path_projet, "dates_pdc.csv"))
data_dem = data_dem %>%  
  mutate(date_demission = as.Date(date_demission, format = "%Y-%m-%d"))
reg_data$demission = FALSE
dems = na.omit(unique(data_dem$date_demission))
for (date in dems) {
  date = as.Date(date, format = "%Y-%m-%d")
  BOOL = (-10 <= date - (reg_data$date)) & ((date - reg_data$date) <= 0)
  reg_data$demission[BOOL] = TRUE
}

#### Regressions ####
# We make a regression with the identified events to test their significance

reg_data2 = reg_data %>% filter(journal==corpus) %>%
  select(date, journal,mois,annee,
         ratio, ratio_desais_trend, y5_trend, weights,gdp_growth,demission,
         Regime, guerre70, WWI, July14,
         ratio_dreyfus, ratio_panama, ratio_greve, ratio_cholera, ratio_guerre, 
         ratio_exposition, ratio_election, ratio_assemblee) %>% 
  mutate(across(!date & !journal & !Regime, as.numeric)) %>% 
  na.omit()

reg_data2[c(3,15:22)] = scale(reg_data2[c(3,15:22)])
reg = lm(ratio ~ I(as.numeric(date)) +
            WWI + July14  + gdp_growth + demission +
           ratio_dreyfus + ratio_greve + ratio_cholera + ratio_exposition + ratio_assemblee*Regime,
         data = reg_data2,weights = weights)
print(summary(reg))
print(corpus)
}
reg_data2$residuals = reg$residuals
plot_ly(reg_data2,x=~date,y=~residuals,mode="lines")

qqnorm(reg$residuals)
qqline(reg$residuals, col = "red", lty = 2)
plot(reg_data2$date, residuals(reg))

# This regression does not take into account the newspapers (may) behave differently : each event (may) has a different effect on each newspaper
# We repeat the regression by interacting each newspaper with each of the events
reg_data2 = reg_data %>% 
  select(date, journal, weights,
         ratio_desais_trend, y5_trend,
         guerre70, Regime, WWI, July14,
         ratio_dreyfus, ratio_panama, ratio_greve, ratio_cholera, ratio_guerre, ratio_exposition, ratio_election, ratio_assemblee,
         gdp_growth, demission) %>% 
  mutate(across(!date & !journal &!Regime, as.numeric)) %>% 
  na.omit()
# Here we name the newspapers we want to test
journaux = unique(reg_data2$journal)
# Here we name the variables we want to test
vars = c("ratio_dreyfus", "ratio_panama", "ratio_greve", "ratio_cholera", 
         "ratio_guerre", "ratio_exposition", "ratio_election", "ratio_assemblee")
# We create the regresion call 
reg_call = "(ratio_desais_trend + y5_trend) ~ Regime + WWI + guerre70"
for (var in vars) {
  for (journal in journaux) {
    # Here we name the pairs that are historically unrelevant
    BOOL_inacc = (var == 'ratio_dreyfus') & (journal == 'huma')
    if (!BOOL_inacc) {
      new_col = paste0(journal, '_', var)
      reg_data2 = reg_data2 %>% 
        mutate(!!as.name(new_col) := 0,)
      BOOL = reg_data2$journal == journal
      reg_data2[BOOL, new_col] = reg_data2[BOOL, var]
      reg_call = paste0(reg_call, " + ", journal, "_", var)
    }
  }
}
# We compute the regression
reg = lm(formula(reg_call),
         data = reg_data2, 
         weights = weights)
# Df to store the results of the for loop
reg_coefs <- data.frame(matrix(NA, nrow = length(vars), ncol = length(journaux)))
rownames(reg_coefs) <- vars
colnames(reg_coefs) <- journaux
reg_lower95 <- data.frame(matrix(NA, nrow = length(vars), ncol = length(journaux)))
rownames(reg_lower95) <- vars
colnames(reg_lower95) <- journaux
reg_upper95 <- data.frame(matrix(NA, nrow = length(vars), ncol = length(journaux)))
rownames(reg_upper95) <- vars
colnames(reg_upper95) <- journaux
journal_ = reg_data2$journal
journal_ = model.matrix( ~ journal_ + 0)
# We loop on the coefficients to retrieve the data
for (var in vars) {
  coefs = coef(reg)
  conf_interval = confint(reg, level = 0.95)
  for (journal in journaux) {
    BOOL_inacc = (var == 'ratio_dreyfus') & (journal == 'huma')
    if (!BOOL_inacc) {
      matching_names1 <- grepl(var, names(coefs))
      matching_names2 <- grepl(journal, names(coefs))
      BOOL = matching_names1 & matching_names2
      value = coef(reg)[BOOL]
      reg_coefs[var, journal] = value
      matching_names1 <- grepl(var, rownames(conf_interval))
      matching_names2 <- grepl(journal, rownames(conf_interval))
      BOOL = matching_names1 & matching_names2
      value = conf_interval[BOOL, '2.5 %']
      reg_lower95[var, journal] = value
      matching_names1 <- grepl(var, rownames(conf_interval))
      matching_names2 <- grepl(journal, rownames(conf_interval))
      BOOL = matching_names1 & matching_names2
      value = conf_interval[BOOL, '97.5 %']
      reg_upper95[var, journal] = value
    }
  }
}
reg_coefs$variable = rownames(reg_coefs)
reg_coefs <- reg_coefs %>%
  pivot_longer(cols = colnames(reg_coefs)[colnames(reg_coefs) !=  'variable'],
               names_to = "journal",
               values_to = "Coef")
reg_lower95$variable = rownames(reg_lower95)
reg_lower95 <- reg_lower95 %>%
  pivot_longer(cols = colnames(reg_lower95)[colnames(reg_lower95) !=  'variable'],
               names_to = "journal",
               values_to = "low")
reg_upper95$variable = rownames(reg_upper95)
reg_upper95 <- reg_upper95 %>%
  pivot_longer(cols = colnames(reg_upper95)[colnames(reg_upper95) !=  'variable'],
               names_to = "journal",
               values_to = "up")
plot_data = reg_coefs %>% left_join(reg_lower95) %>% left_join(reg_upper95)
plot_list = list()

# We plot the data
first = TRUE
for (var in vars) {
  temp = plot_data %>% filter(variable == var)
  p = ggplot(temp, aes(y = journal, x = Coef)) +
    geom_bar(stat = "identity", fill = "blue", alpha = 0.6) +
    geom_errorbar(aes(xmin = low, xmax = up), width = 0.2, color = "red") +
    labs(title = var)
  if (first) {
    p = p + theme(axis.text.y = element_text(angle = 45, hjust = 1))
    first = FALSE
  } else {
    p = p + theme(axis.text.y = element_blank(),
                  axis.title.y = element_blank())
  }
  plot_list = c(plot_list, list(p))
}
grid.arrange(grobs = plot_list, ncol = 8)








