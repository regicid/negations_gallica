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

# Custom path
path_projet = "C:\\Users\\Dell\\Documents\\Projet Benoit"

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
# Rolling sd / mean without outliers
compute_sd_without_outliers <- function(temp, probs) {
  q <- quantile(na.omit(temp), probs = probs)  
  filtered_data <- temp[temp <= q] 
  q <- quantile(na.omit(temp), probs = 1-probs)  
  filtered_data <- temp[temp >= q] 
  sd_without_outliers <- sd(filtered_data, na.rm = TRUE)
  return(sd_without_outliers)
}
compute_mean_without_outliers <- function(temp, probs) {
  q <- quantile(na.omit(temp), probs = probs)  
  filtered_data <- temp[temp <= q]
  q <- quantile(na.omit(temp), probs = 1-probs)  
  filtered_data <- temp[temp >= q] 
  sd_without_outliers <- mean(filtered_data, na.rm = TRUE)
  return(sd_without_outliers)
}

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

# Gallicagram data import
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
# To control for this, we adopt the following strategy :
# 1 - Build observations
#     V(ratio) is observed on observation groups. These groups are based on the journal, the year and ranges/bins of the total variable
#     We only keep groups with at least 30 observations. We need big groups for a reliable measure of the variance, but numerous groups for a reliable number of observations of ratio
# 2 - Build a regression model 
#     V(ratio) = V(Sum(X_i)/T) Where X is the random variable "a written word is 'pas'" and T is the total number of words
#     => V(ratio) = V(X)*T/T^2 = V(X)/T if Xi are identical and independent. We will test this hypothesis later on.
#     => log(V(ratio)) = log(V(X)) - log(T)
#     We cannot compute V(X) = ratio*(1-ratio) because total (the sample size) is often too small. 
#     Based on our data, E(X) is around 0.004. To have a confidence interval of size 0.001 we would need total = 9604000000
#     Thus we resort to econometrics and we compute the model :
#     log(V(ratio)) + log(T) = f(X)
#     We then have : V(ratio) = exp(f(X) - log(T)) = exp(f(X))/T
reg_data = gallicagram_data %>% select(date, journal, ratio, n, total)
reg_data$total_bin_100 = exp(round(log(reg_data$total)*100)/100)
reg_data = reg_data %>%
  group_by(total_bin_100, journal, year(date)) %>%
  summarise(log_variance_bin = log(var(ratio, na.rm = TRUE)),
            group_size = n()) %>% 
  ungroup() %>%
  filter(group_size >= 30) %>% na.omit()
colnames(reg_data)[3] = 'year'
plot(reg_data$log_variance_bin, log(reg_data$total_bin_100))
reg <- lm(log_variance_bin ~ year*journal + log(total_bin_100),
          data = reg_data)
summary(reg)
confint(reg, level = 0.95) # 1 is in the 95% confidence interval of the log(total_bin_100) coeff => the Xi iid hypothesis and the model is validated
plot(residuals(reg))
qqnorm(reg$residuals)
qqline(reg$residuals, col = "red", lty = 2) 
# We can thus force this coefficient to 1 to compute f(X)
fX_reg <- lm(log_variance_bin + log(total_bin_100) ~ year*journal,
             data = reg_data)
summary(fX_reg)
plot(reg_data$year, residuals(fX_reg))
qqnorm(fX_reg$residuals)
qqline(fX_reg$residuals, col = "red", lty = 2) 

# We compute estimated stds and we compute associated weights for each of our journal-date observations
temp = gallicagram_data %>% select(total, journal, date)
temp$year = year(temp$date)
temp = temp[temp$journal %in% unique(fX_reg$model$journal),]
predicted_values <- predict(fX_reg, newdata = temp, type = "response")
temp$var_pred = exp(predicted_values)/temp$total
gallicagram_data = gallicagram_data %>% 
  left_join(temp %>% select(journal, date, var_pred))
gallicagram_data$weights = 1/sqrt(gallicagram_data$var_pred)
plot(density(gallicagram_data$weights, na.rm = TRUE))

# Now, the time series can be standardized
gallicagram_data$ratio = gallicagram_data$ratio / sqrt(gallicagram_data$var_pred) # Divide by estimated standard deviation
# Let's plot this new data
plot_list = list()
for (journal in liste_corpus) {
  BOOL = (gallicagram_data$journal == journal) & (year(gallicagram_data$date) > 1858) & (year(gallicagram_data$date) < 1916)
  data = gallicagram_data[BOOL,]
  new_plot = ggplot(data, aes(x = date)) +
    geom_line(aes(y = ratio))
  new_plot = new_plot +
    labs(title = journal) +
    labs(x = NULL, y = NULL) +
    scale_x_date(
      date_breaks = "1 year",  
      date_labels = "%Y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot_list = c(plot_list, list(new_plot))
}
grid.arrange(grobs = plot_list, ncol = 2)
# As expected, the data is very noisy. We need to remove seasonality if we want to go further

#############################
### SEASONALITY COMPUTING ###
#############################

clean_data = data.frame()
# !! 1 minute * nb of journal to run this loop
for (journal in liste_corpus) {
  
  print(journal) # To keep track
  
  # For loop initialisation
  data = gallicagram_data[gallicagram_data$journal == journal,]
  data = data[order(data$date),]

  # Trend level
  # The variable used to control for the trend level is a moving average centered on date with a width of several years
  # NB : For each trend / seasonality, outliers are removed as the can be considered as unindentified potential events non-related to seasonality 
  moving_average <- rollapply(data$ratio, 
                              width = as.integer(365.25*5),
                              align = "center",
                              FUN = compute_mean_without_outliers, probs = 0.90, 
                              fill = NA)
  data$y5_trend = moving_average
  
  # Month seasonality  
  # The variable used to control for the month seasonality is a moving average centered on date with a width of one month
  window_size = 30
  moving_average <- rollapply(data$ratio, 
                              width = window_size, 
                              align = "center", 
                              FUN = compute_mean_without_outliers, probs = 0.90, 
                              fill = NA)
  temp = data
  temp$s_mois = moving_average
  temp$date = temp$date + years(1)
  temp = temp[c("s_mois", "date")]
  data = data %>% left_join(temp[])
  
  # Week seasonality
  # The variable computes the average of the day effects within the year
  data$jour_s = weekdays(data$date)
  temp1 = data.frame()
  for (day in unique(data$jour_s)) {
    temp2 = data[data$jour_s == day,]
    window_size = as.integer(365.25/7)
    moving_average <- rollapply(temp2$ratio, 
                                width = window_size, 
                                align = "center", 
                                FUN = compute_mean_without_outliers, probs = 0.90, 
                                fill = NA)
    temp2$jour_s = moving_average
    temp2 = temp2[c("date", "jour_s")]
    temp1 = temp1 %>% rbind(temp2)
  }
  data = data %>% select(-jour_s) %>% left_join(temp1, by = "date")
  
  # Exceptional day of the year seasonality
  # We measure it as the multiple year average for the same day of the year
  data$jour_a = paste(month(data$date), day(data$date), sep = "-")
  temp1 = data.frame()
  for (day in unique(data$jour_a)) {
    temp2 = data[data$jour_a == day,]
    window_size = as.integer(365.25*5/365.25)
    moving_average <- rollapply(temp2$ratio, 
                                width = window_size, 
                                align = "center", 
                                FUN = compute_mean_without_outliers, probs = 0.90, 
                                fill = NA)
    temp2$jour_a = moving_average
    temp2 = temp2[c("date", "jour_a")]
    temp1 = temp1 %>% rbind(temp2)
  }
  data = data %>% select(-jour_a) %>% left_join(temp1, by = "date")
  
  data = data %>% 
    # We not remove the trend level from the data
    mutate(ratio_desais = ratio - y5_trend,
           s_mois_desais = s_mois - y5_trend,
           jour_s_desais = jour_s - y5_trend,
           jour_a_desais = jour_a - y5_trend) %>% 
    # We then remove the month seasonality
    mutate(ratio_desais = ratio_desais - s_mois_desais,
           jour_s_desais = jour_s_desais - s_mois_desais,
           jour_a_desais = jour_a_desais - s_mois_desais) %>% 
    # We then remove the week seasonality
    mutate(ratio_desais = ratio_desais - jour_s_desais,
           jour_a_desais = jour_a_desais - jour_s_desais) %>% 
    # We then remove the exceptional day seasonality
    mutate(ratio_desais = ratio_desais - jour_a_desais) %>% 
    select(-s_mois_desais, -jour_s_desais, -jour_a_desais) %>% 
    # We compute a smoothed variable to allow for graphic representation
    mutate(smoothed_ratio_desais_7 = rollmean(ratio_desais, 
                                              k = 7, 
                                              align = "center", 
                                              na.rm = TRUE, fill = NA),
           smoothed_ratio_desais_30 = rollmean(ratio_desais, 
                                               k = 30, 
                                               align = "center", 
                                               na.rm = TRUE, fill = NA))
  
  # End of the loop, data storing
  clean_data = clean_data %>% rbind(data)
}

# Let's plot this new data
plot_list = list()
for (journal in liste_corpus) {
  BOOL = (clean_data$journal == journal) & (year(clean_data$date) > 1858) & (year(clean_data$date) < 1916)
  data = clean_data[BOOL,]
  new_plot = ggplot(data, aes(x = date)) +
    geom_line(aes(y = ratio_desais + y5_trend))
  new_plot = new_plot +
    labs(title = journal) +
    labs(x = NULL, y = NULL) +
    scale_x_date(
      date_breaks = "1 year",  
      date_labels = "%Y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot_list = c(plot_list, list(new_plot))
}
grid.arrange(grobs = plot_list, ncol = 2)
# Now smoothed :
plot_list = list()
for (journal in liste_corpus) {
  BOOL = (clean_data$journal == journal) & (year(clean_data$date) > 1858) & (year(clean_data$date) < 1916)
  data = clean_data[BOOL,]
  new_plot = ggplot(data, aes(x = date)) +
    geom_line(aes(y = smoothed_ratio_desais_30 + y5_trend))
  new_plot = new_plot +
    labs(title = journal) +
    labs(x = NULL, y = NULL) +
    scale_x_date(
      date_breaks = "1 year",  
      date_labels = "%Y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot_list = c(plot_list, list(new_plot))
}
grid.arrange(grobs = plot_list, ncol = 2)
# Under the hypothesis that our data follows a normal law with as mean the moving average
# As our data is standardised, no need to compute a standard dev
clean_data = clean_data %>% 
  mutate(ratio_90 = qnorm(0.80, 
                          mean = y5_trend, 
                          sd = 1),
         ratio_10 = qnorm(0.20, 
                          mean = y5_trend, 
                          sd = 1))
plot_list = list()
for (journal in liste_corpus) {
  BOOL = (clean_data$journal == journal) & (year(clean_data$date) > 1858) & (year(clean_data$date) < 1916)
  data = clean_data[BOOL,]
  new_plot = ggplot(data, aes(x = date)) +
    geom_line(aes(y = smoothed_ratio_desais_30 + y5_trend)) +
    geom_line(aes(y = ratio_90), color = "red", linetype = "dotted") +
    geom_line(aes(y = ratio_10), color = "red", linetype = "dotted")
  new_plot = new_plot +
    labs(title = journal) +
    labs(x = NULL, y = NULL) +
    scale_x_date(
      date_breaks = "1 year",  
      date_labels = "%Y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot_list = c(plot_list, list(new_plot))
}
grid.arrange(grobs = plot_list, ncol = 2)

#############################
### EVENT REGRESSION TEST ###
#############################

# Analysing wether events have an influence on time series
# We will be more precise than in the the above graphic section as we will not resort on smoothing out the data

reg_data = clean_data

# We create the events we want to test
# 2nd Empire
reg_data$EmpII = 0 
BOOL = (as.Date("1854-01-14", format = "%Y-%m-%d") <= reg_data$date) & (reg_data$date <= as.Date("1870-10-27", format = "%Y-%m-%d"))
reg_data$EmpII[BOOL] = 1
# Franco Prussian War
# Commune and political instability
reg_data$guerre70 = 0 
BOOL = (as.Date("1870-07-19", format = "%Y-%m-%d") <= reg_data$date) & (reg_data$date <= as.Date("1871-05-28", format = "%Y-%m-%d"))
reg_data$guerre70[BOOL] = 1
# Third Republic
reg_data$RepIII = 0 
BOOL = (as.Date("1870-09-04", format = "%Y-%m-%d") <= reg_data$date) & (reg_data$date <= as.Date("1940-07-10", format = "%Y-%m-%d"))
reg_data$RepIII[BOOL] = 1
# WWI
reg_data$WWI = 0 
BOOL = (as.Date("1914-07-28", format = "%Y-%m-%d") <= reg_data$date) & (reg_data$date <= as.Date("1918-11-11", format = "%Y-%m-%d"))
reg_data$WWI[BOOL] = 1

# Event proxies
# Some of the events do not have strict timelines and occur with varying intensity
# To take this into account, we create variables that measure this intensity
# The ratio_ variables define the mean occurence of certain words inside the journal_list corpus
# They are then smoothed and set to 0 on irrelevant periods

### Dreyfus ###
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
BOOL = (as.Date("1894-09-01", format = "%Y-%m-%d") <= df_dreyfus$date) & (df_dreyfus$date <= as.Date("1907-06-01", format = "%Y-%m-%d"))
df_dreyfus$ratio_dreyfus[!BOOL] = 0
reg_data = reg_data %>% left_join(df_dreyfus %>% select(date, ratio_dreyfus))

### Panama scandal ###
df_panama = 0
for (corpus in liste_corpus) {
  print(corpus)
  df = data_call(mot = 'panama', debut = time_frame[1], fin = time_frame[2], resolution = "jour", corpus = corpus)
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
BOOL = (as.Date("1889-02-01", format = "%Y-%m-%d") <= df_panama$date) & (df_panama$date <= as.Date("1898-05-30", format = "%Y-%m-%d"))
df_panama$ratio_panama[!BOOL] = 0
reg_data = reg_data %>% left_join(df_panama %>% select(date, ratio_panama))
### Strikes ###
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
reg_data = reg_data %>% left_join(df_greve %>% select(date, ratio_greve))

### Cholera epidemics ###
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
reg_data = reg_data %>% left_join(df_cholera %>% select(date, ratio_cholera))

### War ###
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
reg_data = reg_data %>% left_join(df_guerre %>% select(date, ratio_guerre))

### 1900 Exposition universelle ###
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
reg_data = reg_data %>% left_join(df_expo %>% select(date, ratio_exposition))

### Elections ###
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
reg_data = reg_data %>% left_join(df_election %>% select(date, ratio_election))

### Parlimentary representation ###
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
reg_data = reg_data %>% left_join(df_assemblee %>% select(date, ratio_assemblee))

# Other data sources are retrieved

### French economic growth ###
# Retrieved on https://mail.google.com/mail/u/0/?tab=rm&ogbl#search/eurostar/FMfcgzGtwWHnmzLdXtwbxbGJnWlkKkMQ
data_eco = read.csv(paste0(path_projet, "\\data croissance.csv"))
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
data_dem = read.csv(paste0(path_projet, "\\dates_pdc.csv"))
data_dem = data_dem %>%  
  mutate(date_demission = as.Date(date_demission, format = "%Y-%m-%d"))
reg_data$demission = FALSE
dems = na.omit(unique(data_dem$date_demission))
for (date in dems) {
  date = as.Date(date, format = "%Y-%m-%d")
  BOOL = (-10 <= date - (reg_data$date)) & ((date - reg_data$date) <= 0)
  reg_data$demission[BOOL] = TRUE
}

### Regression ###
# We make a regression with the identified events to test their significance

reg_data2 = reg_data %>% 
  select(date, journal,
         ratio_desais, y5_trend,
         EmpII, guerre70, RepIII, WWI,
         ratio_dreyfus, ratio_panama, ratio_greve, ratio_cholera, ratio_guerre, 
         ratio_exposition, ratio_election, ratio_assemblee,
         gdp_growth, demission) %>% 
  mutate(across(!date & !journal, as.numeric)) %>% 
  na.omit()
reg = lm((ratio_desais + y5_trend) ~ journal +
           EmpII + guerre70 + RepIII + WWI +
           ratio_dreyfus + ratio_panama + ratio_greve + ratio_cholera + ratio_exposition + ratio_assemblee +
           gdp_growth + demission,
         data = reg_data2)
summary(reg)
qqnorm(reg$residuals)
qqline(reg$residuals, col = "red", lty = 2)

# This regression does not take into account the newspapers (may) behave differently : each event (may) has a different effect on each newspaper
# We repeat the regression by interacting each newspaper with each of the events
reg_data2 = reg_data %>% 
  select(date, journal,
         ratio_desais, y5_trend,
         EmpII, guerre70, RepIII, WWI,
         ratio_dreyfus, ratio_panama, ratio_greve, ratio_cholera, ratio_guerre, ratio_exposition, ratio_election, ratio_assemblee,
         gdp_growth, demission) %>% 
  mutate(across(!date & !journal, as.numeric)) %>% 
  na.omit()
# Here we name the newspapers we want to test
journaux = unique(reg_data2$journal)
# Here we name the variables we want to test
vars = c("ratio_dreyfus", "ratio_panama", "ratio_greve", "ratio_cholera", 
         "ratio_guerre", "ratio_exposition", "ratio_election", "ratio_assemblee")
# Here we name the pairs that are historically unrelevant
BOOL1 = (var == 'ratio_dreyfus') & (journal == 'huma')
BOOL_inacc = BOOL1
# We create the regresion call 
reg_call = "(ratio_desais + y5_trend) ~ "
for (var in vars) {
  temp = journal_*reg_data2[,var]
  colnames(temp) = paste0(colnames(temp), '_', var) 
  reg_data2 = reg_data2 %>% cbind(temp)
  for (journal in journaux) {
    if (!BOOL_inacc) {
      reg_call = paste0(reg_call, " + journal_", journal, "_", var)
    }
  }
}
# We compute the regression
reg = lm(formula(reg_call),
         data = reg_data2)
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








