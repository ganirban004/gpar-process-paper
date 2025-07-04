##### IMPORT NECESSARY LIBRARIES #####
library(tidyverse)
library(lubridate)
library(irr)
library(lsr)
library(GoodmanKruskal)
library(MESS)
library(patchwork)





##### IMPORT DATASET AND DEFINE NECESSARY VARIABLES #####
aqi_daily <- read.csv(file.choose(), header = T)  #choose AQI_Dataset_Kolkata.csv
aqi_colors <- c("Good" = "darkgreen", "Satisfactory" = "lightgreen", "Moderate" = "yellow",
                "Poor" = "goldenrod", "Very Poor" = "orange", "Severe" = "red")
aqi_levels <- c("Good", "Satisfactory", "Moderate", "Poor", "Very Poor", "Severe")
aqi_daily <- aqi_daily %>% mutate(AQI_Category = factor(AQI_Category, levels = aqi_levels, ordered = T))





##### DAILY AQI AND AQI CATEGORY TIME SERIES PLOT #####
aqi_full <- aqi_daily %>% mutate(Date_Full = make_date(Year, match(Month, month.name), Date))

### AQI Time Series Plot ###
ggplot(aqi_full, aes(x = Date_Full, y = AQI)) + geom_line(color = "black") + 
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") + labs(x = "Date", y = "AQI") + 
  theme_minimal() + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))



### AQI Category Time Series Plot ###
ggplot(aqi_full, aes(x = Date_Full, y = Category, group = 1)) + geom_line(color = "black") + 
  scale_y_continuous(breaks = 0:5, labels = levels(aqi_daily$AQI_Category)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") + 
  labs(x = "Date", y = "AQI Category", color = "Category") + theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), legend.position = "bottom")



### Rate Evolution Graph ###
rate_data <- aqi_full %>% group_by(Date_Full, AQI_Category) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  complete(Date_Full, AQI_Category, fill = list(count = 0)) %>% arrange(Date_Full) %>% 
  group_by(AQI_Category) %>% arrange(Date_Full) %>% mutate(cum_count = cumsum(count)) %>% ungroup()
ggplot(rate_data, aes(x = Date_Full, y = cum_count, color = AQI_Category)) + geom_line(linewidth = 1) + 
  scale_color_manual(values = aqi_colors) + scale_x_date(date_breaks = "1 year", date_labels = "%Y") + 
  labs(x = "Date", y = "Cumulative Count of Days", color = "AQI Category") + theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), legend.position = "bottom")





##### DISTRIBUTON OF AQI CATEGORIES #####
aqi_marginal <- aqi_daily %>% group_by(AQI_Category) %>% summarise(days = n(), .groups = "drop") %>% 
  mutate(Percentage = days * 100 / sum(days)) %>% select(-days)
ggplot(aqi_marginal, aes(x = AQI_Category, y = Percentage)) + 
  geom_col(color = "black", fill = "skyblue", position = "dodge", width = 0.5) + 
  geom_text(aes(label = paste0(round(Percentage, 2), "%")), position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) + coord_cartesian(ylim = c(0, 50)) +
  labs(x = "Category", y = "Percentage of Days") + theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))





##### AUTO-ASSOCIATION PLOT #####
aqi <- aqi_daily$AQI_Category
lags <- 1:5

### Cohen's Kappa ###
cohen_kappa_lag <- function(series, h) {
  n <- length(series)
  if (n <= h) stop("Lag is too large for the length of the series")
  x <- series[1:(n - h)]
  y <- series[(h + 1):n]
  data_pair <- data.frame(rater1 = x, rater2 = y)
  result <- kappa2(data_pair, weight = "unweighted")
  return(result$value)
}
kappa_vals <- sapply(lags, function(h) cohen_kappa_lag(aqi, h))
cohens_kappa <- ggplot(data.frame(Values = kappa_vals, Lag = lags), aes(x = Lag, y = Values)) + 
  geom_point(shape = 18, color = "black", size = 3) + geom_line(color = "black", linewidth = 1) + 
  labs(title = expression("Cohen's " * kappa), x = "Lag (h)", y = expression(kappa)) + 
  theme_minimal() + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))



### Cramer's Nu ###
cramers_nu_lag <- function(series, h) {
  n <- length(series)
  if (h >= n) stop("Lag is too large for the length of the series.")
  x <- series[1:(n - h)]
  y <- series[(h + 1):n]
  tab <- table(x, y)
  return(cramersV(tab))
}
cramer_vals <- sapply(lags, function(h) cramers_nu_lag(aqi, h))
cramers_nu <- ggplot(data.frame(Values = cramer_vals, Lag = lags), aes(x = Lag, y = Values)) + 
  geom_point(shape = 18, color = "black", size = 3) + geom_line(color = "black", linewidth = 1) + 
  labs(title = expression("Cramer's " * nu), x = "Lag (h)", y = expression(nu)) + theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))



### Goodman Kruskal's Tau ###
gktau_lag <- function(series, h) {
  n <- length(series)
  if (h >= n) stop("Lag is too large for the length of the series.")
  x <- series[1:(n - h)]
  y <- series[(h + 1):n]
  df <- data.frame(X = x, Y = y)
  gk_matrix <- GKtauDataframe(df)
  return(gk_matrix["X", "Y"])
}
tau_vals <- sapply(lags, function(h) gktau_lag(aqi, h))
gk_tau <- ggplot(data.frame(Values = tau_vals, Lag = lags), aes(x = Lag, y = Values)) + 
  geom_point(shape = 18, color = "black", size = 3) + geom_line(color = "black", linewidth = 1) + 
  labs(title = expression("Goodman–Kruskal's " * tau), x = "Lag (h)", y = expression(tau)) + 
  theme_minimal() + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))



### Goodman-Krushkal's Gamma ###
gkgamma_lag <- function(series, h) {
  n <- length(series)
  if (h >= n) stop("Lag is too large for the series length")
  x <- series[1:(n - h)]
  y <- series[(h + 1):n]
  tab <- table(x, y)
  gk_gamma <- gkgamma(tab)
  return(gk_gamma$estimate)
}
gamma_vals <- sapply(lags, function(h) gkgamma_lag(aqi, h))
gk_gamma <- ggplot(data.frame(Values = gamma_vals, Lag = lags), aes(x = Lag, y = Values)) + 
  geom_point(shape = 18, color = "black", size = 3) + geom_line(color = "black", linewidth = 1) + 
  labs(title = expression("Goodman-Krushkal's " * gamma), x = "Lag (h)", y = expression(gamma)) + 
  theme_minimal() + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))



### Pearson Measure X[n]^2 ###
pearson_lag <- function(series, h) {
  n <- length(series)
  if (h >= n) stop("Lag too large for the series length.")
  x <- series[1:(n - h)]
  y <- series[(h + 1):n]
  tab <- table(x, y)
  test <- chisq.test(tab, correct = FALSE)
  return(list(statistic = test$statistic, p_value = test$p.value, df = test$parameter))
}
pearson_vals <- sapply(lags, function(h) pearson_lag(aqi, h)$statistic)
pearson_measure <- ggplot(data.frame(Values = pearson_vals, Lag = lags), aes(x = Lag, y = Values)) + 
  geom_point(shape = 18, color = "black", size = 3) + geom_line(color = "black", linewidth = 1) + 
  labs(title = expression("Pearson Measure " * X[n]^2), x = "Lag (h)", y = expression(X[n]^2)) + 
  theme_minimal() + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))



### Mutual information u ###
mutual_info_lag <- function(series, h) {
  n <- length(series)
  if (h >= n) stop("Lag is too large for the series length.")
  x <- series[1:(n - h)]
  y <- series[(h + 1):n]
  joint_table <- table(x, y)
  joint_prob <- prop.table(joint_table)
  px <- rowSums(joint_prob)
  py <- colSums(joint_prob)
  mi <- 0
  for (i in seq_along(px)) {
    for (j in seq_along(py)) {
      pxy <- joint_prob[i, j]
      if (pxy > 0) {
        mi <- mi + pxy * log(pxy / (px[i] * py[j]))
      }
    }
  }
  return(mi / log(2))  # convert to bits
}
mi_vals <- sapply(lags, function(h) mutual_info_lag(aqi, h))
mi <- ggplot(data.frame(Values = mi_vals, Lag = lags), aes(x = Lag, y = Values)) + 
  geom_point(shape = 18, color = "black", size = 3) + geom_line(color = "black", linewidth = 1) + 
  labs(title = expression("Mutual Information " * u), x = "Lag (h)", y = expression(u)) + 
  theme_minimal() + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))



### Combined Plots ###
print(cohens_kappa | cramers_nu | gk_tau) / (gk_gamma | pearson_measure | mi)





##### SEASONAL PLOT #####
aqi_seasonal <- aqi_daily %>% mutate(Season = case_when(
  Month %in% c("March", "April", "May") ~ "Summer",
  Month %in% c("June", "July", "August", "September") ~ "Monsoon",
  Month %in% c("October", "November") ~ "Autumn",
  Month %in% c("December", "January", "February") ~ "Winter")) %>% 
  mutate(Season = factor(Season, levels = c("Summer", "Monsoon", "Autumn", "Winter"))) %>% 
  group_by(Season, AQI_Category) %>% summarise(days = n(), .groups = "drop") %>% group_by(Season) %>% 
  mutate(Percentage = days * 100 / sum(days)) %>% select(-days) %>% 
  pivot_wider(names_from = AQI_Category, values_from = Percentage, values_fill = 0) %>% 
  pivot_longer(cols = -Season, names_to = "AQI_Category", values_to = "Percentage") %>% 
  mutate(AQI_Category = factor(AQI_Category, levels = aqi_levels)) %>% arrange(Season, AQI_Category)
ggplot(aqi_seasonal, aes(x = Season, y = Percentage, fill = AQI_Category)) +
  geom_col(color = "black", position = "dodge", width = 1) + 
  geom_text(aes(label = paste0(round(Percentage, 2), "%")), position = position_dodge(width = 1), 
            vjust = -0.5, size = 3) + coord_cartesian(ylim = c(0, 75)) +
  labs(x = "Season", y = "Percentage of Days") + scale_fill_manual(values = aqi_colors) + theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), legend.position = "bottom")





##### FESTIVITY PLOT #####
diwali_window <- 7  #take it as one-sided window
prediwali_window <- 15
postdiwali_window <- 15

aqi_festive <- aqi_daily %>% mutate(Date_full = make_date(Year, match(Month, month.name), Date))
diwali_dates <- aqi_festive %>% filter(Is.Diwali == "Diwali") %>% select(Year, Diwali_Date = Date_full)
aqi_festive <- aqi_festive %>% left_join(diwali_dates, by = "Year") %>% 
  mutate(
    Days_from_Diwali = as.integer(Date_full - Diwali_Date), 
    Festivity = case_when(
      Days_from_Diwali >= -(diwali_window + prediwali_window) & Days_from_Diwali < -diwali_window ~ "Pre-Diwali",
      Days_from_Diwali >= -diwali_window  & Days_from_Diwali <= diwali_window ~ "Diwali",
      Days_from_Diwali > diwali_window & Days_from_Diwali <= (diwali_window + postdiwali_window) ~ "Post-Diwali",
      TRUE ~ NA_character_))
aqi_festive <- na.omit(subset(aqi_festive, select = -Is.Diwali)) %>% group_by(Festivity, AQI_Category) %>% 
  summarise(days = n(), .groups = "drop") %>% group_by(Festivity) %>% 
  mutate(Fraction = days / sum(days)) %>% select(-days) %>% 
  pivot_wider(names_from = AQI_Category, values_from = Fraction, values_fill = 0) %>% 
  pivot_longer(cols = -c(Festivity), names_to = "AQI_Category", values_to = "Fraction") %>% 
  mutate(AQI_Category = factor(AQI_Category, levels = aqi_levels, ordered = T)) %>% 
  mutate(Festivity = factor(Festivity, levels = c("Pre-Diwali", "Diwali", "Post-Diwali"))) %>% 
  arrange(Festivity, AQI_Category)
ggplot(aqi_festive, aes(x = Festivity, y = Fraction * 100, fill = AQI_Category)) + 
  geom_col(color = "black", position = "dodge") + coord_cartesian(ylim = c(0, 75)) + 
  geom_text(aes(label = paste0(round(Fraction * 100, 2), "%")), position = position_dodge(width = 1), 
            vjust = -0.5, size = 3) + labs(x = "Period", y = "Percentage of Days") + 
  scale_fill_manual(values = aqi_colors) + theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), legend.position = "bottom")
