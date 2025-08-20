
# LOAD PACKAGES

library("nloptr")
library("tidyverse")

# LOAD FUNCTIONS

source("functions.r")

# ESTIMATE STATIC MODELS

all_frequencies <- c("tick", "01sec", "sec", "10sec", "min", "5min")
all_logdfs <- seq(from = -12, to = 12, length.out = 5001)

file_results <- "freedom_results.rda"

if (file.exists(file_results)) {
  
  load(file_results)
  
} else {
  
  freedom_table <- tibble(frequency = character(), sigma2_log = numeric(), df_log = numeric(), loglik = numeric())
  
  for (f in 1:length(all_frequencies)) {
    
    load(paste0("ibm_2024_", all_frequencies[f], ".rda"))
    
    y <- data$diff
    n <- length(y)
    
    y_table <- table(y)
    
    for (v in 1:length(all_logdfs)) {
      
      est <- est_static(y_table = y_table, df_log = all_logdfs[v])
      
      freedom_table <- freedom_table %>%
        add_row(frequency = all_frequencies[f], sigma2_log = est$sigma2_log, df_log = all_logdfs[v], loglik = est$loglik / n)
      
    }
    
  }
  
  freedom_table <- freedom_table %>%
    group_by(frequency) %>%
    mutate(best = loglik == max(loglik)) %>%
    mutate(local = loglik > lag(loglik) & loglik > lead(loglik)) %>%
    mutate(variance = ifelse(sigma2_log < -740, "Lower Bound", "Optimal")) %>%
    mutate(frequency = case_match(frequency, "tick" ~ "Ultra-High", "01sec" ~ "0.1-Second", "sec" ~ "1-Second", "10sec" ~ "10-Second", "min" ~ "1-Minute", "5min" ~ "5-Minute")) %>%
    mutate(frequency = factor(frequency, levels = c("Ultra-High", "0.1-Second", "1-Second", "10-Second", "1-Minute", "5-Minute")))
  
  save(freedom_table, file = file_results, compress = "xz")
  
}

# DRAW FIGURE 2

cols_rainbow <- c("#f44336", "#ffa500", "#cccc00", "#8fce00", "#6fa8dc", "#8e7cc3")
cols_rainbow_darker <- as.character(sapply(cols_rainbow, function(x) { colorRampPalette(c("#000000", x))(5)[4] }))

ggplot() +
  geom_line(data = freedom_table, aes(x = df_log, y = loglik, color = frequency, linetype = variance)) +
  geom_point(data = freedom_table %>% filter(local), aes(x = df_log, y = loglik, color = frequency)) +
  scale_color_manual(name = "Frequency", values = cols_rainbow) +
  scale_linetype_manual(name = "Scale Parameter", values = c("dashed", "solid")) +
  scale_x_continuous(breaks = seq(from = -10, to = 10, by = 5), minor_breaks = seq(from = -10, to = 10, by = 1)) +
  scale_y_continuous(breaks = c(0, 100, 200, 300), minor_breaks = seq(from = -40, to = 360, length.out = 21)) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(-40, 360)) +
  xlab("ln(ν)") +
  ylab("Log-Likelihood") +
  ggtitle("Log-Likelihood for Various Degrees of Freedom") +
  theme_bw() +
  theme(aspect.ratio = 1)

###
