
# LOAD PACKAGES

library("gasmodel")
library("nloptr")
library("ggpubr")
library("tidyverse")

# LOAD FUNCTIONS

source("functions.r")

# FIT DISTRIBUTION

all_frequencies <- c("tick", "01sec", "sec", "10sec", "min", "5min")
all_dfs <- c(0.01, 1, 2, 4, 8)

y_grid <- seq(from = -12, to = 12, by = 0.001)

file_results <- "distribution_results.rda"

if (file.exists(file_results)) {
  
  load(file_results)
  
} else {
  
  observation_table <- tibble(frequency = character(), y = numeric(), p = numeric())
  distribution_table <- tibble(frequency = character(), y = numeric(), d = numeric(), sigma2 = numeric(), df = numeric())
  
  for (f in 1:length(all_frequencies)) {
    
    load(paste0("ibm_2024_", all_frequencies[f], ".rda"))

    y <- data$diff
    n <- length(y)

    y_table <- table(y)

    prob_tibble <- tibble(y = as.integer(names(y_table)), p = as.numeric(y_table) / n) %>%
      filter(y >= -10 & y <= 10)
    
    observation_table <- observation_table %>%
      add_row(frequency = all_frequencies[f], prob_tibble)

    for (v in 1:length(all_dfs)) {
      
      est <- est_static(y_table = y_table, df_log = log(all_dfs[v]))
      
      dens <- distr_density(y_grid, f = c(0, exp(est$sigma2_log), exp(est$df_log)), distr = "t")
      
      distribution_table <- distribution_table %>%
        add_row(frequency = all_frequencies[f], y = y_grid, d = dens, sigma2 = exp(est$sigma2_log), df = exp(est$df_log))
      
    }
    
  }
  
  distribution_table <- distribution_table %>%
    mutate(label = paste0("ν = ", df, ", σ² = ", format(round(sigma2, 2), nsmall = 2)))
  
  save(observation_table, distribution_table, file = file_results, compress = "xz")
  
}
    
# DRAW FIGURE 1

col_grey <- "#bbbbbb"
cols_rainbow <- c("#f44336", "#ffa500", "#cccc00", "#8fce00", "#6fa8dc", "#8e7cc3")

titles <- paste(c("Ultra-High", "0.1-Second", "1-Second", "10-Second", "1-Minute", "5-Minute"), "Frequency")

figs <- list()

for (f in 1:length(all_frequencies)) {
  
  figs[[f]] <- ggplot() +
    geom_col(data = observation_table %>% filter(frequency == all_frequencies[f]), aes(x = y, y = p), fill = col_grey) +
    geom_line(data = distribution_table %>% filter(frequency == all_frequencies[f]), aes(x = y, y = d, color = label), alpha = 0.8) +
    scale_color_manual(name = "Parameters", values = colorRampPalette(c("#000000", cols_rainbow[f]))(length(all_dfs))) +
    scale_x_continuous(breaks = seq(from = -10, to = 10, by = 5), minor_breaks = seq(from = -10, to = 10, by = 1)) +
    scale_y_continuous(minor_breaks = seq(from = 0, to = 1, length.out = 21)) +
    coord_cartesian(xlim = c(-10, 10), ylim = c(0, 1)) +
    xlab("Price Change") +
    ylab("Probability Mass / Density") +
    ggtitle(titles[f]) +
    theme_bw() +
    theme(aspect.ratio = 1)
  
}

fig <- ggarrange(plotlist = figs, nrow = 3, ncol = 2, align = "hv")
fig <- annotate_figure(fig, top = text_grob("Distribution of Price Changes with Fitted Density of the Student's t-Distribution", size = 14))
fig

###
