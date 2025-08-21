
# LOAD PACKAGES

library("rugarch")
library("fGarch")
library("GAS")
library("gasmodel")
library("tidyverse")

# ESTIMATE GARCH MODELS

load("ibm_2024_5min.rda")

all_frequencies <- c("sec", "min")
all_packages <- c("rugarch", "fGarch", "GAS", "gasmodel")

file_results <- "garch_results.rda"

if (file.exists(file_results)) {
  
  load(file_results)
  
} else {
  
  garch_table <- tibble(frequency = character(), day = as.Date(numeric()), package = character(), parameter = character(), value = numeric())
  
  for (f in 1:length(all_frequencies)) {
    
    load(paste0("ibm_2024_", all_frequencies[f], ".rda"))
    
    all_days <- unique(data$day)
    
    for (d in 1:length(all_days)) {
   
      y <- filter(data, day == all_days[d])$diff
      n <- length(y)
      
      for (p in 1:length(all_packages)) {
        
        message("Frequency: ", all_frequencies[f], ", Day: ", all_days[d], ", Package: ", all_packages[p])
        
        if (all_packages[p] == "rugarch") {
          
          spec <- ugarchspec(variance.model = list(model = "sGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = T), distribution.model = "std")
          est <- ugarchfit(spec = spec, data = y)
          
          mu <- est@fit$solver$sol$pars["mu"]
          omega <- est@fit$solver$sol$pars["omega"]
          alpha <- est@fit$solver$sol$pars["alpha1"]
          phi <- est@fit$solver$sol$pars["beta1"]
          df <- est@fit$solver$sol$pars["shape"]
          loglik <- -min(est@fit$solver$sol$values) / n
          
        } else if (all_packages[p] == "fGarch") {
          
          est <- garchFit(data = y, cond.dist = "std", trace = F)
          
          mu <- coef(est)["mu"]
          omega <- coef(est)["omega"]
          alpha <- coef(est)["alpha1"]
          phi <- coef(est)["beta1"]
          df <- coef(est)["shape"]
          loglik <- -est@fit$llh / n
          
        } else if (all_packages[p] == "GAS") {
          
          spec = UniGASSpec(Dist = "std", GASPar = list(location = F, scale = T, shape = F))
          est <- UniGASFit(GASSpec = spec, data = y)
          
          mu <- coef(est)["kappa1"]
          omega <- coef(est)["kappa2"]
          alpha <- coef(est)["a2"]
          phi <- coef(est)["b2"]
          df <- getFilteredParameters(est)[1, 3]
          loglik <- est@Estimates$IC["llk"] / n
          
        } else if (all_packages[p] == "gasmodel") {
          
          est <- gas(y = y, distr = "t", regress = "sep", par_static = c(T, F, T))
          
          mu <- coef(est)["mean"]
          omega <- coef(est)["log(var)_omega"] * (1 - coef(est)["log(var)_phi1"])
          alpha <- coef(est)["log(var)_alpha1"]
          phi <- coef(est)["log(var)_phi1"]
          df <- coef(est)["df"]
          loglik <- as.numeric(logLik(est)) / n
          
        }
        
        garch_table <- add_row(garch_table, frequency = all_frequencies[f], day = all_days[d], package = all_packages[p], parameter = c("mu", "omega", "alpha", "phi", "df", "loglik"), value = c(mu, omega, alpha, phi, df, loglik))
        
      }
      
    }
    
  }
  
  save(garch_table, file = file_results, compress = "xz")
  
}

# REPORT TABLE 1

garch_table %>%
  mutate(frequency = factor(frequency, levels = c("sec", "min"))) %>%
  mutate(package = factor(package, levels = c("rugarch", "fGarch", "GAS", "gasmodel"))) %>%
  mutate(parameter = factor(parameter, levels = c("mu", "omega", "alpha", "phi", "df", "loglik"))) %>%
  group_by(frequency, package, parameter) %>%
  summarize(value = median(value)) %>%
  pivot_wider(names_from = c(frequency, package), values_from = value) %>%
  arrange(parameter)

###
