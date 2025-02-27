---
title: "Power Experiment 3"
output:
  bookdown::html_document2:
      keep_md: yes
  bookdown::word_document2: default
  bookdown::pdf_document2: default
always_allow_html: true
---




``` r
# load required packages
library("lme4")        # model specification / estimation
library("lmerTest")    # get p-values for mixed models
library("broom.mixed") # extracting data from model fits 
library("tidyverse")   # data wrangling and visualisation
library("patchwork")    # combine plots

# ensure this script returns the same results on each run
set.seed(912783)
```

Statistical power is the probability by which we are able to detect an effect, assuming that it actually exists. The higher you set the benchmark for statistical power, the less likely you are to miss an effect in your sample.

The aim of a power analysis is find out how big a sample you need to achieve the statistical power you want. If there is an effect, more participants will make it more likely to detect it.

To give an overview of this power analysis, we will

1.  Make assumptions about the data generating process (i.e. a model and its parameters).
2.  "Draw" (i.e. simulate) a sample with a given sample size.
3.  Run the statistical analysis we plan to perform on that sample and store the estimates of the model.
4.  Calculate the power: Repeat steps two and three many times. See for how many samples the model recovers a statistically significant estimate for the parameter we set. Statistical power is just the ratio of significant samples among all samples we look at.
5.  Calculate power for different sample sizes: Repeat step four for various sample sizes. See at which sample size we achieve sufficiently high power.

# Functions

## Simulate a sample


``` r
# Parameters for temporary tests
effect_size <- 0.2      # Small effect size (Cohen's d)
sigma <- 1             # Standard deviation of errors

# Function to simulate experiment data
simulate_experiment_data <- function(n, effect_size, sigma) {
  
  # Generate scores for change in competence and trust as standardized scores
  change_competence <- rnorm(n, mean = effect_size, sd = sigma)
  change_trust <- rnorm(n, mean = effect_size, sd = sigma)
  
  # Generate forgetting scores
  forgetting_score <- rnorm(n, mean = -effect_size, sd = sigma)
  
  # Create a data frame
  data <- data.frame(
    participant_id = 1:n,
    discipline = sample(c("archaeology", "entomology"), n, replace = TRUE),
    change_competence = change_competence,  # Post-measure for competence
    change_trust = change_trust,            # Post-measure for trust
    forgetting_score = forgetting_score,    # Independent outcome
    impressiveness = rnorm(n, mean = effect_size, sd = sigma), # more impressiveness expected in experiment
    impressive_forgetting_score = rnorm(n, mean = -effect_size, sd = 1)
  )
  
  return(data)
}

# Function to simulate validation data
simulate_validation_data <- function(n_validation, n_phase1, effect_size, sigma) {
  
  data <- tibble(
    participant_id = 1:n_validation,
    discipline = sample(c("archaeology", "entomology"), n_validation, replace = TRUE),
    condition = sample(c("original", "recalled"), n_validation, replace = TRUE)
  ) %>%
    mutate(
      # assign a text id for recall texts
      text_id = ifelse(condition == "recalled",
                       sample(1:n_phase1, sum(.$condition == "recalled"), replace = TRUE),  # Map to phase 1 participants
                       NA),
      # make a numeric version for effect size calculations
      condition_numeric = ifelse(condition == "original", 0, 1),
      # 0 because no intercept (standardized variables), substract because we expect negative effect sizes for recall
      impressiveness = 0 - effect_size * condition_numeric + rnorm(n_validation, mean = 0, sd = sigma),
      competence = 0 - effect_size * condition_numeric + rnorm(n_validation, mean = 0, sd = sigma),
      trust = 0 - effect_size * condition_numeric + rnorm(n_validation, mean = 0, sd = sigma)
    )
  
  return(data)
}

# test
# experiment_data <- simulate_experiment_data(n = 200, effect_size = effect_size, sigma = sigma)
# validation_data <- simulate_validation_data(n_validation = 100, n_phase1 = 200, effect_size = effect_size, sigma = sigma)
```

## Run models


``` r
run_hypotheses_tests <- function(experiment_data, validation_data) {
  results <- data.frame(
    hypothesis = character(),
    test_type = character(),
    statistic = numeric(),
    p_value = numeric(),
    mean_difference = numeric(),
    stringsAsFactors = FALSE
  )
  
  # H1a: One-sample t-test for change in competence
  if (shapiro.test(experiment_data$change_competence)$p.value > 0.05) {
    h1a <- t.test(experiment_data$change_competence, mu = 0)
    results <- rbind(results, data.frame(
      hypothesis = "H1a",
      test_type = "one sample t-test",
      statistic = h1a$statistic,
      p_value = h1a$p.value,
      mean_difference = h1a$estimate
    ), row.names = NULL)
  } else {
    h1a <- wilcox.test(experiment_data$change_competence, mu = 0, exact = FALSE)
    results <- rbind(results, data.frame(
      hypothesis = "H1a",
      test_type = "Wilcoxon",
      statistic = h1a$statistic,
      p_value = h1a$p.value,
      mean_difference = NA
    ), row.names = NULL)
  }
  
  # H1b: One-sample t-test for change in trust
  if (shapiro.test(experiment_data$change_trust)$p.value > 0.05) {
    h1b <- t.test(experiment_data$change_trust, mu = 0)
    results <- rbind(results, data.frame(
      hypothesis = "H1b",
      test_type = "one sample t-test",
      statistic = h1b$statistic,
      p_value = h1b$p.value,
      mean_difference = h1b$estimate
    ), row.names = NULL)
  } else {
    h1b <- wilcox.test(experiment_data$change_trust, mu = 0, exact = FALSE)
    results <- rbind(results, data.frame(
      hypothesis = "H1b",
      test_type = "Wilcoxon",
      statistic = h1b$statistic,
      p_value = h1b$p.value,
      mean_difference = NA
    ), row.names = NULL)
  }
  
  # H2: Test forgetting score against zero
  if (shapiro.test(experiment_data$forgetting_score)$p.value > 0.05) {
    h2 <- t.test(experiment_data$forgetting_score, mu = 0)
    results <- rbind(results, data.frame(
      hypothesis = "H2",
      test_type = "one sample t-test",
      statistic = h2$statistic,
      p_value = h2$p.value,
      mean_difference = h2$estimate
    ), row.names = NULL)
  } else {
    h2 <- wilcox.test(experiment_data$forgetting_score, mu = 0, exact = FALSE)
    results <- rbind(results, data.frame(
      hypothesis = "H2",
      test_type = "Wilcoxon",
      statistic = h2$statistic,
      p_value = h2$p.value,
      mean_difference = NA
    ), row.names = NULL)
  }
  
  # H3: One-sample test for impressive forgetting score against zero
  if (shapiro.test(experiment_data$impressive_forgetting_score)$p.value > 0.05) {
    h3 <- t.test(experiment_data$impressive_forgetting_score, mu = 0)
    results <- rbind(results, data.frame(
      hypothesis = "H3",
      test_type = "one sample t-test",
      statistic = h3$statistic,
      p_value = h3$p.value,
      mean_difference = h3$estimate
    ), row.names = NULL)
  } else {
    h3 <- wilcox.test(experiment_data$impressive_forgetting_score, mu = 0, exact = FALSE)
    results <- rbind(results, data.frame(
      hypothesis = "H3",
      test_type = "Wilcoxon",
      statistic = h3$statistic,
      p_value = h3$p.value,
      mean_difference = NA
    ), row.names = NULL)
  }
  
  # H4a: Independent t-test for text impressiveness
  h4a <- t.test(
    impressiveness ~ condition,
    data = validation_data
  )
  results <- rbind(results, data.frame(
    hypothesis = "H4a",
    test_type = "independent t-test",
    statistic = h4a$statistic,
    p_value = h4a$p.value,
    mean_difference = h4a$estimate[[1]] - h4a$estimate[[2]]
  ))
  
  # H4b: Independent t-test for competence change
  h4b <- t.test(
    competence ~ condition,
    data = validation_data
  )
  results <- rbind(results, data.frame(
    hypothesis = "H4b",
    test_type = "independent t-test",
    statistic = h4b$statistic,
    p_value = h4b$p.value,
    mean_difference = h4b$estimate[[1]] - h4b$estimate[[2]]
  ))
  
  # H4c: Independent t-test for trust change
  h4c <- t.test(
    trust ~ condition,
    data = validation_data
  )
  results <- rbind(results, data.frame(
    hypothesis = "H4c",
    test_type = "independent t-test",
    statistic = h4c$statistic,
    p_value = h4c$p.value,
    mean_difference = h4c$estimate[[1]] - h4c$estimate[[2]]
  ))
  
  return(results)
}



# test 
# experiment_data <- simulate_experiment_data(n = 200, effect_size = 0.2, sigma = 1)
# validation_data <- simulate_validation_data(n_validation = 100, n_phase1 = 200, effect_size = 0.2, sigma = 1)
# 
# results <- run_hypotheses_tests(experiment_data, validation_data)
# results
```

## Calculate power


``` r
iterate <- function(n_iterations, n, n_validation, effect_size, rho, sigma) {
  
  # create data frame with model results for generated samples
  power_data <- 1:n_iterations %>% 
    purrr::map_df(function(x){
      # this is essentially a for loop - do the following for each 
      # element in 1:iterations
      
      # Generate experiment data
      experiment_data <- simulate_experiment_data(n = n, effect_size = effect_size, sigma = sigma)
      
      # Generate validation data
      validation_data <- simulate_validation_data(n_validation = n_validation, n_phase1 = n, effect_size = effect_size, sigma = sigma)
      
      # Run hypothesis tests
      iteration_results <- run_hypotheses_tests(experiment_data, validation_data)
      
      # Add iteration number
      iteration_results$iteration <- x
      
      # To keep track of progress
      if (x %% 50 == 0) {print(paste("iteration number ", x))}
      
      return(iteration_results)
      
    })
  
  return(power_data)
}

# test
# simulation_results <- iterate(
#   n_iterations = 5,
#   n = 200,
#   n_validation = 100,
#   effect_size = 0.2,
#   rho = 0.2,
#   sigma = 1
# )
# simulation_results
```



``` r
calculate_power <- function(power_data, alpha = 0.05) {
  power_summary <- power_data %>%
    group_by(hypothesis) %>%
    summarize(
      iterations = n(),
      significant = sum(p_value < alpha, na.rm = TRUE),
      power = significant / iterations
    ) %>%
    ungroup()
  
  return(power_summary)
}

# test the function
simulation_results <- iterate(
  n_iterations = 5,
  n = 200,
  n_validation = 100,
  effect_size = 0.2,
  rho = 0.2,
  sigma = 1
)
power <- calculate_power(simulation_results)
power
```

```
## # A tibble: 7 × 4
##   hypothesis iterations significant power
##   <chr>           <int>       <int> <dbl>
## 1 H1a                 5           4   0.8
## 2 H1b                 5           5   1  
## 3 H2                  5           3   0.6
## 4 H3                  5           4   0.8
## 5 H4a                 5           0   0  
## 6 H4b                 5           3   0.6
## 7 H4c                 5           3   0.6
```

## Power by sample size


``` r
power_by_sample_size <- function(sample_sizes, n_iterations, n_validation,
                                 effect_size, rho, sigma, alpha = 0.05) {

    
    # do the `calculate_power()` function for each sample size and store the results
    # in a new variable called `power`
    power <- purrr::map_df(sample_sizes, 
                           function(n){
                             # this is essentially a for loop - 
                             # do the following for each 
                             # element data$n_subj
                             
                             # To keep track of progress
                             print(paste("tested sample size = ", n))
                             
                             # Run the power simulation for the current sample size
                             power_data <- iterate(
                               n_iterations = n_iterations,
                               n = n,
                               n_validation = n_validation,
                               effect_size = effect_size,
                               sigma = sigma
                             )
                             
                             # Calculate power for the current sample size
                             power_summary <- calculate_power(power_data, alpha)
                             
                             # identify respective sample size
                             power_summary$n <- n
                             
                             return(power_summary)
                             
                           })
    
    # add some variables with information about simulation
    data <- power %>% 
      mutate(
        # add identifier variable for number of iterations
        iterations_per_sample = n_iterations)
    
    return(data)
}

# test
# Define parameters as a list
# parameters <- list(
#   sample_sizes = c(50, 100, 200),
#   n_iterations = 100,
#   n_validation = 100,
#   effect_size = 0.2,
#   rho = 0.2,
#   sigma = 1,
#   alpha = 0.05
# )
# 
# # Call the function using do.call
# result <- do.call(power_by_sample_size, parameters)
```


# Simulation



``` r
# list of parameters 
parameters <- list(
  sample_sizes = c(30, 50, 80, 100, 120, 150, 200, 250, 300),
  n_iterations = 1000,
  n_validation = 400,
  effect_size = 0.2,
  rho = 0.2,
  sigma = 1,
  alpha = 0.05
)
```

We run the simulation for different effect sizes (small = 0.2, medium = 0.5, large = 0.8).


``` r
# Define effect sizes to iterate over
effect_sizes <- c(0.2, 0.5, 0.8)

file_name <- "power_simulation.csv"

# only run analysis if a file with that name does not yet exists
if (!file.exists(paste0("data/", file_name))) {
  
  # Run simulation for each effect size
  results <- purrr::map_df(effect_sizes, function(effect_size) {
    # Update the effect size in params list
    parameters$effect_size <- effect_size
    
    # Use do.call to call the power_by_sample_size function with the updated params
    power_data <- do.call(power_by_sample_size, parameters)
    
    # Add the effect size to the results
    power_data$effect_size <- effect_size
    
    return(power_data)
  })
  
  write_csv(results, paste0("data/", file_name))
}
```


``` r
# read simulation data
power <- read_csv("data/power_simulation.csv")
```

```
## Rows: 189 Columns: 7
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr (1): hypothesis
## dbl (6): iterations, significant, power, n, iterations_per_sample, effect_size
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

### Plot results

We can plot the results of this power calculation.


``` r
ggplot(power %>% 
         mutate(effect_size = paste0("d = ", effect_size)), 
       aes(x = n, y = power, color = hypothesis)) +
  geom_line() + 
  geom_point(size = 1.5) + 
  # Add a horizontal line at 90%, our power_threshold
  geom_hline(aes(yintercept = .9), linetype = 'dashed') + 
  # Prettify!
  theme_minimal() + 
  scale_colour_viridis_d(option = "plasma") + 
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, by = 0.1)) + 
  labs(x = 'Sample Size', y = 'Power', 
       title = paste0("Power analysis"),
       caption = paste("iterations per \nsample size =", max(power$iterations))) +
  facet_wrap(~effect_size) 
```

![](power_exp4_files/figure-html/plot-power-1.png)<!-- -->

We can also inspect the precise sample size threshold.


``` r
# accuracy
power %>% 
  filter(effect_size == 0.2) %>%  
  filter(power >= 0.9) %>%
  arrange(n) 
```

```
## # A tibble: 4 × 7
##   hypothesis iterations significant power     n iterations_per_sample
##   <chr>           <dbl>       <dbl> <dbl> <dbl>                 <dbl>
## 1 H1a              1000         945 0.945   300                  1000
## 2 H1b              1000         926 0.926   300                  1000
## 3 H2               1000         934 0.934   300                  1000
## 4 H3               1000         931 0.931   300                  1000
## # ℹ 1 more variable: effect_size <dbl>
```

# Single Study sample

To demonstrate the analytical procedure in the pre-registration, we simulate a single sample. 


``` r
experiment_data <- simulate_experiment_data(n = 200, effect_size = effect_size, sigma = sigma)

write_csv(experiment_data , paste0("data/", "simulated_experiment.csv"))


validation_data <- simulate_validation_data(n_validation = 100, n_phase1 = 200, effect_size = effect_size, sigma = sigma)

write_csv(validation_data, paste0("data/", "simulated_validation.csv"))
```


