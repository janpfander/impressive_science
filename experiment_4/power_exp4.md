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

# Data generating process

## Variables

Outcomes variables: 

- **competence**: “How competent do you think [entomologists] are?” (7-point scale: 1 = exceptionally incompetent, 2 = very incompetent, 3 = quite incompetent, 4 = neither competent nor incompetent, 5 = quite competent, 6 = very competent, 7 = exceptionally competent)

- **trust**:  “How much do you trust the discipline of [entomology]?” (7-point scale: 1 = no trust at all, 2 = almost no trust, 3 = very little trust, 4 = some trust, 5 = quite strong trust, 6 = very strong trust, 7 = absolute trust).

- **impressiveness**: "How impressive did you find the findings of the entomologists described in the text you have read?" (5-point scale: 1 = not impressive at all, 2 = not very impressive, 3 = neither impressive nor unimpressive, 4 = quite impressive, 5 = very impressive)

- **recall**. We divided the text into different information elements. For each participant, we calculate a recall score based on how many of these elements they mention in their open-ended answer. We will code 0 if an element is not mentioned at all, and 1 if it is, not to the exact wording, but containing all the information. In order to have a more fine grained rating of the exactitude with which each point is recalled, we will assign scores of 0.5 if a point is mentioned, but lacks certain key elements. For example, if the original text goes "Some flies can perceive images displayed for just three milliseconds (a thousandth of a second). Since the two vignettes contain a slightly different number of total information elements according to our evaluation grid (13 for archaeology, 11 for entomology), we use a relative measure for the final recall score, namely the share of mentioned elements among all possible elements. 

- **recall of impressing items**. After having given their post evaluation of scientist's competence and trust in the different disciplines, participants will be presented with the evaluation grid for the respective discipline they had been randomized to see. They will be asked to select the 6 (~46% of total information in vignette) most impressive elements in the case of archaeology, and the 5 (~45% of total information in vignette) most impressive elements in the case of entomology. For each participant, we will then compute the score they obtained on the element they subjectively rated to be the most impressive ones. As for the recall variable, to account for the different number of impressive elements between the two vignettes, the final personal perceived loss of impressiveness measure will be a relative one, namely the share of mentioned elements among all impressive elements.

Error:

Errors for the outcomes are generated from a multivariate normal distribution with mean vector $ \mu = (0, 0) $ and covariance matrix. For simplicity, we assume the same errors for all outcomes. 

$$
\Sigma = \begin{bmatrix}
1 & \rho \\
\rho & 1
\end{bmatrix}
$$

## Model

All pre-test measures are assumed to be scaled to have mean 0 and standard deviation 1. We use standardized effect sizes. 


$$
outcome_\text{post} = outcome_\text{pre}  + \beta + \epsilon_{i}
$$

# Functions

## Simulate a sample


``` r
# Load library
library(MASS)
```

```
## 
## Attaching package: 'MASS'
```

```
## The following object is masked from 'package:patchwork':
## 
##     area
```

```
## The following object is masked from 'package:dplyr':
## 
##     select
```

``` r
# Parameters
effect_size <- 0.2      # Small effect size (Cohen's d)
rho <- 0.2             # Weak correlation between pre- and post-scores
sigma <- 1             # Standard deviation of errors

# Function to simulate experiment data
simulate_experiment_data <- function(n, effect_size, rho, sigma) {
  
  # Generate pre-scores as standardized scores
  pre_scores <- rnorm(n, mean = 0, sd = 1)
  
  # Variance-covariance matrix for correlated errors
  cov_mx <- matrix(c(sigma^2, rho, rho, sigma^2), nrow = 2)
  
  # Simulate error terms
  errors <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = cov_mx)
  
  # Calculate post-scores: pre + effect size + error
  post_scores <- pre_scores + effect_size + errors[, 2]
  
  # Create a data frame
  data <- data.frame(
    participant_id = 1:n,
    condition = sample(c("archaeology", "entomology"), n, replace = TRUE),
    pre_competence = pre_scores,
    post_competence = post_scores,
    pre_trust = pre_scores,  # Same pre-scores for trust
    post_trust = post_scores,
    forgetting_score = rnorm(n, mean = -effect_size, sd = 1),  # Independent outcome
    impressiveness = rnorm(n, mean = effect_size, sd = 1), # more impressiveness expected in experiment
    impressive_forgetting_score = rnorm(n, mean = -effect_size, sd = 1)
  )
  
  return(data)
}

# Function to simulate validation data
simulate_validation_data <- function(n_validation, n_phase1, effect_size) {
  
  data <- data.frame(
    participant_id = 1:n_validation,
    condition = sample(c("archaeology", "entomology"), n_validation, replace = TRUE),
    text_id = sample(1:n_phase1, n_validation, replace = TRUE),  # Map to phase 1 participants
    impressiveness = rnorm(n_validation, mean = 0, sd = 1)  # less impressiveness expected in validation study
  )
  return(data)
}

# test
# experiment_data <- simulate_experiment_data(n = 200, effect_size = effect_size, rho = rho, sigma = sigma)
# validation_data <- simulate_validation_data(n_validation = 100, n_phase1 = 200, effect_size = effect_size)
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
  
  # H1a: Paired t-test on competence
  h1a <- t.test(experiment_data$post_competence, experiment_data$pre_competence, paired = TRUE)
  results <- rbind(results, data.frame(
    hypothesis = "H1a",
    test_type = "paired t-test",
    statistic = h1a$statistic,
    p_value = h1a$p.value,
    mean_difference = h1a$estimate
  ), row.names = NULL)
  
  # H1b: Paired t-test on trust
  h1b <- t.test(experiment_data$post_trust, experiment_data$pre_trust, paired = TRUE)
  results <- rbind(results, data.frame(
    hypothesis = "H1b",
    test_type = "paired t-test",
    statistic = h1b$statistic,
    p_value = h1b$p.value,
    mean_difference = h1b$estimate
  ), row.names = NULL)
  
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
  
  # H3a: Independent t-test on impressiveness ratings
  h3a <- t.test(
    x = experiment_data$impressiveness,
    y = validation_data$impressiveness,
    var.equal = TRUE
  )
  results <- rbind(results, data.frame(
    hypothesis = "H3a",
    test_type = "independent t-test",
    statistic = h3a$statistic,
    p_value = h3a$p.value,
    mean_difference = h3a$estimate[1] - h3a$estimate[2]
  ), row.names = NULL)
  
  # H3b: Test impressive forgetting score against zero
  if (shapiro.test(experiment_data$impressive_forgetting_score)$p.value > 0.05) {
    h3b <- t.test(experiment_data$impressive_forgetting_score, mu = 0)
    results <- rbind(results, data.frame(
      hypothesis = "H3b",
      test_type = "one sample t-test",
      statistic = h3b$statistic,
      p_value = h3b$p.value,
      mean_difference = h3b$estimate
    ), row.names = NULL)
  } else {
    h3b <- wilcox.test(experiment_data$impressive_forgetting_score, mu = 0, exact = FALSE)
    results <- rbind(results, data.frame(
      hypothesis = "H3b",
      test_type = "Wilcoxon",
      statistic = h3b$statistic,
      p_value = h3b$p.value,
      mean_difference = NA
    ), row.names = NULL)
  }
  
  return(results)
}


# test 
# experiment_data <- simulate_experiment_data(n = 200, effect_size = 0.2, rho = 0.2, sigma = 1)
# validation_data <- simulate_validation_data(n_validation = 100, n_phase1 = 200, effect_size = 0.2)
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
      experiment_data <- simulate_experiment_data(n = n, effect_size = effect_size, rho = rho, sigma = sigma)
      
      # Generate validation data
      validation_data <- simulate_validation_data(n_validation = n_validation, n_phase1 = n, effect_size = effect_size)
      
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
## # A tibble: 5 × 4
##   hypothesis iterations significant power
##   <chr>           <int>       <int> <dbl>
## 1 H1a                 5           4   0.8
## 2 H1b                 5           4   0.8
## 3 H2                  5           5   1  
## 4 H3a                 5           1   0.2
## 5 H3b                 5           4   0.8
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
                               rho = rho,
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
  n_validation = 200,
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
## Rows: 135 Columns: 7
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
## 1 H1a              1000         941 0.941   300                  1000
## 2 H1b              1000         941 0.941   300                  1000
## 3 H2               1000         924 0.924   300                  1000
## 4 H3b              1000         943 0.943   300                  1000
## # ℹ 1 more variable: effect_size <dbl>
```

# Single Study sample

To demonstrate the analytical procedure in the pre-registration, we simulate a single sample. 


``` r
experiment_data <- simulate_experiment_data(n = 200, effect_size = effect_size, rho = rho, sigma = sigma)

write_csv(experiment_data , paste0("data/", "simulated_experiment.csv"))


validation_data <- simulate_validation_data(n_validation = 100, n_phase1 = 200, effect_size = effect_size)

write_csv(experiment_data , paste0("data/", "simulated_validation.csv"))
```


