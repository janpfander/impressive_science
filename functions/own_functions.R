# round numbers from models
rounded_numbers <- function(x) mutate_if(x, is.numeric, round, 3)

# get model outputs ready for inline text reporting
text_ready <- function(model_output) {
  
  result <- tidy(model_output, conf.int = TRUE)
  
  if ("effect" %in% colnames(result)) {
    # Filter for effect == "fixed" if the variable exists
    result <- result %>% 
      filter(effect == "fixed")
  } 
  
  result <- result %>% 
    # report p.value according to apa standards
    mutate(p.value = case_when(p.value < 0.001 ~ "< .001",
                               TRUE ~ paste0("= ", sprintf("%.3f", p.value))
    )
    ) %>% 
    # all other terms
    rounded_numbers() %>% 
    mutate(ci = glue::glue("[{conf.low}, {conf.high}]")) 
  
  if ("term" %in% colnames(result)) {
    # Filter for effect == "fixed" if the variable exists
    result <- result %>% 
      mutate(
        term = ifelse(term == "(Intercept)", "intercept", term),
        # if there is an interaction (recognized by ":"), name it just interaction
        term = str_replace(term, ".*:.*", "interaction")
      ) %>% 
      split(.$term)
  } 
  
  return(result)
}

# Function for splitting data along several variables (useful for inline reporting)
# taken from here: https://www.tjmahr.com/lists-knitr-secret-weapon/
super_split <- function(.data, ...) {
  dots <- rlang::enquos(...)
  for (var in seq_along(dots)) {
    var_name <- rlang::as_name(dots[[var]])
    .data <- purrr::map_depth(
      .x = .data,
      .depth = var - 1,
      .f = function(xs) split(xs, xs[var_name])
    )
  }
  .data
}


# Function for Exp. 3: 
# to filter data and run mixed model for a given outcome and comparison variable
run_mixed_model <- function(data, outcome) {
  # Filter data to include only the pair we care about
  data_filtered <- data %>% 
    filter(outcome == "knowledge" | outcome == {{outcome}})
  
  # We want to add a binary, numeric variable that identifies whether the outcome is 
  # knowledge or the other outcome of interest.
  data_filtered <- data_filtered %>% 
    mutate(outcome_binary = ifelse(outcome == "knowledge", 1, 0))
  
  # Run the mixed model with random intercepts
  model <- lmer(value ~ treatment * outcome_binary + (1 | id), data = data_filtered)
  
  return(model)
}






