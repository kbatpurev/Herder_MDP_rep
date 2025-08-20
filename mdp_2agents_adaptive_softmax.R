.libPaths()
.libPaths(new="C:/Data/RStudio")
library(tidyverse)
library(Hmisc)
library(dplyr)
library(quantmod) # for findPeaks fxn
library(MDPtoolbox)

# plotting libraries
library(corrplot)
library(RColorBrewer)
library(ggthemes)
library(ggpubr)
library(ggplot2)
library(Cairo)
library(extrafont)
extrafont::loadfonts()
library(kableExtra)
library(truncnorm)


#Set up and parameters
#---------------------------------------------------------------------
# land_grid
land_grid_x <- 5
land_grid_y <- 5
all_cells <- expand.grid(x = 0:(land_grid_x - 1), y = 0:(land_grid_y - 1))

# Actions (stocking rates)
actions <- c("<100", "200-400", "400-600", ">600")
action_income <- c("<100" = 0.8, "200-400" = 2.0, "400-600" = 2.5, ">600" = 3.5)
#action_penalty <- c("<100" = 1.0, "200-400" = 1.0, "400-600" = 1.0, ">600" = 0.7)

# Rangeland states
r_states <- c("intact", "good", "poor", "degraded")
r_code <- setNames(0:3, r_states)

# Weather probabilities
weather_probs <- c("drought" = 0.3, "rainy" = 0.2, "normal" = 0.5)

#Weather, transition and reward functions
#---------------------------------------------------------------------------
draw_weather <- function() {
  sample(names(weather_probs), size=1, prob = weather_probs)
}

transition_rangeland <- function(state, action, weather) {
  if (state == "degraded") return("degraded")  # terminal state
  
  probs <- NULL
  
  # Define probabilities based on action and weather
  if (action == ">600") {
    if (weather == "drought") {
      probs <- c("degraded" = 0.9, "poor" = 0.1)
    } else if (weather == "normal") {
      probs <- c("degraded" = 0.7, "poor" = 0.3)
    } else {  # rainy
      probs <- c("poor" = 0.8, "good" = 0.2)
    }
    
  } else if (action == "400-600") {
    if (weather == "drought") {
      probs <- c("poor" = 0.6, "good" = 0.4)
    } else if (weather == "normal") {
      probs <- c("good" = 0.7, "intact" = 0.3)
    } else {
      probs <- c("good" = 0.4, "intact" = 0.6)
    }
    
  } else if (action == "200-400") {
    if (weather == "drought") {
      probs <- c("good" = 0.3, "intact" = 0.7)
    } else if (weather == "normal") {
      probs <- c("good" = 0.1, "intact" = 0.9)
    } else {
      probs <- c("intact" = 1.0)
    }
    
  } else {  # action == "<100"
    probs <- c("intact" = 1.0)
  }
  
  return(sample(names(probs), 1, prob = probs))
}

compute_reward <- function(action, same_cell, revisit,cell_condition) {
  r <- action_income[action]
  if (same_cell) r <- r * 0.7
  if (revisit) r <- r * 0.7
  
  #A rangeland condition multiplier in reward (with the assumption that herders would benefit from a good quality pasture)
  condition_modifier <- dplyr::case_when(
    cell_condition == "intact"   ~ 1.0,
    cell_condition == "good"     ~ 0.9,
    cell_condition == "poor"     ~ 0.7,
    cell_condition == "degraded" ~ 0.4,
    TRUE ~ 1.0
  )
  
  r <- r * condition_modifier
  return(r)
}

softmax <- function(x, tau = 1.0) {
  exp_x <- exp((x - max(x)) / tau)
  exp_x / sum(exp_x)
}

#Agent Structure
#---------------------------------------------------------------------
create_agent <- function(name, tau = 1.0) {
  list(
    name = name,
    location = sample_n(all_cells, 1),
    last_cell = NA,
    reward_memory = list(),
    memory_count = list(),
    tau = tau
  )
}

#tau is the softmax temperature. 
#Lower tau → more greedy (e.g., 0.1)
#Higher tau → more exploratory (e.g., 2.0)
#tau = 1 is a moderate setting (default in many systems)

select_action <- function(agent) {
  candidates <- expand.grid(x = 0:(land_grid_x - 1), y = 0:(land_grid_y - 1), action = actions)
  candidates$reward <- apply(candidates, 1, function(row) {
    key <- paste(row[["x"]], row[["y"]], row[["action"]], sep = "_")
    if (!is.null(agent$reward_memory[[key]]) && agent$memory_count[[key]] > 0) {
      avg <- agent$reward_memory[[key]] / agent$memory_count[[key]]
    } else {
      avg <- 2 #if an agent hasn't been to a cell before and never tried an action before then they assume reward to be an average of 2
    }
    #this is clearly a naive and potentially overly optimistic approach. Alternative is doing a Bayesian smoothing based on prior knowledge.
    #through a small step like this:
    #prior_mean <- 1.5
    #prior_count <- 1
    #numerator <- ifelse(!is.null(agent$reward_memory[[key]]), agent$reward_memory[[key]], 0) + prior_mean * prior_count
    #denominator <- ifelse(!is.null(agent$memory_count[[key]]), agent$memory_count[[key]], 0) + prior_count
    #avg <- numerator / denominator
    
    # Safe revisit check
    same_as_last_cell <- isTRUE(agent$last_cell["x"] == as.numeric(row[["x"]]) &&
                                  agent$last_cell["y"] == as.numeric(row[["y"]]))
    if (same_as_last_cell) avg <- avg * 0.7
    return(avg)
  })
  
  probs <- softmax(candidates$reward, agent$tau)
  choice <- sample(1:nrow(candidates), 1, prob = probs)
  return(candidates[choice, ])
}


#Simulation loop
#-----------------------------------------------------------------------
simulate_adaptive <- function(steps = 100) {
  land_grid <- all_cells %>% mutate(rangeland = "intact") #assumes that all cells start from an intact condition = which is naive
  agent1 <- create_agent("A1")
  agent2 <- create_agent("A2")
  
  logs <- list()
  weather_trace <- character(steps)  # Store weather for debugging
  
  for (t in 1:steps) {
    weather <- draw_weather()
    weather_trace[t] <- weather
    cat(sprintf("Step %d — Weather: %s\n", t, weather))
    
    a1 <- select_action(agent1)
    a2 <- select_action(agent2)
    
    same_cell <- a1$x == a2$x && a1$y == a2$y
    
    a1_state <- land_grid %>% filter(x == a1$x, y == a1$y) %>% pull(rangeland)
    a2_state <- land_grid %>% filter(x == a2$x, y == a2$y) %>% pull(rangeland)
    
    a1_revisit <- isTRUE(agent1$last_cell["x"] == a1$x &&
                           agent1$last_cell["y"] == a1$y)
    a2_revisit <- isTRUE(agent2$last_cell["x"] == a2$x &&
                           agent2$last_cell["y"] == a2$y)
    
    r1 <- compute_reward(a1$action, same_cell, a1_revisit, a1_state)
    r2 <- compute_reward(a2$action, same_cell, a2_revisit, a2_state)
    
    new1 <- transition_rangeland(a1_state, a1$action, weather)
    new2 <- transition_rangeland(a2_state, a2$action, weather)
    
    land_grid[land_grid$x == a1$x & land_grid$y == a1$y, "rangeland"] <- new1
    land_grid[land_grid$x == a2$x & land_grid$y == a2$y, "rangeland"] <- new2
    
    #agent1 <- update_agent(agent1, a1, r1)
    #agent2 <- update_agent(agent2, a2, r2)
    
    logs[[t]] <- list(
      step = t,
      weather = weather,
      land_grid = land_grid,
      a1 = a1,
      a2 = a2
    )
  }
  
  return(list(logs = logs, weather_history = weather_trace))
}



results <- simulate_adaptive(step=140)

#Visualisation
#--------------------------------------------------------------------------
# Degradation percent over time
degradation_df <- map_dfr(results$logs, ~ .x$land_grid %>%
                            mutate(step = .x$step,
                                   degraded = rangeland == "degraded"))

degradation_summary <- degradation_df %>%
  group_by(step) %>%
  summarise(pct_degraded = mean(degraded) * 100)

ggplot(degradation_summary, aes(x = step, y = pct_degraded)) +
  geom_line(size = 1.2, color = "darkred") +
  labs(title = "Degradation Over Time", y = "% Degraded", x = "Step")

# Policy mix by cell
action_df <- map2_dfr(results$logs, results$weather_history, function(res, weather) {
  tibble(
    step = res$step,
    weather = weather,
    x1 = res$a1$x,
    y1 = res$a1$y,
    a1 = res$a1$action,
    x2 = res$a2$x,
    y2 = res$a2$y,
    a2 = res$a2$action
  )
})
action_long <- action_df %>%
  pivot_longer(
    cols = c(x1, y1, a1, x2, y2, a2),
    names_to = c(".value", "agent"),
    names_pattern = "([a-z]+)([12])"
  ) %>%
  mutate(agent = paste0("A", agent))
action_long <- action_long %>% rename(action = a)

policy_map <- action_long %>%
  group_by(agent, x, y) %>%
  count(action) %>%
  slice_max(n, n = 1, with_ties = FALSE)

ggplot(policy_map, aes(x = x, y = y, fill = action)) +
  geom_tile(color = "white") +
  facet_wrap(~agent) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Dominant Stocking Rate by Agent", x = "land_grid X", y = "land_grid Y")

#Value Iteration
#----------------------------------------------------------------------
initialize_value_function <- function(land_grid, actions, rangeland_states) {
  expand.grid(
    x = unique(land_grid$x),
    y = unique(land_grid$y),
    state = rangeland_states
  ) %>%
    mutate(value = 0)
}


get_possible_transitions <- function(current_state, action, weather) {
  # Return full distribution from the probabilistic version of transition_rangeland()
  if (current_state == "degraded") return(tibble(state = "degraded", prob = 1))
  
  # Example transition table for one action-weather combo
  probs <- switch(action,
                  ">600" = switch(weather,
                                  "drought" = c("degraded" = 0.9, "poor" = 0.1),
                                  "normal"  = c("degraded" = 0.7, "poor" = 0.3),
                                  "rainy"   = c("poor" = 0.8, "good" = 0.2)
                  ),
                  "400-600" = switch(weather,
                                     "drought" = c("poor" = 0.6, "good" = 0.4),
                                     "normal"  = c("good" = 0.7, "intact" = 0.3),
                                     "rainy"   = c("good" = 0.4, "intact" = 0.6)
                  ),
                  "200-400" = switch(weather,
                                     "drought" = c("good" = 0.3, "intact" = 0.7),
                                     "normal"  = c("good" = 0.1, "intact" = 0.9),
                                     "rainy"   = c("intact" = 1.0)
                  ),
                  "<100" = c("intact" = 1.0)
  )
  
  tibble(state = names(probs), prob = probs)
}

expected_value <- function(state, action, value_fn, weather, gamma = 0.95) {
  transitions <- get_possible_transitions(state, action, weather)
  
  # Add rewards (can be adjusted further)
  reward <- compute_reward(action, same_cell = FALSE, revisit = FALSE)  # use average case
  
  transitions %>%
    left_join(value_fn, by = c("state")) %>%
    summarise(ev = sum(prob * (reward + gamma * value))) %>%
    pull(ev)
}

value_iteration <- function(land_grid, actions, rangeland_states, iterations = 20, gamma = 0.95) {
  value_fn <- initialize_value_function(land_grid, actions, rangeland_states)
  
  for (i in 1:iterations) {
    value_fn <- value_fn %>%
      rowwise() %>%
      mutate(value = {
        weather <- draw_weather()  # stochastic weather per step
        max(sapply(actions, function(a) expected_value(state, a, value_fn, weather, gamma)))
      }) %>%
      ungroup()
  }
  
  return(value_fn)
}

extract_policy <- function(value_fn, actions, land_grid, rangeland_states, gamma = 0.95) {
  policy <- value_fn %>%
    rowwise() %>%
    mutate(best_action = {
      weather <- draw_weather()
      action_values <- sapply(actions, function(a) expected_value(state, a, value_fn, weather, gamma))
      actions[which.max(action_values)]
    }) %>%
    ungroup()
  
  return(policy)
}

value_iteration <- function(land_grid, actions, rangeland_states, iterations = 20, gamma = 0.95) {
  value_fn <- initialize_value_function(land_grid, actions, rangeland_states)
  
  for (i in 1:iterations) {
    value_fn <- value_fn %>%
      rowwise() %>%
      mutate(value = {
        weather <- draw_weather()
        max(sapply(actions, function(a) expected_value(state, a, value_fn, weather, gamma)))
      }) %>%
      ungroup()
  }
  
  return(value_fn)
}

land_grid <- expand.grid(
  x = 1:5,
  y = 1:5
) %>%
  mutate(rangeland = "intact")

value_fn <- value_iteration(land_grid, actions, r_states)

policy <- extract_policy(value_fn, actions, land_grid, r_states)

