#_____________________________________________________________________________________________________#
#___ Noise induced order and fitting stochastic models _______________________________________________#
#___ https://royalsocietypublishing.org/doi/10.1098/rstb.2019.0381 ___________________________________#
#___ January 2024 - Marco Fele _______________________________________________________________________#
#_____________________________________________________________________________________________________#

library("tidyverse")

rm(list = ls())

# Functions __________________________________________________________________________________________________________________________________________________________________________________________________________----

# I tried to make this function general 
gillespie_algorithm <- function(states, # vector of size equal to the number of possible states, and entries equal to the number of individuals in that state
                                parameters, # vector of size equal to the number of reactions, and entries equal to the reaction rate
                                omega, # parameter indicating the spatial scale of the system, conventionally equal to the group size for collective decision-making, adjustable for population dynamics
                                stoichiometry_reagents, # matrix of size (rows * columns) number of reactions * number of possible states, with each entry the stoichiometry of the reagents for a given reaction
                                stoichiometry_products) { # matrix of size (rows * columns) number of reactions * number of possible states, with each entry the stoichiometry of the products for a given reaction
  # Useful variables
  n_reactions <- nrow(stoichiometry_reagents)
  n_reagents <- ncol(stoichiometry_reagents)
  
  # Find reaction velocities (or propensity) through the stochastic law of mass action for small population sizes (from: Approximation and inference methods for stochastic biochemical kinetics—a tutorial review)
  velocities <- sapply(1:n_reactions, # return a vector of length equal to the number of reactions, with each entry the reaction velocity
                       function(reaction_number) {
                         reagents_needed <- stoichiometry_reagents[reaction_number, ] 
                         
                         velocity <- 0
                         if(all(states - reagents_needed >= 0)) { # reactions happen (change the value of velocity) only if there are enough reagents (which mean they arrive at zero if you account the reaction happening)
                           velocity <- parameters[reaction_number] * # all this shit comes from me simplifying the stochastic law of mass action, if you understand this kudos
                             1/omega^(sum(reagents_needed) - 1) * # include the resolution of the system
                             prod(c(sapply(1:n_reagents, function(reagent) { # for every reagent find the simplification of the factorial
                                      prod(states[reagent] - 0:(reagents_needed[reagent] - 1)) 
                                    }))[reagents_needed > 0]) #  this is a vector with one entry per reagent, you multiply for every reagent which is involved in the reaction
                         }
                         
                         return(velocity)
                       }) 
  if(all(velocities < 0.00001 & velocities > -0.00001)) return(c(NA, NA, states)) # in case nothing can happen anymore
  
  # Draw waiting time
  waiting_time <- rexp(n = 1, rate = sum(velocities))
  # Draw event
  event <- sample(x = n_reactions, size = 1, prob = velocities)
  # Update states
  new_states <- states - stoichiometry_reagents[event, ] + stoichiometry_products[event, ]
  # Returning event and states is redundant but lets keep it
  return(c(waiting_time, event, new_states))
}

# Parameters ________________________________________________________________________________________________________________________________________________________________________________________________________________----

# Simulation parameters
set.seed(666) 

max_transitions <- 100000
initial_states <- c(25, 25) 
n_replicates <- 20
group_size <- sum(initial_states)
n_species <- length(initial_states)

# Reactions parameters for voter model
parameters <- c(1, 1, 0.01, 0.01, 0.00, 0.00)
names(parameters) <- c("r_a", "r_b", "sigma_a", "sigma_b", "ternary_1", "ternary_2")

# Reactions parameters for higher order interaction model
parameters_ternary <- c(1, 1, 0.01, 0.01, 0.08, 0.08)
names(parameters_ternary) <- c("r_a", "r_b", "sigma_a", "sigma_b", "ternary_1", "ternary_2")

# Stoichiometry matrices
stoichiometry_reagents <- matrix(c(1, 1, # each row is a reaction, each column is a species
                                   1, 1,
                                   1, 0,
                                   0, 1,
                                   2, 1,
                                   1, 2), 
                                 ncol = n_species, 
                                 byrow = T) 
stoichiometry_products <- matrix(c(2, 0, # each row is a reaction, each column is a species
                                   0, 2,
                                   0, 1,
                                   1, 0,
                                   3, 0,
                                   0, 3), 
                                 ncol = n_species,
                                 byrow = T) 

# Voter model ______________________________________________________________________________________________________________________________________________________________________________________----

## Simulation ----

# Create containers to save results
data <- matrix(NA, ncol = 3 + n_species, 
               nrow = max_transitions * n_replicates)
colnames(data) <- c("replicate", "waiting_time", "event", paste("species", 1:n_species, sep = "_"))
state_columns <- 4:ncol(data)

# Stochastic simulation
for(replicate in 1:n_replicates) {
  print(replicate)
  replicate_to_index <- max_transitions * (replicate - 1) # just the way I use to save results when there are multiple replicates
  
  # Set initial conditions
  states <- initial_states # initial conditions
  data[1 + replicate_to_index, ] <- c(replicate, 0, NA, states) # record initial conditions
  
  # Simulate for one replicate
  for(transition in 2:max_transitions) { 
    # Update based on mr Gillespie
    data[transition + replicate_to_index, ] <- c(replicate, 
                                           gillespie_algorithm(states, parameters, omega = group_size,
                                                               stoichiometry_reagents, stoichiometry_products))
    states <- data[transition + replicate_to_index, state_columns]
  }
}

## Visualize results ----

# Modify data
data_m <- data |>
  as.data.frame() |>
  group_by(replicate) |>
  mutate(time = cumsum(waiting_time),
         duration = c(tail(waiting_time, -1), NA)) 

data_l <- data_m |>
  pivot_longer(cols = contains("species"),
               names_to = "species",
               values_to = "number") |>
  group_by(replicate, species) |>
  mutate(proportion = number / group_size) 

# Plot
ggplot(data_l |>
         filter(replicate %in% 1:20, species == "species_1")) +
  geom_step(aes(time, number, color = species, 
                group = interaction(species, replicate))) +
  xlim(c(0, 3000)) +
  facet_wrap(~replicate)

ggplot(data_l |>
         filter(!is.na(duration))) +
  geom_density(aes(number, color = species, weight = duration,
                   group = interaction(species, replicate)), 
               fill = NA)

## Manually calculate autocorrelation function (ACF) (just for fun, ACF is a bit tricky to calculate manually) ----

# Create a time series at one second resolution (I am not sure this is necessary but this is what I did)
data_resolution <- data_m |> 
  mutate(order = species_1,
         time = ceiling(time)) |> # to find state at one second resolution
  select(replicate, time, order) |> 
  group_by(replicate, time) |>
  filter(row_number() == n()) |> # also to find state at one second resolution
  group_by(replicate) |>
  mutate(time_difference = c(diff(time, lag = 1), 1)) |>
  uncount(time_difference) |> # fill in the blanks
  mutate(time = 0:max(time))

# Calculate ACF 

ACF_data <- data_resolution |>
  group_by(replicate) |>
  mutate(avg_order = mean(order),
         deviation = order - avg_order,
         avg_squared_deviation = mean(deviation ^ 2)) |>
  merge(data.frame(tau = seq(0, 300, by = 3))) |>
  mutate(other_time = time + tau) |>
  left_join(data_resolution |> 
              rename(other_order = order),
            by = join_by(replicate == replicate, 
                         other_time == time)) |>
  mutate(other_deviation = other_order - avg_order,
         correlation = deviation * other_deviation / avg_squared_deviation) |>
  group_by(replicate, tau) |>
  summarise(ACF = mean(correlation, na.rm = T),
            sample_size = n()) 

# Plot results 

ggplot(ACF_data) +
  geom_line(aes(tau, ACF, color = replicate, group = replicate)) +
  geom_smooth(aes(tau, ACF), color = "red", linewidth = 2) +
  geom_hline(aes(yintercept = 0)) +
  geom_hline(aes(yintercept = 0.025), lty = "dashed") +
  geom_hline(aes(yintercept = -0.025), lty = "dashed") +
  xlim(c(0, 200))
  
## Fit SDE ----

# Calculate deterministic and stochastic component

data_det_sto <- data_resolution |>
  group_by(replicate) |>
  mutate(order = order / group_size,
        dt_1 = c(tail(order, -1), NA)) |>
  group_by(replicate, order) |>
  mutate(diff = dt_1 - order,
         size = n(),
         first_moment = mean(diff / 1, na.rm = T),
         residual = diff - first_moment * 1) |>
  reframe(first_moment = unique(first_moment, na.rm = T),
          second_moment = mean(residual^2/1, na.rm = T),
          size = unique(size)) |>
  pivot_longer(cols = c("first_moment", "second_moment"),
               names_to = "effect",
               values_to = "value")|>
  group_by(order, effect) |>
  mutate(avg_moment = mean(value, na.rm = T))

# Plot result

ggplot(data_det_sto) + 
  geom_point(aes(order, value, color = replicate, group = replicate)) +
  #geom_line(aes(order, avg_moment), color = "red", linewidth = 2) +
  geom_smooth(aes(order, value, weight = size),
              formula = y ~ splines::bs(x, 4),
              color = "red") +
  facet_wrap(~effect)

# Compare to figure 2: the scale of the plot is different because he uses polarization




# Higher-order interactions model ______________________________________________________________________________________________________________________________________________________________________________________----

## Simulation ----

# Create containers to save results
data_ternary <- matrix(NA, ncol = 3 + n_species, 
               nrow = max_transitions * n_replicates)
colnames(data_ternary) <- c("replicate", "waiting_time", "event", paste("species", 1:n_species, sep = "_"))
state_columns_ternary <- 4:ncol(data_ternary)

# Stochastic simulation
for(replicate in 1:n_replicates) {
  print(replicate)
  replicate_to_index <- max_transitions * (replicate - 1)
  
  # Set initial conditions
  states <- initial_states # initial conditions
  data_ternary[1 + replicate_to_index, ] <- c(replicate, 0, NA, states) # record initial conditions
  
  # Simulate replicate
  for(transition in 2:max_transitions) { 
    data_ternary[transition + replicate_to_index, ] <- c(replicate, 
                                                         gillespie_algorithm(states, parameters_ternary, omega = group_size,
                                                                             stoichiometry_reagents, stoichiometry_products))
    states <- data_ternary[transition + replicate_to_index, state_columns]
  }
}

## Visualize results ----

# Modify data
data_m_ternary <- data_ternary |>
  as.data.frame() |>
  group_by(replicate) |>
  mutate(time = cumsum(waiting_time),
         duration = c(tail(waiting_time, -1), NA)) 

data_l_ternary <- data_m_ternary |>
  pivot_longer(cols = contains("species"),
               names_to = "species",
               values_to = "number") |>
  group_by(replicate, species) |>
  mutate(proportion = number / group_size) 

# Plot
ggplot(data_l_ternary |>
         filter(replicate %in% 1:20, species == "species_1")) +
  geom_step(aes(time, number, color = species, 
                group = interaction(species, replicate))) +
  xlim(c(0, 3000)) +
  facet_wrap(~replicate)

ggplot(data_l_ternary |>
         filter(!is.na(duration))) +
  geom_density(aes(number, color = species, weight = duration,
                   group = interaction(species, replicate)), 
               fill = NA)

## Manually calculate autocorrelation function (ACF) ----

# Create a time series at one second resolution

data_resolution_ternary <- data_m_ternary |> 
  mutate(order = species_1,
         time = ceiling(time)) |> # to find state at one second resolution
  select(replicate, time, order) |> 
  group_by(replicate, time) |>
  filter(row_number() == n()) |> # also to find state at one second resolution
  group_by(replicate) |>
  mutate(time_difference = c(diff(time, lag = 1), 1)) |>
  uncount(time_difference) |> # fill in the blanks
  mutate(time = 0:max(time))

# Calculate ACF 

ACF_data_ternary <- data_resolution_ternary |>
  group_by(replicate) |>
  mutate(avg_order = mean(order),
         deviation = order - avg_order,
         avg_squared_deviation = mean(deviation ^ 2)) |>
  merge(data.frame(tau = seq(0, 300, by = 3))) |>
  mutate(other_time = time + tau) |>
  left_join(data_resolution_ternary |> 
              rename(other_order = order),
            by = join_by(replicate == replicate, 
                         other_time == time)) |>
  mutate(other_deviation = other_order - avg_order,
         correlation = deviation * other_deviation / avg_squared_deviation) |>
  group_by(replicate, tau) |>
  summarise(ACF = mean(correlation, na.rm = T),
            sample_size = n()) # no autocorrelation when the time series ends

# Plot results 

ggplot(ACF_data_ternary) +
  geom_line(aes(tau, ACF, color = replicate, group = replicate)) +
  geom_smooth(aes(tau, ACF), color = "red", linewidth = 2) +
  geom_hline(aes(yintercept = 0)) +
  geom_hline(aes(yintercept = 0.025), lty = "dashed") +
  geom_hline(aes(yintercept = -0.025), lty = "dashed") +
  xlim(c(0, 200))

## Fit SDE ----

# Calculate deterministic and stochastic component

data_det_sto_ternary <- data_resolution_ternary |>
  group_by(replicate) |>
  mutate(order = order / group_size,
         dt_1 = c(tail(order, -1), NA)) |>
  group_by(replicate, order) |>
  mutate(diff = dt_1 - order,
         size = n(),
         first_moment = mean(diff / 1, na.rm = T),
         residual = diff - first_moment * 1) |>
  reframe(first_moment = unique(first_moment, na.rm = T),
          second_moment = mean(residual^2/1, na.rm = T),
          size = unique(size)) |>
  pivot_longer(cols = c("first_moment", "second_moment"),
               names_to = "effect",
               values_to = "value")|>
  group_by(order, effect) |>
  mutate(avg_moment = mean(value, na.rm = T))

# Plot result

ggplot(data_det_sto_ternary) +
  geom_point(aes(order, value, color = replicate, group = replicate)) +
  #geom_line(aes(order, avg_moment), color = "red", linewidth = 2) +
  geom_smooth(aes(order, value, weight = size),
              formula = y ~ splines::bs(x, 4),
              color = "red") +
  facet_wrap(~effect)

# Compare to figure 2: the scale of the plot is different because he uses polarization
