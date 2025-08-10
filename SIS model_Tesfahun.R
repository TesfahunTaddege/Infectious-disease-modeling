# -------------------------------------------------------------
# Project: Infectious Diseases Modeling
# Title: SIS (Susceptible-Infectious-Susceptible) Model Simulation in R
#Author: Tesfahun Taddege (Epidemiologist) ----
#Organization: Amhara Public Health Institute, Bahir Dar, Ethiopia ----
#Date: Jan 23, 2023 ----
#Email: ttaddege@gmail.com---
# This script simulates the SIS model for infectious diseases
# where individuals can become susceptible again after recovery.
# It uses the deSolve package to solve the differential equations
# and generates a plot of the results.
# -------------------------------------------------------------

# 1. Load the required package
if (!requireNamespace("deSolve", quietly = TRUE)) {
  install.packages("deSolve")
}
library(deSolve)

# 2. Define population and initial conditions
population <- 1000
initial_infected_sis <- 1
initial_susceptible_sis <- population - initial_infected_sis

# 3. Define the SIS model equations
sis_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Calculate new infections and recoveries
    new_infections <- beta * S * I / population
    recoveries <- gamma * I
    
    # Differential equations
    dS <- -new_infections + recoveries
    dI <-  new_infections - recoveries
    
    # Return the rates of change
    return(list(c(dS, dI)))
  })
}

# 4. Set the model parameters
# beta: transmission rate
# gamma: recovery rate (1 / duration of infection)
sis_parameters <- c(beta = 0.5, gamma = 0.1)

# 5. Set the initial state of the population
sis_initial_state <- c(S = initial_susceptible_sis, I = initial_infected_sis)

# 6. Define the time sequence
times <- seq(from = 0, to = 150, by = 1) # Simulate for 150 days

# 7. Solve the differential equations
sis_output <- as.data.frame(ode(
  y = sis_initial_state, 
  times = times, 
  func = sis_model, 
  parms = sis_parameters
))
# Convert the output to a data frame for easier handling
# 8. Plot the results and save the file
png("sis_model_plot.png", width = 800, height = 600, res = 100)

plot(sis_output$time, sis_output$S, type = "l", col = "blue", lwd = 2,
     ylim = c(0, population), xlab = "Time (days)", ylab = "Number of Individuals",
     main = "SIS Model Simulation (Reaching Endemic Equilibrium)")
lines(sis_output$time, sis_output$I, type = "l", col = "red", lwd = 2)

legend("right", legend = c("Susceptible", "Infectious"), 
       col = c("blue", "red"), lty = 1, lwd = 2)

dev.off()
# 9. Print a message indicating completion
print("SIS model simulation complete. Plot saved as sis_model_plot.png")