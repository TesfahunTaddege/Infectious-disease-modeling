# =====================================================================
#Project: Infectious Diseases Modeling ----
#Title: SI (Susceptible-Infectious) Model Simulation in R
#Author: Tesfahun Taddege (Epidemiologist) ----
#Organization: Amhara Public Health Institute, Bahir Dar, Ethiopia ----
#Date: Jan 23, 2023 ----
#Email: ttaddege@gmail.com---
#Description: This script simulates the SI model for infectious diseases using R.
#It models the spread of an infectious disease in a closed population
#where individuals can be either susceptible (S) or infectious (I).
# ============================================================

# 1. Load the required package for solving differential equations
if (!requireNamespace("deSolve", quietly = TRUE)) {
  install.packages("deSolve")
}
library(deSolve)

# 2. Define population and initial conditions
population <- 1000
initial_infected_si <- 1
initial_susceptible_si <- population - initial_infected_si

# 3. Define the SI model equations
# This function describes the rate of change between compartments.
si_model <- function(time, state, parameters) {
  # The 'with' function allows us to use variable names directly
  with(as.list(c(state, parameters)), {
    # Calculate the number of new infections
    new_infections <- beta * S * I / population
    
    # Differential equations
    dS <- -new_infections
    dI <-  new_infections
    
    # Return the rates of change
    return(list(c(dS, dI)))
  })
}

# 4. Set the model parameters
# beta: transmission rate (controls how quickly the disease spreads)
si_parameters <- c(beta = 0.5)

# 5. Set the initial state of the population
si_initial_state <- c(S = initial_susceptible_si, I = initial_infected_si)

# 6. Define the time sequence for the simulation
times <- seq(from = 0, to = 50, by = 1) # Simulate for 50 days

# 7. Solve the differential equations using the 'ode' function
si_output <- as.data.frame(ode(
  y = si_initial_state, 
  times = times, 
  func = si_model, 
  parms = si_parameters
))

# 8. Plot the results and save the file
# Create a PNG file
png("si_model_plot.png", width = 800, height = 600, res = 100)

# Create the plot
plot(si_output$time, si_output$S, type = "l", col = "blue", lwd = 2,
     ylim = c(0, population), xlab = "Time (days)", ylab = "Number of Individuals",
     main = "SI Model Simulation")
lines(si_output$time, si_output$I, type = "l", col = "red", lwd = 2)

# Add a legend
legend("right", legend = c("Susceptible", "Infectious"), 
       col = c("blue", "red"), lty = 1, lwd = 2)

# Close the PNG device
dev.off()

print("SI model simulation complete. Plot saved as si_model_plot.png")
