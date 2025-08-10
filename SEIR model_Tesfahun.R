# -------------------------------------------------------------------
#Project: Infectoius Diseases Modeling
#Title: SEIR (Susceptible-Exposed-Infectious-Recovered) Model Simulation in R
##Author: Tesfahun Taddege (Epidemiologist) ----
#Organization: Amhara Public Health Institute, Bahir Dar, Ethiopia ----
#Date: Feb 15, 2023 ----
#Email: ttaddege@gmail.com---
# -------------------------------------------------------------------

# 1. Load the required package
if (!requireNamespace("deSolve", quietly = TRUE)) {
  install.packages("deSolve")
}
library(deSolve)

# 2. Define population and initial conditions
population_seir <- 100000
initial_exposed_seir <- 10 # Start with some exposed individuals
initial_infected_seir <- 0
initial_recovered_seir <- 0
initial_susceptible_seir <- population_seir - initial_exposed_seir

# 3. Define the SEIR model equations
seir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + E + I + R # Total population
    
    # Calculate transitions between compartments
    new_infections <- beta * S * I / N
    progression_to_infectious <- sigma * E
    recoveries <- gamma * I
    
    # Differential equations
    dS <- -new_infections
    dE <-  new_infections - progression_to_infectious
    dI <-  progression_to_infectious - recoveries
    dR <-  recoveries
    
    # Return the rates of change
    return(list(c(dS, dE, dI, dR)))
  })
}

# 4. Set model parameters, including the latent period
gamma_seir <- 1/7      # Recovery rate (infectious period of 7 days)
sigma_seir <- 1/5      # Progression rate (latent period of 5 days)
R0_seir <- 2.5         # Basic Reproduction Number
beta_seir <- R0_seir * gamma_seir # Calculate transmission rate

seir_parameters <- c(beta = beta_seir, gamma = gamma_seir, sigma = sigma_seir)

# 5. Set the initial state of the population
seir_initial_state <- c(
  S = initial_susceptible_seir, 
  E = initial_exposed_seir, 
  I = initial_infected_seir, 
  R = initial_recovered_seir
)

# 6. Define the time sequence
times_seir <- seq(from = 0, to = 200, by = 1) # Simulate for 200 days

# 7. Solve the differential equations
seir_output <- as.data.frame(ode(
  y = seir_initial_state, 
  times = times_seir, 
  func = seir_model, 
  parms = seir_parameters
))

# 8. Plot the results and save the file
png("seir_model_plot.png", width = 800, height = 600, res = 100)

plot(seir_output$time, seir_output$S, type = "l", col = "blue", lwd = 2,
     ylim = c(0, population_seir), xlab = "Time (days)", ylab = "Number of Individuals",
     main = "SEIR Model Simulation (with Latent Period)")
lines(seir_output$time, seir_output$E, type = "l", col = "orange", lwd = 2)
lines(seir_output$time, seir_output$I, type = "l", col = "red", lwd = 2)
lines(seir_output$time, seir_output$R, type = "l", col = "darkgreen", lwd = 2)

legend("topright", legend = c("Susceptible", "Exposed", "Infectious", "Recovered"), 
       col = c("blue", "orange", "red", "darkgreen"), lty = 1, lwd = 2, cex = 0.8)

dev.off()

print("SEIR model simulation complete. Plot saved as seir_model_plot.png")