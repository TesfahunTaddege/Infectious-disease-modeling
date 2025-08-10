# ---------------------------------------------------------
#Project: Infectous Diseases Modeling
#Title: SIR (Susceptible-Infectious-Recovered) Model Simulation in R 
#(Measles Example)
##Author: Tesfahun Taddege (Epidemiologist) ----
#Organization: Amhara Public Health Institute, Bahir Dar, Ethiopia ----
#Date: Feb 15, 2023 ----
#Email: ttaddege@gmail.com---
# ---------------------------------------------------------

# 1. Load the required package
if (!requireNamespace("deSolve", quietly = TRUE)) {
  install.packages("deSolve")
}
library(deSolve)

# 2. Define population and initial conditions for a large-scale outbreak
population_sir <- 100000
initial_infected_sir <- 10
initial_recovered_sir <- 0
initial_susceptible_sir <- population_sir - initial_infected_sir - initial_recovered_sir

# 3. Define the SIR model equations
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + I + R # Total population size
    
    # Calculate new infections and recoveries
    new_infections <- beta * S * I / N
    recoveries <- gamma * I
    
    # Differential equations
    dS <- -new_infections
    dI <-  new_infections - recoveries
    dR <-  recoveries
    
    # Return the rates of change
    return(list(c(dS, dI, dR)))
  })
}

# 4. Set parameters for a highly contagious disease like Measles
R0 <- 15        # Basic Reproduction Number
gamma_sir <- 1/8 # Recovery rate (average infectious period of 8 days)
beta_sir <- R0 * gamma_sir # Calculate transmission rate from R0 and gamma

sir_parameters <- c(beta = beta_sir, gamma = gamma_sir)

# 5. Set the initial state of the population
sir_initial_state <- c(S = initial_susceptible_sir, I = initial_infected_sir, R = initial_recovered_sir)

# 6. Define the time sequence
times_sir <- seq(from = 0, to = 100, by = 1) # Simulate for 100 days

# 7. Solve the differential equations
sir_output <- as.data.frame(ode(
  y = sir_initial_state, 
  times = times_sir, 
  func = sir_model, 
  parms = sir_parameters
))

# 8. Plot the results and save the file
png("sir_model_plot.png", width = 800, height = 600, res = 100)

plot(sir_output$time, sir_output$S, type = "l", col = "blue", lwd = 2,
     ylim = c(0, population_sir), xlab = "Time (days)", ylab = "Number of Individuals",
     main = "SIR Model Simulation for a Measles-like Outbreak")
lines(sir_output$time, sir_output$I, type = "l", col = "red", lwd = 2)
lines(sir_output$time, sir_output$R, type = "l", col = "darkgreen", lwd = 2)

legend("topright", legend = c("Susceptible", "Infectious", "Recovered"), 
       col = c("blue", "red", "darkgreen"), lty = 1, lwd = 2, cex = 0.8)

dev.off()

print("SIR model simulation complete. Plot saved as sir_model_plot.png")
