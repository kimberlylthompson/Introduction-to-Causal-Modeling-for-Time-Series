##################################################################
##         R Practical for Causal Inference Workshop:           ##
##    Introduction to Convergent Cross Mapping and S-map        ##
##                                                              ##
##                  Wednesday, April 16, 2025                   ##
##                                                              ##
##                    Dr. Kimberly Thompson                     ##
##                                                              ##
##     In this script: Simulating Data for CCM and S-map        ##
##                                                              ##
##     R code adapted from:                                     ##
##     Deyle et al. 2016. Tracking and forecasting ecosystem    ##
##        interactions in real time. Proceedings of the Royal   ##
##        Society B. 283.                                       ##
##################################################################

# We will simulate abundance data for 2 predators, 2 consumers, and 
# one shared resource.


###############################################
###                                         ###
###        Step 1: Define Parameters        ###
###             and Functions               ###
###############################################

# Set the random number seed so that our results (in which we will generate
# some random numbers) are reproducible
set.seed(186)

# Define the parameters of the 5-component coupled food chain model 
# We will use the same parameters as Deyle et al. 2016
params <- list(nu1 = 0.1, nu2 = 0.07, lambda1 = 3.2, lambda2 = 2.9,
               C1star = 0.5, C2star = 0.5, mu1 = 0.15, mu2 = 0.15,
               kappa1 = 2.5, kappa2 = 2.0, Rstar = 0.3, k = 1.2)

# nu1, nu2, mu1, mu2: individual metabolic rates
# lambda1, lambda2, kapp1, kapp2: individual ingestion rates
# C1star, C2star, Rstar: half-saturation densities of the consumers, and
#                        resource, respectively

# Function to generate time series of resource
dR <- function(R,C1,C2,P1,P2) R*(1-R/params$k) -
  params$mu1*params$kappa1*(C1*R)/(R+params$Rstar) -
  params$mu2*params$kappa2*(C2*R)/(R+params$Rstar)

# Function to generate time series of consumer 1
dC1 <- function(R,C1,C2,P1,P2) params$mu1*params$kappa1*(C1*R)/(R+params$Rstar) -
  params$nu1*params$lambda1*(P1*C1)/(C1+params$C1star) - params$mu1*C1

# Function to generate time series of consumer 2
dC2 <- function(R,C1,C2,P1,P2) params$mu2*params$kappa2*(C2*R)/(R+params$Rstar) -
  params$nu2*params$lambda2*(P2*C2)/(C2+params$C2star) - params$mu2*C2

# Function to generate time series of predator 1
dP1 <- function(R,C1,C2,P1,P2) params$nu1*params$lambda1*(P1*C1)/
  (C1+params$C1star) - params$nu1*P1

# Function to generate time series of predator 2
dP2 <- function(R,C1,C2,P1,P2) params$nu2*params$lambda2*(P2*C2)/
  (C2+params$C2star) - params$nu2*P2

# Function to concatentate time series of the simulated abundances
dF <- function(R,C1,C2,P1,P2) c(dR(R,C1,C2,P1,P2), dC1(R,C1,C2,P1,P2),
                                dC2(R,C1,C2,P1,P2), dP1(R,C1,C2,P1,P2),
                                dP2(R,C1,C2,P1,P2))



###############################################
###                                         ###
###        Step 2: Simulate Data            ###
###                                         ###
###############################################


# Define simulation parameters
tmax <- 10000     #  Total duration
tau <- 5          #  Sampling interval (we record data every 5 time units)
dt <- 1/100       #  Time step for numerical integration (small for accuracy)
burn <- 200       #  Burn-in period 

# Initialize data storate array with dimensions (tmax/tau rows x 5 columns)
# Columns will store values for 5 variables: R, C1, C2, P1, P2
d <- array(data = 0, dim = c(tmax/tau,5))
colnames(d) <-c('R', 'C1', 'C2', 'P1', 'P2')

# We run the model with initial conditions (R,C1, C2, P1, P2) = (1, 0.5, 0.8, 0.7, 0.8). 
# Note that we start with a burn period to ensure the dynamics of the attractor 
# manifold are captured.

# Set the initial conditions for the 5 variables in the system
# R = 1, C1 = 0.5, C2 = 0.8, P1 = 0.7, P2 = 0.8
Xi <- c(1, 0.5, 0.8, 0.7, 0.8)


### Run the system for the burn-in period without recording data
# This allows transient dynamics to settle before we start collecting data

idex <- 0     # Initialize iteration counter

while(idex < burn * (1/dt)){              # Convert burn time to number of steps
  idex <- idex+1;                         # Increment counter
  Xi <- Xi + dt*do.call(dF,as.list(Xi))
  # The above line does the following:
  # 1. Converts Xi to a list of arguments for function dF
  # 2. Calculates derivates using the functions we defined above
  # 3. Multiplies the results by the time step (this is Euler's method)
  # 4. Adds the results to the current state
}

### Run the main simulation loop

# Store initial conditions (after burn-in) as first row
d[1,] <- Xi

# Reset counters
idex <- 0      # Step counter (in units of dt)
tdex <- 2      # Data storage (next empty row in array 'd')

while(idex < tmax * (1/dt)){     # Run for tmax time units
  
  idex <- idex+1                 # Increment step counter
  
  # Update system using the same numerical integration method as above
  Xi <- Xi + dt*do.call(dF,as.list(Xi))
  
  # Sampling logic: only store data at intervals of 'tau' time units
  # Check if current time is exact multiple of tau AND that we have space in the array.
  # Some reminders of matrix math operators:
  # %% = divide x by y and return the remainder
  # && = evaluates two logical statements and returns a single value if they're both true
  if( (idex*dt) %% tau == 0 && tdex < tmax/tau){
    d[tdex,] <- Xi
    tdex <- tdex+1
  }
}

# Convert to dataframe
d <- as.data.frame(d)

###############################################
###                                         ###
###      Step 3: Save the Simulations       ###
###                                         ###
###############################################

# Saved to Github Data folder
path <- "00_Data/"
write.csv(d, file = "Simulated Abundance Data.csv",
            row.names = FALSE)

