##################################################################
##          R Tutorial for Causal Modeling Workshop:            ##
##    Introduction to Convergent Cross Mapping and S-map        ##
##                                                              ##
##                  Wednesday, April 16, 2025                   ##
##                                                              ##
##                    Dr. Kimberly Thompson                     ##
##                                                              ##
##     In this tutorial: S-mapping                              ##
##                                                              ##
##     R code adapted from:                                     ##
##                                                              ##
##     Deyle, E., May, R., Munch, S., Sugihara, G. 2016.        ##
##        Tracking and forecasting ecosystem interactions in    ##
##        real time. Proc. R. Soc. B 283: 2015258               ##
##                                                              ##                                ##
##                                                              ##
##################################################################

# We again will use simulated data of a community consisting of 
# a shared resource, R, two consumers of this resource: C1 and
# C2, and two predators of these consumers: P1 and P2.

# We will examine how all components of the system affect C2
# through time as a result of the changing state of the system.


###############################################
###                                         ###
###          Install/Load Packages          ###
###                                         ###
###############################################

# List of required packages
packages <- c("ggplot2", "viridis", "tidyverse", "Matrix", "quantreg")

# Install missing packages and load all
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Apply to all packages
invisible(sapply(packages, install_if_missing))

# Call the rEDM library loaded from the earlier session
library(rEDM)


###############################################
###                                         ###
###             Data Loading                ###
###                                         ###
###############################################

# Load the simulated data found in the Github Data folder
path <-     # Fill in the path to the data here
df <- read.csv(paste(path, "Simulated Abundance Data.csv", sep = ""),
                 header = TRUE)


###############################################
###                                         ###
###                S-Mapping                ###
###                                         ###
###############################################

# Smap is basically a weighted linear regression. What does this mean?
# It is an extension of ordinary least squares (OLS) regression. In OLS
# all points (i.e., observations) contribute equally; however as the name
# implies, weighted linear regression applied different weights to points.

# In the case of S-mapping, the weights are assigned based on how close 
# points are in the state space. Note that critically, two points can be 
# separated in time but still be close in location within the state space.

# To solve the weighted linear regression, S-mapping uses singular value
# decomposition (SVD), which is a method for factoring matrices. Details of
# this are beyond the scope of today's workshop, but more information can be 
# found in Deyle et al. 2016.
# 
# Deyle et al. 2016 put solving the weighted linear regression with SVD into
# a single function, which is what we will use.

lm_svdsolve <- function(y, x, ws, subset = seq_along(y)){
  x <- x[subset,]
  y <- y[subset]
  ws <- ws[subset]
  # prepended column of 1s for constant term in linear model
  A <- cbind(1, x) * ws
  A_svd <- svd(A)
  # >>REMOVE SMALL SINGULAR VALUES<<
  s <- A_svd$d
  s_inv <- matrix(0, nrow = dim(x)[2]+1, ncol = dim(x)[2]+1)
  for(i in seq_along(s))
  {
    if(s[i] >= max(s) * 1e-5)
      s_inv[i,i] <- 1/s[i]
  }
  coeff <- A_svd$v %*% s_inv %*% t(A_svd$u) %*% (ws * y)
  coeff <- t(coeff)
  colnames(coeff) <- c("const",colnames(x))
  return(coeff)
}


###############################################
###                                         ###
###        Step 1: Set up the data          ###
###             for S-mapping               ###
###############################################

# Make a matrix of the time series
df <- as.matrix(df)

# we want to focus on the effects on C2 (column 4)
targ_col.test <- 3

# Define the embedding, here values will correspond to the hypothesized effect
# while B.values will correspond to the hypothesized cause
Embedding <- c("R", "C1", "C2", "P1", "P2")
Edim <- length(Embedding)

# generates coeff names
coeff_names.test <- sapply(colnames(df),
                           function(x) paste("d", colnames(df)[targ_col.test],
                                             "d", x, sep = ""))

# Set up the lag: We have the target column, then we add a new column 
# where the first row has the 2nd row's observations, the 3rd row the 
# 2nd row's observations, and so on. Here we only add one lag.
block <- cbind(df[2:dim(df)[1],targ_col.test], 
               df[1:(dim(df)[1]-1),])

# So when you view 'block', you will see a column V1 whose values will 
# be the same as the values in column C2, but starting from the second
# observation. Because of this, you will also see that the length of our
# time series is now no longer 2000, but 1999. This is because we want
# to avoid having an NA value as there would not be a lagged value for 
# the time step prior to time 1.
head(block)


# Similar to the process in CCM, we will standardize the data, including
# the lagged column we just added.

# find the standard deviation of each column in the matrix
norm_consts <- apply(block, 2, function(x) sd(x))

# Standardize the time series in each column so that mean is 0 and SD is 1
block <- as.data.frame(apply(block, 2, function(x) (x-mean(x))/sd(x)))

# Set up additional parameters for Smap. We will use the entire series as the
# library (training) data and the entire series as the prediction (testing)
# data.

# library length = 1:number of rows
lib <- 1:dim(block)[1]

# prediction length = 1:number of rows
pred <- 1:dim(block)[1]

# Set up how to store the output in a blank array
coeff <- array(0, dim = c(length(pred), Edim))

# Assign column names to the blank array
colnames(coeff) <- coeff_names.test

# Convert to a dataframe
coeff <- as.data.frame(coeff)



###############################################
###                                         ###
###        Step 2: Determine the best       ###
###                value of theta           ###
###############################################

# Theta is a parameter that determines the weights given to neighboring points.
# When theta = 0, all points in the state space are treated linearly and given
# the same weight. 
# theta > 0 means that the system is nonlinear and points will 
# be weighted based on their proximity in the state space.

# The below method  for determining theta is based on Deyle et al. 2016's 
# supplement: The simplest way to choose appropriate 
# theta is to examine prediction error as a function of theta. Instead of using
# rho as our metric as we did with CCM, here we will examine
# the mean absolute error between s-map predictions and observations.

# Before estimating theta though, we again need to determine the best embedding
# dimension.

####################################
###                              ###
###  Step 2a: Determine the best ###
###       embedding dimension    ###
####################################

# We will determine the best embedding dimension for the component of the system
# which we are treating as an effect.

# In this case, that is C2.

# Define max E to test.
maxE <- 15

# Initialize a blank vector
rho_C2 <- numeric(maxE)

# Calculate rho for each possible embedding dimension
for (E in 1:maxE) {
  simplex_output <- simplex(block[ , targ_col.test], E = E, tp = 1)
  rho_C2[E] <- simplex_output$rho
}

# Extract the E where rho is at its max
E_C2 <- which.max(rho_C2)


####################################
###                              ###
###  Step 2b: Determine the best ###
###             theta            ###
####################################

# To determine the best theta, we will specify which observations to use
# as training data, and which to use as testing data. Because our results
# could be sensitive to how we divide the data (i.e., what % of data we 
# use for training and testing), we will trial different percent splits
# of our data and then average the resulting thetas for a more robust
# estimate.

# Empty vector to store theta values
theta.vals <- c()

# Test for ideal theta based on different training/testing % splits
for (r in seq(0.2, 0.9, by = 0.1)) {
  
  # Input training data
  E.lib = as.vector(c(1, round(length(block[, 1])*r)))
  
  # Testing data
  E.pred = as.vector(c(round(length(block[, 1])*r) + 1,
                       length(block[ , 1])) )
  
  # if the training testing split results in testing data with only one 
  # observation an error will occur - move to next r if this happens
  if(length(unique(E.pred)) == 1) {
    theta.vals <- c(theta.vals, NA)
    next
  }
  
  t.tmp <- s_map(block[ , targ_col.test], lib = E.lib,
                 pred = E.pred, E = E_C2)
  
  # Fill in the values where MAE is minimized. 
  theta.vals <- c(theta.vals, t.tmp$theta[t.tmp$mae == min(t.tmp$mae)][1])
  
  
  print(r)
  
} # end of r loop 

# Find the average of the theta values.
# Based on common practice,
# if average is under 2, we will use two decimal places, if average is 2 or
# above we will use integers to describe theta.

# In this case, we can see from the vector theta.vals that all the values for
# C2 are 8.

if(mean(theta.vals, na.rm = TRUE) < 2) {
  theta <- round(mean(theta.vals, na.rm = TRUE), digits = 2)
} else {
  if(mean(theta.vals, na.rm = TRUE) >= 2) {
    theta <- round(mean(theta.vals, na.rm = TRUE), digits = 0)
  }
}


###############################################
###                                         ###
###        Step 3: Perform the S-map        ###
###                                         ###
###############################################

# Loop over each prediction point (i.e., all points in 'block') to calculate
# the S-map locally weighted regression for each prediction point
# 1. Calculate the normalized Euclidean distance from the target point to the
#       other points in the state space.
# 2. These distances then become the weights for a weighted linear regression
#       which we solve using singular value decomposition.

# CALCULATE COEFFICIENTS
# We have already defined pred as the number of observations in our 'block'
# object
start.time <- Sys.time() # Takes about 20 seconds
for (ipred in 1:length(pred)){
  
  #target point is excluded from the fitting procedure
  libs = lib[-pred[ipred]]
  
  # CALCULATE WEIGHTS
  # Create a matrix that just repeats the values for the row of the 
  # time series that corresponds to the ipred index (less the column for
  # lagged version of the target column.)
  q <- matrix(as.numeric(block[pred[ipred],2:dim(block)[2]]),
              ncol=Edim, nrow=length(libs), byrow = T)
  
  # calculate euclidean distance for block row indices corresponding
  # to what's in libs
  # libs = row indices minus target row and for columns 2 through 3,
  # so minus the lagged column
  # subtract q (values of the target row), then square each value
  # then find sum of each row
  distances <- sqrt(rowSums((block[libs,2:dim(block)[2]] - q)^2))
  
  # Find the mean of the distances
  dbar <- mean(distances)
  
  # Calculate the weightings, based on theta and the distances
  Ws <- exp(-theta*distances/dbar)
  
  # Solving the regression with SVD
  svd_fit <- lm_svdsolve(block[libs,1], block[libs,2:dim(block)[2]], Ws)
  
  coeff[ipred,] <- svd_fit[-1]
  
  print(ipred)
  
} # end of ipred loop
end.time <- Sys.time()
end.time - start.time


###############################################
###                                         ###
###        Step 4: Plot the results         ###
###                                         ###
###############################################

# We will create a plot that shows each component's effect on C2 (even the 
# effect of C2 on itself)

# Add a time column to the coeff dataframe
coeff$Time <- seq(1, 1999, by = 1)

# Convert the dataframe to long form
coeff_long <- coeff %>%
  pivot_longer(
    cols = -Time,  # Keeps 'Time' column as identifier
    names_to = "Type",  # Stores original column names here
    values_to = "strength"  # Stores values here
  )

# To make the plot easier to evaluate, we will only plot a subset
# of the data, but you can feel free to play around with plotting
# different subsets of the entire time series.

test.plot.data <- coeff_long[coeff_long$Time <= 500, ]

# Create the plot
smap.plot <- ggplot() +
  geom_line(data = test.plot.data, aes(x = Time, y = strength,
                                   color = factor(Type)),
            linewidth = 1) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  scale_x_continuous(name = "Time", lim = c(0, 500),
                     breaks = c(0, 250, 500),
                     labels = c(0, 250, 500)) + 
  scale_y_continuous(name = "Interaction Strength", lim = c(-1.3, 1.6),
                     breaks = c(-1.2, -0.6, 0, 0.6, 1.2),
                     labels = c(-1.2, -0.6, 0, 0.6, 1.2)) +
  scale_color_manual(name = "", values = c("dC2dR" = viridis(100)[25],
                                "dC2dC1" = viridis(100)[45],
                                "dC2dC2" = viridis(100)[65],
                                "dC2dP1" = viridis(100)[85],
                                "dC2dP2" = viridis(100)[100]),
                     labels = c("dC2dR" = "\u2202C2/\u2202R",
                                "dC2dC1" = "\u2202C2/\u2202C1",
                                "dC2dC2" = "\u2202C2/\u2202C2",
                                "dC2dP1" = "\u2202C2/\u2202P1",
                                "dC2dP2" = "\u2202C2/\u2202P2")) +
  theme_minimal() + 
  theme(axis.text.x = element_text(size=22, face="bold")) +
  theme(axis.text.y = element_text(size=22, face="bold")) +
  theme(axis.title.x = element_text(size=22, face="bold", color="gray30")) +
  theme(axis.title.y = element_text(size=22, face="bold", color="gray30")) 
  # theme(legend.position = "none")



#############################################
###                                       ###
###        Try to repeat the analysis     ###
###        by examining the changing      ###
###   interactions for another component  ###
###              of the system            ###
#############################################


