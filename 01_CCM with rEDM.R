##################################################################
##           R Tutorial for Causal Modeling Workshop:           ##
##    Introduction to Convergent Cross Mapping and S-map        ##
##                                                              ##
##                  Wednesday, April 16, 2025                   ##
##                                                              ##
##                    Dr. Kimberly Thompson                     ##
##                                                              ##
##     In this tutorial: CCM with rEDM package                  ##
##                                                              ##
##     R code adapted from:                                     ##
##     Ye, H., Clark, A., Deyle, E., and Sugihara, G. 2019.     ##
##        rEDM: An R package for Empirical Dynamic Modeling     ##
##        and Convergent Cross Mapping. retrieved from          ##
##        https://ha0ye.github.io/rEDM/articles/rEDM.html       ##
##                                                              ##
##     Sugihara, G., Park, J., Deyle, E., Saberski, E.,         ##
##        Smith, C., and Ye, H. 2023. Empirical Dynamic         ##
##        Modeling. retrieved from https://cran.r-project.org   ##
##        /web/packages/rEDM/vignettes/rEDM-tutorial.pdf        ##
##                                                              ##
##     Chang, C., Ushio, M., Hsieh, C. 2017. Empirical dynamic  ##
##        modeling for beginners. Ecological Research 32.       ##
##        10.1007/s11284-017-1469-9                             ##
##                                                              ##
##     Deyle, E., May, R., Munch, S., Sugihara, G. 2016.        ##
##        Tracking and forecasting ecosystem interactions in    ##
##        real time. Proc. R. Soc. B 283: 2015258               ##
##                                                              ##
##     Clark et al. 2015. Spatial convergent cross mapping      ##
##        to detect causal relationships in short time series.  ##
##        Ecology 96(5)                                         ##
##        Society B. 283.                                       ##
##                                                              ##
##################################################################

# The data we will use is a simulated community consisting of 
# a shared resource, R, two consumers of this resource: C1 and
# C2, and two predators of these consumers: P1 and P2.

# We will examine whether there is a causal process between
# Predator 2 and Consumer 2 (testing for bidirectional 
# causality).

# Note that for multivariate time series, it is generally best to 
# standardize the values so that each variable is on a similar 
# scale.


###############################################
###                                         ###
###          Install/Load Packages          ###
###                                         ###
###############################################

# List of required packages
packages <- c("devtools", "multispatialCCM", "tidyverse", "ggplot2",
              "scales", "Metrics", "Kendall", "cocor", "ggridges")

# Install missing packages and load all
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Apply to all packages
invisible(sapply(packages, install_if_missing))


# Install the rEDM package - there are different versions of the 
# package and we all want to be working from the same one:

# So we first check if an existing version of rEDM is installed
# and if so, we remove it
if ("rEDM" %in% rownames(installed.packages())) {
  remove.packages("rEDM")
}

# Then we install rEDM from Github
devtools::install_github("ha0ye/rEDM")

# Load the library
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
###             Data Preparation            ###
###                                         ###
###############################################

# Add a time index: For rEDM, the first column must be Time
df$Time <- seq(1, nrow(df), by = 1)

# Move time to the first column
df <- df[ , c(6, 1:5)]

# Inspect the data
str(df)

# We will scale the variables so that they each have a mean of 0
# and a standard deviation of 1
cols_to_scale <- c("R", "C1", "C2", "P1", "P2")
df[cols_to_scale] <- scale(df[cols_to_scale])

# We can verify the mean and standard deviation by looking at 
# the summary
summary(df)


###############################################
###                                         ###
###      Step 1: Determine Embedding        ###
###             Dimension                   ###
###############################################

# Embedding dimension (E) represents the number of lagged coordinates used 
# to reconstruct the state space. The optimal E can be found using simplex 
# projection, which tests a time series' ability to predict its own dynamics, 
# either with defined training and testing datasets or through leave-one-out 
# cross-validation. We examine the correlation of actual and predicted values
# with the Pearson correlation coefficient, represented by rho. Higher rho 
# values indicate better predictive skill.

# In rEDM, the simplex function automatically performs 
# leave-one-out cross-validation when the library (lib: data to be used for
# training) and prediction (pred: data to be used for testing) parameters
# are not specified.

# In the case of leave-one-out cross validation all points in the time 
# series, except the point to be predicted, are used as training data.

# This step must be performed for both time series when testing for
# bidirectional causality.

# Define maximum E to test
# Theoretically maxE = floor(N/(tau*(E-1))) where N is time series length,
# and tau is the length of the time lag.
# However, practically maxE does not exceed 15 for most time series. 
# The default maxE in rEDM is 10. 
# For short time series maxE <= sqrt(N) where N is again the time series
# length.
# The optimal E is often <= true system dimension + noise

# We will test E values up to 15.
maxE <- 15

# Create empty vectors to store results
rho_P2 <- numeric(maxE)
rho_C2 <- numeric(maxE)

# We will find the best E value for both time series

# Note tp = prediction horizon: Specifies how many steps ahead to forecast
#           (e.g., tp = 3 predicts 3 time steps into the future). It affects
#           prediction targets.
#           Default value is 1, but can be any integer. Positive values 
#           forecast forward; negative values predict backward in time.

# Note tau = time lag: Determines the spacing between coordinates 
#            (e.g., tau = 2 uses every 2nd point). It affects state space 
#            reconstruction.
#            tau is typically 1 (the default in the below function) but can
#            be any positive integer <= time series length

### Example 1: Without any specification of library (training)
###  and prediction (testing) data --> assessment with 
###  leave-one-out cross validation

for (E in 1:maxE) {
  # For predator 2 time series
  simplex_output_A <- rEDM :: simplex(df$P2, E = E, tp = 1)
  rho_P2[E] <- simplex_output_A$rho
  
  # For consumer 2 time series
  simplex_output_B <- simplex(df$C2, E = E, tp = 1)
  rho_C2[E] <- simplex_output_B$rho
}

warnings()
# Note: We receive warnings that there is overlap between our lib and pred
# data (because we did not specify either parameter). As a result, 
# the function will use leave-one-out cross validation to evaluate the 
# results of its predictions.
# Exclusion radius = 0 means that when the function searches for neighbors
# to make a prediction, it excludes the point that is being predicted itself
# from being considered as a neighbor.

# Find optimal embedding dimension (E with maximum rho)
E_P2_loo <- which.max(rho_P2)
E_C2_loo <- which.max(rho_C2)

# Examine the results visually
plot(rho_P2, type = 'l', xlab = "Embedding Dimension", ylab = "Rho Value",
     main = "Ability of P2 ts to predict its own dynamics:\n Leave-One-Out")
abline(v = E_P2_loo, col = "red")

plot(rho_C2, type = 'l', xlab = "Embedding Dimension", ylab = "Rho Value",
     main = "Ability of C2 ts to predict its own dynamics:\n Leave-One-Out")
abline(v = E_C2_loo, col = "red")


### Example 2: specifying library (testing) and prediction (testing) regions

# Create empty vectors to store results
rho_P2_2 <- numeric(maxE)
rho_C2_2 <- numeric(maxE)

# Note: lib/pred can be non-contiguous (e.g., lib = "1 100 151 200")
lib <- c(1,100)   # First 100 points for training / reconstruction
pred <- c(201, 500) # Points 201-500 for testing / prediction

for (E in 1:maxE) {
  # For predator 2 time series
  simplex_output_A2 <- rEDM :: simplex(df$P2, lib = lib, pred = pred,
                                      E = E, tp = 1)
  rho_P2_2[E] <- simplex_output_A2$rho
  
  # For consumer 2 time series
  simplex_output_B2 <- simplex(df$C2, lib = lib, pred = pred,
                               E = E, tp = 1)
  rho_C2_2[E] <- simplex_output_B2$rho
}

# Find optimal embedding dimension (E with maximum rho)
E_P2_lib <- which.max(rho_P2_2)
E_C2_lib <- which.max(rho_C2_2)

# Examine the results visually
plot(rho_P2_2, type = 'l', xlab = "Embedding Dimension", ylab = "Rho Value",
     main = "Ability of P2 ts to predict its own dynamics:\n defining Lib")
abline(v = E_P2_lib, col = "red")

plot(rho_C2_2, type = 'l', xlab = "Embedding Dimension", ylab = "Rho Value",
     main = "Ability of C2 ts to predict its own dynamics:\n defining Lib")
abline(v = E_C2_lib, col = "red")

# In this case the optimal embedding dimension for P2
# remains the same when using leave-one-out cross validation and
# defining the training (library) and testing (prediction) data,
# while we get different results for C2. 

# Note that when defining the training (library) and
# testing (prediction) data that the predictive ability (shown
# by the rho value) decreased steeply as the Embedding Dimension
# increased; whereas in leave-one-out cross validation, the rho values
# were high across the entire range of embedding dimensions.

# This makes sense as we have much more data with which to train the 
# predictions when we use leave-one-out cross validation.

# Before we move on to step 2, we can clean up our workspace and 
# simplify our names for the optimal embedding dimensions of P2 and C2.
# We will proceed with the embedding dimensions determined from 
# specifying the training and testing data.

E_P2 <- E_P2_lib
E_C2 <- E_C2_lib

rm(simplex_output_A, simplex_output_A2, simplex_output_B,
   simplex_output_B2, E, E_C2_lib, E_C2_loo, E_P2_lib, E_P2_loo,
   lib, pred, maxE, rho_P2, rho_P2_2, rho_C2, rho_C2_2,
   cols_to_scale, packages, install_if_missing)



###############################################
###                                         ###
###      Step 2: Test for nonlinearity      ###
###         and stochastic noise            ###
###############################################

# CCM requires coupled nonlinear systems - neither purely random nor 
# completely deterministic.
# This means that while short-term prediction is 
# often possible, information about the state of the system is 
# continually lost over time, hindering long-term prediction.

# In this step, we verify that:
# 1. Predictive skill decays with increasing prediction horizon 
#    (indicating that the system is neither completely random
#    or comletely deterministic)
# 2. The system is nonlinear

# To verify that the predictive skill decays with increasing prediction
# horizon, we again use the simplex function but now we examine 
# multiple tp values for the optimal embedding dimension of each 
# time series that we identified in Step 1.

# Define the library (training) and prediction (testing) ranges
lib <- c(1, 100)
pred <- c(201, 500)

# Test time series P2
P2_test <- simplex(df$P2, 
                  E = E_P2, 
                  tp = 1:10,  # Test prediction horizons 1-10
                  tau = 1,
                  silent = TRUE)  # Suppresses progress messages

# Test time series C2  
C2_test <- simplex(df$C2,
                  E = E_C2,
                  tp = 1:10,
                  tau = 1,
                  silent = TRUE)


# The result is a dataframe that we can use to examine how predictive
# ability (rho) changes as a function of prediction horizon (tp)
plot(P2_test$tp, P2_test$rho, type = 'l', 
     xlab = "Prediction Horizon", ylab = "Rho Value",
     main = "Verification of Prediction Decay for P2 time series")

plot(C2_test$tp, C2_test$rho, type = 'l', 
     xlab = "Prediction Horizon", ylab = "Rho Value",
     main = "Verification of Prediction Decay for C2 time series")

# We see that for both time series, the predictive ability (quantified
# by rho), generally decreases as the prediction horizon increases.

# With this part of Step 2, we have confirmed that the system is neither
# completely stochastic, or completely deterministic.


# To verify that the system is nonlinear we will use the s_map function.
# Note that we will examine S-mapping in more detail in the afternoon
# session. Here though, we will use a small component of the s_map function
# to facilitate verification that our time series meet the requirements
# necessary for applying convergent cross mapping.

# A time series may show predictability even if it is completely stochastic
# because of processes like temporal autocorrelation.
# To ensure that the time series we are investigating are not completely 
# stochastic, we will fit a local linear map to reconstruct the state space.
# What does this mean?

# In the simplex projections we used above, the simplex function made its
# predictions using nearest neighbor interpolation.

# In contrast, we will now use the entire set of points in the state space,
# but with different weightings to create the local space around the point
# we are trying to predict. The weightings will be determined by
# a parameter called theta.

# Suppose we want to predict point q in a state space that is made up of 
# 3 dimensions. If theta = 0 we will weight all other points the in 
# state space equally, so they will all have an equivalent impact on 
# predicting q. This is the same as fitting an autoregressive model
# to the data.

# However, if theta > 0, nearby points receive larger weights, 
# allowing the local linear map to vary in state-space and accommodate 
# nonlinear behavior.

# Here, then, we want to see that values of theta > 0 result in higher
# predictive ability (rho) than does theta = 0.

# P2 time series
smap_output_P2 <- s_map(df$P2, lib, pred, E = E_P2)

# C2 time series
smap_output_C2 <- s_map(df$C2, lib, pred, E = E_C2)

# The result is a dataframe that we can use to examine how predictive
# ability (rho) changes as a function of the weighting given to nearby
# points for the prediction (theta)
plot(smap_output_P2$theta, smap_output_P2$rho, type = 'l', 
     xlab = "Theta (weighting parameter)", ylab = "Rho Value",
     main = "Nonlinearity Test for P2 time series")

plot(smap_output_C2$theta, smap_output_C2$rho, type = 'l', 
     xlab = "Theta (weighting parameter)", ylab = "Rho Value",
     main = "Nonlinearity Test for C2 time series")

# Both plots show that rho increases as theta increases and that the 
# lowest values of rho occur when theta = 0. 

# With this part of Step 2 we have confirmed that the system is in fact
# nonlinear.



###############################################
###                                         ###
###  Step 3: Calculate the two processes'   ###
###     ability to describe each other      ###
###############################################

# Now that we have determined the optimal embedding dimension for each 
# time series, and verified that each display characteristics consistent
# with nonlinearity, we can test for causality between the two series.

# Convergent cross mapping (CCM) is implemented via the ccm function.
# The function expects a single data frame or a matrix with named columns 
# for each time series. The first column must be a time vector. 

# To identify convergence, we compute cross-map skill (Pearson's 
# correlation, rho, between observed and predicted values) over 
# many random subsamples of the time series. These random subsamples
# help us see the confidence intervals of our rho estimates.

# ccm: Run CCM to test bidirectional causality
#   -  lib and pred again define the library (training) and prediction
#      (testing) data which are both specified as intervals. The default is to 
#      use the full time series (with leave-one-out cross validation)
#   -  E uses the optimal embedding dimensions from Step 1 (E_P2, E_C2)
#   -  tau is the time lag to use (defaults to 1)
#   -  tp is the prediction horizon (how far ahead to predict), here it
#      defaults to 0 which means we are predicting to the same time point.
#      So for example with tp = 0 we predict B(t) from A(t) (instantaneous
#      mapping). With tp = 1, we predict B(t+1) from A(t) (forward 
#      prediction).
#   -  norm specifies which method is used in calculating distance from
#      the prediction point to the nearest neighbors. Defaults to 2, which
#      corresponds to Euclidean distance. norm = 1 corresponds to Manhattan
#      distance.
#   -  num_neighbors is the number of nearest neighbors to use. The default
#      value will change based on which distance measure is specified with
#      norm.
#   -  lib_sizes specifies the sequence of subsampled library
#      sizes to test
#   -  num_samples controls how many random samples to take at each lib_size
#   -  random_libs and replace specify how the subsamples from lib_sizes
#      will be generated. Here, setting both to TRUE enables random 
#      sampling with replacement.

# Test P2 --> C2 causality
start.time <- Sys.time()
P2causesC2 <- ccm(
  df, 
  lib = c(1, 500),            # we will not use the full time series
  pred = c(1, 500),           #   to speed up computation time
  E = E_C2,                   # Embedding dimension of effect (C2)
  lib_sizes = seq(10,         # We vary lib size to confirm convergence
                  500, by = 5),  
  num_samples = 100,          # We will take 100 samples
  random_libs = TRUE,         # Randomly sample libraries (like bootstrapping)
  tp = 0,                     # Predict to the same time point
  lib_column = "P2",          # Hypothesized forcing process (P2) 
  target_column = "C2"        # Process being forced (C2)
)
end.time <- Sys.time()
end.time - start.time # Should take only about 20 seconds

# Note that we receive two warnings:
# The first serves as a critical reminder about CCM's unique causal inference
# logic.
# The second is the same warning we saw in Step 1, and is similarly alerting
# us to the fact that predictions were evaluated with leave-one-out cross
# validation since our time points for the two series overlapped.


# Test C2 --> P2 causality
start.time <- Sys.time()
C2causesP2 <- ccm(
  df, 
  lib = c(1, 500),            # we will not use the full time series
  pred = c(1, 500),           #   to speed up computation time
  E = E_P2,                   # Embedding dimension of effect (P2)
  lib_sizes = seq(10,         # We vary lib size to confirm convergence
                  500, by = 5),  
  num_samples = 100,          # We will take 100 samples
  random_libs = TRUE,         # Randomly sample libraries (like bootstrapping)
  tp = 0,                     # Predict to the same time point
  lib_column = "C2",          # Hypothesized forcing process (C2) 
  target_column = "P2"        # Process being forced (P2)
)
end.time <- Sys.time()
end.time - start.time # Should take only about 20 seconds


# The results for each test of causality include correlation coefficients (rho)
# for each lib_size.
# If causality exists we expect that rho will increase as the lib_size
# increases.

# Let's first examine the output provided by the ccm function:
str(P2causesC2)
str(C2causesP2)

# We will create a plot showing how rho changes with library length for
# both causality tests.
# To do this, we first add a type column which will make plotting with 
# ggplot easier
P2causesC2$Type <- "P2causesC2"
C2causesP2$Type <- "C2causesP2"

######################
##  Visualize the   ##
##     results      ##
######################

# 1st calculate the mean and standard deviation for each lib_size
# and combine the dataframes
combined.results <- rbind(
  P2causesC2 %>% 
    group_by(lib_size) %>% 
    summarise(rho.stddev = sd(rho), rho = mean(rho)) %>% 
    mutate(Type = "P2causesC2"),
  C2causesP2 %>% 
    group_by(lib_size) %>% 
    summarise(rho.stddev = sd(rho), rho = mean(rho)) %>% 
    mutate(Type = "C2causesP2")
)

# Create the plot
ccm.plot <- ggplot(data = combined.results, 
                   aes(x = lib_size, y = rho, 
                       ymin = rho - rho.stddev, 
                       ymax = rho + rho.stddev)) +
  geom_line(aes(color = Type)) +
  geom_ribbon(aes(ymin = rho - rho.stddev, 
                  ymax = rho + rho.stddev, 
                  fill = Type), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  theme_bw() +
  theme(axis.text.x = element_text(size=22, face="bold")) +
  theme(axis.text.y = element_text(size=22, face="bold")) +
  theme(axis.title.x = element_text(size=22, face="bold", 
                                    color="gray30")) +
  theme(axis.title.y = element_text(size=22, face="bold", 
                                    color="gray30")) +
  scale_x_continuous(name="Library Length") +
  scale_y_continuous(name="Rho", 
                     limits = c(-0.1, 1),
                     breaks = c(0.0, 0.2, 0.4, 
                                0.6, 0.8, 1.0)) +
  scale_colour_manual(name="", values = c('P2causesC2' = "blueviolet",
                                          'C2causesP2' = "orange1"),
                      labels = c("C2 causes P2",
                                 "P2 causes C2")) +
  scale_fill_manual(name="", values = c('P2causesC2' = "blueviolet",
                                        'C2causesP2' = "orange1"),
                    labels = c("C2 causes P2",
                               "P2 causes C2")) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.2)) +
  # annotate("text", x = 3000, y = 0.32, label = )
  # theme(legend.position = "none") +
  theme(legend.text = element_text(size=16)) +
  guides(colour = guide_legend(override.aes = (list(size=3))))

# View the plot
ccm.plot

# Our results preliminarily indicate bidirectional causality between 
# abundances of Predator 2 and abundances of Consumer 2. 
# Rho values for each increase with increasing library length and converge.

# We can see though that the amount of information transferred (i.e,
# the magnitude of encoding from the causal time series to the effect
# time series) is higher for Predator 2 abundances --> 
# Consumer 2 abundances. 


###############################################
###                                         ###
###    Step 4: Test for significance        ###
###                                         ###
###############################################

# There are two commonly accepted methods for determining whether the 
# causal relationship we identified in Step 3 is significant.

### Method 1: Generating a null distribution of rho values and compare the
###            rho value at the maximal library length from the actual time
###            series against this null distribution.



###########################
###                     ###   
###      P2 --> C2      ###
###                     ###
###########################

# The idea is that we should not be able to predict the time series of P2 
# (the cause) from the time series of the effect (C2) if the values of 
# P2 are shuffled. If we can predict in this case, it would indicate that 
# there is not sufficient evidence of a causal process between P2 --> C2.

# Create matrix of surrogate time series of the hypothesized causal force
surrogate_series <- make_surrogate_data(
  df$P2[1:500],               # the raw time series of the hypothesized 
                              #      causal force
  method = "random_shuffle",  # method for creating the random samples 
  num_surr = 100              # the number of null time series to generate
)

# Convert to a dataframe
surrogate_series <- as.data.frame(surrogate_series)

# Add a time column
surrogate_series$Time <- seq(1, 500, by = 1)

# Add the C2 (hypothesized effect)
surrogate_series$C2 <- df$C2[1:500]

# Reorder the columns so that Time is first
surrogate_series <- surrogate_series[ , c( 1, 2:101, 102)]


# We will repeat steps 1 and 3 for each surrogate time series, now looking at 
# the causal process of surrogate P2 --> C2, but we will automate it to 
# make it more computationally efficient. Note we skip step 2 in the CCM
# process as we already confirmed the two complete time series met the
# requirements.

# Define a surrogate index based on the column numbers
sur.index <- c(2:101)

# Find the best embedding dimension of the effect (C2) as this is what
# is required for the CCM where we test P2 --> C2

# Define maxE
maxE <- 15

# Initialize a blank vector
rho_surrC2 <- numeric(maxE)

# Calculate rho for each possible #
for (E in 1:maxE) {
simplex_output <- simplex(surrogate_series$C2,
                             E = E, tp = 1)
rho_surrC2[E] <- simplex_output$rho
}

# Extract the E where rho is at its max
E_surrC2 <- which.max(rho_surrC2)

# Create a blank vector to store rho values from the surrogate ccm tests
rho.surr <- numeric(length(sur.index))

# Run the CCM for each surrogate time series
start.time <- Sys.time() # Takes about 15 seconds
for(i in 1:length(sur.index)) {

  # Run the ccm
  test <- ccm(
    surrogate_series, 
    lib = c(1, 500),            # we use the full time series
    pred = c(1, 500),          
    E = E_C2,                   # Embedding dimension of effect (C2)
    lib_sizes = seq(400,         # Now we are only concerned with max lib size
                    500, by = 100),  
    num_samples = 50,           # We will take 50 samples to speed up computation
    random_libs = TRUE,         # Randomly sample libraries (like bootstrapping)
    tp = 0,                     # Predict to the same time point
    lib_column = sur.index[i],             # Surrogate series of hypothesized cause (P2) 
    target_column = "C2"        # Process being forced (C2)
  )
  
  # Calculate the mean value of rho at the max library length and add it to 
  # the results vector for the rho from the surrogate tests
  rho.surr[i] <- mean(test$rho[test$lib_size == max(test$lib_size)])
  
  print(i)
}
end.time <- Sys.time()
end.time - start.time

# Assess the significance using the empirical cumulative distribution function
# 1 - ecdf(null distribution column)(actual rho)
rho.signif <- 1 - 
  ecdf(rho.surr)(combined.results$rho[combined.results$lib_size ==
                                        max(combined.results$lib_size) &
                                        combined.results$Type == "P2causesC2"])

# Create dataframe of results from the true time series to facilitate plotting,
# Adding a column for the significance value
true.results <- combined.results %>%
  subset(Type == "P2causesC2") %>%
  filter(lib_size == max(lib_size)) %>%
  mutate(Significance = ifelse(rho.signif < 0.05, "significant", "non"),
         Type = "P2causesC2")

# Make a dataframe of surrogate results
rho.surr.df <- data.frame(rho.surr = rho.surr, Type = "P2causesC2")


###########################
###                     ###   
###      C2 --> P2      ###
###                     ###
###########################

# The idea is that we should not be able to predict the time series of C2 
# (the cause) from the time series of the effect (P2) if the values of 
# C2 are shuffled. If we can predict in this case, it would indicate that 
# there is not sufficient evidence of a causal process between C2 --> P2.

# Create matrix of surrogate time series of the hypothesized causal force
surrogate_series <- make_surrogate_data(
  df$C2[1:500],               # the raw time series of the hypothesized 
  #      causal force
  method = "random_shuffle",  # method for creating the random samples 
  num_surr = 100              # the number of null time series to generate
)

# Convert to a dataframe
surrogate_series <- as.data.frame(surrogate_series)

# Add a time column
surrogate_series$Time <- seq(1, 500, by = 1)

# Add the P2 (hypothesized effect)
surrogate_series$P2 <- df$P2[1:500]

# Reorder the columns so that Time is first
surrogate_series <- surrogate_series[ , c( 1, 2:101, 102)]


# We will repeat steps 1 and 3 for each surrogate time series, now looking at 
# the causal process of surrogate C2 --> P2, but we will automate it to 
# make it more computationally efficient. Note we skip step 2 in the CCM
# process as we already confirmed the two complete time series met the
# requirements.

# Define a surrogate index based on the column numbers
sur.index <- c(2:101)

# Find the best embedding dimension of the effect (P2) as this is what
# is required for the CCM where we test C2 --> P2

# Define maxE
maxE <- 15

# Initialize a blank vector
rho_surrP2 <- numeric(maxE)

# Calculate rho for each possible #
for (E in 1:maxE) {
  simplex_output <- simplex(surrogate_series$P2,
                            E = E, tp = 1)
  rho_surrP2[E] <- simplex_output$rho
}

# Extract the E where rho is at its max
E_surrP2 <- which.max(rho_surrP2)

# Create a blank vector to store rho values from the surrogate ccm tests
rho.surr <- numeric(length(sur.index))

# Run the CCM for each surrogate time series
start.time <- Sys.time() # Takes about 15 seconds
for(i in 1:length(sur.index)) {
  
  # Run the ccm
  test <- ccm(
    surrogate_series, 
    lib = c(1, 500),            # we use the full time series
    pred = c(1, 500),          
    E = E_P2,                   # Embedding dimension of effect (P2)
    lib_sizes = seq(400,         # Now we are only concerned with max lib size
                    500, by = 100),  
    num_samples = 50,           # We will take 50 samples to speed up computation
    random_libs = TRUE,         # Randomly sample libraries (like bootstrapping)
    tp = 0,                     # Predict to the same time point
    lib_column = sur.index[i],             # Surrogate series of hypothesized cause (P2) 
    target_column = "P2"        # Process being forced (C2)
  )
  
  # Calculate the mean value of rho at the max library length and add it to 
  # the results vector for the rho from the surrogate tests
  rho.surr[i] <- mean(test$rho[test$lib_size == max(test$lib_size)])
  
  print(i)
}
end.time <- Sys.time()
end.time - start.time

# Assess the significance using the empirical cumulative distribution function
# 1 - ecdf(null distribution column)(actual rho)
rho.signif <- 1 - 
  ecdf(rho.surr)(combined.results$rho[combined.results$lib_size ==
                                        max(combined.results$lib_size) &
                                        combined.results$Type == "C2causesP2"])

# Add these results to the true results dataframe
true.results <- rbind(true.results,
                      combined.results %>%
                        subset(Type == "C2causesP2") %>%
                        filter(lib_size == max(lib_size)) %>%
                        mutate(Significance = ifelse(rho.signif < 0.05,
                                                     "significant", "non"),
                               Type = "C2causesP2"))

# Add the surrogate results to the dataframe of surrogate results
rho.surr.df <- rbind(rho.surr.df,
                     data.frame(rho.surr = rho.surr, Type = "C2causesP2"))
  

###########################
###                     ###   
###  Graph the results  ###
###                     ###
###########################

significance.plot <- ggplot() +
  geom_density_ridges(data = rho.surr.df,
                      aes(x = rho.surr, y = factor(Type),
                          fill = factor(Type), scale = 0.7,
                      rel_min_height = 0.05),
                      color = "black") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  # True values -  red if significant, gray if not
  geom_point(data = true.results,
             aes(x = rho, y = factor(Type),
                 color = ifelse(Significance == "Significant", "red", "gray50")),
             size = 5, shape = 19) +
  scale_x_continuous(name = "Rho Values", lim = c(-0.20, 0.90),
                     breaks = c(-0.2, 0.0, 0.2, 0.4, 0.6, 0.8),
                     labels = c(-0.2, 0.0, 0.2, 0.4, 0.6, 0.8)) + 
  scale_y_discrete(name = "", labels = c("C2causesP2" = "C2 --> P2",
                                         "P2causesC2" = "P2 --> C2")) +
  scale_fill_manual(values = c("P2causesC2" = "#1f77b4",
                               "C2causesP2" = "#ff7f0e")) +
  theme_minimal() + 
  theme(axis.text.x = element_text(size=22, face="bold")) +
  theme(axis.text.y = element_text(size=22, face="bold")) +
  theme(axis.title.x = element_text(size=22, face="bold", color="gray30")) +
  theme(axis.title.y = element_text(size=22, face="bold", color="gray30")) +
  theme(legend.position = "none")


# We see that for both P2 --> C2 and C2 --> P2, the rho value obtained from the
# true time series far exceeds the rho values obtained from surrogate time series.
# Therefore we can say that there is bidirectional causality between P2 and C2.



### Method 2: Using both Kendall's tau test AND Fisher's Z test.

# This approach has two components. 
# 1. Tests for a significant monotonically (i.e., continuously) increasing
#      trend in rho using the Mann Kendall test. 
# 2. Use Fisher's exact test to test whether the rho value obtained under
#      the maximal library length is significantly higher than the rho value
#      obtained using the minimal library length.

# If both tests are significant then the causal process is deemed significant.

# We will first assess significance for P2 causing C2.
MKP2as.cause <- MannKendall(combined.results$rho[combined.results$Type ==
                                                   "P2causesC2"])
# Inspect the results
MKP2as.cause

# The positive value of tau and the p-value < 0.05 indicate that we have a 
# significant monotonically increasing trend.

# Now we will apply Fisher's Z test to compare the rho values. Note that
# our correlation values are not independent and therefore the version of
# the Fisher's Z test that we use must account for this. Our correlations
# are also overlapping meaning that they have a variable in common.

# Divide combined results by type
P2cC2.results <- combined.results %>%
  filter(Type == "P2causesC2" )


# Note that for overlapping correlations between dependent groups we need
# to specify the correlation between our correlations which sounds very 
# confusing! This is called the intercorrelation coefficient. 
# We don't know this value as each lib_size represents the aggregate results
# of 100 random samples.
# Therefore we will use a conservative estimate of 0.5, indicating that the 
# two correlations might be correlated but they also might not.
FTP2as.cause <- cocor::cocor.dep.groups.overlap(
  r.jk = P2cC2.results$rho[ which.min( P2cC2.results$lib_size ) ], # Cor 1
  r.jh = P2cC2.results$rho[ which.max( P2cC2.results$lib_size ) ], # Cor 2
  r.kh = 0.5,                           # correlation btw the correlations 
  n = nrow(P2cC2.results),              # the 
  var.labels = c("Min Lib", "Max Lib", "Intercorrelation")
)

# Examine the results for the Fisher test which is implemented with the 
# pearson1898 test
FTP2as.cause@pearson1898

# The p-value below 0.05 indicates that the rho value at the maximum library
# length is significantly different than the rho value at the minimum library
# length.

### As a result we can say that there is significant evidence of a causal
### process where P2 affects C2.



# We will do the same for C2 causing P2.
MKC2as.cause <- MannKendall(combined.results$rho[combined.results$Type ==
                                                   "C2causesP2"])

# Inspect the results
MKC2as.cause

# The negative value of tau and the p-value < 0.05 indicate that we have a 
# significant monotonically BUT decreasing trend.

# As we needed a significant increasing trend, we already now that
# C2 --> P2 won't meet the requirement for determining significance.

# But for thoroughness we will still apply Fisher's Z test to compare the 
# rho values. Note that
# our correlation values are not independent and therefore the version of
# the Fisher's Z test that we use must account for this. Our correlations
# are also overlapping meaning that they have a variable in common.

# Divide combined results by type
C2cP2.results <- combined.results %>%
  filter(Type == "C2causesP2" )

# Note that for overlapping correlations between dependent groups we need
# to specify the correlation between our correlations which sounds very 
# confusing! This is called the intercorrelation coefficient. 
# We don't know this value as each lib_size represents the aggregate results
# of 100 random samples.
# Therefore we will use a conservative estimate of 0.5, indicating that the 
# two correlations might be correlated but they also might not.
FTC2as.cause <- cocor::cocor.dep.groups.overlap(
  r.jk = C2cP2.results$rho[ which.min( P2cC2.results$lib_size ) ], # Cor 1
  r.jh = C2cP2.results$rho[ which.max( P2cC2.results$lib_size ) ], # Cor 2
  r.kh = 0.5,                           # correlation btw the correlations 
  n = nrow(C2cP2.results),              # the 
  var.labels = c("Min Lib", "Max Lib", "Intercorrelation")
)

# Examine the results for the Fisher test which is implemented with the 
# pearson1898 test
FTC2as.cause@pearson1898

# The p-value below 0.05 indicates that the rho value at the maximum library
# length is significantly different than the rho value at the minimum library
# length.

### However, since C2 --> P2 did not pass the Mann Kendall test, using this
### method of evaluation (Mann Kendall + Fisher) we can not say that there is 
### significant evidence of a causal process where C2 affects P2.



################################################
###                                          ###
###  Exercise: Try to replicate the analysis ###
###    we did choosing two new components    ###
###     of the system (R, C1, C2, P1, P2)    ###
###                                          ###
################################################
