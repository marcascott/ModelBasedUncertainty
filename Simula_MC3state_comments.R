# =============================================================================
# Simula_MC3state.R
# =============================================================================
#
# Purpose: Simulate longitudinal data from a first-order Markov Chain (MC)
#          model with 3 latent states and 2 time-varying covariates affecting
#          both initial state probabilities and transition probabilities using 
#          LMest package
#          Pennoni F., Pandolfi S. and Bartolucci F. (2025). LMest: An R Package for
#          Estimating Generalized Latent Markov Models. The R Journal, 16(4), 74--101.
#          doi:10.32614/RJ-2024-036
#
# MODEL OVERVIEW:
#   A first-order Markov Chain assumes that the probability of being in a
#   state at time t depends ONLY on the state at time t-1 (the Markov property).
#   Here the model is extended to include covariates:
#     - Initial probabilities Piv(i):   P(state at t=1 | covariates at t=1)
#     - Transition probabilities PI(i,t): P(state at t | state at t-1, covariates)
#
#   The simulation proceeds in two stages:
#     1. Specify parameters (la for initial, psi for transitions)
#     2. Draw state sequences for n individuals over TT time points
#
# DEPENDENCIES:
#   - MASS       : multivariate normal simulation (not used directly but required)
#   - LMest      : Latent Markov model estimation; also provides lmestData(),
#                  matrices2long(), and lmestMc()
#   - draw_mc.R        : custom function to draw state sequences from a MC
#   - comp_Piv.R       : custom function to compute initial state probabilities
#   - comp_PI.R        : custom function to compute transition probability matrices
#   - design_matrices1.R : custom function to build design matrices for the model
#
# OUTPUT:
#   - Simulated dataset (dat1 / dat31_3) with columns: id, time, y1, x1, x2
#   - Empirical transition matrix computed two ways (manual loop and dplyr)
#   - Fitted MC model via lmestMc() for verification
# =============================================================================


# -----------------------------------------------------------------------------
# ENVIRONMENT SETUP
# -----------------------------------------------------------------------------

# Clear the entire R workspace to ensure a clean environment.
# This avoids conflicts with objects left over from previous sessions.
rm(list = ls())

# Load required packages.
# require() is like library() but returns FALSE instead of an error if missing.
require(MASS)    # General-purpose package; provides mvrnorm() and other utilities.
                 # Although not called explicitly here, it is a dependency of LMest.
require(LMest)   # Latent Markov Estimation package.
                 # Provides: lmestData(), matrices2long(), lmestMc() and other
                 # tools for fitting and simulating Markov and Hidden Markov models.

# Load custom helper functions from local source files.
# These scripts must be in the same working directory as this file.
source("functions.R")    # required functions    


# =============================================================================
# SECTION 1: SIMULATION SETTINGS
# =============================================================================

k  <- 3     # Number of latent states in the Markov Chain.
            # States are labelled 1, 2, 3 (e.g., low / medium / high risk).

n  <- 1000  # Number of subjects (individuals) in the simulated panel.

TT <- 31    # Number of time periods observed for each individual.
            # Total observations = n × TT = 31,000 rows.

nc <- 2     # Number of time-varying covariates (x1 and x2).

# model_cov and model_int control which covariates enter which part of the model.
# "all" means ALL covariates affect BOTH initial AND transition probabilities.
# Alternative options (defined in design_matrices1.R) 

model_cov <- "all"   # Covariates affecting TRANSITION probabilities
model_int <- "all"   # Covariates affecting INITIAL state probabilities


# =============================================================================
# SECTION 2: DESIGN MATRICES
# =============================================================================
# The multinomial logit model for a k-state Markov Chain requires design matrices
# that map covariate values to log-odds of state membership.
#
# design_matrices1() returns four matrices:
#
#   G  : design matrix for INITIAL state probabilities (k-1 free parameters per covariate)
#        Rows correspond to the (k-1) free log-odds equations for the initial distribution.
#        The baseline state is dropped for identifiability (here: "central" = state 2).
#
#   Z  : design matrix for TRANSITION probabilities.
#        Encodes all k*(k-1) free log-odds equations for transitions out of each state.
#
#   GG : expanded version of G used internally by comp_Piv().
#
#   DD : difference / contrast matrix used internally for transition parameterisation.
#
# baseline = "central": the MIDDLE state (state 2 of 3) is the reference category.
# This is a modelling choice — changing the baseline changes only coefficient interpretation

out <- design_matrices1(k, nc,
                        baseline   = "central",
                        model_cov  = model_cov,
                        model_int  = model_int)

# Unpack the four design matrices from the output list

G  <- out$G    # Design matrix for initial probabilities
Z  <- out$Z    # Design matrix for transition probabilities
GG <- out$GG   # Expanded initial design matrix
DD <- out$DD   # Contrast matrix for transitions

# Inspect Z dimensions: should be [k*(k-1) * (nc+1)] × [k*(k-1) * (nc+1)]
# The +1 accounts for the intercept column in each logit equation.
dim(Z)

# =============================================================================
# SECTION 3: PARAMETER SPECIFICATION — INITIAL STATE PROBABILITIES
# =============================================================================
# la contains the log-odds parameters for the multinomial logit model of the
# INITIAL state distribution (probability of starting in each state at t=1).
#
# The vector has length (nc+1) * (k-1):
#   - (k-1) = 2 free equations (state 1 vs. baseline, state 3 vs. baseline)
#   - (nc+1) = 3 coefficients per equation (intercept + β_x1 + β_x2)
#
# Values: c(0, 0.5, 1) repeated to fill the vector.
#   - 0   → intercept: equal baseline probability (no preference for any state)
#   - 0.5 → moderate positive effect of first covariate on log-odds of that state
#   - 1   → stronger positive effect of second covariate

la <- array(c(0, 0.5, 1), (nc + 1) * (k - 1))
la   # Print to verify: should show 6 values: 0, 0.5, 1, 0, 0.5, 1


# =============================================================================
# SECTION 4: PARAMETER SPECIFICATION — TRANSITION PROBABILITIES
# =============================================================================
# psi contains the log-odds parameters for the multinomial logit model of the
# TRANSITION probabilities: P(state at t | state at t-1, covariates).
#
# For a k=3 state chain, there are k*(k-1) = 6 free log-odds equations
# (for each of the k=3 origin states, k-1=2 destination states are free).
# Each equation has (nc+1) = 3 coefficients (intercept + 2 covariate effects).
#
#
# STEP 1: Generate random intercepts for each free equation.
# These intercepts are drawn to produce realistic (non-extreme) transition probabilities.
# The formula -log(0.9/0.1 + runif(1)) generates values that
# correspond to transition probabilities in a plausible range.

set.seed(5432)   # Fix the random seed for reproducibility.
                 # Using a fixed seed ensures the same psi is generated every run.

psi <- rep(0, k * (k - 1))   # Initialise psi with zeros (length = 6)

# Fill with random intercepts: one per free log-odds equation.
for (it in 1:(k * (k - 1))) psi[it] <- -log(0.9 / 0.1 + runif(1))
# -log(0.9/0.1 + U) where U ~ Uniform(0,1):
#   0.9/0.1 = 9, so argument to log is in [9, 10]
#   → log values in [log(9), log(10)] ≈ [2.20, 2.30]
#   → negated: approximately [-2.30, -2.20]
# These negative intercepts translate to small baseline transition probabilities
# (i.e., staying in the current state is more likely than switching).

# STEP 2: Define covariate coefficients for transitions.
# bet1 and bet2 are the effects of x1 and x2 on the transition log-odds.
# Both sequences range from positive values, indicating that higher covariate
# values increase the log-odds of transitioning to non-baseline states.
bet1 <- seq(.2, .8, length = 6)   # 6 values from 0.2 to 0.8 (effect of x1)
bet2 <- seq(.7, 1.3, length = 6)  # 6 values from 0.7 to 1.3 (effect of x2)

# STEP 3: Assemble the full psi vector by interleaving intercepts and slopes.
psi <- as.vector(cbind(psi, bet1, bet2))


# =============================================================================
# SECTION 5: SIMULATING TIME-VARYING COVARIATES
# =============================================================================
# Covariates x and y (renamed x1, x2 later) follow AR(1) processes:
#   x[i,t] = 0.5 * x[i,t-1] + ε[i,t],   ε ~ N(0,1)
#   y[i,t] = 0.5 * y[i,t-1] + ε[i,t],   ε ~ N(0,1)
#
# An AR(1) process with autoregressive parameter φ = 0.5:
#   - Mimics realistic temporal dependence in panel data covariates

b <- 30
set.seed(b + 432)   # Separate seed for covariate simulation 
                    # Using a different seed from psi keeps covariate and parameter
                    # simulations independent and reproducible.

# Initialise n × TT matrices for the two covariates (all zeros).
x <- matrix(0, n, TT)
y <- matrix(0, n, TT)

# Draw initial values (t=1) from a standard normal distribution.
# These serve as the starting point of each AR(1) series.
x[, 1] <- rnorm(n)
y[, 1] <- rnorm(n)

# Simulate AR(1) trajectories for t = 2, ..., TT.
# Outer loop: over individuals i (independent AR processes per subject)
# Inner loop: over time points t (sequential dependence within subject)
for (i in 1:n) {
  for (t in 2:TT) {
    x[i, t] <- 0.5 * x[i, t - 1] + rnorm(1)   # AR(1) for covariate 1
    y[i, t] <- 0.5 * y[i, t - 1] + rnorm(1)   # AR(1) for covariate 2
  }
}

# Stack both covariates into a 3D array: XX[i, t, j]
#   - Dimension 1 (i): individual index
#   - Dimension 2 (t): time index
#   - Dimension 3 (j): covariate index (1 = x, 2 = y)
XX <- array(0, c(n, TT, 2))
XX[, , 1] <- as.matrix(x)
XX[, , 2] <- as.matrix(y)


# =============================================================================
# SECTION 6: BUILDING THE MODEL MATRIX
# =============================================================================
# Convert the 3D covariate array to long format and construct a standard
# design matrix (with intercept) as required by comp_Piv() and comp_PI().

# matrices2long() (from LMest) reshapes the n × TT × nc array into a long
# data frame with columns: id, time, Y1, Y2 (one row per individual-timepoint).
res  <- matrices2long(Y = XX)

# Assemble a standard data frame with id, time, and covariate columns.
Xcov <- data.frame(cbind(
  "id"   = res$id,
  "time" = res$time,
  "X1"   = res$Y1,   # Covariate 1 (x) in long format
  "X2"   = res$Y2    # Covariate 2 (y) in long format
))

# model.matrix() creates the full design matrix including the intercept column.
# Formula ~ X1 + X2 → columns: (Intercept), X1, X2
XX  <- model.matrix(as.formula(~ X1 + X2), data = Xcov)

# Extract the time index vector (one entry per row of XX).
tv  <- Xcov[, 2]

# Select only the rows corresponding to t = 1 (the first time point).
# This sub-matrix is used exclusively for the INITIAL state model,
# since initial probabilities are conditioned on covariates at t=1 only.
XX1 <- XX[tv == 1, , drop = FALSE]  


# =============================================================================
# SECTION 7: COMPUTING INITIAL STATE PROBABILITIES
# =============================================================================
# comp_Piv() applies the multinomial logit model for the initial distribution.
# It takes the covariate matrix at t=1 (XX1), the design matrix G, and the
# parameter vector la, and returns the n × k matrix of initial probabilities.
#
# Piv[i, j] = P(state_1 = j | covariates of individual i at t=1)
# Each row sums to 1 (valid probability distribution over k states).
#
# fort = FALSE: use the pure R implementation (not Fortran-compiled code).
# Set fort = TRUE for faster computation with large datasets.

Piv <- comp_Piv(n, k, XX1, G, la, fort = FALSE)


# =============================================================================
# SECTION 8: COMPUTING TRANSITION PROBABILITY MATRICES
# =============================================================================
# comp_PI() computes the transition probabilities for all individuals and
# all time transitions (t-1 → t, for t = 2, ..., TT).
#
# The output PI is an array of dimension n × (TT-1) × k × k:
#   PI[i, t, j, l] = P(state at t+1 = l | state at t = j, covariates of i at t)
#
# Each k × k slice (for fixed i and t) is a valid stochastic matrix:
#   - All entries in [0, 1]
#   - Each row sums to 1
#
# eta = Z %*% psi: the linear predictor vector for all transition log-odds.
#   Z is the design matrix (built in Section 2).
#   psi is the parameter vector (built in Section 4).
#   Their product gives the log-odds of each transition for all individuals/times.


# Compute the linear predictor for transition log-odds.
eta <- Z %*% psi   # Matrix-vector product: [dim(Z)[1] × 1] vector

# Compute the full transition probability array.
PI <- comp_PI(k, n, TT, XX, GG, eta, fort = FALSE)


# =============================================================================
# SECTION 9: DRAWING STATE SEQUENCES (SIMULATION)
# =============================================================================
# draw_mc() uses the initial probabilities (Piv) and transition matrices (PI)
# to simulate a state sequence for each individual over TT time points.
#
# Algorithm (per individual i):
#   1. Draw state at t=1 from Multinomial(1, Piv[i, ])
#   2. For t = 2, ..., TT: draw state at t from Multinomial(1, PI[i, t-1, prev_state, ])
#
# Returns a list with:
#   sim$Y: a data frame in long format with columns id, time, y1 (the simulated state)

sim <- draw_mc(Piv, PI)

# Assemble the final simulated dataset:
#   - sim$Y:        id, time, y1 (simulated state labels 1, 2, or 3)
#   - Xcov[, 3:4]:  x1 and x2 covariate values (columns 3 and 4 of Xcov)
dat1 <- cbind(sim$Y, Xcov[, 3:(2 + nc)])

# Assign meaningful column names for downstream analysis and export.
names(dat1) <- c("id",    # Individual identifier (1 to n)
                 "time",  # Time point (1 to TT)
                 "y1",    # Simulated latent state (1, 2, or 3)
                 "x1",    # Covariate 1 (AR(1) process)
                 "x2")    # Covariate 2 (AR(1) process)
names(dat1)   # Print column names to verify


# =============================================================================
# SECTION 10: PREPARING DATA FOR LMEST ESTIMATION AND DESCRIPTIVE SUMMARY
# =============================================================================
# lmestData() (from LMest) converts a flat data frame into the structured
# format required by LMest's estimation functions.
#
# Arguments:
#   data             : the long-format data frame
#   id               : name of the individual identifier column
#   time             : name of the time column
#   responsesFormula : formula specifying the response (y1) and covariates (x1, x2)
#                      Here y1 is the OBSERVED discrete response (the simulated state).
#                      x1 and x2 are the covariates included in the MC model.

dt <- lmestData(data             = dat1,
                id               = "id",
                time             = "time",
                responsesFormula = y1 ~ x1 + x2)

# summary() with dataSummary = "responses" prints descriptive statistics
# for the response variable y1.
# varType = rep("d", ncol(dt$Y)): specifies that all response variables are
# DISCRETE (categorical), which affects the type of summary shown.
summary(dt,
        dataSummary = "responses",
        varType     = rep("d", ncol(dt$Y)))


# Save the simulated dataset under a versioned name for external use.
# Naming convention: dat[TT]_[k] → dat31_3 = 31 time points, 3 states.
dat31_3 <- dat1

# save(dat31_3, file = "dat31_3.Rdata")
