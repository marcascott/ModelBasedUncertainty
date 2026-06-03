### Functions implemented in the paper:
# Bartolucci, F.; Pandolfi, S.; Pennoni, F. Parsimonious parametrizations of transition matrices850
# of Markov chain and hidden Markov models. Annals of Operations Research 2026. 
# https://link.springer.com/article/10.1007/s10479-025-06986-x
#
# =============================================================================
# design_matrices1.R
# =============================================================================
# Purpose: Construct the design matrices needed to parameterise the
#          multinomial logit models for both the initial state distribution
#          and the transition probabilities of a covariate-dependent
#          first-order Markov Chain.
#
#
#   BASELINE STATE:
#     The multinomial logit model is not identified without fixing one
#     reference category per origin state. Here, 'baseline = "central"'
#     sets the MIDDLE state (state u itself, the diagonal of the transition
#     matrix) as the reference. This is sometimes called the "diagonal
#     reference" parameterisation and is natural for ordinal-type states
#     where staying in the same state is the dominant event.
#
# ARGUMENTS:
#   k          : integer — number of latent states (e.g., 3)
#   nc         : integer — number of time-varying covariates (e.g., 2)
#   baseline   : character — reference category for the multinomial logit.
#                Currently only "central" is implemented, which sets the
#                baseline destination to state u (the origin state itself)
#                for each row u of the transition matrix.
#   model_int  : character — controls which origin states have free intercepts
#                in the transition model. "all" = all k*(k-1) intercepts are
#                free (no equality constraints across origin states).
#   model_cov  : character — controls which origin states share covariate
#                slopes in the transition model. "all" = all k*(k-1) slope
#                vectors are free (separate covariate effects per origin state).
#
# RETURNS:
#   A named list with the following matrices/arrays:
#
#   G  : matrix, dim = k × (k-1)
#        Contrast matrix for the INITIAL distribution.
#        Maps (k-1) free log-odds to all k states by dropping column 1
#        of the k×k identity (the first state is the baseline for G).
#        Used in comp_Piv() as: G %*% Xi %*% la → k-length log-odds vector.
#
#   GG : array, dim = k × (k-1) × k
#        Contrast matrices for the TRANSITION probabilities, one per origin state.
#        GG[,,u]: maps the (k-1) free destination log-odds (for transitions
#        OUT OF state u) to all k states, with state u as the baseline.
#        Used in comp_PI() as: GG[,,u] %*% Xi %*% eta[ind] → k-length log-odds.
#
#   Z  : matrix, dim = k*(k-1) × k*(k-1)*(1+nc)
#        Full design matrix for the TRANSITION model.
#        Encodes the block structure of intercepts (Z1) and covariate slopes (Z2),
#        allowing different parameters per origin state and per covariate.
#
#
#   Z1 : matrix, dim = k*(k-1) × k*(k-1)
#        Design matrix block for INTERCEPTS only.
#        Under model_int = "all": Z1 = I_{k*(k-1)} (identity matrix),
#        meaning each of the k*(k-1) free equations has its own intercept.
#
#   Z2 : matrix, dim = k*(k-1) × k*(k-1)   [only returned when nc > 0]
#        Design matrix block for COVARIATE SLOPES only.
#        Under model_cov = "all": Z2 = I_{k*(k-1)},
#        meaning each of the k*(k-1) free equations has its own slope vector.
# =============================================================================


design_matrices1 <- function(k, nc, baseline, model_int, model_cov) {
  knc <- (k - 1) * (1 + nc)
  
  
  # ===========================================================================
  # MATRIX G — Contrast matrix for the INITIAL distribution
  # ===========================================================================
  # The initial distribution multinomial logit has k states but only (k-1)
  # free log-odds equations (relative to a baseline state).
  #
  #
  # INTERPRETATION:
  # State 1 is the baseline for the initial distribution (its log-odds = 0).
  # G maps a (k-1)-length vector of free log-odds (for states 2, ..., k)
  # to a k-length vector by prepending a 0 for state 1
  #
  # ---------------------------------------------------------------------------
  G <- diag(k)[, -1]
  
  
  # ===========================================================================
  # ARRAY GG — Contrast matrices for the TRANSITION probabilities
  # ===========================================================================
  GG <- array(0, c(k, k - 1, k))
  for (u in 1:k) {
    # -------------------------------------------------------------------------
    # BASELINE = "CENTRAL":
    # For each origin state u, the baseline destination is state u ITSELF
    # (the diagonal element of the transition matrix — "staying in the same
    # state"). This is the "central" or "diagonal reference" parameterisation.
    #
    # -------------------------------------------------------------------------
    if (baseline == "central") GG[, , u] <- diag(k)[, -u]
    
  }
  
  
  # ===========================================================================
  # MATRIX Z1 — Intercept block of the transition design matrix
  # ===========================================================================
  if (model_int == "all") Z1 <- diag(k * (k - 1))
  
  # ===========================================================================
  # FULL DESIGN MATRIX Z — Combines intercepts (Z1) and covariate slopes (Z2)
  # ===========================================================================
  if (nc == 0) Z <- Z1
  
  # ---------------------------------------------------------------------------
  # CASE: WITH COVARIATES (nc > 0)
  # ---------------------------------------------------------------------------
  if (nc > 0) {
    if (model_cov == "all") Z2 <- diag(k * (k - 1))
    # -------------------------------------------------------------------------
    Z <- cbind(Z1 %x% diag(1 + nc)[, 1],      # Intercept selector block
               Z2 %x% diag(1 + nc)[, -1])      # Covariate slope selector block
   out <- list(G  = G,    # k × (k-1): contrast matrix for initial distribution
              Z  = Z,    # k*(k-1) × k*(k-1)*(1+nc): full transition design matrix
              GG = GG,   # k × (k-1) × k: per-origin contrast array for transitions
              Z1 = Z1)   # k*(k-1) × k*(k-1): intercept block
  
  if (nc > 0) out$Z2 <- Z2   # k*(k-1) × k*(k-1): covariate slope block (if applicable)
  
  return(out)
  
}
}
# =============================================================================
# comp_PI.R
# =============================================================================
# Purpose: Compute the array of TRANSITION PROBABILITY MATRICES for a
#           first-order Markov Chain model with covariates
#
#
#   The reference (baseline) destination state is absorbed into the intercept:
#   the log-odds of transitioning to the baseline state is set to 0, so
#   only k-1 free equations are needed per origin state u.
#
#
# ARGUMENTS:
#   k    : integer — number of latent states
#   n    : integer — number of individuals
#   TT   : integer — number of time points
#   XX   : matrix  — design matrix in long format, dim = (n*TT) × (nc+1)
#                    rows are ordered by (individual, time): i=1,t=1; i=1,t=2; ...
#                    columns are: (Intercept), X1, X2, ..., Xnc
#   GG   : array   — contrast/expansion array of dim (k-1) × (k*(k-1)*(1+nc)/k) × k
#                    encodes the multinomial logit contrast structure for each
#                    origin state u (third dimension); produced by design_matrices1()
#   eta  : vector  — linear predictor vector of length k*(k-1)*(1+nc)
#                    computed upstream as Z %*% psi, where Z is the global
#                    transition design matrix and psi is the parameter vector
#   fort : logical — if TRUE, calls a compiled Fortran routine for speed;
#                    if FALSE (default), uses the pure R implementation below
#
# RETURNS:
#   PI   : array of dim k × k × n × TT containing transition probabilities
#
# =============================================================================


comp_PI <- function(k, n, TT, XX, GG, eta, fort = FALSE) {
  
  # ---------------------------------------------------------------------------
  # Derive the number of covariates from the design matrix.
  # XX has (nc + 1) columns: one intercept + nc covariate columns.
  # Subtracting 1 recovers nc (the number of substantive covariates).
  # ---------------------------------------------------------------------------
  nc <- ncol(XX) - 1
  #
  # ===========================================================================
  # BRANCH 1:
  # ===========================================================================
  # Calls a pre-compiled Fortran subroutine "comp_PI" for performance.
  if (fort) {
    
    out <- .Fortran("comp_PI",
                    as.integer(k),          # Number of states
                    as.integer(n),          # Number of individuals
                    as.integer(TT),         # Number of time points
                    as.integer(nc),         # Number of covariates
                    nT  = as.integer(n * TT),          # Total number of rows in XX
                    XX  = XX,                           # Design matrix (passed through)
                    GG,                                 # Contrast array (passed through)
                    eta,                                # Linear predictor vector
                    PI  = array(0, c(k, k, n, TT)))    # Output array (initialised to 0)
    
    PI <- out$PI   # Extract the computed transition probability array
    
  } else {
    
    # =========================================================================
    # BRANCH 2
    # =========================================================================
    # -------------------------------------------------------------------------
    knc <- (k - 1) * (1 + nc)
    
    # Initialise the output array to zero.
    # Dimensions: [origin state u, destination state v, individual i, time t]
    # PI[u, v, i, 1] remains 0 for all u, v, i (no transition at t=1).
    PI <- array(0, c(k, k, n, TT))
    # -------------------------------------------------------------------------
    j <- 0
    
    # -------------------------------------------------------------------------
    # OUTER LOOP: over individuals i = 1, ..., n
    # -------------------------------------------------------------------------
    for (i in 1:n) {
      
      # Increment j to point to row (i, t=1) of XX.
      # This row is consumed here but not used (t=1 has no incoming transition).
      j <- j + 1
      
      # -----------------------------------------------------------------------
      # INNER LOOP: over time points t = 2, ..., TT
      # Transitions are defined only from t=2 onward because a transition
      # requires a previous state (at t-1). At t=1 there is no prior state.
      # -----------------------------------------------------------------------
      for (t in 2:TT) {
        
        # Increment j to point to row (i, t) of XX — the covariate vector at
        # the DESTINATION time point t (covariates are evaluated at time t).
        j <- j + 1
        
        # ---------------------------------------------------------------------
        # BUILD THE INDIVIDUAL-TIME COVARIATE DESIGN BLOCK (Xi):
        #
        # This structure means each free log-odds equation gets its own copy
        # of the covariate vector, allowing separate covariate effects for
        # each destination state within a given origin state.
        # ---------------------------------------------------------------------
        Xi <- diag(k - 1) %x% t(XX[j, ])
        
        # ---------------------------------------------------------------------
        # LOOP OVER ORIGIN STATES u = 1, ..., k:
        # For each origin state, compute the k transition probabilities
        # (the u-th row of the k×k transition matrix for individual i at time t).
        # ---------------------------------------------------------------------
        for (u in 1:k) {
          
          # -------------------------------------------------------------------
          # INDEX BLOCK FOR ORIGIN STATE u:
          # eta is a long vector of ALL transition log-odds parameters,
          # stacked by origin state
          # Each block has length knc = (k-1)*(1+nc).
          # ind selects the knc entries corresponding to origin state u.
          # -------------------------------------------------------------------
          ind <- (u - 1) * knc + (1:knc)
          
          # -------------------------------------------------------------------
          # COMPUTE UNNORMALISED LOG-ODDS (tmp):
          #
          #
          # -------------------------------------------------------------------
          tmp <- exp(c(GG[, , u] %*% Xi %*% eta[ind]))
          
          # -------------------------------------------------------------------
          # SOFTMAX NORMALISATION:
          # Divide each unnormalised odd by the sum across all k states.
          # This converts odds to valid probabilities that sum to 1.
          #
          #
          # -------------------------------------------------------------------
          PI[u, , i, t] <- tmp / sum(tmp)
          
        }  # end loop over origin states u
      }    # end loop over time points t
    }      # end loop over individuals i
    
  }  # end if/else fort
  
  # ---------------------------------------------------------------------------
  # RETURN the completed transition probability array
  # Dimensions: k × k × n × TT
  # ---------------------------------------------------------------------------
  return(PI)
  
}



# =============================================================================
# comp_Piv.R
# =============================================================================
# Purpose: Compute the matrix of INITIAL STATE PROBABILITIES for a
#          covariate-dependent first-order Markov Chain model.
#
#
# ARGUMENTS:
#   n    : integer — number of individuals
#   k    : integer — number of latent states
#   XX1  : matrix  — design matrix at t=1 only, dim = n × (nc+1)
#                    rows are one per individual; columns are (Intercept), X1, ..., Xnc
#                    Produced by subsetting the full XX matrix to rows where t=1.
#   G    : matrix  — contrast/expansion matrix for the initial distribution,
#                    dim = k × ((k-1)*(1+nc))
#                    Maps the (k-1) free log-odds to all k states by inserting
#                    a 0 row for the baseline state. Produced by design_matrices1().
#   la   : vector  — parameter vector for the initial distribution,
#                    length = (k-1) * (1+nc)
#                    Contains intercepts and covariate slopes for each of the
#                    (k-1) free log-odds equations. Set by the user in
#                    Simula_MC3state.R as the 'la' array.
#   fort : logical — if TRUE, calls a compiled Fortran routine for speed;
#                    if FALSE (default), uses the pure R implementation below.
#
# RETURNS:
#   Piv  : matrix of dim n × k containing initial state probabilities.
#          Piv[i, j] = P(S_{i,1} = j | X_{i,1})
#          Each row sums to 1 (valid probability distribution over k states).
#
# =============================================================================


comp_Piv <- function(n, k, XX1, G, la, fort = FALSE) {
  
  # ===========================================================================
  # BRANCH 1
  # ===========================================================================
  #
  # nc1 = ncol(XX1) - 1: the number of covariates (excluding the intercept),
  # derived here because Fortran needs it as an explicit scalar argument
  # (unlike the R branch, which does not use nc1 directly).
  #
  # Piv = matrix(0, n, k): the output matrix, initialised to 0 and passed
  # to Fortran as a writable buffer that will be filled by the subroutine.
  
  if (fort) {
    
    nc1 <- ncol(XX1) - 1   # Number of covariates (intercept column excluded)
    
    out <- .Fortran("comp_Piv",
                    as.integer(n),        # Number of individuals
                    as.integer(k),        # Number of states
                    as.integer(nc1),      # Number of covariates
                    XX1,                  # Design matrix at t=1 (n × (nc1+1))
                    G,                    # Contrast matrix (k × (k-1)*(1+nc1))
                    la,                   # Parameter vector (length (k-1)*(1+nc1))
                    Piv = matrix(0, n, k) # Output: initial probability matrix
    )
    
    Piv <- out$Piv   # Extract the computed initial probability matrix
    
  } else {
    
    # =========================================================================
    # BRANCH 2:  (fort = FALSE)
    # =========================================================================
    # Implements the multinomial logit for the initial distribution in plain R.
    
    # -------------------------------------------------------------------------
    # INITIALISE the output matrix to zero.
    # All n × k entries will be overwritten in the loop below.
    # -------------------------------------------------------------------------
    Piv <- matrix(0, n, k)
    # -------------------------------------------------------------------------
    for (i in 1:n) {
      Xi <- diag(k - 1) %x% t(XX1[i, ])
      tmp <- exp(c(G %*% Xi %*% la))
      Piv[i, ] <- tmp / sum(tmp)
      
    }  # end loop over individuals i
    
  }  # end if/else fort
  
  # ---------------------------------------------------------------------------
  # RETURN the completed initial probability matrix.

  # ---------------------------------------------------------------------------
  return(Piv)
}

# =============================================================================
# draw_mc.R
# =============================================================================
# Purpose: Simulate state sequences from a first-order Markov Chain (MC)
#          with individual- and time-varying transition probabilities.
#
#
#   The simulation draws one realisation per individual by sequentially
#   sampling from these multinomial distributions.
#
# ARGUMENTS:
#   Piv : matrix  — initial state probabilities
#
#   PI  : array   — transition probability matrices
#
# 
#   n  : integer — number of individuals
#   TT : integer — number of time points
#
# RETURNS:
#   A named list with one element:
#     $Y : data frame in long format with columns:
#            id   — individual identifier (1 to n)
#            time — time point (1 to TT)
#            Y1   — simulated state label at that (id, time) (integer, 1 to k)
#
# DEPENDENCIES:
#   matrices2long() from the LMest package — reshapes a wide n × TT matrix
#   into a long data frame with columns id, time, Y1.
# =============================================================================


draw_mc <- function(Piv, PI) {
  
  # ---------------------------------------------------------------------------
  # INITIALISE THE STATE MATRIX
  # Y is an n × TT integer matrix that will hold the simulated state for
  # each individual i at each time point t.
  # Initialised to 0; all entries will be overwritten in the loops below.
  # ---------------------------------------------------------------------------
  Y <- matrix(0, n, TT)
  
  
  # ---------------------------------------------------------------------------
  # STEP 1: DRAW INITIAL STATES (t = 1)
  #
  # For each individual i, draw their starting state from a k-category
  # Multinomial distribution with probability vector Piv[i, ].
  #
  #
  # ---------------------------------------------------------------------------
  for (i in 1:n) Y[i, 1] <- which(rmultinom(1, 1, Piv[i, ]) == 1)
  
  
  # ---------------------------------------------------------------------------
  # STEP 2: DRAW SUBSEQUENT STATES (t = 2, ..., TT)
  #
  # ---------------------------------------------------------------------------
  for (i in 1:n) {
    for (t in 2:TT) {
      Y[i, t] <- which(rmultinom(1, 1, PI[Y[i, t - 1], , i, t]) == 1)
    }
  }
  
  
  # ---------------------------------------------------------------------------
  # STEP 3: RESHAPE FROM WIDE TO LONG FORMAT
  #
  # ---------------------------------------------------------------------------
  res <- matrices2long(Y = Y)
  
  # Assemble the final long-format data frame with standardised column names.
  Y <- cbind("id"   = res$id,     # Individual identifier (1 to n)
             "time" = res$time,   # Time point (1 to TT)
             "Y1"   = res$Y1)     # Simulated state label (1 to k)
  
  
  # ---------------------------------------------------------------------------
  # RETURN a named list containing the simulated dataset
  # ---------------------------------------------------------------------------
  out <- list(Y = Y)
  
  # Note: no explicit return() call — in R, the last evaluated expression
  # in a function body is returned automatically. Here 'out' is returned.
  out
  
}