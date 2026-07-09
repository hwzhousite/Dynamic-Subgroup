# ============================================================
# WHOLE DATA SIMULATION CODE
# ============================================================

# -----------------------------------------------------------------------
# 0.  MARKET FLUCTUATION UTILITIES
#
# All trajectory functions accept an optional `fluct` argument (list) that
# adds market-realistic stochastic fluctuations on top of the smooth
# backbone.  When fluct = NULL the functions reproduce the original smooth
# curves exactly (backward-compatible).
#
# fluct parameters:
#   phi    — AR(1) persistence of the idiosyncratic perturbation (0..1).
#             Higher = smoother wiggles.  Default 0.82 (weekly equity ~).
#   sigma_g — SD of the per-group AR(1) innovation.  Default 0.18.
#   sigma_c — SD of a common-factor AR(1) shock added to all groups
#             (mimics market-wide moves).  Default 0.10.
#   phi_c   — AR(1) persistence of the common factor.  Default 0.90.
#   seed    — optional RNG seed for the fluctuations alone (so the
#             backbone and the fluctuations can be seeded separately).
# -----------------------------------------------------------------------

.ar1_path <- function(m, phi, sigma, seed = NULL) {
  # Generate a zero-mean AR(1) path of length m.
  if (!is.null(seed)) set.seed(seed)
  eps <- rnorm(m, 0, sigma * sqrt(1 - phi^2))   # stationary variance = sigma^2
  x   <- numeric(m)
  x[1] <- rnorm(1, 0, sigma)
  for (tt in seq(2, m)) x[tt] <- phi * x[tt - 1L] + eps[tt]
  x
}

.add_market_fluct <- function(backbone,   # K x m matrix (smooth curves)
                              fluct) {    # list of parameters or NULL
  if (is.null(fluct)) return(backbone)
  K <- nrow(backbone)
  m <- ncol(backbone)

  phi     <- if (!is.null(fluct$phi))     fluct$phi     else 0.82
  sigma_g <- if (!is.null(fluct$sigma_g)) fluct$sigma_g else 0.18
  sigma_c <- if (!is.null(fluct$sigma_c)) fluct$sigma_c else 0.10
  phi_c   <- if (!is.null(fluct$phi_c))   fluct$phi_c   else 0.90
  seed    <- fluct$seed  # may be NULL

  # Common-factor shock (same for all groups — mimics market beta)
  seed_c <- if (!is.null(seed)) seed         else NULL
  seed_g <- if (!is.null(seed)) seed + 1000L else NULL
  common <- .ar1_path(m, phi_c, sigma_c, seed = seed_c)

  # Per-group idiosyncratic shocks
  group_paths <- matrix(NA_real_, K, m)
  for (k in seq_len(K)) {
    sk <- if (!is.null(seed_g)) seed_g + k else NULL
    group_paths[k, ] <- .ar1_path(m, phi, sigma_g, seed = sk)
  }

  backbone + common[col(backbone)] + group_paths
}

# Like .add_market_fluct but multiplies each group's idiosyncratic AR(1)
# shock by a per-group taper vector.  When idio_taper[k, t] = 0 the group
# carries no idiosyncratic noise at time t, guaranteeing that groups with
# identical backbones are EXACTLY identical in the full noisy trajectory.
#
# idio_taper : K x m numeric matrix, values in [0, 1].
#              NULL = all-ones (equivalent to .add_market_fluct).
.add_market_fluct_tapered <- function(backbone, fluct, idio_taper = NULL) {
  if (is.null(fluct)) return(backbone)
  K <- nrow(backbone); m <- ncol(backbone)

  phi     <- if (!is.null(fluct$phi))     fluct$phi     else 0.82
  sigma_g <- if (!is.null(fluct$sigma_g)) fluct$sigma_g else 0.18
  sigma_c <- if (!is.null(fluct$sigma_c)) fluct$sigma_c else 0.10
  phi_c   <- if (!is.null(fluct$phi_c))   fluct$phi_c   else 0.90
  seed    <- fluct$seed

  seed_c  <- if (!is.null(seed)) seed         else NULL
  seed_g  <- if (!is.null(seed)) seed + 1000L else NULL
  common  <- .ar1_path(m, phi_c, sigma_c, seed = seed_c)

  group_paths <- matrix(NA_real_, K, m)
  for (k in seq_len(K)) {
    sk               <- if (!is.null(seed_g)) seed_g + k else NULL
    idio             <- .ar1_path(m, phi, sigma_g, seed = sk)
    tp               <- if (!is.null(idio_taper)) idio_taper[k, ] else rep(1.0, m)
    group_paths[k, ] <- idio * tp
  }

  backbone + common[col(backbone)] + group_paths
}

# -----------------------------------------------------------------------
# 1.  TRUE TRAJECTORY FUNCTIONS
#     Each accepts `fluct = NULL` (smooth) or a list (market-realistic).
# -----------------------------------------------------------------------

gamma_crossing <- function(t, ...) {
  # Sinusoidal common base (the trajectory level at the crossing moment)
  base <- 3.0 + 0.50*sin(2*pi*t + 0.3) + 0.25*sin(5*pi*t  + 0.8) +
          0.12*sin(9*pi*t)

  # Separation amplitude — guaranteed positive so k1 > k2 before crossing
  dev  <- 2.20 + 0.70*sin(2*pi*t + 0.5) + 0.35*sin(4*pi*t  + 1.2) +
          0.18*sin(7*pi*t + 0.3)

  # Crossing sign: +1 (before) → 0 (exact crossing at t=0.50) → −1 (after)
  # Identity: 2*cosine_taper - 1 = cos(pi*(t-t_s)/(t_e-t_s)) on [t_s, t_e]
  tap_sign <- 2 * cosine_taper(t, 0.30, 0.70) - 1

  backbone <- rbind(
    k1 = base + dev * tap_sign,   # starts high, ends low
    k2 = base - dev * tap_sign    # starts low,  ends high
  )

  backbone
}

gamma_merging <- function(t, ...) {
  # Sinusoidal common base (the level both groups converge to once merged)
  base <- 2.5 + 0.45*sin(2*pi*t + 0.2) + 0.22*sin(5*pi*t  + 0.7) +
          0.11*sin(9*pi*t)

  # Per-group deviations: dev1 > 0 and dev2 < 0 always, so k1 > k2 pre-merge
  dev1 <- +1.80 + 0.65*sin(2*pi*t + 0.3) + 0.32*sin(6*pi*t  + 0.5) +
           0.15*sin(10*pi*t + 1.0)
  dev2 <- -1.60 + 0.60*sin(2*pi*t + 1.0) + 0.30*sin(6*pi*t  + 0.4) +
           0.14*sin(9*pi*t  + 0.8)

  # Merge taper: separated before t = 0.35, fully merged by t = 0.65
  tap <- cosine_taper(t, 0.35, 0.65)

  backbone <- rbind(
    k1 = base + dev1 * tap,
    k2 = base + dev2 * tap
  )

  backbone
}

gamma_separated <- function(t, ...) {
  # Sinusoidal common base
  base <- 3.5 + 0.40*sin(2*pi*t + 0.5) + 0.20*sin(6*pi*t  + 1.2) +
          0.10*sin(10*pi*t)

  # Per-group deviations: permanently active (no merge taper)
  # dev1 > 0 and dev2 < 0 always → groups never cross or merge
  dev1 <- +2.00 + 0.65*sin(2*pi*t + 0.2) + 0.32*sin(5*pi*t  + 0.8) +
           0.16*sin(9*pi*t  + 1.3)
  dev2 <- -2.20 + 0.70*sin(2*pi*t + 1.5) + 0.35*sin(4*pi*t  + 0.3) +
           0.18*sin(7*pi*t  + 0.6)

  backbone <- rbind(
    k1 = base + dev1,
    k2 = base + dev2
  )

  backbone
}

cosine_taper <- function(t, t_s, t_e) {
  ifelse(
    t <= t_s, 1,
    ifelse(
      t >= t_e, 0,
      0.5 * (1 + cos(pi * (t - t_s) / (t_e - t_s)))
    )
  )
}

gamma_sequential_merging4 <- function(t, ...) {
  # Sinusoidal common base shared by all four groups once merged
  base <- 3.5 + 0.60 * sin(2*pi*t + 0.4) + 0.30 * sin(6*pi*t  + 1.0) +
          0.15 * sin(12*pi*t)

  # Per-group deviations (active only while the group is separated)
  dev1 <- +2.80 + 0.90 * sin(2*pi*t)       + 0.45 * sin(5*pi*t  + 0.2) +
           0.20 * sin(9*pi*t  + 1.1)
  dev2 <- +0.90 + 0.70 * sin(2*pi*t + 0.5) + 0.35 * sin(7*pi*t  + 0.8) +
           0.18 * sin(11*pi*t + 0.6)
  dev3 <- -0.80 + 0.75 * sin(2*pi*t + 1.0) + 0.40 * sin(5*pi*t  + 1.5) +
           0.20 * sin(10*pi*t + 0.9)
  dev4 <- -2.60 + 0.85 * sin(2*pi*t + 1.8) + 0.42 * sin(6*pi*t  + 0.3) +
           0.22 * sin(8*pi*t  + 2.0)

  # Structural tapers (mirror the three merge windows in the figure)
  tap12 <- cosine_taper(t, 0.20, 0.42)   # k1+k2: fully merged by t = 0.42
  tap34 <- cosine_taper(t, 0.42, 0.65)   # k3+k4: fully merged by t = 0.65
  tapa  <- cosine_taper(t, 0.65, 0.85)   # all four: common base by t = 0.85

  # Backbone — merged pairs are EXACTLY equal once their taper = 0:
  #   k1 = k2  for t >= 0.42  (tap12 = 0)
  #   k3 = k4  for t >= 0.65  (tap34 = 0)
  #   k1=k2=k3=k4 = base  for t >= 0.85  (tapa = 0)
  backbone <- rbind(
    k1 = base + dev1 * tap12 + dev2 * (1 - tap12) * tapa,
    k2 = base + dev2 * tap12 + dev2 * (1 - tap12) * tapa,
    k3 = base + dev3 * tap34 + dev3 * (1 - tap34) * tapa,
    k4 = base + dev4 * tap34 + dev3 * (1 - tap34) * tapa
  )

  backbone
}

# -----------------------------------------------------------------------
# gamma_crossing3: K=3 — two trajectories cross, one permanently separate.
#
# Structure:
#   k1, k2 — crossing pair with sinusoidal backbone; exact crossing
#             (k1 = k2 = base_12) at t = 0.50 via sign-reversing taper
#   k3      — permanently below both k1 and k2 throughout [0, 1];
#             separation guaranteed: max(k3) < min(k1,k2) for all t
# -----------------------------------------------------------------------

gamma_crossing3 <- function(t, ...) {
  # Sinusoidal base shared by the crossing pair
  base_12 <- 3.5 + 0.50*sin(2*pi*t + 0.3) + 0.25*sin(5*pi*t + 0.8) +
             0.12*sin(9*pi*t)

  # Separation amplitude — guaranteed positive (min ≈ 0.67)
  dev_12  <- 1.80 + 0.65*sin(2*pi*t + 0.5) + 0.32*sin(4*pi*t + 1.2) +
             0.16*sin(7*pi*t + 0.3)

  # Crossing sign: +1 → 0 → −1 via smooth cosine over [0.30, 0.70]
  # At midpoint t = 0.50: tap_sign = 0, so k1 = k2 = base_12 exactly
  tap_sign <- 2 * cosine_taper(t, 0.30, 0.70) - 1

  # Separated trajectory — guaranteed below min(k1, k2) for all t
  # Proof: min(base_12 - dev_12) ≥ 2.63 - 2.93 = -0.30
  #        max(k3)               ≈ -1.50 + 0.97  = -0.53 < -0.30  ✓
  k3 <- -1.50 + 0.55*sin(2*pi*t + 1.5) + 0.28*sin(5*pi*t + 0.4) +
         0.14*sin(8*pi*t + 2.0)

  rbind(
    k1 = base_12 + dev_12 * tap_sign,   # starts high, ends low
    k2 = base_12 - dev_12 * tap_sign,   # starts low,  ends high
    k3 = k3                              # always below crossing pair
  )
}

# -----------------------------------------------------------------------
# gamma_branching: K=2 — two trajectories start identical and diverge.
#
# Mirror image of gamma_merging:
#   t < 0.35  : k1 = k2 = base  (EXACTLY identical — not yet branched)
#   0.35–0.65 : smooth cosine divergence via branching taper
#   t > 0.65  : k1 and k2 fully separated (dev1 > 0, dev2 < 0 always)
# -----------------------------------------------------------------------

gamma_branching <- function(t, ...) {
  # Sinusoidal common base (the shared level before branching)
  base <- 3.0 + 0.48*sin(2*pi*t + 0.7) + 0.24*sin(5*pi*t + 1.3) +
          0.12*sin(8*pi*t + 0.5)

  # Per-group deviations — active only after branching
  # dev1 min ≈ 1.70 − 0.62 − 0.31 − 0.15 = 0.62 > 0  (k1 always above base)
  # dev2 max ≈ −1.55 + 0.58 + 0.29 + 0.14 = −0.54 < 0 (k2 always below base)
  dev1 <- +1.70 + 0.62*sin(2*pi*t + 0.4) + 0.31*sin(7*pi*t + 0.9) +
           0.15*sin(11*pi*t + 1.4)
  dev2 <- -1.55 + 0.58*sin(2*pi*t + 1.3) + 0.29*sin(6*pi*t + 0.2) +
           0.14*sin(10*pi*t + 0.7)

  # Branching taper: 0 → 1 (reverse of merge taper)
  # tap = 0 for t <= 0.35 → k1 = k2 = base  (exact identity)
  # tap = 1 for t >= 0.65 → k1 = base+dev1, k2 = base+dev2 (fully separated)
  tap <- 1 - cosine_taper(t, 0.35, 0.65)

  rbind(
    k1 = base + dev1 * tap,   # rises after branching
    k2 = base + dev2 * tap    # falls after branching
  )
}

TRAJ_FNS <- list(
  crossing             = gamma_crossing,
  merging              = gamma_merging,
  separated            = gamma_separated,
  sequential_merging4  = gamma_sequential_merging4,
  crossing3            = gamma_crossing3,
  branching            = gamma_branching
)

# -----------------------------
# 2. TRAJECTORY RESOLVER
# -----------------------------

resolve_trajectory <- function(traj_type,
                               t_seq,
                               K,
                               p,
                               trajectory_args = list()) {
  if (!traj_type %in% names(TRAJ_FNS)) {
    stop(
      "Unknown traj_type = '", traj_type, "'. Available options are: ",
      paste(names(TRAJ_FNS), collapse = ", ")
    )
  }
  
  G_raw <- do.call(
    TRAJ_FNS[[traj_type]],
    c(list(t = t_seq), trajectory_args)
  )
  
  if (!is.matrix(G_raw)) {
    stop("Trajectory function must return a K x m matrix.")
  }
  
  if (nrow(G_raw) != K) {
    stop(
      "Trajectory '", traj_type, "' returns K = ", nrow(G_raw),
      " groups, but simulate_data() was called with K = ", K, "."
    )
  }
  
  if (ncol(G_raw) != length(t_seq)) {
    stop("Trajectory must have length(t_seq) columns.")
  }
  
  G <- array(NA_real_, dim = c(K, length(t_seq), p))
  
  for (j in seq_len(p)) {
    G[, , j] <- G_raw * sqrt(j)
  }
  
  G
}

# -----------------------------
# 3. DYNAMIC SUBGROUP PATHS
# -----------------------------

simulate_group_paths <- function(n,
                                 m,
                                 K,
                                 group_probs = rep(1 / K, K),
                                 class_change = FALSE,
                                 change_time = NULL,
                                 change_prob = 0.30,
                                 transition_matrix = NULL,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (is.null(change_time)) {
    change_time <- floor(m / 2)
  }
  
  if (length(group_probs) != K) {
    stop("group_probs must have length K.")
  }
  
  if (abs(sum(group_probs) - 1) > 1e-8) {
    stop("group_probs must sum to 1.")
  }
  
  if (change_time < 1L || change_time > m) {
    stop("change_time must be between 1 and m.")
  }
  
  g0 <- sample(
    seq_len(K),
    size = n,
    replace = TRUE,
    prob = group_probs
  )
  
  group_mat <- matrix(NA_integer_, nrow = n, ncol = m)
  group_mat[, ] <- g0
  
  if (!isTRUE(class_change)) {
    return(group_mat)
  }
  
  if (is.null(transition_matrix)) {
    transition_matrix <- matrix(
      change_prob / (K - 1),
      nrow = K,
      ncol = K
    )
    diag(transition_matrix) <- 1 - change_prob
  }
  
  if (!is.matrix(transition_matrix)) {
    stop("transition_matrix must be a K x K matrix.")
  }
  
  if (!all(dim(transition_matrix) == c(K, K))) {
    stop("transition_matrix must have dimension K x K.")
  }
  
  if (any(transition_matrix < 0)) {
    stop("transition_matrix cannot contain negative probabilities.")
  }
  
  row_sums <- rowSums(transition_matrix)
  
  if (any(abs(row_sums - 1) > 1e-8)) {
    stop("Each row of transition_matrix must sum to 1.")
  }
  
  if (change_time > 1L) {
    group_mat[, 1:(change_time - 1L)] <- g0
  }
  
  for (i in seq_len(n)) {
    current_group <- g0[i]
    
    for (tt in change_time:m) {
      current_group <- sample(
        seq_len(K),
        size = 1L,
        prob = transition_matrix[current_group, ]
      )
      
      group_mat[i, tt] <- current_group
    }
  }
  
  group_mat
}

# -----------------------------
# 4. MISSING DATA INJECTION
# -----------------------------

inject_missing <- function(dat,
                           miss_mech = c(
                             "trunc_start",    # deterministic left truncation
                             "random_entry",   # random late entry
                             "random_drop",    # random scattered dropout
                             "random_exit",    # random early exit
                             "random_both"     # random entry + exit
                           ),
                           trunc_prop = 0.25,
                           miss_prop  = 0.25,
                           seed = NULL) {
  miss_mech <- match.arg(miss_mech)
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(dat$Y)
  m <- ncol(dat$Y)
  
  has_X <- !is.null(dat$X)
  has_Z <- !is.null(dat$Z)
  
  n_trunc <- max(1L, floor(trunc_prop * m))
  exit_earliest <- max(1L, floor((1 - trunc_prop) * m))
  
  for (i in seq_len(n)) {
    mask_t <- switch(
      miss_mech,

      # Deterministic: mask the first n_trunc time points for every subject
      trunc_start = seq_len(n_trunc),

      random_entry = {
        entry_i <- sample.int(n_trunc, size = 1L)
        if (entry_i <= 1L) integer(0L) else seq_len(entry_i - 1L)
      },

      random_drop = which(runif(m) < miss_prop),

      random_exit = {
        exit_i <- sample.int(m - exit_earliest + 1L, size = 1L) +
          exit_earliest - 1L
        if (exit_i >= m) integer(0L) else seq.int(exit_i + 1L, m)
      },

      # BUG-1 FIX: use explicit left/right masks with early-exit guards
      random_both = {
        entry_i <- sample.int(n_trunc, size = 1L)
        exit_i  <- sample.int(m - exit_earliest + 1L, size = 1L) +
          exit_earliest - 1L
        # ensure at least one observed point
        if (exit_i <= entry_i) exit_i <- min(entry_i + 1L, m)
        left_mask  <- if (entry_i <= 1L) integer(0L) else seq_len(entry_i - 1L)
        right_mask <- if (exit_i  >= m)  integer(0L) else seq.int(exit_i + 1L, m)
        union(left_mask, right_mask)
      }
    )
    
    if (length(mask_t) == 0L) next
    
    dat$Y[i, mask_t] <- NA_real_
    
    if (has_X) {
      dat$X[i, mask_t, ] <- NA_real_
    }
    
    if (has_Z) {
      dat$Z[i, mask_t, ] <- NA_real_
    }
  }
  
  dat
}

# -----------------------------
# 5. UNIFIED DATA SIMULATION
# -----------------------------

simulate_data <- function(n,
                          m = 24L,
                          K = 2L,
                          p = 1L,
                          r = 0L,
                          sigma = 1,
                          traj_type = c(
                            "crossing",
                            "merging",
                            "separated",
                            "sequential_merging4",
                            "crossing3",
                            "branching"
                          ),
                          trajectory_args = list(),

                          # fluct: AR(1) market fluctuation parameters passed
                          # into trajectory functions.  NULL = smooth curves.
                          # list(phi, sigma_g, sigma_c, phi_c, seed)
                          fluct = NULL,

                          # garch_alpha: GARCH(1,1) weight on lagged squared
                          # residual for volatility clustering in Y.
                          # 0 = homoskedastic.  ~0.15 mimics equity data.
                          garch_alpha = 0,

                          alpha_true = NULL,
                          group_probs = rep(1 / K, K),
                          t_seq = seq(0, 1, length.out = m),
                          x_mean = 1,
                          x_sd = 0.1,
                          z_mean = 0,
                          z_sd = 1,
                          seed = NULL,
                          
                          # subgroup transition
                          class_change = FALSE,
                          change_time = NULL,
                          change_prob = 0.20,
                          transition_matrix = NULL,
                          group_seed = NULL,
                          
                          # missing data
                          missing = FALSE,
                          miss_mech = c(
                            "trunc_start",
                            "random_entry",
                            "random_drop",
                            "random_exit",
                            "random_both"
                          ),
                          trunc_prop = 0.25,
                          miss_prop = 0.20,
                          miss_seed = NULL) {
  traj_type <- match.arg(traj_type,
    choices = c("crossing", "merging", "separated",
                "sequential_merging4", "crossing3", "branching"))
  miss_mech <- match.arg(miss_mech)
  
  if (!is.null(seed)) set.seed(seed)
  
  stopifnot(n > 0, K > 0, p > 0, r >= 0, sigma > 0)
  stopifnot(length(group_probs) == K)
  stopifnot(abs(sum(group_probs) - 1) < 1e-8)
  
  m <- length(t_seq)
  
  G <- resolve_trajectory(
    traj_type = traj_type,
    t_seq = t_seq,
    K = K,
    p = p,
    trajectory_args = trajectory_args
  )
  
  true_groups <- simulate_group_paths(
    n = n,
    m = m,
    K = K,
    group_probs = group_probs,
    class_change = class_change,
    change_time = change_time,
    change_prob = change_prob,
    transition_matrix = transition_matrix,
    seed = group_seed
  )
  
  X <- array(
    rnorm(n * m * p, mean = x_mean, sd = x_sd),
    dim = c(n, m, p)
  )
  
  if (r > 0L) {
    Z <- array(
      rnorm(n * m * r, mean = z_mean, sd = z_sd),
      dim = c(n, m, r)
    )
    
    if (is.null(alpha_true)) {
      alpha_true <- matrix(1, nrow = m, ncol = r)
    } else if (is.vector(alpha_true) && length(alpha_true) == r) {
      alpha_true <- matrix(
        rep(alpha_true, each = m),
        nrow = m,
        ncol = r
      )
    }
    
    stopifnot(
      is.matrix(alpha_true),
      nrow(alpha_true) == m,
      ncol(alpha_true) == r
    )
  } else {
    Z <- NULL
    alpha_true <- NULL
  }
  
  beta_true <- array(NA_real_, dim = c(n, m, p))
  Y <- matrix(NA_real_, n, m)

  # GARCH(1,1) volatility: h_t = sigma^2*(1 - garch_alpha) + garch_alpha*eps_{t-1}^2
  # When garch_alpha = 0 this collapses to constant sigma (original behaviour).
  garch_beta0 <- sigma^2 * (1 - min(garch_alpha, 0.9))

  for (i in seq_len(n)) {
    h_t   <- sigma^2   # conditional variance at t=1
    eps_prev <- 0      # lagged residual

    for (tt in seq_len(m)) {
      k <- true_groups[i, tt]

      beta_it <- G[k, tt, ]
      beta_true[i, tt, ] <- beta_it

      mu_it <- sum(X[i, tt, ] * beta_it)
      if (!is.null(Z))
        mu_it <- mu_it + sum(Z[i, tt, ] * alpha_true[tt, ])

      # GARCH update: h_t = omega + alpha * eps_{t-1}^2  (simplified, beta=0)
      if (tt > 1L)
        h_t <- garch_beta0 + garch_alpha * eps_prev^2

      eps_it <- rnorm(1, 0, sqrt(h_t))
      Y[i, tt] <- mu_it + eps_it
      eps_prev  <- eps_it
    }
  }
  
  if (K == 2L) {
    gap <- G[1, , 1] - G[2, , 1]
    snr <- mean(gap^2) / (4 * sigma^2)
  } else {
    signal_t <- sapply(seq_len(m), function(tt) {
      stats::var(as.vector(G[, tt, ]))
    })
    snr <- mean(signal_t) / sigma^2
  }
  
  dat <- list(
    Y = Y,
    X = X,
    Z = Z,
    alpha_true = alpha_true,
    beta_true = beta_true,
    true_groups = true_groups,
    initial_groups = true_groups[, 1],
    G = G,
    snr = snr,
    t_seq = t_seq,
    traj_type = traj_type,
    class_change = class_change,
    change_time = change_time,
    change_prob = change_prob,
    transition_matrix = transition_matrix,
    missing = FALSE,
    miss_mech = NULL,
    trunc_prop = NULL,
    miss_prop = NULL
  )
  
  if (isTRUE(missing)) {
    dat <- inject_missing(
      dat = dat,
      miss_mech = miss_mech,
      trunc_prop = trunc_prop,
      miss_prop = miss_prop,
      seed = miss_seed
    )
    
    dat$missing <- TRUE
    dat$miss_mech <- miss_mech
    dat$trunc_prop <- trunc_prop
    dat$miss_prop <- miss_prop
  }
  
  dat
}


