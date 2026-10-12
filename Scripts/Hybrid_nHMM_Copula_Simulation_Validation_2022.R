

# PART 1
##########################################################################################################
## Hybrid nHMM – State-conditional margins (zero-inflated + ECDF+GPD) – R-vine copula (zone-wise)     ##
## Simulation for 2022 (MJJA), diagnostics (PIT, Rosenblatt, pairwise GoF), and monthly validation.    ##
## Author: (MAURICIO HERRERA-MARÍN)                                                                                       ##
##########################################################################################################
## Run from the repository root, e.g.:
## Rscript Scripts/Hybrid_nHMM_Copula_Simulation_Validation_2022.R
## The input CSV is shipped inside the ZIP under Data/.
base_path <- normalizePath(getwd(), mustWork = TRUE)
data_dir <- file.path(base_path, "Data")
data_zip_path <- file.path(data_dir, "mapa_estaciones_precipitaciones_clean.csv.zip")
stations_csv_path <- file.path(data_dir, "mapa_estaciones_precipitaciones_clean.csv")
if (!dir.exists(data_dir)) {
  stop("Data/ not found. Run the script from the root of the GitHub repository.")
}
if (!file.exists(stations_csv_path)) {
  if (!file.exists(data_zip_path)) {
    stop("Missing input ZIP: ", data_zip_path)
  }
  zip_entries <- utils::unzip(data_zip_path, list = TRUE)$Name
  csv_entries <- zip_entries[grepl("\\.csv$", zip_entries, ignore.case = TRUE)]
  if (length(csv_entries) != 1L) {
    stop("Expected exactly one CSV in ", data_zip_path,
         "; found ", length(csv_entries), ". Check the archive manually.")
  }
  utils::unzip(data_zip_path, files = csv_entries, exdir = data_dir,
               junkpaths = TRUE, overwrite = FALSE)
  extracted_csv <- file.path(data_dir, basename(csv_entries))
  if (!file.exists(stations_csv_path) && file.exists(extracted_csv)) {
    if (!file.copy(extracted_csv, stations_csv_path)) {
      stop("Unable to copy extracted CSV to expected filename: ", stations_csv_path)
    }
  }
}
if (!file.exists(stations_csv_path)) {
  stop("Input CSV was not found after extraction: ", stations_csv_path)
}
datos_base_l <- read.csv(stations_csv_path)

## ------------------------------ Libraries -------------------------------------------------------------
pacman::p_load(
  dplyr, tidyr, purrr, tibble, data.table,
  lubridate, ggplot2, stringr,
  VineCopula, copula, evd,
  parallel, goftest, ADGofTest
)

## ------------------------------ Theming (optional) ---------------------------------------------------
tema_box <- theme_bw() +
  theme(
    plot.title   = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 13),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey92")
  )

## ------------------------------ Parameters -----------------------------------------------------------
set.seed(42)

selected_zone  <- 5            # <- change zone id here
M_per_day      <- 100          # number of simulations per day (your setting)
train_end      <- as.Date("2021-08-31")
months_MJJA    <- 5:8
test_year      <- 2022

## Vine fitting options (expanded families + truncation to avoid overfit in high d)
vine_familyset <- NA           # NA => allow all (incl. rotated BB families) in VineCopula
vine_selcrit   <- "BIC"        # penalize complexity
vine_truncate  <- 2            # truncate vine at Tree 2 (robust in d~50–60)
min_rows_vine  <- 200          # minimum rows per state to fit vine

## Tail model options
u_star         <- 0.90         # body/upper-tail threshold for positive precipitation
min_excess_gpd <- 25           # fit GPD only if enough exceedances

## ------------------------------ Expect data ----------------------------------------------------------
stopifnot(all(c("dates","codigo_estacion","precipitation","labels","States","medoids") %in% names(datos_base_l)))
datos_base_l <- datos_base_l %>%
  mutate(dates = as.Date(dates),
         Year  = year(dates),
         is_train = dates <= train_end)

## ------------------------------ Select zone & medoid -------------------------------------------------
df_zone <- datos_base_l %>% filter(labels == selected_zone)
medoid_station <- df_zone %>% filter(medoids == 1) %>% pull(codigo_estacion) %>% unique()
if(length(medoid_station) != 1) stop("Medoid not found or not unique in selected zone.")

## ------------------------------ State-conditional zero-inflated margins ------------------------------
# We model each (station, state) margin F_{m|s} as:
# P(R=0) = p0; for R>0 use ECDF in the body and GPD in the upper tail above the u* quantile.
# This design aligns with PIT diagnostics and improves extremes.

fit_margin_state <- function(x, u_star = 0.90, min_excess = 25){
  x <- x[is.finite(x)]
  if(length(x) < 40L){
    # too few samples => fallback: point-mass at 0 + ECDF on positives only
    p0   <- mean(x <= 0, na.rm = TRUE)
    xpos <- sort(x[x > 0])
    Fe   <- if(length(xpos) > 0) stats::ecdf(xpos) else NULL
    return(list(p0 = p0, Fe = Fe, x_sorted = xpos,
                thr = Inf, Fu = 1, gpd = NULL))
  }
  p0   <- mean(x <= 0, na.rm = TRUE)
  xpos <- sort(x[x > 0])
  if(length(xpos) < 10L){
    Fe <- NULL
    return(list(p0 = p0, Fe = Fe, x_sorted = xpos,
                thr = Inf, Fu = 1, gpd = NULL))
  }
  Fe   <- stats::ecdf(xpos)
  thr  <- stats::quantile(xpos, probs = u_star, na.rm = TRUE, type = 7)
  Fu   <- Fe(thr)  # CDF at threshold within positive support
  exc  <- xpos[xpos > thr] - thr
  gpd  <- NULL
  if(length(exc) >= min_excess){
    gpd <- tryCatch(evd::fpot(xpos, threshold = thr, model = "gpd", std.err = FALSE),
                    error = function(e) NULL)
  }
  list(p0 = p0, Fe = Fe, x_sorted = xpos,
       thr = thr, Fu = Fu, gpd = gpd)
}

# Hybrid CDF (randomized at zero for PIT stability)
F_hurdle <- function(r, model){
  if(is.null(model)) return(NA_real_)
  p0 <- model$p0
  if(r <= 0){
    # randomized PIT within [0, p0]; for reproducibility use deterministic jitter if needed
    return(runif(1, min = 0, max = p0))
  }
  Fe <- model$Fe
  if(is.null(Fe)) return(p0)  # no positive support info
  if(r <= model$thr || !is.finite(model$thr) || is.null(model$gpd)){
    return(p0 + (1 - p0) * Fe(r))
  }else{
    # upper tail via GPD
    p_ex  <- evd::pgpd(r - model$thr, loc = 0,
                       scale = model$gpd$estimate["scale"],
                       shape = model$gpd$estimate["shape"])
    return(p0 + (1 - p0) * (model$Fu + (1 - model$Fu)*p_ex))
  }
}

# Hybrid quantile for simulation (invert CDF approximately)
Q_hurdle <- function(u, model){
  if(is.null(model)) return(NA_real_)
  p0 <- model$p0
  if(u <= p0) return(0)  # atom at zero
  if(is.null(model$Fe) || length(model$x_sorted) == 0) return(0)
  u_pos <- (u - p0) / (1 - p0)
  # if below tail threshold, use empirical quantile of positive support
  if(u_pos <= model$Fu || is.null(model$gpd) || !is.finite(model$thr)){
    return(as.numeric(stats::quantile(model$x_sorted, probs = u_pos, type = 7)))
  }
  # in the upper tail: invert GPD
  p_tail <- (u_pos - model$Fu) / (1 - model$Fu)
  par <- model$gpd$estimate
  q_ex <- evd::qgpd(p_tail, loc = 0, scale = par["scale"], shape = par["shape"])
  return(as.numeric(model$thr + q_ex))
}

## Fit (station, state) margins on TRAIN
marginals_tbl <- datos_base_l %>%
  filter(is_train) %>%
  group_by(codigo_estacion, States) %>%
  summarise(model = list(fit_margin_state(precipitation, u_star, min_excess_gpd)),
            .groups = "drop")

## ------------------------------ Build U-matrices per state (TRAIN) -----------------------------------
# We compute U_{t,m} = F_{m|s}(R_{t,m}) using the state-conditional margin model.
# These matrices are used to fit a vine per state within the selected zone.

calc_U_zone_state <- function(df_zone, state_id, marginals_tbl){
  df <- df_zone %>%
    filter(is_train, States == state_id) %>%
    select(dates, codigo_estacion, States, precipitation)
  if(nrow(df) == 0) return(NULL)
  
  dfU <- df %>%
    left_join(marginals_tbl, by = c("codigo_estacion","States")) %>%
    mutate(U = purrr::map2_dbl(precipitation, model, ~ F_hurdle(.x, .y))) %>%
    mutate(codigo_estacion = as.character(codigo_estacion)) %>%
    select(dates, codigo_estacion, U) %>%
    tidyr::pivot_wider(names_from = codigo_estacion, values_from = U) %>%
    arrange(dates) %>%
    suppressWarnings(tidyr::drop_na())
  
  if(nrow(dfU) == 0) return(NULL)
  list(U = as.matrix(dfU %>% select(-dates)),
       station_ids = colnames(dfU %>% select(-dates)))
}

## ------------------------------ Fit a vine per state in the zone -------------------------------------
vine_models <- vector("list", length = 5)
names(vine_models) <- as.character(sort(unique(df_zone$States)))

for(st in sort(unique(df_zone$States))){
  Uz <- calc_U_zone_state(df_zone, st, marginals_tbl)
  if(is.null(Uz) || nrow(Uz$U) < min_rows_vine){
    vine_models[[as.character(st)]] <- NULL
    next
  }
  message(sprintf("Fitting vine: Zone %s, State %s  (n=%d, d=%d)",
                  selected_zone, st, nrow(Uz$U), ncol(Uz$U)))
  vine_fit <- tryCatch(
    RVineStructureSelect(Uz$U,
                         familyset     = vine_familyset,   # all families incl. rotations
                         selectioncrit = vine_selcrit,
                         trunclevel    = vine_truncate),
    error = function(e) { message("  -> vine fit failed: ", e$message); NULL }
  )
  vine_models[[as.character(st)]] <- list(vine = vine_fit,
                                          stations = Uz$station_ids)
}

## ------------------------------ Goodness-of-fit diagnostics (FIXED) ----------------------------------
## (i) Rosenblatt (multivariate PIT) + univariate EDF tests
## (ii) Pairwise GoF for first tree (Tree 1) using the vine's internal matrices

# Helper: extract Tree 1 pairs (edges) from an RVine object via its matrices
extract_T1_pairs <- function(rvine_obj){
  if(is.null(rvine_obj) || is.null(rvine_obj$Matrix)) return(NULL)
  d <- ncol(rvine_obj$Matrix)
  # For Tree 1, VineCopula encodes pairs using last column + diagonal of Matrix
  tibble::tibble(
    i      = rvine_obj$Matrix[d, (d-1):1],
    j      = diag(rvine_obj$Matrix)[(d-1):1],
    family = rvine_obj$family[(d-1):1, d],
    par    = rvine_obj$par[(d-1):1, d],
    par2   = rvine_obj$par2[(d-1):1, d]
  )
}

gof_results <- purrr::imap_dfr(vine_models, function(vm, st){
  # vm is list(vine=<RVineMatrix>, stations=<chr>)
  if(is.null(vm) || is.null(vm$vine)) return(tibble(State = as.integer(st)))
  
  # Build U for this state on TRAIN (same margins you fitted)
  Uz <- calc_U_zone_state(df_zone, as.integer(st), marginals_tbl)
  if(is.null(Uz) || nrow(Uz$U) == 0) return(tibble(State = as.integer(st)))
  
  # (i) Rosenblatt transform (multivariate PIT) and univariate uniformity (first components)
  Rz <- tryCatch(VineCopula::RVinePIT(Uz$U, vm$vine), error = function(e) NULL)
  uni_tab <- if(!is.null(Rz)){
    ktest <- min(5L, ncol(Rz))
    purrr::map_dfr(seq_len(ktest), function(j){
      u <- Rz[, j]
      tibble::tibble(
        comp = j,
        ks_p  = ks.test(u, "punif")$p.value,
        cvm_p = goftest::cvm.test(u, "punif")$p.value
      )
    })
  } else NULL
  
  # (ii) First-tree (Tree 1) pairs via matrices + parametric White GoF
  pairs <- extract_T1_pairs(vm$vine)
  pair_p_min <- NA_real_
  if(!is.null(pairs) && nrow(pairs) > 0){
    pvals <- purrr::pmap_dbl(pairs, function(i, j, family, par, par2){
      ui <- Uz$U[, as.integer(i)]
      uj <- Uz$U[, as.integer(j)]
      # Parametric GoF for the fitted family/params (fast B; increase to 1000 for final runs)
      tt <- tryCatch(
        VineCopula::BiCopGofTest(ui, uj, family = family, par = par, par2 = par2,
                                 method = "White", B = 200),
        error = function(e) NULL
      )
      if(is.null(tt)) NA_real_ else tt$p.value
    })
    pair_p_min <- suppressWarnings(min(pvals, na.rm = TRUE))
    if(!is.finite(pair_p_min)) pair_p_min <- NA_real_
  }
  
  tibble::tibble(
    State       = as.integer(st),
    n           = nrow(Uz$U),
    d           = ncol(Uz$U),
    pit_min_ks  = if(is.null(uni_tab)) NA_real_ else min(uni_tab$ks_p,  na.rm = TRUE),
    pit_min_cvm = if(is.null(uni_tab)) NA_real_ else min(uni_tab$cvm_p, na.rm = TRUE),
    pair_T1_min_p = pair_p_min
  )
})

print(gof_results)

## ------------------------------ Daily simulation for 2022 (medoid) -----------------------------------
# We simulate from the state-specific vine and invert with the (medoid, state) margin model.

# Helper: get margin model for (medoid, state)
get_medoid_model <- function(st){
  mm <- marginals_tbl %>%
    filter(codigo_estacion == medoid_station, States == st) %>%
    pull(model)
  if(length(mm) == 0) NULL else mm[[1]]
}

# Build daily simulation tibble for MJJA 2022
dates_2022 <- df_zone %>%
  filter(Year == test_year, month(dates) %in% months_MJJA,
         codigo_estacion == medoid_station) %>%
  pull(dates) %>% unique() %>% sort()

sim_list <- vector("list", length(dates_2022))
for(i in seq_along(dates_2022)){
  d0 <- dates_2022[i]
  st <- df_zone %>%
    filter(dates == d0, codigo_estacion == medoid_station) %>%
    pull(States) %>% unique()
  if(length(st) == 0) next
  
  vm <- vine_models[[as.character(st)]]
  med_model <- get_medoid_model(st)
  if(is.null(med_model)){
    sim_list[[i]] <- tibble(dates = d0, simulated = rep(NA_real_, M_per_day))
    next
  }
  
  # draw U for the medoid coordinate
  if(!is.null(vm) && !is.null(vm$vine)){
    # sample full d-dim vector then pick medoid's coordinate if present
    Udraw <- RVineSim(M_per_day, vm$vine)
    # map medoid to column index used in vine fit (if present)
    med_col <- match(as.character(medoid_station), vm$stations)
    u_vec <- if(!is.na(med_col)) Udraw[, med_col] else stats::runif(M_per_day)
  }else{
    # fallback: independent Uniform(0,1) if vine not available for this state
    u_vec <- stats::runif(M_per_day)
  }
  
  # invert to precipitation using the medoid state-conditional margin
  r_vec <- vapply(u_vec, Q_hurdle, numeric(1), model = med_model)
  sim_list[[i]] <- tibble(dates = rep(d0, M_per_day), simulated = r_vec)
}
sim_results <- bind_rows(sim_list) %>% arrange(dates)

## ------------------------------ Observations (medoid) and summaries ----------------------------------
obs_medoid <- df_zone %>%
  filter(codigo_estacion == medoid_station) %>%
  select(dates, precip_orig = precipitation)

# Daily 90% interval + mean of simulations
sim_daily_summary <- sim_results %>%
  group_by(dates) %>%
  summarise(
    q05  = quantile(simulated, 0.05, na.rm = TRUE),
    q95  = quantile(simulated, 0.95, na.rm = TRUE),
    mean = mean(simulated, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(obs_medoid, by = "dates")

# Plot: daily ribbon (MJJA 2022)
p_daily <- sim_daily_summary %>%
  filter(year(dates) == test_year, month(dates) %in% months_MJJA) %>%
  ggplot(aes(x = dates)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "skyblue", alpha = 0.35) +
  geom_line(aes(y = mean), color = "steelblue") +
  geom_line(aes(y = precip_orig), color = "orange") +
  labs(title = sprintf("Observed vs simulated precipitation – Zone %s medoid (MJJA %d)", 
                       selected_zone, test_year),
       x = NULL, y = "Precipitation (mm)") +
  theme_bw()
print(p_daily)

## ------------------------------ Monthly boxplots (MJJA 2022) -----------------------------------------
sim_month_2022 <- sim_results %>%
  mutate(month = factor(month(dates),
                        levels = months_MJJA,
                        labels = c("May","Jun","Jul","Aug")))

obs_month_2022 <- obs_medoid %>%
  filter(year(dates) == test_year, month(dates) %in% months_MJJA) %>%
  mutate(month = factor(month(dates), levels = months_MJJA,
                        labels = c("May","Jun","Jul","Aug"))) %>%
  group_by(month) %>%
  summarise(obs_mean = mean(precip_orig, na.rm = TRUE), .groups = "drop")

p_month <- ggplot(sim_month_2022, aes(x = month, y = simulated)) +
  geom_boxplot(fill = "grey85", color = "grey30", outlier.size = 0.4) +
  geom_point(data = obs_month_2022, aes(y = obs_mean), color = "orange", size = 3) +
  geom_line(data = obs_month_2022, aes(y = obs_mean, group = 1), color = "orange", linewidth = 0.9) +
  labs(title = sprintf("Monthly distributions of daily simulations – Zone %s medoid (%d)", 
                       selected_zone, test_year),
       x = "Month (MJJA)", y = "Precipitation (mm)") +
  tema_box
print(p_month)

## ------------------------------ Validation metrics (AW, POC, AAD) ------------------------------------
# Definitions used here (daily, aggregated over MJJA 2022):
#  - AW: mean width of the central 90% interval, E[q95 - q05]
#  - POC: % of observed daily values falling inside [q05, q95]
#  - AAD: mean asymmetry degree, E[ |(q95 + q05)/2 - median| / ((q95 - q05)/2) ]
#         (bounded in [0,∞); 0 means perfectly symmetric)

# Compute daily quantiles incl. median
sim_q_daily <- sim_results %>%
  filter(year(dates) == test_year, month(dates) %in% months_MJJA) %>%
  group_by(dates) %>%
  summarise(
    q05 = quantile(simulated, 0.05, na.rm = TRUE),
    q50 = quantile(simulated, 0.50, na.rm = TRUE),
    q95 = quantile(simulated, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(obs_medoid, by = "dates")

AW  <- mean(sim_q_daily$q95 - sim_q_daily$q05, na.rm = TRUE)
POC <- mean((sim_q_daily$precip_orig >= sim_q_daily$q05) &
              (sim_q_daily$precip_orig <= sim_q_daily$q95), na.rm = TRUE) * 100
AAD <- mean(abs(((sim_q_daily$q95 + sim_q_daily$q05)/2 - sim_q_daily$q50) /
                  pmax((sim_q_daily$q95 - sim_q_daily$q05)/2, 1e-6)), na.rm = TRUE)

metrics_tbl <- tibble(Zone = selected_zone, AW_mm = AW, POC_pct = POC, AAD = AAD)
print(metrics_tbl)

## ------------------------------ Optional: conditional PIT histograms (TRAIN) --------------------------
# Emission diagnostics per state (using the state-conditional margin above)
pit_train_zone <- df_zone %>%
  filter(is_train) %>%
  left_join(marginals_tbl, by = c("codigo_estacion","States")) %>%
  mutate(U_PIT = purrr::map2_dbl(precipitation, model, ~ F_hurdle(.x, .y)))

ggplot(pit_train_zone, aes(x = pmin(pmax(U_PIT, 0), 1))) +
  geom_histogram(bins = 20, color = "grey25", fill = "steelblue") +
  facet_wrap(~ States, ncol = 3) +
  labs(title = sprintf("Conditional PIT by HMM state (TRAIN) – Zone %s", selected_zone),
       x = expression(U[PIT]), y = "Count") +
  theme_bw()

##########################################################################################################
## Notes:
## 1) This script uses state-conditional margins with a zero-inflated hurdle and ECDF+GPD tail, which
##    directly addresses the PIT distortions observed in earlier diagnostics.
## 2) The vine family set is expanded (NA => all families incl. rotations) and truncated at Tree 2 to
##    reduce overfitting/instability in high dimension; selection is BIC-based by default.
## 3) The Rosenblatt-based p-values (pit_min_ks/cvm) and pairwise GoF (Tree-1) provide interpretable,
##    scalable goodness-of-fit checks for the dependence model.
## 4) Daily simulations for 2022 are performed by state using the fitted vine (when available) and are
##    inverted with the medoid’s (state-conditional) margin; M_per_day = 100 as requested.
## 5) AW/POC/AAD metrics are computed on daily scale over MJJA 2022; you can also aggregate by month.
##########################################################################################################


########################################################################################################
## BLOCK V — Multivariate verification where dependence matters
## Goals:
##  (A) Compute Energy Score (ES) and Variogram Score (VS) on a FIXED station set per state (d >= 2).
##  (B) Show skill on spatial aggregates: areal rainfall CRPS (dependence affects the aggregate).
##  (C) Event-based skill where joint behaviour matters: 
##      - ">= k wet stations" (Brier)
##      - joint exceedance in pairs (Brier; tail coherence)
##
## Assumptions:
##  - You already ran BLOCK 1 (hybrid margins + vine fitting), so the following exist:
##    datos_base_l, df_zone, marginals_tbl, vine_models, selected_zone, months_MJJA, test_year,
##    and the functions Q_hurdle(), fit_margin_state(), plus packages loaded there.
##  - HMM state S_t is common across stations at a given date (usual in nHMM workflows).
##  - We evaluate MJJA in 'test_year' only, for the selected zone.
########################################################################################################

pacman::p_load(scoringRules, forecast, dplyr, tidyr, purrr, tibble, lubridate, ggplot2)

## ------------------------------ Helper: fixed station set per state -----------------------------------
# Build a fixed set of stations for each state 's' to keep dimension d constant across days in TEST.
# We intersect: stations in the state-specific vine, stations with a fitted margin (station,state),
# and stations sufficiently complete in TEST days of that state (>= min_days_frac of the state's dates).

get_fixed_station_set <- function(state_id, min_days_frac = 0.95){
  # (i) vine stations (order used when fitting the vine)
  vm <- vine_models[[as.character(state_id)]]
  if (is.null(vm) || is.null(vm$vine) || is.null(vm$stations)) return(character(0))
  S_vine <- as.character(vm$stations)
  
  # (ii) stations with (station,state) margin
  S_marg <- marginals_tbl %>% 
    filter(States == state_id) %>% 
    pull(codigo_estacion) %>% 
    as.character() %>% unique()
  
  # (iii) stations well observed in TEST days of this state
  test_dates_s <- df_zone %>% 
    filter(Year == test_year, month(dates) %in% months_MJJA, States == state_id) %>%
    distinct(dates) %>% pull()
  if (length(test_dates_s) == 0) return(character(0))
  
  keep_tbl <- df_zone %>% 
    filter(Year == test_year, month(dates) %in% months_MJJA, States == state_id) %>%
    group_by(codigo_estacion) %>%
    summarise(n_days = sum(is.finite(precipitation)), .groups = "drop") %>%
    mutate(frac = n_days / length(test_dates_s)) %>%
    filter(frac >= min_days_frac)
  S_keep <- as.character(keep_tbl$codigo_estacion)
  
  # Final intersection. Keep the vine order.
  S <- S_vine[S_vine %in% intersect(S_marg, S_keep)]
  S
}

## ------------------------------ Helper: state-wise model table for a station set ----------------------
# Return a table (station,state) -> hybrid margin model in the exact order of 'station_order'.

get_models_for_state <- function(state_id, station_order){
  out <- marginals_tbl %>%
    filter(States == state_id, codigo_estacion %in% station_order) %>%
    mutate(codigo_estacion = as.character(codigo_estacion)) %>%
    arrange(match(codigo_estacion, station_order))
  # Safety check: ensure 1 model per station
  if (nrow(out) != length(station_order)) return(NULL)
  out
}

## ------------------------------ Helper: build observation vector per day ------------------------------
# Return y (numeric vector of length d) for the given date/state and station order; NULL if any missing.

get_y_vec <- function(date, state_id, station_order){
  day <- df_zone %>% 
    filter(dates == date, States == state_id, codigo_estacion %in% station_order) %>%
    mutate(codigo_estacion = as.character(codigo_estacion)) %>%
    arrange(match(codigo_estacion, station_order))
  if (nrow(day) != length(station_order)) return(NULL)
  y <- as.numeric(day$precipitation)
  if (any(!is.finite(y))) return(NULL)
  y
}

## ------------------------------ Helper: simulate M×d ensemble under vine / independence ---------------
# Both return a numeric matrix with nrow = M (members), ncol = d (in 'station_order').

simulate_ensemble_vine <- function(state_id, station_order, M = 500){
  vm <- vine_models[[as.character(state_id)]]
  if (is.null(vm) || is.null(vm$vine)) return(NULL)
  # Draw U on the full vine then select columns that match 'station_order'
  U_full <- VineCopula::RVineSim(M, vm$vine)            # matrix [M × d_vine]
  # Map vine columns to station names; the fitting code saved them in 'stations'
  vine_names <- as.character(vm$stations)
  idx <- match(station_order, vine_names)
  if (any(is.na(idx))) return(NULL)
  U <- U_full[, idx, drop = FALSE]                      # [M × d]
  # Map to precipitation using hybrid margins per station
  mods <- get_models_for_state(state_id, station_order)
  if (is.null(mods)) return(NULL)
  out <- matrix(NA_real_, nrow = M, ncol = length(station_order))
  for (j in seq_along(station_order)){
    mod_j <- mods$model[[j]]
    out[, j] <- vapply(U[, j], Q_hurdle, numeric(1), model = mod_j)
  }
  colnames(out) <- station_order
  out
}

simulate_ensemble_indep <- function(state_id, station_order, M = 500){
  d <- length(station_order)
  if (d == 0) return(NULL)
  U <- matrix(stats::runif(M * d), nrow = M, ncol = d)
  mods <- get_models_for_state(state_id, station_order)
  if (is.null(mods)) return(NULL)
  out <- matrix(NA_real_, nrow = M, ncol = d)
  for (j in seq_len(d)){
    mod_j <- mods$model[[j]]
    out[, j] <- vapply(U[, j], Q_hurdle, numeric(1), model = mod_j)
  }
  colnames(out) <- station_order
  out
}

## ------------------------------ Helper: day-level scores (ES/VS + areal CRPS + events) ----------------
# Returns a one-row tibble of metrics for the date; NULL if cannot be computed.

score_one_day <- function(date, state_id, station_order, M = 500, alpha_vs = 0.5,
                          k_wet_frac = 0.2, thr_joint = 10){
  # observation vector (d >= 2 required for multivariate ES/VS)
  y <- get_y_vec(date, state_id, station_order)
  if (is.null(y)) return(NULL)
  d <- length(y)
  if (d < 2) return(NULL)
  
  # ensembles (both must be fully numeric with same dimensions)
  ens_v <- simulate_ensemble_vine(state_id, station_order, M = M)
  ens_i <- simulate_ensemble_indep(state_id, station_order, M = M)
  if (is.null(ens_v) || is.null(ens_i)) return(NULL)
  if (!is.numeric(ens_v) || !is.numeric(ens_i)) return(NULL)
  if (!all(dim(ens_v) == dim(ens_i))) return(NULL)
  if (ncol(ens_v) != d) return(NULL)
  
  # ---- Multivariate scores (scoringRules expects: y length d, dat [M × d]) ----
  ES_v <- tryCatch(scoringRules::es_sample(y = y, dat = ens_v), error = function(e) NA_real_)
  ES_i <- tryCatch(scoringRules::es_sample(y = y, dat = ens_i), error = function(e) NA_real_)
  W    <- matrix(1, d, d); diag(W) <- 0
  VS_v <- tryCatch(scoringRules::vs_sample(y = y, dat = ens_v, w = W, p = alpha_vs), error = function(e) NA_real_)
  VS_i <- tryCatch(scoringRules::vs_sample(y = y, dat = ens_i, w = W, p = alpha_vs), error = function(e) NA_real_)
  
  # ---- Areal rainfall (spatial mean) – univariate CRPS where dependence matters ----
  y_areal   <- mean(y)
  draw_areal_v <- rowMeans(ens_v)
  draw_areal_i <- rowMeans(ens_i)
  CRPS_areal_v <- tryCatch(scoringRules::crps_sample(y_areal, dat = as.numeric(draw_areal_v)),
                           error = function(e) NA_real_)
  CRPS_areal_i <- tryCatch(scoringRules::crps_sample(y_areal, dat = as.numeric(draw_areal_i)),
                           error = function(e) NA_real_)
  
  # ---- Event: ">= k wet stations" (k based on fraction of d) – Brier score ----
  k <- max(2L, floor(k_wet_frac * d))
  obs_k  <- as.integer(sum(y > 0) >= k)
  pr_k_v <- mean(rowSums(ens_v > 0) >= k)
  pr_k_i <- mean(rowSums(ens_i > 0) >= k)
  BS_k_v <- (pr_k_v - obs_k)^2
  BS_k_i <- (pr_k_i - obs_k)^2
  
  # ---- Event: pairwise joint exceedance (>= thr_joint mm) – average Brier over first P pairs ----
  # Choose up to 10 lexicographic pairs to keep it light and consistent across days.
  Pairs <- t(utils::combn(seq_len(d), 2))
  P <- min(nrow(Pairs), 10L)
  if (P > 0){
    Pairs <- Pairs[seq_len(P), , drop = FALSE]
    obs_pair <- apply(Pairs, 1, function(ix) as.integer(all(y[ix] > thr_joint)))
    pr_pair_v <- apply(Pairs, 1, function(ix) mean(ens_v[, ix[1]] > thr_joint & ens_v[, ix[2]] > thr_joint))
    pr_pair_i <- apply(Pairs, 1, function(ix) mean(ens_i[, ix[1]] > thr_joint & ens_i[, ix[2]] > thr_joint))
    BS_pair_v <- mean((pr_pair_v - obs_pair)^2)
    BS_pair_i <- mean((pr_pair_i - obs_pair)^2)
  } else {
    BS_pair_v <- NA_real_; BS_pair_i <- NA_real_
  }
  
  tibble(
    date = as.Date(date), state = as.integer(state_id), d = d,
    ES_v = ES_v, ES_i = ES_i, VS_v = VS_v, VS_i = VS_i,
    CRPSA_v = CRPS_areal_v, CRPSA_i = CRPS_areal_i,
    BS_k_v = BS_k_v, BS_k_i = BS_k_i,
    BS_pair_v = BS_pair_v, BS_pair_i = BS_pair_i
  )
}

## ------------------------------ Main loop: compute scores over MJJA-TEST ------------------------------
verify_zone_multivar <- function(M = 500, alpha_vs = 0.5, min_days_frac = 0.95,
                                 k_wet_frac = 0.2, thr_joint = 10, verbose = TRUE){
  
  # Candidate days in TEST (MJJA) for the selected zone
  cal <- df_zone %>%
    filter(Year == test_year, month(dates) %in% months_MJJA) %>%
    distinct(dates, States) %>%
    arrange(dates)
  
  if (nrow(cal) == 0){
    message("No TEST dates in zone.")
    return(NULL)
  }
  
  # Build fixed station sets per state
  states_in_cal <- sort(unique(cal$States))
  S_fixed <- purrr::map(states_in_cal, ~get_fixed_station_set(.x, min_days_frac = min_days_frac))
  names(S_fixed) <- as.character(states_in_cal)
  
  if (verbose){
    cat(sprintf("Zone %s — fixed station sets per state (d>=2 required):\n", selected_zone))
    for (s in states_in_cal){
      cat(sprintf("  State %s: d_fixed = %d\n", s, length(S_fixed[[as.character(s)]])))
    }
  }
  
  # Score day by day (skip if d<2, missing obs, or ensembles cannot be built)
  res <- vector("list", nrow(cal))
  pb  <- txtProgressBar(min = 0, max = nrow(cal), style = 3)
  for (ii in seq_len(nrow(cal))){
    dt <- cal$dates[ii]; st <- cal$States[ii]
    S  <- S_fixed[[as.character(st)]]
    if (length(S) >= 2){
      res[[ii]] <- score_one_day(dt, st, S, M = M, alpha_vs = alpha_vs,
                                 k_wet_frac = k_wet_frac, thr_joint = thr_joint)
    } else {
      res[[ii]] <- NULL
    }
    setTxtProgressBar(pb, ii)
  }
  close(pb)
  
  out <- bind_rows(res)
  if (nrow(out) == 0){
    message("No valid multivariate days (likely d<2 after fixing station sets or missing data).")
    return(NULL)
  }
  out
}

## ------------------------------ Run and summarise ------------------------------------------------------
set.seed(123)
scores_mv <- verify_zone_multivar(
  M = 500,            # ensemble size
  alpha_vs = 0.5,     # VS power
  min_days_frac = 0.95,
  k_wet_frac = 0.2,
  thr_joint = 10,
  verbose = TRUE
)

if (is.null(scores_mv)) {
  stop("Verification produced no rows. Check d_fixed per state and data completeness.")
}

# Basic summaries (vine better if ES/VS/CRPSA smaller; BS smaller is better)
summary_tbl <- scores_mv %>%
  summarise(
    n_days   = n(),
    ES_v     = mean(ES_v,     na.rm = TRUE), ES_i     = mean(ES_i,     na.rm = TRUE),
    VS_v     = mean(VS_v,     na.rm = TRUE), VS_i     = mean(VS_i,     na.rm = TRUE),
    CRPSA_v  = mean(CRPSA_v,  na.rm = TRUE), CRPSA_i  = mean(CRPSA_i,  na.rm = TRUE),
    BS_k_v   = mean(BS_k_v,   na.rm = TRUE), BS_k_i   = mean(BS_k_i,   na.rm = TRUE),
    BS_pair_v= mean(BS_pair_v,na.rm = TRUE), BS_pair_i= mean(BS_pair_i,na.rm = TRUE)
  )
print(summary_tbl)

# Diebold–Mariano on areal CRPS (lower is better for vine; H=1 one-step ahead)
dm_CRPSA_p <- tryCatch(
  forecast::dm.test(e1 = scores_mv$CRPSA_v, e2 = scores_mv$CRPSA_i, alternative = "less", h = 1)$p.value,
  error = function(e) NA_real_
)
cat(sprintf("DM p-value (vine better on areal CRPS if p < 0.05): %s\n",
            ifelse(is.na(dm_CRPSA_p), "NA", format(dm_CRPSA_p, digits = 3))))

## ------------------------------ Visual diagnostics ----------------------------------------------------
# (1) Time series of score gains (indep − vine) – positive values favour vine
p_es <- ggplot(scores_mv, aes(x = date, y = ES_i - ES_v)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() +
  labs(title = sprintf("Daily ES gain: independence − vine (Zone %s)", selected_zone),
       y = "ΔES (indep − vine)", x = NULL) + theme_bw()
print(p_es)

p_vs <- ggplot(scores_mv, aes(x = date, y = VS_i - VS_v)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() +
  labs(title = sprintf("Daily VS gain: independence − vine (Zone %s)", selected_zone),
       y = "ΔVS (indep − vine)", x = NULL) + theme_bw()
print(p_vs)

p_crpsA <- ggplot(scores_mv, aes(x = date, y = CRPSA_i - CRPSA_v)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() +
  labs(title = sprintf("Daily areal-CRPS gain: independence − vine (Zone %s)", selected_zone),
       y = "ΔCRPS_areal (indep − vine)", x = NULL) + theme_bw()
print(p_crpsA)

# (2) Boxplots of gains across MJJA
g_box <- scores_mv %>%
  transmute(`ΔES` = ES_i - ES_v,
            `ΔVS` = VS_i - VS_v,
            `ΔCRPS_areal` = CRPSA_i - CRPSA_v) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "gain")
p_box <- ggplot(g_box, aes(x = metric, y = gain)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(fill = "grey85") +
  labs(title = sprintf("Score gains (independence − vine), Zone %s", selected_zone),
       x = NULL, y = "gain (>0 favours vine)") + theme_bw()
print(p_box)

# (3) Event Brier gains
p_bs <- scores_mv %>%
  transmute(date,
            `ΔBS (>=k wet)` = BS_k_i - BS_k_v,
            `ΔBS (pair>=thr)` = BS_pair_i - BS_pair_v) %>%
  pivot_longer(-date, names_to = "event", values_to = "gain") %>%
  ggplot(aes(x = date, y = gain)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() + facet_wrap(~event, ncol = 1, scales = "free_y") +
  labs(title = sprintf("Brier gains: independence − vine (Zone %s)", selected_zone),
       x = NULL, y = "ΔBrier (>0 favours vine)") + theme_bw()
print(p_bs)

########################################################################################################
## Notes:
## - The crucial fix vs previous attempts is the *fixed station set per state* (d constant, d>=2),
##   guaranteeing ES/VS are well-defined and dimension-matched for every scored day.
## - Areal CRPS uses the spatial mean of ensemble draws; this is where dependence changes the 
##   predictive distribution of the aggregate and can reveal vine skill.
## - Event-based Brier scores target joint behaviour (co-occurrence of wet days, joint exceedance).
## - If 'scores_mv' is empty, usually either no state reaches d>=2 at the required completeness
##   or the vine had not been fitted for that state in the zone.
########################################################################################################



## ------------------------------ Figure export helpers -----------------------------------------------
# Create output dirs
out_dir_main <- file.path(getwd(), "fig_main")       # figures for the paper
out_dir_supp <- file.path(getwd(), "fig_supp")       # figures for the supplement
dir.create(out_dir_main, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_supp, showWarnings = FALSE, recursive = TRUE)

# Generic saver: writes PDF (vector) + PNG (hi-res)
save_fig <- function(plot_obj, base_name, path_main = out_dir_main, path_supp = out_dir_supp,
                     width = 9, height = 6, dpi = 700, to_main = TRUE, to_supp = TRUE){
  if (inherits(plot_obj, "ggplot")) {
    if (to_main) {
      ggsave(filename = file.path(path_main, paste0(base_name, ".pdf")),
             plot = plot_obj, width = width, height = height, units = "in", device = cairo_pdf)
      ggsave(filename = file.path(path_main, paste0(base_name, ".png")),
             plot = plot_obj, width = width, height = height, units = "in", dpi = dpi)
    }
    if (to_supp) {
      ggsave(filename = file.path(path_supp, paste0(base_name, ".pdf")),
             plot = plot_obj, width = width, height = height, units = "in", device = cairo_pdf)
      ggsave(filename = file.path(path_supp, paste0(base_name, ".png")),
             plot = plot_obj, width = width, height = height, units = "in", dpi = dpi)
    }
  } else {
    message("Object is not a ggplot; skipping: ", base_name)
  }
}

## ------------------------------ Capture PIT figure to an object -------------------------------------
# In tu bloque 1, el PIT se dibuja “en línea”. Guardémoslo en un objeto para exportarlo.
p_pit <- ggplot(pit_train_zone, aes(x = pmin(pmax(U_PIT, 0), 1))) +
  geom_histogram(bins = 20, color = "grey25", fill = "steelblue") +
  facet_wrap(~ States, ncol = 3) +
  labs(title = sprintf("Conditional PIT by HMM state (TRAIN) – Zone %s", selected_zone),
       x = expression(U[PIT]), y = "Count") +
  theme_bw()

## ------------------------------ Save main figures ----------------------------------------------------
# Paper-ready (elige las que van al cuerpo principal)
save_fig(p_daily,  base_name = sprintf("Z%02d_MJJA%d_daily_ribbon", selected_zone, test_year))
save_fig(p_month,  base_name = sprintf("Z%02d_MJJA%d_monthly_boxplots", selected_zone, test_year))
save_fig(p_pit,    base_name = sprintf("Z%02d_train_conditional_PIT", selected_zone, test_year),
         width = 10, height = 7)

# Si ejecutaste el bloque multivariante y tienes estas figuras:
if (exists("p_es"))     save_fig(p_es,     base_name = sprintf("Z%02d_ES_gain", selected_zone), height = 4.5)
if (exists("p_vs"))     save_fig(p_vs,     base_name = sprintf("Z%02d_VS_gain", selected_zone), height = 4.5)
if (exists("p_crpsA"))  save_fig(p_crpsA,  base_name = sprintf("Z%02d_CRPS_areal_gain", selected_zone), height = 4.5)
if (exists("p_box"))    save_fig(p_box,    base_name = sprintf("Z%02d_gain_boxplots", selected_zone), height = 5.5)
if (exists("p_bs"))     save_fig(p_bs,     base_name = sprintf("Z%02d_Brier_gains", selected_zone), height = 5.5)

## ------------------------------ (Opcional) Exportar tablas clave ------------------------------------
# GoF por estado y resumen multivariante a CSV (suplemento)
if (exists("gof_results")) {
  readr::write_csv(gof_results,
                   file.path(out_dir_supp, sprintf("Z%02d_gof_results.csv", selected_zone)))
}
if (exists("scores_mv")) {
  readr::write_csv(scores_mv,
                   file.path(out_dir_supp, sprintf("Z%02d_scores_multivar_MJJA%d.csv", selected_zone, test_year)))
}


# PART 2


##########################################################################################################
## Hybrid nHMM — State-conditional margins (zero-inflated + ECDF + GPD) — R-vine copula (zone-wise)
## Simulation for 2022 (MJJA), diagnostics (PIT, Rosenblatt, pairwise GoF), and multivariate validation.
## Author: (MAURICIO HERRERA-MARÍN)  
##########################################################################################################
## NOTE:
##   This script assumes that data objects and helper functions used downstream are available.
##   If you have a preparatory script (e.g., "PART 1),
##   run it beforehand to load any missing data/functions.

## ============================== 0) Setup & Libraries ===================================================

## Run from the repository root, e.g.:
## Rscript Scripts/Hybrid_nHMM_Copula_Simulation_Validation_2022.R
## The input CSV is shipped inside the ZIP under Data/.
base_path <- normalizePath(getwd(), mustWork = TRUE)
data_dir <- file.path(base_path, "Data")
data_zip_path <- file.path(data_dir, "mapa_estaciones_precipitaciones_clean.csv.zip")
stations_csv_path <- file.path(data_dir, "mapa_estaciones_precipitaciones_clean.csv")
if (!dir.exists(data_dir)) {
  stop("Data/ not found. Run the script from the root of the GitHub repository.")
}
if (!file.exists(stations_csv_path)) {
  if (!file.exists(data_zip_path)) {
    stop("Missing input ZIP: ", data_zip_path)
  }
  zip_entries <- utils::unzip(data_zip_path, list = TRUE)$Name
  csv_entries <- zip_entries[grepl("\\.csv$", zip_entries, ignore.case = TRUE)]
  if (length(csv_entries) != 1L) {
    stop("Expected exactly one CSV in ", data_zip_path,
         "; found ", length(csv_entries), ". Check the archive manually.")
  }
  utils::unzip(data_zip_path, files = csv_entries, exdir = data_dir,
               junkpaths = TRUE, overwrite = FALSE)
  extracted_csv <- file.path(data_dir, basename(csv_entries))
  if (!file.exists(stations_csv_path) && file.exists(extracted_csv)) {
    if (!file.copy(extracted_csv, stations_csv_path)) {
      stop("Unable to copy extracted CSV to expected filename: ", stations_csv_path)
    }
  }
}
if (!file.exists(stations_csv_path)) {
  stop("Input CSV was not found after extraction: ", stations_csv_path)
}
datos_base_l <- read.csv(stations_csv_path)

## Packages
pacman::p_load(
  dplyr, tidyr, purrr, tibble, data.table,
  lubridate, ggplot2, stringr,
  VineCopula, copula, evd,
  parallel, goftest, ADGofTest,
  scoringRules, forecast, readr
)

## Plot theming
tema_box <- theme_bw() +
  theme(
    plot.title       = element_text(size = 16, face = "bold"),
    axis.title       = element_text(size = 13),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey92")
  )

## Output directories (main and supplement) + legacy fallback "supp_figures" for compatibility
out_dir_main <- file.path(getwd(), "fig_main")
out_dir_supp <- file.path(getwd(), "fig_supp")
dir.create(out_dir_main, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_supp, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(getwd(), "supp_figures"), showWarnings = FALSE, recursive = TRUE)  # legacy name

## Generic figure saver: PDF (vector) + PNG (hi-res). Uses cairo if available.
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
save_fig <- function(plot_obj, base_name, path_main = out_dir_main, path_supp = out_dir_supp,
                     width = 9, height = 6, dpi = 700, to_main = TRUE, to_supp = TRUE){
  if (!inherits(plot_obj, "ggplot")) {
    message("Object is not a ggplot; skipping: ", base_name)
    return(invisible(NULL))
  }
  if (to_main) {
    ggsave(filename = file.path(path_main, paste0(base_name, ".pdf")),
           plot = plot_obj, width = width, height = height, units = "in", device = pdf_device)
    ggsave(filename = file.path(path_main, paste0(base_name, ".png")),
           plot = plot_obj, width = width, height = height, units = "in", dpi = dpi)
  }
  if (to_supp) {
    ggsave(filename = file.path(path_supp, paste0(base_name, ".pdf")),
           plot = plot_obj, width = width, height = height, units = "in", device = pdf_device)
    ggsave(filename = file.path(path_supp, paste0(base_name, ".png")),
           plot = plot_obj, width = width, height = height, units = "in", dpi = dpi)
  }
}

## Small utilities
clamp01 <- function(u, eps = 1e-12) pmin(pmax(u, eps), 1 - eps)

## Core run configuration
selected_zone  <- 5            # <- change zone id here
M_per_day      <- 300          # daily simulations for medoid visuals (100 recommended)
train_end      <- as.Date("2021-08-31")
months_MJJA    <- 5:8
test_year      <- 2022

## Vine fitting options
vine_familyset <- NA           # NA => allow all (including rotations)
vine_selcrit   <- "BIC"        # penalize complexity
vine_truncate  <- 2            # truncate vine at Tree 2 (robust for high-dim)
min_rows_vine  <- 200          # minimum rows per state to fit vine

## Tail model options
u_star         <- 0.90         # GPD threshold on positive precipitation
min_excess_gpd <- 25           # minimum exceedances to fit GPD

## ============================== 1) Parameters & Data Preflight ========================================

set.seed(42)

verify_zone_multivar <- function(M = 500, alpha_vs = 0.5, min_days_frac = 0.95,
                                 k_wet_frac = 0.2, thr_joint = 10, verbose = TRUE){
  
  # Candidate TEST days (MJJA) for the selected zone
  cal <- df_zone %>%
    filter(Year == test_year, month(dates) %in% months_MJJA) %>%
    distinct(dates, States) %>%
    arrange(dates)
  
  if (nrow(cal) == 0){
    message("No TEST dates in zone.")
    return(NULL)
  }
  
  # Build fixed station sets per state
  states_in_cal <- sort(unique(cal$States))
  S_fixed <- purrr::map(states_in_cal, ~get_fixed_station_set(.x, min_days_frac = min_days_frac))
  names(S_fixed) <- as.character(states_in_cal)
  
  if (verbose){
    cat(sprintf("Zone %s — fixed station sets per state (d>=2 required):\n", selected_zone))
    for (s in states_in_cal){
      cat(sprintf("  State %s: d_fixed = %d\n", s, length(S_fixed[[as.character(s)]])))
    }
  }
  
  # Score day-by-day (skip if d<2, missing obs, or ensembles cannot be built)
  res <- vector("list", nrow(cal))
  pb  <- txtProgressBar(min = 0, max = nrow(cal), style = 3)
  for (ii in seq_len(nrow(cal))){
    dt <- cal$dates[ii]; st <- cal$States[ii]
    S  <- S_fixed[[as.character(st)]]
    if (length(S) >= 2){
      res[[ii]] <- score_one_day(dt, st, S, M = M, alpha_vs = alpha_vs,
                                 k_wet_frac = k_wet_frac, thr_joint = thr_joint)
    } else {
      res[[ii]] <- NULL
    }
    setTxtProgressBar(pb, ii)
  }
  close(pb)
  
  out <- bind_rows(res)
  if (nrow(out) == 0){
    message("No valid multivariate days (likely d<2 after fixing station sets or missing data).")
    return(NULL)
  }
  out
}

## Subset to selected zone & detect medoid
df_zone <- datos_base_l %>% filter(labels == selected_zone)
medoid_station <- df_zone %>% filter(medoids == 1) %>% pull(codigo_estacion) %>% unique()
if (length(medoid_station) != 1) stop("Medoid not found or not unique in selected zone.")

## Multivariate scores for TEST (zone-level verification)
scores_mv <- verify_zone_multivar(
  M = 500,            # ensemble size
  alpha_vs = 0.5,     # VS power
  min_days_frac = 0.95,
  k_wet_frac = 0.2,
  thr_joint = 10,
  verbose = TRUE
)

## Required data columns check
stopifnot(all(c("dates","codigo_estacion","precipitation","labels","States","medoids") %in% names(datos_base_l)))

## Prepare base data frame with helpers
datos_base_l <- datos_base_l %>%
  mutate(dates = as.Date(dates),
         Year  = year(dates),
         is_train = dates <= train_end)

## Select zone & medoid (recompute after adding helpers)
df_zone <- datos_base_l %>% filter(labels == selected_zone)
medoid_station <- df_zone %>% filter(medoids == 1) %>% pull(codigo_estacion) %>% unique()
if (length(medoid_station) != 1) stop("Medoid not found or not unique in selected zone.")

## ============================== 2) State-conditional Margin Model =====================================

# Fit per (station, state) margin: zero atom + ECDF body + GPD tail
fit_margin_state <- function(x, u_star = 0.90, min_excess = 25){
  x <- x[is.finite(x)]
  if (length(x) < 40L){
    p0   <- mean(x <= 0, na.rm = TRUE)
    xpos <- sort(x[x > 0])
    Fe   <- if (length(xpos) > 0) stats::ecdf(xpos) else NULL
    return(list(p0 = p0, Fe = Fe, x_sorted = xpos,
                thr = Inf, Fu = 1, gpd = NULL))
  }
  p0   <- mean(x <= 0, na.rm = TRUE)
  xpos <- sort(x[x > 0])
  if (length(xpos) < 10L){
    Fe <- NULL
    return(list(p0 = p0, Fe = Fe, x_sorted = xpos,
                thr = Inf, Fu = 1, gpd = NULL))
  }
  Fe   <- stats::ecdf(xpos)
  thr  <- stats::quantile(xpos, probs = u_star, na.rm = TRUE, type = 7)
  Fu   <- Fe(thr)
  exc  <- xpos[xpos > thr] - thr
  gpd  <- NULL
  if (length(exc) >= min_excess){
    gpd <- tryCatch(evd::fpot(xpos, threshold = thr, model = "gpd", std.err = FALSE),
                    error = function(e) NULL)
  }
  list(p0 = p0, Fe = Fe, x_sorted = xpos,
       thr = thr, Fu = Fu, gpd = gpd)
}

# Hybrid CDF (randomized PIT at zero)
F_hurdle <- function(r, model){
  if (is.null(model)) return(NA_real_)
  p0 <- model$p0
  if (r <= 0){
    return(runif(1, min = 0, max = p0))
  }
  Fe <- model$Fe
  if (is.null(Fe)) return(p0)
  if (r <= model$thr || !is.finite(model$thr) || is.null(model$gpd)){
    return(p0 + (1 - p0) * Fe(r))
  } else {
    p_ex  <- evd::pgpd(r - model$thr, loc = 0,
                       scale = model$gpd$estimate["scale"],
                       shape = model$gpd$estimate["shape"])
    return(p0 + (1 - p0) * (model$Fu + (1 - model$Fu) * p_ex))
  }
}

# Hybrid quantile (inverse CDF approximation)
Q_hurdle <- function(u, model){
  if (is.null(model)) return(NA_real_)
  u <- clamp01(u)  # avoid 0/1 boundary issues
  p0 <- model$p0
  if (u <= p0) return(0)
  if (is.null(model$Fe) || length(model$x_sorted) == 0) return(0)
  u_pos <- (u - p0) / (1 - p0)
  if (u_pos <= model$Fu || is.null(model$gpd) || !is.finite(model$thr)){
    return(as.numeric(stats::quantile(model$x_sorted, probs = u_pos, type = 7)))
  }
  p_tail <- (u_pos - model$Fu) / (1 - model$Fu)
  par <- model$gpd$estimate
  q_ex <- evd::qgpd(p_tail, loc = 0, scale = par["scale"], shape = par["shape"])
  return(as.numeric(model$thr + q_ex))
}

## Fit margins on TRAIN
marginals_tbl <- datos_base_l %>%
  filter(is_train) %>%
  group_by(codigo_estacion, States) %>%
  summarise(model = list(fit_margin_state(precipitation, u_star, min_excess_gpd)),
            .groups = "drop")

## ============================== 3) Build U-matrices & Fit Vine per State ==============================

calc_U_zone_state <- function(df_zone, state_id, marginals_tbl){
  df <- df_zone %>%
    filter(is_train, States == state_id) %>%
    select(dates, codigo_estacion, States, precipitation)
  if (nrow(df) == 0) return(NULL)
  
  dfU <- df %>%
    left_join(marginals_tbl, by = c("codigo_estacion","States")) %>%
    mutate(U = purrr::map2_dbl(precipitation, model, ~ F_hurdle(.x, .y))) %>%
    mutate(codigo_estacion = as.character(codigo_estacion)) %>%
    select(dates, codigo_estacion, U) %>%
    tidyr::pivot_wider(names_from = codigo_estacion, values_from = U) %>%
    arrange(dates) %>%
    suppressWarnings(tidyr::drop_na())
  
  if (nrow(dfU) == 0) return(NULL)
  list(U = as.matrix(dfU %>% select(-dates)),
       station_ids = colnames(dfU %>% select(-dates)))
}

vine_models <- vector("list", length = 5)
names(vine_models) <- as.character(sort(unique(df_zone$States)))

for (st in sort(unique(df_zone$States))){
  Uz <- calc_U_zone_state(df_zone, st, marginals_tbl)
  if (is.null(Uz) || nrow(Uz$U) < min_rows_vine){
    vine_models[[as.character(st)]] <- NULL
    next
  }
  message(sprintf("Fitting vine: Zone %s, State %s  (n=%d, d=%d)",
                  selected_zone, st, nrow(Uz$U), ncol(Uz$U)))
  vine_fit <- tryCatch(
    RVineStructureSelect(Uz$U,
                         familyset     = vine_familyset,
                         selectioncrit = vine_selcrit,
                         trunclevel    = vine_truncate),
    error = function(e) { message("  -> vine fit failed: ", e$message); NULL }
  )
  vine_models[[as.character(st)]] <- list(vine = vine_fit,
                                          stations = Uz$station_ids)
}

## ============================== 4) Vine GoF Diagnostics (TRAIN) =======================================

extract_T1_pairs <- function(rvine_obj){
  if (is.null(rvine_obj) || is.null(rvine_obj$Matrix)) return(NULL)
  d <- ncol(rvine_obj$Matrix)
  tibble::tibble(
    i      = rvine_obj$Matrix[d, (d-1):1],
    j      = diag(rvine_obj$Matrix)[(d-1):1],
    family = rvine_obj$family[(d-1):1, d],
    par    = rvine_obj$par[(d-1):1, d],
    par2   = rvine_obj$par2[(d-1):1, d]
  )
}

gof_results <- purrr::imap_dfr(vine_models, function(vm, st){
  if (is.null(vm) || is.null(vm$vine)) return(tibble(State = as.integer(st)))
  
  Uz <- calc_U_zone_state(df_zone, as.integer(st), marginals_tbl)
  if (is.null(Uz) || nrow(Uz$U) == 0) return(tibble(State = as.integer(st)))
  
  Rz <- tryCatch(VineCopula::RVinePIT(Uz$U, vm$vine), error = function(e) NULL)
  uni_tab <- if (!is.null(Rz)){
    ktest <- min(5L, ncol(Rz))
    purrr::map_dfr(seq_len(ktest), function(j){
      u <- Rz[, j]
      tibble::tibble(
        comp  = j,
        ks_p  = ks.test(u, "punif")$p.value,
        cvm_p = goftest::cvm.test(u, "punif")$p.value
      )
    })
  } else NULL
  
  pairs <- extract_T1_pairs(vm$vine)
  pair_p_min <- NA_real_
  if (!is.null(pairs) && nrow(pairs) > 0){
    pvals <- purrr::pmap_dbl(pairs, function(i, j, family, par, par2){
      ui <- Uz$U[, as.integer(i)]
      uj <- Uz$U[, as.integer(j)]
      tt <- tryCatch(
        VineCopula::BiCopGofTest(ui, uj, family = family, par = par, par2 = par2,
                                 method = "White", B = 200),
        error = function(e) NULL
      )
      if (is.null(tt)) NA_real_ else tt$p.value
    })
    pair_p_min <- suppressWarnings(min(pvals, na.rm = TRUE))
    if (!is.finite(pair_p_min)) pair_p_min <- NA_real_
  }
  
  tibble::tibble(
    State         = as.integer(st),
    n             = nrow(Uz$U),
    d             = ncol(Uz$U),
    pit_min_ks    = if (is.null(uni_tab)) NA_real_ else min(uni_tab$ks_p,  na.rm = TRUE),
    pit_min_cvm   = if (is.null(uni_tab)) NA_real_ else min(uni_tab$cvm_p, na.rm = TRUE),
    pair_T1_min_p = pair_p_min
  )
})

print(gof_results)



## ============================== 5) Medoid Simulation & Univariate Checks ==============================

# Helper: retrieve the fitted marginal model for the zone medoid at a given state
get_medoid_model <- function(st){
  mm <- marginals_tbl %>%
    filter(codigo_estacion == medoid_station, States == st) %>%
    pull(model)
  if (length(mm) == 0) NULL else mm[[1]]
}

# MJJA dates in test_year available at the medoid
dates_2022 <- df_zone %>%
  filter(Year == test_year, month(dates) %in% months_MJJA,
         codigo_estacion == medoid_station) %>%
  pull(dates) %>% unique() %>% sort()

# Daily medoid simulations (state-conditional, via vine if available; otherwise iid uniforms)
sim_list <- vector("list", length(dates_2022))
for (i in seq_along(dates_2022)){
  d0 <- dates_2022[i]
  st <- df_zone %>%
    filter(dates == d0, codigo_estacion == medoid_station) %>%
    pull(States) %>% unique()
  if (length(st) == 0) next
  
  vm <- vine_models[[as.character(st)]]
  med_model <- get_medoid_model(st)
  if (is.null(med_model)){
    sim_list[[i]] <- tibble(dates = d0, simulated = rep(NA_real_, M_per_day))
    next
  }
  
  if (!is.null(vm) && !is.null(vm$vine)){
    Udraw   <- RVineSim(M_per_day, vm$vine)
    med_col <- match(as.character(medoid_station), vm$stations)
    u_vec   <- if (!is.na(med_col)) Udraw[, med_col] else stats::runif(M_per_day)
  } else {
    u_vec <- stats::runif(M_per_day)
  }
  
  u_vec <- clamp01(u_vec) # safety
  r_vec <- vapply(u_vec, Q_hurdle, numeric(1), model = med_model)
  sim_list[[i]] <- tibble(dates = rep(d0, M_per_day), simulated = r_vec)
}
sim_results <- bind_rows(sim_list) %>% arrange(dates)

# Observed series at the medoid
obs_medoid <- df_zone %>%
  filter(codigo_estacion == medoid_station) %>%
  select(dates, precip_orig = precipitation)

## Daily ribbon (90% band) vs observed at the medoid
sim_daily_summary <- sim_results %>%
  group_by(dates) %>%
  summarise(
    q05  = quantile(simulated, 0.05, na.rm = TRUE),
    q95  = quantile(simulated, 0.95, na.rm = TRUE),
    mean = mean(simulated, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(obs_medoid, by = "dates")

p_daily <- sim_daily_summary %>%
  filter(year(dates) == test_year, month(dates) %in% months_MJJA) %>%
  ggplot(aes(x = dates)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "skyblue", alpha = 0.35) +
  geom_line(aes(y = mean), color = "steelblue", linewidth = 0.7) +
  geom_line(aes(y = precip_orig), color = "orange", linewidth = 1) +
  labs(title = sprintf("Observed vs simulated precipitation — Zone %s medoid (MJJA %d)",
                       selected_zone, test_year),
       x = NULL, y = "Precipitation (mm)") +
  theme_bw()
print(p_daily)

## Monthly distributions of daily simulations vs monthly observed mean (medoid)
sim_month_2022 <- sim_results %>%
  mutate(month = factor(month(dates),
                        levels = months_MJJA,
                        labels = c("May","Jun","Jul","Aug")))

obs_month_2022 <- obs_medoid %>%
  filter(year(dates) == test_year, month(dates) %in% months_MJJA) %>%
  mutate(month = factor(month(dates), levels = months_MJJA,
                        labels = c("May","Jun","Jul","Aug"))) %>%
  group_by(month) %>%
  summarise(obs_mean = mean(precip_orig, na.rm = TRUE), .groups = "drop")

p_month <- ggplot(sim_month_2022, aes(x = month, y = simulated)) +
  geom_boxplot(fill = "grey85", color = "grey30", outlier.size = 0.4) +
  geom_point(data = obs_month_2022, aes(y = obs_mean), color = "orange", size = 3) +
  geom_line(data = obs_month_2022, aes(y = obs_mean, group = 1), color = "orange", linewidth = 0.9) +
  labs(title = sprintf("Monthly distributions of daily simulations — Zone %s medoid (%d)",
                       selected_zone, test_year),
       x = "Month (MJJA)", y = "Precipitation (mm)") +
  tema_box
print(p_month)

## Simple interval diagnostics at the medoid (AW, POC, AAD)
sim_q_daily <- sim_results %>%
  filter(year(dates) == test_year, month(dates) %in% months_MJJA) %>%
  group_by(dates) %>%
  summarise(
    q05 = quantile(simulated, 0.05, na.rm = TRUE),
    q50 = quantile(simulated, 0.50, na.rm = TRUE),
    q95 = quantile(simulated, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(obs_medoid, by = "dates")

AW  <- mean(sim_q_daily$q95 - sim_q_daily$q05, na.rm = TRUE)
POC <- mean((sim_q_daily$precip_orig >= sim_q_daily$q05) &
              (sim_q_daily$precip_orig <= sim_q_daily$q95), na.rm = TRUE) * 100
AAD <- mean(abs(((sim_q_daily$q95 + sim_q_daily$q05) / 2 - sim_q_daily$q50) /
                  pmax((sim_q_daily$q95 - sim_q_daily$q05) / 2, 1e-6)), na.rm = TRUE)

metrics_tbl <- tibble(Zone = selected_zone, AW_mm = AW, POC_pct = POC, AAD = AAD)
print(metrics_tbl)

## PIT (TRAIN) by HMM state — single faceted figure
pit_train_zone <- df_zone %>%
  filter(is_train) %>%
  left_join(marginals_tbl, by = c("codigo_estacion","States")) %>%
  mutate(U_PIT = purrr::map2_dbl(precipitation, model, ~ F_hurdle(.x, .y)))

p_pit <- ggplot(pit_train_zone, aes(x = pmin(pmax(U_PIT, 0), 1))) +
  geom_histogram(bins = 20, color = "grey25", fill = "steelblue") +
  facet_wrap(~ States, ncol = 3) +
  labs(title = sprintf("Conditional PIT by HMM state (TRAIN) — Zone %s", selected_zone),
       x = expression(U[PIT]), y = "Count") +
  theme_bw()
print(p_pit)

## ============================== 6) Multivariate Verification (ES / VS / CRPS_areal / Events) =====================

## --- Helpers to construct fixed station sets (constant d per state in TEST) ---
get_fixed_station_set <- function(state_id, min_days_frac = 0.95){
  vm <- vine_models[[as.character(state_id)]]
  if (is.null(vm) || is.null(vm$vine) || is.null(vm$stations)) return(character(0))
  S_vine <- as.character(vm$stations)
  
  S_marg <- marginals_tbl %>%
    filter(States == state_id) %>%
    pull(codigo_estacion) %>%
    as.character() %>% unique()
  
  test_dates_s <- df_zone %>%
    filter(Year == test_year, month(dates) %in% months_MJJA, States == state_id) %>%
    distinct(dates) %>% pull()
  if (length(test_dates_s) == 0) return(character(0))
  
  keep_tbl <- df_zone %>%
    filter(Year == test_year, month(dates) %in% months_MJJA, States == state_id) %>%
    group_by(codigo_estacion) %>%
    summarise(n_days = sum(is.finite(precipitation)), .groups = "drop") %>%
    mutate(frac = n_days / length(test_dates_s)) %>%
    filter(frac >= min_days_frac)
  S_keep <- as.character(keep_tbl$codigo_estacion)
  
  S <- S_vine[S_vine %in% intersect(S_marg, S_keep)]
  S
}

# Align marginal models with a given station ordering (state-specific)
get_models_for_state <- function(state_id, station_order){
  out <- marginals_tbl %>%
    filter(States == state_id, codigo_estacion %in% station_order) %>%
    mutate(codigo_estacion = as.character(codigo_estacion)) %>%
    arrange(match(codigo_estacion, station_order))
  if (nrow(out) != length(station_order)) return(NULL)
  out
}

# Extract observed vector y(date, state, station_order)
get_y_vec <- function(date, state_id, station_order){
  day <- df_zone %>%
    filter(dates == date, States == state_id, codigo_estacion %in% station_order) %>%
    mutate(codigo_estacion = as.character(codigo_estacion)) %>%
    arrange(match(codigo_estacion, station_order))
  if (nrow(day) != length(station_order)) return(NULL)
  y <- as.numeric(day$precipitation)
  if (any(!is.finite(y))) return(NULL)
  y
}

## --- Ensemble simulators (uniform clamping for stability) ---
simulate_ensemble_vine <- function(state_id, station_order, M = 500){
  vm <- vine_models[[as.character(state_id)]]
  if (is.null(vm) || is.null(vm$vine)) return(NULL)
  U_full    <- VineCopula::RVineSim(M, vm$vine)    # [M × d_vine]
  vine_names <- as.character(vm$stations)
  idx       <- match(station_order, vine_names)
  if (any(is.na(idx))) return(NULL)
  U <- U_full[, idx, drop = FALSE]                 # [M × d]
  U <- apply(U, 2, clamp01)
  mods <- get_models_for_state(state_id, station_order)
  if (is.null(mods)) return(NULL)
  out <- matrix(NA_real_, nrow = M, ncol = length(station_order))
  for (j in seq_along(station_order)){
    out[, j] <- vapply(U[, j], Q_hurdle, numeric(1), model = mods$model[[j]])
  }
  colnames(out) <- station_order
  out
}

simulate_ensemble_indep <- function(state_id, station_order, M = 500){
  d <- length(station_order)
  if (d == 0) return(NULL)
  U <- matrix(stats::runif(M * d), nrow = M, ncol = d)
  U <- apply(U, 2, clamp01)
  mods <- get_models_for_state(state_id, station_order)
  if (is.null(mods)) return(NULL)
  out <- matrix(NA_real_, nrow = M, ncol = d)
  for (j in seq_len(d)){
    out[, j] <- vapply(U[, j], Q_hurdle, numeric(1), model = mods$model[[j]])
  }
  colnames(out) <- station_order
  out
}

## --- Robust ES/VS evaluator (auto-detect orientation: M×d vs d×M) ---
es_vs_from_ens <- function(y, ens, alpha_vs = 0.5){
  y   <- as.numeric(y)
  ens <- as.matrix(ens)
  ok  <- stats::complete.cases(as.data.frame(ens))
  ens <- ens[ok, , drop = FALSE]
  if (nrow(ens) < 2L) return(c(ES = NA_real_, VS = NA_real_))
  d <- length(y); if (d < 2L) return(c(ES = NA_real_, VS = NA_real_))
  
  candidates <- list()
  if (ncol(ens) == d) candidates <- c(candidates, list(ens))     # M×d
  if (nrow(ens) == d) candidates <- c(candidates, list(t(ens)))  # d×M -> M×d
  if (length(candidates) == 0) return(c(ES = NA_real_, VS = NA_real_))
  
  W <- matrix(1, d, d); diag(W) <- 0
  for (dat in candidates){
    ES <- tryCatch(as.numeric(scoringRules::es_sample(y = y, dat = dat)),
                   error = function(e) NA_real_)
    VS <- tryCatch(as.numeric(scoringRules::vs_sample(y = y, dat = dat, w = W, p = alpha_vs)),
                   error = function(e) NA_real_)
    if (is.finite(ES) && is.finite(VS)) return(c(ES = ES, VS = VS))
  }
  c(ES = NA_real_, VS = NA_real_)
}

## --- Day-level multivariate scores and event metrics ---
score_one_day <- function(date, state_id, station_order, M = 500, alpha_vs = 0.5,
                          k_wet_frac = 0.2, thr_joint = 10){
  
  y <- get_y_vec(date, state_id, station_order)
  if (is.null(y)) return(NULL)
  d <- length(y); if (d < 2) return(NULL)
  
  ens_v <- simulate_ensemble_vine(state_id, station_order, M = M)
  ens_i <- simulate_ensemble_indep(state_id, station_order, M = M)
  if (is.null(ens_v) || is.null(ens_i)) return(NULL)
  
  # Keep rows valid in both ensembles
  cc_v <- stats::complete.cases(as.data.frame(ens_v))
  cc_i <- stats::complete.cases(as.data.frame(ens_i))
  cc   <- cc_v & cc_i
  ens_v <- ens_v[cc, , drop = FALSE]
  ens_i <- ens_i[cc, , drop = FALSE]
  if (nrow(ens_v) < 2L || nrow(ens_i) < 2L) return(NULL)
  
  if (!is.numeric(ens_v) || !is.numeric(ens_i)) return(NULL)
  if (!all(dim(ens_v) == dim(ens_i))) return(NULL)
  if (ncol(ens_v) != d || any(!is.finite(y))) return(NULL)
  
  mv_v <- es_vs_from_ens(y, ens_v, alpha_vs = alpha_vs)
  mv_i <- es_vs_from_ens(y, ens_i, alpha_vs = alpha_vs)
  ES_v <- mv_v["ES"]; ES_i <- mv_i["ES"]
  VS_v <- mv_v["VS"]; VS_i <- mv_i["VS"]
  
  # Areal (mean) CRPS from ensemble draws
  y_areal      <- mean(y)
  draw_areal_v <- rowMeans(ens_v)
  draw_areal_i <- rowMeans(ens_i)
  
  CRPS_areal_v <- tryCatch(
    scoringRules::crps_sample(y_areal, dat = as.numeric(draw_areal_v)),
    error = function(e) NA_real_
  )
  CRPS_areal_i <- tryCatch(
    scoringRules::crps_sample(y_areal, dat = as.numeric(draw_areal_i)),
    error = function(e) NA_real_
  )
  
  # Event 1: ">= k wet stations" (k based on fraction of d) — Brier score
  k      <- max(2L, floor(k_wet_frac * d))
  obs_k  <- as.integer(sum(y > 0) >= k)
  pr_k_v <- mean(rowSums(ens_v > 0) >= k)
  pr_k_i <- mean(rowSums(ens_i > 0) >= k)
  BS_k_v <- (pr_k_v - obs_k)^2
  BS_k_i <- (pr_k_i - obs_k)^2
  
  # Event 2: pairwise joint exceedance (>= thr_joint mm) — average Brier over first P pairs
  Pairs <- t(utils::combn(seq_len(d), 2))
  P <- min(nrow(Pairs), 10L)
  if (P > 0) {
    Pairs <- Pairs[seq_len(P), , drop = FALSE]
    obs_pair  <- apply(Pairs, 1, function(ix) as.integer(all(y[ix] > thr_joint)))
    pr_pair_v <- apply(Pairs, 1, function(ix)
      mean(ens_v[, ix[1]] > thr_joint & ens_v[, ix[2]] > thr_joint))
    pr_pair_i <- apply(Pairs, 1, function(ix)
      mean(ens_i[, ix[1]] > thr_joint & ens_i[, ix[2]] > thr_joint))
    BS_pair_v <- mean((pr_pair_v - obs_pair)^2)
    BS_pair_i <- mean((pr_pair_i - obs_pair)^2)
  } else {
    BS_pair_v <- NA_real_; BS_pair_i <- NA_real_
  }
  
  tibble(
    date = as.Date(date), state = as.integer(state_id), d = d,
    ES_v = ES_v, ES_i = ES_i, VS_v = VS_v, VS_i = VS_i,
    CRPSA_v = CRPS_areal_v, CRPSA_i = CRPS_areal_i,
    BS_k_v = BS_k_v, BS_k_i = BS_k_i,
    BS_pair_v = BS_pair_v, BS_pair_i = BS_pair_i
  )
}

# ========================= Save figures (paper + supplement) =========================
# This block assumes you already created:
#   - out_dir_main <- file.path(getwd(), "fig_main")
#   - out_dir_supp <- file.path(getwd(), "fig_supp")
#   - save_fig() helper (writes PDF+PNG; see earlier block)
# If not available, the guard below (re)defines them quickly.

if (!exists("out_dir_main")) out_dir_main <- file.path(getwd(), "fig_main")
if (!exists("out_dir_supp")) out_dir_supp <- file.path(getwd(), "fig_supp")
dir.create(out_dir_main, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_supp, showWarnings = FALSE, recursive = TRUE)

if (!exists("save_fig")) {
  save_fig <- function(plot_obj, base_name,
                       path_main = out_dir_main, path_supp = out_dir_supp,
                       width = 9, height = 6, dpi = 700,
                       to_main = TRUE, to_supp = TRUE) {
    if (!inherits(plot_obj, "ggplot")) return(invisible(NULL))
    if (to_main) {
      ggsave(file.path(path_main, paste0(base_name, ".pdf")),
             plot = plot_obj, width = width, height = height, units = "in",
             device = cairo_pdf)
      ggsave(file.path(path_main, paste0(base_name, ".png")),
             plot = plot_obj, width = width, height = height, units = "in", dpi = dpi)
    }
    if (to_supp) {
      ggsave(file.path(path_supp, paste0(base_name, ".pdf")),
             plot = plot_obj, width = width, height = height, units = "in",
             device = cairo_pdf)
      ggsave(file.path(path_supp, paste0(base_name, ".png")),
             plot = plot_obj, width = width, height = height, units = "in", dpi = dpi)
    }
  }
}

# ---------------- Paper-ready (send only to main) ----------------
if (exists("p_daily")) {
  save_fig(p_daily,
           base_name = sprintf("Z%02d_MJJA%d_daily_ribbon", selected_zone, test_year),
           to_main = TRUE, to_supp = FALSE)
}
if (exists("p_month")) {
  save_fig(p_month,
           base_name = sprintf("Z%02d_MJJA%d_monthly_boxplots", selected_zone, test_year),
           to_main = TRUE, to_supp = FALSE)
}

# ---------------- Supplement (S-section naming; supplement only) ----------------
# S1 — Conditional PIT (faceted by state)
if (exists("p_pit")) {
  save_fig(p_pit,
           base_name = sprintf("S1_PIT_Z%02d", selected_zone),  # e.g., S1_PIT_Z04.pdf/png
           width = 10, height = 7,
           to_main = FALSE, to_supp = TRUE)
}

# S3 — Multivariate gains (ES / VS / CRPS_areal)
if (exists("p_es")) {
  save_fig(p_es,
           base_name = sprintf("S3_ES_gain_Z%02d", selected_zone),
           height = 4.5, to_main = FALSE, to_supp = TRUE)
}
if (exists("p_vs")) {
  save_fig(p_vs,
           base_name = sprintf("S3_VS_gain_Z%02d", selected_zone),
           height = 4.5, to_main = FALSE, to_supp = TRUE)
}
if (exists("p_crpsA")) {
  save_fig(p_crpsA,
           base_name = sprintf("S3_CRPS_areal_gain_Z%02d", selected_zone),
           height = 4.5, to_main = FALSE, to_supp = TRUE)
}
if (exists("p_box")) {
  save_fig(p_box,
           base_name = sprintf("S3_gain_boxplots_Z%02d", selected_zone),
           height = 5.5, to_main = FALSE, to_supp = TRUE)
}

# S4 — Event Brier gains
if (exists("p_bs")) {
  save_fig(p_bs,
           base_name = sprintf("S4_Brier_gains_Z%02d", selected_zone),
           height = 5.5, to_main = FALSE, to_supp = TRUE)
}

# ---------------- Key tables (supplement CSV exports) ----------------
if (exists("gof_results")) {
  readr::write_csv(
    gof_results,
    file.path(out_dir_supp, sprintf("S2_Z%02d_gof_results_tree1_rosenblatt.csv", selected_zone))
  )
}
if (exists("scores_mv")) {
  readr::write_csv(
    scores_mv,
    file.path(out_dir_supp, sprintf("S3_Z%02d_scores_multivar_MJJA%d.csv", selected_zone, test_year))
  )
}
if (exists("summary_tbl")) {
  readr::write_csv(
    summary_tbl,
    file.path(out_dir_supp, sprintf("S3_Z%02d_summary_multivar_MJJA%d.csv", selected_zone, test_year))
  )
}
if (exists("metrics_tbl")) {
  readr::write_csv(
    metrics_tbl,
    file.path(out_dir_supp, sprintf("S0_Z%02d_medoid_interval_metrics_MJJA%d.csv", selected_zone, test_year))
  )
}


# ========================= Save figures (paper + supplement) =========================
# This block assumes you already created:
#   - out_dir_main <- file.path(getwd(), "fig_main")
#   - out_dir_supp <- file.path(getwd(), "fig_supp")
#   - save_fig() helper (writes PDF+PNG; see earlier block)
# If not available, the guard below (re)defines them quickly.

if (!exists("out_dir_main")) out_dir_main <- file.path(getwd(), "fig_main")
if (!exists("out_dir_supp")) out_dir_supp <- file.path(getwd(), "fig_supp")
dir.create(out_dir_main, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_supp, showWarnings = FALSE, recursive = TRUE)

if (!exists("save_fig")) {
  save_fig <- function(plot_obj, base_name,
                       path_main = out_dir_main, path_supp = out_dir_supp,
                       width = 9, height = 6, dpi = 700,
                       to_main = TRUE, to_supp = TRUE) {
    if (!inherits(plot_obj, "ggplot")) return(invisible(NULL))
    if (to_main) {
      ggsave(file.path(path_main, paste0(base_name, ".pdf")),
             plot = plot_obj, width = width, height = height, units = "in",
             device = cairo_pdf)
      ggsave(file.path(path_main, paste0(base_name, ".png")),
             plot = plot_obj, width = width, height = height, units = "in", dpi = dpi)
    }
    if (to_supp) {
      ggsave(file.path(path_supp, paste0(base_name, ".pdf")),
             plot = plot_obj, width = width, height = height, units = "in",
             device = cairo_pdf)
      ggsave(file.path(path_supp, paste0(base_name, ".png")),
             plot = plot_obj, width = width, height = height, units = "in", dpi = dpi)
    }
  }
}

# ---------------- Paper-ready (send only to main) ----------------
if (exists("p_daily")) {
  save_fig(p_daily,
           base_name = sprintf("Z%02d_MJJA%d_daily_ribbon", selected_zone, test_year),
           to_main = TRUE, to_supp = FALSE)
}
if (exists("p_month")) {
  save_fig(p_month,
           base_name = sprintf("Z%02d_MJJA%d_monthly_boxplots", selected_zone, test_year),
           to_main = TRUE, to_supp = FALSE)
}

# ---------------- Supplement (S-section naming; supplement only) ----------------
# S1 — Conditional PIT (faceted by state)
if (exists("p_pit")) {
  save_fig(p_pit,
           base_name = sprintf("S1_PIT_Z%02d", selected_zone),  # e.g., S1_PIT_Z04.pdf/png
           width = 10, height = 7,
           to_main = FALSE, to_supp = TRUE)
}

# S3 — Multivariate gains (ES / VS / CRPS_areal)
if (exists("p_es")) {
  save_fig(p_es,
           base_name = sprintf("S3_ES_gain_Z%02d", selected_zone),
           height = 4.5, to_main = FALSE, to_supp = TRUE)
}
if (exists("p_vs")) {
  save_fig(p_vs,
           base_name = sprintf("S3_VS_gain_Z%02d", selected_zone),
           height = 4.5, to_main = FALSE, to_supp = TRUE)
}
if (exists("p_crpsA")) {
  save_fig(p_crpsA,
           base_name = sprintf("S3_CRPS_areal_gain_Z%02d", selected_zone),
           height = 4.5, to_main = FALSE, to_supp = TRUE)
}
if (exists("p_box")) {
  save_fig(p_box,
           base_name = sprintf("S3_gain_boxplots_Z%02d", selected_zone),
           height = 5.5, to_main = FALSE, to_supp = TRUE)
}

# S4 — Event Brier gains
if (exists("p_bs")) {
  save_fig(p_bs,
           base_name = sprintf("S4_Brier_gains_Z%02d", selected_zone),
           height = 5.5, to_main = FALSE, to_supp = TRUE)
}

# ---------------- Key tables (supplement CSV exports) ----------------
if (exists("gof_results")) {
  readr::write_csv(
    gof_results,
    file.path(out_dir_supp, sprintf("S2_Z%02d_gof_results_tree1_rosenblatt.csv", selected_zone))
  )
}
if (exists("scores_mv")) {
  readr::write_csv(
    scores_mv,
    file.path(out_dir_supp, sprintf("S3_Z%02d_scores_multivar_MJJA%d.csv", selected_zone, test_year))
  )
}
if (exists("summary_tbl")) {
  readr::write_csv(
    summary_tbl,
    file.path(out_dir_supp, sprintf("S3_Z%02d_summary_multivar_MJJA%d.csv", selected_zone, test_year))
  )
}
if (exists("metrics_tbl")) {
  readr::write_csv(
    metrics_tbl,
    file.path(out_dir_supp, sprintf("S0_Z%02d_medoid_interval_metrics_MJJA%d.csv", selected_zone, test_year))
  )
}


###
#######

# ======================= Robust plots for ES/VS/CRPSA/Brier gains =======================
# Prereqs: scores_mv, selected_zone, (optional) save_fig(), out_dir_supp/out_dir_main

# 1) Ensure score columns are numeric (tibbles can keep list-cols)
num_cols <- c("ES_v","ES_i","VS_v","VS_i","CRPSA_v","CRPSA_i",
              "BS_k_v","BS_k_i","BS_pair_v","BS_pair_i")
for (cn in num_cols) {
  if (is.list(scores_mv[[cn]])) scores_mv[[cn]] <- unlist(scores_mv[[cn]])
  scores_mv[[cn]] <- suppressWarnings(as.numeric(scores_mv[[cn]]))
}
# Coerce dates
if (!inherits(scores_mv$date, "Date")) {
  scores_mv$date <- as.Date(scores_mv$date)
}

# 2) Build gains with finite filtering and safe fallbacks
library(dplyr); library(tidyr); library(ggplot2)

gains <- scores_mv %>%
  mutate(
    dES      = ES_i     - ES_v,
    dVS      = VS_i     - VS_v,
    dCRPSA   = CRPSA_i  - CRPSA_v,
    dBS_k    = BS_k_i   - BS_k_v,
    dBS_pair = BS_pair_i- BS_pair_v
  ) %>%
  arrange(date)

n_finite_ES    <- sum(is.finite(gains$dES))
n_finite_VS    <- sum(is.finite(gains$dVS))
n_finite_CRPSA <- sum(is.finite(gains$dCRPSA))
n_finite_BS_k  <- sum(is.finite(gains$dBS_k))
n_finite_BS_p  <- sum(is.finite(gains$dBS_pair))

# Helper: fallback plot text if no finite values
fallback_plot <- function(title_txt, ylab_txt){
  ggplot() +
    annotate("text", x = 0, y = 0, label = "No finite values to plot\n(check ES/VS computation)",
             size = 4) +
    labs(title = title_txt, y = ylab_txt, x = NULL) +
    theme_void()
}

# 3) Time-series plots (with fallback if empty)

# ES gain
p_es <- if (n_finite_ES > 0) {
  gains %>% filter(is.finite(dES)) %>%
    ggplot(aes(x = date, y = dES)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line() + geom_point(size = 0.7) +
    labs(title = sprintf("Daily ES gain: independence − vine (Zone %s)", selected_zone),
         y = "\u0394ES (indep − vine)", x = NULL) +
    theme_bw()
} else {
  fallback_plot(sprintf("Daily ES gain: independence − vine (Zone %s)", selected_zone),
                "\u0394ES (indep − vine)")
}

# VS gain
p_vs <- if (n_finite_VS > 0) {
  gains %>% filter(is.finite(dVS)) %>%
    ggplot(aes(x = date, y = dVS)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line() + geom_point(size = 0.7) +
    labs(title = sprintf("Daily VS gain: independence − vine (Zone %s)", selected_zone),
         y = "\u0394VS (indep − vine)", x = NULL) +
    theme_bw()
} else {
  fallback_plot(sprintf("Daily VS gain: independence − vine (Zone %s)", selected_zone),
                "\u0394VS (indep − vine)")
}

# Areal CRPS gain
p_crpsA <- if (n_finite_CRPSA > 0) {
  gains %>% filter(is.finite(dCRPSA)) %>%
    ggplot(aes(x = date, y = dCRPSA)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line() + geom_point(size = 0.7) +
    labs(title = sprintf("Daily areal-CRPS gain: independence − vine (Zone %s)", selected_zone),
         y = "\u0394CRPS_areal (indep − vine)", x = NULL) +
    theme_bw()
} else {
  fallback_plot(sprintf("Daily areal-CRPS gain: independence − vine (Zone %s)", selected_zone),
                "\u0394CRPS_areal (indep − vine)")
}

# 4) Boxplots of gains across MJJA (only if any finite values exist)
g_box_df <- tibble::tibble(
  `ΔES`         = gains$dES,
  `ΔVS`         = gains$dVS,
  `ΔCRPS_areal` = gains$dCRPSA
) %>%
  tidyr::pivot_longer(everything(), names_to = "metric", values_to = "gain") %>%
  filter(is.finite(gain))

p_box <- if (nrow(g_box_df) > 0) {
  ggplot(g_box_df, aes(x = metric, y = gain)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_boxplot(fill = "grey85") +
    labs(title = sprintf("Score gains (independence − vine), Zone %s", selected_zone),
         x = NULL, y = "gain (>0 favours vine)") +
    theme_bw()
} else {
  fallback_plot(sprintf("Score gains (independence − vine), Zone %s", selected_zone),
                "gain (>0 favours vine)")
}

# 5) Event Brier gains (>=k wet; pair >= threshold)
brier_long <- gains %>%
  transmute(date,
            `ΔBS (>=k wet)`   = dBS_k,
            `ΔBS (pair>=thr)` = dBS_pair) %>%
  pivot_longer(-date, names_to = "event", values_to = "gain") %>%
  filter(is.finite(gain))

p_bs <- if (nrow(brier_long) > 0) {
  ggplot(brier_long, aes(x = date, y = gain)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line() + geom_point(size = 0.7) +
    facet_wrap(~event, ncol = 1, scales = "free_y") +
    labs(title = sprintf("Brier gains: independence − vine (Zone %s)", selected_zone),
         x = NULL, y = "ΔBrier (>0 favours vine)") +
    theme_bw()
} else {
  fallback_plot(sprintf("Brier gains: independence − vine (Zone %s)", selected_zone),
                "ΔBrier (>0 favours vine)")
}

# 6) Print to screen
print(p_es); print(p_vs); print(p_crpsA); print(p_box); print(p_bs)

# 7) Save (uses save_fig() if present; otherwise falls back to ggsave)
save_or_ggsave <- function(p, base_name, width = 8, height = 4.5) {
  if (exists("save_fig")) {
    save_fig(p, base_name = base_name, height = height, width = width)
  } else {
    dir.create("fig_supp", showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(file.path("fig_supp", paste0(base_name, ".pdf")),
                    p, width = width, height = height, device = cairo_pdf)
    ggplot2::ggsave(file.path("fig_supp", paste0(base_name, ".png")),
                    p, width = width, height = height, dpi = 700)
  }
}

save_or_ggsave(p_es,     sprintf("Z%02d_daily_ES_gain",       selected_zone))
save_or_ggsave(p_vs,     sprintf("Z%02d_daily_VS_gain",       selected_zone))
save_or_ggsave(p_crpsA,  sprintf("Z%02d_daily_CRPSA_gain",    selected_zone))
save_or_ggsave(p_box,    sprintf("Z%02d_gain_boxplots",       selected_zone), height = 5.5)
save_or_ggsave(p_bs,     sprintf("Z%02d_Brier_gains",         selected_zone), height = 5.5)
# ======================================================================================


####
#########

# ========= S1.2 — KS / CvM p-values for (multivariate) PIT via Rosenblatt =========
# Requirements: vine_models, calc_U_zone_state(), df_zone, marginals_tbl, selected_zone
# Packages: dplyr, tibble, purrr, VineCopula, goftest

build_S1_pvals_table <- function(max_components = 5) {
  states <- sort(unique(df_zone$States))
  out <- purrr::map_dfr(states, function(st) {
    vm <- vine_models[[as.character(st)]]
    if (is.null(vm) || is.null(vm$vine)) {
      return(tibble::tibble(
        Zone = selected_zone, State = st, Comp = NA_integer_,
        KS_p = NA_real_, CvM_p = NA_real_, n_days = NA_integer_,
        Notes = "No vine fitted for this state"
      ))
    }
    Uz <- calc_U_zone_state(df_zone, st, marginals_tbl)
    if (is.null(Uz) || nrow(Uz$U) == 0) {
      return(tibble::tibble(
        Zone = selected_zone, State = st, Comp = NA_integer_,
        KS_p = NA_real_, CvM_p = NA_real_, n_days = 0L,
        Notes = "No training rows for this state in zone"
      ))
    }
    # Rosenblatt transform under the fitted vine (multivariate PIT)
    Rz <- tryCatch(VineCopula::RVinePIT(Uz$U, vm$vine), error = function(e) NULL)
    if (is.null(Rz)) {
      return(tibble::tibble(
        Zone = selected_zone, State = st, Comp = NA_integer_,
        KS_p = NA_real_, CvM_p = NA_real_, n_days = nrow(Uz$U),
        Notes = "RVinePIT failed"
      ))
    }
    ktest <- min(max_components, ncol(Rz))
    tibble::tibble(
      Zone   = selected_zone,
      State  = st,
      Comp   = seq_len(ktest),
      KS_p   = purrr::map_dbl(seq_len(ktest), ~ ks.test(Rz[, .x], "punif")$p.value),
      CvM_p  = purrr::map_dbl(seq_len(ktest), ~ goftest::cvm.test(Rz[, .x], "punif")$p.value),
      n_days = nrow(Rz),
      Notes  = sprintf("Tree trunc=%s; sel=%s", as.character(vine_truncate), vine_selcrit)
    )
  })
  out
}

S1_pvals_tbl <- build_S1_pvals_table(max_components = 5)

# Save to CSV (to import in LaTeX with siunitx/csvsimple or similar)
readr::write_csv(
  S1_pvals_tbl,
  file.path(out_dir_supp, sprintf("S1_Z%02d_PIT_KS_CvM_pvalues.csv", selected_zone))
)

# (Optional) print LaTeX rows for manual pasting
apply(S1_pvals_tbl, 1, function(r) {
  cat(sprintf("%s & %s & %s & %s & %s & %s & %s \\\\\n",
              r[["Zone"]], r[["State"]], r[["Comp"]],
              format(as.numeric(r[["KS_p"]]), digits = 3, scientific = TRUE),
              format(as.numeric(r[["CvM_p"]]), digits = 3, scientific = TRUE),
              r[["n_days"]], r[["Notes"]]))
})


# ===================================================================
# S3.1 — Tail stability (mean-excess) and GPD parameter table
# Requires: df_zone, marginals_tbl, medoid_station, selected_zone, out_dir_supp, save_fig()
# ===================================================================

library(dplyr); library(tidyr); library(ggplot2); library(purrr)

# Mean-excess curve: E[X - u | X > u]
mean_excess_curve <- function(xpos, grid = NULL){
  xpos <- sort(xpos[xpos > 0])
  if (length(xpos) < 30) return(tibble(u = numeric(0), e = numeric(0)))
  if (is.null(grid)) {
    # Robust threshold grid (avoid extreme upper tail instability)
    probs <- seq(0.60, 0.98, by = 0.02)
    grid  <- as.numeric(quantile(xpos, probs = probs, names = FALSE, type = 7))
    grid  <- unique(pmin(grid, max(xpos) - .Machine$double.eps))
  }
  e_u <- map_dbl(grid, function(u){
    exc <- xpos[xpos > u] - u
    if (length(exc) < 5) return(NA_real_)
    mean(exc)
  })
  tibble(u = grid, e = e_u)
}

states_z <- sort(unique(df_zone$States))

# --- (A) Mean-excess figure per state (medoid) ---
me_df <- map_dfr(states_z, function(st){
  x <- df_zone %>%
    filter(is_train, States == st, codigo_estacion == medoid_station, precipitation > 0) %>%
    pull(precipitation)
  if (length(x) == 0) return(tibble(State = st, u = numeric(0), e = numeric(0)))
  me  <- mean_excess_curve(x)
  if (nrow(me) == 0) return(tibble(State = st, u = numeric(0), e = numeric(0)))
  me$State <- st
  me
})

thr_tbl <- marginals_tbl %>%
  filter(codigo_estacion == medoid_station, States %in% states_z) %>%
  transmute(State = States, thr = purrr::map_dbl(model, ~ .x$thr),
            n_exc = purrr::map_int(model, ~ sum(.x$x_sorted > .x$thr)))

p_me <- ggplot(me_df, aes(x = u, y = e)) +
  geom_line() +
  geom_vline(data = thr_tbl, aes(xintercept = thr), linetype = 2) +
  facet_wrap(~ State, scales = "free", ncol = 3) +
  labs(title = sprintf("Tail-stability (mean excess) — Zone %s (medoid)", selected_zone),
       x = "Threshold u (mm)", y = expression(E[X - u ~"|"~ X > u])) +
  theme_bw()

# Save with the exact name referenced by your LaTeX
save_fig(p_me,
         base_name = sprintf("S3_mean_excess_Z%02d", selected_zone),
         width = 10, height = 7, to_main = FALSE, to_supp = TRUE)

# --- (B) GPD parameter table (medoid) per state with bootstrap CIs ---
gpd_boot_ci <- function(xpos, thr, B = 400, conf = 0.95){
  xpos <- xpos[xpos > 0]
  exc  <- xpos[xpos > thr] - thr
  if (length(exc) < 25) return(c(Scale_L = NA, Scale_U = NA, Shape_L = NA, Shape_U = NA))
  set.seed(123)
  sc <- sh <- numeric(B)
  for (b in seq_len(B)){
    exc_b <- sample(exc, replace = TRUE)
    x_b   <- c(xpos[xpos <= thr], thr + exc_b)
    fit_b <- tryCatch(evd::fpot(x_b, threshold = thr, model = "gpd", std.err = FALSE),
                      error = function(e) NULL)
    if (!is.null(fit_b)) {
      sc[b] <- unname(fit_b$estimate["scale"])
      sh[b] <- unname(fit_b$estimate["shape"])
    } else {
      sc[b] <- NA; sh[b] <- NA
    }
  }
  sc <- sc[is.finite(sc)]; sh <- sh[is.finite(sh)]
  if (length(sc) < 50 || length(sh) < 50)
    return(c(Scale_L = NA, Scale_U = NA, Shape_L = NA, Shape_U = NA))
  a <- (1 - conf)/2; b <- 1 - a
  c(Scale_L = quantile(sc, a), Scale_U = quantile(sc, b),
    Shape_L = quantile(sh, a), Shape_U = quantile(sh, b))
}

gpd_tbl <- map_dfr(states_z, function(st){
  md <- marginals_tbl %>% filter(codigo_estacion == medoid_station, States == st)
  if (nrow(md) == 0) {
    return(tibble(State = st, Threshold_mm = NA_real_, Scale = NA_real_,
                  Scale_L = NA_real_, Scale_U = NA_real_,
                  Shape = NA_real_, Shape_L = NA_real_, Shape_U = NA_real_,
                  n_exc = NA_integer_))
  }
  m    <- md$model[[1]]
  thr  <- as.numeric(m$thr)
  xpos <- m$x_sorted
  nexc <- sum(xpos > thr)
  if (is.null(m$gpd)) {
    tibble(State = st, Threshold_mm = thr, Scale = NA_real_,
           Scale_L = NA_real_, Scale_U = NA_real_,
           Shape = NA_real_, Shape_L = NA_real_, Shape_U = NA_real_,
           n_exc = nexc)
  } else {
    est   <- m$gpd$estimate
    ci    <- gpd_boot_ci(xpos, thr, B = 400, conf = 0.95)
    tibble(State = st,
           Threshold_mm = thr,
           Scale  = unname(est["scale"]),
           Scale_L = ci["Scale_L"], Scale_U = ci["Scale_U"],
           Shape  = unname(est["shape"]),
           Shape_L = ci["Shape_L"], Shape_U = ci["Shape_U"],
           n_exc = nexc)
  }
})

readr::write_csv(
  gpd_tbl,
  file.path(out_dir_supp, sprintf("S3_Z%02d_gpd_params_medoid.csv", selected_zone))
)

# ==================================================================
# S4 — Vine-copula selection & goodness-of-fit
# ==================================================================

# ---- Helpers -----------------------------------------------------
fam_is_arch <- function(f) {
  f %in% c(
    3,4,5,6,7,8,9,10,         # Clayton, Gumbel, Frank, Joe, BB1, BB6, BB7, BB8
    13,14,16,17,18,19,20,     # rotations
    23,24,26,27,28,29,30,
    33,34,36,37,38,39,40
  )
}
detect_trunc <- function(rv) {
  if (exists("vine_truncate")) return(vine_truncate)
  d <- ncol(rv$Matrix)
  max_tree <- 0
  for (t in 1:(d-1)) {
    col <- d - t + 1
    fam_col <- rv$family[(d-1):(d-t), col]
    if (any(fam_col != 0, na.rm = TRUE)) max_tree <- t
  }
  max_tree
}
num <- function(x) ifelse(is.finite(x), x, NA_real_)
txt <- function(x) {
  x <- as.character(x)
  ifelse(is.na(x) | x=="NA", "", x)
}

# ---- S4.1: Structure summary and family mix ----------------------
summarize_vine_structure <- function(selected_zone) {
  states <- sort(unique(df_zone$States))
  purrr::map_dfr(states, function(st) {
    vm <- vine_models[[as.character(st)]]
    if (is.null(vm) || is.null(vm$vine)) {
      return(tibble::tibble(
        Zone = selected_zone, State = st, Type = "(none)", Trunc = "",
        pct_Gaussian = "", pct_Arch = "", Notes = "no vine"
      ))
    }
    rv  <- vm$vine
    fam <- as.vector(rv$family)
    fam <- fam[fam != 0]                 # exclude independences
    nE  <- length(fam)
    pctG <- if (nE > 0) mean(fam == 1)*100 else NA_real_  # Gaussian share
    pctA <- if (nE > 0) mean(fam_is_arch(fam))*100 else NA_real_
    Uz <- calc_U_zone_state(df_zone, st, marginals_tbl)
    n_rows <- if (is.null(Uz)) NA_integer_ else nrow(Uz$U)
    n_cols <- if (is.null(Uz)) NA_integer_ else ncol(Uz$U)
    tibble::tibble(
      Zone = selected_zone,
      State = st,
      Type = "R-vine",
      Trunc = detect_trunc(rv),
      pct_Gaussian = round(pctG, 1),
      pct_Arch     = round(pctA, 1),
      Notes = sprintf("n=%s, d=%s, edges=%d",
                      txt(n_rows), txt(n_cols), nE)
    )
  }) %>%
    dplyr::mutate(
      dplyr::across(c(pct_Gaussian, pct_Arch), ~ ifelse(is.finite(.), ., NA_real_)),
      dplyr::across(everything(), ~ ifelse(is.na(.), "", .))
    )
}

S4_struct <- summarize_vine_structure(selected_zone)

readr::write_csv(
  S4_struct,
  file.path(out_dir_supp, sprintf("S4_Z%02d_vine_structure.csv", selected_zone))
)

# Optional: LaTeX rows
for (i in seq_len(nrow(S4_struct))) {
  r <- S4_struct[i, ]
  cat(sprintf("%s & %s & %s & %s & %s & %s & %s \\\\\n",
              r$Zone, r$State, r$Type, r$Trunc,
              r$pct_Gaussian, r$pct_Arch, r$Notes))
}

# ---- S4.2: White tests on Tree 1 --------------------------------
extract_T1_pairs <- function(rvine_obj){
  if (is.null(rvine_obj) || is.null(rvine_obj$Matrix)) return(NULL)
  d <- ncol(rvine_obj$Matrix)
  tibble::tibble(
    i      = rvine_obj$Matrix[d, (d-1):1],
    j      = diag(rvine_obj$Matrix)[(d-1):1],
    family = rvine_obj$family[(d-1):1, d],
    par    = rvine_obj$par[(d-1):1, d],
    par2   = rvine_obj$par2[(d-1):1, d]
  )
}

white_tbl <- purrr::imap_dfr(vine_models, function(vm, st_chr){
  st <- as.integer(st_chr)
  Uz <- calc_U_zone_state(df_zone, st, marginals_tbl)
  if (is.null(vm) || is.null(vm$vine) || is.null(Uz) || nrow(Uz$U) == 0) return(NULL)
  pairs <- extract_T1_pairs(vm$vine)
  if (is.null(pairs) || nrow(pairs) == 0) return(NULL)
  purrr::pmap_dfr(pairs, function(i, j, family, par, par2){
    ui <- Uz$U[, as.integer(i)]
    uj <- Uz$U[, as.integer(j)]
    wt <- tryCatch(
      VineCopula::BiCopGofTest(ui, uj, family = family, par = par, par2 = par2,
                               method = "White", B = 200),
      error = function(e) NULL
    )
    pval <- if (is.null(wt)) NA_real_ else wt$p.value
    tibble::tibble(
      Zone   = selected_zone,
      State  = st,
      Pair   = sprintf("(%d,%d)", as.integer(i), as.integer(j)),
      Family = tryCatch(VineCopula::BiCopName(family = family, short = TRUE),
                        error = function(e) as.character(family)),
      p_value = pval,
      n      = nrow(Uz$U)
    )
  })
})

if (is.null(white_tbl) || nrow(white_tbl) == 0) {
  white_tbl <- tibble::tibble(
    Zone=selected_zone, State=integer(), Pair=character(),
    Family=character(), p_value=double(), n=integer()
  )
}

# Export detailed and per-state-minimum CSVs
white_export <- white_tbl %>%
  dplyr::mutate(
    p_str = dplyr::if_else(is.finite(p_value),
                           formatC(p_value, format = "e", digits = 2),
                           "")
  ) %>%
  dplyr::select(Zone, State, Pair, Family, p_value = p_str, n)

readr::write_csv(
  white_export,
  file.path(out_dir_supp, sprintf("S4_Z%02d_white_T1.csv", selected_zone))
)

white_min <- white_tbl %>%
  dplyr::group_by(Zone, State) %>%
  dplyr::summarise(min_p = suppressWarnings(min(p_value, na.rm = TRUE)),
                   n = dplyr::first(n), .groups = "drop") %>%
  dplyr::mutate(min_p = ifelse(is.finite(min_p),
                               formatC(min_p, format = "e", digits = 2), ""))
readr::write_csv(
  white_min,
  file.path(out_dir_supp, sprintf("S4_Z%02d_white_T1_min.csv", selected_zone))
)

# Optional: LaTeX rows (detailed)
for (i in seq_len(nrow(white_export))) {
  r <- white_export[i, ]
  cat(sprintf("%s & %s & %s & %s & %s & %s \\\\\n",
              r$Zone, r$State, r$Pair, r$Family, r$p_value, r$n))
}

# -----------------------------------------------------------------
# S5 — Multivariate scoring (ES/VS)
# -----------------------------------------------------------------

# ---- S5.1 Station-set coverage per state (TEST period) ----
S5_cov <- scores_mv %>%
  dplyr::mutate(zone = selected_zone) %>%
  dplyr::group_by(zone, state) %>%
  dplyr::summarise(d = max(d, na.rm = TRUE),
                   days_scored = dplyr::n(),
                   .groups = "drop") %>%
  dplyr::arrange(state)

readr::write_csv(
  S5_cov,
  file.path(out_dir_supp, sprintf("S5_Z%02d_station_set_coverage.csv", selected_zone))
)

# Optional: LaTeX rows
apply(S5_cov, 1, function(r){
  cat(sprintf("%d & %d & %d & %d \\\\\n",
              as.integer(r[["zone"]]),
              as.integer(r[["state"]]),
              as.integer(r[["d"]]),
              as.integer(r[["days_scored"]])))
})

# ---- S5.2 Mean ES/VS and Diebold–Mariano p-values (vine vs independence) ----
# Lower scores are better; H1: vine < indep (one-sided DM, "less")
dm_safe <- function(x, y, h = 1L){
  p <- NA_real_
  out <- tryCatch(forecast::dm.test(x, y, alternative = "less", h = h),
                  error = function(e) NULL)
  if (!is.null(out) && is.finite(out$p.value)) p <- out$p.value
  p
}

S5_dm <- scores_mv %>%
  dplyr::mutate(zone = selected_zone) %>%
  dplyr::group_by(zone) %>%
  dplyr::summarise(
    ES_vine   = mean(ES_v, na.rm = TRUE),
    ES_indep  = mean(ES_i, na.rm = TRUE),
    pDM_ES    = dm_safe(ES_v, ES_i, h = 1L),
    VS_vine   = mean(VS_v, na.rm = TRUE),
    VS_indep  = mean(VS_i, na.rm = TRUE),
    pDM_VS    = dm_safe(VS_v, VS_i, h = 1L),
    .groups = "drop"
  )

readr::write_csv(
  S5_dm,
  file.path(out_dir_supp, sprintf("S5_Z%02d_ESVS_DM.csv", selected_zone))
)

# Optional: LaTeX rows
apply(S5_dm, 1, function(r){
  cat(sprintf("%d & %.3f & %.3f & %s & %.3f & %.3f & %s \\\\\n",
              as.integer(r[["zone"]]),
              as.numeric(r[["ES_vine"]]),
              as.numeric(r[["ES_indep"]]),
              ifelse(is.finite(as.numeric(r[["pDM_ES"]])) ,
                     format(as.numeric(r[["pDM_ES"]]), digits = 3, scientific = TRUE), "---"),
              as.numeric(r[["VS_vine"]]),
              as.numeric(r[["VS_indep"]]),
              ifelse(is.finite(as.numeric(r[["pDM_VS"]])) ,
                     format(as.numeric(r[["pDM_VS"]]), digits = 3, scientific = TRUE), "---")))
})

# =================================================================
# S6.1 — Areal CRPS (vine vs independence) + DM("less")
# =================================================================

# Ensure numeric vectors
num_cols <- c("CRPSA_v","CRPSA_i","BS_k_v","BS_k_i","BS_pair_v","BS_pair_i")
for (cn in num_cols) {
  if (is.list(scores_mv[[cn]])) scores_mv[[cn]] <- unlist(scores_mv[[cn]])
  scores_mv[[cn]] <- suppressWarnings(as.numeric(scores_mv[[cn]]))
}
if (!inherits(scores_mv$date, "Date")) scores_mv$date <- as.Date(scores_mv$date)

safe_mean <- function(x) { x <- x[is.finite(x)]; if (length(x)==0) return(NA_real_); mean(x) }
fmt_num   <- function(x) if (is.finite(x)) sprintf("%.3f", x) else "---"
fmt_p     <- function(p) if (is.finite(p)) format(p, digits=3, scientific=TRUE) else "---"

dm_less <- function(x, y, h=1L) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 10L) return(NA_real_)
  out <- tryCatch(forecast::dm.test(x, y, alternative="less", h=h),
                  error = function(e) NULL)
  if (is.null(out) || !is.finite(out$p.value)) NA_real_ else out$p.value
}

S6_crps <- scores_mv %>%
  dplyr::mutate(zone = selected_zone) %>%
  dplyr::group_by(zone) %>%
  dplyr::summarise(
    CRPS_vine  = safe_mean(CRPSA_v),
    CRPS_indep = safe_mean(CRPSA_i),
    pDM_CRPS   = dm_less(CRPSA_v, CRPSA_i, h=1L),
    .groups="drop"
  )

# CSV export
readr::write_csv(S6_crps, file.path(out_dir_supp, sprintf("S6_Z%02d_areal_CRPS_DM.csv", selected_zone)))

# Optional: LaTeX rows
if (nrow(S6_crps)==0) {
  cat(sprintf("%d & --- & --- & --- \\\\\n", selected_zone))
} else {
  for (i in seq_len(nrow(S6_crps))) {
    cat(sprintf("%d & %s & %s & %s \\\\\n",
                as.integer(S6_crps$zone[i]),
                fmt_num(S6_crps$CRPS_vine[i]),
                fmt_num(S6_crps$CRPS_indep[i]),
                fmt_p(S6_crps$pDM_CRPS[i])))
  }
}

# Optional: save S6 figure from existing plot
if (exists("p_crpsA")) {
  if (exists("save_fig")) {
    save_fig(p_crpsA, base_name=sprintf("S6_areal_CRPS_gains_Z%02d", selected_zone),
             height=4.5, to_main=FALSE, to_supp=TRUE)
  } else {
    dir.create("fig_supp", showWarnings=FALSE, recursive=TRUE)
    ggplot2::ggsave(file.path("fig_supp", sprintf("S6_areal_CRPS_gains_Z%02d.pdf", selected_zone)),
                    p_crpsA, width=8, height=4.5, device=cairo_pdf)
  }
}

# =================================================================
# S6.2 — Brier scores for events (>=k wet) and (pair >= threshold) + DM
# =================================================================

S6_brier <- scores_mv %>%
  dplyr::mutate(zone = selected_zone) %>%
  dplyr::group_by(zone) %>%
  dplyr::summarise(
    BS_k_vine    = safe_mean(BS_k_v),
    BS_k_indep   = safe_mean(BS_k_i),
    pDM_BS_k     = dm_less(BS_k_v, BS_k_i, h=1L),
    BS_pair_vine = safe_mean(BS_pair_v),
    BS_pair_indep= safe_mean(BS_pair_i),
    pDM_BS_pair  = dm_less(BS_pair_v, BS_pair_i, h=1L),
    .groups="drop"
  )

readr::write_csv(S6_brier, file.path(out_dir_supp, sprintf("S6_Z%02d_Brier_events_DM.csv", selected_zone)))

# Optional: LaTeX rows
if (nrow(S6_brier)==0) {
  cat(sprintf("%d & --- & --- & --- & --- & --- & --- \\\\\n", selected_zone))
} else {
  for (i in seq_len(nrow(S6_brier))) {
    cat(sprintf("%d & %s & %s & %s & %s & %s & %s \\\\\n",
                as.integer(S6_brier$zone[i]),
                fmt_num(S6_brier$BS_k_vine[i]),   fmt_num(S6_brier$BS_k_indep[i]),   fmt_p(S6_brier$pDM_BS_k[i]),
                fmt_num(S6_brier$BS_pair_vine[i]),fmt_num(S6_brier$BS_pair_indep[i]),fmt_p(S6_brier$pDM_BS_pair[i])))
  }
}

# Optional: save S6 Brier figure (reuses p_bs if available)
if (exists("p_bs")) {
  if (exists("save_fig")) {
    save_fig(p_bs, base_name=sprintf("S6_Brier_events_Z%02d", selected_zone),
             height=5.0, to_main=FALSE, to_supp=TRUE)
  } else {
    dir.create("fig_supp", showWarnings=FALSE, recursive=TRUE)
    ggplot2::ggsave(file.path("fig_supp", sprintf("S6_Brier_events_Z%02d.pdf", selected_zone)),
                    p_bs, width=8, height=5.0, device=cairo_pdf)
  }
}

# =================================================================
# S8 — Out-of-sample verification
# MJJA 2022 (current) and rolling-origin skeleton 2012–2021
# Requires: datos_base_l, df_zone, marginals_tbl, vine_models,
#           simulate_ensemble_vine(), simulate_ensemble_indep(),
#           get_fixed_station_set(), get_y_vec(), clamp01, Q_hurdle, save_fig, etc.
# =================================================================

pacman::p_load(patchwork)

# ---------- S8.1 PIT (medoid) and areal-CRPS distribution (MJJA test_year) ----------
# Uses: fixed station set per state for areal-CRPS; medoid PIT from the VINE ensemble.

compute_S8_for_zone_year <- function(zone_id = selected_zone, year_eval = test_year,
                                     months_eval = 5:8, M = 500,
                                     alpha_vs = 0.5, k_wet_frac = 0.2, thr_joint = 10) {
  
  # Subset zone
  df_z <- datos_base_l %>% dplyr::filter(labels == zone_id)
  
  # Zone medoid
  medoid_station_z <- df_z %>% dplyr::filter(medoids == 1) %>%
    dplyr::pull(codigo_estacion) %>% unique()
  if (length(medoid_station_z) != 1) return(NULL)
  
  # MJJA dates for evaluation year
  dates_eval <- df_z %>%
    dplyr::filter(lubridate::year(dates) == year_eval,
                  lubridate::month(dates) %in% months_eval) %>%
    dplyr::pull(dates) %>% unique() %>% sort()
  
  # Loop days: state (medoid), M×d simulations (VINE/INDEP), areal metrics, and medoid PIT
  pit_vals <- c()
  crpsA_v  <- c()
  crpsA_i  <- c()
  
  for (d0 in dates_eval) {
    # Medoid state at day d0
    st0 <- df_z %>%
      dplyr::filter(dates == d0, codigo_estacion == medoid_station_z) %>%
      dplyr::pull(States) %>% unique()
    if (length(st0) == 0) next
    
    # Fixed station set for state st0
    station_set <- get_fixed_station_set(st0, min_days_frac = 0.95)
    if (length(station_set) < 2) next
    
    # Observations for that set
    y <- get_y_vec(d0, st0, station_set)
    if (is.null(y) || any(!is.finite(y))) next
    y_areal <- mean(y)
    
    # Ensembles
    ens_v <- simulate_ensemble_vine(st0, station_set, M = M)
    ens_i <- simulate_ensemble_indep(st0, station_set, M = M)
    if (is.null(ens_v) || is.null(ens_i)) next
    
    # Areal draws and CRPS
    draw_areal_v <- rowMeans(ens_v)
    draw_areal_i <- rowMeans(ens_i)
    crps_v <- tryCatch(scoringRules::crps_sample(y_areal, dat = as.numeric(draw_areal_v)),
                       error = function(e) NA_real_)
    crps_i <- tryCatch(scoringRules::crps_sample(y_areal, dat = as.numeric(draw_areal_i)),
                       error = function(e) NA_real_)
    crpsA_v <- c(crpsA_v, crps_v)
    crpsA_i <- c(crpsA_i, crps_i)
    
    # Medoid PIT using VINE ensemble (univariate, ECDF of medoid ensemble)
    vm <- vine_models[[as.character(st0)]]
    med_col <- if (!is.null(vm) && !is.null(vm$stations)) {
      match(as.character(medoid_station_z), vm$stations)
    } else NA_integer_
    if (is.na(med_col)) {
      # Fallback: simulate U ~ Unif for the medoid and invert via its marginal
      u_med  <- stats::runif(M)
      med_mod <- marginals_tbl %>% dplyr::filter(codigo_estacion == medoid_station_z,
                                                 States == st0) %>% dplyr::pull(model)
      if (length(med_mod) == 1) {
        x_med_draw <- vapply(u_med, Q_hurdle, numeric(1), model = med_mod[[1]])
      } else {
        x_med_draw <- rep(NA_real_, M)
      }
    } else {
      x_med_draw <- ens_v[, med_col]
    }
    # Observed medoid y
    y_med <- df_z %>% dplyr::filter(dates == d0, codigo_estacion == medoid_station_z) %>%
      dplyr::pull(precipitation) %>% as.numeric()
    if (length(y_med) == 1 && is.finite(y_med) && any(is.finite(x_med_draw))) {
      Fhat <- stats::ecdf(x_med_draw[is.finite(x_med_draw)])
      pit_vals <- c(pit_vals, Fhat(y_med))
    }
  }
  
  # Figures: PIT histogram and CRPS_areal distribution (violin/box)
  df_crps <- tibble::tibble(
    model = rep(c("vine","indep"),
                times = c(length(crpsA_v), length(crpsA_i))),
    CRPS  = c(crpsA_v, crpsA_i)
  ) %>% dplyr::filter(is.finite(CRPS))
  
  p_pit <- ggplot2::ggplot(data.frame(PIT = pit_vals), aes(x = PIT)) +
    geom_histogram(bins = 20, color = "grey25", fill = "steelblue") +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = sprintf("Zone %02d — PIT (medoid, MJJA %d)", zone_id, year_eval),
         x = "PIT", y = "Count") + theme_bw()
  
  p_crps <- ggplot2::ggplot(df_crps, aes(x = model, y = CRPS)) +
    geom_violin(fill = "grey85", color = "grey35") +
    geom_boxplot(width = 0.25, outlier.size = 0.4) +
    labs(title = sprintf("Zone %02d — Areal CRPS (MJJA %d)", zone_id, year_eval),
         x = NULL, y = "CRPS (areal)") + theme_bw()
  
  p_panel <- p_pit + p_crps + patchwork::plot_layout(ncol = 2)
  save_fig(p_panel,
           base_name = "S8_pit_crps_panels",  # -> fig_supp/S8_pit_crps_panels.pdf/png
           to_main = FALSE, to_supp = TRUE, width = 11, height = 4.8)
  
  # --------- Areal pinball loss for tau in {0.90, 0.95} ----------
  pinball_loss <- function(y, q, tau) (tau - as.numeric(y < q)) * (y - q)
  
  rows_pb <- list()
  for (d0 in dates_eval) {
    st0 <- df_z %>% dplyr::filter(dates == d0, codigo_estacion == medoid_station_z) %>%
      dplyr::pull(States) %>% unique()
    if (length(st0) == 0) next
    
    station_set <- get_fixed_station_set(st0, min_days_frac = 0.95)
    if (length(station_set) < 2) next
    
    y <- get_y_vec(d0, st0, station_set)
    if (is.null(y) || any(!is.finite(y))) next
    y_areal <- mean(y)
    
    ens_v <- simulate_ensemble_vine(st0, station_set, M = M)
    ens_i <- simulate_ensemble_indep(st0, station_set, M = M)
    if (is.null(ens_v) || is.null(ens_i)) next
    
    draw_areal_v <- rowMeans(ens_v)
    draw_areal_i <- rowMeans(ens_i)
    
    for (tau in c(0.90, 0.95)) {
      qv <- stats::quantile(draw_areal_v, tau, na.rm = TRUE, names = FALSE)
      qi <- stats::quantile(draw_areal_i, tau, na.rm = TRUE, names = FALSE)
      rows_pb[[length(rows_pb) + 1]] <- tibble::tibble(
        zone = zone_id, date = d0, tau = tau,
        PL_vine  = pinball_loss(y_areal, qv, tau),
        PL_indep = pinball_loss(y_areal, qi, tau)
      )
    }
  }
  pinball_tbl <- dplyr::bind_rows(rows_pb)
  
  # MJJA averages per zone
  pinball_zone <- pinball_tbl %>%
    dplyr::group_by(zone, tau) %>%
    dplyr::summarise(PL_vine  = mean(PL_vine,  na.rm = TRUE),
                     PL_indep = mean(PL_indep, na.rm = TRUE),
                     .groups = "drop")
  
  readr::write_csv(pinball_zone,
                   file.path(out_dir_supp, sprintf("S8_Z%02d_pinball_MJJA%d.csv", zone_id, year_eval)))
  
  # Optional: LaTeX lines for \ref{tab:S8-pinball} (rounded to 3 decimals)
  cat("% --- S8 pinball (areal) LaTeX rows ---\n")
  for (tau in c(0.90, 0.95)) {
    rr <- pinball_zone %>% dplyr::filter(tau == !!tau)
    if (nrow(rr) == 1) {
      cat(sprintf("%d & PL$_{%.2f}$ vine: %.3f & PL$_{%.2f}$ indep: %.3f \\\\\n",
                  rr$zone, tau, rr$PL_vine, tau, rr$PL_indep))
    }
  }
  
  invisible(list(pit = pit_vals, crps = df_crps, pinball_zone = pinball_zone))
}

# Run for current zone/year
S8_res <- compute_S8_for_zone_year(zone_id = selected_zone, year_eval = test_year)

# ---------- (Optional) Rolling-origin 2012–2021 skeleton ----------
# Note: for true rolling-origin, refit marginals_tbl and vine_models
# with train_end = 31-Aug-(year-1) before simulating MJJA 'year'.

# years_hind <- 2012:2021
# zones_all  <- sort(unique(datos_base_l$labels))
# res_hind   <- list()
# for (yy in years_hind) {
#   train_end_y <- as.Date(sprintf("%d-08-31", yy - 1))
#   # 1) Recompute is_train for that cut:
#   datos_base_l <- datos_base_l %>%
#     dplyr::mutate(is_train = dates <= train_end_y)
#   # 2) Refit marginals_tbl and vine_models (reuse your code from above).
#   # 3) For each zone, call compute_S8_for_zone_year(zone_id=zz, year_eval=yy).
# }

# =================================================================
# S9 — Sensitivity panels
# Requires: datos_base_l, df_zone, selected_zone, test_year, months_MJJA,
#           clamp01(), es_vs_from_ens(), save_fig(), and dplyr/ggplot2/tidyr/purrr/readr
# =================================================================

library(dplyr); library(tidyr); library(purrr); library(ggplot2)

# 1) Sensitivity grid
family_opts <- list(
  "Gaussian"        = c(1),               # Gaussian only
  "Gauss+Arch"      = c(1,3,4,5,6),       # Gaussian + Archimedean (Clayton, Gumbel, Frank, Joe)
  "Gauss+t+Arch"    = c(1,2,3,4,5,6)      # add t
)
trunc_vals  <- c(2,3,4)
u_vals      <- c(0.85, 0.90, 0.95)

M_per_day_s9 <- 200   # ensemble size for POC/diagnostics
alpha_vs     <- 0.5
k_wet_frac   <- 0.2
thr_pair     <- 10

# 2) Parameterized helpers (same design, explicit arguments)
fit_margin_state_u <- function(x, u_star = 0.90, min_excess = 25){
  x <- x[is.finite(x)]
  if(length(x) < 40L){
    p0   <- mean(x <= 0, na.rm = TRUE)
    xpos <- sort(x[x > 0]); Fe <- if(length(xpos)>0) stats::ecdf(xpos) else NULL
    return(list(p0=p0, Fe=Fe, x_sorted=xpos, thr=Inf, Fu=1, gpd=NULL))
  }
  p0   <- mean(x <= 0, na.rm = TRUE)
  xpos <- sort(x[x > 0]); Fe <- if(length(xpos)>0) stats::ecdf(xpos) else NULL
  if(is.null(Fe) || length(xpos)<10L) return(list(p0=p0, Fe=NULL, x_sorted=xpos, thr=Inf, Fu=1, gpd=NULL))
  thr <- stats::quantile(xpos, probs=u_star, na.rm=TRUE, type=7)
  Fu  <- Fe(thr); exc <- xpos[xpos>thr] - thr
  gpd <- NULL
  if(length(exc) >= min_excess){
    gpd <- tryCatch(evd::fpot(xpos, threshold=thr, model="gpd", std.err=FALSE), error=function(e) NULL)
  }
  list(p0=p0, Fe=Fe, x_sorted=xpos, thr=thr, Fu=Fu, gpd=gpd)
}

build_marginals_tbl <- function(u_star, min_excess = min_excess_gpd){
  datos_base_l %>%
    filter(is_train) %>%
    group_by(codigo_estacion, States) %>%
    summarise(model = list(fit_margin_state_u(precipitation, u_star=u_star, min_excess=min_excess)),
              .groups="drop")
}

calc_U_zone_state_by <- function(df_zone, state_id, marginals_tbl){
  df <- df_zone %>%
    filter(is_train, States==state_id) %>%
    select(dates, codigo_estacion, States, precipitation)
  if(nrow(df)==0) return(NULL)
  dfU <- df %>%
    left_join(marginals_tbl, by=c("codigo_estacion","States")) %>%
    mutate(U = purrr::map2_dbl(precipitation, model, ~ F_hurdle(.x, .y))) %>%
    mutate(codigo_estacion = as.character(codigo_estacion)) %>%
    select(dates,codigo_estacion,U) %>%
    tidyr::pivot_wider(names_from=codigo_estacion, values_from=U) %>%
    arrange(dates) %>% suppressWarnings(tidyr::drop_na())
  if(nrow(dfU)==0) return(NULL)
  list(U = as.matrix(dfU %>% select(-dates)),
       station_ids = colnames(dfU %>% select(-dates)))
}

fit_vines_by <- function(marginals_tbl, familyset, trunclevel, min_rows = min_rows_vine){
  out <- vector("list", length = 5)
  names(out) <- as.character(sort(unique(df_zone$States)))
  for(st in sort(unique(df_zone$States))){
    Uz <- calc_U_zone_state_by(df_zone, st, marginals_tbl)
    if(is.null(Uz) || nrow(Uz$U) < min_rows){ out[[as.character(st)]] <- NULL; next }
    vfit <- tryCatch(
      VineCopula::RVineStructureSelect(Uz$U, familyset = familyset,
                                       selectioncrit = vine_selcrit,
                                       trunclevel = trunclevel),
      error=function(e) NULL
    )
    out[[as.character(st)]] <- list(vine=vfit, stations = Uz$station_ids)
  }
  out
}

get_fixed_station_set_by <- function(vine_models, marginals_tbl, state_id, min_days_frac=0.95){
  vm <- vine_models[[as.character(state_id)]]
  if (is.null(vm) || is.null(vm$vine) || is.null(vm$stations)) return(character(0))
  S_vine <- as.character(vm$stations)
  S_marg <- marginals_tbl %>% filter(States==state_id) %>% pull(codigo_estacion) %>% as.character() %>% unique()
  test_dates_s <- df_zone %>%
    filter(Year==test_year, month(dates)%in%months_MJJA, States==state_id) %>%
    distinct(dates) %>% pull()
  if(length(test_dates_s)==0) return(character(0))
  keep_tbl <- df_zone %>%
    filter(Year==test_year, month(dates)%in%months_MJJA, States==state_id) %>%
    group_by(codigo_estacion) %>%
    summarise(n_days = sum(is.finite(precipitation)), .groups="drop") %>%
    mutate(frac = n_days/length(test_dates_s)) %>% filter(frac>=min_days_frac)
  S_keep <- as.character(keep_tbl$codigo_estacion)
  S_vine[S_vine %in% intersect(S_marg, S_keep)]
}

get_models_for_state_by <- function(marginals_tbl, state_id, station_order){
  out <- marginals_tbl %>%
    filter(States==state_id, codigo_estacion %in% station_order) %>%
    mutate(codigo_estacion = as.character(codigo_estacion)) %>%
    arrange(match(codigo_estacion, station_order))
  if(nrow(out)!=length(station_order)) return(NULL)
  out
}

simulate_ensemble_vine_by <- function(vine_models, marginals_tbl, state_id, station_order, M=500){
  vm <- vine_models[[as.character(state_id)]]
  if(is.null(vm) || is.null(vm$vine)) return(NULL)
  U_full <- VineCopula::RVineSim(M, vm$vine)
  vine_names <- as.character(vm$stations)
  idx <- match(station_order, vine_names)
  if(any(is.na(idx))) return(NULL)
  U <- U_full[, idx, drop=FALSE]; U <- apply(U, 2, clamp01)
  mods <- get_models_for_state_by(marginals_tbl, state_id, station_order); if(is.null(mods)) return(NULL)
  out <- matrix(NA_real_, nrow=M, ncol=length(station_order))
  for(j in seq_along(station_order)){
    out[,j] <- vapply(U[,j], Q_hurdle, numeric(1), model = mods$model[[j]])
  }
  colnames(out) <- station_order
  out
}

simulate_ensemble_indep_by <- function(marginals_tbl, state_id, station_order, M=500){
  d <- length(station_order); if(d==0) return(NULL)
  U <- matrix(runif(M*d), nrow=M, ncol=d); U <- apply(U, 2, clamp01)
  mods <- get_models_for_state_by(marginals_tbl, state_id, station_order); if(is.null(mods)) return(NULL)
  out <- matrix(NA_real_, nrow=M, ncol=d)
  for(j in seq_len(d)){
    out[,j] <- vapply(U[,j], Q_hurdle, numeric(1), model = mods$model[[j]])
  }
  colnames(out) <- station_order
  out
}

score_one_day_by <- function(date, state_id, station_order, vine_models, marginals_tbl,
                             M=400, alpha_vs=0.5, k_wet_frac=0.2, thr_joint=10){
  day <- df_zone %>%
    filter(dates==date, States==state_id, codigo_estacion %in% station_order) %>%
    mutate(codigo_estacion = as.character(codigo_estacion)) %>%
    arrange(match(codigo_estacion, station_order))
  if(nrow(day)!=length(station_order)) return(NULL)
  y <- as.numeric(day$precipitation); if(any(!is.finite(y))) return(NULL)
  ens_v <- simulate_ensemble_vine_by(vine_models, marginals_tbl, state_id, station_order, M)
  ens_i <- simulate_ensemble_indep_by(marginals_tbl, state_id, station_order, M)
  if(is.null(ens_v) || is.null(ens_i)) return(NULL)
  cc <- complete.cases(as.data.frame(ens_v)) & complete.cases(as.data.frame(ens_i))
  ens_v <- ens_v[cc,,drop=FALSE]; ens_i <- ens_i[cc,,drop=FALSE]
  if(nrow(ens_v)<2L || nrow(ens_i)<2L) return(NULL)
  
  mv_v <- es_vs_from_ens(y, ens_v, alpha_vs = alpha_vs)
  mv_i <- es_vs_from_ens(y, ens_i, alpha_vs = alpha_vs)
  ES_v <- mv_v["ES"]; ES_i <- mv_i["ES"]
  VS_v <- mv_v["VS"]; VS_i <- mv_i["VS"]
  
  y_areal <- mean(y)
  CRPS_v  <- tryCatch(scoringRules::crps_sample(y_areal, dat = rowMeans(ens_v)), error=function(e) NA_real_)
  CRPS_i  <- tryCatch(scoringRules::crps_sample(y_areal, dat = rowMeans(ens_i)), error=function(e) NA_real_)
  
  tibble(ES_gain = ES_i - ES_v, VS_gain = VS_i - VS_v, CRPS_gain = CRPS_i - CRPS_v)
}

# 3) Evaluate sensitivity across the full grid
sens_results <- list()

for(flab in names(family_opts)){
  famset <- family_opts[[flab]]
  for(trunclevel in trunc_vals){
    for(u_star_cur in u_vals){
      
      # (i) Margins with threshold u*
      marg_tbl <- build_marginals_tbl(u_star = u_star_cur, min_excess = min_excess_gpd)
      
      # (ii) Vines per state with {familyset, trunclevel}
      vines_cur <- fit_vines_by(marg_tbl, familyset = famset, trunclevel = trunclevel, min_rows = min_rows_vine)
      
      # (iii) POC (medoid, MJJA test)
      medoid_station <- df_zone %>% filter(medoids==1) %>% pull(codigo_estacion) %>% unique()
      dates_test <- df_zone %>%
        filter(Year==test_year, month(dates)%in%months_MJJA, codigo_estacion==medoid_station) %>%
        pull(dates) %>% unique() %>% sort()
      
      sim_list <- vector("list", length(dates_test))
      for(i in seq_along(dates_test)){
        d0 <- dates_test[i]
        st <- df_zone %>% filter(dates==d0, codigo_estacion==medoid_station) %>% pull(States) %>% unique()
        med_model <- marg_tbl %>% filter(codigo_estacion==medoid_station, States==st) %>% pull(model)
        med_model <- if(length(med_model)==0) NULL else med_model[[1]]
        vm <- vines_cur[[as.character(st)]]
        
        if(is.null(med_model)){
          sim_list[[i]] <- tibble(dates=d0, sim=rep(NA_real_, M_per_day_s9))
          next
        }
        if(!is.null(vm) && !is.null(vm$vine)){
          Udraw <- VineCopula::RVineSim(M_per_day_s9, vm$vine)
          mcol  <- match(as.character(medoid_station), vm$stations)
          u_vec <- if(!is.na(mcol)) Udraw[,mcol] else runif(M_per_day_s9)
        } else {
          u_vec <- runif(M_per_day_s9)
        }
        u_vec <- clamp01(u_vec)
        r_vec <- vapply(u_vec, Q_hurdle, numeric(1), model = med_model)
        sim_list[[i]] <- tibble(dates=rep(d0, M_per_day_s9), sim=r_vec)
      }
      sim_df <- bind_rows(sim_list) %>% arrange(dates)
      obs_df <- df_zone %>% filter(codigo_estacion==medoid_station) %>%
        select(dates, precip_orig = precipitation)
      
      poc_tbl <- sim_df %>%
        group_by(dates) %>%
        summarise(q05 = quantile(sim,0.05,na.rm=TRUE),
                  q95 = quantile(sim,0.95,na.rm=TRUE), .groups="drop") %>%
        left_join(obs_df, by="dates") %>%
        filter(year(dates)==test_year, month(dates)%in%months_MJJA)
      
      POC <- mean((poc_tbl$precip_orig >= poc_tbl$q05) & (poc_tbl$precip_orig <= poc_tbl$q95),
                  na.rm=TRUE) * 100
      
      # (iv) Multivariate ES/VS/CRPS (MJJA averages over states with d >= 2)
      gains_all <- list()
      for(st in sort(unique(df_zone$States))){
        S <- get_fixed_station_set_by(vines_cur, marg_tbl, state_id = st, min_days_frac = 0.95)
        if(length(S) < 2) next
        dates_s <- df_zone %>%
          filter(Year==test_year, month(dates)%in%months_MJJA, States==st) %>%
          distinct(dates) %>% pull()
        day_gains <- purrr::map_dfr(dates_s, ~ score_one_day_by(.x, st, S, vines_cur, marg_tbl,
                                                                M=300, alpha_vs=alpha_vs,
                                                                k_wet_frac=k_wet_frac, thr_joint=thr_pair))
        if(nrow(day_gains)>0) gains_all[[length(gains_all)+1]] <- day_gains
      }
      gains_df <- if(length(gains_all)==0) tibble(ES_gain=NA_real_, VS_gain=NA_real_, CRPS_gain=NA_real_) else bind_rows(gains_all)
      
      sens_results[[length(sens_results)+1]] <- tibble(
        family = flab, trunc = trunclevel, u_star = u_star_cur,
        POC = POC,
        dES = mean(gains_df$ES_gain, na.rm=TRUE),
        dVS = mean(gains_df$VS_gain, na.rm=TRUE),
        dCRPSA = mean(gains_df$CRPS_gain, na.rm=TRUE)
      )
    }
  }
}

sens_summary <- bind_rows(sens_results)

# 4) Sensitivity panel plot
lab_metric <- c(POC = "POC (% in 90% band)",
                dES = "\u0394ES (indep - vine)",
                dVS = "\u0394VS (indep - vine)",
                dCRPSA = "\u0394CRPS_areal (indep - vine)")

plot_s9 <- sens_summary %>%
  pivot_longer(cols = c(POC, dES, dVS, dCRPSA), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("POC","dES","dVS","dCRPSA"),
                         labels = lab_metric)) %>%
  ggplot(aes(x = u_star, y = value, group = factor(trunc), shape = factor(trunc))) +
  geom_hline(data = ~ subset(.x, grepl("^\u0394", as.character(metric))), yintercept = 0, linetype = 2) +
  geom_line() + geom_point(size = 2) +
  facet_grid(metric ~ family, scales = "free_y") +
  scale_x_continuous(breaks = u_vals) +
  labs(x = expression(u^"*"), y = NULL, shape = "trunc") +
  theme_bw() +
  theme(strip.text = element_text(face="bold"))

print(plot_s9)

# 5) Save to supplement
save_fig(plot_s9, base_name = "S9_sensitivity_panels", to_main = FALSE, to_supp = TRUE,
         width = 10, height = 10)



