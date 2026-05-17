# =============================================================================
# PBPK_MRL_cattle_Adult.R  -- INTEGRATED PIPELINE
#
# Combines two stages in a single self-contained R script:
#
#   PART A: PBPK model for iAs in beef cattle (Hung 2021 control dose)
#           - Deterministic ODE + 10000-run population Monte Carlo
#           - Output: Transfer Factors (TF) for muscle/liver/kidney
#
#   PART B: Build TF_cattle distribution (mean + 95% CI) from PBPK output
#
#   PART C: MRL risk assessment (Adult-only, lifetime = 70 years)
#           - 10000-trial Monte Carlo cancer risk model
#           - Output: Maximum Residue Limit (MRL) of iAs in cattle feed
#
# Adult-only design (no age stratification):
#   ED_age = c(Adult = 70)
#   LT     = 70 years
#   Data is pre-filtered to agegroup = "Adult" so no ED warnings appear.
#
# All outputs land in:
#   <this-script-dir>/output/
# =============================================================================


## 0. Libraries and output paths -----------------------------------------------

suppressPackageStartupMessages({
  # PBPK
  library(deSolve)
  library(tidyverse)
  library(parallel)
  library(gridExtra)
  library(grid)
  # MRL
  library(readxl)
  library(MASS)
  library(truncnorm)
  library(scales)
})

# Resolve dplyr conflicts (MASS::select shadows dplyr::select)
select <- dplyr::select
filter <- dplyr::filter

# --- Path resolution (portable) ----------------------------------------------
# The script auto-locates the repository root assuming this script lives under
# <repo>/scripts/  and data sits at  <repo>/data/  -- output goes to <repo>/output/.
# Works for: source(), Rscript, and RStudio.

resolve_script_dir <- function() {
  # 1. Rscript
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(file_arg) > 0) return(normalizePath(dirname(file_arg)))
  # 2. RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    path <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                     error = function(e) "")
    if (nzchar(path)) return(normalizePath(dirname(path)))
  }
  # 3. Fallback: current working directory
  normalizePath(getwd())
}

SCRIPT_DIR <- resolve_script_dir()
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, ".."))
OUT_DIR    <- file.path(REPO_ROOT, "output")
EXCEL_PATH <- file.path(REPO_ROOT, "data", "food_intake_iAs.xlsx")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_path <- function(filename) file.path(OUT_DIR, filename)

cat("INTEGRATED PIPELINE: PBPK -> TF -> MRL  (Adult only)\n")
cat("=====================================================\n")
cat("Output directory:", OUT_DIR, "\n\n")


# =============================================================================
# =============================================================================
#                            PART A: PBPK MODEL
# =============================================================================
# =============================================================================


## A1. PBPK parameters ---------------------------------------------------------

parms <- c(

  # --- Physiological (Lin 2020) ---
  BW  = 621, QCC = 5.45,

  # Fractional blood flows
  FQ_muscle = 0.28, FQ_kid = 0.11, FQ_hep = 0.07, FQ_gi = 0.38, FQ_rest = 0.16,

  # Fractional tissue volumes
  FV_lung = 0.0085, FV_muscle = 0.3610, FV_kid = 0.0021, FV_liv = 0.0122,
  FV_gi = 0.0751, FV_vb = 0.029925, FV_ab = 0.009975, FV_rest = 0.5012,

  MW = 75,

  # Tissue/blood partition coefficients (Yu values)
  P_AsIII_lung   = 4.15, P_AsV_lung   = 4.15, P_MMA_lung   = 1.800, P_DMA_lung   = 2.075,
  P_AsIII_muscle = 2.60, P_AsV_muscle = 2.60, P_MMA_muscle = 1.800, P_DMA_muscle = 2.800,
  P_AsIII_kid    = 4.15, P_AsV_kid    = 4.15, P_MMA_kid    = 1.800, P_DMA_kid    = 2.075,
  P_AsIII_liv    = 5.30, P_AsV_liv    = 5.30, P_MMA_liv    = 2.350, P_DMA_liv    = 2.650,
  P_AsIII_gi     = 2.80, P_AsV_gi     = 2.80, P_MMA_gi     = 1.200, P_DMA_gi     = 1.400,
  P_AsIII_rest   = 2.60, P_AsV_rest   = 2.60, P_MMA_rest   = 1.800, P_DMA_rest   = 2.800,

  # Absorption rate constants (hr-1)
  Ka_AsIII = 1.23e+00, Ka_AsV = 9.25e-01, Ka_MMA = 2.16e+00, Ka_DMA = 2.16e+00,

  # Oxidation / reduction
  K_red_AsV_to_AsIII = 7.04e+00, K_ox_AsIII_to_AsV = 9.41e+00,

  # Methylation - liver
  Vmax_AsIII_to_MMA_liv = 1.60e+02, Km_AsIII_to_MMA_liv = 100,
  Vmax_AsIII_to_DMA_liv = 3.21e+02, Km_AsIII_to_DMA_liv = 100,
  Vmax_MMA_to_DMA_liv   = 2.28e+02, Km_MMA_to_DMA_liv   = 100,

  # Methylation - kidney
  Vmax_AsIII_to_MMA_kid = 1.07e+02, Km_AsIII_to_MMA_kid = 100,
  Vmax_AsIII_to_DMA_kid = 1.43e+02, Km_AsIII_to_DMA_kid = 100,
  Vmax_MMA_to_DMA_kid   = 7.13e+01, Km_MMA_to_DMA_kid   = 100,

  # Urinary, biliary, faecal excretion
  k_urine_AsIII = 2.16e+01, k_urine_AsV = 2.16e+01,
  k_urine_MMA   = 9.25e+01, k_urine_DMA = 4.01e+01,
  eF_AsV = 6.43e-03, eB_AsV = 9.21e-02,

  # --- Oral doses (ug/day) - Hung 2021 iAs_control raw data ---
  # AsIII: 9420 g/day * 14.3 ug/kg DM / 1000 = 134.71 ug/day
  # AsV  : 9420 * 89.0 / 1000 = 838.38 ug/day
  # MMA  : 9420 * 14.4 / 1000 = 135.65 ug/day
  # DMA  : 9420 *  1.8 / 1000 =  16.96 ug/day
  PdoseC_AsIII = 134.71, PdoseC_AsV = 838.38, PdoseC_MMA = 135.65, PdoseC_DMA = 16.96,

  # Simulation schedule
  TSTOP = 200, tlen = 24, FREQ = 24,

  # Feed exposure (BTF/TF normalization)
  feed_intake_kg   = 9.420,   # kg DM/day
  C_feed_iAs_ug_kg = 103.3      # AsIII (14.3) + AsV (89), Hung 2021 control
)


## A2. PBPK helper functions ---------------------------------------------------

pulse_fn <- function(t, period, width) as.numeric(t %% period < width)

rtnorm <- function(n, mean, sd, lo, hi) {
  p_lo <- pnorm(lo, mean, sd)
  p_hi <- pnorm(hi, mean, sd)
  qnorm(runif(n, p_lo, p_hi), mean, sd)
}

compute_derived <- function(p) {
  BW  <- unname(p["BW"]); QCC <- unname(p["QCC"])
  Q_lung   <- QCC * BW
  Q_adjust <- unname(p["FQ_muscle"] + p["FQ_kid"] + p["FQ_hep"] +
                     p["FQ_gi"]    + p["FQ_rest"])
  Q_muscle <- Q_lung * unname(p["FQ_muscle"]) / Q_adjust
  Q_kid    <- Q_lung * unname(p["FQ_kid"])    / Q_adjust
  Q_hep    <- Q_lung * unname(p["FQ_hep"])    / Q_adjust
  Q_gi     <- Q_lung * unname(p["FQ_gi"])     / Q_adjust
  Q_rest   <- Q_lung * unname(p["FQ_rest"])   / Q_adjust
  Q_liv    <- Q_hep + Q_gi

  FV_sum <- unname(p["FV_lung"] + p["FV_muscle"] + p["FV_kid"] + p["FV_liv"] +
                   p["FV_gi"]   + p["FV_vb"]     + p["FV_ab"]  + p["FV_rest"])
  V_adj  <- 1 / FV_sum
  V_lung   <- BW * unname(p["FV_lung"])   * V_adj
  V_muscle <- BW * unname(p["FV_muscle"]) * V_adj
  V_kid    <- BW * unname(p["FV_kid"])    * V_adj
  V_liv    <- BW * unname(p["FV_liv"])    * V_adj
  V_gi     <- BW * unname(p["FV_gi"])     * V_adj
  V_vb     <- BW * unname(p["FV_vb"])     * V_adj
  V_ab     <- BW * unname(p["FV_ab"])     * V_adj
  V_rest   <- BW * unname(p["FV_rest"])   * V_adj

  c(Q_lung = Q_lung, Q_muscle = Q_muscle, Q_kid = Q_kid, Q_hep = Q_hep,
    Q_gi   = Q_gi,   Q_rest   = Q_rest,   Q_liv = Q_liv,
    V_lung = V_lung, V_muscle = V_muscle, V_kid = V_kid, V_liv = V_liv,
    V_gi   = V_gi,   V_vb     = V_vb,     V_ab  = V_ab,  V_rest = V_rest)
}


## A3. PBPK state variable names (52) -----------------------------------------

state_names <- c(
  "Dose_rate_AsIII", "Dose_rate_AsV", "Dose_rate_MMA", "Dose_rate_DMA",
  "AAO_AsIII", "AAO_AsV", "AAO_MMA", "AAO_DMA",
  "MET_AsIII_to_MMA_kid", "MET_AsIII_to_DMA_kid", "MET_MMA_to_DMA_kid",
  "MET_AsIII_to_MMA_liv", "MET_AsIII_to_DMA_liv", "MET_MMA_to_DMA_liv",
  "AMT_AsIII_lung",   "AMT_AsV_lung",   "AMT_MMA_lung",   "AMT_DMA_lung",
  "AMT_AsIII_muscle", "AMT_AsV_muscle", "AMT_MMA_muscle", "AMT_DMA_muscle",
  "AMT_AsIII_kid",    "AMT_AsV_kid",    "AMT_MMA_kid",    "AMT_DMA_kid",
  "AMT_AsIII_liv",    "AMT_AsV_liv",    "AMT_MMA_liv",    "AMT_DMA_liv",
  "AMT_AsIII_gi",     "AMT_AsV_gi",     "AMT_MMA_gi",     "AMT_DMA_gi",
  "AMT_AsIII_rest",   "AMT_AsV_rest",   "AMT_MMA_rest",   "AMT_DMA_rest",
  "AMT_AsIII_vb",     "AMT_AsV_vb",     "AMT_MMA_vb",     "AMT_DMA_vb",
  "AMT_AsIII_ab",     "AMT_AsV_ab",     "AMT_MMA_ab",     "AMT_DMA_ab",
  "AMT_AsIII_urine", "AMT_AsV_urine", "AMT_MMA_urine", "AMT_DMA_urine",
  "AMT_AsV_bile", "AMT_AsV_faecal"
)
stopifnot(length(state_names) == 52)


## A4. PBPK ODE function -------------------------------------------------------

pbpk_cattle_ode <- function(t, y, parms) {
  with(as.list(c(y, parms)), {

    if (t < TSTOP) {
      Exposure    <- pulse_fn(t, FREQ, tlen)
      ROral_AsIII <- (PdoseC_AsIII / MW / tlen) * Exposure
      ROral_AsV   <- (PdoseC_AsV   / MW / tlen) * Exposure
      ROral_MMA   <- (PdoseC_MMA   / MW / tlen) * Exposure
      ROral_DMA   <- (PdoseC_DMA   / MW / tlen) * Exposure
    } else {
      ROral_AsIII <- 0; ROral_AsV <- 0; ROral_MMA <- 0; ROral_DMA <- 0
    }

    dDose_rate_AsIII <- ROral_AsIII - Dose_rate_AsIII * Ka_AsIII
    dDose_rate_AsV   <- ROral_AsV   - Dose_rate_AsV   * Ka_AsV
    dDose_rate_MMA   <- ROral_MMA   - Dose_rate_MMA   * Ka_MMA
    dDose_rate_DMA   <- ROral_DMA   - Dose_rate_DMA   * Ka_DMA

    dAAO_AsIII <- Dose_rate_AsIII * Ka_AsIII
    dAAO_AsV   <- Dose_rate_AsV   * Ka_AsV
    dAAO_MMA   <- Dose_rate_MMA   * Ka_MMA
    dAAO_DMA   <- Dose_rate_DMA   * Ka_DMA

    # Concentrations
    C_AsIII_lung <- AMT_AsIII_lung / V_lung; Ca_AsIII_lung <- C_AsIII_lung / P_AsIII_lung
    C_AsV_lung   <- AMT_AsV_lung   / V_lung; Ca_AsV_lung   <- C_AsV_lung   / P_AsV_lung
    C_MMA_lung   <- AMT_MMA_lung   / V_lung; Ca_MMA_lung   <- C_MMA_lung   / P_MMA_lung
    C_DMA_lung   <- AMT_DMA_lung   / V_lung; Ca_DMA_lung   <- C_DMA_lung   / P_DMA_lung

    C_AsIII_muscle <- AMT_AsIII_muscle / V_muscle; CV_AsIII_muscle <- C_AsIII_muscle / P_AsIII_muscle
    C_AsV_muscle   <- AMT_AsV_muscle   / V_muscle; CV_AsV_muscle   <- C_AsV_muscle   / P_AsV_muscle
    C_MMA_muscle   <- AMT_MMA_muscle   / V_muscle; CV_MMA_muscle   <- C_MMA_muscle   / P_MMA_muscle
    C_DMA_muscle   <- AMT_DMA_muscle   / V_muscle; CV_DMA_muscle   <- C_DMA_muscle   / P_DMA_muscle

    C_AsIII_kid <- AMT_AsIII_kid / V_kid; CV_AsIII_kid <- C_AsIII_kid / P_AsIII_kid
    C_AsV_kid   <- AMT_AsV_kid   / V_kid; CV_AsV_kid   <- C_AsV_kid   / P_AsV_kid
    C_MMA_kid   <- AMT_MMA_kid   / V_kid; CV_MMA_kid   <- C_MMA_kid   / P_MMA_kid
    C_DMA_kid   <- AMT_DMA_kid   / V_kid; CV_DMA_kid   <- C_DMA_kid   / P_DMA_kid

    C_AsIII_liv <- AMT_AsIII_liv / V_liv; CV_AsIII_liv <- C_AsIII_liv / P_AsIII_liv
    C_AsV_liv   <- AMT_AsV_liv   / V_liv; CV_AsV_liv   <- C_AsV_liv   / P_AsV_liv
    C_MMA_liv   <- AMT_MMA_liv   / V_liv; CV_MMA_liv   <- C_MMA_liv   / P_MMA_liv
    C_DMA_liv   <- AMT_DMA_liv   / V_liv; CV_DMA_liv   <- C_DMA_liv   / P_DMA_liv

    C_AsIII_gi <- AMT_AsIII_gi / V_gi; CV_AsIII_gi <- C_AsIII_gi / P_AsIII_gi
    C_AsV_gi   <- AMT_AsV_gi   / V_gi; CV_AsV_gi   <- C_AsV_gi   / P_AsV_gi
    C_MMA_gi   <- AMT_MMA_gi   / V_gi; CV_MMA_gi   <- C_MMA_gi   / P_MMA_gi
    C_DMA_gi   <- AMT_DMA_gi   / V_gi; CV_DMA_gi   <- C_DMA_gi   / P_DMA_gi

    C_AsIII_rest <- AMT_AsIII_rest / V_rest; CV_AsIII_rest <- C_AsIII_rest / P_AsIII_rest
    C_AsV_rest   <- AMT_AsV_rest   / V_rest; CV_AsV_rest   <- C_AsV_rest   / P_AsV_rest
    C_MMA_rest   <- AMT_MMA_rest   / V_rest; CV_MMA_rest   <- C_MMA_rest   / P_MMA_rest
    C_DMA_rest   <- AMT_DMA_rest   / V_rest; CV_DMA_rest   <- C_DMA_rest   / P_DMA_rest

    C_AsIII_vb <- AMT_AsIII_vb / V_vb; C_AsV_vb <- AMT_AsV_vb / V_vb
    C_MMA_vb   <- AMT_MMA_vb   / V_vb; C_DMA_vb <- AMT_DMA_vb / V_vb
    C_AsIII_ab <- AMT_AsIII_ab / V_ab; C_AsV_ab <- AMT_AsV_ab / V_ab
    C_MMA_ab   <- AMT_MMA_ab   / V_ab; C_DMA_ab <- AMT_DMA_ab / V_ab

    # Methylation fluxes
    RMET_AsIII_to_MMA_kid <- Vmax_AsIII_to_MMA_kid * C_AsIII_kid / (Km_AsIII_to_MMA_kid + C_AsIII_kid)
    RMET_AsIII_to_DMA_kid <- Vmax_AsIII_to_DMA_kid * C_AsIII_kid / (Km_AsIII_to_DMA_kid + C_AsIII_kid)
    RMET_MMA_to_DMA_kid   <- Vmax_MMA_to_DMA_kid   * C_MMA_kid   / (Km_MMA_to_DMA_kid   + C_MMA_kid)
    RMET_AsIII_to_MMA_liv <- Vmax_AsIII_to_MMA_liv * C_AsIII_liv / (Km_AsIII_to_MMA_liv + C_AsIII_liv)
    RMET_AsIII_to_DMA_liv <- Vmax_AsIII_to_DMA_liv * C_AsIII_liv / (Km_AsIII_to_DMA_liv + C_AsIII_liv)
    RMET_MMA_to_DMA_liv   <- Vmax_MMA_to_DMA_liv   * C_MMA_liv   / (Km_MMA_to_DMA_liv   + C_MMA_liv)

    dMET_AsIII_to_MMA_kid <- RMET_AsIII_to_MMA_kid
    dMET_AsIII_to_DMA_kid <- RMET_AsIII_to_DMA_kid
    dMET_MMA_to_DMA_kid   <- RMET_MMA_to_DMA_kid
    dMET_AsIII_to_MMA_liv <- RMET_AsIII_to_MMA_liv
    dMET_AsIII_to_DMA_liv <- RMET_AsIII_to_DMA_liv
    dMET_MMA_to_DMA_liv   <- RMET_MMA_to_DMA_liv

    # Tissue mass balance
    dAMT_AsIII_lung <- Q_lung*C_AsIII_vb - Q_lung*Ca_AsIII_lung +
      K_red_AsV_to_AsIII*AMT_AsV_lung - K_ox_AsIII_to_AsV*AMT_AsIII_lung
    dAMT_AsV_lung   <- Q_lung*C_AsV_vb   - Q_lung*Ca_AsV_lung   -
      K_red_AsV_to_AsIII*AMT_AsV_lung + K_ox_AsIII_to_AsV*AMT_AsIII_lung
    dAMT_MMA_lung   <- Q_lung*C_MMA_vb - Q_lung*Ca_MMA_lung
    dAMT_DMA_lung   <- Q_lung*C_DMA_vb - Q_lung*Ca_DMA_lung

    dAMT_AsIII_muscle <- Q_muscle*C_AsIII_ab - Q_muscle*CV_AsIII_muscle +
      K_red_AsV_to_AsIII*AMT_AsV_muscle - K_ox_AsIII_to_AsV*AMT_AsIII_muscle
    dAMT_AsV_muscle   <- Q_muscle*C_AsV_ab   - Q_muscle*CV_AsV_muscle   -
      K_red_AsV_to_AsIII*AMT_AsV_muscle + K_ox_AsIII_to_AsV*AMT_AsIII_muscle
    dAMT_MMA_muscle   <- Q_muscle*C_MMA_ab - Q_muscle*CV_MMA_muscle
    dAMT_DMA_muscle   <- Q_muscle*C_DMA_ab - Q_muscle*CV_DMA_muscle

    dAMT_AsIII_kid <- Q_kid*C_AsIII_ab - Q_kid*CV_AsIII_kid +
      K_red_AsV_to_AsIII*AMT_AsV_kid - K_ox_AsIII_to_AsV*AMT_AsIII_kid -
      RMET_AsIII_to_MMA_kid - RMET_AsIII_to_DMA_kid -
      k_urine_AsIII*(AMT_AsIII_kid / P_AsIII_kid)
    dAMT_AsV_kid   <- Q_kid*C_AsV_ab - Q_kid*CV_AsV_kid -
      K_red_AsV_to_AsIII*AMT_AsV_kid + K_ox_AsIII_to_AsV*AMT_AsIII_kid -
      k_urine_AsV*(AMT_AsV_kid / P_AsV_kid)
    dAMT_MMA_kid   <- Q_kid*C_MMA_ab - Q_kid*CV_MMA_kid +
      RMET_AsIII_to_MMA_kid - RMET_MMA_to_DMA_kid -
      k_urine_MMA*(AMT_MMA_kid / P_MMA_kid)
    dAMT_DMA_kid   <- Q_kid*C_DMA_ab - Q_kid*CV_DMA_kid +
      RMET_AsIII_to_DMA_kid + RMET_MMA_to_DMA_kid -
      k_urine_DMA*(AMT_DMA_kid / P_DMA_kid)

    dAMT_AsIII_liv <- Q_hep*C_AsIII_ab - Q_liv*CV_AsIII_liv + Q_gi*CV_AsIII_gi +
      K_red_AsV_to_AsIII*AMT_AsV_liv - K_ox_AsIII_to_AsV*AMT_AsIII_liv -
      RMET_AsIII_to_MMA_liv - RMET_AsIII_to_DMA_liv
    dAMT_AsV_liv   <- Q_hep*C_AsV_ab - Q_liv*CV_AsV_liv + Q_gi*CV_AsV_gi -
      K_red_AsV_to_AsIII*AMT_AsV_liv + K_ox_AsIII_to_AsV*AMT_AsIII_liv -
      eB_AsV*(AMT_AsV_liv / P_AsV_liv)
    dAMT_MMA_liv   <- Q_hep*C_MMA_ab - Q_liv*CV_MMA_liv + Q_gi*CV_MMA_gi +
      RMET_AsIII_to_MMA_liv - RMET_MMA_to_DMA_liv
    dAMT_DMA_liv   <- Q_hep*C_DMA_ab - Q_liv*CV_DMA_liv + Q_gi*CV_DMA_gi +
      RMET_AsIII_to_DMA_liv + RMET_MMA_to_DMA_liv

    dAMT_AsIII_gi <- Q_gi*C_AsIII_ab - Q_gi*CV_AsIII_gi +
      K_red_AsV_to_AsIII*AMT_AsV_gi - K_ox_AsIII_to_AsV*AMT_AsIII_gi +
      Dose_rate_AsIII*Ka_AsIII
    dAMT_AsV_gi   <- Q_gi*C_AsV_ab - Q_gi*CV_AsV_gi -
      K_red_AsV_to_AsIII*AMT_AsV_gi + K_ox_AsIII_to_AsV*AMT_AsIII_gi -
      eF_AsV*(AMT_AsV_gi / P_AsV_gi) + Dose_rate_AsV*Ka_AsV
    dAMT_MMA_gi   <- Q_gi*C_MMA_ab - Q_gi*CV_MMA_gi + Dose_rate_MMA*Ka_MMA
    dAMT_DMA_gi   <- Q_gi*C_DMA_ab - Q_gi*CV_DMA_gi + Dose_rate_DMA*Ka_DMA

    dAMT_AsIII_rest <- Q_rest*C_AsIII_ab - Q_rest*CV_AsIII_rest +
      K_red_AsV_to_AsIII*AMT_AsV_rest - K_ox_AsIII_to_AsV*AMT_AsIII_rest
    dAMT_AsV_rest   <- Q_rest*C_AsV_ab - Q_rest*CV_AsV_rest -
      K_red_AsV_to_AsIII*AMT_AsV_rest + K_ox_AsIII_to_AsV*AMT_AsIII_rest
    dAMT_MMA_rest   <- Q_rest*C_MMA_ab - Q_rest*CV_MMA_rest
    dAMT_DMA_rest   <- Q_rest*C_DMA_ab - Q_rest*CV_DMA_rest

    dAMT_AsIII_vb <- Q_muscle*CV_AsIII_muscle + Q_kid*CV_AsIII_kid +
      Q_liv*CV_AsIII_liv + Q_rest*CV_AsIII_rest - Q_lung*C_AsIII_vb
    dAMT_AsV_vb   <- Q_muscle*CV_AsV_muscle   + Q_kid*CV_AsV_kid   +
      Q_liv*CV_AsV_liv   + Q_rest*CV_AsV_rest   - Q_lung*C_AsV_vb
    dAMT_MMA_vb   <- Q_muscle*CV_MMA_muscle   + Q_kid*CV_MMA_kid   +
      Q_liv*CV_MMA_liv   + Q_rest*CV_MMA_rest   - Q_lung*C_MMA_vb
    dAMT_DMA_vb   <- Q_muscle*CV_DMA_muscle   + Q_kid*CV_DMA_kid   +
      Q_liv*CV_DMA_liv   + Q_rest*CV_DMA_rest   - Q_lung*C_DMA_vb

    Q_out <- Q_muscle + Q_kid + Q_hep + Q_gi + Q_rest
    dAMT_AsIII_ab <- Q_lung*Ca_AsIII_lung - Q_out*C_AsIII_ab +
      K_red_AsV_to_AsIII*AMT_AsV_ab - K_ox_AsIII_to_AsV*AMT_AsIII_ab
    dAMT_AsV_ab   <- Q_lung*Ca_AsV_lung   - Q_out*C_AsV_ab   -
      K_red_AsV_to_AsIII*AMT_AsV_ab + K_ox_AsIII_to_AsV*AMT_AsIII_ab
    dAMT_MMA_ab   <- Q_lung*Ca_MMA_lung - Q_out*C_MMA_ab
    dAMT_DMA_ab   <- Q_lung*Ca_DMA_lung - Q_out*C_DMA_ab

    dAMT_AsIII_urine <- k_urine_AsIII*(AMT_AsIII_kid / P_AsIII_kid)
    dAMT_AsV_urine   <- k_urine_AsV  *(AMT_AsV_kid   / P_AsV_kid)
    dAMT_MMA_urine   <- k_urine_MMA  *(AMT_MMA_kid   / P_MMA_kid)
    dAMT_DMA_urine   <- k_urine_DMA  *(AMT_DMA_kid   / P_DMA_kid)
    dAMT_AsV_bile    <- eB_AsV       *(AMT_AsV_liv   / P_AsV_liv)
    dAMT_AsV_faecal  <- eF_AsV       *(AMT_AsV_gi    / P_AsV_gi)

    list(c(
      dDose_rate_AsIII, dDose_rate_AsV, dDose_rate_MMA, dDose_rate_DMA,
      dAAO_AsIII, dAAO_AsV, dAAO_MMA, dAAO_DMA,
      dMET_AsIII_to_MMA_kid, dMET_AsIII_to_DMA_kid, dMET_MMA_to_DMA_kid,
      dMET_AsIII_to_MMA_liv, dMET_AsIII_to_DMA_liv, dMET_MMA_to_DMA_liv,
      dAMT_AsIII_lung,   dAMT_AsV_lung,   dAMT_MMA_lung,   dAMT_DMA_lung,
      dAMT_AsIII_muscle, dAMT_AsV_muscle, dAMT_MMA_muscle, dAMT_DMA_muscle,
      dAMT_AsIII_kid,    dAMT_AsV_kid,    dAMT_MMA_kid,    dAMT_DMA_kid,
      dAMT_AsIII_liv,    dAMT_AsV_liv,    dAMT_MMA_liv,    dAMT_DMA_liv,
      dAMT_AsIII_gi,     dAMT_AsV_gi,     dAMT_MMA_gi,     dAMT_DMA_gi,
      dAMT_AsIII_rest,   dAMT_AsV_rest,   dAMT_MMA_rest,   dAMT_DMA_rest,
      dAMT_AsIII_vb,     dAMT_AsV_vb,     dAMT_MMA_vb,     dAMT_DMA_vb,
      dAMT_AsIII_ab,     dAMT_AsV_ab,     dAMT_MMA_ab,     dAMT_DMA_ab,
      dAMT_AsIII_urine, dAMT_AsV_urine, dAMT_MMA_urine, dAMT_DMA_urine,
      dAMT_AsV_bile, dAMT_AsV_faecal
    ))
  })
}


## A5. Deterministic PBPK run --------------------------------------------------

cat("\n[PART A] PBPK MODEL\n")
cat("-------------------\n")
cat("A1. Running deterministic PBPK...\n")

derived_det <- compute_derived(parms)
y0          <- setNames(rep(0, length(state_names)), state_names)
times       <- seq(0, parms["TSTOP"], by = 1)

out <- ode(
  y     = y0,    times  = times,
  func  = pbpk_cattle_ode,
  parms = c(parms, derived_det),
  method = "lsoda"
)
out_df <- as.data.frame(out)
cat(sprintf("    Done: %d time points, %d state variables\n",
            nrow(out_df), ncol(out_df) - 1))


## A6. Mass balance assertions -------------------------------------------------

cat("A2. Checking mass balance...\n")

# Concentrations + body burden for the assertion below
out_df <- with(as.list(c(parms, derived_det)), {
  out_df |> mutate(
    Day = time / 24,
    C_AsIII_muscle = AMT_AsIII_muscle / V_muscle,
    C_AsV_muscle   = AMT_AsV_muscle   / V_muscle,
    C_iAs_muscle   = C_AsIII_muscle + C_AsV_muscle,
    C_AsIII_kid    = AMT_AsIII_kid    / V_kid,
    C_AsV_kid      = AMT_AsV_kid      / V_kid,
    C_iAs_kid      = C_AsIII_kid + C_AsV_kid,
    C_AsIII_liv    = AMT_AsIII_liv    / V_liv,
    C_AsV_liv      = AMT_AsV_liv      / V_liv,
    C_iAs_liv      = C_AsIII_liv + C_AsV_liv,
    AMT_AsIII = AMT_AsIII_ab + AMT_AsIII_vb + AMT_AsIII_lung + AMT_AsIII_liv +
                AMT_AsIII_gi + AMT_AsIII_kid + AMT_AsIII_muscle + AMT_AsIII_rest,
    AMT_AsV   = AMT_AsV_ab   + AMT_AsV_vb   + AMT_AsV_lung   + AMT_AsV_liv +
                AMT_AsV_gi   + AMT_AsV_kid   + AMT_AsV_muscle + AMT_AsV_rest,
    AMT_MMA   = AMT_MMA_ab   + AMT_MMA_vb   + AMT_MMA_lung   + AMT_MMA_liv +
                AMT_MMA_gi   + AMT_MMA_kid   + AMT_MMA_muscle + AMT_MMA_rest,
    AMT_DMA   = AMT_DMA_ab   + AMT_DMA_vb   + AMT_DMA_lung   + AMT_DMA_liv +
                AMT_DMA_gi   + AMT_DMA_kid   + AMT_DMA_muscle + AMT_DMA_rest,
    AMT_urine = AMT_AsIII_urine + AMT_AsV_urine + AMT_MMA_urine + AMT_DMA_urine,
    AMT_excretion = AMT_urine + AMT_AsV_bile + AMT_AsV_faecal,
    AAO_total = AAO_AsIII + AAO_AsV + AAO_MMA + AAO_DMA,
    BAL = AAO_total - (AMT_AsIII + AMT_AsV + AMT_MMA + AMT_DMA + AMT_excretion)
  )
})

bal_A_pct <- max(abs(out_df$BAL)) / max(out_df$AAO_total) * 100

daily_input_umol <- with(as.list(parms),
  (PdoseC_AsIII + PdoseC_AsV + PdoseC_MMA + PdoseC_DMA) / MW)
total_fed_t <- (daily_input_umol / 24) * pmin(out_df$time, unname(parms["TSTOP"]))
depot_sum   <- with(out_df, Dose_rate_AsIII + Dose_rate_AsV +
                            Dose_rate_MMA   + Dose_rate_DMA)
total_accounted <- with(out_df,
  AMT_AsIII + AMT_AsV + AMT_MMA + AMT_DMA + AMT_excretion) + depot_sum
bal_B_pct <- max(abs(total_fed_t - total_accounted)) /
             total_fed_t[length(total_fed_t)] * 100

cat(sprintf("    (A) max |BAL|/AAO    = %.4e %% (tol 1e-08)\n", bal_A_pct))
cat(sprintf("    (B) max |fed-acct|/f = %.4e %% (tol 1e-03)\n", bal_B_pct))
stopifnot("Mass balance (A) violated" = bal_A_pct < 1e-8)
stopifnot("Mass balance (B) violated" = bal_B_pct < 1e-3)
cat("    PASS\n")


## A7. Population PBPK setup ---------------------------------------------------

cat("A3. Population PBPK setup...\n")

pop_dist <- list(
  BW        = c(mean = 621,       sd = 6.21,       lo = 602.37,    hi = 639.63),
  QCC       = c(mean = 5.45,      sd = 0.545,      lo = 3.815,     hi = 7.085),
  FQ_muscle = c(mean = 0.28,      sd = 0.09,       lo = 0.01,      hi = 0.55),
  FQ_kid    = c(mean = 0.11,      sd = 0.011,      lo = 0.077,     hi = 0.143),
  FQ_hep    = c(mean = 0.07,      sd = 0.007,      lo = 0.049,     hi = 0.091),
  FQ_gi     = c(mean = 0.38,      sd = 0.038,      lo = 0.266,     hi = 0.494),
  FQ_rest   = c(mean = 0.16,      sd = 0.016,      lo = 0.0001,    hi = 0.598),
  FV_lung   = c(mean = 0.0085,    sd = 0.0025,     lo = 0.001,     hi = 0.016),
  FV_muscle = c(mean = 0.3610,    sd = 0.1137,     lo = 0.0199,    hi = 0.7021),
  FV_kid    = c(mean = 0.0021,    sd = 0.0005,     lo = 0.0006,    hi = 0.0036),
  FV_liv    = c(mean = 0.0122,    sd = 0.0018,     lo = 0.0068,    hi = 0.0176),
  FV_gi     = c(mean = 0.0751,    sd = 0.00751,    lo = 0.05257,   hi = 0.09763),
  FV_vb     = c(mean = 0.029925,  sd = 0.0029925,  lo = 0.0209475, hi = 0.0389025),
  FV_ab     = c(mean = 0.009975,  sd = 0.0009975,  lo = 0.0069825, hi = 0.0129675),
  FV_rest   = c(mean = 0.5012,    sd = 0.05012,    lo = 0.1112,    hi = 0.8912)
)

# CANONICAL run size - LOCKED for manuscript reproducibility.
# All reported numbers and figures derive from this exact configuration:
#   PBPK population = 10000 | PBPK seed = 4729
#   MRL Monte Carlo = 10000 | MRL  seed = 123  (see config below)
n_runs_pbpk <- 10000L
cat("    n_runs =", n_runs_pbpk, "(locked canonical)\n")

set.seed(4729)
parm_samples <- do.call(cbind, lapply(pop_dist, function(d) {
  rtnorm(n_runs_pbpk, d["mean"], d["sd"], d["lo"], d["hi"])
}))
colnames(parm_samples) <- names(pop_dist)
parms_fixed <- parms[!names(parms) %in% names(pop_dist)]


## A8. Single-run wrapper ------------------------------------------------------

run_one_pop <- function(i) {
  p_i   <- parm_samples[i, ]
  d_i   <- compute_derived(p_i)
  all_p <- c(parms_fixed, p_i, d_i)
  tryCatch({
    out_i <- deSolve::ode(
      y = y0, times = times, func = pbpk_cattle_ode,
      parms = all_p, method = "lsoda"
    )
    df_i <- as.data.frame(out_i)
    data.frame(
      run          = i,
      time         = df_i$time,
      C_iAs_muscle = (df_i$AMT_AsIII_muscle + df_i$AMT_AsV_muscle) / all_p["V_muscle"],
      C_iAs_liv    = (df_i$AMT_AsIII_liv    + df_i$AMT_AsV_liv)    / all_p["V_liv"],
      C_iAs_kid    = (df_i$AMT_AsIII_kid    + df_i$AMT_AsV_kid)    / all_p["V_kid"]
    )
  }, error = function(e) {
    warning(sprintf("Run %d failed: %s", i, e$message)); NULL
  })
}


## A9. Run population in parallel ----------------------------------------------

n_cores <- max(1L, detectCores() - 1L)
cat(sprintf("A4. Population MC: %d runs on %d cores...\n", n_runs_pbpk, n_cores))

t_start <- proc.time()
if (.Platform$OS.type == "windows") {
  cl <- makeCluster(n_cores)
  clusterExport(cl, varlist = c("pbpk_cattle_ode", "pulse_fn", "compute_derived",
                                "y0", "times", "parms_fixed", "parm_samples"))
  invisible(clusterEvalQ(cl, library(deSolve)))
  pop_raw <- parLapply(cl, seq_len(n_runs_pbpk), run_one_pop)
  stopCluster(cl)
} else {
  pop_raw <- mclapply(seq_len(n_runs_pbpk), run_one_pop, mc.cores = n_cores)
}
elapsed <- (proc.time() - t_start)["elapsed"]
n_ok <- sum(!sapply(pop_raw, is.null))
cat(sprintf("    Done: %d/%d runs in %.1f s\n", n_ok, n_runs_pbpk, elapsed))

pop_df <- do.call(rbind, Filter(Negate(is.null), pop_raw))


## A10. Convert to ug/kg and compute BTF/TF -----------------------------------

MW_As            <- unname(parms["MW"])
feed_intake_kg   <- unname(parms["feed_intake_kg"])
C_feed_iAs_ug_kg <- unname(parms["C_feed_iAs_ug_kg"])

pop_df_conv <- pop_df |>
  mutate(
    Day             = time / 24,
    C_iAs_muscle_ug = C_iAs_muscle * MW_As,
    C_iAs_liv_ug    = C_iAs_liv    * MW_As,
    C_iAs_kid_ug    = C_iAs_kid    * MW_As,
    BTF_muscle      = C_iAs_muscle_ug / C_feed_iAs_ug_kg,
    BTF_liv         = C_iAs_liv_ug    / C_feed_iAs_ug_kg,
    BTF_kid         = C_iAs_kid_ug    / C_feed_iAs_ug_kg,
    TF_muscle       = BTF_muscle / feed_intake_kg,
    TF_liv          = BTF_liv    / feed_intake_kg,
    TF_kid          = BTF_kid    / feed_intake_kg
  )


## A11. PBPK plots -------------------------------------------------------------

cat("A5. Saving PBPK plots and CSVs...\n")

# Plot: Deterministic species in liver/kidney
det_long <- out_df |>
  select(Day,
         Kid_AsIII = C_AsIII_kid, Kid_AsV = C_AsV_kid,
         Liv_AsIII = C_AsIII_liv, Liv_AsV = C_AsV_liv,
         Mus_AsIII = C_AsIII_muscle, Mus_AsV = C_AsV_muscle) |>
  pivot_longer(-Day, names_to = c("Tissue", "Species"), names_sep = "_",
               values_to = "Conc") |>
  mutate(Tissue = recode(Tissue, Kid = "Kidney", Liv = "Liver", Mus = "Muscle"))

p_det <- ggplot(det_long, aes(x = Day, y = Conc, colour = Species)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~Tissue, scales = "free_y") +
  labs(title = "PBPK deterministic: AsIII and AsV in cattle tissues",
       x = "Time (days)", y = expression(C~(mu*mol~L^{-1}))) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(out_path("A1_pbpk_deterministic_species.png"),
       p_det, width = 10, height = 5, dpi = 300)

# Mass balance plot - fractional error (% of cumulative absorbed iAs)
# Plotting the relative error (rather than absolute umol) shows directly that
# the error stays far below the 1% regulatory criterion. t = 0 is dropped to
# avoid 0/0 (nothing absorbed yet).
out_df$BAL_pct <- ifelse(out_df$AAO_total > 0,
                         out_df$BAL / out_df$AAO_total * 100, NA_real_)
max_bal_pct <- max(abs(out_df$BAL_pct), na.rm = TRUE)

p_bal <- ggplot(subset(out_df, !is.na(BAL_pct)), aes(Day, BAL_pct)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(colour = "firebrick", linewidth = 0.8) +
  labs(title    = "PBPK mass-balance fractional error",
       subtitle = sprintf("max |error| = %.2e %%   (regulatory criterion: < 1%%)",
                           max_bal_pct),
       x = "Time (days)",
       y = "Mass-balance error (% of cumulative absorbed iAs)") +
  theme_bw(base_size = 11)
ggsave(out_path("A2_pbpk_mass_balance.png"),
       p_bal, width = 8, height = 4, dpi = 300)

# Population ribbon plot (ug/kg ww)
pct_fns_pbpk <- list(
  p05 = \(x) quantile(x, 0.05), p25 = \(x) quantile(x, 0.25),
  p50 = \(x) quantile(x, 0.50), p75 = \(x) quantile(x, 0.75),
  p95 = \(x) quantile(x, 0.95)
)

pct_pop <- pop_df_conv |>
  group_by(time, Day) |>
  summarise(across(c(C_iAs_muscle_ug, C_iAs_liv_ug, C_iAs_kid_ug),
                   pct_fns_pbpk, .names = "{.col}_{.fn}"),
            .groups = "drop")

pop_rib <- pct_pop |>
  select(Day,
         Muscle_p05 = C_iAs_muscle_ug_p05, Muscle_p25 = C_iAs_muscle_ug_p25,
         Muscle_p50 = C_iAs_muscle_ug_p50, Muscle_p75 = C_iAs_muscle_ug_p75,
         Muscle_p95 = C_iAs_muscle_ug_p95,
         Liver_p05  = C_iAs_liv_ug_p05,    Liver_p25  = C_iAs_liv_ug_p25,
         Liver_p50  = C_iAs_liv_ug_p50,    Liver_p75  = C_iAs_liv_ug_p75,
         Liver_p95  = C_iAs_liv_ug_p95,
         Kidney_p05 = C_iAs_kid_ug_p05,    Kidney_p25 = C_iAs_kid_ug_p25,
         Kidney_p50 = C_iAs_kid_ug_p50,    Kidney_p75 = C_iAs_kid_ug_p75,
         Kidney_p95 = C_iAs_kid_ug_p95) |>
  pivot_longer(-Day, names_to = c("Tissue", "Pct"), names_sep = "_",
               values_to = "Conc") |>
  pivot_wider(names_from = Pct, values_from = Conc) |>
  mutate(Tissue = factor(Tissue, levels = c("Kidney", "Liver", "Muscle")))

p_pop <- ggplot(pop_rib, aes(x = Day)) +
  geom_ribbon(aes(ymin = p05, ymax = p95), fill = "#4393c3", alpha = 0.18) +
  geom_ribbon(aes(ymin = p25, ymax = p75), fill = "#4393c3", alpha = 0.38) +
  geom_line(aes(y = p50), colour = "steelblue4", linewidth = 0.9) +
  facet_wrap(~Tissue, scales = "free_y", ncol = 1) +
  labs(title = sprintf("Population PBPK iAs (n = %d, Hung 2021 dose)", n_runs_pbpk),
       x = "Time (days)", y = expression(C[iAs]~(mu*g~kg^{-1}~ww))) +
  theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"))
ggsave(out_path("A3_pbpk_population_ribbon.png"),
       p_pop, width = 8, height = 10, dpi = 300)


## A12. Publication-grade tissue distribution figures ------------------------
# Two figures using Okabe-Ito colorblind-safe palette, saved as TIFF (300 dpi,
# LZW compression) for journal submission plus PNG for quick viewing.
#   A4: Time-course (median + 95% CI ribbon, single panel, all tissues overlaid)
#   A5: Steady-state distribution (violin + jitter + boxplot)
# Concentrations are taken from C_iAs_*_ug (already in ug/kg ww using MW=75
# consistent with the rest of the script).

cat("A6. Publication-grade tissue plots (TIFF + PNG)...\n")

TISSUE_COLORS <- c(
  "Liver"                = "#D55E00",
  "Kidney"               = "#E69F00",
  "Muscle - Other offal" = "#0072B2"
)
TISSUE_LEVELS <- c("Liver", "Kidney", "Muscle - Other offal")

# --- Long-format population data, with negative-artifact filter --------------
# Tiny negative values can appear at solver transients; report and remove them.
n_neg <- sum(pop_df_conv$C_iAs_muscle_ug < 0 |
             pop_df_conv$C_iAs_liv_ug    < 0 |
             pop_df_conv$C_iAs_kid_ug    < 0, na.rm = TRUE)
cat(sprintf("    Numerical artifacts removed: %d / %d rows (%.4f%%)\n",
            n_neg, nrow(pop_df_conv), 100 * n_neg / nrow(pop_df_conv)))
if (n_neg > 0) {
  t_neg <- pop_df_conv |>
    filter(C_iAs_muscle_ug < 0 | C_iAs_liv_ug < 0 | C_iAs_kid_ug < 0) |>
    pull(time) |> range()
  cat(sprintf("    Range of t with negatives: %.2f - %.2f hours\n",
              t_neg[1], t_neg[2]))
}

raw_long <- pop_df_conv |>
  filter(C_iAs_muscle_ug >= 0, C_iAs_liv_ug >= 0, C_iAs_kid_ug >= 0) |>
  select(t = time,
         `Muscle - Other offal` = C_iAs_muscle_ug,
         Liver                  = C_iAs_liv_ug,
         Kidney                 = C_iAs_kid_ug) |>
  pivot_longer(-t, names_to = "Tissue", values_to = "Concentration") |>
  mutate(Tissue = factor(Tissue, levels = TISSUE_LEVELS))

# --- A4: Publication time-course (median + 95% CI ribbon) -------------------
ts_cattle <- raw_long |>
  group_by(t, Tissue) |>
  summarise(
    Median = median(Concentration, na.rm = TRUE),
    Low    = quantile(Concentration, 0.025, na.rm = TRUE),
    High   = quantile(Concentration, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

p_tc <- ggplot(ts_cattle, aes(x = t, group = Tissue)) +
  geom_ribbon(aes(ymin = Low, ymax = High, fill = Tissue),
              alpha = 0.15, color = NA) +
  geom_line(aes(y = Median, color = Tissue, linetype = Tissue),
            linewidth = 1) +
  scale_color_manual(values = TISSUE_COLORS, name = "Tissue") +
  scale_fill_manual(values = TISSUE_COLORS, name = "Tissue") +
  scale_linetype_manual(values = c("solid", "longdash", "solid"),
                        name = "Tissue") +
  scale_x_continuous(
    name   = "Exposure time (hours)",
    breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_y_continuous(
    name   = expression(paste("iAs concentration (", mu, "g/kg)")),
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, NA)
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    axis.title  = element_text(color = "black", size = 14),
    axis.text   = element_text(color = "black", size = 12),
    axis.line   = element_line(color = "black", linewidth = 0.5),
    plot.margin = margin(t = 15, r = 20, b = 10, l = 10, unit = "pt"),
    legend.position   = c(0.20, 0.85),
    legend.background = element_rect(fill = "white", color = "black",
                                     linewidth = 0.3),
    legend.title      = element_blank(),
    legend.text       = element_text(size = 11)
  )

ggsave(out_path("A4_Cattle_Tissue_TimeCourse.tiff"), p_tc,
       width = 8, height = 6, dpi = 300, compression = "lzw")
ggsave(out_path("A4_Cattle_Tissue_TimeCourse.png"),  p_tc,
       width = 8, height = 6, dpi = 300)

# --- A5: Publication steady-state distribution (violin + jitter + boxplot) --
ss_cattle <- raw_long |> filter(t == max(t, na.rm = TRUE))

ss_summary <- ss_cattle |>
  group_by(Tissue) |>
  summarise(
    n      = n(),
    Mean   = mean(Concentration),
    Median = median(Concentration),
    SD     = sd(Concentration),
    P2.5   = quantile(Concentration, 0.025),
    P97.5  = quantile(Concentration, 0.975),
    .groups = "drop"
  )

write.csv(ss_cattle,  out_path("A_SteadyState_PopDistribution.csv"),
          row.names = FALSE)
write.csv(ss_summary, out_path("A_SteadyState_PopSummary.csv"),
          row.names = FALSE)

cat("    Steady-state summary:\n")
print(ss_summary, digits = 4)

set.seed(42)   # reproducible jitter positions
p_dist <- ggplot(ss_cattle, aes(x = Tissue, y = Concentration, fill = Tissue)) +
  geom_violin(trim = FALSE, alpha = 0.35, color = NA) +
  geom_jitter(width = 0.08, alpha = 0.15, size = 0.6,
              color = "gray25", shape = 16) +
  geom_boxplot(width = 0.18, fill = "white", color = "black", alpha = 0.85,
               outlier.shape = NA, linewidth = 0.5) +
  scale_fill_manual(values = TISSUE_COLORS) +
  scale_y_continuous(
    name = expression(paste("Steady-state iAs concentration (",
                            mu, "g/kg)")),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  labs(x = "Cattle tissue compartments") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title      = element_text(color = "black", size = 14),
    axis.text.x     = element_text(color = "black", size = 13, face = "bold"),
    axis.text.y     = element_text(color = "black", size = 12),
    axis.line       = element_line(color = "black", linewidth = 0.5),
    plot.margin     = margin(t = 15, r = 20, b = 10, l = 10, unit = "pt")
  )

ggsave(out_path("A5_Cattle_SteadyState_Distribution.tiff"), p_dist,
       width = 8, height = 6, dpi = 300, compression = "lzw")
ggsave(out_path("A5_Cattle_SteadyState_Distribution.png"),  p_dist,
       width = 8, height = 6, dpi = 300)


## A13. Save PBPK CSVs --------------------------------------------------------

write.csv(out_df, out_path("A_pbpk_deterministic.csv"), row.names = FALSE)

pop_export <- pop_df_conv |>
  select(run, time, Day, C_iAs_muscle, C_iAs_liv, C_iAs_kid,
         C_iAs_muscle_ug, C_iAs_liv_ug, C_iAs_kid_ug,
         BTF_muscle, BTF_liv, BTF_kid, TF_muscle, TF_liv, TF_kid)
write.csv(pop_export, out_path(sprintf("A_pbpk_pop_%druns.csv", n_runs_pbpk)),
          row.names = FALSE)


# =============================================================================
# =============================================================================
#               PART B: BUILD TF_cattle FROM PBPK OUTPUT (in memory)
# =============================================================================
# =============================================================================

cat("\n[PART B] TF EXTRACTION\n")
cat("----------------------\n")

ss_pbpk <- pop_df_conv |> filter(time == max(time))

TF_cattle <- list(
  meat = list(
    mean = mean(ss_pbpk$TF_muscle),
    q025 = quantile(ss_pbpk$TF_muscle, 0.025, names = FALSE),
    q500 = quantile(ss_pbpk$TF_muscle, 0.500, names = FALSE),
    q975 = quantile(ss_pbpk$TF_muscle, 0.975, names = FALSE)
  ),
  liver = list(
    mean = mean(ss_pbpk$TF_liv),
    q025 = quantile(ss_pbpk$TF_liv, 0.025, names = FALSE),
    q500 = quantile(ss_pbpk$TF_liv, 0.500, names = FALSE),
    q975 = quantile(ss_pbpk$TF_liv, 0.975, names = FALSE)
  ),
  kidney = list(
    mean = mean(ss_pbpk$TF_kid),
    q025 = quantile(ss_pbpk$TF_kid, 0.025, names = FALSE),
    q500 = quantile(ss_pbpk$TF_kid, 0.500, names = FALSE),
    q975 = quantile(ss_pbpk$TF_kid, 0.975, names = FALSE)
  )
)
# 'others' = remaining tissues compartment of PBPK; shares partition
# coefficients with muscle -> TF = TF_meat
TF_cattle$others <- TF_cattle$meat

cat("TF (d/kg) from PBPK population output:\n")
for (tis in c("meat", "liver", "kidney", "others")) {
  tfi <- TF_cattle[[tis]]
  cat(sprintf("  %-8s: mean=%.6f, 95%% CI=[%.6f, %.6f]\n",
              tis, tfi$mean, tfi$q025, tfi$q975))
}

# Also persist a small TF_summary.csv for audit
tf_summary <- data.frame(
  tissue = c("meat", "liver", "kidney"),
  mean   = c(TF_cattle$meat$mean,  TF_cattle$liver$mean,  TF_cattle$kidney$mean),
  q025   = c(TF_cattle$meat$q025,  TF_cattle$liver$q025,  TF_cattle$kidney$q025),
  q500   = c(TF_cattle$meat$q500,  TF_cattle$liver$q500,  TF_cattle$kidney$q500),
  q975   = c(TF_cattle$meat$q975,  TF_cattle$liver$q975,  TF_cattle$kidney$q975),
  n_runs = n_runs_pbpk,
  units  = "d/kg",
  source = sprintf("Integrated pipeline run %s", format(Sys.time(), "%Y-%m-%d %H:%M"))
)
write.csv(tf_summary, out_path("B_TF_summary.csv"), row.names = FALSE)


# =============================================================================
# =============================================================================
#                   PART C: MRL RISK ASSESSMENT (Adult-only)
# =============================================================================
# =============================================================================

cat("\n[PART C] MRL RISK ASSESSMENT (Adult only)\n")
cat("-----------------------------------------\n")


## C1. Configuration ----------------------------------------------------------

config <- list(
  n_sim              = 10000,
  seed               = 123,
  acceptable_risk    = 1e-6,
  confidence_level   = 0.95,
  target_correlation = 0.5,
  bioavail_min       = 0.70,
  bioavail_max       = 0.90,
  use_tf_uncertainty = TRUE
)
set.seed(config$seed)

ED_age <- list("Adult" = 70)
LT     <- 70

cat("C1. Config:\n")
cat("    Trials             :", config$n_sim, "\n")
cat("    Acceptable risk    :", config$acceptable_risk, "(1 in 1,000,000)\n")
cat("    BW-FI correlation  :", config$target_correlation, "\n")
cat("    Bioavailability    : Uniform(", config$bioavail_min, ",",
    config$bioavail_max, ")\n")
cat("    ED_age = Adult     :", ED_age$Adult, "years\n")
cat("    Lifetime LT        :", LT, "years\n")
cat("    Weight             :", ED_age$Adult / LT, "(= 1 for Adult-only)\n")


## C2. Load Excel data and filter to Adult only -------------------------------

cat("C2. Loading Excel data...\n")
stopifnot("Excel file not found" = file.exists(EXCEL_PATH))

food_intake_data <- read_excel(EXCEL_PATH, sheet = "food intake")
bodyweight_data  <- read_excel(EXCEL_PATH, sheet = "bodyweight")
csf_data         <- read_excel(EXCEL_PATH, sheet = "csf")

# Adult-only filter (silences ED_age warnings later)
bodyweight_data  <- bodyweight_data  |> filter(agegroup == "Adult")
food_intake_data <- food_intake_data |> filter(agegroup == "Adult")

stopifnot("No Adult rows in bodyweight sheet"  = nrow(bodyweight_data)  > 0,
          "No Adult rows in food_intake sheet" = nrow(food_intake_data) > 0)

cat("    Food intake (Adult): ", nrow(food_intake_data), "records\n")
cat("    Bodyweight  (Adult): ", nrow(bodyweight_data),  "records\n")
cat("    CSF              : ", nrow(csf_data),           "records\n")


## C3. Cattle data preparation ------------------------------------------------

tissues     <- c("meat", "liver", "kidney", "others")
age_groups  <- "Adult"
populations <- c("general public", "consumer only")

cattle_data <- food_intake_data |>
  filter(food == "cattle") |>
  mutate(
    mean = as.numeric(mean),
    sd   = as.numeric(sd),
    min  = as.numeric(min),
    max  = as.numeric(max),
    sd   = ifelse(sd == 0 & mean > 0, 0.1 * mean, sd),
    min  = ifelse(is.na(min), pmax(0, mean - 3 * sd), min),
    max  = ifelse(is.na(max), mean + 3 * sd, max)
  ) |>
  filter(!is.na(mean), !is.na(sd), sd > 0, mean >= 0, max > min)
cat("    Cattle records   : ", nrow(cattle_data), "\n")


## C4. Other parameters --------------------------------------------------------

feeding_rate_params <- list(mean = 9.42, sd = 0.94)   # Waegeneers et al. 2011
csf_skin    <- csf_data$csf[csf_data$cancer == "skin"]
csf_lung    <- csf_data$csf[csf_data$cancer == "lung"]
csf_bladder <- csf_data$csf[csf_data$cancer == "bladder"]

cat("    CSF [(ug/kg-day)^-1]: skin=", csf_skin,
    "lung=", csf_lung, "bladder=", csf_bladder, "\n")


## C5. Sampling helpers --------------------------------------------------------

sample_tf <- function(tissue_name) {
  tf <- TF_cattle[[tissue_name]]
  est_sd <- (tf$q975 - tf$q025) / 3.92   # 95% CI width = 3.92 * SD
  rtruncnorm(1, a = tf$q025, b = tf$q975, mean = tf$mean, sd = est_sd)
}

sample_bw <- function(bw_mean, bw_sd, bw_min, bw_max) {
  if (bw_sd <= 0) return(bw_mean)
  rtruncnorm(1, a = bw_min, b = bw_max, mean = bw_mean, sd = bw_sd)
}

sample_fi_conditional <- function(bw_sampled, bw_mean, bw_sd,
                                  fi_mean, fi_sd, fi_min, fi_max,
                                  target_corr = 0.5) {
  if (bw_sd <= 0 || fi_sd <= 0) return(fi_mean)
  z_bw    <- (bw_sampled - bw_mean) / bw_sd
  epsilon <- rnorm(1)
  z_fi    <- target_corr * z_bw + sqrt(1 - target_corr^2) * epsilon
  fi_raw  <- fi_mean + z_fi * fi_sd
  max(fi_min, min(fi_max, fi_raw))
}


## C6. MC simulation engine ----------------------------------------------------

run_simulation <- function(population_type, n_trials) {

  cat("    Simulating:", population_type, "...\n")

  pop_data <- cattle_data |> filter(population == population_type)
  if (nrow(pop_data) == 0) stop("No data for: ", population_type)

  MRL_skin_vec    <- numeric(n_trials)
  MRL_lung_vec    <- numeric(n_trials)
  MRL_bladder_vec <- numeric(n_trials)
  MRL_unified_vec <- numeric(n_trials)
  total_ef_vec    <- numeric(n_trials)
  bioavail_vec    <- numeric(n_trials)
  feeding_rate_vec <- numeric(n_trials)

  # Diagnostic records for post-simulation validation (one row per trial x tissue)
  n_rec      <- n_trials * length(tissues)
  rec_trial  <- integer(n_rec)
  rec_tissue <- character(n_rec)
  rec_bw     <- numeric(n_rec)
  rec_fi     <- numeric(n_rec)
  rec_tf     <- numeric(n_rec)
  rec_i      <- 0L

  # Pre-extract BW and FI parameters
  bw_params <- list()
  for (ag in age_groups) {
    bw_row <- bodyweight_data[bodyweight_data$agegroup == ag, ]
    if (nrow(bw_row) == 0) next
    bw_params[[ag]] <- list(
      mean = as.numeric(bw_row$mean[1]), sd = as.numeric(bw_row$sd[1]),
      min  = as.numeric(bw_row$min[1]),  max = as.numeric(bw_row$max[1])
    )
  }
  fi_params <- list()
  for (ag in age_groups) {
    fi_params[[ag]] <- list()
    for (tis in tissues) {
      fi_row <- pop_data |> filter(agegroup == ag, tissue == !!tis)
      if (nrow(fi_row) == 0) next
      fi_m <- as.numeric(fi_row$mean[1]); fi_s <- as.numeric(fi_row$sd[1])
      if (is.na(fi_m) || is.na(fi_s)) next
      fi_params[[ag]][[tis]] <- list(
        mean = fi_m, sd = fi_s,
        min  = as.numeric(fi_row$min[1]), max = as.numeric(fi_row$max[1])
      )
    }
  }

  for (trial in 1:n_trials) {

    feeding_rate    <- rtruncnorm(1, a = 0,
                                  mean = feeding_rate_params$mean,
                                  sd   = feeding_rate_params$sd)
    bioavailability <- runif(1, config$bioavail_min, config$bioavail_max)

    bw_sampled <- list()
    for (ag in age_groups) {
      bp <- bw_params[[ag]]; if (is.null(bp)) next
      bw_sampled[[ag]] <- sample_bw(bp$mean, bp$sd, bp$min, bp$max)
    }

    total_ef <- 0
    for (tissue in tissues) {

      tf_val <- if (config$use_tf_uncertainty) sample_tf(tissue) else TF_cattle[[tissue]]$mean

      tissue_age_sum <- 0
      fi_recorded    <- NA_real_
      for (age_group in age_groups) {
        ed_val <- ED_age[[age_group]]
        if (is.null(ed_val) || ed_val == 0) next
        bw_val <- bw_sampled[[age_group]]; if (is.null(bw_val)) next
        fp     <- fi_params[[age_group]][[tissue]]; if (is.null(fp)) next
        bp     <- bw_params[[age_group]]

        fi_val <- sample_fi_conditional(
          bw_sampled = bw_val, bw_mean = bp$mean, bw_sd = bp$sd,
          fi_mean = fp$mean, fi_sd = fp$sd, fi_min = fp$min, fi_max = fp$max,
          target_corr = config$target_correlation
        )
        # Age-weighted: (IR_age / 1000) * ED_age / (BW_age * LT)
        tissue_age_sum <- tissue_age_sum + (fi_val / 1000) * ed_val / (bw_val * LT)
        fi_recorded    <- fi_val   # capture for validation (Adult-only)
      }

      # Diagnostic record, one row per trial x tissue. TF and BW are always
      # sampled; FI is NA for tissues with no consumption data (liver and
      # kidney have zero reported cattle intake in the Taiwan dataset).
      rec_i <- rec_i + 1L
      rec_trial[rec_i]  <- trial
      rec_tissue[rec_i] <- tissue
      rec_bw[rec_i]     <- bw_sampled[[age_groups[1]]]
      rec_fi[rec_i]     <- fi_recorded
      rec_tf[rec_i]     <- tf_val

      total_ef <- total_ef + tf_val * feeding_rate * bioavailability * tissue_age_sum
    }

    if (total_ef > 0) {
      MRL_skin_vec[trial]    <- config$acceptable_risk / (csf_skin    * total_ef)
      MRL_lung_vec[trial]    <- config$acceptable_risk / (csf_lung    * total_ef)
      MRL_bladder_vec[trial] <- config$acceptable_risk / (csf_bladder * total_ef)
      MRL_unified_vec[trial] <- min(MRL_skin_vec[trial], MRL_lung_vec[trial], MRL_bladder_vec[trial])
    } else {
      MRL_skin_vec[trial] <- NA; MRL_lung_vec[trial] <- NA
      MRL_bladder_vec[trial] <- NA; MRL_unified_vec[trial] <- NA
    }
    total_ef_vec[trial]     <- total_ef
    bioavail_vec[trial]     <- bioavailability
    feeding_rate_vec[trial] <- feeding_rate
  }

  results <- data.frame(
    trial        = 1:n_trials,
    MRL_skin     = MRL_skin_vec,
    MRL_lung     = MRL_lung_vec,
    MRL_bladder  = MRL_bladder_vec,
    MRL_unified  = MRL_unified_vec,
    total_ef     = total_ef_vec,
    bioavail     = bioavail_vec,
    feeding_rate = feeding_rate_vec
  )
  # Attach trial x tissue diagnostic records (BW, FI, TF) for validation
  attr(results, "diag") <- data.frame(
    trial  = rec_trial[seq_len(rec_i)],
    tissue = rec_tissue[seq_len(rec_i)],
    bw     = rec_bw[seq_len(rec_i)],
    fi     = rec_fi[seq_len(rec_i)],
    tf     = rec_tf[seq_len(rec_i)]
  )
  results
}


## C7. Run MC simulations -----------------------------------------------------

cat("C3. Running Monte Carlo simulations (", config$n_sim, "trials each)...\n")
t_start_mrl <- proc.time()
results_general  <- run_simulation("general public", config$n_sim)
results_consumer <- run_simulation("consumer only",  config$n_sim)
elapsed_mrl <- (proc.time() - t_start_mrl)["elapsed"]
cat(sprintf("    Done in %.1f s\n", elapsed_mrl))


## C8. Statistical analysis ---------------------------------------------------

cat("C4. Statistical analysis...\n")

calc_summary <- function(results, pop_name) {
  MRL <- results$MRL_unified[!is.na(results$MRL_unified)]
  if (length(MRL) == 0) return(NULL)
  data.frame(
    Population = pop_name, N_Valid = length(MRL),
    Mean = mean(MRL), Median = median(MRL),
    SD = sd(MRL), CV = sd(MRL) / mean(MRL),
    Min = min(MRL), Max = max(MRL),
    P01 = quantile(MRL, 0.01), P05 = quantile(MRL, 0.05),
    P10 = quantile(MRL, 0.10), P50 = quantile(MRL, 0.50),
    P95 = quantile(MRL, 0.95), P99 = quantile(MRL, 0.99)
  )
}

summary_general  <- calc_summary(results_general,  "General Public")
summary_consumer <- calc_summary(results_consumer, "Consumer Only")
combined_summary <- rbind(summary_general, summary_consumer)

cat("\n  Summary (MRL_unified, ug/kg):\n")
print(combined_summary[, c("Population", "N_Valid", "Mean", "Median", "P05", "P95", "CV")],
      digits = 5, row.names = FALSE)


## C9. Regulatory recommendation ----------------------------------------------

p5_general  <- summary_general$P05
p5_consumer <- summary_consumer$P05
final_MRL   <- min(p5_general, p5_consumer)
protecting  <- if (p5_general < p5_consumer) "General Public" else "Consumer Only"

cat("\n  Regulatory recommendation:\n")
cat(sprintf("    P5 General Public: %8.4f ug/kg\n", p5_general))
cat(sprintf("    P5 Consumer Only : %8.4f ug/kg\n", p5_consumer))
cat(sprintf("    FINAL MRL        : %8.4f ug/kg  (protecting %s)\n",
            final_MRL, protecting))


## C10. Sensitivity analysis ---------------------------------------------------

cat("C5. Sensitivity analysis (general public)...\n")
sens_data <- results_general |>
  select(MRL_unified, total_ef, bioavail) |>
  filter(!is.na(MRL_unified))

corr_ef  <- cor(sens_data$MRL_unified, sens_data$total_ef, method = "spearman")
corr_bio <- cor(sens_data$MRL_unified, sens_data$bioavail, method = "spearman")
cat(sprintf("    Spearman r(MRL, total_ef)  = %.4f\n", corr_ef))
cat(sprintf("    Spearman r(MRL, bioavail) = %.4f\n", corr_bio))

sensitivity_results <- data.frame(
  Parameter   = c("Total Exposure Factor", "Bioavailability"),
  Correlation = c(corr_ef, corr_bio)
) |> mutate(Abs_Corr = abs(Correlation)) |> arrange(desc(Abs_Corr))


## C11. Post-simulation validation --------------------------------------------
# Verifies the mathematical / structural integrity of the Monte Carlo run.
# Two categories:
#   (1) Pearson r(BW, FI) per tissue vs the target correlation (0.50)
#   (2) Distributional integrity: a one-sample Kolmogorov-Smirnov test
#       comparing each Monte Carlo-sampled input to the distribution it was
#       DRAWN from (its intended target) - truncated-normal for TF and feeding
#       rate, uniform for bioavailability. This confirms the sampler reproduced
#       the intended distributions, rather than testing against plain normality
#       (TF is a truncated normal by construction).
# Performed on the general-public run.

cat("C6. Post-simulation validation (general public)...\n")
diag_gp <- attr(results_general, "diag")

# (1) BW-FI Pearson correlation per tissue ------------------------------------
# Liver and kidney have zero reported cattle intake (FI = NA) -> no correlation.
val_corr <- do.call(rbind, lapply(tissues, function(tis) {
  s  <- diag_gp[diag_gp$tissue == tis, ]
  ok <- is.finite(s$bw) & is.finite(s$fi)
  data.frame(tissue = tis, n = sum(ok),
             r_pearson = if (sum(ok) >= 3)
               cor(s$bw[ok], s$fi[ok], method = "pearson") else NA_real_,
             stringsAsFactors = FALSE)
}))
cat(sprintf("    (1) BW-FI Pearson correlation (target r = %.2f):\n",
            config$target_correlation))
for (i in seq_len(nrow(val_corr))) {
  if (is.na(val_corr$r_pearson[i]))
    cat(sprintf("        %-8s no consumption data (FI = 0)\n",
                val_corr$tissue[i]))
  else
    cat(sprintf("        %-8s r = %.3f  (n = %d)\n",
                val_corr$tissue[i], val_corr$r_pearson[i], val_corr$n[i]))
}

# (2) Distributional integrity vs intended target ----------------------------
# One-sample KS test against the exact distribution each input was drawn from.
# The KS D statistic (max CDF deviation) is the primary metric; at n = 10000
# the p-value is sensitive to ordinary sampling noise and is secondary.
ks_target <- function(x, cdf, ...) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(data.frame(D = NA_real_, p_value = NA_real_))
  k <- suppressWarnings(ks.test(x, cdf, ...))
  data.frame(D = unname(k$statistic), p_value = k$p.value)
}

# (2a) TF vs intended truncated-normal target (a = q025, b = q975)
val_tf <- do.call(rbind, lapply(tissues, function(tis) {
  tf  <- TF_cattle[[tis]]
  esd <- (tf$q975 - tf$q025) / 3.92          # same SD that sample_tf() uses
  k   <- ks_target(diag_gp$tf[diag_gp$tissue == tis], "ptruncnorm",
                   a = tf$q025, b = tf$q975, mean = tf$mean, sd = esd)
  data.frame(variable = paste0("TF_", tis), target = "truncated normal",
             D = k$D, p_value = k$p_value, stringsAsFactors = FALSE)
}))

# (2b) Feeding rate vs intended truncated-normal target (a = 0)
k_fr <- ks_target(results_general$feeding_rate, "ptruncnorm",
                  a = 0, b = Inf, mean = feeding_rate_params$mean,
                  sd = feeding_rate_params$sd)
val_fr <- data.frame(variable = "Feeding_rate", target = "truncated normal",
                     D = k_fr$D, p_value = k_fr$p_value, stringsAsFactors = FALSE)

# (2c) Bioavailability vs intended uniform target
k_bio <- ks_target(results_general$bioavail, "punif",
                   min = config$bioavail_min, max = config$bioavail_max)
val_bio <- data.frame(variable = "Bioavailability", target = "uniform",
                      D = k_bio$D, p_value = k_bio$p_value, stringsAsFactors = FALSE)

val_dist <- rbind(val_tf, val_fr, val_bio)
cat("    (2) Distributional integrity - KS test vs intended target (small D = match):\n")
for (i in seq_len(nrow(val_dist)))
  cat(sprintf("        %-16s vs %-18s D = %.4f, p = %.4f\n",
              val_dist$variable[i], val_dist$target[i],
              val_dist$D[i], val_dist$p_value[i]))
cat(sprintf("        -> max D = %.4f across all sampled inputs\n",
            max(val_dist$D, na.rm = TRUE)))

# Feeding-rate target check (mean / SD)
cat(sprintf("        Feeding rate: mean = %.3f (target %.2f), SD = %.3f (target %.2f)\n",
            mean(results_general$feeding_rate), feeding_rate_params$mean,
            sd(results_general$feeding_rate),   feeding_rate_params$sd))

# Export validation summary --------------------------------------------------
val_out <- rbind(
  data.frame(Category = "BW-FI correlation", Variable = val_corr$tissue,
             Metric = "Pearson r", Value = round(val_corr$r_pearson, 4),
             P_value = NA_real_,
             Reference = sprintf("target r = %.2f", config$target_correlation),
             stringsAsFactors = FALSE),
  data.frame(Category = "Distributional integrity", Variable = val_dist$variable,
             Metric = paste0("KS D vs ", val_dist$target),
             Value = round(val_dist$D, 4),
             P_value = round(val_dist$p_value, 4),
             Reference = "small D = matches intended target",
             stringsAsFactors = FALSE)
)
write.csv(val_out, out_path("C_validation_summary.csv"), row.names = FALSE)

# Validation figure: BW-FI scatter (only tissues with consumption data) -------
diag_fi    <- diag_gp[is.finite(diag_gp$fi), ]
fi_tissues <- intersect(tissues, unique(diag_fi$tissue))
diag_fi$tissue <- factor(diag_fi$tissue, levels = fi_tissues)
lab_df <- data.frame(
  tissue = factor(val_corr$tissue[!is.na(val_corr$r_pearson)], levels = fi_tissues),
  label  = sprintf("r = %.3f", val_corr$r_pearson[!is.na(val_corr$r_pearson)])
)
p_val <- ggplot(diag_fi, aes(bw, fi)) +
  geom_point(alpha = 0.08, size = 0.4, colour = "#4393c3") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "firebrick", linewidth = 0.8) +
  geom_text(data = lab_df, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.15, vjust = 1.6, size = 3.8, fontface = "bold",
            inherit.aes = FALSE) +
  facet_wrap(~tissue, scales = "free_y") +
  labs(title = sprintf("BW-FI correlation by tissue (target r = %.2f)",
                       config$target_correlation),
       subtitle = "General public; red line = linear fit. Liver and kidney omitted (zero reported intake).",
       x = "Body weight (kg)", y = "Food intake (g/day)") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))
ggsave(out_path("C5_BW_FI_correlation.png"), p_val,
       width = 9, height = 4.5, dpi = 300)


## C12. MRL plots --------------------------------------------------------------

cat("C7. Creating MRL plots...\n")

plot_data <- rbind(
  data.frame(Population = "General Public", MRL_unified = results_general$MRL_unified),
  data.frame(Population = "Consumer Only",  MRL_unified = results_consumer$MRL_unified)
) |> filter(!is.na(MRL_unified))

# P1: Distribution histogram
pC1 <- ggplot(plot_data, aes(x = MRL_unified, fill = Population)) +
  geom_histogram(alpha = 0.7, bins = 50, position = "identity",
                 color = "black", linewidth = 0.2) +
  # extra right-side expansion + plot margin so the last (largest) axis
  # label is not clipped at the image edge
  scale_x_log10(labels = scales::comma_format(),
                expand = expansion(mult = c(0.03, 0.12))) +
  geom_vline(xintercept = final_MRL, color = "red",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = final_MRL * 1.2, y = Inf,
           label = paste("Proposed MRL:", round(final_MRL, 2), "ug/kg"),
           vjust = 1.5, hjust = 0, color = "red", size = 4, fontface = "bold") +
  labs(x = "Unified MRL (ug/kg, log-scale)", y = "Frequency", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        plot.margin = margin(t = 10, r = 28, b = 10, l = 10, unit = "pt"))
ggsave(out_path("C1_MRL_distribution.png"),
       pC1, width = 8, height = 6, dpi = 300, bg = "white")

# P2: Boxplot
pC2 <- ggplot(plot_data, aes(Population, MRL_unified, fill = Population)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  scale_y_log10(labels = scales::comma_format()) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "red") +
  geom_hline(yintercept = final_MRL, color = "red", linetype = "dashed") +
  labs(title = "MRL comparison (cattle, Adult only)",
       subtitle = "Red diamond = mean | red line = final recommendation",
       x = NULL, y = "MRL (ug/kg, log-scale)") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(out_path("C2_MRL_boxplot.png"),
       pC2, width = 8, height = 6, dpi = 300, bg = "white")

# P3: Cancer endpoint comparison (general public)
endpoint_data <- data.frame(
  Cancer_Type = rep(c("Skin", "Lung", "Bladder", "Unified"),
                    each = nrow(results_general)),
  MRL_Value   = c(results_general$MRL_skin, results_general$MRL_lung,
                  results_general$MRL_bladder, results_general$MRL_unified)
) |> filter(!is.na(MRL_Value))

pC3 <- ggplot(endpoint_data, aes(Cancer_Type, MRL_Value, fill = Cancer_Type)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_log10(labels = scales::comma_format()) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(title = "MRL by cancer endpoint (General Public, cattle)",
       subtitle = "Unified = min(Skin, Lung, Bladder)",
       x = "Endpoint", y = "MRL (ug/kg, log-scale)") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(out_path("C3_MRL_by_endpoint.png"),
       pC3, width = 8, height = 6, dpi = 300, bg = "white")

# P4: Sensitivity tornado
pC4 <- ggplot(sensitivity_results,
              aes(reorder(Parameter, Abs_Corr), Correlation)) +
  geom_col(aes(fill = Correlation > 0), width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "#2196F3", "FALSE" = "#F44336")) +
  coord_flip() +
  labs(title = "Sensitivity tornado",
       subtitle = "Spearman rank correlation with MRL_unified",
       x = NULL, y = "Rank correlation") +
  theme_minimal(base_size = 12)
ggsave(out_path("C4_sensitivity_tornado.png"),
       pC4, width = 8, height = 4, dpi = 300, bg = "white")


## C13. Export CSV tables ------------------------------------------------------

cat("C8. Exporting tables and summary...\n")

write.csv(results_general,  out_path("C_results_general_public.csv"), row.names = FALSE)
write.csv(results_consumer, out_path("C_results_consumer_only.csv"),  row.names = FALSE)
write.csv(combined_summary, out_path("C_summary_statistics.csv"),     row.names = FALSE)

regulatory_summary <- data.frame(
  Parameter = c("General_Public_P5", "Consumer_Only_P5", "Final_MRL",
                "Acceptable_Risk", "Lifetime_LT", "Simulation_Trials",
                "TF_meat_mean", "TF_liver_mean", "TF_kidney_mean", "PBPK_n_runs"),
  Value     = c(round(p5_general, 6), round(p5_consumer, 6), round(final_MRL, 6),
                config$acceptable_risk, LT, config$n_sim,
                TF_cattle$meat$mean, TF_cattle$liver$mean, TF_cattle$kidney$mean,
                n_runs_pbpk),
  Units     = c("ug/kg", "ug/kg", "ug/kg", "unitless", "years", "trials",
                "d/kg", "d/kg", "d/kg", "runs")
)
write.csv(regulatory_summary, out_path("C_regulatory_summary.csv"), row.names = FALSE)

# Table 4 publication
calc_endpoint_stats <- function(values) {
  v <- values[!is.na(values)]
  list(mean = mean(v), p5 = quantile(v, 0.05), p95 = quantile(v, 0.95))
}
gen_skin    <- calc_endpoint_stats(results_general$MRL_skin)
gen_lung    <- calc_endpoint_stats(results_general$MRL_lung)
gen_bladder <- calc_endpoint_stats(results_general$MRL_bladder)
con_skin    <- calc_endpoint_stats(results_consumer$MRL_skin)
con_lung    <- calc_endpoint_stats(results_consumer$MRL_lung)
con_bladder <- calc_endpoint_stats(results_consumer$MRL_bladder)

table4 <- data.frame(
  Cancer_endpoint = rep(c("Skin", "Lung", "Bladder"), 2),
  Population      = c(rep("General public", 3), rep("Consumer only", 3)),
  Mean = c(gen_skin$mean, gen_lung$mean, gen_bladder$mean,
           con_skin$mean, con_lung$mean, con_bladder$mean),
  P5   = c(gen_skin$p5, gen_lung$p5, gen_bladder$p5,
           con_skin$p5, con_lung$p5, con_bladder$p5),
  P95  = c(gen_skin$p95, gen_lung$p95, gen_bladder$p95,
           con_skin$p95, con_lung$p95, con_bladder$p95)
)
write.csv(table4, out_path("C_Table4_MRL_by_endpoint.csv"), row.names = FALSE)


## C14. Graphical abstract ----------------------------------------------------
# Five-panel summary banner: Feed -> Population PBPK -> Transfer Factors ->
# Monte Carlo -> Proposed MRL. All numbers are pulled live from the pipeline
# (n_runs_pbpk, TF_cattle, config$n_sim, final_MRL) so the figure always
# matches the canonical run.

cat("C9. Building graphical abstract...\n")

# --- palette ---
ga_bg <- "#F7F9FC"; ga_feed_c <- "#E8F0FE"; ga_pbpk_c <- "#D4E8D1"
ga_tf_c <- "#FFF3E0"; ga_mc_c <- "#F3E5F5"; ga_mrl_c <- "#FFEBEE"
ga_accent <- "#1565C0"; ga_green <- "#2E7D32"; ga_orange <- "#E65100"
ga_purple <- "#6A1B9A"; ga_red <- "#C62828"; ga_border <- "#546E7A"
ga_text <- "#212121"; ga_sub <- "#666666"

# --- helpers ---
ga_canvas <- function(bg = "white") {
  ggplot() +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 10), expand = FALSE) +
    theme_void() +
    theme(plot.background  = element_rect(fill = bg, color = NA),
          panel.background = element_rect(fill = bg, color = NA),
          plot.margin = margin(0, 0, 0, 0))
}
ga_ellipse <- function(cx, cy, a, b, n = 100) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + a * cos(theta), y = cy + b * sin(theta))
}
ga_arrow <- function() {
  ga_canvas(ga_bg) +
    annotate("segment", x = 2, xend = 8, y = 5, yend = 5,
             arrow = arrow(length = unit(0.15, "inches"), type = "closed"),
             color = ga_border, linewidth = 1.2)
}

# --- live numbers ---
ga_n_pbpk    <- format(n_runs_pbpk, big.mark = ",")
ga_n_mc      <- format(config$n_sim, big.mark = ",")
ga_feed_rate <- unname(parms["feed_intake_kg"])
ga_mrl      <- round(final_MRL, 1)
# Current international tAs feed MRLs span 2,000-30,000 ug/kg
# (CFIA 2015; EU 2019; FSANZ 2001; U.S. NRC 2005)
ga_fold_lo  <- round(2000  / ga_mrl)
ga_fold_hi  <- round(30000 / ga_mrl)

# --- Panel 1: Feed input ---
p_feed <- ga_canvas(ga_feed_c) +
  annotate("rect", xmin = 0.2, xmax = 9.8, ymin = 0.2, ymax = 9.8,
           fill = NA, color = ga_border, linewidth = 0.6) +
  annotate("text", x = 5, y = 9.0, label = "FEED INPUT",
           size = 3.2, fontface = "bold", color = ga_accent) +
  annotate("rect", xmin = 2.8, xmax = 7.2, ymin = 3.0, ymax = 6.3,
           fill = "#8D6E63", color = "#6D4C41", linewidth = 0.5) +
  annotate("polygon", x = c(3.1, 6.9, 6.5, 3.5), y = c(6.3, 6.3, 7.0, 7.0),
           fill = "#A1887F", color = "#6D4C41", linewidth = 0.5) +
  annotate("text", x = 5, y = 5.3, label = "iAs", size = 4.5,
           fontface = "bold", color = "white") +
  annotate("text", x = 5, y = 4.0, label = "FEED", size = 2.8,
           fontface = "bold", color = "#D7CCC8") +
  annotate("text", x = 5, y = 1.5,
           label = sprintf("Feeding rate = %.2f kg DM/day", ga_feed_rate),
           size = 2.4, color = ga_sub)

# --- Panel 2: Population PBPK ---
ga_body   <- ga_ellipse(5.0, 5.8, 3.4, 2.0)
ga_head   <- ga_ellipse(8.5, 6.8, 1.1, 0.9)
ga_liver  <- ga_ellipse(4.2, 6.2, 1.0, 0.60)
ga_kidney <- ga_ellipse(5.5, 5.8, 0.8, 0.50)
ga_offal  <- ga_ellipse(6.7, 5.2, 0.8, 0.55)
ga_muscle <- ga_ellipse(3.3, 4.8, 0.9, 0.5)

p_pbpk <- ga_canvas(ga_pbpk_c) +
  annotate("rect", xmin = 0.1, xmax = 9.9, ymin = 0.1, ymax = 9.9,
           fill = NA, color = ga_border, linewidth = 0.6) +
  annotate("text", x = 5, y = 9.2, label = "POPULATION PBPK",
           size = 3.2, fontface = "bold", color = ga_green) +
  annotate("text", x = 5, y = 8.5,
           label = sprintf("n = %s virtual cattle", ga_n_pbpk),
           size = 2.5, color = ga_sub) +
  geom_polygon(data = ga_body,   aes(x, y), fill = "#E0E0E0", color = "#888888", linewidth = 0.5) +
  geom_polygon(data = ga_head,   aes(x, y), fill = "#EEEEEE", color = "#888888", linewidth = 0.5) +
  annotate("segment", x = c(2.8, 3.8, 6.2, 7.2), xend = c(2.8, 3.8, 6.2, 7.2),
           y = rep(3.8, 4), yend = rep(2.8, 4), color = "#888888", linewidth = 0.8) +
  geom_polygon(data = ga_liver,  aes(x, y), fill = "#A5D6A7", color = "#2E7D32", linewidth = 0.5) +
  geom_polygon(data = ga_kidney, aes(x, y), fill = "#FFCC80", color = "#E65100", linewidth = 0.5) +
  geom_polygon(data = ga_offal,  aes(x, y), fill = "#BCAAA4", color = "#5D4037", linewidth = 0.5) +
  geom_polygon(data = ga_muscle, aes(x, y), fill = "#EF9A9A", color = "#C62828", linewidth = 0.5) +
  annotate("text", x = 4.2, y = 6.2, label = "Liver",  size = 1.8, fontface = "bold", color = "#1B5E20") +
  annotate("text", x = 5.5, y = 5.8, label = "Kidney", size = 1.6, fontface = "bold", color = "#BF360C") +
  annotate("text", x = 6.7, y = 5.2, label = "Other\noffal", size = 1.5, fontface = "bold",
           color = "#3E2723", lineheight = 0.8) +
  annotate("text", x = 3.3, y = 4.8, label = "Muscle", size = 1.6, fontface = "bold", color = "#B71C1C") +
  annotate("text", x = 5, y = 2.2,
           label = expression(paste("Allometric scaling (BW"^"0.75", ")")),
           size = 2.5, color = ga_text) +
  annotate("text", x = 5, y = 1.5, label = "Michaelis-Menten methylation",
           size = 2.5, color = ga_text) +
  annotate("text", x = 5, y = 0.8, label = "Steady state: ~200 hours",
           size = 2.5, color = ga_text)

# --- Panel 3: Transfer factors (live TF_cattle values) ---
ga_tf <- data.frame(
  tissue = factor(c("Muscle", "Other\noffal", "Kidney", "Liver"),
                  levels = c("Muscle", "Other\noffal", "Kidney", "Liver")),
  TF = c(TF_cattle$meat$mean, TF_cattle$others$mean,
         TF_cattle$kidney$mean, TF_cattle$liver$mean)
)
ga_tf_max <- max(ga_tf$TF)

p_tf <- ggplot(ga_tf, aes(x = tissue, y = TF)) +
  geom_col(aes(fill = tissue, color = tissue), width = 0.7, linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.4f", TF), color = tissue),
            angle = 90, hjust = -0.2, vjust = 0.5, size = 2.2,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("Muscle" = "#EF9A9A", "Other\noffal" = "#BCAAA4",
                               "Kidney" = "#FFCC80", "Liver" = "#A5D6A7")) +
  scale_color_manual(values = c("Muscle" = "#C62828", "Other\noffal" = "#5D4037",
                                "Kidney" = "#E65100", "Liver" = "#2E7D32")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.45)), limits = c(0, NA)) +
  labs(title = "TRANSFER FACTORS",
       subtitle = expression(paste("TF"["iAs/iAs"], " (d/kg)")),
       caption = "Liver > Kidney > Offal ≈ Muscle\nFirst iAs-specific TFs") +
  theme_void() +
  theme(
    plot.background  = element_rect(fill = ga_tf_c, color = ga_border, linewidth = 0.6),
    panel.background = element_rect(fill = ga_tf_c, color = NA),
    plot.title    = element_text(hjust = 0.5, size = 9.0, face = "bold",
                                 color = ga_orange, margin = margin(t = 6, b = 0)),
    plot.subtitle = element_text(hjust = 0.5, size = 7.5, color = ga_sub,
                                 margin = margin(t = 2, b = 6)),
    plot.caption  = element_text(hjust = 0.5, size = 6.0, color = ga_text,
                                 lineheight = 1.3, margin = margin(t = 5, b = 2)),
    axis.text.x   = element_text(size = 6.5, color = ga_text,
                                 margin = margin(t = 3), lineheight = 0.8),
    axis.line.x   = element_line(color = "#AAAAAA", linewidth = 0.4),
    legend.position = "none",
    plot.margin   = margin(3, 8, 3, 8)
  ) +
  coord_cartesian(clip = "off")

# --- Panel 4: Monte Carlo ---
ga_mc_x  <- seq(0, 10, length.out = 200)
ga_mc_df <- data.frame(x = ga_mc_x, y = dnorm(ga_mc_x, mean = 5.5, sd = 1.5))
ga_p5_x  <- qnorm(0.05, mean = 5.5, sd = 1.5)

p_mc <- ga_canvas(ga_mc_c) +
  annotate("rect", xmin = 0.1, xmax = 9.9, ymin = 0.1, ymax = 9.9,
           fill = NA, color = ga_border, linewidth = 0.6) +
  annotate("text", x = 5, y = 9.2, label = "MONTE CARLO",
           size = 3.2, fontface = "bold", color = ga_purple) +
  annotate("text", x = 5, y = 8.5, label = sprintf("%s iterations", ga_n_mc),
           size = 2.5, color = ga_sub) +
  geom_ribbon(data = ga_mc_df, aes(x = x, ymin = 3.5, ymax = 3.5 + y * 15),
              fill = "#CE93D8", alpha = 0.5) +
  geom_line(data = ga_mc_df, aes(x = x, y = 3.5 + y * 15),
            color = ga_purple, linewidth = 0.7) +
  annotate("segment", x = ga_p5_x, xend = ga_p5_x, y = 3.2, yend = 7.5,
           color = ga_red, linewidth = 0.7, linetype = "dashed") +
  annotate("text", x = ga_p5_x, y = 7.9, label = expression(P[5]),
           size = 3, fontface = "bold", color = ga_red) +
  annotate("text", x = 5, y = 2.8, label = "The consumer-only",
           size = 2.5, color = ga_text) +
  annotate("text", x = 5, y = 2.1, label = "Cattle products consumption",
           size = 2.5, color = ga_text) +
  annotate("text", x = 5, y = 1.4, label = "Skin cancer endpoint",
           size = 2.5, color = ga_text) +
  annotate("text", x = 5, y = 0.6, label = expression(paste("Target: CR = ", 10^-6)),
           size = 2.6, fontface = "bold", color = ga_red)

# --- Panel 5: Proposed MRL ---
p_mrl <- ga_canvas(ga_mrl_c) +
  annotate("rect", xmin = 0.2, xmax = 9.8, ymin = 0.2, ymax = 9.8,
           fill = NA, color = ga_red, linewidth = 1.2) +
  annotate("text", x = 5, y = 9.0, label = "PROPOSED MRL",
           size = 3.2, fontface = "bold", color = ga_red) +
  annotate("text", x = 5, y = 7.0, label = sprintf("%.1f", ga_mrl),
           size = 13, fontface = "bold", color = ga_red) +
  annotate("text", x = 5, y = 5.6, label = "µg/kg",
           size = 5, fontface = "bold", color = ga_red) +
  annotate("text", x = 5, y = 4.8, label = "iAs in cattle feed",
           size = 2.5, color = ga_sub) +
  annotate("segment", x = 1.5, xend = 8.5, y = 4.0, yend = 4.0,
           color = "#CCCCCC", linewidth = 0.4) +
  annotate("text", x = 5, y = 3.3, label = "vs. Current tAs feed MRLs",
           size = 2.6, fontface = "bold", color = ga_text) +
  annotate("text", x = 5, y = 2.2, label = "2,000–30,000 µg/kg",
           size = 2.5, color = ga_sub) +
  annotate("text", x = 5, y = 1.2,
           label = sprintf("%s–%s× higher",
                           format(ga_fold_lo, big.mark = ","),
                           format(ga_fold_hi, big.mark = ",")),
           size = 2.6, fontface = "bold", color = ga_red)

# --- Title / footer banners ---
p_ga_title <- ggplot() +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 3), expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = ga_accent, color = NA)) +
  annotate("text", x = 5, y = 2.0,
           label = "Health-protective feed MRLs for inorganic arsenic in cattle",
           size = 5.2, fontface = "bold", color = "white") +
  annotate("text", x = 5, y = 0.8,
           label = "A population PBPK-probabilistic risk assessment framework",
           size = 4.2, fontface = "italic", color = "#BBDEFB")

p_ga_foot <- ggplot() +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 3), expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = "#ECEFF1", color = NA)) +
  annotate("text", x = 5, y = 2.0,
           label = "Feed-to-tissue iAs transfer pathway -> Probabilistic cancer risk assessment -> Health-based feed MRL",
           size = 3.2, fontface = "italic", color = ga_border) +
  annotate("text", x = 5, y = 0.8,
           label = "The framework for deriving iAs-specific MRLs in cattle feed",
           size = 3, fontface = "bold", color = ga_accent)

# --- assemble ---
ga_middle <- arrangeGrob(
  ggplotGrob(p_feed), ggplotGrob(ga_arrow()), ggplotGrob(p_pbpk),
  ggplotGrob(ga_arrow()), ggplotGrob(p_tf), ggplotGrob(ga_arrow()),
  ggplotGrob(p_mc), ggplotGrob(ga_arrow()), ggplotGrob(p_mrl),
  ncol = 9, widths = c(3, 1, 3, 1, 3, 1, 3, 1, 3))

ga_final <- arrangeGrob(ggplotGrob(p_ga_title), ga_middle, ggplotGrob(p_ga_foot),
                        nrow = 3, heights = c(1, 6, 1))

ggsave(out_path("C6_graphical_abstract.png"), ga_final,
       width = 2656, height = 1062, units = "px", dpi = 300, bg = "white")


## C15. Final console summary --------------------------------------------------

cat("\n", strrep("=", 65), "\n", sep = "")
cat("FINAL INTEGRATED-PIPELINE SUMMARY (Adult-only)\n")
cat(strrep("=", 65), "\n", sep = "")

cat("\nPART A - PBPK\n")
cat(sprintf("  Dose (Hung 2021)    : AsIII %.2f, AsV %.2f, MMA %.2f, DMA %.2f ug/day\n",
            parms["PdoseC_AsIII"], parms["PdoseC_AsV"],
            parms["PdoseC_MMA"], parms["PdoseC_DMA"]))
cat(sprintf("  Population runs     : %d (steady state at t = %d h)\n",
            n_runs_pbpk, unname(parms["TSTOP"])))

cat("\nPART B - TF (d/kg)\n")
for (tis in c("meat", "liver", "kidney")) {
  tfi <- TF_cattle[[tis]]
  cat(sprintf("  TF_%-7s: mean = %.6f | 95%% CI [%.6f, %.6f]\n",
              tis, tfi$mean, tfi$q025, tfi$q975))
}

cat("\nPART C - MRL\n")
cat(sprintf("  Feed: C_iAs = %.0f ug/kg DM | FR = %.2f kg/day | LT = %d y | Adult only\n",
            C_feed_iAs_ug_kg, feed_intake_kg, LT))
cat(sprintf("  General Public P5: %.4f ug/kg\n", p5_general))
cat(sprintf("  Consumer Only P5 : %.4f ug/kg\n", p5_consumer))
cat(sprintf("  FINAL MRL        : %.4f ug/kg  (protecting %s)\n",
            final_MRL, protecting))
cat(sprintf("  Risk level       : %.0e (1 in 1,000,000)\n",
            config$acceptable_risk))

cat("\nAll outputs in:", OUT_DIR, "\n")
cat(strrep("=", 65), "\n", sep = "")
