# PEACE-3 data preparation for MAIC/ML-NMR
# Input: PEACE-3 IPD data transfer 12May2025
# Output: IPD with TTE outcomes (OS, PFS) and baseline characteristics
#         Age (cont), ECOG (cat), WHO-PS (cat), PSA (cont & cat)
# Author: David Aceituno
# Date: 2025-11-14

# 1. Packages --------------------------------------------------------------
pkgs <- c("tidyverse","haven","survival","here", "gtsummary")
lapply(pkgs, library, character.only = T)

# 2. Global settings and working directory -------------------------------------
peace3_directory <- here("data", "data transfer 12May2025")
results_directory <- here("results")
source(here("R", "utils.R"))

# 3. Import data ---------------------------------------------------------------
adsl <- read_sas(paste0(peace3_directory, "/adsl.sas7bdat")) 
adtte <- read_sas(paste0(peace3_directory, "/adtte.sas7bdat")) 

# 4. Create TTE variables  --------------------------------------
death_df <- adtte |> 
  dplyr::filter(PARAMCD == "DEATH") |>
  dplyr::transmute(USUBJID,
                   time_death  = as.numeric(AVAL/(365.25)),
                   event_death = as.integer(1 - CNSR))

pfs_df <- adtte |>
  dplyr::filter(PARAMCD == "RPFS") |>
  dplyr::transmute(USUBJID,
                   time_pfs  = as.numeric(AVAL/(365.25)),
                   event_pfs = as.integer(1 - CNSR))

tte_all <- full_join(death_df, pfs_df, by = "USUBJID")


# Quick check  
nrow(tte_all) # 413 participants (check with Howard why there are 33 fewer patients than Tombal 2023)
sum(tte_all$event_death) # 219 events 
sum(tte_all$event_pfs) # 275 events 

# 5. Baseline characteristics ----------------------------------------------------
# According to the SAP: age, WHO-PS, Gleason, PSA
peace3_bsl <- merge(tte_all, adsl, by = "USUBJID", all.x = TRUE)

baseline_all <- peace3_bsl |> 
  mutate(
    treatment = TRT01P,
    who = BASEWHO,
    age = AGE,
    gleason_m8 = if_else(GLEASONN >= 8, ">=8", "<8", missing = NA),
    psa =  PSABL) |> 
  select(USUBJID, treatment, time_death, event_death, event_pfs, time_pfs, 
         age, who, gleason_m8, psa)  |> 
  distinct(USUBJID, .keep_all = TRUE)

# Create table by treatment
baseline_all |> select(treatment, age, who, gleason_m8, psa) |> 
  tbl_summary(by = treatment,
              statistic = list(
                age ~ "{mean} ({sd})",                          
                psa ~ "{median} ({p25}, {p75})",                
                all_categorical() ~ "{n} ({p}%)"),
              label = list(
                age        ~ "Age, years",
                who        ~ "WHO performance status",
                gleason_m8 ~ "Gleason score ≥8",
                psa        ~ "PSA, µg/L"), missing_text = "Missing") |> add_overall(last = TRUE)

# 6. Export datasets ----------------------------------------------------
write_csv(baseline_all, here(results_directory, "PEACE3_ipd_for_MAIC.csv"))
saveRDS(baseline_all, here("data", "peace3_ipd_for_MAIC.RDS"))
