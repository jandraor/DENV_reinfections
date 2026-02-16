source("./R_scripts/CPC.R")
source("./R_scripts/NMC.R")
source("./R_scripts/KFCS.R")

CPC_df <- CPC_get_tidy_data()

n_CPC <- CPC_df$subject_id |> unique() |> length()

NMC_df <- NMC_get_placebo_data()

n_NMC <- NMC_df$subjectNo |> unique() |> length()

KFCS_df <- KFCS_get_titre_data()

n_KFCS <- KFCS_df$subjectNo |> unique() |> length()

n_CPC + n_NMC + n_KFCS

#-------------------------------------------------------------------------------
CPC_range_df <- CPC_df |> group_by(subject_id) |>
  mutate(min_date = min(collection_date),
         max_date = max(collection_date)) |>
  select(subject_id, min_date, max_date) |>
  ungroup() |>
  unique() |>
  mutate(range = time_length(interval(min_date, max_date), "years"))

CPC_person_years <- sum(CPC_range_df$range)

CPC_range_df$range |> mean() |> round(1)

NMC_range_df <- NMC_df |> group_by(subjectNo) |>
  mutate(min_date = min(collected_date),
         max_date = max(collected_date)) |>
  select(subjectNo, min_date, max_date) |>
  ungroup() |>
  unique() |>
  mutate(range = time_length(interval(min_date, max_date), "years"))

max_follow_up <- max(NMC_range_df$range)

person_years <- sum(NMC_range_df$range)

NMC_range_df$range |> mean() |> round(1)

KFCS_range_df <- KFCS_df |>
  group_by(subjectNo) |>
  mutate(min_date = min(collection_date),
         max_date = max(collection_date)) |>
  select(subjectNo, min_date, max_date) |>
  ungroup() |>
  unique() |>
  mutate(range = time_length(interval(min_date, max_date), "years"))

max_follow_up <- max(KFCS_range_df$range)

KFCS_range_df$range |> mean() |> round(1)

KFCS_person_years <- sum(KFCS_range_df$range)

cat("\n KFCS person years: ", KFCS_person_years)

# Age inclusion-----------------------------------------------------------------
NMC_range_inclusion <- NMC_df$age_cyd_inclusion |> range()

CPC_range_inclusion <- CPC_df |> group_by(subject_id) |>
  filter(age == min(age)) |> ungroup() |>
  pull(age) |> range()

KFCS_enrolment <- KFCS_df |> group_by(subjectNo) |>
  filter(age == min(age)) |> ungroup()

KFCS_enrolment |>
  mutate(serostatus = ifelse(log2_mean == 0, "seronegative", "seropositive")) |>
  group_by(serostatus) |> summarise(n = n()) |> mutate(pct = n / sum(n))

KFCS_range_inclusion <- KFCS_enrolment |> pull(age) |> range()

#-------------------------------------------------------------------------------

plac_df |> select(subjectNo, starts_with("log2_D")) |>
  pivot_longer(-subjectNo) |>
  mutate(below_LOD = ifelse(value == 0, TRUE, FALSE)) |>
  group_by(below_LOD) |>
  summarise(n = n()) |>
  mutate(pct = n / sum(n))

#-------------------------------------------------------------------------------

n_all_5 <- KFCS_df |> filter(log2_mean == 0) |> nrow()

total_meas <- nrow(KFCS_df)

n_all_5 / total_meas

#-------------------------------------------------------------------------------

KFCS_df |> select(subjectNo, starts_with("log2_D")) |>
  pivot_longer(-subjectNo) |>
  mutate(below_LOD = ifelse(value == 0, TRUE, FALSE)) |>
  group_by(below_LOD) |>
  summarise(n = n()) |>
  mutate(pct = n / sum(n))

