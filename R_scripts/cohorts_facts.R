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

#-------------------------------------------------------------------------------

HS_infections <- read_csv("./data/NMC/Salje_et_al_2021/infection_simulations.csv") |>
  filter(vaccine == 0) |>
  group_by(simID) |>
  summarise(n_infections = n())


HS_summary <- HS_infections |>
  summarise(
    q2.5  = quantile(n_infections, 0.025),
    mean  = mean(n_infections),
    q97.5 = quantile(n_infections, 0.975))

finalSubjectData <- read.table("./data/NMC/Salje_et_al_2021/SubjectData.extra_v2.txt")

colnames(finalSubjectData)<-c("newID","sero","vaccine","NoAssays",
                              "NoSympSerKnown","NoSympSerUnknown",
                              "delayStartStudy","delayLastBlood",
                              "delayStartVac3","delayStartVac1",
                              "delayLastAnnualDraw","SeroStatB1","seroY1")

HS_plac_individuals <-  finalSubjectData |>  filter(vaccine == 0) |>
  mutate(plac_id = row_number(),
         range_days = (delayLastBlood - delayStartVac3),
         range      =  range_days / 365)

HS_person_years <- sum(HS_plac_individuals$range)

HS_incidence <- HS_summary$mean / HS_person_years

plac_df <- NMC_get_placebo_data()

anchor_df <- plac_df |> filter(visit_no == 6) |>
  select(subjectNo, days_bleed) |>
  arrange(subjectNo) |>
  mutate(plac_id = row_number()) |>
  left_join(HS_plac_individuals[, c("plac_id", "range_days")],
            by = "plac_id") |>
  mutate(final_draw = days_bleed + range_days) |>
  rename(start_draw = days_bleed)

range_df <- plac_df |> select(subjectNo, days_bleed) |>
  left_join(anchor_df[, c("subjectNo", "start_draw", "final_draw")],
            by = "subjectNo") |>
  filter(days_bleed >= start_draw & days_bleed <= final_draw) |>
  group_by(subjectNo) |>
  summarise(range = (max(days_bleed) - min(days_bleed)) / 365)

person_years <- sum(range_df$range)

subset_inf_df <- NMC_get_infection_df() |>
  select(subjectNo, days_bleed, is_inf) |>
  left_join(anchor_df[, c("subjectNo", "start_draw", "final_draw")],
            by = "subjectNo") |>
  filter(days_bleed > start_draw & days_bleed <= final_draw)

n_inf <- sum(subset_inf_df$is_inf)
