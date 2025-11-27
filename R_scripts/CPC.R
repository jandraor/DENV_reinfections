read_CPC_data <- function()
{
  raw_data <- read_csv("./data/CPC/CPC Specimen List_Cohort_Household_3Sep14_TT_16JUN16.csv",
                       show_col_types = FALSE)

  raw_data |>
    rename(subject_id          = `COHORT SUBJECT NUMBER`,
           Dengue_Nested_PCR   = `Dengue Nested PCR`,
           date_acute_blood    = `DATE OF ACUTE BLOOD`,
           Y01_collection_date = `Y01 date of Collection`,
           Y02_collection_date = `Y02 Date of Collection`,
           Y03_collection_date = `Y03 Date of Collection`)
}

CPC_get_tidy_data <- function()
{
  clean_data <- read_CPC_data()

  collection_dates_df <- clean_data |>
    select(subject_id, contains("collection_date")) |>
    pivot_longer(-subject_id, values_to = "collection_date") |>
    mutate(time = as.numeric(str_sub(name, 2, 3)),
           collection_date = dmy(collection_date)) |>
    select(-name) |>
    mutate(key = paste(subject_id, time, sep = "_")) |>
    select(key, collection_date) |>
    unique()

  clean_data |>
    select(subject_id,
           age_yr,
           contains("HAI_D")) |>
    pivot_longer(matches("HAI_D"), values_to = "HAI") |>
    filter(!is.na(HAI)) |>
    mutate(HAI       = ifelse(HAI == "<10", 5, as.numeric(HAI)),
           log_2_HAI = log2_transform(HAI),
           log_HAI   = log(HAI),
           time      = as.numeric(str_sub(name, 2, 3)),
           age       = age_yr + time - 1,
           key       = paste(subject_id, time, sep = "_")) |>
    left_join(collection_dates_df, by = "key") |>
    group_by(subject_id, age, collection_date) |>
    summarise(log2_mean = mean(log_2_HAI),
              log_mean  = mean(log_HAI),
              .groups    = "drop")
}

CPC_get_infection_df <- function()
{
  CPC_infection_df <- CPC_get_tidy_data() |>
    group_by(subject_id) |>
    arrange(subject_id, age) |>
    mutate(lag_log2_titre   = lag(log2_mean),
           lag_log_titre    = lag(log_mean),
           ratio_log2       = ifelse(log2_mean == 0 & lag_log2_titre == 0,
                                     0, log2_mean / lag_log2_titre),
           is_titre_inf     = ifelse(ratio_log2 >= 1.6, TRUE, FALSE),
           start_interval   = lag(collection_date),
           end_interval     = collection_date) |>
    filter(!is.na(lag_log2_titre)) |>
    ungroup()

  symp_df <- CPC_get_symp_infections()

  df_list <- split(CPC_infection_df, CPC_infection_df$subject_id)

  CPC_infection_df <- imap_dfr(df_list, \(df, s_id) {

    subject_PCR <- symp_df |> filter(subject_id == s_id)

    df$PCR <- FALSE

    if(nrow(subject_PCR) > 0)
    {
      PCR_dates <- subject_PCR$date_acute_blood

      df$PCR <- sapply(
        seq_len(nrow(df)),
        function(i) any(PCR_dates >= df$start_interval[i] &
                          PCR_dates <= df$end_interval[i]))

    }

    df
  }) |>
    mutate(is_inf = is_titre_inf | PCR)


  CPC_infection_df
}

CPC_get_prob_inf_at_age <- function()
{
  infection_df <- CPC_get_infection_df()

  CPC_age_inf <- infection_df |> group_by(age) |>
    summarise(n_infections  = sum(is_inf),
              n_individuals = n()) |>
    mutate(pct = n_infections / n_individuals) |>
    filter(!is.na(age))
}

CPC_get_symp_infections <- function()
{
  read_CPC_data() |>
    select(subject_id, date_acute_blood, Dengue_Nested_PCR) |>
    filter(str_detect(Dengue_Nested_PCR, "DEN")) |>
    mutate(date_acute_blood = dmy(date_acute_blood))
}

CPC_get_prob_symp_inf <- function()
{
  infection_df <- CPC_get_infection_df()

  infection_df |> group_by(age) |>
    summarise(n_infections  = sum(is_inf),
              n_symp        = sum(PCR))
}
