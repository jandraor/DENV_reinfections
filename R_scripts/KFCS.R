dmy2 <- function(date_str, threshold = 2025) {
  date <- dmy(date_str)               # parse normally
  year(date) <- ifelse(year(date) > threshold, year(date) - 100, year(date))
  date
}

# Impute value: Value assigned to measurements <10. In this dataset, values are
#   already imputed to 5.
KFCS_read_data <- function(imputed_value = NULL)
{
  KFCS_df <- read_csv("./data/KFCS/KFCS_longData.csv",
           show_col_types = FALSE) |>
    mutate(
      start_interval = dmy(START_dateSpecimenCollected),
      end_interval   = dmy(END_dateSpecimenCollected),
      dateBirth2 = ifelse(str_detect(dateBirth, "dd"),
                          str_glue("01-01-{str_sub(dateBirth, -4, -1)}"),
                          dateBirth)) |>
    mutate(dateBirth2 = dmy2(dateBirth2)) |>
    mutate(
      across(
        .cols = matches("HAI"),
        .fns = ~na_if(., 0)))

  if(!is.null(imputed_value))
  {
    KFCS_df <- KFCS_df |>
      mutate(across(matches("HAI"), ~ if_else(.x == 5, imputed_value, .x)))
  }

  KFCS_df <- KFCS_df |>
    mutate(
      across(
        .cols = matches("HAI"),
        .fns = list(log2 = log2_transform, log = log),
        .names = "{.fn}_{.col}"))

  KFCS_df <- KFCS_df |>
    mutate(log2_mean_start_DEN =
             rowMeans(select(KFCS_df, starts_with("log2_START_HAI_DEN")),
                      na.rm = TRUE),
           log_mean_start_DEN =
             rowMeans(select(KFCS_df, starts_with("log_START_HAI_DEN")),
                      na.rm = TRUE),
           log2_mean_end_DEN =
             rowMeans(select(KFCS_df, starts_with("log2_END_HAI_DEN")),
                      na.rm = TRUE),
           log_mean_end_DEN =
             rowMeans(select(KFCS_df, starts_with("log_END_HAI_DEN")),
                      na.rm = TRUE))

  KFCS_seroneg_ids <- KFCS_df |>
    select(subjectNo, start_interval, log2_mean_start_DEN) |>
    group_by(subjectNo) |>
    filter(start_interval == min(start_interval),
           log2_mean_start_DEN == 0) |> pull(subjectNo) |> unique()

  KFCS_df <- KFCS_df |>
    mutate(serostatus = ifelse(subjectNo %in% KFCS_seroneg_ids,
                               "seronegative", "seropositive"))

  KFCS_df
}

KFCS_get_infection_df <- function(imputed_val = NULL)
{
  KFCS_df <- KFCS_read_data(imputed_val)

  infection_df <- KFCS_df |>
    select(subjectNo, serostatus, dateBirth2, start_interval,
           log2_mean_start_DEN, end_interval, log2_mean_end_DEN,
           log_mean_start_DEN, log_mean_end_DEN, xgBoostPred) |>
    mutate(age_start    = as.numeric(start_interval - dateBirth2) / 365.25,
           age_end      = as.numeric(end_interval - dateBirth2) / 365.25,
           age_round    = round(age_end),
           is_titre_inf = ifelse(xgBoostPred < 0.8, FALSE, TRUE),
           ratio_log2   = ifelse(log2_mean_end_DEN == 0 & log2_mean_start_DEN == 0,
                                     0, log2_mean_end_DEN / log2_mean_start_DEN),
           is_titre_alt_inf = ifelse(ratio_log2 >= 1.6, TRUE, FALSE))


  KFCS_PCR <-  KFCS_get_symp_infections()

  df_list <- split(infection_df, infection_df$subjectNo)


  infection_df <- imap_dfr(df_list, \(df, subject_id) {

    subject_PCR <- KFCS_PCR |> filter(subjectNo == subject_id)

    df$PCR <- FALSE

    df$serotype <- "Subclinical"

    if(nrow(subject_PCR) > 0)
    {
      PCR_dates <- subject_PCR$dateEvaluationA1

      df$PCR <- sapply(
        seq_len(nrow(df)),
        function(i) any(PCR_dates >= df$start_interval[i] &
                          PCR_dates <= df$end_interval[i]))

      df[df$PCR, "serotype"] <- subject_PCR$serotype
    }

    df
  }) |>
    mutate(is_inf = is_titre_inf | PCR,
           is_alt_inf = is_titre_alt_inf | PCR)
}

KFCS_get_prob_inf_at_age <- function()
{
  infection_df <- KFCS_get_infection_df()


  KFCS_age_inf <- infection_df |> group_by(age_round) |>
    summarise(n_infections  = sum(is_inf),
              n_individuals = n()) |>
    mutate(pct = n_infections / n_individuals) |>
    filter(!is.na(age_round))
}

KFCS_plot_prob_inf_age <- function()
{
  KFCS_age_inf <- KFCS_get_prob_inf_at_age()

  ggplot(KFCS_age_inf |> filter(age_round >=2), aes(age_round, pct)) +
    geom_line() +
    theme_classic()
}

KFCS_get_symp_infections <- function()
{
  KFCS_PCR <- read_csv("./data/KFCS/Analysis_Illness_20240722.csv",
                       show_col_types = FALSE) |>
    filter(pcrDengueTypeDEN1 == 1 |
             pcrDengueTypeDEN2 == 1 |
             pcrDengueTypeDEN3 == 1 |
             pcrDengueTypeDEN4 == 1 |
             pcrResult == "Dengue") |>
    select(subjectNo, dateEvaluationA1, pcrResult, contains("pcrDengueType")) |>
    mutate(DEN_sum = rowSums(across(starts_with("pcrDengueType")))) |>
    mutate(dateEvaluationA1 = ymd(dateEvaluationA1),
           serotype = case_when(
             pcrDengueTypeDEN1 == 1 & DEN_sum == 1 ~ "DENV-1",
             pcrDengueTypeDEN2 == 1 & DEN_sum == 1 ~ "DENV-2",
             pcrDengueTypeDEN3 == 1 & DEN_sum == 1 ~ "DENV-3",
             pcrDengueTypeDEN4 == 1 & DEN_sum == 1 ~ "DENV-4",
             TRUE ~ "Untyped/Mixed")) |>
    select(-contains("pcrDengueType"))
}

KFCS_get_rise_df <- function()
{
  PCR_df <- read_csv("./data/KFCS/Analysis_Illness_20240722.csv",
                     show_col_types = FALSE) |>
    filter(pcrResult == "Dengue")

  df <- PCR_df |>
    select(subjectNo, id, haiBA1_sampleDate, contains("haiBA1_D")) |>
    filter(haiBA1_sampleDate != "NULL") |>
    mutate(
      across(
        .cols = matches("haiBA1_D"),
        .fns = as.numeric)) |>
    mutate(
      across(
        .cols = matches("haiBA1_D"),
        .fns = list(log2 = log2_transform),
        .names = "{.fn}_{.col}"))

  acute_df <- df |>
    mutate(log2_mean =
             rowMeans(select(df, starts_with("log2_haiBA1_D")),
                      na.rm = TRUE)) |>
    rename(sample_date = haiBA1_sampleDate) |>
    select(subjectNo, id, sample_date, log2_mean) |>
    mutate(sample_type = "A1")

  df <- PCR_df |>
    select(subjectNo, id, haiBC1_sampleDate, contains("haiBC1_D")) |>
    filter(haiBC1_sampleDate != "NULL") |>
    mutate(
      across(
        .cols = matches("haiBC1_D"),
        .fns = as.numeric)) |>
    mutate(
      across(
        .cols = matches("haiBC1_D"),
        .fns = list(log2 = log2_transform),
        .names = "{.fn}_{.col}"))

  conv_1_df <- df |>
    mutate(log2_mean =
             rowMeans(select(df, starts_with("log2_haiBC1_D")),
                      na.rm = TRUE)) |>
    rename(sample_date = haiBC1_sampleDate) |>
    select(subjectNo, id, sample_date, log2_mean) |>
    mutate(sample_type = "C1")

  df <- PCR_df |>
    select(subjectNo, id, haiBC2_sampleDate, contains("haiBC2_D")) |>
    filter(haiBC2_sampleDate != "NULL") |>
    mutate(
      across(
        .cols = matches("haiBC2_D"),
        .fns = as.numeric)) |>
    mutate(
      across(
        .cols = matches("haiBC2_D"),
        .fns = list(log2 = log2_transform),
        .names = "{.fn}_{.col}"))

  conv_2_df <- df |>
    mutate(log2_mean =
             rowMeans(select(df, starts_with("log2_haiBC2_D")),
                      na.rm = TRUE)) |>
    rename(sample_date = haiBC2_sampleDate) |>
    select(subjectNo, id, sample_date, log2_mean) |>
    mutate(sample_type = "C2")


  sample_df <- bind_rows(acute_df, conv_1_df, conv_2_df) |>
    mutate(sample_date = ymd(sample_date))

  titre_df <- KFCS_get_titre_data()

  df_list <- split(sample_df, sample_df$id)

  rise_df <- map_df(df_list, \(df){

    peak_df <- df |> filter(log2_mean == max(log2_mean)) |>
      filter(sample_date == min(sample_date))

    acute_df <- df |> filter(sample_type == "A1")

    infection_time <- acute_df$sample_date

    s_id <- unique(df$subjectNo)

    s_df <- titre_df |> filter(subjectNo == s_id)

    if(nrow(s_df) == 0) return(NULL)

    baseline_df <- s_df |>
      filter(collection_date < infection_time) |>
      filter(collection_date == max(collection_date))

    baseline_date <- baseline_df$collection_date

    gap <- as.numeric(infection_time - baseline_date)

    if(gap > 365) return(NULL)


    data.frame(subject_id     = s_id ,
               baseline       = baseline_df$log2_mean,
               baseline_time  = baseline_date,
               peak           = peak_df$log2_mean,
               peak_time      = peak_df$sample_date,
               infection_time = infection_time)  |>
      mutate(rise   = peak - baseline,
             delta_t = peak_time - baseline_time,
             serostatus = ifelse(baseline == 0, "Seronegative", "Seropositive"))
  })

  rise_df
}

KFCS_get_PCR_trajectories <- function()
{
  rise_df <- KFCS_get_rise_df() |>
    mutate(id = row_number())

  rise_list <- split(rise_df, rise_df$id)

  infection_df <- KFCS_get_infection_df() |>
    select(subjectNo, start_interval, log2_mean_start_DEN, end_interval,
           log2_mean_end_DEN, is_inf)

  map_df(rise_list, \(df) {

    s_id      <- df$subject_id
    peak_time <- df$peak_time
    id        <- df$id

    s_df <- infection_df |> filter(subjectNo == s_id) |>
      mutate(is_peak_here = peak_time > start_interval &
               peak_time <= end_interval) |>
      filter(end_interval > peak_time) |>
      mutate(inf_counter = cumsum(is_inf)) |>
      filter(inf_counter == 1)


    if(nrow(s_df) == 0) return(NULL)

    second_meas <- s_df |> slice(1) |>
      pull(end_interval)

    second_titre <- s_df |> slice(1) |>
      pull(log2_mean_end_DEN)

    rest_df <- s_df |> slice(-1)

    data.frame(id = id,
               subjectNo       = s_id,
               collection_date = c(peak_time, second_meas, rest_df$end_interval),
               log2_mean       = c(df$peak, second_titre, rest_df$log2_mean_end_DEN)) |>
      mutate(years_post_peak = as.numeric(collection_date - min(collection_date)) / 365)

  })

}

KFCS_get_prob_symp_inf <- function()
{
  infection_df <- KFCS_get_infection_df()

  infection_df |> group_by(age_round) |>
    summarise(n_infections  = sum(is_inf),
              n_symp        = sum(PCR),
              n_individuals = n()) |>
    filter(!is.na(age_round))
}

KFCS_estimate_deltas <- function(df)
{
  # to remove short-term dynamics
  df <- df |>
    mutate(prev_meas_inf = lag(is_inf, default = FALSE)) |>
    filter(!is_inf) |> filter(!prev_meas_inf)

  n_rows <- nrow(df)

  if(n_rows == 0) return(NULL)

  df <- df |> arrange(start_interval) |>
    rename(start_titre = log2_mean_start_DEN,
           end_titre   = log2_mean_end_DEN) |>
    mutate(delta_time  = time_length(end_interval - start_interval, "years"),
           delta_titre = end_titre  - start_titre) |>
    select(inf_id, start_interval, start_titre, end_interval, end_titre,
           delta_time, delta_titre)

  if(n_rows == 1) return(df)

  combos <- t(combn(seq_len(n_rows), 2))

  df |>
    bind_rows(
      data.frame(
        inf_id         = unique(df$inf_id),
        start_interval = df$start_interval[combos[, 1]],
        start_titre    = df$start_titre[combos[, 1]],
        end_interval   = df$end_interval[combos[, 2]],
        end_titre      = df$end_titre[combos[, 2]]) |>
        mutate(
          delta_time  = time_length(end_interval - start_interval, "years"),
          delta_titre = end_titre  - start_titre))
}

KFCS_estimate_deltas_PCR <- function(df) {

  df <- df |> filter(!PCR) # to remove short-term dynamics

  n_rows <- nrow(df)

  if(n_rows == 0) return(NULL)

  df <- df |> arrange(start_interval) |>
    rename(start_titre = log2_mean_start_DEN,
           end_titre   = log2_mean_end_DEN) |>
    mutate(delta_time  = time_length(end_interval - start_interval, "years"),
           delta_titre = end_titre  - start_titre) |>
    select(inf_id, start_interval, start_titre, end_interval, end_titre,
           delta_time, delta_titre)

  if(n_rows == 1) return(df)

  combos <- t(combn(seq_len(n_rows), 2))

  df |>
    bind_rows(
      data.frame(
        inf_id         = unique(df$inf_id),
        start_interval = df$start_interval[combos[, 1]],
        start_titre    = df$start_titre[combos[, 1]],
        end_interval   = df$end_interval[combos[, 2]],
        end_titre      = df$end_titre[combos[, 2]]) |>
        mutate(
          delta_time  = time_length(end_interval - start_interval, "years"),
          delta_titre = end_titre  - start_titre))
}

KFCS_get_titre_data <- function()
{
  KFCS_df <- KFCS_read_data() |>
    select(subjectNo, dateBirth2, log2_mean_start_DEN, log2_mean_end_DEN,
           log_mean_start_DEN, log_mean_end_DEN, start_interval, end_interval,
           starts_with("log2_START_HAI_DEN"),
           starts_with("log2_END_HAI_DEN"))

  df_list <- split(KFCS_df, KFCS_df$subjectNo)

  map_df(df_list, \(df) {

    first_row_df <- df |> slice(1)

    rest_df <- df |> slice(-1)

    collection_dates <- c(first_row_df$start_interval[1],
                          first_row_df$end_interval[1],
                          rest_df$end_interval)

    log2_means <- c(first_row_df$log2_mean_start_DEN[1],
                    first_row_df$log2_mean_end_DEN[1],
                    rest_df$log2_mean_end_DEN)

    log2_D1_vals <- c(first_row_df$log2_START_HAI_DEN1[1],
                      first_row_df$log2_END_HAI_DEN1[1],
                      rest_df$log2_END_HAI_DEN1)

    log2_D2_vals <- c(first_row_df$log2_START_HAI_DEN2[1],
                      first_row_df$log2_END_HAI_DEN2[1],
                      rest_df$log2_END_HAI_DEN2)

    log2_D3_vals <- c(first_row_df$log2_START_HAI_DEN3[1],
                      first_row_df$log2_END_HAI_DEN3[1],
                      rest_df$log2_END_HAI_DEN3)

    log2_D4_vals <- c(first_row_df$log2_START_HAI_DEN4[1],
                      first_row_df$log2_END_HAI_DEN4[1],
                      rest_df$log2_END_HAI_DEN4)

    data.frame(subjectNo       = unique(df$subjectNo),
               dateBirth2      = unique(df$dateBirth2),
               collection_date = collection_dates,
               log2_mean       = log2_means,
               log2_D1         = log2_D1_vals,
               log2_D2         = log2_D2_vals,
               log2_D3         = log2_D3_vals,
               log2_D4         = log2_D4_vals) |>
      mutate(age = round(as.numeric(collection_date - dateBirth2) / 365.25))
  })
}

KFCS_get_decay_rates <- function()
{
  KFCS_decay <- KFCS_get_infection_df() |>
    group_by(subjectNo) |>
    arrange(subjectNo, start_interval) |>
    mutate(inf_idx = cumsum(is_inf)) |>
    filter(!(serostatus == "seronegative" & inf_idx == 0)) |>
    mutate(prev_meas_inf = lag(is_inf, default = FALSE)) |>
    ungroup() |>
    filter(!is_inf) |>
    filter(!prev_meas_inf) |>
    mutate(titre_diff = log_mean_end_DEN - log_mean_start_DEN,
           time_diff  = as.numeric(end_interval - start_interval) / 365.25,
           slope      = titre_diff / time_diff)|>
    rename(age = age_round) |>
    filter(age  >= 2,
           between(time_diff, 0.9, 2))

  KFCS_decay
}
