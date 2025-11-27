dmy2 <- function(date_str, threshold = 2025) {
  date <- dmy(date_str)               # parse normally
  year(date) <- ifelse(year(date) > threshold, year(date) - 100, year(date))
  date
}

KFCS_read_data <- function()
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
        .fns = ~na_if(., 0))
    ) |>
    mutate(
      across(
        .cols = matches("HAI"),
        .fns = list(log2 = log2_transform, log = log),
        .names = "{.fn}_{.col}"))

  KFCS_df |>
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
}

KFCS_get_infection_df <- function()
{
  KFCS_df <- KFCS_read_data()

  infection_df <- KFCS_df |>
    select(subjectNo, dateBirth2, start_interval, log2_mean_start_DEN,
           end_interval, log2_mean_end_DEN,
           log_mean_start_DEN, log_mean_end_DEN, xgBoostPred) |>
    mutate(age_start    = as.numeric(start_interval - dateBirth2) / 365.25,
           age_end      = as.numeric(end_interval - dateBirth2) / 365.25,
           age_round    = round((age_start + age_end) / 2),
           is_titre_inf = ifelse(xgBoostPred < 0.8, FALSE, TRUE))

  KFCS_PCR <-  KFCS_get_symp_infections()

  df_list <- split(infection_df, infection_df$subjectNo)


  infection_df <- imap_dfr(df_list, \(df, subject_id) {

    subject_PCR <- KFCS_PCR |> filter(subjectNo == subject_id)

    df$PCR <- FALSE

    if(nrow(subject_PCR) > 0)
    {
      PCR_dates <- subject_PCR$dateEvaluationA1

      df$PCR <- sapply(
        seq_len(nrow(df)),
        function(i) any(PCR_dates >= df$start_interval[i] &
                          PCR_dates <= df$end_interval[i]))

    }

    df
  }) |>
    mutate(is_inf = is_titre_inf | PCR)
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
    select(subjectNo, dateEvaluationA1, pcrResult) |>
    mutate(dateEvaluationA1 = ymd(dateEvaluationA1)) |>
    filter(pcrResult == "Dengue")
}

KFCS_get_prob_symp_inf <- function()
{
  infection_df <- KFCS_get_infection_df()

  infection_df |> group_by(age_round) |>
    summarise(n_infections  = sum(is_inf),
              n_symp        = sum(PCR)) |>
    filter(!is.na(age_round))
}

KFCS_estimate_deltas <- function(df)
{
  df <- df |> filter(!is_inf) # to remove short-term dynamics

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
           log_mean_start_DEN, log_mean_end_DEN, start_interval, end_interval)

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

    data.frame(subjectNo       = unique(df$subjectNo),
               dateBirth2      = unique(df$dateBirth2),
               collection_date = collection_dates,
               log2_mean       = log2_means) |>
      mutate(age = round(as.numeric(collection_date - dateBirth2) / 365.25))
  })
}
