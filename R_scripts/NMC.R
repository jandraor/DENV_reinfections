format_NMC_titre <- function(raw_df)
{
  titre_df <- raw_df |>
    mutate(
      across(
        .cols = matches("^D\\d"),
        .fns  = ~ suppressWarnings(if_else(.x == "<10.0", 5, as.numeric(.x)))
      )) |>
    mutate(across(
      .cols = matches("^D\\d"),
      .fns = log2_transform,
      .names = "log2_{.col}")) |>
    mutate(across(
      .cols = matches("^D\\d"),
      .fns = log,
      .names = "log_{.col}")) |>
    mutate(collected_date = dmy(collected_date),
           date_vaccine1  = dmy(date_vaccine1),
           days_bleed     = ifelse(is.na(days_bleed) & timepoint == "BL0004",
                                   730, days_bleed),
           # dplyr::if_else() keeps the date format
           collected_date = dplyr::if_else(
             is.na(collected_date) & timepoint == "BL0004",
             date_vaccine1 + days(days_bleed),
             collected_date),
           collected_year = year(collected_date))

  titre_df <- titre_df |>
    mutate(log2_mean = rowMeans(select(titre_df, starts_with("log2_D")),
                                na.rm = TRUE),
           log_mean  = rowMeans(select(titre_df, starts_with("log_D")),
                                na.rm = TRUE))
}

NMC_get_placebo_data <- function()
{
  raw_data       <- read_csv("./data/NMC/NMC_stats_AFRIMS_PRNT-2025-01-15.csv",
                             show_col_types = FALSE)

  titre_df <- format_NMC_titre(raw_data)

  seroneg_ids <- titre_df |>  group_by(subjectNo) |>
    filter(days_bleed == min(days_bleed)) |>
    filter(log2_mean == 0, vac_group == "Placebo") |> pull(subjectNo)

  plac_df <- titre_df |> filter(vac_group == "Placebo") |>
    mutate(age = round(age_cyd_inclusion + days_bleed / 365, 0),
           serostatus = ifelse(subjectNo %in% seroneg_ids,
                               "seronegative",
                               "seropositive"))

  plac_df
}

NMC_get_infection_df <- function(cut_off)
{
  fp <- str_glue("./data/NMC/NMC_infection_{cut_off}.rds")

  if(!file.exists(fp))
  {
    plac_df <- NMC_get_placebo_data()

    plac_ids <- unique(plac_df$subjectNo)

    df_list <- split(plac_df, plac_df$subjectNo)

    plac_inf <- NMC_get_symptomatic_infections(plac_ids)

    infection_df <- imap_dfr(df_list, \(f_df, id) {

      PCR_df <- plac_inf |>
        filter(subject_no == id)

      detect_infections(f_df |> rename(titre = log2_mean),
                        PCR_df, cutoff = cut_off) |>
        remove_multiple_measurements_in_a_year()
    }) |> select(subjectNo, serostatus, collected_year, age, days_bleed, titre,
                 log_mean, PCR_infection, titre_infection, infection, serotype,
                 contains("log2_D")) |>
      mutate(key = paste(subjectNo, collected_year, sep = "_"))

    # Individuals for whom it is not possible to determine whether there was
    # an infection during the first year because there was only one measurement
    keys_first_year <- plac_df |> group_by(subjectNo) |>
      filter(collected_year == min(collected_year)) |>
      group_by(subjectNo, collected_year) |>
      count() |> arrange(desc(n)) |> ungroup() |>
      filter(n == 1) |>
      mutate(key = paste(subjectNo, collected_year, sep = "_")) |>
      pull(key)

    infection_df <- infection_df |>
      mutate(is_inf = ifelse(key %in% keys_first_year,
                             NA, infection))

    saveRDS(infection_df, fp)
  } else infection_df <- readRDS(fp)

  infection_df
}

NMC_get_prob_inf_at_age <- function()
{
  infection_df <- NMC_get_infection_df(1.18)

  infection_df |>
    filter(!is.na(is_inf)) |> group_by(age) |>
    summarise(n_infections  = sum(is_inf),
              n_individuals = n()) |>
    mutate(pct = n_infections / n_individuals)
}

NMC_get_prob_symp_inf <- function()
{
  infection_df <- NMC_get_infection_df(1.18)

  infection_df |>
    filter(!is.na(is_inf)) |> group_by(age) |>
    summarise(n_infections  = sum(is_inf),
              n_symp        = sum(PCR_infection),
              n_individuals = n())
}

filter_short_term_dynamics <- function(df, plac_inf, cutoff_st = 365)
{
  first_meas_time <- df$days_bleed[[1]]

  if(df$infection[[1]] == 1)
  {
    df <- df |> filter((days_bleed - first_meas_time) > cutoff_st)
  } else
  {
    subject_id <- unique(df$subjectNo)

    subject_PCR_df <- plac_inf |> filter(subject_no == subject_id)

    if(nrow(subject_PCR_df) > 0)
    {
      # time since confirmation
      time_since_conf <- first_meas_time - subject_PCR_df$dengue_days_pd1

      time_since_conf <- suppressWarnings(min(time_since_conf[time_since_conf >= 0]))

      if(is.finite(time_since_conf) && time_since_conf <= 365)
      {
        df <- df |>
          filter((days_bleed - first_meas_time + time_since_conf) > cutoff_st)
      }

    }
  }

  df
}

estimate_delta_times <- function(df)
{
  n_meas <- nrow(df)

  if(n_meas == 1) return(NULL)

  combos <- t(combn(seq_len(n_meas), 2))

  data.frame(
    inf_id = unique(df$inf_id),
    days_bleed_1 = df$days_bleed[combos[, 1]],
    titre_1      = df$titre[combos[, 1]],
    days_bleed2 = df$days_bleed[combos[, 2]],
    titre_2      = df$titre[combos[, 2]]) |>
    mutate(delta_time = days_bleed2 - days_bleed_1,
           delta_titre = titre_2 - titre_1)

}

pairwise_delta_by_serotype <- function(df)
{
  n_meas <- nrow(df)

  if(n_meas == 1) return(NULL)

  combos <- t(combn(seq_len(n_meas), 2))

  data.frame(
    inf_id = unique(df$inf_id),
    days_bleed_1 = df$days_bleed[combos[, 1]],
    titre_D1_1   = df$log2_D1_NT[combos[, 1]],
    titre_D2_1   = df$log2_D2_NT[combos[, 1]],
    titre_D3_1   = df$log2_D3_NT[combos[, 1]],
    titre_D4_1   = df$log2_D4_NT[combos[, 1]],
    days_bleed_2 = df$days_bleed[combos[, 2]],
    titre_D1_2   = df$log2_D1_NT[combos[, 2]],
    titre_D2_2   = df$log2_D2_NT[combos[, 2]],
    titre_D3_2   = df$log2_D3_NT[combos[, 2]],
    titre_D4_2   = df$log2_D4_NT[combos[, 2]]) |>
    mutate(delta_time     = days_bleed_2 - days_bleed_1,
           delta_titre_D1 = titre_D1_2 - titre_D1_1,
           delta_titre_D2 = titre_D2_2 - titre_D2_1,
           delta_titre_D3 = titre_D3_2 - titre_D3_1,
           delta_titre_D4 = titre_D4_2 - titre_D4_1)
}

NMC_get_symptomatic_infections <- function(plac_ids)
{
  PCR_infections <- read_csv("./data/NMC/NMC_PCR_infections_manual.csv",
                             show_col_types = FALSE) |>
    arrange(subject_no)

  plac_inf <- PCR_infections |> filter(subject_no %in% plac_ids)
}

NMC_get_annual_pairs <- function()
{
  cut_off <- 1.18

  plac_df <- NMC_get_placebo_data()

  plac_ids <- unique(plac_df$subjectNo)

  df_list <- split(plac_df, plac_df$subjectNo)

  plac_inf <- NMC_get_symptomatic_infections(plac_ids)

  infection_detection_df <- plac_df |> filter(visit_no >= 6 ) |>
    arrange(subjectNo, days_bleed) |>
    group_by(subjectNo) |>
    mutate(prev_titre_log2          = lag(log2_mean),
           titre_diff               = log2_mean - prev_titre_log2,
           is_titre_inf             = ifelse(titre_diff > cut_off, 1, 0),
           previous_collection_date = lag(collected_date)) |>
    ungroup()

  df_list <- split(infection_detection_df, infection_detection_df$subjectNo)

  PCR_df <- NMC_get_symptomatic_infections(plac_ids)

  infection_detection_df <- imap_dfr(df_list, \(df, subject_id) {

    subject_PCR <- PCR_df  |> filter(subject_no == subject_id)

    df$PCR <- FALSE

    df$serotype <- "Subclinical"

    if(nrow(subject_PCR) > 0)
    {
      vacc_date <- df$date_vaccine1[1]

      first_meas <- min(df$collected_date)

      subject_PCR <- subject_PCR |>
        mutate(inf_date = vacc_date + dengue_days_pd1) |>
        filter(inf_date >= first_meas)

      PCR_dates <- subject_PCR$inf_date

      df$PCR <- sapply(
        seq_len(nrow(df)),
        function(i) any(PCR_dates >= df$previous_collection_date[i] &
                          PCR_dates <= df$collected_date[i]))

      df[df$PCR, "serotype"] <- subject_PCR$serotype
    }

    df
  }) |>
    mutate(is_inf = is_titre_inf | PCR) |>
    filter(!is.na(is_titre_inf)) |>
    mutate(result = case_when(is_inf & PCR ~ "PCR Conf",
                              is_inf & !PCR ~ "Subclinical",
                              !is_inf ~ "No infection"))
}

NMC_get_binned_decay <- function(cut_off_val)
{
  infection_detection_df <- NMC_get_infection_df(cut_off_val)

  infection_df <- infection_detection_df |>
    mutate(is_inf = ifelse(is.na(is_inf), 0, is_inf)) |>
    group_by(subjectNo) |>
    mutate(inf_idx = cumsum(is_inf),
           inf_id  = paste(subjectNo, inf_idx, sep = "_")) |>
    ungroup()

  df_list <- split(infection_df, infection_df$inf_id)

  plac_ids <- unique(infection_detection_df$subjectNo)

  plac_inf <- NMC_get_symptomatic_infections(plac_ids)

  delta_df <- map_df(df_list, \(df) {

    df <- filter_short_term_dynamics(df, plac_inf, 365)

    if(nrow(df) <= 1) return(NULL)

    estimate_delta_times(df)
  }) |> mutate(delta_year = delta_time / 365,
               bin_delta  = round(delta_year, 0))

  NMC_mean_estimates <- delta_df |> add_bootstrap_CI() |>
    mutate(cohort = "NMC")

  NMC_mean_estimates
}

NMC_get_PCR_rise <- function()
{
  plac_df <- NMC_get_placebo_data()

  plac_ids <- unique(plac_df$subjectNo)

  plac_inf <- NMC_get_symptomatic_infections(plac_ids)

  tol <- 365 # tolerance in days

  plac_infections_list <- transpose(plac_inf)

  result_df <- imap_dfr(plac_infections_list, \(inf_obj, i) {

    subject_id <- inf_obj$subject_no
    inf_time   <- inf_obj$dengue_days_pd1

    subject_titre <- plac_df |>
      filter(subjectNo == subject_id) |>
      arrange(days_bleed) |>
      select(subjectNo, timepoint, days_bleed, visit_no, log2_mean)

    # Number of illness investigation blood draws
    n_illness_inv <- str_detect(subject_titre$timepoint, "BLA|BLC") |> sum()

    if(n_illness_inv == 0)
    {
      draws_time     <- subject_titre$days_bleed
      meas_after_inf <- draws_time[draws_time - inf_time > 0]

      next_meas <- min(meas_after_inf)

      if(next_meas - inf_time > tol)
      {
        cat("\n-------Rejection I---------")
        cat("\n Subject id: ", subject_id)
        cat("\n Infection time: ", inf_time)
        cat("\n Next meas after infection: ", next_meas)
        cat("\n Time to next meas: ", next_meas - inf_time)
        return(NULL)
      }


      rise_df <- data.frame(subject_id     = subject_id,
                            baseline       = NA,
                            baseline_time  = NA,
                            peak           = NA,
                            peak_time      = inf_time,
                            infection_time = inf_time) |>
        mutate(rise     = NA,
               delta_t = NA,
               serostatus = NA)
      return(rise_df)
    }

    titre_below_30 <- subject_titre |> filter((days_bleed - inf_time) < 30) |>
      arrange(days_bleed)

    if(nrow(titre_below_30) == 0)
    {
      draws_time     <- subject_titre$days_bleed
      meas_after_inf <- draws_time[draws_time - inf_time > 0]

      next_meas <- min(meas_after_inf)

      if(next_meas - inf_time > tol)
      {
        cat("\n-------Rejection II---------")
        cat("\n Subject id: ", subject_id)
        cat("\n Infection time: ", inf_time)
        cat("\n Next meas after infection: ", next_meas)
        cat("\n Time to next meas: ", next_meas - inf_time)
        return(NULL)
      }

      rise_df <- data.frame(subject_id     = subject_id,
                            baseline       = NA,
                            baseline_time  = NA,
                            peak           = NA,
                            peak_time      = inf_time,
                            infection_time = inf_time) |>
        mutate(rise       = NA,
               delta_t    = NA,
               serostatus = NA)

      return(rise_df)
    }

    peak_df <- titre_below_30 |> filter(log2_mean == max(log2_mean)) |>
      select(subjectNo, days_bleed, log2_mean, timepoint)

    if(!str_detect(peak_df$timepoint, "BLA|BLC"))
    {
      draws_time     <- subject_titre$days_bleed
      meas_after_inf <- draws_time[draws_time - inf_time > 0]
      next_meas      <- min(meas_after_inf)

      if(next_meas - inf_time > tol)
      {
        cat("\n-------Rejection III---------")
        cat("\n Subject id: ", subject_id)
        cat("\n Infection time: ", inf_time)
        cat("\n Next meas after infection: ", next_meas)
        cat("\n Time to next meas: ", next_meas - inf_time)
        return(NULL)
      }

      rise_df <- data.frame(subject_id     = subject_id,
                            baseline       = NA,
                            baseline_time  = NA,
                            peak           = NA,
                            peak_time      = inf_time,
                            infection_time = inf_time) |>
        mutate(rise       = NA,
               delta_t    = NA,
               serostatus = NA)

      return(rise_df)
    }

    peak_time <- peak_df$days_bleed

    baseline_df <- titre_below_30 |> filter(days_bleed < peak_time,
                                            peak_time - days_bleed < tol) |>
      filter(log2_mean == min(log2_mean))

    rise_df <- data.frame(subject_id     = subject_id,
                          baseline       = baseline_df$log2_mean,
                          baseline_time  = baseline_df$days_bleed,
                          peak           = peak_df$log2_mean,
                          peak_time      = peak_time,
                          infection_time = inf_time) |>
      mutate(rise   = peak - baseline,
             delta_t = peak_time - baseline_time,
             serostatus = ifelse(baseline == 0, "Seronegative", "Seropositive"))

  })
}

NMC_get_decay_rates <- function()
{
  infection_detection_df <- NMC_get_infection_df(1.18)

  infection_df <- infection_detection_df |>
    mutate(is_inf = ifelse(is.na(is_inf), 0, is_inf)) |>
    group_by(subjectNo) |>
    mutate(inf_idx = cumsum(is_inf),
           inf_id  = paste(subjectNo, inf_idx, sep = "_")) |>
    ungroup() |>
    # Removes measurements from seronegative individuals before their first infection.
    filter(!(serostatus == "seronegative" & inf_idx == 0))

  df_list <- split(infection_df, infection_df$inf_id)

  df_list <- unname(df_list)

  decay_df <- imap_dfr(df_list, \(df, i) {

    df <- filter_short_term_dynamics(df, plac_inf)

    if(nrow(df) <= 1) return(NULL)

    df |>
      mutate(prev_log_mean  = lag(log_mean),
             titre_diff     = log_mean - prev_log_mean,
             time_diff      = (days_bleed - lag(days_bleed)) / 365,
             slope          = titre_diff / time_diff) |>
      filter(!is.na(slope),  between(time_diff, 0.9, 2))
  })
}

NMC_get_measurements_summary <- function()
{
  fn <- "./data/NMC/NMC_measurements_summary.rds"

  if(!file.exists(fn))
  {
    df <- NMC_get_placebo_data() |> select(starts_with("log2_D"))

    means <- colMeans(df, na.rm = TRUE)
    sds   <- sapply(df, sd, na.rm = TRUE)

    summary_list <- list(means = means,
                         sds   = sds)

    saveRDS(summary_list, fn)

  } else summary_list <- readRDS(fn)

  summary_list
}
