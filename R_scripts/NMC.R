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

  titre_df  |> filter(vac_group == "Placebo") |>
    mutate(age = round(age_cyd_inclusion + days_bleed / 365, 0))
}

NMC_get_infection_df <- function()
{
  plac_df <- NMC_get_placebo_data()

  plac_ids <- unique(plac_df$subjectNo)

  df_list <- split(plac_df, plac_df$subjectNo)

  plac_inf <- NMC_get_symptomatic_infections(plac_ids)


  infection_df <- imap_dfr(df_list, \(f_df, id) {

    PCR_df <- plac_inf |>
      filter(subject_no == id)

    detect_infections(f_df |> rename(titre = log2_mean),
                      PCR_df, cutoff = 1.18) |>
      remove_multiple_measurements_in_a_year()
  }) |> select(subjectNo, collected_year, age, days_bleed, titre, log_mean,
               PCR_infection, titre_infection, infection, serotype) |>
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

 infection_df |>
   mutate(is_inf = ifelse(key %in% keys_first_year,
                          NA, infection))
}

NMC_get_prob_inf_at_age <- function()
{
  infection_df <- NMC_get_infection_df()

  infection_df |>
    filter(!is.na(is_inf)) |> group_by(age) |>
    summarise(n_infections  = sum(is_inf),
              n_individuals = n()) |>
    mutate(pct = n_infections / n_individuals)
}

NMC_get_prob_symp_inf <- function()
{
  infection_df <- NMC_get_infection_df()

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

NMC_get_symptomatic_infections <- function(plac_ids)
{
  PCR_infections <- read_csv("./data/NMC/NMC_PCR_infections_manual.csv",
                             show_col_types = FALSE) |>
    arrange(subject_no)

  plac_inf <- PCR_infections |> filter(subject_no %in% plac_ids)
}

