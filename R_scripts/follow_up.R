detect_infections <- function(df, PCRdata, cutoff) {

  vacc_date             <- df$date_vaccine1[1]
  first_collection_date <- min(df$collected_date)

  # Inferred infections
  df <- df |>
    arrange(days_bleed) |>
    mutate(delta_t = as.numeric(days_bleed - lag(days_bleed))) |>
    mutate(titre_infection = ifelse(titre > lag(titre, default = Inf) + cutoff,
                                    1, 0)) |>
    # This handles measurements taken on the rise
    mutate(next_delta_t    = lead(delta_t, default = Inf),
           next_titre      = lead(titre),
           titre_infection = ifelse(
             next_delta_t < 60 &
               (next_titre > lag(titre, default = Inf) + cutoff), 1,
             titre_infection)) |>
    select(-next_delta_t, -next_titre) |>
    mutate(consecutive_infection = ifelse(titre_infection == 1 &
                                          lag(titre_infection == 1), 1, 0)) |>
    mutate(titre_infection = ifelse(titre_infection == 1 &
                                      consecutive_infection == 1 &
                                      delta_t < 180, 0, titre_infection)) |>
    select(-consecutive_infection)

  #-----------------------------------------------------------------------------

  min_days_bleed <- min(df$days_bleed)

  PCRdata <- PCRdata |> filter(dengue_days_pd1 >= min_days_bleed)

  df$PCR_infection <- 0
  df$serotype      <- "Subclinical"

  PCR_inf <- nrow(PCRdata)

  if(PCR_inf > 0)  {

    PCRdata <- PCRdata |> mutate(inf_date = vacc_date + dengue_days_pd1)

    df <-  df |> mutate(cd_previous = lag(collected_date,
                                          default = first_collection_date))

    for(i in 1:PCR_inf) {

      inf_date <- PCRdata$inf_date[i]

      # Position
      pos <- which(inf_date > df$cd_previous & inf_date <= df$collected_date)

      if(length(pos) == 0) next

      titre_inf_prev_pos <- df$titre_infection[pos - 1]

      # How many days ago was the previous measurement
      time_prev_meas <- df$delta_t[pos]

      # Handling measurements on the rise
      if(titre_inf_prev_pos == 1 & time_prev_meas < 180) pos <- pos - 1

      df$PCR_infection[pos] <- 1
      df$serotype[pos]      <- PCRdata$serotype[i]
    }
  }

  df <- df |>
    mutate(infection = ifelse(PCR_infection == 1 | titre_infection == 1, 1, 0))
}

remove_multiple_measurements_in_a_year <- function(df) {

  df <- df |>
    mutate(collected_year  = year(collected_date))

  df_list <- split(df, df$collected_year)

  map_df(df_list, \(df_year) {

    if(nrow(df_year) == 1) return(df_year)

    if(nrow(df_year) > 1) {

      infections <- pull(df_year, infection)

      if(sum(infections) > 1) {
        msg <- "More than one infection was detected in a year."
        stop(msg ,call. = FALSE)
      }

      if(sum(infections) == 0) {
        df_year <- df_year |> filter(collected_date == max(collected_date))
      }

      if(sum(infections) == 1) {
        pos_inf  <- which(infections == 1)
        inf_date <- df_year$collected_date[pos_inf]
        df_year <- df_year |> filter(collected_date == inf_date)
      }

      return(df_year)

    }

    msg <- "remove_multiple_measurements_in_a_year() can't handle this data.frame"

    stop(msg ,call. = FALSE)
  })
}
