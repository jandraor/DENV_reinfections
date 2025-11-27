estimate_rise <- function(baseline_titre, rise_options) {
  max(0, rise_options$intercept + rise_options$slope * max(0, baseline_titre))
}

simulate_titres_updating_decay <- function(subject_df, stop_time,
                                           long_decay_rate)
{
  s_id      <- unique(subject_df$subject_id)

  inf_times <- round(subject_df$time)

  n_inf <- length(inf_times)

  titre_mat <- matrix(NA, nrow =  n_inf, ncol = stop_time + 1,
                      dimnames = list(NULL, 0:stop_time))

  for(inf_idx in seq_along(inf_times))
  {
    inf_time <- inf_times[[inf_idx]]

    if(inf_idx == 1) A0 <- inv_log2_transform(5)

    if(inf_idx > 1)
    {
      baseline_prev_inf <- titre_mat[1:(inf_idx - 1), inf_time + 1]

      baseline_log2 <- sum(log2_transform(baseline_prev_inf))

      # rise_options <- list(intercept = 10.33,
      #                      slope     = -0.93)

      rise_options <- list(intercept = 5.1,
                           slope     = -0.74)

      rise <- estimate_rise(baseline_log2, rise_options)

      A0 <- inv_log2_transform(rise)
    }

    start_period <- inf_time

    end_period <- stop_time

    sim_length <- end_period - start_period

    sim_titre <- A0 * exp(-long_decay_rate[[inf_idx]] * (0:sim_length))

    floor <- 5

    idx_from_now_on <- (inf_time + 1):(end_period + 1)

    titre_mat[inf_idx, idx_from_now_on] <- pmax(sim_titre, floor)

    if(inf_idx > 1)
    {
      for(j in 1:((inf_idx - 1)))
      {
        current_titre <- titre_mat[j, inf_time + 1]
        updated_titre <- current_titre * exp(-long_decay_rate[[inf_idx]] * (0:sim_length))

        titre_mat[j, idx_from_now_on] <- pmax(updated_titre, floor)
      }
    }
  }

  as.data.frame(titre_mat) |>
    mutate(inf_idx = seq_along(inf_times),
           .before = everything()) |>
    pivot_longer(-inf_idx, names_to = "day",
                 values_to = "titre") |>
    mutate(day        = as.numeric(day),
           subject_id = s_id)
}

simulate_titres_scenarios <- function(damp_df)
{
  lapply(1:3, \(sce_id) {

    fn_inf <- str_glue("./saved_objects/sim_infections/inf_{sce_id}.rds")
    inf_df <- readRDS(fn_inf)

    fn <- str_glue("./saved_objects/sim_titres/sce_{sce_id}.rds")

    if(!file.exists(fn))
    {
      inf_df <- inf_df |> mutate(chunk_id = (subject_id - 1) %/% 1000)

      chunk_list <- split(inf_df, inf_df$chunk_id)

      sim_pack <- map_df(chunk_list, simulate_titres_from_histories,
                         damp_df = damp_df)

      saveRDS(sim_pack, fn)
    } else sim_pack<- readRDS(fn)

    sim_pack
  })
}

simulate_titres_from_histories <- function(infection_history_df, damp_df) {

  df_list <- split(infection_history_df, infection_history_df$subject_id)

  meas_events <- 365 * seq(1:80)

  sce_lbl <- unique(infection_history_df$scenario_lbl)

  sim_decay_rate <- damp_df$decay_rate

  plan(multisession, workers = future::availableCores() - 1)

  titres_df <- future_map_dfr(
    df_list,
    simulate_titres_updating_decay,
    stop_time = 365 * 80,
    long_decay_rate = sim_decay_rate / 365) |>
    filter(day %in% meas_events) |>
    mutate(scenario_lbl = sce_lbl)

  titres_df
}
