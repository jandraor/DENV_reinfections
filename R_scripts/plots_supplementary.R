plot_S1 <- function()
{
  pairs_df <- NMC_get_annual_pairs()

  pairs_df$result <- factor(pairs_df$result,
                            levels = c("No infection", "PCR Conf", "Subclinical"))

  overall_mean <- mean(pairs_df$titre_diff)

  decay_mean <- pairs_df |> filter(titre_diff < 0) |> pull(titre_diff) |> mean()

  PCR_mean <- pairs_df |> filter(result == "PCR Conf") |> pull(titre_diff) |>
    mean()

  subcl_mean <- pairs_df |> filter(result == "Subclinical") |>
    pull(titre_diff) |> mean()

  ggplot(pairs_df, aes(titre_diff, fill = result)) +
    geom_histogram(colour = "white",
                   position = position_stack(reverse = TRUE)) +
    geom_vline(xintercept = 1.18) +
    geom_vline(xintercept = overall_mean, linetype = "dashed") +
    geom_vline(xintercept = decay_mean, linetype = "dotdash") +
    geom_vline(xintercept = PCR_mean, linetype = "dotted",
               colour = "#E67700") +
    geom_vline(xintercept = subcl_mean, linetype = "dotted",
               colour = "#506d7E") +
    scale_fill_manual(values = c("No infection" = "#B3BEC6",
                                 "Subclinical"  = "#506d7E",
                                 "PCR Conf"     = "#E67700")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x    = "Titer difference between sequential annual blood draws",
         y    = "Number of pairs",
         fill = "") +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.8,0.8))
}

plot_S3 <- function()
{
  infection_detection_df <- NMC_get_infection_df()

  plac_ids <- unique(infection_detection_df$subjectNo)

  plac_inf <- NMC_get_symptomatic_infections(plac_ids)

  infection_df <- infection_detection_df |>
    mutate(is_inf = ifelse(is.na(is_inf), 0, is_inf)) |>
    group_by(subjectNo) |>
    mutate(inf_idx = cumsum(is_inf),
           inf_id  = paste(subjectNo, inf_idx, sep = "_")) |>
    ungroup()

  df_list <- split(infection_df, infection_df$inf_id)

  delta_df <- map_df(df_list, \(df) {

    df <- filter_short_term_dynamics(df, plac_inf, 365)

    if(nrow(df) <= 1) return(NULL)

    pairwise_delta_by_serotype(df)
  }) |> mutate(delta_year = delta_time / 365,
               bin_delta  = round(delta_year, 0)) |>
    select(inf_id, contains("delta"), -delta_time, -delta_year) |>
    pivot_longer(c(-inf_id, -bin_delta),
                 values_to = "delta_titre",
                 names_to = "serotype") |>
    mutate(serotype = str_replace(serotype, "delta_titre_D", "DENV")) |>
    filter(!is.na(delta_titre))

  df_list <- split(delta_df, delta_df$serotype)

  delta_summary_df <- map_df(df_list, \(df) {

    summary_df <- df |> add_bootstrap_CI() |>
      mutate(serotype = unique(df$serotype)) |>
      filter(n > 30)

    fit <- lm(mean ~ 0 + bin_delta, data = summary_df)

    slope <- fit$coefficients[[1]]

    summary_df <- summary_df |> mutate(slope = round(slope, 2))
  })

  ggplot(delta_summary_df, aes(bin_delta, mean)) +
    geom_smooth(method = "lm", formula = y ~ 0 + x,
                aes(group = serotype, fill = serotype, colour = serotype),
                fullrange = TRUE,
                se = FALSE) +
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5, x = bin_delta,
                      colour = serotype),
                  width = 0) +
    geom_point(aes(colour = serotype, x = bin_delta, y = mean),
               shape = 15, size = 2) +
   expand_limits(x = 0, y = 0) +
   scale_x_continuous(limits = c(0, 8)) +
   scale_y_continuous(limits = c(-2, NA)) +
    scale_colour_manual(values = c("DENV1" = DENV_1_4[[1]],
                                   "DENV2" = DENV_1_4[[2]],
                                   "DENV3" = DENV_1_4[[3]],
                                   "DENV4" = DENV_1_4[[4]])) +
    scale_fill_manual(values = c("DENV1" = DENV_1_4[[1]],
                                 "DENV2" = DENV_1_4[[2]],
                                 "DENV3" = DENV_1_4[[3]],
                                 "DENV4" = DENV_1_4[[4]])) +
    labs(x = "Years between blood draws",
         y = "Titre difference (log2)",
         fill   = "",
         colour = "") +
    theme(legend.position = "inside",
          legend.position.inside = c(0.35, 0.35))
}
