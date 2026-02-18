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

plot_S3_S4 <- function()
{
  summary_list <- NMC_get_measurements_summary()

  infection_detection_df <- NMC_get_infection_df(cut_off = 1.18) |>
    mutate(std_log2_D1 = log2_D1_NT - summary_list$means[["log2_D1_NT"]],
           std_log2_D2 = log2_D2_NT - summary_list$means[["log2_D2_NT"]],
           std_log2_D3 = log2_D3_NT - summary_list$means[["log2_D3_NT"]],
           std_log2_D4 = log2_D4_NT - summary_list$means[["log2_D4_NT"]],
           is_inf = ifelse(is.na(is_inf), 0, is_inf))

  inf_detection_all_df <- infection_detection_df |>
    select(subjectNo, days_bleed, is_inf, contains("log2_D"))

  infection_df <- inf_detection_all_df |> filter(is_inf == 1)

  ids <- infection_df |> group_by(subjectNo) |>
    summarise(n_inf = n()) |>
    filter(n_inf < 4) |>  # there are few individuals with >=4 infections
    pull(subjectNo)

  blood_after_df <- infection_df |> filter(subjectNo %in% ids) |>
    arrange(subjectNo, days_bleed) |>
    group_by(subjectNo) |>
    mutate(seq_inf = row_number()) |>
    select(subjectNo, seq_inf, contains("log2_D")) |>
    ungroup()

  blood_after_df$max_sero <- max.col(blood_after_df[, c(3:6)])

  blood_after_df$max_sero_std <- max.col(blood_after_df[, c(7:10)])

  df_list <- split(blood_after_df, blood_after_df$subjectNo)

  blood_after_df <- map_df(df_list, \(df) {

    df$n_distinct     <- NA
    df$n_distinct_std <- NA

    for(i in seq_len(nrow(df)))
    {
      df[i, "n_distinct"]     <- df$max_sero[1:i] |> unique() |> length()
      df[i, "n_distinct_std"] <- df$max_sero_std[1:i] |> unique() |> length()
    }

    df
  })

  boot_ci <- function(x, B = 1000) {
    boots <- replicate(B, mean(sample(x, replace = TRUE)))
    quantile(boots, c(0.025, 0.975))
  }

  summary_df <- blood_after_df |>
    group_by(seq_inf) |>
    summarise(
      mean      = mean(n_distinct),
      lower     = boot_ci(n_distinct)[1],
      upper     = boot_ci(n_distinct)[2],
      mean_std  = mean(n_distinct_std),
      lower_std = boot_ci(n_distinct_std)[1],
      upper_std = boot_ci(n_distinct_std)[2],
      .groups = "drop")

  fit <- lm(I(mean - 1) ~ 0 + I(seq_inf - 1), data = summary_df)

  p <- predict(fit, newdata = data.frame(seq_inf = c(1,2,3)))

  summary_df$pred_mean <- p + 1

  fit2 <- lm(I(mean_std - 1) ~ 0 + I(seq_inf - 1), data = summary_df)

  p2 <- predict(fit2, newdata = data.frame(seq_inf = c(1,2,3)))

  summary_df$pred_mean_std <- p2 + 1

  g_S2 <- ggplot(summary_df, aes(seq_inf)) +
    geom_line(aes(y = pred_mean, linetype = "Unadjusted"), colour = "grey60") +
    geom_line(aes(y = pred_mean_std, linetype = "Assay adjusted"),
              colour = NMC_all) +
    geom_pointrange(aes(y = mean_std, ymin = lower_std, ymax = upper_std),
                    colour = NMC_all) +
    scale_y_continuous(limits = c(0.91, 3),
                       breaks = 1:3) +
    scale_x_continuous(breaks = 1:3) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    geom_hline(yintercept = 1, linetype = "dotted", colour = "grey50") +
    annotate("text", label = "Original antigenic sin",
             x = 2, y = 0.925) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    annotate("text", x = 2, y = 2.125,
             label = "Independece across serotypes", angle = 45) +
    coord_equal(xlim = c(1, NA)) +
    labs(x = "Number of infection",
         y = "# of serotypes with the highest titer",
         linetype = NULL) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.3, 0.8),
      legend.key.width = unit(1.5, "cm")) +
    guides(linetype = guide_legend(override.aes = list(linewidth = 1.2)))

  inf_detect_seroneg_df <- infection_detection_df |>
    filter(serostatus == "seronegative")

  inf_seroneg_df <- inf_detect_seroneg_df |> filter(is_inf == 1)

  seroneg_ids <- inf_seroneg_df |> group_by(subjectNo) |>
    summarise(n_inf = n()) |>
    filter(n_inf < 4) |>  # there are few individuals with >=4 infections
    pull(subjectNo)

  bld_after_seroneg_df <- blood_after_df |>
    filter(subjectNo %in% seroneg_ids)

  summary_seroneg_df <- bld_after_seroneg_df |>
    group_by(seq_inf) |>
    summarise(
      mean_std  = mean(n_distinct_std),
      lower_std = boot_ci(n_distinct_std)[1],
      upper_std = boot_ci(n_distinct_std)[2],
      .groups = "drop")

  fit3 <- lm(I(mean_std - 1) ~ 0 + I(seq_inf - 1), data = summary_seroneg_df )

  p3 <- predict(fit3, newdata = data.frame(seq_inf = c(1,2,3)))

  summary_seroneg_df$pred_mean_std <- p3 + 1


  g_S3 <- ggplot(summary_df, aes(seq_inf)) +
    geom_line(aes(y = pred_mean_std, linetype = "All individuals"),
              colour = "grey50",
              position = position_nudge(x = 0.02)) +
    geom_line(data = summary_seroneg_df,
              aes(y = pred_mean_std, linetype = "Seronegative"),
              colour = NMC_all,
              position = position_nudge(x = -0.02)) +
    geom_pointrange(aes(y = mean_std, ymin = lower_std, ymax = upper_std),
                    colour = "grey50",
                    position = position_nudge(x = 0.02)) +
    geom_pointrange(data = summary_seroneg_df,
                    aes(y = mean_std, ymin = lower_std, ymax = upper_std),
                    colour = NMC_all,
                    position = position_nudge(x = -0.02)) +
    scale_y_continuous(limits = c(0.91, 3),
                       breaks = 1:3) +
    scale_x_continuous(breaks = 1:3) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    coord_equal(xlim = c(1, NA)) +
    labs(x = "Number of infection",
         y = "# of serotypes with the highest titer",
         linetype = NULL) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.3, 0.6),
      legend.key.width = unit(1.5, "cm")) +
    guides(linetype = guide_legend(override.aes = list(linewidth = 1.2)))

  list(S2 = g_S2,
       S3 = g_S3)
}

plot_S5 <- function()
{

  raw <- read_excel("./data/NMC/NMC Laboratory testing result to Henrik 23JUNE23.xlsx")

  Bleeds <- colnames(raw) |>
    grep(pattern = '^BL', value = T) |>
    str_extract('^BL[0-9A-Z\\-]+') |>
    unique()


  processTiters = function(x) {
    mult = str_extract(x, '>|<')
    ifelse( is.na(mult)
            , 1
            , c('<' = 1/2, '>' = 2)[mult]
    ) *
      as.numeric(gsub('^[^0-9]', '', x))
  }

  dat <-
    Bleeds |>
    lapply(function(bleed_name){
      out <-
        raw |>
        select(starts_with(bleed_name)) |>
        rename_all(function(x) gsub(paste0('^',bleed_name,' '), '', x)) |>
        pivot_longer(cols = c(-No,-CollectionDate),names_to = "assay", values_to = "titer") |>
        mutate(
          CollectionDate = as.Date(CollectionDate)
          , virus = str_extract(assay, '^[A-Z1-4]+')
          , assay = str_extract(assay, '[A-Z]+$')
          , titer = processTiters(titer)
        )
      out = split(out, out$assay)
      inner_join(out$HAI, out$PRNT,
                 by = c('No', 'virus', 'CollectionDate'),
                 suffix = c('.hai', '.prnt'))
    })


  # GMT across DENV titres
  dat.gmt <-
    dat |>
    lapply(function(x) {
      x |>
        filter(grepl('^D[1-4]', virus)) |>
        group_by(No) |>
        summarise(
          logGMT.hai  = mean(1 + log2(titer.hai / 10)),
          logGMT.prnt = mean(1 + log2(titer.prnt / 10)))
    }) |> do.call(what = rbind)

  # fit linear model to translate PRNT to HAI
  fit <- glm(logGMT.hai ~ logGMT.prnt, data = dat.gmt)
  summary(fit)

  cor.test(dat.gmt$logGMT.hai, dat.gmt$logGMT.prnt, method = 'pearson')

  dat.gmt |>
    ggplot(aes(x = logGMT.prnt, y = logGMT.hai)) +
    geom_point(size = 1, shape = 1, alpha = 0.2, colour = NMC_all) +
    geom_smooth(method = 'lm', colour = NMC_all, fill = NMC_all) +
    coord_fixed(ratio = 1) +
    # actual means
    geom_pointrange(data =
                      dat.gmt |>
                      mutate(
                        prntBucket = cut(logGMT.prnt,
                                         breaks = c(-Inf,0:10,Inf))) |>
                      group_by(prntBucket) |>
                      summarise(
                        N        = n(),
                        avg.prnt = mean(logGMT.prnt),
                        avg.hai  = mean(logGMT.hai),
                        sd.prnt  = sd(logGMT.prnt),
                        sd.hai   = sd(logGMT.hai)),
                    aes(x = avg.prnt, y = avg.hai,
                        ymin = avg.hai - 1.96 * sd.hai / sqrt(N),
                        ymax = avg.hai + 1.96 * sd.hai / sqrt(N)),
                    alpha = 0.5, colour = NMC_all) +
    labs(
      x = expression(PRNT~(log[2]^"*")),
      y = expression(HAI~(log[2]^"*")))
}

plot_S6 <- function()
{
  rise_df <- NMC_get_PCR_rise()

  ggplot(rise_df, aes(rise)) +
    geom_histogram(colour = "white",binwidth = 1, fill = NMC_PCR) +
    scale_x_continuous(n.breaks = 6) +
    labs(y = "Frequency",
         x = expression("Rise in titers (" * log[2]^"*" * " scale)"))
}

plot_S7 <- function()
{
  infection_detection_df <- NMC_get_infection_df(cut_off = 1.18)

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
         y = expression("Titer difference (" * log[2]^"*" * ")"),
         fill   = "",
         colour = "") +
    theme(legend.position = "inside",
          legend.position.inside = c(0.35, 0.35))
}

plot_S8 <- function()
{
  sens_df <- map_df(c(1, 1.18, 1.5), \(cutoff) {

    mean_df <- NMC_get_binned_decay(cutoff) |>
      mutate(cutoff_val = as.factor(cutoff))

    mean_df
  }) |> filter(n > 30)

  df_list <- split(sens_df, sens_df$cutoff_val)

  slope_df <- imap_dfr(df_list, \(df, co) {

    fit <- lm(mean ~ 0 + bin_delta, data = df)

    ci <- confint(fit, level = 0.95)

    data.frame(
      cut_off = co,
      slope   = coef(fit)[["bin_delta"]],
      lower   = ci["bin_delta", 1],
      upper = ci["bin_delta", 2])
  })

  ggplot(sens_df, aes(bin_delta, mean)) +
    geom_smooth(method = "lm", formula = y ~ 0 + x,
                aes(group = cutoff_val, linetype = cutoff_val),
                colour = NMC_all,
                fullrange = TRUE, se  = FALSE) +
    scale_linetype_manual(values = c("dotted", "solid", "dashed")) +
    expand_limits(x = 0, y = 0) +
    labs(x = "Years between blood draws",
         y = expression("Titer difference (" * log[2]^"*" * ")"),
         linetype = "Cutoff") +
    scale_x_continuous(limits = c(0, NA)) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.3, 0.3),
      legend.key.width = unit(1.5, "cm")) +
    guides(linetype = guide_legend(override.aes = list(linewidth = 1.2)))

}

plot_S9 <- function()
{
  df <- map_df(c(2, 5, 8), \(imputed_val) {

    KFCS_infection_detection_df <- KFCS_get_infection_df(imputed_val)

    KFCS_infection_df <- KFCS_infection_detection_df |>
      group_by(subjectNo) |>
      arrange(subjectNo, start_interval) |>
      mutate(inf_idx = cumsum(is_inf),
             inf_id  = paste(subjectNo, inf_idx, sep = "_")) |>
      ungroup()

    df_list <- split(KFCS_infection_df, KFCS_infection_df$inf_id)

    KFCS_delta_df <- map_df(df_list, KFCS_estimate_deltas) |>
      mutate(bin_delta    = round(delta_time, 0))

    KFCS_mean_estimates <- KFCS_delta_df |> add_bootstrap_CI() |>
      mutate(cohort = "KFCS",
             imputed_val = as.factor(imputed_val)) |>
      filter(bin_delta > 0)

    fit <- lm(mean ~ 0 + bin_delta, data = KFCS_mean_estimates)
    print(fit)

    print(round(confint(fit), 2))

    KFCS_mean_estimates
  })

  ggplot(df, aes(bin_delta, mean)) +
    geom_smooth(method = "lm", formula = y ~ 0 + x,
                aes(group = imputed_val,
                    linetype = imputed_val),
                alpha = 0.1,
                fullrange = TRUE,
                se = FALSE,
                colour = KFCS_clr) +
    expand_limits(x = 0, y = 0) +
    scale_linetype_manual(values = c("dotted", "solid", "dashed")) +
    labs(x = "Years between blood draws",
         y = expression("Titer difference (" * log[2]^"*" * ")"),
         linetype = "Imputed value") +
    scale_x_continuous(limits = c(0, NA)) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.3, 0.3),
      legend.key.width = unit(1.5, "cm")) +
    guides(linetype = guide_legend(override.aes = list(linewidth = 1.2)))
}
