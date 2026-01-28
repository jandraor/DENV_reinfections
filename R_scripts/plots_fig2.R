plot_2A <- function(df)
{
  ggplot(df, aes(collected_date, log2_mean)) +
    geom_line(aes(group = subjectNo), colour = NMC_all, alpha = 0.1) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_date(date_breaks = "3 years",
                 date_labels = "%Y") +
    geom_smooth(se = FALSE, colour = NMC_all) +
    labs(x = "Year", y = "Titre (log2)",
         subtitle = "NMC")
}

plot_2B <- function(df, df2)
{
  ggplot(df |> filter(!is.na(serostatus)), aes(baseline, rise)) +
    geom_point(colour = NMC_PCR, aes(shape = serostatus)) +
    labs(x     = "Baseline (log2)", y = "Rise (log2)",
         shape = NULL) +
    scale_y_continuous(limits = c(0, NA)) +
    geom_line(data = df2, aes(x = baseline, y = rise), colour = "grey50") +
    theme(legend.text = element_text(color = "grey55"),
          legend.position = c(0.4, 0.25),
          legend.direction = "vertical",
          legend.margin = margin(t = 3, r = 3, b = 3, l = 3),
          plot.caption = element_text(colour = "grey50"),
          legend.background = element_rect(colour = "grey80", linewidth = 0.1))
}

plot_2C <- function(df)
{
  ggplot(df, aes(years_post_peak, titre)) +
    geom_line(aes(group = inf_id), colour = NMC_PCR, alpha = 0.5) +
    geom_smooth(se = FALSE, colour = NMC_PCR, linewidth = 1.5) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = "Years post-infection", y = "Titre (log2)")
}

# df: Symptomatic and subclinical infections have been removed
# df2: Only symptomatic infection have been removed
plot_2D <- function(df)
{
  df <- df |> filter(bin_delta < 9)

  ggplot(df, aes(bin_delta, mean)) +
    geom_smooth(method = "lm", formula = y ~ 0 + x,
                aes(linetype = "Detected", group = cohort, fill = cohort,
                    colour = cohort),
                alpha = 0.1,
                fullrange = TRUE) +
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5, colour = cohort),
                  width = 0,
                  position = position_nudge(
                    x = ifelse(df$cohort == "KFCS", 0.1, 0))) +
    geom_point(shape = 15,
               aes(colour = cohort),
               position = position_nudge(x = ifelse(df$cohort == "KFCS",
                                                    0.1, 0)),
               size = 0.5) +
    expand_limits(x = 0, y = 0) +
    scale_colour_manual(values = c(KFCS_clr, NMC_all)) +
    labs(x = "Years between blood draws",
         y = "Titre difference (log2)",
         linetype = "Subclinical infections") +
    scale_x_continuous(limits = c(0, NA)) +
    theme(legend.position = "none")
}

plot_2E <- function(df, sim_df)
{
  ggplot(df, aes(age, MA)) +
    geom_errorbar(aes(ymin = lower, ymax = upper, colour = cohort),
                  width = 0, alpha = 0.25) +
    geom_point(shape = 20, aes(colour =  cohort),
               size = 0.1, show.legend = FALSE) +
    geom_line(data = sim_df, aes(y = mean_decay, group = cohort,
                                 colour = cohort), linetype = "dashed") +
    labs(x = "Age", y = "Decay rate", colour = "Cohort") +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_colour_manual(values = c("KFCS" = KFCS_clr, "NMC" = NMC_all)) +
    theme(legend.position  = "none",
          legend.direction = "horizontal")
}

plot_2F <- function(pred_df, data_df)
{
  ggplot(pred_df, aes(age, q50)) +
    geom_line(linewidth = 0.5, aes(colour = cohort)) +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = cohort), alpha = 0.5) +
    geom_point(data = data_df, aes(y = mean, colour = cohort), shape = 1,
               size = 0.5) +
    scale_colour_manual(values = c("NMC" = NMC_all, "KFCS" = KFCS_clr,
                                   "CPC" = CPC_clr)) +
    scale_fill_manual(values = c("NMC" = NMC_all, "KFCS" = KFCS_clr,
                                   "CPC" = CPC_clr)) +
    scale_y_continuous(limits = c(NA, 12)) +
    scale_x_continuous(limits = c(NA, 95)) +
    labs(x = "Age", y = "Mean titre (log2)") +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank())
}

