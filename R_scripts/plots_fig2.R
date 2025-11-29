# df: Symptomatic and subclinical infections have been removed
# df2: Only symptomatic infection have been removed
plot_2C <- function(df, df2, NMC_all, NMC_PCR)
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
    geom_smooth(data = df2 |> filter(bin_delta < 9),
                method = "lm", formula = y ~ 0 + x,
                aes(linetype = "Undetected", group = cohort, colour = cohort),
                alpha = 0.1,
                fullrange = TRUE, se = FALSE) +
    scale_colour_manual(values = c(KFCS_clr, NMC_all)) +
    labs(x = "Years between blood draws",
         y = "Titre difference (log2)",
         linetype = "Subclinical infections") +
    theme_classic() +
    scale_x_continuous(limits = c(0, NA)) +
    theme(axis.line = element_line(colour = "grey75"),
          legend.position = "none",
          legend.title = element_text(size = 8, colour = "grey65"),
          legend.text  = element_text(size = 7, colour = "grey60"),
          legend.margin = margin(t = 2, r = 2, b = 2, l = 2))
}

plot_2D <- function(df, KFCS_clr, NMC_clr)
{
  ggplot(df, aes(age, MA)) +
    geom_errorbar(aes(ymin = lower, ymax = upper, colour = cohort),
                  width = 0, alpha = 0.25) +
    geom_point(shape = 20, aes(colour =  cohort),
               size = 0.1, show.legend = FALSE) +
    labs(x = "Age", y = "Decay rate", colour = "Cohort") +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_colour_manual(values = c("KFCS" = KFCS_clr, "NMC" = NMC_clr)) +
    theme_classic() +
    theme(axis.line = element_line(colour = "grey75"),
          legend.position = "none",
          legend.direction = "horizontal")
}

plot_2E <- function(CPC_titre, NMC_titre, KFCS_titre,
                    CPC_avg_titre, NMC_avg_titre, KFCS_avg_titre,
                    CPC_clr, NMC_clr, KFCS_clr)
{
  ggplot(CPC_titre, aes(age)) +
    geom_point(shape = 1, aes(y = log2_mean, colour = "CPC"), alpha = 0.1) +
    geom_point(data = NMC_titre, aes(x = age, y = log2_mean * 0.72,
                                   colour = "NMC"), shape = 1,
               alpha = 0.02) +
    geom_point(data = KFCS_titre_df, aes(x = age, y = log2_mean, colour = "KFCS"),
               shape = 1, alpha = 0.02) +
    geom_line(data = CPC_avg_titre,
              aes(x = age, y = mean, colour = "CPC"), linewidth = 1.5) +
    geom_line(data = NMC_avg_titre, aes(x = age, y = mean * 0.72,
                                        colour = "NMC"), linewidth = 1.5) +
    geom_line(data = KFCS_avg_titre,
              aes(x = age, y = mean, colour = "KFCS"), linewidth = 1.5) +
    labs(y = "Titre (log2)", x = "Age", colour = "Cohort") +
    scale_colour_manual(values = c(CPC_clr, KFCS_clr, NMC_all)) +
    theme_classic() +
    theme(axis.line = element_line(colour = "grey75"),
          legend.position = "bottom",
          legend.direction = "horizontal")
}

plot_2F_1 <- function(sim_df, data_df, CPC_clr, KFCS_clr)
{
  ggplot(sim_df, aes(age_bin)) +
    geom_line(aes(y = q50, colour = cohort, group = cohort), alpha = 0.75,
              linetype = "11") +
    geom_errorbar(data = data_df,
                  aes(ymin = lower, ymax = upper,  colour = cohort),
                  width = 0.1, position = position_nudge(
                    x = ifelse(data_df$cohort == "NMC", 0.2, 0)
                  )) +
    geom_point(data = data_df, aes(y = pct, colour = cohort),
               size = 0.5, position = position_nudge(
                 x = ifelse(data_df$cohort == "NMC", 0.2, 0)
               )) +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = cohort,
                    group = cohort), alpha = 0.25) +
    facet_wrap(~country, nrow = 1) +
    scale_colour_manual(values = c(CPC_clr, KFCS_clr, NMC_all)) +
    scale_fill_manual(values = c(CPC_clr, KFCS_clr)) +
    scale_x_discrete(labels = labels) +
    scale_y_continuous(limits = c(0, NA), labels = label_percent()) +
    labs(x      = "Age group", y = "Infections",
         colour = "Cohort") +
    theme_classic() +
    theme(axis.line = element_line(colour = "grey75"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     size = 6),
          plot.caption = element_text(colour = "grey65"),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title = element_text(colour = "grey50"))

}

plot_2F_2 <- function(df, extra_cohort)
{
  ggplot(df, aes(x = age_bin, y = value)) +
    geom_bar(
      stat = "identity",
      position = "stack",
      aes(fill = key)) +
    scale_y_continuous(labels = percent_format()) +
    scale_x_discrete(labels = labels) +
    scale_fill_manual(values = c("#F1EADB", CPC_clr, "#F1EADB", KFCS_clr)) +
    geom_point(data = extra_cohort, aes(age_bin, symptomatic),
               colour = NMC_all) +
    labs(x = "Age group", y = "Symptomatic", fill = NULL) +
    guides(fill = guide_legend(keywidth = 0.3, keyheight = 0.3)) +
    facet_wrap(~country, nrow = 1) +
    theme_classic() +
    theme(axis.line = element_line(colour =  "grey65"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     size = 6),
          axis.text.y = element_text(size = 8),
          axis.title = element_text(size = 10, colour = "grey50"),
          legend.position = "none",
          legend.direction = "vertical",
          legend.text = element_text(size = 5),
          legend.margin = margin(0, 0, 0, 0),
          strip.background = element_blank(),
          strip.text = element_blank())
}
