plot_1A <- function(df)
{
  ggplot(df, aes(x = collected_year, y = n_inf, fill = serotype)) +
    geom_bar(stat = "identity") +
    labs(x = "Year", y = "Number of infections", subtitle = "NMC",
         fill = NULL) +
    scale_fill_manual(values = c(DENV_1_4, untyped, subclinical)) +
    theme(legend.position = c(0.85, 0.8),
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = 5, colour = "grey25"))
}

plot_1B <- function(df)
{
  ggplot(df, aes(x = collected_year, y = n_inf, fill = serotype)) +
    geom_bar(stat = "identity") +
    scale_x_continuous(breaks = 2013:2014) +
    labs(x = "Year", y = NULL, subtitle = "CPC",
         fill = NULL) +
    scale_fill_manual(values = c("DENV-1" = DENV_1_4[1],
                                 "DENV-2" = DENV_1_4[2],
                                 "DENV-3" = DENV_1_4[3],
                                 "DENV-4" = DENV_1_4[4],
                                 "Subclinical" = subclinical)) +
    theme(legend.position = "none")
}

plot_1C <- function(df)
{
  ggplot(df, aes(x = collected_year, y = n_inf, fill = serotype)) +
    geom_bar(stat = "identity") +
    labs(x = "Year", y = NULL, subtitle = "KFCS",
         fill = NULL) +
    scale_fill_manual(values = c(DENV_1_4, untyped, subclinical)) +
    theme(legend.position = "none",
          legend.text = element_text(size = 8, colour = "grey25"))
}

plot_1D <- function(sim_df, data_df)
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
    scale_colour_manual(values = c("CPC" = CPC_clr, "NMC" = NMC_all)) +
    scale_fill_manual(values = c("CPC" = CPC_clr, "NMC" = NMC_all)) +
    scale_x_discrete(labels = labels) +
    scale_y_continuous(limits = c(0, NA), labels = label_percent()) +
    labs(x      = "Age group", y = "Infections",
         colour = "Cohort", subtitle = "CPC/NMC") +
    theme(axis.line = element_line(colour = "grey75"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.caption = element_text(colour = "grey65"),
          legend.position = "none")

}

plot_1E <- function(sim_df, data_df)
{
  ggplot(sim_df, aes(age_bin)) +
    geom_line(aes(y = q50, colour = cohort, group = cohort), alpha = 0.75,
              linetype = "11") +
    geom_errorbar(data = data_df,
                  aes(ymin = lower, ymax = upper,  colour = cohort),
                  width = 0.1) +
    geom_point(data = data_df, aes(y = pct, colour = cohort),
               size = 0.5) +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = cohort,
                    group = cohort), alpha = 0.25) +
    scale_colour_manual(values = c("KFCS" = KFCS_clr)) +
    scale_fill_manual(values = KFCS_clr) +
    scale_x_discrete(labels = labels) +
    scale_y_continuous(limits = c(0, NA), labels = label_percent()) +
    labs(x      = "Age group", y = "Infections",
         colour = "Cohort", subtitle = "KFCS") +
    theme(axis.line = element_line(colour = "grey75"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.caption = element_text(colour = "grey65"),
          legend.position = "none")

}

plot_1F <- function(df)
{
  df <- df |>
    mutate(undef = ifelse(RR == 0, TRUE, FALSE),
           RR_plot = ifelse(undef, 0.02, RR)) # 0.02 arbitrary for plotting

  ggplot(df, aes(age_bin, RR_plot)) +
    geom_point(aes(color = cohort, shape = undef),
               position = position_dodge(width = 0.6)) +
    scale_y_log10() +
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5, group = cohort,
                      colour = cohort),
                  position = position_dodge(width = 0.6),
                  width = 0) +
    scale_x_discrete(labels = labels) +
    scale_colour_manual(values = c("KFCS" = KFCS_clr, "CPC" = CPC_clr,
                                   "NMC" = NMC_all)) +

    scale_shape_manual(values = c(`FALSE` = 16, `TRUE`  = 1)) +
    labs(x = "Age group",
         y = "Relative risk (vs 0–9)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     size = 6),
          legend.position = "none")

}
