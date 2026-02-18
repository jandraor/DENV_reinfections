clrs_1A_to_C <- c("DENV-1"        = DENV_1_4[1],
                  "DENV-2"        = DENV_1_4[2],
                  "DENV-3"        = DENV_1_4[3],
                  "DENV-4"        = DENV_1_4[4],
                  "Untyped/Mixed" = untyped,
                  "Subclinical"   = subclinical)

clrs_1D_to_F <- c("CPC" = CPC_clr, "NMC" = NMC_all, "KFCS" = KFCS_clr)

leg_text_size <- 14

plot_1A <- function(df)
{
  ggplot(df, aes(x = collected_year, y = n_inf, fill = serotype)) +
    geom_bar(stat = "identity", colour = NA) +
    labs(x = "Year", y = "Number of infections", subtitle = "NMC",
         fill = NULL) +
    scale_fill_manual(values = clrs_1A_to_C,
                      guide  = guide_legend(nrow = 1)) +
    scale_y_continuous(limits = c(0, NA)) +
    theme(legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = leg_text_size, colour = "grey25"),
          legend.position = "bottom")
}

plot_1B <- function(df)
{
  # For CPC, samples were collected in the first trimester of every year.
  #  Therefore, it is more likely than infections occurred in the year before.
  ggplot(df, aes(x = collected_year - 1, y = n_inf, fill = serotype)) +
    geom_bar(stat = "identity") +
    scale_x_continuous(breaks = 2012:2013) +
    labs(x = "Year", y = NULL, subtitle = "CPC",
         fill = NULL) +
    scale_y_continuous(breaks = c(0, 70, 140)) +
    scale_fill_manual(values = clrs_1A_to_C, guide  = guide_legend(nrow = 1)) +
    theme(legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = leg_text_size, colour = "grey25"),
          legend.position = "bottom")
}

plot_1C <- function(df)
{
  ggplot(df, aes(x = collected_year, y = n_inf, fill = serotype)) +
    geom_bar(stat = "identity") +
    labs(x = "Year", y = NULL, subtitle = "KFCS",
         fill = NULL) +
    scale_fill_manual(values = clrs_1A_to_C,
                      guide  = guide_legend(nrow = 1)) +
    theme(legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = leg_text_size, colour = "grey25"),
          legend.position = "bottom")
}

plot_1D <- function(sim_df, data_df)
{
  ggplot(sim_df, aes(age_bin)) +
    geom_line(aes(y = q50, colour = cohort, group = cohort), alpha = 0.75,
              linetype = "11", show.legend = FALSE) +
    geom_errorbar(data = data_df,
                  aes(ymin = lower, ymax = upper,  colour = cohort),
                  width = 0.1, position = position_nudge(
                    x = ifelse(data_df$cohort == "NMC", 0.2, 0)),
                  show.legend = FALSE) +
    geom_point(data = data_df, aes(y = pct, colour = cohort),
               size = 0.5, position = position_nudge(
                 x = ifelse(data_df$cohort == "NMC", 0.2, 0)
               )) +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = cohort,
                    group = cohort), alpha = 0.25, show.legend = FALSE) +
    scale_colour_manual(values = clrs_1D_to_F) +
    scale_fill_manual(values = clrs_1D_to_F) +
    scale_x_discrete(labels = labels) +
    scale_y_continuous(limits = c(0, 0.3)) +
    labs(x      = "Age group", y = "Probability of infection",
         colour = "Cohort", subtitle = "CPC/NMC") +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.caption = element_text(colour = "grey65"),
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = leg_text_size, colour = "grey25"),
          legend.position = "bottom")

}

plot_1E <- function(sim_df, data_df)
{
  ggplot(sim_df, aes(age_bin)) +
    geom_line(aes(y = q50, colour = cohort, group = cohort), alpha = 0.75,
              linetype = "11", show.legend = FALSE) +
    geom_errorbar(data = data_df,
                  aes(ymin = lower, ymax = upper,  colour = cohort),
                  width = 0.1, show.legend = FALSE) +
    geom_point(data = data_df, aes(y = pct, colour = cohort),
               size = 0.5) +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = cohort,
                    group = cohort), alpha = 0.25, show.legend = FALSE) +
    scale_colour_manual(values = clrs_1D_to_F) +
    scale_fill_manual(values = KFCS_clr) +
    scale_x_discrete(labels = labels) +
    scale_y_continuous(limits = c(0, 0.3)) +
    labs(x      = "Age group", y = "Probability of infection",
         colour = "Cohort", subtitle = "KFCS") +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    theme(axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.caption    = element_text(colour = "grey65"),
          legend.key.size = unit(0.4, "cm"),
          legend.text     = element_text(size = leg_text_size,
                                         colour = "grey25"),
          legend.position = "bottom")

}

plot_1F <- function(df)
{
  ggplot(prob_symp_age_df, aes(age_bin, pct_symp_given_age,
                               colour = cohort, group = cohort)) +
    geom_linerange(aes(ymin = 0, ymax = upper),
                   position = position_dodge(width = 0.6),
                   linewidth = 0.4, alpha = 0.25, show.legend = FALSE) +
    geom_point(position = position_dodge(width = 0.6), size = 1.5) +
    scale_y_continuous(trans = "log1p" ) +
    scale_colour_manual(values = clrs_1D_to_F) +
    labs(x = "Age group", y = "Prob symptomatic infection",
         subtitle = "All cohorts",
         colour = "Cohort") +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.key.size = unit(0.4, "cm"),
      legend.text     = element_text(size = leg_text_size,
                                     colour = "grey25"),
      legend.position = "bottom")
}
