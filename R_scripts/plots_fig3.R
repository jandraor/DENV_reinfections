plot_3A <- function(df, label_df, clrs, linetype_vctr)
{
  ggplot(df, aes(age, pct)) +
    geom_line(aes(colour = as.factor(lbl), group = as.factor(lbl),
                  linetype = as.factor(lbl)), linewidth =  1.2) +
    geom_text(data = label_df, aes(label = label, x = x, y = y,
                                  colour = as.factor(sce_id)),
              parse = TRUE) +
    scale_y_continuous(limits = c(0, NA), labels = label_percent()) +
    labs(y = "% infections", x = "Age", colour = NULL,
         linetype = NULL) +
    scale_colour_manual(values = clrs) +
    scale_linetype_manual(values = linetype_vctr) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.line = element_line(colour = "grey65"),
          axis.title = element_text(colour = "grey55"),
          legend.position = "none")
}

plot_3B <- function(df, clr_vctr, linetype_vct)
{
  ggplot(df, aes(age, meas)) +
    geom_line(aes(age, mean, colour = scenario_lbl, linetype = scenario_lbl),
              linewidth =  1.2) +
    scale_colour_manual(values = c(clr_dgm, clr_reinf, clr_reinf, "#C7CEAC")) +
    scale_linetype_manual(values = c("solid", "11", "solid", "solid")) +
    scale_y_continuous(limits = c(0, NA)) +
    guides(colour   = guide_legend(direction = "vertical")) +
    theme_classic() +
    labs(x = "Age", y = "Mean titre(log2)", colour = NULL,
         linetype = NULL) +
    theme(axis.line = element_line(colour = "grey65"),
          axis.title = element_text(size = 10, colour = "grey55"),
          legend.position = "none")
}

plot_3C <- function(df, clr_vctr, linetype_vct)
{
  ggplot(df, aes(age, mean)) +
    geom_line(aes(group    = as.factor(scenario_lbl),
                  colour   = as.factor(scenario_lbl),
                  linetype = as.factor(scenario_lbl)),
              linewidth = 1.2) +
    labs(x = "Age", y = "Mean cumulative infections") +
    scale_colour_manual(values = clr_vctr) +
    scale_linetype_manual(values = linetype_vct) +
    theme_classic() +
    theme(legend.position = "none",
          axis.line = element_line(colour = "grey75"),
          axis.title = element_text(colour = "grey55"))
}

plot_3D <- function(df)
{
  my_colours <- c("grey85", viridis(10))

  ggplot(df, aes(x = age, y = frac, fill = factor(n_inf_group),
                 colour = factor(n_inf_group))) +
    geom_area(position = "stack") +
    scale_fill_manual(values = my_colours, name = "Cumulative infections") +
    scale_colour_manual(values = my_colours, name = "Cumulative infections") +
    labs(x = "Age", y = "Density", subtitle = parse(text = sce_names[[1]])) +
    scale_y_continuous(labels = label_percent()) +
    facet_wrap(~lbl, labeller = label_parsed) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_classic() +
    theme(legend.position = "bottom",
          axis.line = element_line(colour = "grey70"),
          legend.key.size = unit(0.25, "cm"),
          legend.margin = margin(t = -5, b = -2.5),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title = element_text(colour = "grey55"))
}
