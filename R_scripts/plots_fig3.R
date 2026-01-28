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
    theme(legend.position = "none")
}

plot_3B <- function(df, clr_vctr, linetype_vct)
{
  ggplot(df, aes(age, mean_titre)) +
    geom_line(aes(colour = as.factor(scenario),
                  linetype = as.factor(scenario)),
              linewidth =  1.2) +
    scale_colour_manual(values = clr_vctr) +
    scale_linetype_manual(values = c("solid", "11", "solid", "solid")) +
    scale_y_continuous(limits = c(0, NA)) +
    guides(colour   = guide_legend(direction = "vertical")) +
    labs(x = "Age", y = "Mean titre(log2)", colour = NULL,
         linetype = NULL) +
    theme(legend.position = "none")
}

plot_3C <- function(df, clr_vctr, linetype_vct)
{
  ggplot(df, aes(age, mean)) +
    geom_line(aes(group    = as.factor(scenario),
                  colour   = as.factor(scenario),
                  linetype = as.factor(scenario)),
              linewidth = 1.2) +
    labs(x = "Age", y = "Mean cumulative infections") +
    scale_colour_manual(values = clr_vctr) +
    scale_linetype_manual(values = linetype_vct) +
    theme_classic() +
    theme(legend.position = "none",
          axis.line = element_line(colour = "grey75"),
          axis.title = element_text(colour = "grey55"))
}

plot_3D <- function(df, n_ind, sce_name)
{
  frac_df <- df |>
    mutate(bin_inf = ifelse(n_inf_cum < 10, as.character(n_inf_cum), "10+"),
           bin_inf = factor(bin_inf, levels = c(0:9, "10+"))) |>
    group_by(age, bin_inf) |>
    summarise(frac = sum(count) / n_ind,
              .groups = "drop")

  my_colours <- c("grey85", viridis(10))

  ggplot(frac_df, aes(x = age, y = frac, fill = factor(bin_inf),
                 colour = factor(bin_inf))) +
    geom_area(position = "stack") +
    scale_fill_manual(values = my_colours, name = "Cumulative infections") +
    scale_colour_manual(values = my_colours, name = "Cumulative infections") +
    labs(x = "Age", y = "Proportion", subtitle = parse(text = sce_names)) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.2, "cm"),
          legend.text = element_text(size = 6),
          legend.margin = margin(t = -5, b = -2.5),
          legend.title = element_text(size = 7.5),
          legend.spacing.x = unit(0.01, "cm"))
}
