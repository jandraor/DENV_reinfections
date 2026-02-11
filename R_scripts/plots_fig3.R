plot_3A <- function(df, label_df, clrs, linetype_vctr)
{
  ggplot(df, aes(age, pct)) +
    geom_line(aes(colour = as.factor(lbl), group = as.factor(lbl),
                  linetype = as.factor(lbl)), linewidth =  1.2) +
    geom_text(data = label_df, aes(label = label, x = x, y = y,
                                  colour = as.factor(sce_id)),
              fontface = "bold") +
    scale_y_continuous(limits = c(0, NA)) +
    labs(y = "Probability of infection", x = "Age", colour = NULL,
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
    labs(x = "Age", y = "Mean titer(log2)", colour = NULL,
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
    theme(legend.position = "none",
          axis.line = element_line(colour = "grey75"),
          axis.title = element_text(colour = "grey55"))
}

plot_3D <- function(df, n_ind, sce_name)
{
  frac_df <- df |>
    mutate(bin_inf = ifelse(n_inf_cum < 9, as.character(n_inf_cum), "9+"),
           bin_inf = factor(bin_inf, levels = c(0:8, "9+"))) |>
    group_by(age, bin_inf) |>
    summarise(frac = sum(count) / n_ind,
              .groups = "drop")

  label_df <- data.frame(x = c(4, 9.5, 15, 22, 33, 46, 57, 65, 71, 77),
                         y = c(0.95, 0.85, 0.75, 0.65, 0.55, 0.43, 0.33, 0.22, 0.13, 0.05),
                         label = unique(frac_df$bin_inf))

  my_colours <- c("grey85", viridis(10))

  ggplot(frac_df, aes(x = age, y = frac)) +
    geom_area(position = "stack", aes(colour = factor(bin_inf),
                                      fill = factor(bin_inf))) +
    scale_fill_manual(values = my_colours, name = "Cumulative infections") +
    scale_colour_manual(values = my_colours, name = "Cumulative infections") +
    geom_label(
      data = label_df,
      aes(x = x, y = y, label = label),
      colour = "black",
      label.size = 0,
      fill = alpha("white", 0.05),
      size = 4,
      fontface = "bold",
      show.legend = FALSE) +
    labs(x = "Age", y = "Proportion", subtitle = sce_name) +
    theme(legend.position = "none",
          legend.key.size = unit(0.2, "cm"),
          legend.text = element_text(size = 6),
          legend.margin = margin(t = -5, b = -2.5),
          legend.title = element_text(size = 7.5),
          legend.spacing.x = unit(0.01, "cm"))
}
