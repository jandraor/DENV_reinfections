plot_profile <- function(df, cutoff_value, x_pos, y_pos, ci_vals, par_name)
{
  ci_lbl <- str_glue("95% CI: {ci_vals[[1]]} - {ci_vals[[2]]}")

  ggplot(df, aes(value, ll)) +
    geom_point(colour = "steelblue") +
    geom_hline(yintercept = cutoff_value, linetype = "dashed") +
    annotate("text", x = x_pos, y = y_pos, label = ci_lbl, colour = "steelblue",
             size = 2.5) +
    labs(title = parse(text = par_name),
         y = "Log-lik", x = "Parameter value") +
    geom_smooth(se = FALSE)
}
