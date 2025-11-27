calculate_pct_infection_by_age <- function(df)
{
  df |> group_by(age_sample) |>
    summarise(n_infections  = sum(infection),
              n_individuals = n()) |>
    mutate(pct = n_infections / n_individuals)
}

group_infections_in_bins <- function(df, breaks, labels)
{
  df |> mutate(age_bin = cut(age,
                             breaks = breaks,
                             levels = levels)) |>
    group_by(age_bin) |>
    summarise(n_infections  = sum(n_infections),
              n_individuals = sum(n_individuals)) |>
    mutate(pct = n_infections / n_individuals)

}

group_symp_infection_in_bins <- function(df, breaks, labels)
{
  df |> mutate(age_bin = cut(age,
                             breaks = breaks,
                             levels = levels)) |>
    group_by(age_bin) |>
    summarise(n_infections = sum(n_infections),
              n_symp       = sum(n_symp)) |>
    mutate(pct_symp = n_symp / n_infections)


}

add_binomial_CI <- function(df)
{
  df <- df |>
    rowwise() |>
    mutate(ci_inf = list(binom.confint(n_infections, n_individuals,
                                       methods = "wilson"))) |>
    mutate(lower = ci_inf$lower,
           upper = ci_inf$upper) |>
    select(-ci_inf) |>
    ungroup()

  df
}

add_bootstrap_CI <- function(df, cohort_name)
{
  df |>
    group_by(bin_delta) |>
    summarise(
      boot_obj = list(boot(delta_titre, boot_mean, R = 1000)),
      .groups = "drop"
    ) |>
    rowwise() |>
    mutate(
      ci    = list(boot.ci(boot_obj, type = "perc")$percent[4:5]),
      mean  = boot_obj$t0,
      q2.5  = ci[[1]],
      q97.5 = ci[[2]]) |>
    mutate(cohort = cohort_name)
}


