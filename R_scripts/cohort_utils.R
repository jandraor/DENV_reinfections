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
  df <- df |> mutate(age_bin = cut(age,
                             breaks = breaks,
                             labels = labels,
                             right  = FALSE)) |>
    group_by(age_bin) |>
    summarise(n_infections  = sum(n_infections),
              n_symp        = sum(n_symp),
              n_individuals = sum(n_individuals)) |>
    mutate(pct_symp = n_symp / n_infections,
           pct_symp_given_age = n_symp / n_individuals)

  df
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

add_binomial_symp_CI <- function(df)
{
  df <- df |>
    rowwise() |>
    mutate(ci_inf = list(binom.confint(n_symp, n_individuals,
                                       methods = "wilson"))) |>
    mutate(lower = ci_inf$lower,
           upper = ci_inf$upper) |>
    select(-ci_inf) |>
    ungroup()

  df
}

add_bootstrap_CI <- function(df)
{
  df |>
    group_by(bin_delta) |>
    summarise(
      boot_obj = list(boot(delta_titre, boot_mean, R = 1000)),
      n        = n(), .groups = "drop") |>
    rowwise() |>
    mutate(
      ci    = list(boot.ci(boot_obj, type = "perc")$percent[4:5]),
      mean  = boot_obj$t0,
      q2.5  = ci[[1]],
      q97.5 = ci[[2]])
}

estimate_symptomatic_RR <- function(df_list, a_ref, n_ref)
{
  map_df(df_list, \(df) {

    if(nrow(df) == 0) return(NULL)

    CI_obj <- wald_rr(
      a     = df$n_symp,
      n     = df$n_individuals,
      a_ref = a_ref,
      n_ref = n_ref)

    if(df$n_symp == 0)
    {
      return(data.frame(age_bin = unique(df$age_bin),
                        q2.5    = NA,
                        RR      = 0,
                        q97.5   = NA))
    }

    data.frame(age_bin = unique(df$age_bin),
               q2.5    = CI_obj$lower,
               RR      = CI_obj$RR,
               q97.5   = CI_obj$upper)

  })
}

# a      = symptomatic cases in group
# n      = total infections in group
# a_ref  = symptomatic cases in reference group
# n_ref  = total infections in reference group
wald_rr <- function(a, n, a_ref, n_ref, alpha = 0.05)
{
  a_cc     <- a
  n_cc     <- n
  a_ref_cc <- a_ref
  n_ref_cc <- n_ref

  # Relative risk
  rr <- (a_cc / n_cc) / (a_ref_cc / n_ref_cc)

  # Standard error of log(RR)
  se_log_rr <- sqrt(
    (1 / a_cc - 1 / n_cc) +
      (1 / a_ref_cc - 1 / n_ref_cc)
  )

  z <- qnorm(1 - alpha / 2)

  # Confidence interval
  lower <- exp(log(rr) - z * se_log_rr)
  upper <- exp(log(rr) + z * se_log_rr)

  list(RR = rr, lower = lower, upper = upper)
}
