simulate_two_cohorts <- function(iter_val, param_obj, cohort_df1, cohort_df2,
                                 breaks, labels, ll_val)
{
  sim_1 <- simulate_DENV_infections_cohort(
    lambda_serotype = param_obj$lambda_1,
    loss_rate       = param_obj$rho,
    cohort_df       = cohort_df1)

  age_inf_1 <- calculate_pct_infection_by_age(sim_1) |>
    rename(age = age_sample) |>
    group_infections_in_bins(breaks, labels) |>
    mutate(iter    = iter_val,
           log_lik = ll_val)

  sim_2 <- simulate_DENV_infections_cohort(
    lambda_serotype = param_obj$lambda_2,
    loss_rate       = param_obj$rho,
    cohort_df       = cohort_df2)

  age_inf_2 <- calculate_pct_infection_by_age(sim_2) |>
    rename(age = age_sample) |>
    group_infections_in_bins(breaks, labels) |>
    mutate(iter    = iter_val,
           log_lik = ll_val)

  list(cohort_1 = age_inf_1, cohort_2 = age_inf_2)
}
