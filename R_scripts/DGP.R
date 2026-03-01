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

estimate_avg_titre_by_age <- function(log_first_peak, decay_rate,
                                      phi, beta, inf_times_list, final_age)
{
  titre_mat <- simulate_DENV_long_decay_titres(
    inf_times_list = inf_times_list,
    log_first_peak = log_first_peak,
    decay_rate_vec = decay_rate,
    phi            = phi,
    beta           = beta,
    final_age      = final_age)

  titre_mat[titre_mat < 10] <- 5

  log2_titre_mat <- log2_transform(titre_mat)

  avg_titre <- colMeans(log2_titre_mat)

  avg_titre
}

get_decay_rate <- function() 0.2 * exp(-0.5*(0:3))
