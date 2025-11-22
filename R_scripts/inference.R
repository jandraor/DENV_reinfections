approximate_prob_inf_at_age <- function(lambda_serotype,
                                        loss_rate,
                                        n_individuals,
                                        end_year)
{
  set.seed(1853)
  inf_df <- simulate_DENV_infections_since_birth(lambda_serotype,
                                                 loss_rate,
                                                 n_individuals, end_year)

  inf_df |> group_by(age) |> count() |>
    mutate(pct = round(n / n_individuals, 2))
}

loglik_fun <- function(pars, data_df, final_age)
{
  lambda_val    <- inv.logit(pars[[1]])
  loss_rate_val <- inv.logit(pars[[2]])

  age_prob_df <- approximate_prob_inf_at_age(lambda_val,
                                             loss_rate_val,
                                             n_individuals = 1e5,
                                             final_age)
  -sum(dbinom(x    = data_df$n_infections,
              size = data_df$n_individuals,
              prob = age_prob_df$pct[-1],log = TRUE))

}

loglik_2_cohorts <- function(pars, data_df1, data_df2,
                             final_age_1, final_age_2,
                             n_particles)
{
  lambda_1      <- boot::inv.logit(pars[[1]])
  lambda_2      <- boot::inv.logit(pars[[2]])
  loss_rate_val <- boot::inv.logit(pars[[3]])

  age_prob_df1 <- approximate_prob_inf_at_age(lambda_1,
                                              loss_rate_val,
                                              n_individuals = n_particles,
                                              final_age_1)

  ll_1 <- -sum(dbinom(x    = data_df1$n_infections,
                      size = data_df1$n_individuals,
                      prob = age_prob_df1$pct[-1],
                      log = TRUE))

  age_prob_df2 <- approximate_prob_inf_at_age(lambda_2,
                                              loss_rate_val,
                                              n_individuals = n_particles,
                                              final_age_2)

  ll_2 <- -sum(dbinom(x   = data_df2$n_infections,
                     size = data_df2$n_individuals,
                     prob = age_prob_df2$pct[-1],
                     log = TRUE))

  ll_1 + ll_2
}

profile_ll_fun <- function(estimated_pars, data_df1, data_df2,
                           final_age_1, final_age_2,
                           n_particles, fixed_par, fixed_pos)
{
  pars <- append(estimated_pars, fixed_par, after = fixed_pos - 1)

  loglik_2_cohorts(pars,
                   data_df1,
                   data_df2,
                   final_age_1,
                   final_age_2,
                   n_particles)
}

maximise_over_profile <- function(par_vals, par_name, MLE_vals, fixed_pos,
                                  data_df1, data_df2, final_age_1, final_age_2)
{
  n_cores <- future::availableCores() - 1

  plan(multisession, workers = n_cores)

  with_progress(
  {
    p <- progressor(steps = length(par_vals))

    profile_df <- future_imap_dfr(par_vals, run_profile_optim,
                                  par_name    = par_name,
                                  MLE_vals    = MLE_vals,
                                  fixed_pos   = fixed_pos,
                                  data_df1    = data_df1,
                                  data_df2    = data_df2,
                                  final_age_1 = final_age_1,
                                  final_age_2 = final_age_2,
                                  p           = p)
  })

  plan(sequential)
  gc()

  profile_df
}

run_profile_optim <- function(par_value, iter, par_name, MLE_vals, fixed_pos,
                              data_df1, data_df2, final_age_1, final_age_2,
                              p = NULL)
{
  fn <- str_glue("./saved_objects/inference/profile_{par_name}/iter_{iter}.rds")

  if(!file.exists(fn))
  {
    estimated_pars <- logit(MLE_vals[-fixed_pos])

    inits <- append(MLE_vals[-fixed_pos], par_value,
                    after = fixed_pos - 1)

    optim_result <- optim(estimated_pars,
                          method      = "Nelder-Mead",
                          fn          = profile_ll_fun,
                          data_df1    = data_df1,
                          data_df2    = data_df2,
                          final_age_1 = final_age_1,
                          final_age_2 = final_age_2,
                          n_particles = 1e5,
                          fixed_par   = logit(par_value),
                          fixed_pos   = fixed_pos)

    estimates <- append(optim_result$par, logit(par_value),
                        after = fixed_pos - 1)

    df <- data.frame(
      lambda_1_init = inits[[1]],
      lambda_2_init = inits[[2]],
      rho_init      = inits[[3]],
      lambda_1_MLE  = inv.logit(estimates[[1]]),
      lambda_2_MLE  = inv.logit(estimates[[2]]),
      rho_MLE       = inv.logit(estimates[[3]]),
      ll            = -optim_result$value,
      convergence   = optim_result$convergence)

    saveRDS(df, fn)
  } else df <- readRDS(fn)

  if (!is.null(p)) p() # Progressor

  df
}


