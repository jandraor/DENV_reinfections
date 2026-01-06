profile_fun_2 <- function(estimated_pars, titre_data_list, age_inf_data_list,
                          final_age_vctr, n_indiv, fixed_par, fixed_pos,
                          seed_vec)
{

  pars <- append(estimated_pars, fixed_par, after = fixed_pos - 1)

  log_lik_titre_prob_inf(pars, titre_data_list, age_inf_data_list,
                         final_age_vctr, n_indiv, seed_vec)
}

optimise_over_fixed_value <- function(starting_point, id, titre_data_list,
                                      age_inf_data_list, final_age_vctr,
                                      fixed_pos, n_indiv)
{
  set.seed(1742)
  seed_vec <- sample.int(1e7, 10)

  fn <- str_glue("./saved_objects/inference/two_datasets/profile_{fixed_pos}/iter_{id}.rds")

  if(!file.exists(fn))
  {
    fixed_val     <- starting_point[[fixed_pos]]
    unc_fixed_val <- link_funs[[fixed_pos]](fixed_val)

    starting_point[[fixed_pos]] <- NULL

    par_vctr <- unlist(starting_point, use.names = TRUE)

    names_vctr <- names(par_vctr)

    #Unconstrained (unc)
    unc_par_vctr <- map_dbl(names_vctr, \(nm) link_funs[[nm]](par_vctr[[nm]]))

    names(unc_par_vctr) <- names_vctr

    res <- nloptr(
      x0                = unc_par_vctr,
      eval_f            = profile_fun_2,
      titre_data_list   = titre_data_list,
      age_inf_data_list = age_inf_data_list,
      final_age_vctr    = final_age_vctr,
      n_indiv           = n_indiv,
      fixed_par         = unc_fixed_val,
      fixed_pos         = fixed_pos,
      seed_vec          = seed_vec,
      opts = list(algorithm = "NLOPT_LN_SBPLX",
                  maxeval     = 3500,
                  xtol_rel    = 1e-8,
                  ftol_rel    = 1e-10,
                  print_level = 0))

    saveRDS(res, fn)

  } else res <- readRDS(fn)

  res
}

construct_profile <- function(starting_points, fixed_pos, titre_data_list,
                              age_inf_data_list, final_age_vctr, n_indiv)
{
  future_imap(
    starting_points,
    optimise_over_fixed_value,
    titre_data_list   = titre_data_list,
    age_inf_data_list = age_inf_data_list,
    final_age_vctr    = final_age_vctr,
    fixed_pos         = fixed_pos,
    n_indiv           = n_indiv)
}

format_profile_df <- function(obj_list, vals)
{
  param_profile_df <- map_df(obj_list, \(prof_obj) {

    sol        <- prof_obj$solution
    names(sol) <- setdiff(colnames, param)
    ll_val     <- -prof_obj$objective

    for(nm in names(sol))
    {
      sol[[nm]] <- inverse_link_funs[[nm]](sol[[nm]])
    }

    df_sol <- sol |> t() |> as.data.frame() |>
      mutate(ll = ll_val)
  }) |>
    mutate(value = vals) |>
    group_by(value) |>
    filter(ll == max(ll)) |>
    ungroup() |>
    rename("{param}" := value)

  param_profile_df <- param_profile_df[, c(colnames, "ll")]

  param_profile_df
}
