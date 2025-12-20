profile_fun_2 <- function(estimated_pars, titre_data_list, age_inf_data_list,
                          final_age_vctr, n_indiv, fixed_par, fixed_pos)
{

  pars <- append(estimated_pars, fixed_par, after = fixed_pos - 1)

  log_lik_titre_prob_inf(pars, titre_data_list, age_inf_data_list,
                         final_age_vctr, n_indiv)
}

optimise_over_fixed_value <- function(fixed_val, id, MLE_vals,
                                      titre_data_list, age_inf_data_list,
                                      final_age_vctr, fixed_pos,
                                      n_indiv)
{
  fn <- str_glue("./saved_objects/inference/two_datasets/profile_{fixed_pos}/iter_{id}.rds")

  if(!file.exists(fn))
  {
    unc_fixed_val <- function_list[[fixed_pos]](fixed_val)


    MLE_vals[[fixed_pos]] <- NULL

    par_vctr <- unlist(MLE_vals, use.names = TRUE)

    names_vctr <- names(par_vctr)

    #Unconstrained (unc)
    unc_par_vctr <- map_dbl(names_vctr, \(nm) function_list[[nm]](par_vctr[[nm]]))

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
      opts = list(algorithm = "NLOPT_LN_SBPLX",
                  maxeval     = 20000,
                  xtol_rel    = 1e-8,
                  ftol_rel    = 1e-10,
                  print_level = 0))

    saveRDS(res, fn)

  } else res <- readRDS(fn)

  res
}

construct_profile <- function(fixed_vals, fixed_pos, MLE_vals, titre_data_list,
                              age_inf_data_list, final_age_vctr, n_indiv)
{
  future_imap(
    fixed_vals,
    optimise_over_fixed_value,
    MLE_vals = MLE_vals,
    titre_data_list = titre_data_list,
    age_inf_data_list = age_inf_data_list,
    final_age_vctr = final_age_vctr,
    fixed_pos  = fixed_pos,
    n_indiv    = n_indiv)
}
