profile_fun <- function(estimated_pars, age_inf_data_list,
                          final_age_vctr, n_indiv, fixed_par, fixed_pos,
                          seed_vec)
{

  pars <- append(estimated_pars, fixed_par, after = fixed_pos - 1)

  log_lik_prob_inf(pars, age_inf_data_list, final_age_vctr, n_indiv, seed_vec)
}

optimise_over_fixed_value <- function(starting_point, id, age_inf_data_list,
                                      final_age_vctr, fixed_pos, n_indiv)
{
  set.seed(1544)
  seed_vec <- sample.int(1e7, 10)

  fn <- str_glue("./saved_objects/inference/one_dataset/profile_{fixed_pos}/iter_{id}.rds")

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
      eval_f            = profile_fun,
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
profile_fun_2 <- function(estimated_pars, titre_data_list, age_inf_data_list,
                          final_age_vctr, n_indiv, fixed_par, fixed_pos,
                          seed_vec)
{

  pars <- append(estimated_pars, fixed_par, after = fixed_pos - 1)

  log_lik_titre_prob_inf(pars, titre_data_list, age_inf_data_list,
                         final_age_vctr, n_indiv, seed_vec)
}

optimise_over_fixed_value_2 <- function(starting_point, id, titre_data_list,
                                        age_inf_data_list, final_age_vctr,
                                        fixed_pos, n_indiv, ds)
{
  set.seed(1742)
  seed_vec <- sample.int(1e7, 10)

  fn <- str_glue("./saved_objects/inference/{ds}/profile_{fixed_pos}/iter_{id}.rds")

  if(!file.exists(fn))
  {
    fixed_name    <- names(starting_point)[[fixed_pos]]
    fixed_val     <- starting_point[[fixed_pos]]
    unc_fixed_val <- link_funs[[fixed_name]](fixed_val)

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

    dir.create(dirname(fn), recursive = TRUE, showWarnings = FALSE)
    saveRDS(res, fn)

    if (file.exists(fn)) {
      message(str_glue("Finished optimisation {id} (fixed_pos = {fixed_pos}) -> saved to {fn}"))
    } else {
      warning(str_glue("Optimisation {id} finished but file was NOT written: {fn}"))
    }

  } else res <- readRDS(fn)

  res
}

construct_profile <- function(starting_points, fixed_pos, titre_data_list,
                              age_inf_data_list, final_age_vctr, n_indiv,
                              ds, computation = "parallel")
{
  if(ds == "one_dataset")
  {
    profile_res <- future_imap(
      starting_points,
      optimise_over_fixed_value,
      age_inf_data_list = age_inf_data_list,
      final_age_vctr    = final_age_vctr,
      fixed_pos         = fixed_pos,
      n_indiv           = n_indiv)
  }

  if(ds == "two_datasets")
  {
    profile_res <- future_imap(
      starting_points,
      optimise_over_fixed_value_2,
      titre_data_list   = titre_data_list,
      age_inf_data_list = age_inf_data_list,
      final_age_vctr    = final_age_vctr,
      fixed_pos         = fixed_pos,
      n_indiv           = n_indiv,
      ds                = ds)
  }

  if(ds == "alternative")
  {
    arg_list <- list(.x = starting_points,
                     .f = optimise_over_fixed_value_2,
                     titre_data_list = titre_data_list,
                     age_inf_data_list = age_inf_data_list,
                     final_age_vctr    = final_age_vctr,
                     fixed_pos         = fixed_pos,
                     n_indiv           = n_indiv,
                     ds                = ds)

    if(computation == "parallel")
    {
      profile_res <- do.call(future_imap, arg_list)
    }

    if(computation == "sequential")
    {
      profile_res <- do.call(imap, arg_list)
    }
  }

  profile_res
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

alternative_profile <- function(param,
                                titre_data_list,
                                age_inf_data_list,
                                final_age_vctr,
                                n_indiv,
                                box,
                                computation = "parallel")
{
  arg_list <- list(titre_data_list   = titre_data_list,
                   age_inf_data_list = age_inf_data_list,
                   final_age_vctr    = final_age_vctr,
                   n_indiv           = n_indiv,
                   ds                = "alternative",
                   computation       = computation)

  param_names <- c("lambda_1", "lambda_2", "log_A0", "phi", "sd")

  if (!param %in% param_names) {
    stop(paste0(param," is an invalid parameter name."))
  }

  if(param == "lambda_1")
  {
    grid_vals <- seq(0.038, 0.04, by = 0.00001)

    set.seed(1656)

    profile_design(
      lambda_1  = grid_vals,
      lower     = box[1, setdiff(param_names, param)],
      upper     = box[2, setdiff(param_names, param)],
      nprof     = 20,
      type      = "sobol") -> guesses_df
  }

  if(param == "lambda_2")
  {
    grid_vals <- seq(0.05, 0.07, by = 0.00025)

    set.seed(0954)

    guesses_df <- profile_design(
      lambda_2  = grid_vals,
      lower = box[1, setdiff(param_names, param)],
      upper = box[2, setdiff(param_names, param)],
      nprof = 20, type = "sobol")
  }

  if(param == "log_A0")
  {
    grid_vals <- seq(1.6, 2.0, by = 0.0025)

    set.seed(0836)

    guesses_df <- profile_design(
      log_A0 = grid_vals,
      lower  = box[1, setdiff(param_names, param)],
      upper  = box[2, setdiff(param_names, param)],
      nprof  = 20, type = "sobol")
  }

  if(param == "phi")
  {
    grid_vals <- seq(6.6, 7.2, by = 0.01)

    set.seed(2109)

    guesses_df <- profile_design(
      phi    = grid_vals,
      lower  = box[1, setdiff(param_names, param)],
      upper  = box[2, setdiff(param_names, param)],
      nprof  = 20, type = "sobol")
  }

  if(param == "sd")
  {
    grid_vals <- seq(2.2, 2.9, by = 0.01)

    set.seed(2118)

    guesses_df <- profile_design(
      sd     = grid_vals,
      lower  = box[1, setdiff(param_names, param)],
      upper  = box[2, setdiff(param_names, param)],
      nprof  = 20, type = "sobol")
  }

  message(str_glue("...Constructing profile for {param}..."))

  guesses_df <- guesses_df[, param_names]

  # starting points(sp)
  sp_list <- transpose(guesses_df)

  arg_list$fixed_pos       <- which(param_names %in% param)
  arg_list$starting_points <- sp_list

  prof_objs <- do.call(construct_profile, arg_list)

  list(prof_objs   = prof_objs,
       grid_values = guesses_df[[param]])
}
