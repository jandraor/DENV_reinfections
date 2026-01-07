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
    mutate(pct = n / n_individuals)
}

find_MLE <- function(start_list, age_inf_data_list, final_age_vctr, n_indiv)
{
  set.seed(1422)
  seed_vec <- sample.int(1e7, 10)

  future_map(start_list, \(start_obj) {

    iter_id <- start_obj$iter_id

    fn <- str_glue("./saved_objects/inference/one_dataset/MLE/opt_{iter_id}.rds")

    if(!file.exists(fn))
    {
      start_obj[["iter_id"]] <- NULL

      par_vctr <- unlist(start_obj, use.names = TRUE)

      res <- nloptr(
        x0                = par_vctr,
        eval_f            = log_lik_prob_inf,
        age_inf_data_list = age_inf_data_list,
        final_age_vctr    = final_age_vctr,
        n_indiv           = n_indiv,
        seed_vec          = seed_vec,
        opts = list(algorithm = "NLOPT_LN_SBPLX",
                    maxeval   = 3000,
                    xtol_rel  = 1e-8,
                    ftol_rel  = 1e-10,
                    print_level = 0))

      saveRDS(res, fn)
    } else res <- readRDS(fn)

    res
  })

}


find_MLE_2 <- function(start_list, titre_data_list, age_inf_data_list,
                     final_age_vctr, n_indiv)
{
  set.seed(1742)
  seed_vec <- sample.int(1e7, 10)

  future_map(start_list, \(start_obj) {

    iter_id <- start_obj$iter_id

    fn <- str_glue("./saved_objects/inference/two_datasets/MLE/opt_{iter_id}.rds")

    if(!file.exists(fn))
    {
      message(str_glue("Starting iter: {iter_id}"))

      start_obj[["iter_id"]] <- NULL

      par_vctr <- unlist(start_obj, use.names = TRUE)

      res <- nloptr(
        x0                = par_vctr,
        eval_f            = log_lik_titre_prob_inf,
        titre_data_list   = titre_data_list,
        age_inf_data_list = age_inf_data_list,
        final_age_vctr    = final_age_vctr,
        n_indiv           = n_indiv,
        seed_vec          = seed_vec,
        opts = list(algorithm = "NLOPT_LN_SBPLX",
                    maxeval   = 3000,
                    xtol_rel  = 1e-8,
                    ftol_rel  = 1e-10,
                    print_level = 0))

      message(str_glue("Finished iter: {iter_id}"))

      saveRDS(res, fn)
    } else res <- readRDS(fn)

    res
  })
}

link_funs <- list(
  lambda_1 = logit,
  lambda_2 = logit,
  rho      = logit,
  log_A0   = log,
  phi      = log,
  sd       = log)

inverse_link_funs <- list(
  lambda_1 = inv.logit,
  lambda_2 = inv.logit,
  rho      = inv.logit,
  log_A0   = exp,
  phi      = exp,
  sd       = exp)


get_starting_points <- function(ds)
{
  if(ds == "one_dataset")
  {
    set.seed(07012026)

    start_df <- sobol_design(
      lower = c("lambda_1" = 0.01,
                "lambda_2" = 0.01,
                "rho"      = 0.001),
      upper = c("lambda_1" = 0.25,
                "lambda_2" = 0.25,
                "rho"      = 0.25),
      nseq =  200)


  }

  if(ds == "two_datasets")
  {
    set.seed(111020225)

    start_df <- sobol_design(
      lower = c("lambda_1" = 0.01,
                "lambda_2" = 0.01,
                "rho"      = 0.001,
                "log_A0"   = 0.1,
                "phi"      = 1,
                "sd"       = 0.01),
      upper = c("lambda_1" = 0.25,
                "lambda_2" = 0.25,
                "rho"      = 0.25,
                "log_A0"   = 4,
                "phi"      = 10,
                "sd"       = 10),
      nseq =  200)
  }

  par_names <- colnames(start_df)


  for (nm in  par_names)
  {
    start_df[[nm]] <- link_funs[[nm]](start_df[[nm]])
  }

  start_df <- start_df |> mutate(iter_id = row_number(),
                                 .before = everything())

  start_list <- transpose(start_df)

  start_list
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

log_lik_prob_inf <- function(pars, age_inf_data_list, final_age_vctr, n_indiv,
                             seed_vec)
{
  # Infection parameters-------------------
  lambda_1 <- inv.logit(pars[[1]])
  lambda_2 <- inv.logit(pars[[2]])
  lambdas  <- c(lambda_1, lambda_2)
  rho      <- inv.logit(pars[[3]])

  # cat("\n---------------")
  # cat("\n lambda 1: ", lambdas[[1]])
  # cat("\n lambda 2: ", lambdas[[2]])
  # cat("\n rho: ", rho)

  n_replicates <- length(seed_vec)

  ll_vals <- map_dbl(seq_len(n_replicates), \(k) {

    set.seed(seed_vec[k])

    inf_list <- sapply(1:2, \(cohort_idx) {

      inf_df <- simulate_DENV_infections_since_birth(
        lambda_serotype   = lambdas[[cohort_idx]],
        loss_rate         = rho,
        final_age         = final_age_vctr[[cohort_idx]],
        n_individuals     = n_indiv)

      age_inf_data <- age_inf_data_list[[cohort_idx]]

      age_tab      <- tabulate(inf_df$age + 1, nbins = max(final_age_vctr[[cohort_idx]]) + 1)
      prob_inf_age <- age_tab / n_indiv

      dbinom(x    = age_inf_data$n_infections,
             size = age_inf_data$n_individuals,
             prob = prob_inf_age[-c(1, 2)],
             log  = TRUE) |> sum()
    }) |> sum()
  })

  mean_ll <- -mean(ll_vals) # negative log-likelihood
  mcse    <- sd(ll_vals) / sqrt(length(ll_vals))

  # Attach MCSE as an attribute
  attr(mean_ll, "MCSE") <- mcse

  mean_ll
}

log_lik_titre_prob_inf <- function(pars, titre_data_list, age_inf_data_list,
                                   final_age_vctr, n_indiv, seed_vec)
{
  # Infection parameters-------------------
  lambda_1 <- inv.logit(pars[[1]])
  lambda_2 <- inv.logit(pars[[2]])
  lambdas  <- c(lambda_1, lambda_2)
  rho      <- inv.logit(pars[[3]])

  # Decay rate dynamics--------------------
  decay_rate <- get_decay_rate()


  # Peak dynamics--------------------------
  log_A0  <- exp(pars[[4]])
  phi     <- exp(pars[[5]])
  beta    <- 1

  sd <- exp(pars[[6]])

  # cat("\n---------------")
  # cat("\n lambda 1: ", lambdas[[1]])
  # cat("\n lambda 2: ", lambdas[[2]])
  # cat("\n rho: ", rho)
  # cat("\n log A0: ", log_A0 )
  # cat("\n phi: ", phi)
  # cat("\n sd: ", sd)

  n_replicates <- length(seed_vec)

  ll_vals <- map_dbl(seq_len(n_replicates), \(k) {

    set.seed(seed_vec[k])

    inf_list <- lapply(1:2, \(cohort_idx) {

      simulate_DENV_infections_since_birth(
        lambda_serotype   = lambdas[[cohort_idx]],
        loss_rate         = rho,
        final_age         = final_age_vctr[[cohort_idx]],
        n_individuals     = n_indiv)
    })

    ll_age_inf <- sapply(1:2, \(cohort_idx) {

      inf_df       <- inf_list[[cohort_idx]]
      age_inf_data <- age_inf_data_list[[cohort_idx]]

      age_tab <- tabulate(inf_df$age + 1, nbins = max(final_age_vctr[[cohort_idx]]) + 1)
      prob_inf_age <- age_tab / n_indiv

      dbinom(x    = age_inf_data$n_infections,
             size = age_inf_data$n_individuals,
             prob = prob_inf_age[-c(1, 2)],
             log  = TRUE) |> sum()
    }) |> sum()

    ll_titre <- sapply(1:2, \(cohort_idx) {

      inf_df <- inf_list[[cohort_idx]]

      if(nrow(inf_df) == 0) return(-Inf)

      ids <- seq_len(n_indiv)

      inf_df <- inf_df[order(inf_df$subject_id, inf_df$age), ]

      inf_times_list <- split(inf_df$age,
                              factor(inf_df$subject_id, levels = ids))

      avg_titre_vctr <- estimate_avg_titre_by_age(
        log_first_peak = log_A0,
        decay_rate     = decay_rate,
        phi            = phi,
        beta           = beta,
        inf_times_list = inf_times_list,
        final_age      = final_age_vctr[[cohort_idx]])

      titre_df <- titre_data_list[[cohort_idx]]

      dnorm(x    = titre_df$mean,
            mean = avg_titre_vctr[-1],
            sd   = sd / sqrt(titre_df$n), log = TRUE) |>
      sum()
    }) |> sum()

    ll <- ll_age_inf + ll_titre

    ll
  })

  #cat("\n Log lik: ", -mean(ll_vals))

  # Compute mean log-likelihood and MCSE
  mean_ll <- -mean(ll_vals)            # negative log-likelihood
  mcse    <- sd(ll_vals) / sqrt(length(ll_vals))

  # Attach MCSE as an attribute
  attr(mean_ll, "MCSE") <- mcse

  mean_ll
}

get_loglik_2 <- function()
{
  fldr <- "./saved_objects/inference/two_datasets/MLE"

  files <- list.files(path = fldr, pattern = "^opt")

  sol_df <- map_df(files, \(fn) {

    fp <- file.path(fldr, fn)

    res <- readRDS(fp)

    sol <- res$solution

    status <- res$status

    if(status == 5) return (NULL)

    names(sol) <- c("lambda_1",
                    "lambda_2",
                    "rho" ,
                    "log_A0" ,
                    "phi" ,
                    "sd")

    for (nm in names(inverse_link_funs))
    {
      sol[[nm]] <- inverse_link_funs[[nm]](sol[[nm]])
    }

    as.data.frame(as.list(sol)) |>
      mutate(ll = -res$objective)
  })

  sol_df
}

link_pars <- function(param_obj)
{
  for(nm in names(param_obj))
  {
    param_obj[[nm]] <- link_funs[[nm]](param_obj[[nm]])
  }

  param_obj
}

get_decay_rate <- function() 0.3 * exp(-0.5*(0:3))

source("./R_scripts/inference_data.R")
source("./R_scripts/inference_profile.R")
