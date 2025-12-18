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



find_MLE_2 <- function(start_list, titre_data_list, age_inf_data_list,
                     final_age_vctr)
{
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
        n_indiv           = 1e4,
        opts = list(algorithm = "NLOPT_LN_SBPLX",
                    maxeval   = 20000,
                    xtol_rel  = 1e-8,
                    ftol_rel  = 1e-10))

      message(str_glue("Finished iter: {iter_id}"))

      saveRDS(res, fn)
    } else res <- readRDS(fn)

    res
  })
}

function_list <- list(
  lambda_1 = logit,
  lambda_2 = logit,
  rho      = logit,
  log_A0   = log,
  phi      = log,
  sd_1     = log,
  sd_2     = log)

inverse_link_funs <- list(
  lambda_1 = inv.logit,
  lambda_2 = inv.logit,
  rho      = inv.logit,
  log_A0   = exp,
  phi      = exp,
  sd_1     = exp,
  sd_2     = exp)


get_starting_points <- function()
{
  set.seed(111020225)

  start_df <- sobol_design(
    lower = c("lambda_1" = 0.01,
              "lambda_2" = 0.01,
              "rho"      = 0.001,
              "log_A0"   = 0.1,
              "phi"      = 1,
              "sd_1"     = 0.01,
              "sd_2"     = 0.01),
    upper = c("lambda_1" = 0.25,
              "lambda_2" = 0.25,
              "rho"      = 0.25,
              "log_A0"   = 4,
              "phi"      = 10,
              "sd_1"     = 10,
              "sd_2"     = 10),
    nseq =  200)

  for (nm in names(function_list))
  {
    start_df[[nm]] <- function_list[[nm]](start_df[[nm]])
  }

  start_df <- start_df |> mutate(iter_id = row_number(),
                                 .before = everything())

  start_list <- transpose(start_df)
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

log_lik_titre_prob_inf <- function(pars, titre_data_list, age_inf_data_list,
                                   final_age_vctr, n_indiv)
{
  # Infection parameters-------------------
  lambda_1 <- inv.logit(pars[[1]])
  lambda_2 <- inv.logit(pars[[2]])
  lambdas  <- c(lambda_1, lambda_2)
  rho      <- inv.logit(pars[[3]])

  # Decay rate dynamics--------------------
  decay_rate <- 0.2 * exp(-0.5*(0:3))


  # Peak dynamics--------------------------
  log_A0  <- exp(pars[[4]])
  phi     <- exp(pars[[5]])
  beta    <- 1

  sd_vals <- exp(pars[6:7])

  # cat("\n---------------")
  # cat("\n log A0: ", log_A0 )
  # cat("\n Decay rate: ", decay_rate)
  # cat("\n phi: ", phi)
  # cat("\n sd_val: ", sd_vals[[1]])
  # cat("\n sd_val2: ", sd_vals[[2]])
  # cat("\n lambda 1: ", lambdas[[1]])
  # cat("\n lambda 2: ", lambdas[[2]])
  # cat("\n rho: ", rho)
  # cat("\n r_1: ", r_1)
  #cat("\n alpha: ", alpha)
  #cat("\n beta: ", beta)

  set.seed(1150)

  inf_list <- lapply(1:2, \(cohort_idx) {

    simulate_DENV_infections_since_birth(
      lambda_serotype   = lambdas[[cohort_idx]],
      loss_rate         = rho,
      final_age         = final_age_vctr[[cohort_idx]],
      n_individuals     = n_indiv) |>
      rename(subject_id = infected_ind)
  })

  ll_age_inf <- sapply(1:2, \(cohort_idx) {

    inf_df       <- inf_list[[cohort_idx]]
    age_inf_data <- age_inf_data_list[[cohort_idx]]

    age_tab <- tabulate(inf_df$age + 1, nbins = max(final_age_vctr[[cohort_idx]]) + 1)
    prob_inf_age <- age_tab / n_indiv

    dbinom(x    = age_inf_data$n_infections,
           size = age_inf_data$n_individuals,
           prob = prob_inf_age[-1],
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
          sd   = sd_vals[[cohort_idx]] / sqrt(titre_df$n), log = TRUE) |>
      sum()
  }) |> sum()

  ll <- ll_age_inf + ll_titre


  # cat("\n ll_age_inf:", ll_age_inf)
  # cat("\n log_prior: ", log_prior)
  # cat("\n log_lik: ", ll_age_inf + ll_titre)
  # cat("\n---------------")

  -ll
}

get_MLE_2 <- function()
{
  fldr <- "./saved_objects/inference/two_datasets/MLE"

  files <- list.files(path = fldr, pattern = "^opt")

  sol_df <- map_df(files, \(fn) {

    fp <- file.path(fldr, fn)

    res <- readRDS(fp)

    sol <- res$solution

    names(sol) <- c("lambda_1",
                    "lambda_2",
                    "rho" ,
                    "log_A0" ,
                    "phi" ,
                    "sd_1" ,
                    "sd_2")

    for (nm in names(inverse_link_funs))
    {
      sol[[nm]] <- inverse_link_funs[[nm]](sol[[nm]])
    }

    as.data.frame(as.list(sol)) |>
      mutate(ll = -res$objective)
  })

  sol_df |> filter(ll == max(ll)) |> as.list()
}

source("./R_scripts/inference_profile.R")
