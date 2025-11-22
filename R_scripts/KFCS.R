dmy2 <- function(date_str, threshold = 2025) {
  date <- dmy(date_str)               # parse normally
  year(date) <- ifelse(year(date) > threshold, year(date) - 100, year(date))
  date
}

KFCS_read_data <- function()
{
  read_csv("./data/KFCS/KFCS_longData.csv") |> 
    mutate(
      start_interval = dmy(START_dateSpecimenCollected),
      end_interval   = dmy(END_dateSpecimenCollected),
      dateBirth2 = ifelse(str_detect(dateBirth, "dd"),
                          str_glue("01-01-{str_sub(dateBirth, -4, -1)}"), 
                          dateBirth)) |> 
    mutate(dateBirth2 = dmy2(dateBirth2))
}

KFCS_get_infection_df <- function()
{
  KFCS_df <- KFCS_read_data()
  
  infection_df <- KFCS_df |> 
    select(subjectNo, dateBirth2, start_interval, end_interval, xgBoostPred) |> 
    mutate(age_start    = as.numeric(start_interval - dateBirth2) / 365.25,
           age_end      = as.numeric(end_interval - dateBirth2) / 365.25,
           age_round    = round((age_start + age_end) / 2),
           is_titre_inf = ifelse(xgBoostPred < 0.8, FALSE, TRUE))
  
  KFCS_PCR <-  KFCS_get_symp_infections()
  
  df_list <- split(infection_df, infection_df$subjectNo)
  
  
  infection_df <- imap_dfr(df_list, \(df, subject_id) {
    
    subject_PCR <- KFCS_PCR |> filter(subjectNo == subject_id)
    
    df$PCR <- FALSE
    
    if(nrow(subject_PCR) > 0)
    {
      PCR_dates <- subject_PCR$dateEvaluationA1
      
      df$PCR <- sapply(
        seq_len(nrow(df)), 
        function(i) any(PCR_dates >= df$start_interval[i] & 
                          PCR_dates <= df$end_interval[i]))
      
    }
    
    df
  }) |> 
    mutate(is_inf = is_titre_inf | PCR)
}

KFCS_get_prob_inf_at_age <- function()
{
  infection_df <- KFCS_get_infection_df()

  
  KFCS_age_inf <- infection_df |> group_by(age_round) |> 
    summarise(n_infections  = sum(is_inf),
              n_individuals = n()) |> 
    mutate(pct = n_infections / n_individuals) |> 
    filter(!is.na(age_round))
}

KFCS_plot_prob_inf_age <- function()
{
  KFCS_age_inf <- KFCS_get_prob_inf_at_age()
  
  ggplot(KFCS_age_inf |> filter(age_round >=2), aes(age_round, pct)) +
    geom_line() + 
    theme_classic()
}

KFCS_get_symp_infections <- function()
{
  KFCS_PCR <- read_csv("./data/KFCS/Analysis_Illness_20240722.csv") |> 
    select(subjectNo, dateEvaluationA1, pcrResult) |> 
    mutate(dateEvaluationA1 = ymd(dateEvaluationA1)) |> 
    filter(pcrResult == "Dengue")
}










