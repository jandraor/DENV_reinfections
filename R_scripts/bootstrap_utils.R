boot_mean <- function(data_vct, indices) {
  resampled_data_vct <- data_vct[indices] # resample data by indices
  return(mean(resampled_data_vct))        # compute mean of resample
}
