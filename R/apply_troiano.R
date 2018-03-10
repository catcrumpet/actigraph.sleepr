#' Apply the Troiano algorithm
#'
#' The Troiano algorithm detects periods of non-wear in activity data from an ActiGraph device. Such intervals are likely to represent invalid data and therefore should be excluded from downstream analysis. The algorithm formalizes a technique used to analyze the 2003-2004 NHANES data; the original SAS source code can be found at \url{http://riskfactor.cancer.gov/tools/nhanes_pam/}.
#' @param agdb A \code{tibble} (\code{tbl}) of activity data (at least) an \code{epochlength} attribute.
#' @param activity_threshold Highest activity per minute level to be considered "zero"; an epoch with activity exceeding the threshold is considered a "spike". The default threshold is 0.
#' @param min_period_len Minimum number of consecutive "zero" minutes to start a non-wear period. The default is 60.
#' @param max_nonzero_count Epochs with activity greater than \code{max_nonzero_count} per minute are labeled as "zero". The default is Inf.
#' @param spike_tolerance Also known as artifactual movement interval. At most \code{spike_tolerance} "nonzero" minutes can occur in sequence during a non-wear period without interrupting it. The default is 2.
#' @param spike_stoplevel An activity spike that exceeds \code{spike_stoplevel} per minute counts ends a non-wear period, even if the spike tolerance has not been reached. The default is 100.
#' @param use_magnitude Logical. If true, the magnitude of the axis vector (e.g., axis1, axis2, axis3) is used to measure activity; otherwise the axis1 value is used. The default is FALSE.
#' @param endat_nnz_seq Logical. If true, a non-wear period ends with a run of nonzero minutes that is longer than \code{spike_tolerance}. This corresponds to the option \emph{"Require consecutive epochs outside of the activity threshold"} in ActiLife's Wear Time Validation menu. The default is TRUE.
#' @details
#' The Troiano algorithm specifies that a non-wear period starts with \code{min_period_len} consecutive minutes of "zero" activity and ends with more than \code{spike_tolerance} minutes of "nonzero" activity.
#'
#' This implementation of the algorithm has been adapted to epoch lengths that are divisors of 60. Function parameters are denominated in minutes and then converted to appropriate values in epochs. For example, when using data with 30-second epochs, setting \code{activity_threshold} to 100 is interpreted as 100 counts per minute, which would be converted to 50 per 30-second epoch, and specifying a \code{spike_tolerance} of 10 is interpreted as 10 minutes, which would be converted to 20 30-second epochs. This conversion does not necessarily produce equivalent results between different epoch lengths. Two sequential 30-second epochs may have activity counts of 51 and 0, the first of which would be above a \code{spike_stoplevel} of 100 counts per minute (50 counts per 30 seconds). However, an equivalent 60-second epoch would have an activity count of 51, which would fall below the same threshold.
#' @return A summary \code{tibble} of the detected non-wear periods. If the activity data is grouped, then non-wear periods are detected separately for each group.
#' @references RP Troiano, D Berrigan, KW Dodd, LC Mâsse, T Tilert and M McDowell. Physical activity in the united states measured by accelerometer. \emph{Medicine & Science in Sports & Exercise}, 40(1):181–188, 2008.
#' @references ActiLife 6 User's Manual by the ActiGraph Software Department. 04/03/2012.
#' @seealso \code{\link{apply_choi}}, \code{\link{collapse_epochs}}
#' @examples
#' library("dplyr")
#' data("gtxplus1day")
#'
#' gtxplus1day %>%
#'   collapse_epochs(60) %>%
#'   apply_troiano()
#' @export

apply_troiano <- function(agdb,
                          activity_threshold = 0,
                          min_period_len = 60,
                          max_nonzero_count = Inf,
                          spike_tolerance = 2,
                          spike_stoplevel = 100,
                          use_magnitude = FALSE,
                          endat_nnz_seq = TRUE) {

  check_args_nonwear_periods(agdb, use_magnitude)

  epoch_len <- attr(agdb, "epochlength")
  epochs_per_min <- 60L / epoch_len

  activity_threshold <- activity_threshold / epochs_per_min
  min_period_len <- min_period_len * epochs_per_min
  max_nonzero_count <- max_nonzero_count / epochs_per_min
  spike_tolerance <- spike_tolerance * epochs_per_min
  spike_stoplevel <- spike_stoplevel / epochs_per_min

  if (endat_nnz_seq) {
    nonwear <- agdb %>%
      do(apply_troiano_seq_(.data,
                            epoch_len,
                            activity_threshold,
                            min_period_len,
                            max_nonzero_count,
                            spike_tolerance,
                            spike_stoplevel,
                            use_magnitude))
  } else {
    nonwear <- agdb %>%
      do(apply_troiano_nonseq_(.data,
                               epoch_len,
                               activity_threshold,
                               min_period_len,
                               max_nonzero_count,
                               spike_tolerance,
                               spike_stoplevel,
                               use_magnitude))
  }

  nonwear <-
    structure(nonwear,
              class = c("tbl_period", "tbl_df", "tbl", "data.frame"),
              nonwear_algorithm = "Troiano",
              min_period_len = min_period_len,
              max_nonzero_count = max_nonzero_count,
              spike_tolerance = spike_tolerance,
              spike_stoplevel = spike_stoplevel,
              activity_threshold = activity_threshold,
              endat_nnz_seq = endat_nnz_seq,
              use_magnitude = use_magnitude,
              epochlength = attr(agdb, "epochlength"))

  if (is.grouped_df(agdb))
    nonwear <- nonwear %>% group_by(!!! groups(agdb))

  nonwear
}

apply_troiano_seq_ <- function(data,
                               epoch_len,
                               activity_threshold,
                               min_period_len,
                               max_nonzero_count,
                               spike_tolerance,
                               spike_stoplevel,
                               use_magnitude) {

  data %>%
    add_magnitude() %>%
    mutate(count = if (use_magnitude) magnitude else axis1,
           wear = if_else(count <= activity_threshold | count > max_nonzero_count, 0L, 1L),
           wear = if_else(count > spike_stoplevel, 2L, wear)) %>%
    group_by(rleid = rleid(wear)) %>%
    summarise(wear = first(wear),
              timestamp = first(timestamp),
              length = n()) %>%
    mutate(wear = if_else(wear == 1L &
                            lead(wear, default = 1L) == 0L &
                            length <= spike_tolerance,
                          NA_integer_, wear),
           # Since `na.locf` can't impute leading NAs, fill in those with 1s
           wear = if_else(row_number() == 1 & is.na(wear), 1L, wear),
           # Fill in NAs with the most recent zero/nonzero wear state
           wear = na.locf(wear)) %>%
    group_by(rleid = rleid(wear)) %>%
    summarise(wear = first(wear),
              timestamp = first(timestamp),
              length = sum(length)) %>%
    filter(wear == 0L, length >= min_period_len) %>%
    rename(period_start = timestamp) %>%
    mutate(period_end = period_start + seconds(length * epoch_len)) %>%
    select(period_start, period_end, length)
}

apply_troiano_nonseq_ <- function(data,
                                  epoch_len,
                                  activity_threshold,
                                  min_period_len,
                                  max_nonzero_count,
                                  spike_tolerance,
                                  spike_stoplevel,
                                  use_magnitude) {

  data %>%
    add_magnitude() %>%
    mutate(count = if (use_magnitude) magnitude else as.numeric(axis1),
           count = if_else(count > max_nonzero_count, 0, count),
           length = wle(count, activity_threshold, spike_tolerance, spike_stoplevel)) %>%
    filter(length >= min_period_len) %>%
    rename(period_start = timestamp) %>%
    select(period_start, length) %>%
    mutate(period_end = period_start + seconds(length * epoch_len)) %>%
    mutate(a = time_length(period_start - first(period_start), "second"),
           b = time_length(period_end - first(period_start), "second")) %>%
    # Remove periods which overlap with previous periods
    filter(overlap(a, b)) %>%
    select(period_start, period_end, length)
}
