#' Apply the Choi algorithm
#'
#' The Choi algorithm detects periods of non-wear in activity data from an ActiGraph device. Such intervals are likely to represent invalid data and therefore should be excluded from downstream analysis.
#' @param agdb A \code{tibble} (\code{tbl}) of activity data (at least) an \code{epochlength} attribute. The epoch length must be an exact divisor of 60.
#' @inheritParams apply_troiano
#' @param min_period_len Minimum number of consecutive "zero" minutes to start a non-wear period. The default is 90.
#' @param min_window_len The minimum number of consecutive "zero" minutes immediately preceding and following a spike of artifactual movement. The default is 30.
#' @details
#' The Choi algorithm extends the Troiano algorithm by requiring that short spikes of artifactual movement during a non-wear period are preceded and followed by \code{min_window_len} consecutive "zero" epochs.
#'
#' This implementation of the algorithm has been adapted to epoch lengths that are divisors of 60.
#' @return A summary \code{tibble} of the detected non-wear periods. If the activity data is grouped, then non-wear periods are detected separately for each group.
#' @references L Choi, Z Liu, CE Matthews and MS Buchowski. Validation of accelerometer wear and nonwear time classification algorithm. \emph{Medicine & Science in Sports & Exercise}, 43(2):357–364, 2011.
#' @references ActiLife 6 User's Manual by the ActiGraph Software Department. 04/03/2012.
#' @seealso \code{\link{apply_troiano}}, \code{\link{collapse_epochs}}
#' @examples
#' library("dplyr")
#' data("gtxplus1day")
#'
#' gtxplus1day %>%
#'   collapse_epochs(60) %>%
#'   apply_choi()
#' @export

apply_choi <- function(agdb,
                       min_period_len = 90,
                       min_window_len = 30,
                       spike_tolerance = 2,
                       use_magnitude = FALSE) {

  check_args_nonwear_periods(agdb, use_magnitude)
  stopifnot(min_window_len >= spike_tolerance)

  epoch_len <- attr(agdb, "epochlength")
  epochs_per_min <- 60L / epoch_len

  min_period_len <- min_period_len * epochs_per_min
  min_window_len <- min_window_len * epochs_per_min
  spike_tolerance <- spike_tolerance * epochs_per_min

  nonwear <- agdb %>%
    do(apply_choi_(.data,
                   epoch_len,
                   min_period_len,
                   min_window_len,
                   spike_tolerance,
                   use_magnitude))
  nonwear <-
    structure(nonwear,
              class = c("tbl_period", "tbl_df", "tbl", "data.frame"),
              nonwear_algorithm = "Choi",
              min_period_len = min_period_len,
              min_window_len = min_window_len,
              spike_tolerance = spike_tolerance,
              use_magnitude = use_magnitude,
              epochlength = attr(agdb, "epochlength"))

  if (is.grouped_df(agdb))
    nonwear <- nonwear %>% group_by(!!! groups(agdb))

  nonwear
}

apply_choi_ <- function(data,
                        epoch_len,
                        min_period_len,
                        min_window_len,
                        spike_tolerance,
                        use_magnitude) {
  data %>%
    add_magnitude() %>%
    mutate(count = if (use_magnitude) magnitude else axis1,
           wear = if_else(count == 0, 0L, 1L)) %>%
    group_by(rleid = rleid(wear)) %>%
    summarise(wear = first(wear), timestamp = first(timestamp), length = n()) %>%
    # Let (spike, zero, zero, spike) -> (spike of length 4)
    # as long as (zero, zero) is shorter than spike_tolerance
    mutate(wear = if_else(wear == 0L & length < spike_tolerance, 1L, wear)) %>%
    group_by(rleid = rleid(wear)) %>%
    summarise(wear = first(wear), timestamp = first(timestamp), length = sum(length)) %>%
    # Ignore artifactual movement intervals
    mutate(wear =
             if_else(wear == 1L &
                       length <= spike_tolerance &
                       lead(length, default = 0L) >= min_window_len &
                       lag(length, default = 0L) >= min_window_len,
                     0L, wear)) %>%
    group_by(rleid = rleid(wear)) %>%
    summarise(wear = first(wear), timestamp = first(timestamp), length = sum(length)) %>%
    filter(wear == 0L,
           # TODO: Filtering if the row_number is 1 or n(),
           # regardless of the period length, means that
           # the initial and final non-wear periods can be shorter.
           length >= min_period_len # | row_number() %in% c(1, n())
           ) %>%
    rename(period_start = timestamp) %>%
    mutate(period_end = period_start + seconds(length * epoch_len)) %>%
    select(period_start, period_end, length)
}
