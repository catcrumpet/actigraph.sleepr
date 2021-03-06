check_epochlen_is_60 <- function(agdb, algorithm) {
  if (attr(agdb, "epochlength") != 60)
    stop("Epochs should have length 60s to apply ", algorithm, ". Epochs can ",
         "be aggregated with `collapse_epochs`.")
}

check_epochlen_divisor_60 <- function(agdb) {
  if (60 %% attr(agdb, "epochlength"))
    stop("Epochs should be exact divisors of 60.")
}

check_no_missing_timestamps <- function(agdb) {
  if (has_missing_epochs(agdb))
    stop("Missing timestamps. Epochs should be evenly spaced from ",
         "first(timestamp) to last(timestamp).")
}

check_no_missing_counts <- function(agdb, var) {
  if (anyNA(agdb[[var]]))
    stop("Missing ", var, " counts. These can be imputed with `impute_epochs`.")
}

check_no_missing_state <- function(agdb) {
  if (!exists("sleep", agdb))
    stop("Missing asleep/awake (0/1) indicator column. These states can be ",
         "inferred with `apply_sadeh` or `apply_cole_kripke`.")
  if (anyNA(agdb$sleep))
    stop("Missing asleep/awake values.")
}

check_has_variable <- function(agdb, var) {
  if (!exists(var, where = agdb))
    stop("tbl_agd does not have variable ", var)
}

check_args_sleep_scores <- function(agdb, algorithm) {
  check_epochlen_is_60(agdb, algorithm)
  check_no_missing_timestamps(agdb)
  check_no_missing_counts(agdb, "axis1")
}

check_args_sleep_periods <- function(agdb, algorithm) {
  check_epochlen_is_60(agdb, algorithm)
  check_no_missing_timestamps(agdb)
  check_no_missing_state(agdb)
}

check_args_nonwear_periods <- function(agdb, use_magnitude) {
  check_epochlen_divisor_60(agdb)
  check_no_missing_timestamps(agdb)
  check_no_missing_counts(agdb, "axis1")
  if (use_magnitude) {
    check_no_missing_counts(agdb, "axis2")
    check_no_missing_counts(agdb, "axis3")
  }
}

check_args_filter <- function(agdb, var) {
  check_has_variable(agdb, var)
  check_no_missing_timestamps(agdb)
  check_no_missing_counts(agdb, var)
}

check_args_collapse_method <- function(agdb, epoch_len_out) {
  if (epoch_len_out %% attr(agdb, "epochlength"))
    stop("Output epoch length is not an exact multiple of input epoch length.")
  check_no_missing_timestamps(agdb)
  check_no_missing_counts(agdb, "axis1")
}

check_args_pa_scores <- function(agdb, age, algorithm) {
  check_epochlen_divisor_60(agdb)
  check_no_missing_timestamps(agdb)
  check_no_missing_counts(agdb, "axis1")
  if (length(age) == 0 | is.na(age))
    stop("Age is missing.")
  if (age < 6)
    stop("Age cannot be less than 6.")
}
