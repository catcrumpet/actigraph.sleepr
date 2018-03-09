#' Re-integrate epochs
#'
#' Collapse post-filtered activity counts into larger epoch "buckets".
#' @param agdb A \code{tibble} (\code{tbl}) of activity data (at least) an \code{epochlength} attribute.
#' @param epoch_len_out Output (longer) epoch length in seconds, must be exact multiple of the input epoch length.
#' @param use_incomplete logical. Set to \code{TRUE} to follow ActiLife convention, which collapses all observed epochs even if they are incomplete.
#' @return A \code{tibble} (\code{tbl}) of activity data collapsed into epochs of specified length.
#' @references ActiLife 6 User's Manual by the ActiGraph Software Department. 04/03/2012.
#' @details
#' Output epochs start at the first timestamp of the input data.
#'
#' @examples
#' library("dplyr")
#' data("gtxplus1day")
#'
#' gtxplus1day %>%
#'   collapse_epochs(30)
#'
#' gtxplus1day %>%
#'   collapse_epochs(60)
#' @export

collapse_epochs <- function(agdb, epoch_len_out,
                            use_incomplete = TRUE) {

  check_args_collapse_method(agdb, epoch_len_out)
  collapse_factor <- epoch_len_out / attr(agdb, "epochlength")
  if (collapse_factor == 1) return(agdb)

  agdb <- agdb %>%
    do(collapse_epochs_(.data, epoch_len_out, collapse_factor, use_incomplete))

  attr(agdb, "epochlength") <- epoch_len_out
  attr(agdb, "epochcount") <- nrow(agdb)

  agdb
}

collapse_epochs_ <- function(data, epoch_len_out, collapse_factor, use_incomplete) {

  # Exclude lux which is summarised by `floor(mean(lux))`
  selected <- intersect(colnames(data),
                        c("axis1", "axis2", "axis3", "steps",
                          "inclineoff", "inclinestanding",
                          "inclinesitting", "inclinelying"))

  data <- data %>%
    select_at(vars("timestamp", selected)) %>%
    mutate(timestamp = first(timestamp) + seconds(floor(time_length(timestamp - first(timestamp), "second") / epoch_len_out) * epoch_len_out)) %>%
    mutate(n = 1L) %>%
    group_by(timestamp) %>%
    summarise_all(sum)

  if (!use_incomplete) {
    data <- data %>% filter(n == collapse_factor)
  }

  data %>% select(-n)
}
