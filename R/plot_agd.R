#' Plot activity values
#'
#' Plot a time series of activity values (by default, the counts on the vertical axis \emph{axis1}).
#' @param agdb A \code{tibble} (\code{tbl}) of activity data (at least) an \code{epochlength} attribute.
#' @param var The activity variable to plot. The default is axis1.
#' @param color Activity line color.
#' @param nrow,ncol Number of rows and columns. Relevant only if the activity data is grouped.
#' @examples
#' library("dplyr")
#' data("gtxplus1day")
#' data <- gtxplus1day %>%
#'   collapse_epochs(60) %>%
#'   apply_cole_kripke()
#'
#' plot_activity(data, color = "gray")
#' plot_activity(data, color = "sleep")
#' @export
plot_activity <- function(agdb, var = "axis1", color = "black",
                          nrow = NULL, ncol = NULL) {
  selected_vars <- intersect(names(agdb), color)
  if (length(selected_vars)) is.var <- TRUE else is.var <- FALSE
  if (is.var) {
    p <- ggplot(agdb, aes_string("timestamp", var,
                                 color = color, fill = color)) +
      geom_bar(stat = "identity")
  } else {
    p <- ggplot(agdb, aes_string("timestamp", var)) +
      geom_bar(stat = "identity", color = color, fill = color)
  }
  if (is.grouped_df(agdb)) {
    p <- p +
      facet_wrap(as.character(groups(agdb)),
                 nrow = nrow, ncol = ncol,
                 scales = "free_x")
  }
  p + theme_light() + labs(x = "time", y = var) + guides(fill = "none")
}
#' Plot activity and periods
#'
#' Plot activity values as a time series and periods as polygons.
#' @inheritParams plot_activity
#' @param periods A \code{tibble} of periods with at least two columns \code{start_var} and \code{end_var}.
#' @param act_var The activity variable to plot. The default is axis1.
#' @param start_var The column which indicates when the time periods start.
#' @param end_var The column which indicates when the time periods end.
#' @param fill Polygon fill color.
#' @examples
#' library("dplyr")
#' library("lubridate")
#' data("gtxplus1day")
#'
#' # Detect sleep periods using Sadeh as the sleep/awake algorithm
#' # and Tudor-Locke as the sleep period algorithm
#' periods_sleep <- gtxplus1day %>%
#'   collapse_epochs(60) %>%
#'   apply_cole_kripke() %>%
#'   apply_tudor_locke(min_sleep_period = 60)
#'
#' plot_activity_period(gtxplus1day, periods_sleep, "axis1",
#'                      "in_bed_time", "out_bed_time")
#' @export
plot_activity_period <- function(agdb, periods, act_var,
                                 start_var, end_var,
                                 color = "black", fill = "#525252",
                                 ncol = NULL, nrow = NULL) {
  stopifnot(identical(groups(agdb), groups(periods)))
  height <- max(agdb[[act_var]], 1000)
  polys <- periods %>%
    # Convert the timestamps to strings and back to datetime objects
    # to avoid error about "attributes not identical across variables"
    mutate_at(c(start_var, end_var),
              function(x) strftime(x, tz = tz(x), usetz = TRUE)) %>%
    construct_period_polys(start_var, end_var, height) %>%
    mutate_at(vars(x), ymd_hms)
  p <- plot_activity(agdb, var = act_var, color = color,
                     nrow = nrow, ncol = ncol)
  if (nrow(polys)) {
    p <- p +
      geom_polygon(data = polys, aes(x = x, y = y, group = id),
                   fill = fill, alpha = 0.2,
                   show.legend = FALSE, inherit.aes = FALSE)
  }
  p
}
#' Transform periods to polygons
#'
#' Transform periods to polygons, so that they can be overlaid on an activity plot.
#' @inheritParams plot_activity_period
#' @param height The height of the rectangles.
#' @export
construct_period_polys <- function(periods, start_var, end_var, height) {
  periods %>%
    do(construct_period_polys_(., start_var, end_var, height))
}
construct_period_polys_ <- function(periods, start_var, end_var, height) {
  nr <- nrow(periods)
  df <- periods %>%
    mutate(id = row_number()) %>%
    select_("id", start_var, end_var) %>%
    gather(key, x, - id) %>% select(- key) %>%
    arrange(id, x)
  df[rep(seq_len(2 * nr), each = 2), ] %>%
    mutate(y = rep(c(height, 0, 0, height), times = nr))
}