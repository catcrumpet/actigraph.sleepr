#' Categorize activity into physical activity categories using Troiano cutpoints
#'
#' @param agdb A \code{tibble} (\code{tbl}) of activity data (at least) an \code{epochlength} attribute.
#' @param age The age of the participant in whole years. This must be a value that is greater than or equal to 6 years. The function will first check the \code{agdb} for an appropriate age attribute, if it is not available, then an age must be explicitly specified.
#' @return A \code{tibble} (\code{tbl}) of activity data. A new column \code{pa} of ordinal factor indicates the physical activity category for each epoch with labels "sed", "lig", "mod", "vig", and "ext".
#' @examples
#' library("dplyr")
#' data("gtxplus1day")
#'
#' gtxplus1day %>%
#'   add_pa_troiano()
#' @export

add_pa_troiano <- function(agdb, age = NA) {
  if (is.na(age))
    age <- attr(agdb, "age")

  age <- as.integer(age)

  check_args_pa_scores(agdb, age, "Troiano")

  attr(agdb, "pa_cutpoints") <- "Troiano"
  agdb %>% do(add_pa_troiano_(.data, age))
}

add_pa_troiano_ <- function(data, age) {
  cuts_vec <- case_when(age == 6 ~ c(1400, 3758),
                        age == 7 ~ c(1515, 3947),
                        age == 8 ~ c(1638, 4147),
                        age == 9 ~ c(1770, 4360),
                        age == 10 ~ c(1910, 4588),
                        age == 11 ~ c(2059, 4832),
                        age == 12 ~ c(2220, 5094),
                        age == 13 ~ c(2393, 5375),
                        age == 14 ~ c(2580, 5679),
                        age == 15 ~ c(2781, 6007),
                        age == 16 ~ c(3000, 6363),
                        age == 17 ~ c(3239, 6751),
                        age >= 18 ~ c(2020, 5999)) %>%
                        {c(0, 100, ., 16000, Inf)}

  data %>%
    mutate(pa = cut(axis1 * (60L / attr(data, "epochlength")),
                    breaks = cuts_vec,
                    labels = c("sed", "lig", "mod", "vig", "ext"),
                    include.lowest = TRUE,
                    right = FALSE,
                    ordered_result = TRUE))
}
