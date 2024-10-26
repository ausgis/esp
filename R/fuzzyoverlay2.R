#' Spatial fuzzy overlay between variables pairwise
#'
#' @param formula A formula.
#' @param data A `data.frame` or `tibble` of discretized data.
#' @param method (optional) Spatial overlay method. One of `and`, `or`, `intersection`.
#' Default is `and`.
#'
#' @return A list
#' \describe{
#' \item{\code{overlay}}{overlay results between pairs of variables}
#' \item{\code{variable}}{pairwise interacting variable}
#' }
#' @export
#'
#' @examples
#' sim = tibble::tibble(y = stats::runif(7,0,10),
#'                      x1 = c(1,rep(2,3),rep(3,3)),
#'                      x2 = c(rep(1,2),rep(2,2),rep(3,3)),
#'                      x3 = c(rep(1,3),rep(2,2),rep(3,2)))
#' fo1 = fuzzyoverlay2(y ~ .,data = sim, method = 'and')
#' fo1
#' fo2 = fuzzyoverlay2(y ~ .,data = sim, method = 'or')
#' fo2
#' fo3 = fuzzyoverlay2(y ~ .,data = sim, method = 'intersection')
#' fo3
#'
fuzzyoverlay2 = \(formula, data, method = "and"){
  if (!(method %in% c("and","or","intersection"))){
    stop("`method` must be one of `and`,`or` or `intersection`!")
  }
  formula = stats::as.formula(formula)
  formula.vars = all.vars(formula)
  if (formula.vars[2] != "."){
    data = dplyr::select(data,dplyr::all_of(formula.vars))
  }
  yname = formula.vars[1]
  xname = colnames(data)[-which(colnames(data) == yname)]
  xinteract = utils::combn(xname,2,simplify = FALSE)
  variable1 = purrr::map_chr(seq_along(xinteract), \(.x) xinteract[[.x]][1])
  variable2 = purrr::map_chr(seq_along(xinteract), \(.x) xinteract[[.x]][2])

  if (method == "intersection"){
    suppressMessages({res = purrr::map2_dfc(variable1,variable2, \(.v1,.v2) {
      dti = dplyr::select(data, dplyr::all_of(c(.v1,.v2)))
      resout = purrr::reduce(dti,paste,sep = '_')
      resout = as.integer(as.factor(resout))
    })})
  } else {
    suppressMessages({res = purrr::map2_dfc(variable1,variable2, \(.v1,.v2) {
      dti = dplyr::select(data, dplyr::all_of(c(yname,.v1,.v2)))
      resout = sdsfun::fuzzyoverlay(paste0(yname, " ~ ."), dti,method)
      resout = as.integer(as.factor(resout))
    })})
  }
  names(res) = paste0("xi",seq_along(variable1))
  IntersectionSymbol = rawToChar(as.raw(c(0x20, 0xE2, 0x88, 0xA9, 0x20)))
  variable = paste0(variable1,IntersectionSymbol,variable2)
  res = list("overlay" = res, "variable" = variable)
  return(res)
}
