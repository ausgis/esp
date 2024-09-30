fuzzyoverlay2 = \(formula, data, method = "and"){
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

  suppressWarnings({res = purrr::map2_dfc(variable1,variable2, \(.v1,.v2) {
    dti = dplyr::select(data, dplyr::all_of(c(yname,.v1,.v2)))
    resout = sdsfun::fuzzyoverlay(paste0(yname, " ~ ."), dti,method)
    resout = as.integer(as.factor(resout))
  })})
  names(resout) = paste0("xinteract",seq_along(variable1))
  return(resout)
}
