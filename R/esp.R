esp = \(formula, data, discvar = NULL, discnum = 10,
        overlaymethod = 'and', alpha = 0.95) {
  formula = stats::as.formula(formula)
  formula.vars = all.vars(formula)
  if (formula.vars[2] != "."){
    data = dplyr::select(data,dplyr::all_of(formula.vars))
  }
  yname = formula.vars[1]
  geom = sf::st_geometry(data)
  gname = sdsfun::sf_geometry_name(data)
  xname = colnames(data)[-which(colnames(data) %in% c(yname,gname))]

  if (is.null(discvar)) {
    xdiscname = xname
    xundiscname = NULL
  } else {
    xdiscname = discvar
    xundiscname = xname[-which(xname %in% discvar)]
  }
  discdf = dplyr::select(data,dplyr::all_of(c(yname,xdiscname)))

}

