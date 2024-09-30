esp = \(formula, data, discvar = NULL, discnum = 3:8,
        alpha = 0.75, bw = "AIC", adaptive = TRUE,
        kernel = "gaussian", overlay = 'and', ...) {
  formula = stats::as.formula(formula)
  formula.vars = all.vars(formula)
  if (formula.vars[2] != "."){
    data = dplyr::select(data,dplyr::all_of(formula.vars))
  }
  yname = formula.vars[1]
  yvec = data[, yname, drop = TRUE]
  geom = sf::st_geometry(data)
  gdist = sdsfun::sf_distance_matrix(data)
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
  if (!is.null(xundiscname)) {
    undiscdf = data |>
      sf::st_drop_geometry() |>
      tibble::as_tibble() |>
      dplyr::select(dplyr::all_of(xundiscname)) |>
      dplyr::mutate(dplyr::across(dplyr::everything(),\(x){
        if (inherits(x,"factor")){
          x = as.integer(x)
        } else if (inherits(x,'character')) {
          x = as.integer(as.factor(x))
        }
        return(x)
      }))

    names(undiscdf) = paste0('x',seq(length(xdiscname) + 1,by = 1,
                                     length.out = length(xundiscname) + 1))
  } else {
    undiscdf = NULL
  }

  gwrcoefs = esp::gwr_betas(paste0(yname," ~ ."),discdf,bw,adaptive,kernel)
  discdf = purrr::map2_dfr(gwrcoefs, names(gwrcoefs), \(.coef,.name) {
    D0 = stats::dist(.coef)
    resh = ClustGeo::hclustgeo(D0,as.dist(gdist),alpha)
    resdisc = tibble::as_tibble(stats::cutree(resh,discnum))
    names(resdisc) = paste0("disc_",discnum)
    resdisc = dplyr::mutate(resdisc,xname = .name)
    resdisc = tibble::rowid_to_column(resdisc)
    return(resdisc)
  }) |>
    tidyr::pivot_longer(
      cols = dplyr::starts_with("disc_"),
      names_to = "discnum",
      names_prefix = "disc_",
      values_to = "disc",
      values_drop_na = TRUE
  ) |>
    dplyr::group_split(discnum) |>
    purrr::map(\(.df) {
      .res = .df |>
        dplyr::select(-discnum) |>
        tidyr::pivot_wider(names_from = xname,
                           values_from = disc) |>
        dplyr::select(-rowid)
      names(.res) = paste0('x',seq_along(.res))
      if (!is.null(undiscdf)){.res = dplyr::bind_cols(.res,undiscdf)}
      return(.res)
  })

  discdf = purrr::map(discdf, \(.df) {
    .res = sf::st_set_geometry(.df,geom)
    return(.res)
  })
  gwrcoefs = purrr::map(discdf, \(.df) {
    esp::gwr_betas("y~.",.df,bw,adaptive,kernel)
  })


  return(discdf)
}

