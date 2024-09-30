.gwr_hclust_disc = \(formula, data, listw = NULL, overlay = 'and',
                     discvar = NULL, discnum = 3:8, alpha = 0.75,
                     bw = "AIC", adaptive = TRUE, kernel = "gaussian",
                     model = 'lag', cores = 1, ...) {
  doclust = FALSE
  if (inherits(cores, "cluster")) {
    doclust = TRUE
  } else if (cores > 1) {
    doclust = TRUE
    cores = parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cores), add = TRUE)
  }

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

  if (is.null(listw)) {
    listw = spdep::nb2listw(sdsfun::spdep_nb(data), style = "W", zero.policy = TRUE)
  }

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
  gwr_hclust = \(n,discnum,alpha,...) {
    D0 = stats::dist(gwrcoefs[,n,drop = TRUE])
    resh = ClustGeo::hclustgeo(D0,stats::as.dist(gdist),alpha,...)
    resdisc = tibble::as_tibble(stats::cutree(resh,discnum))
    names(resdisc) = paste0("disc_",discnum)
    resdisc = dplyr::mutate(resdisc,xname = names(gwrcoefs)[n])
    resdisc = tibble::rowid_to_column(resdisc)
    return(resdisc)
  }

  if (doclust) {
    out_g = parallel::parLapply(cores, seq_along(gwrcoefs), gwr_hclust,
                                discnum = discnum, alpha = alpha, ...)
    discdf = tibble::as_tibble(do.call(rbind, out_g))
  } else {
    discdf = purrr::map_dfr(seq_along(gwrcoefs), gwr_hclust,
                            discnum = discnum, alpha = alpha, ...)
  }

  discdf = discdf |>
    tidyr::pivot_longer(
      cols = dplyr::starts_with("disc_"),
      names_to = "discnum",
      names_prefix = "disc_",
      values_to = "disc",
      values_drop_na = TRUE
    ) |>
    dplyr::group_split(discnum)

  bind_discdf = \(n){
    .res = discdf[[n]] |>
      dplyr::select(-discnum) |>
      tidyr::pivot_wider(names_from = xname,
                         values_from = disc) |>
      dplyr::select(-rowid)
    names(.res) = paste0('x',seq_along(.res))
    if (!is.null(undiscdf)){.res = dplyr::bind_cols(.res,undiscdf)}
    fdf = dplyr::bind_cols(tibble::tibble(y = yvec),.res)
    fdfres = esp::fuzzyoverlay2("y ~ .",fdf,overlay)[[1]]
    res = dplyr::bind_cols(.res,fdfres)
    return(res)
  }
  if (doclust) {
    discdf = parallel::parLapply(cores, seq_along(discdf), bind_discdf)
  } else {
    discdf = purrr::map(seq_along(discdf), bind_discdf)
  }

  do_dummy = \(n) {
    .df = sdsfun::dummy_tbl(discdf[[n]])
    .res = dplyr::bind_cols(tibble::tibble(y = yvec),.df)
    # .res = sf::st_set_geometry(.res,geom)
    return(.res)
  }
  if (doclust) {
    discsf = parallel::parLapply(cores, seq_along(discdf), do_dummy)
  } else {
    discsf = purrr::map(seq_along(discdf), do_dummy)
  }

  run_slm = \(n,listw,model){
    suppressWarnings({if (model == "lag") {
      g = spatialreg::lagsarlm("y ~ .",discsf[[n]],listw,Durbin =  FALSE)
    } else {
      g = spatialreg::errorsarlm("y ~ .",discsf[[n]],listw,Durbin =  FALSE)
    }})
    return(stats::coef(g))
  }
  if (doclust) {
    sarcoef = parallel::parLapply(cores, seq_along(discsf), run_slm,
                                  listw = listw, model = model)
  } else {
    sarcoef = purrr::map(seq_along(discsf), run_slm,
                        listw = listw, model = model)
  }

  res = list("coef" = sarcoef,
             "disc" = discdf,
             "dummy" = discsf)
  return(res)
}

