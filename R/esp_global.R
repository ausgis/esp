#' @title global stratified power
#' @description
#' Equivalent geographical detector q-statistic under a spatial linear regression framework.
#'
#' @param formula A formula
#' @param data An `sf` object of observation data. Please note that the column names of the independent
#' variables should not be `all` or `none`.
#' @param listw (optional) A `listw`. See `spdep::mat2listw()` and `spdep::nb2listw()` for details.
#' @param discvar (optional) Name of continuous variable columns that need to be discretized. Noted that
#' when `formula` has `discvar`, `data` must have these columns. Default is `all`, which means all independent
#' variables are used as `discvar`. When `discvar` is set to `none`, all independent variables do not need to
#' be discretized.
#' @param discnum (optional) Number of discretization. Default all will use `3:8`.
#' @param model (optional) The type of linear model used, with the spatial lag model as default.
#' @param Durbin (optional) Whether to consider spatial Durbin terms, default is `false`.
#' @param overlay (optional) Spatial overlay method. One of `and`, `or`, `intersection`. Default is `and`.
#' @param alpha (optional) Controlling the strength of spatial soft constraints, the larger the `alpha`,
#' the stronger the spatial soft constraint. Default is `0.75`.
#' @param bw (optional) The bandwidth used in selecting models. The optimal bandwidth can be
#' selected using one of two methods: `AIC`, and `CV`. Default will use `AIC`.
#' @param adaptive (optional) Whether the bandwidth value is adaptive or not. Default is `TRUE`.
#' @param kernel (optional) Kernel function. Default is `gaussian`.
#' @param cores (optional) Positive integer (default is 1). When cores are greater than 1, use
#' multi-core parallel computing.
#' @param ... (optional) Other arguments passed to `ClustGeo::hclustgeo()`.
#'
#' @return A list.
#' \describe{
#' \item{\code{factor}}{factor detection result}
#' \item{\code{interaction}}{Interactive detection results}
#' \item{\code{pred}}{predicted values of response variables required for the calculation of q-values}
#' \item{\code{disc}}{spatial discretization results}
#' \item{\code{y}}{values corresponding to the response variable}
#' \item{\code{xvar}}{all the explanatory variables}
#' \item{\code{allvar}}{all explanatory variables and their pairwise interaction variables.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' depression = system.file('extdata/Depression.csv',package = 'gdverse') |>
#'   readr::read_csv() |>
#'   sf::st_as_sf(coords = c('X','Y'), crs = 4326)
#' g = esp_global(Depression_prevelence ~ ., data = depression, cores = 6)
#' g$factor
#' }
esp_global = \(formula, data, listw = NULL, discvar = "all", discnum = 3:8,
               model = 'lag', Durbin = FALSE, overlay = 'and', alpha = 0.75,
               bw = "AIC", adaptive = TRUE, kernel = "gaussian", cores = 1, ...) {
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

  if (c("all","none") %in% xname) {
    stop("The column names of the independent variables should not be `all` or `none`.")
  }

  if (is.null(listw)) {
    listw = spdep::nb2listw(sdsfun::spdep_nb(data), style = "W", zero.policy = TRUE)
  }

  if (discvar == "all"){
    xdiscname = xname
    xundiscname = NULL
  } else if (discvar == "none"){
    xdiscname = NULL
    xundiscname = xname
  } else {
    xdiscname = discvar
    xundiscname = xname[-which(xname %in% discvar)]
  }

  if (!is.null(xdiscname)) {
    discdf = dplyr::select(data,dplyr::all_of(c(yname,xdiscname)))
    if (!is.null(xundiscname)) {
      xvarname = c(xdiscname,xundiscname)
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

      names(undiscdf) = paste0('x',seq(length(xdiscname) + 1, by = 1,
                                       length.out = length(xundiscname) + 1))
    } else {
      xvarname = xdiscname
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
  } else {
    discnum = 0
    xvarname = xundiscname
    discdf = dplyr::select(data,dplyr::all_of(c(yname,xvarname)))
  }

  get_slm = \(n,listw,model,Durbin){
    slmvar = names(discdf[[n]])
    dummydf = sdsfun::dummy_tbl(discdf[[n]])
    slmlevelvar = names(dummydf)
    slmx = sapply(slmvar, function(x) {
      matched = grep(paste0("^", x), slmlevelvar, value = TRUE)
      res = paste(matched, collapse = "+")
      return(unname(res))
    })
    dummydf = dplyr::bind_cols(tibble::tibble(y = yvec),dummydf)

    suppressMessages({slm_res = purrr::map_dfc(slmx, \(.varname) {
      slmformula = paste0("y ~ ",.varname)
      suppressWarnings({if (model == "lag") {
        g = spatialreg::lagsarlm(slmformula, dummydf, listw,
                                 Durbin = Durbin,zero.policy = TRUE)
      } else if (model == "error") {
        g = spatialreg::errorsarlm(slmformula, dummydf, listw,
                                   Durbin = Durbin,zero.policy = TRUE)
      } else {
        g = stats::lm(slmformula, dummydf)
      }})

      aicv = stats::AIC(g)
      loglikv = as.numeric(stats::logLik(g))
      if (model == "ols") {
        fity = g$fitted.values
        g = summary(g)
        pv = suppressWarnings(stats::pf(g$fstatistic[1],g$fstatistic[2],
                                        g$fstatistic[3],lower.tail = FALSE))
      } else {
        fity = as.numeric(stats::predict(g, pred.type = 'TC', listw = listw, re.form = NA))
        g = summary(g)
        pv = as.numeric(g$LR1$p.value)
      }
      return(list("pred" = fity, "pvalue" = pv,
                  "AIC" = aicv, "LogLik" = loglikv))
    })})
    return(slm_res)
  }
  if (doclust) {
    slmres = parallel::parLapply(cores, seq_along(discdf), get_slm, listw = listw,
                                 model = model, Durbin = Durbin)
  } else {
    slmres = purrr::map(seq_along(discdf), get_slm, listw = listw,
                        model = model, Durbin = Durbin)
  }

  xinteract = utils::combn(xvarname,2,simplify = FALSE)
  variable1 = purrr::map_chr(seq_along(xinteract), \(.x) xinteract[[.x]][1])
  variable2 = purrr::map_chr(seq_along(xinteract), \(.x) xinteract[[.x]][2])
  IntersectionSymbol = rawToChar(as.raw(c(0x20, 0xE2, 0x88, 0xA9, 0x20)))
  Interactname = paste0(variable1,IntersectionSymbol,variable2)
  allvarname = c(xvarname,Interactname)

  aicv = purrr::map(slmres, \(.df){
    .res = dplyr::slice(dplyr::select(.df,dplyr::starts_with("AIC")),1)
    names(.res) = allvarname
    return(.res)
  })
  loglikv = purrr::map(slmres, \(.df){
    .res = dplyr::slice(dplyr::select(.df,dplyr::starts_with("LogLik")),1)
    names(.res) = allvarname
    return(.res)
  })
  y_pred = purrr::map(slmres, \(.df){
    .res = dplyr::select(.df,dplyr::starts_with("pred"))
    names(.res) = allvarname
    return(.res)
  })
  pv = purrr::map(slmres, \(.df){
    .res = dplyr::slice(dplyr::select(.df,dplyr::starts_with("pvalue")),1)
    names(.res) = allvarname
    return(.res)
  })
  qv = purrr::map(seq_along(y_pred),\(n){
    qvalue = SLMQ(as.matrix(y_pred[[n]]),yvec)
    resqv = tibble::tibble(Variable = names(pv[[n]]),
                           Qvalue = qvalue,
                           Pvalue = as.numeric(pv[[n]]),
                           AIC = as.numeric(aicv[[n]]),
                           LogLik = as.numeric(loglikv[[n]]),
                           DiscNum = discnum[n])
    return(resqv)
  })

  fdv = purrr::map_dfr(seq_along(qv),\(n) qv[[n]][seq_along(xvarname),])
  idv = purrr::map_dfr(seq_along(qv),\(n) {
    qv_disc = qv[[n]][,"Qvalue",drop = TRUE]
    names(qv_disc) = allvarname
    idtype = purrr::pmap_chr(list(qv12 = qv_disc[Interactname],
                                  qv1 = qv_disc[variable1],
                                  qv2 = qv_disc[variable2]),
                             InteractionType)
    return(tibble::tibble(Variable = Interactname,
                          Interaction = idtype,
                          Qv1 = qv_disc[variable1],
                          Qv2 = qv_disc[variable2],
                          Qv12 = qv_disc[Interactname],
                          DiscNum = discnum[n]))
  })

  res = list("factor" = fdv,
             "interaction" = idv,
             "pred" = y_pred,
             "disc" = discdf,
             "y" = yvec,
             "xvar" = xvarname,
             "allvar" = allvarname)
  return(res)
}

