#' @title Spatially Explicit Stratified Power (SESP) Model
#' @description
#' Equivalent geographical detector \eqn{q}-values under a spatial linear regression framework.
#' @note
#' Note that when the number of continuous independent variables is small (three or fewer), the built-in spatial explicit
#' discretization in `sesp` may overestimate the variable q value. In such cases, it is recommended to discretize the
#' variables beforehand and then input them into `sesp` for computation.
#'
#' @param formula A formula for enhanced stratified power model.
#' @param data An `sf` object of observation data. Please note that the column names of the independent
#' variables should not be `all` or `none`.
#' @param listw (optional) A `listw` object. See `spdep::mat2listw()` and `spdep::nb2listw()` for details.
#' @param discvar (optional) Name of continuous variable columns that need to be discretized. Noted that
#' when `formula` has `discvar`, `data` must have these columns. Default is `all`, which means all independent
#' variables are used as `discvar`. When `discvar` is set to `none`, all independent variables do not need to
#' be discretized.
#' @param discnum (optional) Number of discretization. Default all will use `3:8`.
#' @param model (optional) The type of linear model used, default is `ols`. The `model` value must be any of
#' `ols`, `gwr`, `lag` or `error`.
#' @param durbin (optional) Whether to consider spatial durbin terms, default is `false`.
#' @param overlay (optional) Spatial overlay method. One of `and`, `or`, `intersection`. Default is `and`.
#' @param alpha (optional) Controlling the strength of spatial soft constraints, the larger the `alpha`,
#' the stronger the spatial soft constraint. Default is `0.5`.
#' @param bw (optional) The bandwidth used in selecting models. The optimal bandwidth can be
#' selected using one of two methods: `AIC`, and `CV`. Default is `AIC`.
#' @param adaptive (optional) Whether the bandwidth value is adaptive or not. Default is `TRUE`.
#' @param kernel (optional) Kernel function. Default is `gaussian`.
#' @param increase_rate (optional) The critical increase rate of the number of discretization.
#' Default is `5%`.
#' @param cores (optional) Positive integer (default is 1). When cores are greater than 1, use
#' multi-core parallel computing.
#' @param ... (optional) Other arguments passed to `ClustGeo::hclustgeo()`.
#'
#' @return A list.
#' \describe{
#' \item{\code{factor}}{global factor detection result}
#' \item{\code{interaction}}{global interactive detection results}
#' \item{\code{optdisc}}{independent variable optimal spatial discretization}
#' \item{\code{allfactor}}{factor detection results corresponding to different numbers of discreteization}
#' \item{\code{model}}{regression model used to estimate equivalent q values}
#' }
#' @export
#'
#' @examples
#' NTDs = sf::st_as_sf(gdverse::NTDs, coords = c('X','Y'))
#' g = sesp(incidence ~ ., data = NTDs, discvar = 'none',
#'          model = 'ols', overlay = 'intersection', cores = 1)
#' g
#'
sesp = \(formula, data, listw = NULL, discvar = "all", discnum = 3:8,
         model = 'ols', durbin = FALSE, overlay = 'and', alpha = 0.5, bw = "AIC",
         adaptive = TRUE, kernel = "gaussian", increase_rate = 0.05, cores = 1, ...) {
  if (!(model %in% c("ols","gwr","lag","error"))){
    stop("`model` must be one of `ols`,`gwr`,`lag` or `error`!")
  }
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
  data = sdsfun::tbl_all2int(data)
  yname = formula.vars[1]
  yvec = data[, yname, drop = TRUE]
  geom = sf::st_geometry(data)
  gdist = sdsfun::sf_distance_matrix(data)
  gname = sdsfun::sf_geometry_name(data)
  xname = colnames(data)[-which(colnames(data) %in% c(yname,gname))]

  if (any(c("all","none") %in% xname)) {
    stop("The column names of the independent variables should not be `all` or `none`.")
  }
  if (length(xname) <= 1) {
    stop("At least two independent variables are required.")
  }

  if (is.null(listw)){
    globallw = spdep::nb2listw(sdsfun::spdep_nb(data), style = "W", zero.policy = TRUE)
  } else {
    globallw = listw
  }

  if (any(discvar == "all")){
    xdiscname = xname
    xundiscname = NULL
  } else if (any(discvar == "none")){
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

    gwrcoefs = sesp::gwr_betas(paste0(yname," ~ ."),discdf,bw,adaptive,kernel)
    gwr_hclust = \(n,discnum,alpha,...) {
      moran_g = sdsfun::moran_test(
        sf::st_set_geometry(tibble::tibble(x = gwrcoefs[,n,drop = TRUE]),geom)
      )
      moran_v = dplyr::pull(moran_g$result,2)
      moran_p = dplyr::pull(moran_g$result,6)
      se_alpha = dplyr::if_else(moran_v>=alpha&moran_p<=0.05,
                                moran_v,alpha,missing = alpha)
      D0 = stats::dist(gwrcoefs[,n,drop = TRUE])
      resh = ClustGeo::hclustgeo(D0,stats::as.dist(gdist),se_alpha,...)
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
      if (!is.null(undiscdf)){.res = dplyr::bind_cols(.res,undiscdf)}
      names(.res) = paste0('x',seq_along(.res))
      return(.res)
    }
    if (doclust) {
      discdf = parallel::parLapply(cores, seq_along(discdf), bind_discdf)
    } else {
      discdf = purrr::map(seq_along(discdf), bind_discdf)
    }
  } else {
    discnum = 0
    xvarname = xundiscname
    discdf = data |>
      sf::st_drop_geometry() |>
      dplyr::select(dplyr::all_of(xvarname))
    names(discdf) = paste0('x',seq_along(discdf))
    discdf = list(discdf)
  }

  get_slm = \(n,listw,model,durbin,bw,adaptive,kernel){
    slmvar = names(discdf[[n]])
    dummydf = sdsfun::dummy_tbl(discdf[[n]])
    slmlevelvar = names(dummydf)
    slmx = sapply(slmvar, function(x) {
      matched = grep(paste0("^", x, "_"), slmlevelvar, value = TRUE)
      res = paste(matched, collapse = "+")
      return(unname(res))
    })
    dummydf = dplyr::bind_cols(tibble::tibble(y = yvec),dummydf)

    suppressMessages({slm_res = purrr::map_dfc(slmx, \(.varname) {
      slmformula = paste0("y ~ ",.varname)
      suppressWarnings({if (model == "lag") {
        g = spatialreg::lagsarlm(slmformula, dummydf, listw,
                                 Durbin = durbin,zero.policy = TRUE)
      } else if (model == "error") {
        g = spatialreg::errorsarlm(slmformula, dummydf, listw,
                                   Durbin = durbin,zero.policy = TRUE)
      } else if (model == "gwr") {
        dummydf = sf::st_set_geometry(dummydf,geom)
        g = GWmodel3::gwr_basic(
          slmformula, dummydf, bw = bw, adaptive = adaptive, kernel = kernel
        )
      } else {
        g = stats::lm(slmformula, dummydf)
      }})

      if (model != "gwr") {
        aicv = stats::AIC(g)
        bicv = stats::AIC(g)
        loglikv = as.numeric(stats::logLik(g))
        if (model == 'error'){
          fity = as.numeric(stats::predict(g, pred.type = 'trend', listw = listw, re.form = NA))
        } else {
          fity = as.numeric(stats::predict(g, pred.type = 'TC', listw = listw, re.form = NA))
        }
        return(list("pred" = fity, "AIC" = aicv,
                    "BIC" = bicv, "LogLik" = loglikv))
      } else {
        aicv = g$diagnostic$AIC
        fity = stats::predict(g,dummydf)$yhat
        return(list("pred" = fity, "AIC" = aicv))
      }
    })})
    return(slm_res)
  }
  if (doclust) {
    slmres = parallel::parLapply(cores, seq_along(discdf), get_slm, listw = globallw, model = model,
                                 durbin = durbin, bw = bw, adaptive = adaptive, kernel = kernel)
  } else {
    slmres = purrr::map(seq_along(discdf), get_slm, listw = globallw, model = model,
                        durbin = durbin, bw = bw, adaptive = adaptive, kernel = kernel)
  }

  xinteract = utils::combn(xvarname,2,simplify = FALSE)
  variable1 = purrr::map_chr(seq_along(xinteract), \(.x) xinteract[[.x]][1])
  variable2 = purrr::map_chr(seq_along(xinteract), \(.x) xinteract[[.x]][2])
  IntersectionSymbol = rawToChar(as.raw(c(0x20, 0xE2, 0x88, 0xA9, 0x20)))
  Interactname = paste0(variable1,IntersectionSymbol,variable2)
  allvarname = c(xvarname,Interactname)

  y_pred = purrr::map(slmres, \(.df){
    .res = dplyr::select(.df,dplyr::starts_with("pred"))
    names(.res) = xvarname
    return(.res)
  })
  aicv = purrr::map(slmres, \(.df){
    .res = dplyr::slice(dplyr::select(.df,dplyr::starts_with("AIC")),1)
    names(.res) = xvarname
    return(.res)
  })
  if (model != "gwr") {
    bicv = purrr::map(slmres, \(.df){
      .res = dplyr::slice(dplyr::select(.df,dplyr::starts_with("BIC")),1)
      names(.res) = xvarname
      return(.res)
    })
    loglikv = purrr::map(slmres, \(.df){
      .res = dplyr::slice(dplyr::select(.df,dplyr::starts_with("LogLik")),1)
      names(.res) = xvarname
      return(.res)
    })
    qv = purrr::map(seq_along(y_pred),\(n){
      qvalue = SLMQ(yvec,as.matrix(y_pred[[n]]))
      resqv = tibble::tibble(Variable = names(aicv[[n]]),
                             Qvalue = qvalue,
                             AIC = as.numeric(aicv[[n]]),
                             BIC = as.numeric(bicv[[n]]),
                             LogLik = as.numeric(loglikv[[n]]),
                             DiscNum = discnum[n])
    return(resqv)
    })} else {
      qv = purrr::map(seq_along(y_pred),\(n){
        qvalue = SLMQ(yvec,as.matrix(y_pred[[n]]))
        resqv = tibble::tibble(Variable = names(aicv[[n]]),
                               Qvalue = qvalue,
                               AIC = as.numeric(aicv[[n]]),
                               DiscNum = discnum[n])
        return(resqv)
      })
  }

  fdv = purrr::list_rbind(qv)
  fdv$Variable = factor(fdv$Variable,levels = xvarname)

  if (length(y_pred) > 1){
    suppressWarnings({opt_discnum = dplyr::group_split(fdv,Variable) |>
      purrr::map_dbl(\(.df) sdsfun::loess_optnum(.df$Qvalue, .df$DiscNum,
                                                 increase_rate = increase_rate)[1])})
    opt_discdf = purrr::map_dfc(seq_along(opt_discnum),
                                \(.n) {dn = which(discnum == opt_discnum[.n])
                                return(dplyr::select(discdf[[dn]],
                                                     dplyr::all_of(paste0("x",.n))))
                                })
    opt_fdv = dplyr::group_split(fdv,Variable) |>
      purrr::map2_dfr(opt_discnum,
                      \(.qv,.discn) dplyr::filter(.qv,DiscNum == .discn)) |>
      dplyr::select(-DiscNum)
  } else {
    opt_discnum = NULL
    opt_discdf = discdf[[1]]
    opt_fdv = dplyr::select(fdv,-DiscNum)
  }

  fdf = dplyr::bind_cols(tibble::tibble(y = yvec),opt_discdf)
  fdfres = sesp::fuzzyoverlay2("y ~ .",fdf,overlay)[[1]]
  idvdf = dplyr::bind_cols(opt_discdf,fdfres)
  idvvar = names(idvdf)
  dummyidvdf = sdsfun::dummy_tbl(idvdf)
  idvlevelvar = names(dummyidvdf)
  dummyidvdf = dplyr::bind_cols(tibble::tibble(y = yvec),dummyidvdf)

  get_idv = \(svar,listw,model,durbin,bw,adaptive,kernel){
    slmx = sapply(svar, function(x) {
      matched = grep(paste0("^", x, "_"), idvlevelvar, value = TRUE)
      res = paste(matched, collapse = "+")
      return(unname(res))
    })

    suppressMessages({slm_res = purrr::map_dfc(slmx, \(.varname) {
      slmformula = paste0("y ~ ",.varname)
      suppressWarnings({if (model == "lag") {
        g = spatialreg::lagsarlm(slmformula, dummyidvdf, listw,
                                 Durbin = durbin,zero.policy = TRUE)
      } else if (model == "error") {
        g = spatialreg::errorsarlm(slmformula, dummyidvdf, listw,
                                   Durbin = durbin,zero.policy = TRUE)
      } else if (model == "gwr") {
        dummyidvdf = sf::st_set_geometry(dummyidvdf,geom)
        g = GWmodel3::gwr_basic(
          slmformula, dummyidvdf, bw = bw, adaptive = adaptive, kernel = kernel
        )
      } else {
        g = stats::lm(slmformula, dummyidvdf)
      }})

      if (model != "gwr") {
        aicv = stats::AIC(g)
        bicv = stats::AIC(g)
        loglikv = as.numeric(stats::logLik(g))
        if (model == 'error'){
          fity = as.numeric(stats::predict(g, pred.type = 'trend', listw = listw, re.form = NA))
        } else {
          fity = as.numeric(stats::predict(g, pred.type = 'TC', listw = listw, re.form = NA))
        }
        return(list("pred" = fity, "AIC" = aicv,
                    "BIC" = bicv, "LogLik" = loglikv))
      } else {
        aicv = g$diagnostic$AIC
        fity = stats::predict(g,dummyidvdf)$yhat
        return(list("pred" = fity, "AIC" = aicv))
      }
    })})
    return(slm_res)
  }
  if (doclust) {
    slmres = parallel::parLapply(cores, idvvar, get_idv, listw = globallw, model = model,
                                 durbin = durbin, bw = bw, adaptive = adaptive, kernel = kernel)
  } else {
    slmres = purrr::map(idvvar, get_idv, listw = globallw, model = model,
                        durbin = durbin, bw = bw, adaptive = adaptive, kernel = kernel)
  }

  suppressMessages({y_pred = purrr::map_dfc(slmres,
                          \(.df) dplyr::select(.df,dplyr::starts_with("pred")))})
  qv_disc = SLMQ(yvec,as.matrix(y_pred))
  names(qv_disc) = allvarname
  idtype = purrr::pmap_chr(list(qv12 = qv_disc[Interactname],
                                qv1 = qv_disc[variable1],
                                qv2 = qv_disc[variable2]),
                           InteractionType)
  opt_idv = tibble::tibble(Variable = Interactname,
                           Interaction = idtype,
                           Qv1 = qv_disc[variable1],
                           Qv2 = qv_disc[variable2],
                           Qv12 = qv_disc[Interactname],
                           Variable1 = variable1,
                           Variable2 = variable2)

  opt_fdv$Variable = as.character(opt_fdv$Variable)
  fdv$Variable = as.character(fdv$Variable)
  names(opt_discdf) = xvarname
  res = list("factor" = opt_fdv,
             "interaction" = opt_idv,
             "optdisc" = opt_discdf,
             "allfactor" = fdv,
             "model" = SLMUsed(model,durbin))
  class(res) = "sespm"
  return(res)
}

#' @title print sesp model
#' @export
#' @noRd
print.sespm = \(x, ...) {
  cat("***          Spatially Explicit Stratified Power     \n")
  cat(paste0("\n Q values are estimated using *",x$model,"* \n"))
  cat("\n -------------- Global Power of Determinant : ------------\n")
  PrintGlobalQ(utils::head(x$factor,10))
  cat("\n -------------  Global Variable Interaction : ------------\n")
  PrintGlobalQ(utils::head(dplyr::select(x$interaction,1:2),10))
  cat("\n! Only the top ten items of global scale are displayed.")
  cat("\n! The others can be accessed through specific subsets.\n")
}

#' @title plot sesp model
#' @export
#' @noRd
plot.sespm = \(x, slicenum = 2, ...) {
  g_factor = x$factor |>
    dplyr::arrange(dplyr::desc(Qvalue)) |>
    dplyr::mutate(Variable = forcats::fct_reorder(Variable, Qvalue, .desc = TRUE)) |>
    dplyr::mutate(Variable_col = c("first",rep("others",times = nrow(x$factor)-1))) |>
    dplyr::mutate(Qvtext = paste0(sprintf("%4.2f", Qvalue * 100), "%"))
  fig1 = ggplot2::ggplot(g_factor,
                         ggplot2::aes(x = Qvalue, y = Variable,
                                      fill = Variable_col)) +
    ggplot2::geom_col() +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::scale_fill_manual(breaks = c("first", "others"),
                               values = c("#DE3533","#808080")) +
    ggplot2::geom_text(data = dplyr::slice(g_factor, seq(1,slicenum)),
                       ggplot2::aes(label = Qvtext), hjust = 1.25, color = "black",
                       family = 'serif', fontface = "bold") +
    ggplot2::geom_text(data = dplyr::slice(g_factor, -seq(1,slicenum)),
                       ggplot2::aes(label = Qvtext),
                       hjust = -0.1, color = "black", fontface = "bold") +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(family = "serif",
                                                       face = "bold.italic"),
                   legend.position = "off")

  gv1 = dplyr::count(x$interaction,Variable1)
  gv2 = dplyr::count(x$interaction,Variable2)
  gv = gv1 |>
    dplyr::left_join(gv2,by = 'n') |>
    dplyr::arrange(dplyr::desc(n))
  g_interaction = dplyr::mutate(x$interaction,
                      Variable1 = factor(Variable1,levels = gv$Variable1),
                      Variable2 = factor(Variable2,levels = rev(gv$Variable2))
                  )
  fig2 = ggplot2::ggplot(g_interaction,
                         ggplot2::aes(x = Variable1, y = Variable2,
                                      size = Qv12, color = Interaction)) +
    ggplot2::geom_point(alpha = 1) +
    ggplot2::scale_size(range = c(1,10),guide = 'none') +
    ggplot2::scale_color_manual(values = c("Enhance, nonlinear" = "#EA4848",
                                           "Independent" = "#E08338",
                                           "Enhance, bi-" = "#F2C55E",
                                           "Weaken, uni-" = "#6EE9EF",
                                           "Weaken, nonlinear" = "#558DE8")) +
    ggplot2::labs(x = "", y = "", size = "", color = "") +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(family = 'serif',face = "bold.italic"),
      axis.text.x = ggplot2::element_text(family = 'serif',face = "bold.italic"),
      legend.position = "inside",
      legend.justification = c('right','bottom'),
      legend.background = ggplot2::element_rect(fill = 'transparent')
    )

  fig_p = patchwork::wrap_plots(fig1,fig2, ncol = 2)
  return(fig_p)
}
