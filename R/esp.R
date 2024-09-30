esp = \(formula, data, listw = NULL, overlay = 'and',
        discvar = NULL, discnum = 3:8, alpha = 0.75,
        bw = "AIC", adaptive = TRUE, kernel = "gaussian",
        model = 'lag', cores = 1, ...) {
  g = esp:::.gwr_hclust_disc(formula, data, listw, overlay, discvar, discnum,
                             alpha, bw, adaptive, kernel, model, cores, ...)
  sarcoef = g[[1]]
  discdf = g[[2]]
  y = g[[3]]
  xvar = g[[4]]
  allvar = g[[5]]
}

