library(ggplot2)
library(dplyr)
library(tidyr)

autopick_threshold_gg <- function(
    x,
    kin,
    fitrange_PLOD,
    FPtol_pairs,
    use4th,
    selecto = c("ML","paranoid"),
    NVAR = 10,
    plot_bins = NULL,
    shading_density = 10,     # kept for signature parity; not used in ggplot shading
    want_all_results = FALSE,
    ...                       # ignored for plotting; ggplot handles theme/labels internally
){
  stopifnot(length(fitrange_PLOD) == 2, all(is.finite(fitrange_PLOD)))
  fitrange_PLOD <- sort(fitrange_PLOD)
  min_PLOD <- fitrange_PLOD[1]
  max_PLOD <- fitrange_PLOD[2]
  
  li <- if (x %is.a% "snpgeno") x$locinfo else x
  
  E_HSP <- sum(li$E_HSP)
  E_UP  <- sum(li$E_UP)
  
  E2 <- E_HSP
  E3 <- (E_UP + E_HSP) / 2
  E4 <- (3*E_UP + 1*E_HSP) / 4
  
  if (E3 > min_PLOD) {
    warning(sprintf('fitrange_PLOD goes below 3rd-order mean, which is %5.2f; probably a bad idea', E3))
  }
  if (E2 > max_PLOD) {
    stop(sprintf('fitrange_PLOD exceeds 2nd-order mean, which is %5.2f; noooo!', E2))
  }
  
  kin_PLOD  <- kin$PLOD %such.that% (. %in.range% fitrange_PLOD)
  overmean  <- kin_PLOD %such.that% (. > E2)
  emp_V_HSP <- mean(sqr(overmean - E2))
  
  vpk <- var_PLOD_kin(x, emp_V_HSP, n_meio = c(3,4))
  
  CDF <- function(plod, P234, SD234){
    P234[1]*pnorm(plod, mean=E2, sd=SD234[1]) +
      P234[2]*pnorm(plod, mean=E3, sd=SD234[2]) +
      P234[3]*pnorm(plod, mean=E4, sd=SD234[3])
  }
  mixlglk3 <- function(P3, SD234){
    P234 <<- c(1-P3, P3, 0)
    pdf  <<- (1-P3)*pdf2 + P3*pdf3
    Pr_in_range <<- diff(CDF(fitrange_PLOD, P234, SD234))
    pdf <- pdf / Pr_in_range
    sum(log(pdf))
  }
  mixlglk4 <- function(Ppar, SD234){
    P3 <<- Ppar[1]
    P4 <<- (1-P3) * Ppar[2]
    P2 <<- 1 - P3 - P4
    pdf <<- P2*pdf2 + P3*pdf3 + P4*pdf4
    Pr_in_range <<- diff(CDF(fitrange_PLOD, c(P2,P3,P4), SD234))
    pdf <- pdf / Pr_in_range
    sum(log(pdf))
  }
  
  set_stuff <- function(nlocal = sys.parent()) mlocal({
    Nall    <- Nobs / Pr_in_range
    N3      <- Nall * P3
    thresh  <- qnorm(FPtol_pairs / N3, mean = E3, sd = SD3, lower.tail = FALSE)
    Pr_FNeg <- pnorm(thresh, mean = E2, sd = SD2)
  })
  
  Nobs <- length(kin_PLOD)
  SD2  <- sqrt(emp_V_HSP)
  pdf2 <- dnorm(kin_PLOD, mean = E2, sd = SD2)
  pdf  <- 0*pdf2
  Pr_in_range <- P2 <- P3 <- P4 <- (-999)
  
  V3_range <- sort(vpk[,"M3"])
  SD3i <- seq(from = sqrt(V3_range[1]), to = sqrt(V3_range[2]), length = NVAR)
  
  V4_range <- sort(vpk[,"M4"])
  SD4i <- seq(from = sqrt(V4_range[1]), to = sqrt(V4_range[2]), length = NVAR)
  
  allvals3 <- allvals4 <- NULL
  
  for (ivar in seq_along(SD3i)) {
    SD3 <- SD3i[ivar]
    SD4 <- SD4i[ivar]
    this_SD234 <- c(SD2, SD3, SD4)
    
    pdf3 <- dnorm(kin_PLOD, mean = E3, sd = SD3)
    pdf4 <- dnorm(kin_PLOD, mean = E4, sd = SD4)
    
    silly <- c(0.01, 1-0.01)
    bestio <- optimize(mixlglk3, silly, maximum = TRUE, SD234 = this_SD234)
    P3 <- bestio$maximum
    if (min(abs(P3 - silly)) < 0.01) {
      warning(sprintf('Silly "best" fit for SD3=%5.1f', SD3))
    }
    lglk <- mixlglk3(P3, this_SD234)
    set_stuff()
    
    newvals3 <- returnList(SD2, SD3, SD4 = 0.1, P2 = 1-P3, P3, P4 = 0, Nall, lglk, thresh, Pr_FNeg)
    allvals3 <- if (is.null(allvals3)) data.frame(newvals3) else rbind(allvals3, newvals3)
    
    if (use4th) {
      startio <- c(0.1, 0.7)
      bestio <- nlminb(startio, NEG(mixlglk4),
                       lower = c(0, 1.5/(1+1.5)), upper = c(1, 5/(1+5)),
                       SD234 = this_SD234)
      lglk <- mixlglk4(bestio$par, this_SD234)
      set_stuff()
      newvals4 <- returnList(SD2, SD3, SD4, P2, P3, P4, N3, Nall, lglk, thresh, Pr_FNeg)
      allvals4 <- if (is.null(allvals4)) data.frame(newvals4) else rbind(allvals4, newvals4)
    }
  }
  
  selecto <- match.arg(selecto)
  allvals  <- if (use4th) allvals4 else allvals3
  i_highest <- which.max(allvals$thresh)
  i_fittest <- which.max(allvals$lglk)
  picki     <- if (selecto == "paranoid") i_highest else i_fittest
  
  # --------------------- ggplot "Graph?" section ---------------------
  if (!is.null(plot_bins)) {
    # colours to match your scheme elsewhere
    kincols <- c(HSP = "skyblue", HTP = "grey40", HCP = "orange")
    
    # histogram breaks/mids like base code
    brks <- seq(from = min(kin$PLOD, na.rm = TRUE),
                to   = max_PLOD, by = plot_bins)
    mids <- brks[-length(brks)] + diff(brks)[1]/2
    
    df_hist <- tibble(PLOD = kin$PLOD) |> filter(PLOD < max_PLOD)
    
    vals <- as.list(allvals[picki, ])
    SD2 <- vals$SD2; SD3 <- vals$SD3; SD4 <- vals$SD4
    P2  <- vals$P2;  P3  <- vals$P3;  P4  <- if (isTRUE(use4th)) vals$P4 else 0
    Nall <- vals$Nall; thresh <- vals$thresh
    
    # composite counts
    comp_counts <- Nall * diff(
      (P2 * pnorm(brks, mean = E2, sd = SD2)) +
        (P3 * pnorm(brks, mean = E3, sd = SD3)) +
        (P4 * pnorm(brks, mean = E4, sd = SD4))
    )
    comp_df <- tibble(mids = mids, count = comp_counts, what = "Overall")
    
    # components (HSP = 2nd order, HTP = 3rd order, HCP = 4th if used)
    comp_list <- list(
      HSP = diff(pnorm(brks, mean = E2, sd = SD2)) * Nall * P2,
      HTP = diff(pnorm(brks, mean = E3, sd = SD3)) * Nall * P3
    )
    if (isTRUE(use4th)) comp_list$HCP <- diff(pnorm(brks, mean = E4, sd = SD4)) * Nall * P4
    
    comp_parts <- bind_rows(lapply(names(comp_list), function(nm) {
      tibble(mids = mids, count = comp_list[[nm]], what = nm)
    }))
    
    shade_df <- tibble(xmin = -Inf, xmax = min_PLOD, ymin = -Inf, ymax = Inf)
    
    p <- ggplot(df_hist, aes(PLOD)) +
      geom_histogram(breaks = brks, fill = "grey85", colour = "black") +
      geom_rect(data = shade_df, inherit.aes = FALSE,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "white", alpha = 0.7) +
      geom_vline(xintercept = thresh, linetype = 2, linewidth = 1, colour = "black") +
      geom_line(data = comp_df, aes(x = mids, y = count), linewidth = 1.2, colour = "black") +
      geom_line(data = comp_parts,
                aes(x = mids, y = count, colour = what),
                linetype = 2, linewidth = 1) +
      scale_colour_manual(values = kincols,
                          breaks = intersect(names(kincols), comp_parts$what),
                          name = NULL) +
      labs(
        title = sprintf("Autothresh: %s, #FP = %.1f", selecto, FPtol_pairs),
        x = "PLOD", y = "Count"
      ) +
      # >>> add small right margin so last label isn't clipped <<<
      scale_x_continuous(
        limits = c(min(kin$PLOD, na.rm = TRUE), max_PLOD),
        expand = expansion(mult = c(0, 0.05))
      ) +
      theme_bw()+
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12)
      ) +
      guides(colour = guide_legend(override.aes = list(linetype = 2)))
    
    #print(p)
  }
  # -------------------------------------------------------------------
  
  # Return the threshold object just like original
  flatto <- function(matto){
    m <- c(matto)
    names(m) <- outer(rownames(matto), colnames(matto),
                      function(x,y) sprintf('%s_%s', x, y))
    m
  }
  
  threshold <- allvals$thresh[picki]
  attributes(threshold) <- c(
    returnList(fitrange_PLOD, selecto, use4th, FPtol_pairs),
    info = list(c(
      flatto(vpk[,-1]),
      vpk@info %without.name% "n_meio",
      unlist(returnList(E2, E3, E4)),
      unlist(allvals[picki, ] %without.name% "lglk")
    ))
  )
  
  if (want_all_results){
    MAX_lglk <- max(allvals$lglk)
    allvals$lglk  <- allvals$lglk  - MAX_lglk
    allvals3$lglk <- allvals3$lglk - MAX_lglk
    if (use4th) allvals4$lglk <- allvals4$lglk - MAX_lglk
    threshold@allvals3 <- allvals3
    threshold@allvals4 <- allvals4
  }
  if (exists("p")) attr(threshold, "plot") <- p
  
  return(threshold)
}
