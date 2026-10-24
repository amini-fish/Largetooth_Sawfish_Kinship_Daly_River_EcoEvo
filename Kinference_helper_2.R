hsp_histo_gg  <- function(
    kinpairs,
    lb,
    ub = NULL,
    fullsib_cut,
    bin = 5,
    HSPmean = TRUE,
    HSPdist = TRUE,
    POPmean = TRUE,
    FSPmean = TRUE,
    UPmean  = FALSE,   # often far to the left
    main = ""
){
  # --- Input checks
  if (!("PLOD" %in% names(kinpairs)))
    stop("`kinpairs` must contain a 'PLOD' column.")
  
  PLOD <- suppressWarnings(as.numeric(kinpairs$PLOD))
  PLOD <- PLOD[is.finite(PLOD)]
  if (length(PLOD) == 0) stop("No finite PLOD values found.")
  
  if (is.null(ub)) ub <- max(PLOD, na.rm = TRUE) + bin
  if (is.null(lb)) stop("Please provide `lb` (lower bound).")
  
  # --- Restrict to plotting window, set breaks/mids like base::hist
  brks <- seq(lb, ub, by = bin)
  if (length(brks) < 2) stop("Choose a larger (ub - lb) or smaller `bin`.")
  mids <- brks[-length(brks)] + diff(brks)[1] / 2
  
  df <- tibble(PLOD = PLOD) |> filter(PLOD >= lb, PLOD <= ub)
  df$bin <- cut(df$PLOD, breaks = brks, include.lowest = TRUE, right = TRUE)
  obs_counts <- as.integer(tabulate(df$bin, nbins = length(brks) - 1))
  obs_df <- tibble(mids = mids, counts = obs_counts)
  
  # --- Kinship colours (consistent with your other plot)
  kincols <- c(
    POP = "orange",
    HSP = "skyblue",
    FSP = "navyblue",
    UP  = "black"
  )
  
  # --- Which mean lines to show
  has_attr <- function(nm) !is.null(attr(kinpairs, nm))
  get_mean <- function(nm) attr(kinpairs, nm)
  
  want_means <- c(POP = POPmean, HSP = HSPmean, FSP = FSPmean, UP = UPmean)
  vlines <- bind_rows(
    lapply(names(want_means)[want_means], function(k) {
      nm <- paste0("mean_", k)
      if (has_attr(nm))
        tibble(kin = k, x = get_mean(nm), col = kincols[k])
      else NULL
    })
  )
  
  # --- Expected HSP distribution (points + line)
  hsp_exp_df <- NULL
  if (isTRUE(HSPdist) && isTRUE(HSPmean) && has_attr("mean_HSP")) {
    E.hsp <- get_mean("mean_HSP")
    sel <- df$PLOD > E.hsp & df$PLOD < fullsib_cut
    V.hsp <- mean((df$PLOD[sel] - E.hsp)^2, na.rm = TRUE)
    
    if (is.finite(V.hsp) && V.hsp > 0) {
      total_in_window <- sum(sel, na.rm = TRUE)
      exp_num <- 2 * total_in_window *
        (pnorm(brks[-1],  mean = E.hsp, sd = sqrt(V.hsp)) -
           pnorm(brks[-length(brks)], mean = E.hsp, sd = sqrt(V.hsp)))
      
      hsp_exp_df <- tibble(mids = mids, exp = exp_num)
    }
  }
  
  # --- Build plot
  p <- ggplot(df, aes(PLOD)) +
    geom_histogram(breaks = brks,
                   fill = "grey80", colour = "black") +
    labs(title = main, x = "PLOD", y = "Count") +
    scale_x_continuous(limits = c(lb, ub),
                       expand = expansion(mult = c(0, 0.01))) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 14), 
      axis.text = element_text(size = 12),
      legend.position = "none",
    )
  
  # --- Mean lines
  if (!is.null(vlines) && nrow(vlines)) {
    p <- p + geom_vline(data = vlines,
                        aes(xintercept = x, colour = kin),
                        linewidth = 1) +
      scale_colour_manual(values = kincols,
                          breaks = intersect(names(kincols), vlines$kin),
                          name = "Kin Type") +
      guides(colour = guide_legend(title = NULL))
  }
  
  # --- Expected HSP overlay
  if (!is.null(hsp_exp_df)) {
    p <- p +
      geom_point(data = hsp_exp_df, aes(x = mids, y = exp),
                 colour = kincols["HSP"], size = 2) +
      geom_line(data = hsp_exp_df, aes(x = mids, y = exp),
                colour = kincols["HSP"], linewidth = 1)
  }
  
  return(p)
}

