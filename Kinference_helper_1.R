library(ggplot2)
library(dplyr)

plot_histoPLOD_gg <- function(hsps,
                              UP = TRUE, HSP = TRUE, POP = TRUE, FSP = TRUE,
                              showUP = c(SPA = TRUE, Normal = FALSE),
                              main = "") {
  
  x <- hsps@bins
  y <- hsps@n_PLODs_in_bin
  binmids <- hsps@bins + (hsps@bins[2] - hsps@bins[1]) / 2
  df <- data.frame(PLOD = x[-1], logFreq = log10(pmax(head(y, -1), 0.1)))
  
  kincols <- c(UP = "black", POP = "orange", HSP = "skyblue", FSP = "navyblue")
  
  p <- ggplot(df, aes(PLOD, logFreq)) +
    geom_step(colour = "black", linewidth = 0.7) +
    labs(title = main, x = "PLOD", y = expression(log[10]*"(Frequency)")) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 14), 
      axis.text = element_text(size = 12),
      legend.position = "none",
    ) +
    scale_y_continuous(limits = c(-1,3.5))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank())
  
  shown_kinships <- c()
  if (UP)  shown_kinships <- c(shown_kinships, "UP")
  if (HSP) shown_kinships <- c(shown_kinships, "HSP")
  if (POP) shown_kinships <- c(shown_kinships, "POP")
  if (FSP) shown_kinships <- c(shown_kinships, "FSP")
  
  for (kin in shown_kinships) {
    mean_val <- attr(hsps, paste0("mean_", kin))
    if (!is.null(mean_val))
      p <- p + geom_vline(xintercept = mean_val, colour = kincols[kin], linewidth = 1)
  }
  
  if (isTRUE(showUP["SPA"])) {
    spa_y <- log10(diff(hsps@binprobs) * sum(y[binmids < 0]))
    spa_df <- data.frame(x = binmids, logFreq = spa_y)
    p <- p + geom_line(data = spa_df, aes(x, logFreq),
                       colour = kincols["UP"], linewidth = 1)
  }
  
  if (isTRUE(showUP["Normal"])) {
    norm_y <- log10(diff(c(0, pnorm(binmids,
                                    mean = hsps@mean_UP,
                                    sd = sqrt(hsps@var_UP)) * sum(y[binmids < 0]))))
    norm_df <- data.frame(x = x, logFreq = norm_y)
    p <- p + geom_line(data = norm_df, aes(x, logFreq),
                       colour = kincols["UP"], linetype = "dashed", linewidth = 1)
  }
  
  # single plot output (no second '+ p')
  return(p)
}

