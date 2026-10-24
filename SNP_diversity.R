suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(dartRverse)
  library(dartR.sexlinked)
  library(hierfstat)
  library(HardyWeinberg)
  library(devtools)
  library(ggpubr)
  library(patchwork)
})

## ---- Use your original working directory (unchanged) ----
setwd("C:/Users/samue/Desktop/Honours/analysis")

## -------------------- 7) Diversity summaries & plots (your outputs preserved) --------------------
gl <- get(load("C:/Users/samue/Desktop/PMD_clean.Rdata")); gl
gl.report.diversity(gl)

## By cohort (year_caught)
pop(gl) <- gl@other$ind.metrics$year_caught

pophet <- gl.report.heterozygosity(gl, error.bar = "SE"); pophet

## Also compute per-pop again (seems duplicated in your script; we keep behaviour)
het <- gl.report.heterozygosity(gl, error.bar = "SE")
het <- rbind(het, pophet)

## Preserve your custom rownames order
rownames(het) <- c("2011", "2012", "2017", "2018", "2024", "All")

library(forcats)
het <- het %>%
  mutate(pop = fct_recode(pop, "ALL" = "pop1"))
head(het)

## Ho
ho <- ggplot(data = het, aes(x = pop, y = Ho)) +
  geom_bar(aes(fill = pop), colour = "black", stat = "identity") +
  geom_errorbar(aes(ymin = Ho - HoSE, ymax = Ho + HoSE), width = 0.3) +
  theme_bw() +
  scale_fill_brewer(palette = "YlOrRd") +
  xlab("Cohort") + ylab("Ho") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.position = "none") +
  labs(fill = "Ho")
print(ho)

## uHe
uhe <- ggplot(data = het, aes(x = pop, y = uHe)) +
  geom_bar(aes(fill = pop), colour = "black", stat = "identity") +
  geom_errorbar(aes(ymin = uHe - uHeSE, ymax = uHe + uHeSE), width = 0.3) +
  theme_bw() +
  scale_fill_brewer(palette = "YlOrRd") +
  xlab("Cohort") + ylab("uHE") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.position = "none") +
  labs(fill = "uHe")
print(uhe)

## FIS
fis <- ggplot(data = het, aes(x = pop, y = FIS)) +
  geom_bar(aes(fill = pop), colour = "black", stat = "identity") +
  geom_errorbar(aes(ymin = FIS - FISSE, ymax = FIS + FISSE), width = 0.3) +
  theme_bw() +
  scale_fill_brewer(palette = "YlOrRd") +
  xlab("Cohort") + ylab("FIS") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none") +
  labs(fill = "FIS")
print(fis)

library(ggpubr)
het_plot <- ggarrange(ho, uhe, fis, ncol = 1, nrow = 3, align = "v")
print(het_plot)

## ---- SAVE NAME PRESERVED ----
devEMF::emf("C:/Users/samue/Desktop/het_plot.emf", width = 10, height = 12)
print(het_plot)
dev.off()

## -------------------- 8) Allelic richness (hierfstat) --------------------
gl <- get(load("C:/Users/samue/Desktop/PMD_clean.Rdata")); gl
pop(gl) <- gl@other$ind.metrics$year_caught
gi <- gl2gi(gl)
hfstat <- genind2hierfstat(gi)

ar <- allelic.richness(hfstat, diploid = TRUE)

## Your summary calcs kept (though n changed to 3464 once; we keep your code)
summary(ar$Ar)
sd(ar$Ar$`1`)/sqrt(3646)  # as in your script

ar <- as.data.frame(ar$Ar)
mean.ar <- colMeans(ar); mean.ar
sd.ar <- lapply(ar, stats::sd); sd.ar
n <- 3464
se.ar <- mapply(function(s, n) s / sqrt(n), sd.ar, n); se.ar

library(tidyverse)
ar <- ar %>%
  tibble::rownames_to_column("ID") %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "Ar")

AR_plot <- ggplot(data = ar, aes(x = Year, y = Ar, fill = as.factor(Year))) +
  geom_boxplot(colour = "black") +
  scale_fill_brewer(palette = "YlOrRd") +
  xlab("Cohort") + ylab("Allelic Richness") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none") +
  labs(fill = "Year")
print(AR_plot)

het_plot <- ggarrange(ho, uhe, fis, AR_plot,
                      ncol = 2, nrow = 2, align = "v",
                      labels = c("A","B","C","D"))
print(het_plot)

## ---- SAVE NAME PRESERVED ----
devEMF::emf("C:/Users/samue/Desktop/diversity_plot.emf", width = 10, height = 6)
print(het_plot)
dev.off()

## -------------------- 9) Optional: per-pop FIS summary --------------------
fis_results <- hierfstat::basic.stats(gi)
ind_Fis <- rowMeans(fis_results$Fis, na.rm = TRUE)  ## proxy; interpret cautiously
ind_Fis

## (inbreedR loaded in your original, but not used—keeping optional)
# install.packages("inbreedR"); library(inbreedR)

## End
