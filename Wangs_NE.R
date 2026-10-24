### Implementing Wang's Sibship Ne (2009)

## Need Q_FS #fs/n2

#Nb = 1/Q_HS + 2*Q_FS

## Load the data 

library(dplyr)
library(ggplot2)
library(ggnewscale)
library(RColorBrewer)
library(ggpubr)

setwd("C:/Users/samue/Desktop/Research/Honours/Chapter_2_Relatedness_paper")

#pairs <- read.csv("PMD_ckmrsim_complete.csv")
pairs <- read.csv("C:/Users/samue/Desktop/Research/Honours/Chapter_2_Relatedness_paper/Kinference/kinference_all.csv")

pairs <- pairs %>%
  select(id_1, id_2, rel)

glimpse(pairs)

meta <- read.csv("C:/Users/samue/Desktop/Research/Honours/Chapter_2_Relatedness_paper/Kinference/metadata.csv")

##

calc_Qs_by_year <- function(pairs, meta) {
  # Clean and standardise names
  pairs$id1 <- as.character(pairs$id_1)
  pairs$id2 <- as.character(pairs$id_2)
  pairs$rel <- as.character(pairs$rel)
  meta$id   <- as.character(meta$id)
  
  # Helper to compute QFS/QHS for one cohort
  calc_Qs <- function(ids, pairs) {
    sub_pairs <- subset(pairs, id_1 %in% ids & id_2 %in% ids)
    n <- length(unique(ids))
    tot_pairs <- choose(n, 2)
    QFS <- sum(sub_pairs$rel == "full-sibs") / tot_pairs
    QHS <- sum(sub_pairs$rel == "half-sibs") / tot_pairs
    data.frame(N = n, QFS = QFS, QHS = QHS)
  }
  
  # Loop over each unique year (not using `by()`)
  years <- sort(unique(meta$year_caught))
  out <- lapply(years, function(yr) {
    ids <- meta$id[meta$year_caught == yr]
    cbind(year = yr, calc_Qs(ids, pairs))
  })
  
  do.call(rbind, out)
}

qs_by_year <- calc_Qs_by_year(pairs, meta)
print(qs_by_year)


Nb_2012 <- 1/(0.003333333 + (2*0.02333333)); Nb_2012

qs_by_year$Nb <- 4 / (2*qs_by_year$QFS) + qs_by_year$QHS

qs_by_year

## SENb = sqrt(var(Nb)) -> 95% CI = Nb +- 1.96 x SE Nb 

qs_by_year <- within(qs_by_year, {
  P <- choose(N, 2)
  Nb <- 4 / (2*QFS) + QHS
  
  dQFS <- -2 / (QHS + 2 * QFS)^2
  dQHS <- -1 / (QHS + 2 * QFS)^2
  varQFS <- QFS * (1 - QFS) / P
  varQHS <- QHS * (1 - QHS) / P
  
  varNb <- dQFS^2 * varQFS + dQHS^2 * varQHS
  seNb  <- sqrt(varNb)
  Nb_L  <- Nb - 1.96 * seNb
  Nb_U  <- Nb + 1.96 * seNb
})

qs_by_year

qs_by_year <- qs_by_year %>%
  filter(!Nb_U %in% "NaN")

qs_by_year

Ne_overall <- length(qs_by_year$Nb) / sum(1 / qs_by_year$Nb)
Ne_overall

# Combine all cohorts into one analysis
n_total <- length(unique(c(pairs$id_1, pairs$id_2)))
tot_pairs <- choose(n_total, 2)

QFS_all <- sum(pairs$rel == "full-sibs") / tot_pairs
QHS_all <- sum(pairs$rel == "half-sibs") / tot_pairs

Ne_overall <- 4 / (2* QFS_all) + QHS_all
Ne_overall

P <- tot_pairs
dQFS <- -2 / (QHS_all + 2*QFS_all)^2
dQHS <- -1 / (QHS_all + 2*QFS_all)^2
varQFS <- QFS_all*(1-QFS_all)/P
varQHS <- QHS_all*(1-QHS_all)/P
varNe <- dQFS^2*varQFS + dQHS^2*varQHS
seNe  <- sqrt(varNe)
Ne_L  <- Ne_overall - 1.96 * seNe
Ne_U  <- Ne_overall + 1.96 * seNe

round(c(Ne = Ne_overall, SE = seNe, CI_L = Ne_L, CI_U = Ne_U), 2)
