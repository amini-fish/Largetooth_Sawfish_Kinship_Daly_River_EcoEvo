# ------------------------------------------------------------
# Test: is the proportion of kin pairs (kin==1) higher within
#       than between communities? Uses χ² vs permutation null.
#
# pairs_df  : data.frame with columns id1, id2, kin (0/1)
# groups    : named vector id -> community  OR data.frame(id, community)
# nperm     : permutations of community labels (sizes preserved)
# ------------------------------------------------------------

setwd("C:/Users/samue/Desktop/Research/Honours")

data = read.csv("PMD_ckmrsim_complete.csv")
kinf_data = read.csv("C:/Users/samue/Desktop/Research/Honours/Chapter_2_Relatedness_paper/Kinference/kinference_all.csv")
meta = read.csv("C:/Users/samue/Desktop/Research/Honours/Kyne_Daly River Pristis pristis data_SamAminiHonours.csv")

glimpse(meta)
glimpse(data)

pairs_df = data %>%
  dplyr::select(id_1, id_2, Relatives) %>% 
  dplyr::rename(c(id1 = id_1, id2 = id_2, kin = Relatives)) %>%
  mutate(
    id1 = gsub("_", "", id1),   # remove underscore from IDs
    id2 = gsub("_", "", id2), 
    kin = as.integer(kin))

head(pairs_df)
head(meta)

kinf_data$kin <- ifelse(kinf_data$rel == "unrelated", 0, 1)
head(kinf_data)

kinf_data <- kinf_data %>%
  dplyr::mutate(
    id1 = gsub("_", "", id_1),   # remove underscore from IDs
    id2 = gsub("_", "", id_2)) %>%
  left_join(meta %>% 
              dplyr::select(Genetic_code, Year),
            by = c("id1" = "Genetic_code")) %>%
  dplyr::rename(year_1 = Year) %>%
  left_join(meta %>% dplyr::select(Genetic_code, Year),
            by = c("id2" = "Genetic_code")) %>%
  dplyr::rename(year_2 = Year)

kinf_pairs = kinf_data %>%
  dplyr::select(id_1, id_2, kin) %>% 
  dplyr::rename(c(id1 = id_1, id2 = id_2)) %>%
  mutate(
    id1 = gsub("_", "", id1),   # remove underscore from IDs
    id2 = gsub("_", "", id2), 
    kin = as.integer(kin))


kinf_2012 = kinf_data %>%
  filter(year_1 == "2012" & year_2 == "2012") %>%
  dplyr::select(id_1, id_2, kin) %>%
  dplyr::rename(c(id1 = id_1, id2 = id_2)) %>%
  mutate(
    id1 = gsub("_", "", id1),   # remove underscore from IDs
    id2 = gsub("_", "", id2), 
    kin = as.integer(kin))
  
kinf_2012
  
## 2012 only model

pairs_2012 = data %>%
  filter(year_1 == "2012" & year_2 == "2012") %>%
  dplyr::select(id_1, id_2, Relatives) %>%
  dplyr::rename(c(id1 = id_1, id2 = id_2, kin = Relatives)) %>%
  mutate(
    id1 = gsub("_", "", id1),   # remove underscore from IDs
    id2 = gsub("_", "", id2), 
    kin = as.integer(kin))


head(pairs_2012)

## 2017 only model (not worth it I reckon...)

pairs_2017 = data %>%
  filter(year_1 == "2017" & year_2 == "2017") %>%
  dplyr::select(id_1, id_2, Relatives) %>%
  dplyr::rename(c(id1 = id_1, id2 = id_2, kin = Relatives)) %>%
  mutate(
    id1 = gsub("_", "", id1),   # remove underscore from IDs
    id2 = gsub("_", "", id2), 
    kin = as.integer(kin))

head(pairs_2017)

## Load individual metadata 
meta <- read.csv("C:/Users/samue/Desktop/Research/Honours/Kyne_Daly River Pristis pristis data_SamAminiHonours.csv")

glimpse(meta)

### Billabong 
groups_df = meta %>%
  dplyr::select(Genetic_code,Site_ID_simple) %>%
  dplyr::rename(id = Genetic_code) %>%
  dplyr::rename(community = Site_ID_simple)

# OR 

### Year

groups_df = meta %>%
  dplyr::select(Genetic_code,Year) %>%
  dplyr::rename(id = Genetic_code) %>%
  dplyr::rename(community = Year)

head(groups_df)

## Run the permutatio test...

test_kin_within_vs_between_bin <- function(pairs_df,
                                           groups,
                                           nperm = 1000,
                                           seed = 42,
                                           yates_correction = FALSE) {
  stopifnot(all(c("id1","id2","kin") %in% names(pairs_df)))
  if (!all(pairs_df$kin %in% c(0,1))) stop("'kin' must be 0/1.")
  set.seed(seed)
  
  # Standardise groups to named vector
  if (is.data.frame(groups)) {
    stopifnot(all(c("id","community") %in% names(groups)))
    grp_vec <- setNames(as.character(groups$community), groups$id)
  } else {
    grp_vec <- groups
  }
  
  
  
  # Keep pairs with labels for both IDs
  keep <- pairs_df$id1 %in% names(grp_vec) & pairs_df$id2 %in% names(grp_vec)
  dat  <- pairs_df[keep, , drop = FALSE]
  if (nrow(dat) == 0) stop("No pairs with both group labels found.")
  
  g1 <- grp_vec[dat$id1]
  g2 <- grp_vec[dat$id2]
  
  within  <- g1 == g2
  is_kin  <- dat$kin == 1L
  
  # Observed 2x2 table: rows {within, between} x cols {kin, non-kin}
  obs_tab <- matrix(c(sum(within & is_kin),
                      sum(within & !is_kin),
                      sum(!within & is_kin),
                      sum(!within & !is_kin)),
                    nrow = 2, byrow = TRUE,
                    dimnames = list(c("within","between"), c("kin","nonkin")))
  
  # Observed χ² and effect (diff in proportions)
  obs_chi <- suppressWarnings(chisq.test(obs_tab, correct = yates_correction)$statistic[[1]])
  p_within  <- obs_tab["within","kin"]  / sum(obs_tab["within", ])
  p_between <- obs_tab["between","kin"] / sum(obs_tab["between", ])
  obs_diff  <- p_within - p_between
  
  # Prep for fast permutations (shuffle labels across individuals)
  uniq_ids <- sort(unique(c(dat$id1, dat$id2)))
  idx1     <- match(dat$id1, uniq_ids)
  idx2     <- match(dat$id2, uniq_ids)
  base_grp <- grp_vec[uniq_ids]
  
  chi_perm  <- numeric(nperm)
  diff_perm <- numeric(nperm)
  
  for (b in seq_len(nperm)) {
    shuf_grp <- sample(base_grp, length(base_grp), replace = FALSE)
    within_b <- shuf_grp[idx1] == shuf_grp[idx2]
    
    tab_b <- matrix(c(sum(within_b & is_kin),
                      sum(within_b & !is_kin),
                      sum(!within_b & is_kin),
                      sum(!within_b & !is_kin)),
                    nrow = 2, byrow = TRUE)
    
    chi_perm[b] <- suppressWarnings(chisq.test(tab_b, correct = yates_correction)$statistic[[1]])
    
    pw  <- tab_b[1,1] / sum(tab_b[1,])
    pb  <- tab_b[2,1] / sum(tab_b[2,])
    diff_perm[b] <- pw - pb
  }
  
  # One-sided p (alternative: within > between)
  p_chi  <- (sum(chi_perm >= obs_chi) + 1) / (nperm + 1)
  p_diff <- (sum(diff_perm >= obs_diff) + 1) / (nperm + 1)
  
  list(
    observed_table           = obs_tab,
    observed_chisq           = obs_chi,
    observed_diff_prop       = obs_diff,
    perm_chisq               = chi_perm,
    perm_diff_prop           = diff_perm,
    pvalue_chisq_one_sided   = p_chi,
    pvalue_diff_one_sided    = p_diff,
    details = list(n_pairs = nrow(dat), nperm = nperm, yates_correction = yates_correction)
  )
} 

# Ensure consistent naming and types
pairs_df <- pairs_df %>%
  dplyr::mutate(
    id1 = trimws(as.character(id1)),
    id2 = trimws(as.character(id2)),
    kin = as.integer(kin)
  )

# groups_df should already have: id (individual), billabong (community)
groups_df <- groups_df %>%
  dplyr::mutate(
    id = trimws(as.character(id)),
    community = trimws(as.character(community))
  )

head(groups_df)

# Keep only IDs that appear in your kin-pair table
ids_in_pairs <- sort(unique(c(pairs_df$id1, pairs_df$id2)))

groups_df <- groups_df %>%
  dplyr::filter(id %in% ids_in_pairs, !is.na(community), nzchar(community))

head(groups_df)

# Create named vector: id → billabong
groups_vec <- setNames(groups_df$community, groups_df$id)

##-------------------------------------------------
## Run the test
##-------------------------------------------------

# groups_df: columns id, community
res <- test_kin_within_vs_between_bin(pairs_df = pairs_2012, #pairs_2012
                                      groups   = groups_df,  # or a named vector id->community
                                      nperm    = 10000)

res$observed_table
res$observed_diff_prop # small difference of 0.01 (1%) - matches with GLMM
res$pvalue_chisq_one_sided # highly significant 0.01 (< 0.05*)

res

perm_p_ci <- function(b, m, conf = 0.95) {
  p <- (b + 1)/(m + 1)
  se <- sqrt(p*(1 - p)/m)
  z  <- qnorm( (1 + conf)/2 )
  c(p = p, lower = max(0, p - z*se), upper = min(1, p + z*se))
}

# Example: b = 0 exceedances out of m = 10000
perm_p_ci(b = 0, m = 10000)

## VIsualise the null distributions for each model 

perm_df <- data.frame(
  perm_stat = res$perm_diff_prop,  # or res$perm_chisq
  observed = res$observed_diff_prop
)

loc_plot <- ggplot(perm_df, aes(x = perm_stat)) +
  geom_density(bins = 10, fill = "grey80", colour = "black") +
  geom_vline(aes(xintercept = observed),
             colour = "red", linetype = "solid", linewidth = 1) +
  labs(x = "Permuted difference in proportion (Δp)",
       y = "Frequency",
       title = "Permutation null distribution for capture location",
       subtitle = sprintf("Observed Δp = %.3f", res$observed_diff_prop)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14))

loc_plot


# groups_df: columns id, community
res_year <- test_kin_within_vs_between_bin(pairs_df = pairs_df,
                                           groups   = groups_df,
                                           nperm    = 10000)

## VIsualise the null distributions for each model 

year_perm_df <- data.frame(
  perm_stat = res_year$perm_diff_prop,  # or res$perm_chisq
  observed = res_year$observed_diff_prop
)

year_plot <- ggplot(year_perm_df, aes(x = perm_stat)) +
  geom_density(fill = "grey80", colour = "black") +
  geom_vline(aes(xintercept = observed),
             colour = "red", linetype = "solid", linewidth = 1) +
  labs(x = "Permuted difference in proportion (Δp)",
       y = "Frequency",
       title = "Permutation null distribution for capture year",
       subtitle = sprintf("Observed Δp = %.3f", res$observed_diff_prop)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14))

print(year_plot)

library(ggpubr)

perm_combined <- ggarrange(year_plot, 
                           loc_plot, 
                           ncol = 2, 
                           nrow = 1, 
                           labels = c("A", "B"))

print(perm_combined)

devEMF::emf("C:/Users/samue/Desktop/permutation_results.emf", width = 10, height = 6)  
print(perm_combined)
dev.off()
