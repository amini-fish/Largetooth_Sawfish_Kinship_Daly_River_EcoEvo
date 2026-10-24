setwd("C:/Users/samue/Desktop/Honours/analysis")

### INSTALL WHERE NEEDED (assumed starting from new PC for reporducing) ### 

install.packages("hierfstat")
devtools::install_version("ggplot2", "3.4.4")
#install.packages("related")
install.packages("related", repos="http://R-Forge.R-project.org")
install.packages("parrallel")

### LOAD DESIRED PACKAGES INTO R ###

library(dartRverse)
library(ggplot2)
library(hierfstat)
library(dplyr)
library(devtools)
library(viridis)
library(related)
library(ggpubr)

## Not run: 
#---Read data into R---#

gl <- get(load("C:/Users/samue/Desktop/Honours/PMD_clean.Rdata")); gl

input <- dartR.base::gl2related(gl, save = F)

#---Calculate Relatedness---#

output <- coancestry(genotype.data = input, 
                     trioml = 1,
                     allow.inbreeding = T,
                     rng.seed = 42)

#---View Point Estimates---#
output$relatedness

## End(Not run)	

### Save trioML results into a new dataframe

trioML <- data.frame(output$relatedness[, c(2,3,5)])

mean <- mean(trioML$trioml)
sd <- sd(trioML$trioml)
med <- median(trioML$trioml)

trioML$rel <- trioML$trioml/2
summary(trioML$rel)

head(trioML)



rel.hist <- ggplot(trioML, aes(x = as.numeric(rel))) +
  geom_histogram(
    bins = 50,
    colour = "black",
    fill = "darkgrey"
  ) +
  geom_vline(xintercept = mean, linewidth = 1.1, colour = "green") +
  geom_vline(xintercept = c(0.125, 0.25),
             linewidth = 1.1,
             colour = c("skyblue", "navyblue")) +
  geom_vline(xintercept = c(0.092, 0.158, 0.204, 0.296),
             colour = "black", linewidth = 1, linetype = "dashed") +
  scale_x_continuous(
    name = "Relatedness",
    limits = c(-0.01, 0.325),
    breaks = seq(0.0, 0.325, by = 0.05),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 1000),
    expand = c(0, 0),
    breaks = seq(0, 1000, by = 200)
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14)
  ) 

print(rel.hist)

devEMF::emf("C:/Users/samue/Desktop/tioML_hist.emf", width = 10, height = 6)
print(rel.hist)
dev.off()


## Extract the indiviudals with relatedness above x and y -  need to find a solid guideline for this assignment not arbitrary numbers.


ind_metadata <- read.csv("C:/Users/samue/Desktop/10658_prelim_report/Pristis_metadata_SRA.csv")

glimpse(ind_metadata)

ind_metadata <- ind_metadata %>%
  mutate(year = year(dmy(Date)))

colnames(trioML) <- c("id_1", "id_2", "rel")

trio_rel <- trioML %>%
  left_join(ind_metadata, by = c("id_1" = "id"))

names(trio_rel)

# 1) Build pre-suffixed copies of the metadata
meta1 <- ind_metadata %>%
  transmute(
    id,
    birth_year_1  = birth_year,
    billabong_1   = Site_ID,
    capture_year1 = year_caught
  )

meta2 <- ind_metadata %>%
  transmute(
    id,
    birth_year_2  = birth_year,
    billabong_2   = Site_ID,
    capture_year2 = year_caught
  )

head(trioML)

# 2) Join without any post-join renaming
trio_rel2 <- trioML %>%
  left_join(meta1, by = c("id_1" = "id")) %>%
  left_join(meta2, by = c("id_2" = "id")) %>%
  mutate(
    birth_year_diff  = abs(birth_year_1 - birth_year_2),
    Within_Billabong = as.integer(billabong_1 == billabong_2),
    Within_Cohort    = as.integer(capture_year1 == capture_year2)
  ) %>%
  # NOTE: your object has 'rel' (not 'value'), per names(trio_rel)
  dplyr::select(
    id_1, id_2, rel,
    billabong_1, billabong_2,
    birth_year_1, birth_year_2,
    capture_year1, capture_year2,
    birth_year_diff, Within_Billabong, Within_Cohort
  )

head(trio_rel2)


trio_rel3 <- trio_rel2 %>%
  unite("year_caught_both",
        as.character("capture_year1"),as.character("capture_year2"), 
        sep = ":",
        remove = FALSE, na.rm = FALSE) %>%
  mutate(year_caught_both=factor(year_caught_both))

## Do the same for Billabong

trio_rel3 <- trio_rel3 %>% 
  unite("billabong_both",
        as.character("billabong_1"),as.character("billabong_2"), 
        sep = ":",
        remove = FALSE, na.rm = FALSE) %>%
  mutate(billabong_both=factor(billabong_both))


head(trio_rel3)

levels(trio_rel3$billabong_both)

## Read to run the boxplot

# Start from your object
df <- trio_rel3 %>%
  mutate(rel = rel / 2) %>%
  filter(!is.na(billabong_both)) %>%
  # 1) Clean strings & fix typos
  mutate(
    billabong_both = str_squish(billabong_both),
    billabong_both = str_replace_all(billabong_both, fixed("Tymualagun"), "Tyumalagun"), 
    year_caught_both = str_squish(year_caught_both),
    year_caught_both = sapply(strsplit(year_caught_both, ":"), function(v) paste(sort(v), collapse = ":"))) %>%
  # 2) Canonicalise pair order: "A:B" sorted alphabetically
  mutate(
    billabong_both = sapply(strsplit(billabong_both, ":", fixed = TRUE), function(v) {
      v <- trimws(v)
      if (length(v) < 2 || anyNA(v)) return(NA_character_)
      paste(sort(v), collapse = ":")
    })
  ) %>%
  tidyr::drop_na(billabong_both)

##

# Helper: "B:A" or "2014_2012" -> "A:B"
canon_pair <- function(x) {
  x <- as.character(x)
  x <- gsub("_", ":", x, fixed = TRUE)
  vapply(strsplit(x, ":", fixed = TRUE), function(v) {
    v <- trimws(v)
    if (length(v) < 2 || anyNA(v)) return(NA_character_)
    paste(sort(v), collapse = ":")
  }, FUN.VALUE = character(1))
}

# 1) Clean + canonicalise both columns
df <- trio_rel3 %>%
  mutate(
    billabong_both   = str_squish(billabong_both),
    billabong_both   = gsub("Tymualagun", "Tyumalagun", billabong_both, fixed = TRUE),
    billabong_both   = canon_pair(billabong_both),
    
    year_caught_both = str_squish(year_caught_both),
    year_caught_both = canon_pair(year_caught_both)
  ) %>%
  tidyr::drop_na(billabong_both, year_caught_both)

# If you want one row per unique combo, uncomment:
# %>% distinct(billabong_both, year_caught_both, .keep_all = TRUE)

# 2) Build ordering for billabong pairs (same location first, then freq, then median rel)
ord <- df %>%
  group_by(billabong_both) %>%
  reframe(
    n_pairs  = n(),
    med_rel  = median(rel, na.rm = TRUE),
    same_loc = grepl("^(.*):\\1$", billabong_both),
    .groups = "drop"
  ) %>%
  arrange(desc(same_loc), desc(n_pairs), med_rel) %>%
  pull(billabong_both) %>%
  unique()

# 3) Apply factor levels
df$billabong_both <- factor(df$billabong_both, levels = ord)

# (Optional) inspect
levels(df$billabong_both)
df$year_caught_both <- factor(df$year_caught_both)

levels(df$year_caught_both)

write.csv(df, "C:/Users/samue/Desktop/trioml_cormlpe.csv")

df <- read.csv("C:/Users/samue/Desktop/Honours/trioml_cormlpe.csv")

# 4) Plot (labels on 2 lines)
rel_bp <- ggplot(df, aes(x = billabong_both, y = rel)) +
  geom_hline(yintercept = c(0.092, 0.158, 0.204, 0.296),
             colour = "black", linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = 0.125, linewidth = 1.1, colour = "skyblue") +
  geom_hline(yintercept = 0.250, linewidth = 1.1, colour = "navyblue") +
  geom_jitter(aes(fill = as.factor(Within_Cohort)),
              width = 0.15, height = 0,
              size = 3.2, shape = 21, colour = "black", alpha = 0.85) +
  scale_fill_manual(values = c("0" = "darkgrey", "1" = "black")) +
  scale_x_discrete(labels = function(x) gsub(":", "\n", x)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
    legend.position = "none",
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 14)
  ) +
  xlab("Pairwise Location") +
  ylab("Relatedness")

rel_bp

####

head(df)



df$billabong_both <- recode_factor(df$billabong_both, 
                                   "Daly River:Tymualagun" = "Tyumalagun:Daly River",
                                   "Daly River:Dalagarr" = "Dalagarr:Daly River", 
                                   "Wunullen:Dalagarr" = "Dalagarr:Wunullen", 
                                   "Tyumalagun:Dalagarr" = "Dalagarr:Tyumalagun")

levels(df$billabong_both)

##### Prelim visualisation #####

rel_bp <- ggplot(data = df, aes(x = billabong_both, y = rel)) +
  
  geom_hline(yintercept = 0.092, col = "black", linewidth = 1, linetype="dashed") +
  geom_hline(yintercept = 0.158, col = "black", linewidth = 1,  linetype="dashed") +
  geom_hline(yintercept = 0.204, col = "black", linewidth = 1, linetype="dashed") +
  geom_hline(yintercept = 0.296, col = "black", linewidth = 1,  linetype="dashed") +
  geom_hline(yintercept = 0.125, linewidth = 1.1, col = "skyblue") +
  geom_hline(yintercept = 0.250, linewidth = 1.1, col = "navyblue") + 
  
  geom_jitter(aes(fill = as.factor(Within_Cohort)),
              size = 3.5, 
              col = "black", 
              shape = 21, 
              alpha = 0.9) + 
  scale_fill_manual(values = c("0" = "darkgrey", "1" = "black")) +
  
  theme_bw()+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none", 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  
  xlab("Pairwise Location") +
  ylab("Relatedness")

plot(rel_bp)

devEMF::emf("C:/Users/samue/Desktop/emibd9_bp.emf", width = 10, height = 6)
print(rel_bp)
dev.off()

## Combined plot

rel_combined <- ggarrange(rel.hist, 
                      rel_bp,
                      ncol = 1, 
                      nrow = 2, 
                      align = "v", 
                      labels = c("A", "B"), 
                      legend = "none")

print(rel_combined)

devEMF::emf("C:/Users/samue/Desktop/trioML_plot.emf", width = 11, height = 12)  
print(rel_combined)
dev.off()























