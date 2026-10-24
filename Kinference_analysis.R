
## Load packages 

library(adegenet)
library(gdsfmt)
library(SNPRelate)
library(gbasics)
library(kinference)
library(dartRverse)
library(adegenet)
library(ggplot2)
library(hierfstat)
library(dplyr)
library(devtools)
library(fpeek)
library(mvbutils)
library(atease)
library(vecless)
library(gridExtra)

## Load relevant datasets + make relevant datasets 

citation("dartRverse")

gl <- get(load("C:/Users/samue/Desktop/Honours/kinference_analysis/daly_geno_raw.Rdata")); gl

#View(as.matrix(gl))

gl <- gl.filter.monomorphs(gl, v = 5); gl

gl <- gl.filter.callrate(gl, method = "ind", t = 0.7); gl

gl <- gl.filter.monomorphs(gl, v = 5); gl


#gl.smearplot(gl)

## Our GL object has:
# 114 inds
# 4843 loci

# 1 pop

setwd("C:/Users/samue/Desktop/Honours/kinference_analysis")

gl2vcf(gl, plink.bin.path = "C:/Users/samue/Desktop/Honours/kinference_analysis/plink",
       outfile = "PMD_geno.vcf",
       outpath = "C:/Users/samue/Desktop/Honours/kinference_analysis")

df <- read_vcf2snpgeno(vfilename = "C:/Users/samue/Desktop/Honours/kinference_analysis/PMD_geno.vcf")

head(df)

## OR 

gl2gds(gl, 
        outfile = "PMD_gds.gds", 
        outpath = "C:/Users/samue/Desktop/Honours/kinference_analysis",
        v = 5)

snpgdsOpen(filename = "C:/Users/samue/Desktop/Honours/kinference_analysis/PMD_gds.gds", allow.duplicate = TRUE)

gds <- "C:/Users/samue/Desktop/Honours/kinference_analysis/PMD_gds.gds"

g <- SNPRelate::snpgdsOpen(gds, allow.duplicate = TRUE)
on.exit(SNPRelate::snpgdsClose(g), add = TRUE)

snpg <- gbasics::read_snpgds2snpgeno(filename = gds, 
                    locusID = "snp.id",
                    sampleID = "sample.id",
                    infoFrame = NULL,
                    locinfoFields = c(
                                      "snp.position", 
                                      "snp.chromosome",
                                      "snp.allele"))

## assuming snpg@diplos == genotypes4_ambig

## (genotypes4_ambig is defined by define_genotypes() )

## Skip me

df@diplos

snpg1 <- df

#snpg1[ df == "AB"] <- "BBO" ## recode 'AB' to 'OO'
#snpg1[ df == "BBO"] <- "AB" ## recode 'OO' to 'BBO'
#snpg1[ df == "AAO"] <- "AB" ## ...and so on
#snpg1[ df == "BBO"] <- "AAO"

## ...and then proceed with analysis using 'snpg1'
  
## Load to be compatible with kinference

head(snpg1)

snpg1 <- est_ALF_ABO_quick(snpg1)


saveRDS(df, "C:/Users/samue/Desktop/Honours/kinference_analysis/PMD_geno.snpgeno")

## Begin QA/QC

df <- snpg1 

pvals <- kinference::check6and4(df, thresh_pchisq_6and4 = c(1e-4, 1e-6), show6 = FALSE)

## Remove the worst loci p < 1x10^-10

df_1 <- df[, pvals$pval4 > 1e-8] 
str(df_1) ## 3336 loci retained from ~5000

pvals_1 <- check6and4(df_1, thresh_pchisq_6and4 = c(1e-4, 1e-6), show6 = FALSE) # this seems to look best? 

## Check for duplicated samples 

dups <- find_duplicates(df_1, max_diff_loci = 50) # none that are alarming yay!

## Check for fish from different pops/species

ilglk <- ilglk_geno(df_1, list(nclass = 100, xlim = c(-2, 2))) ## do I even need a cutoff?!

## Estimate power to find kin i.e., how well will this set of markers work 

df_2 <- kinference::kin_power(df_1, k = 0.5) 

#head(df_2$locinfo)

## Examine the heterozygosities of our inds: 
## Importnant for nulls + contamination (low + high)

# Contam 

hetz_rich <- hetzminoo_fancy(df_2, "rich", showPlot = FALSE, hist_pars = NULL)

# Degrdaded samples 

hetz_poor <- kinference::hetzminoo_fancy(df_2, "poor", showPlot = FALSE)

## Quick clean up 

hetzminoo <- data.frame(rich = hetz_rich, 
                        poor = hetz_poor)

hetzminoo # has both rich and poor

rich_plot <- ggplot(data = hetzminoo) + 
  geom_histogram(aes(x = rich), col = "black", fill = "grey") +
  geom_vline(xintercept = 0.2, col = "black", linewidth = 1, linetype="dashed") +
  xlab("Rich")+
  ylab("Frequency")+
  theme_bw()

print(rich_plot)

poor_plot <- ggplot(data = hetzminoo) + 
  geom_histogram(aes(x = poor), col = "black", fill = "pink") +
  geom_vline(xintercept = 0.2, col = "black", linewidth = 1, linetype="dashed") +
  xlab("Poor") +
  ylab("Frequency")+
  theme_bw()

hetz_plot <- ggarrange(rich_plot, 
                       poor_plot, 
                       nrow = 2, 
                       ncol = 1, 
                       labels = c("A", "B"), 
                       legend = "none")

print(hetz_plot)

devEMF::emf("C:/Users/samue/Desktop/kinference_supps/hetz_polot.emf", width = 12, height = 8)  
print(hetz_plot)
dev.off()


df_3 <- df_2[(hetz_rich > 0.2)]; str(df_3)

hetz_poor <- hetzminoo_fancy(df_3, "poor")

# sanity check # sanity checkdf_2 

## Revisit our allele frequencies

pvals_2 <- check6and4(df_3, thresh_pchisq_6and4 = c(1e-4, 1e-6), show6 = FALSE)

## It appears we have only removed our good loci? 

lociPvals <- data.frame(pvals = c(pvals_2$pval4, pvals_1$pval4),
                        stage = c(rep("after", length(pvals_2$pval4)),
                                  rep("before", length(pvals_1$pval4)))
)

lm1 <- lm(pvals ~ stage, data = lociPvals)
summary(lm1) 

## The p-values improved by aout 3.3%

###-------------------------------
### Kin Finding 
###-------------------------------

## Kin finding workflow ##

# Whilst unlikely, check data for POPS first

pops <- find_POPs(df_3, limit_pairs = choose(nrow(df_3), 2), keep_thresh = 1, nbins = 100)

hist(pops$wpsex, breaks = 30, xlim = c(0.018, 0.15), ylim = c(0, 20), xlab = "wpsex")

## There's a bit of a dip around 0.04, but it's not _clear_ separation
abline(v = 0.04, lty = 2, col = "red")
box()

with(pops, plot.default(nABOO ~ wpsex))
abline(v = 0.04, col = "red", lty = 2)

splitPOPs <- split_FSPs_from_POPs(df_3, pops[pops$wpsex < 0.055, ], gerr = 0.01)

str(splitPOPs)

ggplot() + 
  geom_histogram(data = splitPOPs, aes(x = PLOD_FP), bins = 50, fill = "grey", col = "black") + 
  geom_vline(xintercept = splitPOPs@E_FSP, col = "orange") + 
  geom_vline(xintercept = splitPOPs@E_POP, col = "navyblue") + 
  theme_bw()

splitPOPs@E_POP

POPs <- splitPOPs[splitPOPs$PLOD_FP < -50,] ## our first set of definite kin pairs
POPs$ij <- paste(POPs$i, POPs$j) ## a pair ID

POPs

## find HSPs-or-stronger using find_HSPs

HSPs_1 <- find_HSPs(df_3, keep_thresh = -20, limit_pairs = choose(nrow(df_3),2))

PLOD_loghisto(HSPs_1)

PLOD_plot <- plot_histoPLOD_gg(HSPs_1, UP = TRUE, HSP = TRUE, POP = TRUE, FSP = TRUE,
                  showUP = c(SPA = TRUE, Normal = FALSE))


print(PLOD_plot)

devEMF::emf("C:/Users/samue/Desktop/kinference_PLOD_plot.emf", width = 10, height = 6)  
print(PLOD_plot)
dev.off()

HSP_histo(HSPs_1, fullsib_cut = 300, lb = -10, bin = 10)

HSP_plod_plot <- hsp_histo_gg(HSPs_1,   # data.frame with column 'PLOD'
             lb            = -10,
             ub            = 500,
             fullsib_cut   = 300,         # threshold used for HSP variance window
             bin           = 10,
             HSPmean = TRUE, HSPdist = TRUE,
             POPmean = TRUE, FSPmean = TRUE, UPmean = FALSE)

print(HSP_plod_plot)

devEMF::emf("C:/Users/samue/Desktop/HSP_PLOD_plot.emf", width = 10, height = 6)  
print(HSP_plod_plot)
dev.off()

## The hist looks weird asf with high freq of HSP? 

POPsOrFSPs <- HSPs_1[HSPs_1$PLOD > 300,]

POPsOrFSPs$ij <- paste(POPsOrFSPs$i, POPsOrFSPs$j) ## a pair ID to match the one in POPs

POPsOrFSPs

FSPs <- POPsOrFSPs[!POPsOrFSPs$ij %in% POPs$ij,] ## find our full-sibling pairs

FSPs

nrow(FSPs)

## This should give us all the FSPs in our data? 

HSPs_1$ij <- paste(HSPs_1$i, HSPs_1$j)

HSPs_2 <- HSPs_1[!HSPs_1$ij %in% c(POPs$ij, FSPs$ij),]

## plot the empirical PLOD scores with the expected distribution
## of the PLOD for half-sibling pairs

thresh <- autopick_threshold(df_3, kin = HSPs_2,
                             fitrange_PLOD = c(0, 250),
                             FPtol_pairs = 1, 
                             use4th = TRUE, 
                             selecto = "ML", 
                             plot_bins = 5)

thresh

    # <- the ggplot object
thresh_gg <- autopick_threshold_gg(df_3, kin = HSPs_2,
                             fitrange_PLOD = c(0, 250),
                             FPtol_pairs = 1, 
                             use4th = TRUE, 
                             selecto = "ML", 
                             plot_bins = 5)


p <- attr(thresh_gg, "plot")

devEMF::emf("C:/Users/samue/Desktop/Autothreshold_PLOD_plot.emf", width = 10, height = 6)  
print(p)
dev.off()

HSPs <- HSPs_2[HSPs_2$PLOD > thresh,]

fpos_rate <- thresh@FPtol_pairs / nrow(HSPs); fpos_rate ## approx 2% of our 'HSPs' are false-positives
fneg_rate <- thresh@info["Pr_FNeg"]; fneg_rate ## we're missing about 0.5% of HSPs to false-negatives
nrow(HSPs) ## 52 inferred half-sibling pairs

## Add intividual metadata to the kin pairs ##

POPs # This is a false positive
nrow(FSPs)
nrow(HSPs)

kinf_kin <- rbind(FSPs, HSPs); kinf_kin
dim(kinf_kin)

kinf_kin$rel <- ifelse(kinf_kin$PLOD > 300, "full-sibs", "half-sibs") # this needs to be adjusted to fit the datasets...

kinf_kin

ind_names <- df_3$info; ind_names

ind_names <- ind_names %>%
  mutate(Our_sample = sub("^pop1_", "", Our_sample)) %>%
  mutate(row_id = row_number())

ind_names

kinf_kin <- kinf_kin %>%
  left_join(ind_names, by = c("i" = "row_id" )) %>%
  rename(ind_i = Our_sample)%>%
  left_join(ind_names, by = c("j" = "row_id" )) %>%
  rename(ind_j = Our_sample)

kinf_kin

ind_metadata <- read.csv("C:/Users/samue/Desktop/10658_prelim_report/Pristis_metadata_SRA.csv")

head(ind_metadata)

kinf_kin2 <- kinf_kin %>%
  left_join(ind_metadata, by = c("ind_i" = "id")) %>%
  rename(birth_year_1 = birth_year, billabong_1 = Site_ID, capture_year1 = year_caught) %>%
  left_join(ind_metadata, by = c("ind_j" = "id")) %>%
  rename(birth_year_2 = birth_year, billabong_2 = Site_ID, capture_year2 = year_caught) %>%
  mutate(
    birth_year_diff = abs(birth_year_1 - birth_year_2),  # Absolute age difference
    pw_loc = as.integer(billabong_1 == billabong_2), # returns true false as 1,0
    pw_yr = as.integer(capture_year1 == capture_year2)) %>% #returns true false as 1,0
  rename(id_1 = ind_i, id_2 = ind_j) %>%
  dplyr::select(id_1, id_2, rel,billabong_1, billabong_2, birth_year_1, birth_year_2, capture_year1, capture_year2, birth_year_diff, pw_loc, pw_yr)
  

kin_pairs <- tibble(kinf_kin2)

print(kin_pairs)

table(kin_pairs$rel, kin_pairs$pw_yr)

## Need to add all other kin pairs to the dataset and then export for visualisation

## Use the ckmrsim kin pair dataset (with all combos) to get the non-kin

all_pairs <- read.csv("C:/Users/samue/Desktop/PMD_ckmrsim_complete.csv")

dim(all_pairs)

all_pairs <- all_pairs %>%
  dplyr::select(id_1, id_2, rel, pw_loc, pw_yr)

## now we need 

dim(kinf_kin2)

kinf_kin <- kinf_kin2 %>%
  dplyr::select(id_1, id_2, rel, pw_loc, pw_yr)


names(kinf_kin)
names(all_pairs)

write.csv(kin_pairs, "C:/Users/samue/Desktop/kinference_kin.csv")

head(df_3)
glimpse(df_3)

more_kin <- pops %>%
  left_join(ind_names, by = c("i" = "row_id" )) %>%
  rename(ind_i = Our_sample)%>%
  left_join(ind_names, by = c("j" = "row_id" )) %>%
  rename(ind_j = Our_sample)

more_kin

## 1. Standardise pair order in both data frames

pairs_df <- more_kin %>%
  rename(id_1 = ind_i, id_2 = ind_j) %>%
  mutate(
    id_1 = str_trim(as.character(id_1)),
    id_2 = str_trim(as.character(id_2)),
    a = pmin(id_1, id_2),
    b = pmax(id_1, id_2)
  )

kinf_kin <- kinf_kin %>%
  mutate(
    id_1 = str_trim(as.character(id_1)),
    id_2 = str_trim(as.character(id_2)),
    a = pmin(id_1, id_2),
    b = pmax(id_1, id_2)
  ) %>%
  distinct(a, b, .keep_all = TRUE)   # <- prevent many-to-many joins

head(pairs_df)
head(kinf_kin)

## 2. Join and fill in missing with "unrelated"

out <- pairs_df %>%
  left_join(
    kinf_kin %>% select(a, b, rel),
    by = c("a","b")) %>%
  mutate(rel = ifelse(is.na(rel), "unrelated", rel)) %>%
  dplyr::select(id_1 = a, id_2 = b, rel)


out <- pairs_df %>%
  dplyr::left_join(kinf_kin[, c("a","b","rel")], by = c("a","b")) %>%
  dplyr::mutate(rel = ifelse(is.na(rel), "unrelated", rel)) %>%
  dplyr::select(id_1 = a, id_2 = b, rel)

head(out)
table(out$rel)


write.csv(out, "C:/Users/samue/Desktop/kinference_all.csv")
