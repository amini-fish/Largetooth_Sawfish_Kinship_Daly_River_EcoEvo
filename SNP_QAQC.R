setwd("C:/Users/samue/Desktop/Honours/analysis")

#Install some packages

install_github("green-striped-gecko/dartR.sexlinked@dev")
#install.packages("hierfstat")
#devtools::install_version("ggplot2", "3.4.4")

install.packages("HardyWeinberg")

## Load all required packages
library(dartR.sexlinked)
library(dartRverse)
library(ggplot2)
library(hierfstat)
library(dplyr)
library(devtools)
library(HardyWeinberg)

data <- "C:/Users/samue/Desktop/10658_prelim_report/Report_DPri25-10658_SNP.csv"
meta <- "C:/Users/samue/Desktop/10658_prelim_report/Pristis_metadata_SRA.csv"

## Load all genotype data and metadata

data.gl <- dartR.base::gl.read.dart(filename = data, 
                                    ind.metafile = meta,
                                    nas = "-", 
                                    topskip = 6, 
                                    probar = T); data.gl


save(data.gl, file = "C:/Users/samue/Desktop/daly_geno_raw.Rdata")


gl <- get(load("C:/Users/samue/Desktop/daly_geno_raw.Rdata")); gl

gl <- gl.filter.allna(gl, recalc = T, v = 5); gl

unknowns <- gl@other$ind.metrics[gl@other$ind.metrics$sex == "Unknown",1]; unknowns

gl_sexy <- gl.drop.ind(gl, ind.list = unknowns, recalc = T, v = 5)

gl_sexy@other$ind.metrics$sex <- droplevels(gl_sexy@other$ind.metrics$sex)

sex_linked <- gl.keep.sexlinked(gl, 
                    system = "xy",
                    ncores = 10, 
                    plot.theme = theme_bw())

gl.infer.sex(sex_linked, system = "xy", seed = 42)


gl.smearplot(data.gl) 

data.gl <- gl.keep.pop(data.gl, pop.list = "Daly")

data.gl <- gl.filter.monomorphs(data.gl, v = 5); data.gl

## Save unfiltered data 

#save(data.gl, file = "pristis_geno_raw.RData") 

###### FILTERING #####

## Get rid of sex linked loci 

## Get rid of unreliable loci

gl.report.reproducibility(data.gl)

data.gl <- dartR.base::gl.filter.reproducibility(data.gl, threshold = 0.98)

## Callrate 

dartR.base::gl.report.callrate(data.gl)

data.gl <- dartR.base::gl.filter.callrate(data.gl, method = "loc", threshold = 0.85)

data.gl

#Get rid of low and super high read depth loci

dartR.base::gl.report.rdepth(data.gl)

data.gl <- dartR.base::gl.filter.rdepth(data.gl, lower = 10, upper = 150)

data.gl

## Remove secondary loci - important to keep unlinked loci

data.gl <- dartR.base::gl.filter.secondaries(data.gl, method = "best") 

data.gl@other$loc.metrics

#Very low filter – this is only to get rid of your really bad individuals

dartR.base::gl.report.callrate(data.gl, method = "ind")

data.gl <- dartR.base::gl.filter.callrate(data.gl, method = "ind", threshold = 0.8)

data.gl

data.gl <- dartR.base::gl.filter.monomorphs(data.gl)

data.gl

#Always run this after removing individuals – removes loci that are no longer variable

data.gl <- dartR.base::gl.filter.hwe(data.gl)

data.gl

gl.report.maf(data.gl)

data.gl <- gl.filter.maf(data.gl, t = 0.025, v = 5)


data.gl


## remove evidence of DNA contamination ## Important for kin finding 

dartR.base::gl.report.heterozygosity(data.gl, method = "ind")

data.gl <- dartR.base::gl.filter.heterozygosity(data.gl,t.lower = 0.225,  t.upper = 0.5) ## Visualise our cleaned data 

gl.smearplot(data.gl)

data.gl

####################################################

data.gl <- gl.filter.allna(data.gl)

data.gl

data.gl <- gl.filter.monomorphs(data.gl); data.gl

gl.report.ld(data.gl, probar = T, nchunks = 8)

data.gl
### SAVE DALY RIVER INDS ###


save(data.gl, file = "C:/Users/samue/Desktop/PMD_clean.Rdata")

## LD??

### Check our filtering steps ###

data.gl@other$history

### SAVE CLEANED SNP DATA ###

save(data.gl, file = "pristis_geno_cleaned.Rdata")

?gl.report.heterozygosity

gl <- get(load("C:/Users/samue/Desktop/PMD_clean.Rdata")); gl

gl.report.diversity(gl)

pophet <- gl.report.heterozygosity(gl, 
                               error.bar = "SE")

pophet

pop(gl) <- gl@other$ind.metrics$year_caught

het <- gl.report.heterozygosity(gl, 
                         error.bar = "SE")


head(het)

het <- rbind(het, pophet)
rownames(het)<- c("2011", "2012", "2017", "2018", "2024", "All")

library(forcats)

het <- het %>%
  mutate(pop = fct_recode(pop, 
                          "ALL" = "pop1"))
  
head(het)

ho <- ggplot(data = het, aes(x = pop, y = Ho)) +
  geom_bar(aes(fill = pop), colour = "black" ,
           stat = "identity") +
  geom_errorbar(aes(ymin = Ho - HoSE, ymax = Ho + HoSE),
                width = 0.3) + 
  theme_bw() + 
  scale_fill_brewer(palette = "YlOrRd") + 
  xlab("Cohort") + 
  ylab("Ho") + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14), 
        legend.position = "none") + 
  labs(fill = "Ho")

print(ho)

uhe <- ggplot(data = het, aes(x = pop, y = uHe)) +
  geom_bar(aes(fill = pop),colour = "black" ,
           stat = "identity") +
  geom_errorbar(aes(ymin = uHe - uHeSE, ymax = uHe + uHeSE),
                width = 0.3) + 
  theme_bw() + 
  scale_fill_brewer(palette = "YlOrRd") + 
  xlab("Cohort") + 
  ylab("uHE") + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14), 
        legend.position = "none") + 
  labs(fill = "uHe")

print(uhe)

fis <- ggplot(data = het, aes(x = pop, y = FIS)) +
  geom_bar(aes(fill = pop), colour = "black",
             stat = "identity") +
  geom_errorbar(aes(ymin = FIS - FISSE, ymax = FIS + FISSE),
                width = 0.3) + 
  theme_bw() + 
  scale_fill_brewer(palette = "YlOrRd") + 
  xlab("Cohort") + 
  ylab("FIS") + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14),
        legend.position = "none") + 
  labs(fill = "FIS")

print(fis)

library(patchwork)
library(ggpubr)

het_plot <- ggarrange(ho , uhe, fis, ncol = 1, nrow = 3, align = "v")

print(het_plot)

devEMF::emf("C:/Users/samue/Desktop/het_plot.emf", width = 10, height = 12)  
print(het_plot)
dev.off()

## Allelic richness using Heirfstat

gl <- get(load("C:/Users/samue/Desktop/PMD_clean.Rdata")); gl
pop(gl) <- gl@other$ind.metrics$year_caught
gi <- gl2gi(gl)

hfstat <- genind2hierfstat(gi)

#calculate allelic richness

ar <- allelic.richness(hfstat, diploid = TRUE)

summary(ar$Ar)

sd(ar$Ar$`1`)/sqrt(3646)

head(ar$Ar)

ar <- as.data.frame(ar$Ar)

ar

mean.ar <- colMeans(ar); mean.ar
mean.ar

sd.ar <- lapply(ar, stats::sd)
sd.ar

head(ar)

n <- 3464

se.ar <- mapply(function(s, n) s / sqrt(n), sd.ar, n)
se.ar


library(dplyr)
library(tidyverse)

ar <- ar %>%
  rownames_to_column("ID") %>%
  pivot_longer(
    cols = starts_with("20"), 
    names_to = "Year", 
    values_to = "Ar"
  )


library(wesanderson)

colo <- c("skyblue", "#74c476",  "#FDD835", "orange", "lightpink")

AR_plot <- ggplot(data = ar, aes(x = Year, y = Ar, fill = as.factor(Year))) +
  geom_boxplot(colour = "black") +
  scale_fill_brewer(palette = "YlOrRd") +
  xlab("Cohort") +
  ylab("Allelic Richness") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14), 
        legend.position = "none") +

  labs(fill = "Year")

print(AR_plot)

het_plot <- ggarrange(ho, 
                      uhe, 
                      fis, 
                      AR_plot, 
                      ncol = 2, 
                      nrow = 2, 
                      align = "v", 
                      labels = c("A", "B", "C", "D"))

print(het_plot)

devEMF::emf("C:/Users/samue/Desktop/diversity_plot.emf", width = 10, height = 6)  
print(het_plot)
dev.off()

fis_results <- hierfstat::basic.stats(gi)

# Extract individual Fis (per locus averaged across loci)
ind_Fis <- rowMeans(fis_results$Fis, na.rm = TRUE)

ind_Fis

install.packages("inbreedR")

library(inbreedR)



####

ggarrange(p1, p2, p3, p4,
          
          ncol = 2, nrow = 2)


## ============================================================
## Pristis SNP pipeline (dartR): read → QC/filter → save → diversity
## Consolidated, annotated, and preserving EXACT save filenames
## ============================================================

## ---- (Optional) one-time installs ----
# install.packages(c("HardyWeinberg","ggpubr","patchwork"))
# devtools::install_github("green-striped-gecko/dartR.sexlinked@dev")
# install.packages("hierfstat")

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

## ---- Original input paths (unchanged) ----
data <- "C:/Users/samue/Desktop/10658_prelim_report/Report_DPri25-10658_SNP.csv"
meta <- "C:/Users/samue/Desktop/10658_prelim_report/Pristis_metadata_SRA.csv"

## -------------------- 1) Load genotypes + metadata --------------------
data.gl <- dartR.base::gl.read.dart(
  filename     = data,
  ind.metafile = meta,
  nas          = "-",
  topskip      = 6,
  probar       = TRUE
); data.gl

## ---- SAVE NAME PRESERVED ----
save(data.gl, file = "C:/Users/samue/Desktop/daly_geno_raw.Rdata")

## For convenience, keep a working copy called `gl`
gl <- get(load("C:/Users/samue/Desktop/daly_geno_raw.Rdata")); gl

## -------------------- 2) Initial tidy --------------------
gl <- gl.filter.allna(gl, recalc = TRUE, v = 5); gl

## Drop unknown-sex individuals only if present (you had this)
if ("sex" %in% names(gl@other$ind.metrics)) {
  unknowns <- gl@other$ind.metrics[gl@other$ind.metrics$sex == "Unknown", 1]
  if (length(unknowns)) {
    gl <- gl.drop.ind(gl, ind.list = unknowns, recalc = TRUE, v = 5)
    gl@other$ind.metrics$sex <- droplevels(gl@other$ind.metrics$sex)
  }
}

## -------------------- 3) Sex-linked loci: IDENTIFY then REMOVE --------------------
## (You kept sex-linked earlier by mistake; here we drop them from gl)
sex_linked <- tryCatch(
  gl.keep.sexlinked(gl, system = "xy", ncores = max(1, parallel::detectCores()-2), plot.theme = theme_bw()),
  error = function(e) NULL
)

if (!is.null(sex_linked)) {
  gl <- gl.drop.loc(gl, loc.list = locNames(sex_linked), recalc = TRUE, v = 3)
}

## -------------------- 5) Locus-level QC --------------------

dartR.base::gl.report.callrate(gl)  ## locus
gl <- dartR.base::gl.filter.callrate(gl, method = "loc", threshold = 0.85); gl

dartR.base::gl.report.rdepth(gl)
gl <- dartR.base::gl.filter.rdepth(gl, lower = 10, upper = 150);gl

## -------------------- 6) Individual-level QC --------------------
dartR.base::gl.report.callrate(gl, method = "ind")
gl <- dartR.base::gl.filter.callrate(gl, method = "ind", threshold = 0.8); gl

## Remove monomorphs after dropping poor individuals
gl <- dartR.base::gl.filter.monomorphs(gl); gl

## HWE
gl <- dartR.base::gl.filter.hwe(gl,  alpha = 0.05); gl

## MAF
gl.report.maf(gl)
gl <- gl.filter.maf(gl, t = 0.025, v = 5); gl

## Heterozygosity (contamination guardrails you used)
dartR.base::gl.report.heterozygosity(gl, method = "ind")
gl <- dartR.base::gl.filter.heterozygosity(gl, t.lower = 0.225, t.upper = 0.5)

## Visualise cleaned data (as you had)
gl.smearplot(gl)

## Final tidy
gl <- gl.filter.allna(gl)
gl <- gl.filter.monomorphs(gl); gl

## ---- SAVE NAMES PRESERVED ----
## “SAVE DALY RIVER INDS”
save(gl, file = "C:/Users/samue/Desktop/PMD_clean.Rdata")

## History of filters
gl@other$history

## “SAVE CLEANED SNP DATA”
save(gl, file = "pristis_geno_cleaned.Rdata")

