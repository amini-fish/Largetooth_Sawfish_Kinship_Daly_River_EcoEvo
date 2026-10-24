#  mtDNA pipeline: read DArT sheets → FASTA → concatenate → hap --------


# Load packages -----------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(janitor)
  library(ggplot2)
  library(Biostrings)
  library(DECIPHER)
  library(apex)
  library(ape)
  library(pegas)
  library(RColorBrewer)
})

# Paths -------------------------------------------------------------------
base_dir      <- "C:/Users/samue/Desktop/Research/Honours/Chapter_2_Relatedness_paper/mtDNA_data"
dart_report   <- "C:/Users/samue/Desktop/Research/Honours/Chapter_2_Relatedness_paper/mtDNA_data/Report-DPri25-10659.xlsx"
out_dir       <- base_dir
meta_csv      <- "C:/Users/samue/Desktop/Research/Honours/10658_prelim_report/Pristis_metadata_SRA.csv"
kin_pairs_csv <- "C:/Users/samue/Desktop/Research/Honours/kinference_all.csv"


## Sheets → FASTA outputs
targets <- list(
  CTYB    = file.path(out_dir, "cytb.fasta")
  #ND4     = file.path(out_dir, "nd4.fasta"),
  #CONTROL = file.path(out_dir, "control.fasta")
)

targets


# Read in mtDNA data ------------------------------------------------------

read_mt_sheet_to_fasta <- function(dart_report, sheet, fasta_out, qc_hist = FALSE) {
  stopifnot(file.exists(dart_report))
  
  df <- read_excel(dart_report, sheet = sheet) %>% 
    clean_names()
  
  if (!all(c("sequence", "genotype_name") %in% names(df))) {
    stop(sprintf("Sheet '%s' missing required columns.", sheet))
  }
  
  df$sequence <- sub("-+$", "", df$sequence)
  
  if (qc_hist && "total_counts_ont_all_amplicons" %in% names(df)) {
    ggplot(df, aes(total_counts_ont_all_amplicons)) +
      geom_histogram(bins = 50) +
      labs(title = paste("QC:", sheet))
  }
  
  fa <- DNAStringSet(stats::setNames(df$sequence, df$genotype_name))
  writeXStringSet(fa, fasta_out)
  invisible(fasta_out)
}

## Process all sheets
invisible(lapply(names(targets), function(s) {
  read_mt_sheet_to_fasta(dart_report, 
                         s, 
                         targets[[s]], 
                         qc_hist = FALSE)
}))


# Read in FASTA files and concatenate -------------------------------------

fasta_files <- unname(unlist(targets))
stopifnot(all(file.exists(fasta_files)))

all_fasta <- apex::read.multiFASTA(fasta_files, add.gaps = TRUE)
concat    <- apex::concatenate(all_fasta, exclude.gaps = TRUE)
stopifnot(inherits(concat, "DNAbin"))


# Build haplotypes  -------------------------------------------------------
haps <- pegas::haplotype(concat)
net  <- pegas::haploNet(haps)


# Input metadata to inform grouping ---------------------------------------
meta <- read.csv(meta_csv, stringsAsFactors = FALSE) %>%
  mutate(id = as.character(id))

rownames(concat) <- as.character(rownames(concat))

keep_ids <- intersect(meta$id, rownames(concat))
concat   <- concat[keep_ids, , drop = FALSE]
meta     <- meta[match(keep_ids, meta$id), ]

grp <- factor(meta$year_caught)
hf  <- pegas::haploFreq(concat, grp)

pal <- brewer.pal(max(3, nlevels(grp)), "YlOrRd")

plot(net,
     pie = hf,
     bg = pal,
     size = sqrt(attr(net, "freq")),
     col = "black",
     col.link = "black",
     show.mutation = 1,
     labels = TRUE)
legend("topleft", legend = levels(grp), fill = pal, bty = "n")


# Assign individuals to haplotypes ----------------------------------------
make_haplotype_map <- function(concat_dnabin, haps_obj) {
  idx <- attr(haps_obj, "index")
  ids <- rownames(concat_dnabin)
  
  hap_names <- rownames(haps_obj)
  if (is.null(hap_names) || any(hap_names == "")) {
    hap_names <- paste0("H", seq_along(idx))
  }
  
  hap_df <- data.frame(
    haplotype = rep(hap_names, lengths(idx)),
    row_idx   = unlist(idx),
    stringsAsFactors = FALSE
  )
  
  hap_df$individual <- ids[hap_df$row_idx]
  hap_df$row_idx <- NULL
  
  id2hap <- stats::setNames(hap_df$haplotype, hap_df$individual)
  list(hap_df = hap_df, id2hap = id2hap)
}

map    <- make_haplotype_map(concat, haps)
hap_df <- map$hap_df
id2hap <- map$id2hap


# Assign kin to maternal haplotypes ---------------------------------------
kin_pairs <- read.csv(kin_pairs_csv, stringsAsFactors = FALSE) %>%
  mutate(id_1 = as.character(id_1),
         id_2 = as.character(id_2))

kin_annot <- kin_pairs %>%
  mutate(
    hap_1 = unname(id2hap[id_1]),
    hap_2 = unname(id2hap[id_2]),
    same_haplotype = (!is.na(hap_1) & !is.na(hap_2) & hap_1 == hap_2)
  )

print(table(kin_annot$rel, kin_annot$same_haplotype))

missing_ids <- setdiff(unique(c(kin_annot$id_1, kin_annot$id_2)), names(id2hap))
message("Missing IDs in mtDNA alignment: ", length(missing_ids))

## Add in cohorts to see distribution of maternal paternal kin over time

kin_annot <- kin_annot %>%
  left_join(
    meta %>% select(id, year_caught),
    by = c("id_1" = "id")
  ) %>%
  dplyr::rename(year_caught_1 = year_caught)

kin_annot <- kin_annot %>%
  left_join(
    meta %>% select(id, year_caught),
    by = c("id_2" = "id")
  ) %>%
  dplyr::rename(year_caught_2 = year_caught)

kin_annot

write.csv(kin_annot, file.path(out_dir, "kinference_mtDNA_results.csv"), row.names = FALSE)


# mtDNA diversity indices and general population genetics -----------------
nuc_div <- pegas::nuc.div(concat, variance = T); nuc_div
se_nuc <- sqrt(1.431903e-08); se_nuc

hap_div <- pegas::hap.div(concat, variance = T); hap_div
hap_se <- sqrt(0.0008175514); hap_se

haplotypes_all <- pegas::haplotype(concat)

## By cohort 

meta <- read.csv("C:/Users/samue/Desktop/Research/Honours/10658_prelim_report/Pristis_metadata_SRA.csv")
meta$id <- as.character(meta$id)

rownames(concat) <- as.character(rownames(concat))

# Match order and keep only overlapping IDs
keep <- intersect(meta$id, rownames(dna))
concat <- concat[keep, ]
meta <- meta[match(keep, meta$id), ]

# Helper function for mtdna to consolidate outputs 

mtDNA_stats <- function(dna_sub) {
  haps <- pegas::haplotype(dna_sub)   # compute haplotypes for this subset only
  data.frame(
    N          = nrow(dna_sub),
    Haplotypes = nrow(haps),
    Hd         = pegas::hap.div(dna_sub, variance = T),
    Pi         = pegas::nuc.div(dna_sub, variance = T),
    Theta_S    = pegas::theta.s(dna_sub),
    Tajimas_D  = pegas::tajima.test(dna_sub)$D
  )
}

# Loop it by year

year_stats <- meta %>%
  split(.$year_caught) %>%
  lapply(function(df) {
    ids <- df$id
    # Subset the DNAbin correctly (drop = FALSE keeps matrix structure)
    dna_sub <- concat[rownames(concat) %in% ids, , drop = FALSE]
    mtDNA_stats(dna_sub)
  }) %>%
  bind_rows(.id = "year_caught")

year_stats

## Outputs (she ain't pretty)

h_se11 <- sqrt(0.002940365); h_se12
n_se11 <- sqrt(1.312654e-06); n_se12

h_se12 <- sqrt(0.002940365); h_se12
n_se12 <- sqrt(1.312654e-06); n_se12

h_se17 <- sqrt(0.005029106); h_se17
n_se17 <- sqrt(3.008383e-07); n_se17

h_se18 <- sqrt(0.004475992); h_se18
n_se18 <- sqrt(8.297118e-07); n_se18

h_se24 <- sqrt(0.002095053); h_se24
n_se24 <- sqrt(5.440222e-07); n_se24


# Random match probability and such... ------------------------------------

## random match probability 
counts <- c(2, 64, 40, 3, 15, 1, 2, 1) #manual and gross but it works...be sure to adjust
N <- sum(counts)

# Probability two random individuals share a haplotype
P_match <- sum(counts * (counts - 1)) / (N * (N - 1))
P_match


hap_counts <- as.integer(attr(haps, "freq"))
N_total    <- sum(hap_counts)

P_match <- if (N_total > 1) {
  sum(hap_counts * (hap_counts - 1)) / (N_total * (N_total - 1))
} else NA

message(sprintf("Random match probability (shared haplotype): %.6f", P_match))


# Clean network plot ------------------------------------------------------

plot(net,
     pie = hf,
     bg = pal,
     size = sqrt(attr(net, "freq")),
     col = "black",
     col.link = "black",
     show.mutation = 2,
     threshold = c(1, 10),
     cex = 1.6,
     lwd = 1.8,
     labels = TRUE,
     directed = TRUE)

legend("topleft", legend = levels(grp), fill = pal, bty = "n", cex = 0.9)

