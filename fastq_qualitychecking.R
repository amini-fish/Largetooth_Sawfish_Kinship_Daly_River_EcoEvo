if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ShortRead")
install.packages("fastqcr")

## load packages ## 

library(ShortRead)
library(ggplot2)

## Load our fastq files as raw reads to asses the per. base pair quality

file  = "C:/Users/samue/Desktop/Honours/Chapter_2_Relatedness_paper/mtDNA_data/fastq_files/4457856.FASTQ.gz"

fq <- readFastq(file)

## Some quick eda 

reads = sread(fq)
head(reads)

summary(reads) ## approximately 10065 bp 

## Now we want to focus on the read lengths more 

widths = as.data.frame(reads@ranges@width); head(widths)

summary(widths$`reads@ranges@width`)

ggplot(data = widths) +
  geom_histogram(aes(x = reads@ranges@width), 
                 colour = "black") +
  theme_bw()

## graph the quality scores 

quals = quality(fq)

head(quals) # letters = good, symbols = bad

# now we want the average quality score

numqscores = as(quals, 'matrix')
avgscore = rowMeans(numqscores, na.rm = T)
avgscore = as.data.frame(avgscore)

ggplot(avgscore) +
  geom_histogram(aes(x = avgscore), 
                 col = 'black') +
  theme_bw()

## to read multiple files I guess you can do this: 

fastq_files <- list.files("C:/Users/samue/Desktop/Honours/Chapter_2_Relatedness_paper/mtDNA_data/fastq_files/", pattern = ".FASTQ.gz", full.names = TRUE)

fastq_files

all_reads <- lapply(fastq_files, readFastq); all_reads

## Try fastqc 

?fastqc_install()

fastqcr::fastqc(fq.dir = "C:/Users/samue/Desktop/Honours/Chapter_2_Relatedness_paper/mtDNA_data/fastq_files", # FASTQ files directory
                qc.dir = "C:/Users/samue/Desktop/Honours/Chapter_2_Relatedness_paper/mtDNA_data/fastqc_results", # Results direcory
                threads = 10)

### Straight from bioconductor website

fls <- dir("C:/Users/samue/Desktop/Honours/Chapter_2_Relatedness_paper/mtDNA_data/fastq_files", 
           ".FASTQ.gz", full=TRUE)

qaSummary <- qa(fls, type="fastq")

qaSummary

head(qaSummary[["readCounts"]])

summary(qaSummary[["readCounts"]])

hist(qaSummary[["readCounts"]]$read, breaks = 50, main = "Read counts per sample", xlab = "Reads")

head(qaSummary[["baseCalls"]])

summary(qaSummary[["readQualityScore"]])

#browseURL(report(qaSummary))

summary(qaSummary[["perCycle"]])

df_q <- qaSummary[["readQualityScore"]]

library(ggplot2)

ggplot(df_q, aes(x = quality, y = density)) +
  geom_line() +
  theme_bw() +
  labs(title = "Per-read average quality distribution",
       x = "Average Phred quality per read",
       y = "Density")


rc <- qaSummary[["readCounts"]]
ggplot(rc, aes(x = read)) + 
  geom_histogram(bins = 50) +
  theme_bw() +
  labs(x = "Reads per file", y = "Number of samples")

base_df <- as.data.frame(qaSummary[["baseCalls"]])
head(base_df)


fq <- ShortRead::readFastq(fls[1])
table(ShortRead::width(fq))

citation("pegas")
citation("hierfstat")
