### Can we estimate the effective population size? - in the face of pop structure 

# Pierre will likely have something to say about this...

library(dartRverse)
library(ggplot2)
library(hierfstat)
library(dplyr)
library(devtools)
library(dartR.base)
library(dartR.captive)
library(gplots)
library(graph4lg)
library(viridis)
library(ggraph)
library(hierfstat)

## Load genetic data

gl <- get(load("C:/Users/samue/Desktop/PMD_clean.Rdata")); gl

table(gl@other$ind.metrics$year_caught)

## we use 2012, 2017 and 2024

Ne <- dartR.popgen::gl.LDNe(gl, 
                             mating = "random",
                             singleton.rm = TRUE,
                             Waples.correction = "nChromosomes",
                             Waples.correction.value = 46,
                             outfile = "PMD_LDne.txt",
                             neest.path = "C:/Users/samue/Desktop",
                             plot_theme = theme_bw())


Ne

## 2012 

gl2 <- gl

pop(gl2) <- gl2@other$ind.metrics$year_caught

gl_2012 <- gl.keep.pop(gl2, pop.list = "2012"); gl_2012

gl_2012

gl_2012 <- gl.filter.monomorphs(gl_2012); gl_2012

nes_2012 <- dartR.popgen::gl.LDNe(gl_2012, 
                                  mating = "random",
                                  singleton.rm = TRUE,
                                  Waples.correction = "nChromosomes",
                                  Waples.correction.value = 46,
                                  outfile = "PMD_LDne_17.txt",
                                  neest.path = "C:/Users/samue/Desktop",
                                  plot_theme = theme_bw())

nes_2012

## 2017 only 

pop(gl) <- gl@other$ind.metrics$year_caught

gl_2017 <- gl.keep.pop(gl, pop.list = "2017"); gl_2017 

gl_2017 <- gl.filter.monomorphs(gl_2017); gl_2017

nes_2017 <- dartR.popgen::gl.LDNe(gl_2017, 
                             mating = "random",
                             singleton.rm = TRUE,
                             Waples.correction = "nChromosomes",
                             Waples.correction.value = 46,
                             outfile = "PMD_LDne_17.txt",
                             neest.path = "C:/Users/samue/Desktop",
                             plot_theme = theme_bw())
nes_2017

## 2018 

gl_2018 <- gl.keep.pop(gl2, pop.list = "2018"); gl_2018

gl_2018 <- gl.filter.monomorphs(gl_2018); gl_2018

nes_2018 <- dartR.popgen::gl.LDNe(gl_2018, 
                                  mating = "random",
                                  singleton.rm = TRUE,
                                  Waples.correction = "nChromosomes",
                                  Waples.correction.value = 46,
                                  naive = T,
                                  outfile = "PMD_LDne_17.txt",
                                  neest.path = "C:/Users/samue/Desktop",
                                  plot_theme = theme_bw())

nes_2018

## 2024 

gl_2024 <- gl.keep.pop(gl2, pop.list = "2024"); gl_2024

gl_2024 <- gl.filter.monomorphs(gl_2024); gl_2024

nes_2024 <- dartR.popgen::gl.LDNe(gl_2024, 
                                  mating = "random",
                                  singleton.rm = TRUE,
                                  Waples.correction = "nChromosomes",
                                  Waples.correction.value = 46,
                                  outfile = "PMD_LDne_17.txt",
                                  neest.path = "C:/Users/samue/Desktop",
                                  plot_theme = theme_bw())

nes_2024

### COLONY Prep

t3 <- gl
t3$other$ind.metrics$offspring <- "yes"
t3$other$ind.metrics$mother <- "no"
t3$other$ind.metrics$father <- "no"
t3$other$ind.metrics$id <- indNames(t3)


gl2colony(gl = t3,
          filename_out = "C:/Users/samue/Desktop/colony_all.txt",
          probability_father = 0.5,
          probability_mother = 0.5,
          seed = 42,                # Seed for random number generator
          update_allele_freq = 0,     # 0/1=Not updating/updating allele frequency
          di_mono_ecious = 2,         # 2/1=Dioecious/Monoecious species
          inbreed = 0,                # 0/1=no inbreeding/inbreeding
          haplodiploid = 0,           # 0/1=Diploid species/HaploDiploid species
          polygamy_male = 0,          # 0/1=Polygamy/Monogamy for males
          polygamy_female = 0,        # 0/1=Polygamy/Monogamy for females
          clone_inference = 1,        # 0/1=Clone inference =No/Yes
          scale_shibship = 1,         # 0/1=Scale full sibship=No/Yes
          sibship_prior = 0,          # 0/1/2/3/4=No/Weak/Medium/Strong/Optimal sibship prior
          known_allele_freq = 0,      # 0/1=Unknown/Known population allele frequency
          num_runs = 10,               # Number of runs
          length_run = 3,             # 1/2/3/4=short/medium/long/very long run
          monitor_method = 0,         # 0/1=Monitor method by Iterate#/Time in second
          monitor_interval = 10000,   # Monitor interval in Iterate# / in seconds
          windows_gui = 1,            # 0/1=No/Yes for run with Windows GUI
          likelihood = 1,             # 0/1/2=PairLikelihood score/Fulllikelihood/FPLS
          precision_fl = 2,           # 0/1/2/3=Low/Medium/High/Very high precision with Full-likelihood
          marker_id = 'mk@',          # Marker Ids
          marker_type = '0@',         # Marker types, 0/1=Codominant/Dominant
          allelic_dropout = '0.000@', # Allelic dropout rate at each locus
          other_typ_err = '0.05@',    # Other typing error rate at each locus
          paternity_exclusion_threshold = '0 0',
          maternity_exclusion_threshold = '0 0',
          paternal_sibship = 0,
          maternal_sibship = 0,
          excluded_paternity = 0,
          excluded_maternity = 0,
          excluded_paternal_sibships = 0,
          excluded_maternity_sibships = 0
)
