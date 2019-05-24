library(dplyr)


# sample IDs --------------------------------------------------------------

gen_samples <- function(n, snps) {
  rep(1:n, each = snps)
}


# SNP IDs -----------------------------------------------------------------

gen_snps <- function(n, samples) {
  rep(1:n, times = samples)
}


# SNP alleles -------------------------------------------------------------

gen_alleles <- function(n) {
  as.integer(rbinom(n, size = 1, prob = 0.5))
}


# chromosomal crossover ---------------------------------------------------

crossover <- function(alleles1, alleles2) {
  len <- length(alleles1)
  chiasma <- sample.int(len - 1, 1)
  tibble(
    a1 = c(alleles1[1:chiasma], alleles2[(chiasma + 1):len]),
    a2 = c(alleles2[1:chiasma], alleles1[(chiasma + 1):len])
  )
}


# maternal allele count ---------------------------------------------------

maternal_allele_n <- function(n, coverage, fetal_fraction, vmr) {
  maternal_fraction <- 1 - fetal_fraction
  mean <- coverage * maternal_fraction / 2
  dispersion_param <- coverage / (vmr - 1)  # ReadDepth (http://dx.doi.org/10.1371/journal.pone.0016327)
  rnbinom(n, mu = mean, size = dispersion_param)
}


# fetal allele count ------------------------------------------------------

fetal_allele_n <- function(n, coverage, fetal_fraction, vmr) {
  mean <- coverage * fetal_fraction / 2
  dispersion_param <- coverage / (vmr - 1)  # ReadDepth (http://dx.doi.org/10.1371/journal.pone.0016327)
  rnbinom(n, mu = mean, size = dispersion_param)
}


# genotypes ---------------------------------------------------------------

gen_genotypes <- function(samples, snps) {
  gens <- tibble(
    sample = gen_samples(samples, snps),  # sample ID
    snp = gen_snps(snps, samples),  # SNP ID
    ma1 = gen_alleles(samples * snps),  # maternal allele1
    ma2 = gen_alleles(samples * snps),  # maternal allele2
    pa1 = gen_alleles(samples * snps),  # paternal allele1
    pa2 = gen_alleles(samples * snps)  # paternal allele2
    )

  gens <- bind_cols(
    gens,
    gens %>%
    group_by(sample) %>%
      do(crossover(.$ma1, .$ma2)) %>%
      ungroup() %>%
      select(
        fa_mat = a1,  # maternally inherited fetal allele1
        fa_mat_tri = a2  # maternally inherited fetal allele2 (trisomy)
      ),
    gens %>%
    group_by(sample) %>%
      do(crossover(.$pa1, .$pa2)) %>%
      ungroup() %>%
      select(
        fa_pat = a1,  # paternally inherited fetal allele1
        fa_pat_tri = a2  # paternally inherited fetal allele2 (trisomy)
      )
  )
  gens
}


# allele counts -----------------------------------------------------------

gen_counts <- function(gens, coverage, fetal_fraction, vmr) {
  n <- nrow(gens)
  tibble(
    sample = gens$sample,  # sample ID
    snp = gens$snp,  # SNP ID
    ma1_n = maternal_allele_n(n, coverage, fetal_fraction, vmr),  # maternal allele1 count
    ma2_n = maternal_allele_n(n, coverage, fetal_fraction, vmr),  # maternal allele2 count
    fa_mat_n = fetal_allele_n(n, coverage, fetal_fraction, vmr),  # maternally inherited fetal allele count
    fa_pat_n = fetal_allele_n(n, coverage, fetal_fraction, vmr),  # paternally inherited fetal allele count
    fa_tri_n = fetal_allele_n(n, coverage, fetal_fraction, vmr)  # trisomic fetal allele count
  )
}

# normal ------------------------------------------------------------------

sim_euploidy <- function(gens, counts) {
  tibble(
    sample = gens$sample,  # sample ID
    snp = gens$snp,  # SNP ID
    a0_n =  # wild-type allele count
      if_else(gens$ma1 == 0, counts$ma1_n, 0) +
      if_else(gens$ma2 == 0, counts$ma2_n, 0) +
      if_else(gens$fa_mat == 0, counts$fa_mat_n, 0) +
      if_else(gens$fa_pat == 0, counts$fa_pat_n, 0),
    a1_n =  # alternative allele count
      if_else(gens$ma1 == 1, counts$ma1_n, 0) +
      if_else(gens$ma2 == 1, counts$ma2_n, 0) +
      if_else(gens$fa_mat == 1, counts$fa_mat_n, 0) +
      if_else(gens$fa_pat == 1, counts$fa_pat_n, 0),
    mat_gen = gens$ma1 + gens$ma2,  # maternal genotype
    fet_gen = gens$fa_mat + gens$fa_pat,  # fetal genotype
    major_count = pmax(a0_n, a1_n),  # major allele count
    minor_count = pmin(a0_n, a1_n),  # minor allele count
    total_count = a0_n + a1_n,  # total allele count
    allelic_ratio = major_count / minor_count
  )
}


# maternal trisomy --------------------------------------------------------

sim_maternal_trisomy <- function(gens, counts) {
  tibble(
    sample = gens$sample,  # sample ID
    snp = gens$snp,  # SNP ID
    a0_n =  # wild-type allele count
      if_else(gens$ma1 == 0, counts$ma1_n, 0) +
      if_else(gens$ma2 == 0, counts$ma2_n, 0) +
      if_else(gens$fa_mat == 0, counts$fa_mat_n, 0) +
      if_else(gens$fa_pat == 0, counts$fa_pat_n, 0) +
      if_else(gens$fa_mat_tri == 0, counts$fa_tri_n, 0),
    a1_n =  # alternative allele count
      if_else(gens$ma1 == 1, counts$ma1_n, 0) +
      if_else(gens$ma2 == 1, counts$ma2_n, 0) +
      if_else(gens$fa_mat == 1, counts$fa_mat_n, 0) +
      if_else(gens$fa_pat == 1, counts$fa_pat_n, 0) +
      if_else(gens$fa_mat_tri == 1, counts$fa_tri_n, 0),
    mat_gen = gens$ma1 + gens$ma2,  # maternal genotype
    fet_gen = gens$fa_mat + gens$fa_pat + gens$fa_mat_tri,  # fetal genotype
    major_count = pmax(a0_n, a1_n),  # major allele count
    minor_count = pmin(a0_n, a1_n),  # minor allele count
    total_count = a0_n + a1_n,  # total allele count
    allelic_ratio = major_count / minor_count
  )
}


# paternal trisomy --------------------------------------------------------

sim_paternal_trisomy <- function(gens, counts) {
  tibble(
    sample = gens$sample,  # sample ID
    snp = gens$snp,  # SNP ID
    a0_n =  # wild-type allele count
      if_else(gens$ma1 == 0, counts$ma1_n, 0) +
      if_else(gens$ma2 == 0, counts$ma2_n, 0) +
      if_else(gens$fa_mat == 0, counts$fa_mat_n, 0) +
      if_else(gens$fa_pat == 0, counts$fa_pat_n, 0) +
      if_else(gens$fa_pat_tri == 0, counts$fa_tri_n, 0),
    a1_n =  # alternative allele count
      if_else(gens$ma1 == 1, counts$ma1_n, 0) +
      if_else(gens$ma2 == 1, counts$ma2_n, 0) +
      if_else(gens$fa_mat == 1, counts$fa_mat_n, 0) +
      if_else(gens$fa_pat == 1, counts$fa_pat_n, 0) +
      if_else(gens$fa_pat_tri == 1, counts$fa_tri_n, 0),
    mat_gen = gens$ma1 + gens$ma2,  # maternal genotype
    fet_gen = gens$fa_mat + gens$fa_pat + gens$fa_pat_tri,  # fetal genotype
    major_count = pmax(a0_n, a1_n),  # major allele count
    minor_count = pmin(a0_n, a1_n),  # minor allele count
    total_count = a0_n + a1_n,  # total allele count
    allelic_ratio = major_count / minor_count
  )
}
