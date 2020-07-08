source("sim.R")


library(arrow)
library(glue)


samples <- 100
snps <- 1000  # SNPs per sample
covs <- seq(500, 15000, by = 500)  # mean coverage
ffs <- seq(0.01, 0.2, by = 0.01)  # fetal fraction
vmr <- 3  # variance-to-mean ratio

set.seed(42)

gens <- gen_genotypes(samples, snps)

set_path <- "train_data"
dir.create(set_path)
write_feather(gens, glue("{set_path}/gens.feather"))

for (cov in covs) {
  cat("coverage:", cov, "\n")
  for (ff in ffs) {
    cat("fetal fraction:", ff, "\n")
    counts <- gen_counts(gens, cov, ff, vmr)

    # write to file
    file_path <- glue("{set_path}/coverage/{cov}/fetal_fraction/{ff}")
    dir.create(file_path, recursive = T)
    write_feather(counts, glue("{file_path}/counts.feather"))
    write_feather(sim_euploidy(gens, counts),
                  glue("{file_path}/euploidy.feather"))
    write_feather(sim_maternal_trisomy(gens, counts),
                  glue("{file_path}/maternal_trisomy.feather"))
    write_feather(sim_paternal_trisomy(gens, counts),
                  glue("{file_path}/paternal_trisomy.feather"))
  }
}
