source("sim.R")


library(dplyr)
library(feather)


samples <- 10000
snps <- 1000  # SNPs per sample
covs <- seq(500, 15000, by = 500)  # mean coverage
ffs <- seq(0.01, 0.2, by = 0.01)  # fetal fraction
vmr <- 3  # variance-to-mean ratio

set.seed(42)

gens <- gen_genotypes(samples, snps)

path <- "data"
dir.create(path)
write_feather(gens, sprintf("%s/gens.feather", path))

for (cov in covs) {
  cat("coverage:", cov, "\n")
  for (ff in ffs) {
    cat("fetal fraction:", ff, "\n")
    counts <- gen_counts(gens, cov, ff, vmr)

    # write data
    out_path <- sprintf("%s/coverage/%s/fetal_fraction/%s", path, cov, ff)
	dir.create(out_path, recursive = T)
    write_feather(counts, sprintf("%s/counts.feather", out_path))
    write_feather(sim_euploidy(gens, counts), sprintf("%s/euploidy.feather", out_path))
    write_feather(sim_maternal_trisomy(gens, counts), sprintf("%s/maternal_trisomy.feather", out_path))
    write_feather(sim_paternal_trisomy(gens, counts), sprintf("%s/paternal_trisomy.feather", out_path))
  }
}
