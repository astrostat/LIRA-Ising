req_pkgs <- c(
  "ggplot2", "reshape2", "RColorBrewer", "foreach", "parallel", "doMC",
  "dplyr", "Brobdingnag", "FITSio", "igraph", "coda", "gmp", "MASS"
)
missing <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(paste("Missing required R packages:", paste(missing, collapse = ", ")))
}

suppressPackageStartupMessages({
  library(foreach)
  library(doMC)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  library(dplyr)
  library(Brobdingnag)
  library(gmp)
  library(igraph)
})

source_files <- list.files("R", full.names = TRUE)
invisible(lapply(source_files, source))

set.seed(123)

n <- 8
img_true <- rectBase(
  bkg_param = 2,
  max_param = 8,
  n_img = n,
  width = 4,
  height = 4
)

n_draws <- 20
lambda <- t(replicate(n_draws, c(img_true + matrix(rnorm(n^2, 0, 1), n, n))))

n_g <- n^2 + 1
g_file <- tempfile(fileext = ".txt")
writeLines(paste(rep("1", n_g), collapse = ","), g_file)

G <- loadPartition(n = n, g_file = g_file)

ising <- isingGibbs(
  lambda = lambda,
  G = G,
  init_iter = 20,
  burn_iter = 5,
  beta_niter = 5,
  init_seed = 123,
  ncores = 1
)

z_max <- getBound(Ziter = ising$ising_array, lambda = lambda, param = ising$param)
print(length(z_max))
