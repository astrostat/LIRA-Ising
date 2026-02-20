# LIRA-Ising: Estimating the boundaries of irregularly shaped sources

Bayesian workflow to identify the boundaries of irregularly shaped sources. First, LIRA (pylira) is run to obtain posterior multiscale counts that account for known structures, background, and the telescope PSF. This package then uses an Ising model to group pixels with similar intensities into background and source regions, and extracts a boundary from the most likely pixel aggregation.

## Paper

This repository reproduces the method described in:
https://iopscience.iop.org/article/10.3847/1538-4365/ae3716

## What this repo provides

- LIRA post-processing to reshape `.out` + `.param` files into an iteration-by-pixel matrix.
- Ising Gibbs sampler with Swendsen-Wang updates to infer pixel assignments.
- Boundary estimation (`getBound`) and optional genetic algorithm alternative (`gaBound`).
- Plot helpers to visualize images and boundaries.

## Requirements

- R with packages: `ggplot2`, `reshape2`, `RColorBrewer`, `foreach`, `parallel`, `doMC`, `dplyr`, `Brobdingnag`, `FITSio`, `igraph`, `coda`, `gmp`, `MASS`.
- LIRA (pylira) output files: an `.out` file with posterior draws and a `.param` file with MCMC diagnostics.
- A partition function file for the image size `n` (downloadable from Paul Beale's Ising model tables):
	https://spot.colorado.edu/~beale/index1.html

## Prerequisite: pylira

This repo expects you to run LIRA using `pylira` first. The upstream project is:
https://github.com/astrostat/pylira

If you have cloned it into `pylira/` (for example at the same level as this repo), follow the pylira README to install its dependencies and run the LIRA model. At a minimum, you will need a working Python environment and the dependencies listed in pylira's documentation.

### What pylira produces

LIRA produces at least two files used here:

- `*.out`: posterior draws of the multiscale counts (one draw per iteration, all pixels concatenated).
- `*.param`: MCMC diagnostics including `logPost`, `expectedMSCounts`, and possibly `bkgScale`.

### What LIRA-Ising consumes

This package consumes those LIRA outputs via `liraPost()`:

- `out` is the `*.out` file.
- `param` is the `*.param` file.

`liraPost()` returns a matrix with one row per LIRA draw and one column per pixel. That matrix (`lambda`) is the primary input to `isingGibbs()` and `getBound()`.

## Installation

From the repo root:

```bash
R CMD INSTALL .
```

Or from R:

```r
install.packages(c(
	"ggplot2", "reshape2", "RColorBrewer", "foreach", "parallel", "doMC",
	"dplyr", "Brobdingnag", "FITSio", "igraph", "coda", "gmp", "MASS"
))
install.packages("devtools")
devtools::install_local(".")
```

## Quickstart

```r
library(LIsegmentation)

# 1) Read LIRA output
lambda <- liraPost(
	out = "path/to/lira_output.out",
	param = "path/to/lira_output.param",
	burn = 100,
	plot = TRUE
)

# 2) Load partition function for image size n x n
n <- sqrt(ncol(lambda))
G <- loadPartition(n = n, g_file = "path/to/beale_g_table.txt")

# 3) Run Ising model MCMC
ising <- isingGibbs(
	lambda = lambda,
	G = G,
	init_iter = 500,
	burn_iter = 50,
	ncores = 4
)

# 4) Extract a boundary (MAP over an ad hoc candidate set)
z_max <- getBound(Ziter = ising$ising_array, lambda = lambda, param = ising$param)
z_img <- array(z_max, dim = c(n, n))

# 5) Optional: probability map of source membership
p_img <- cAvgIsing(ising$ising_array)

# 6) Plot (img can be observed counts or a posterior mean image)
img_mean <- array(colMeans(lambda), dim = c(n, n))
plotSource(img = img_mean, bound = z_img, title = "LIRA-Ising boundary")
```

## Toy example (no pylira)

This example builds a synthetic image, simulates LIRA-like draws, and runs the
Ising boundary estimation. It demonstrates the complete workflow and API without requiring LIRA output files.

A complete standalone script is available in [lira_ising_toy.R](lira_ising_toy.R). You can run it directly:

```bash
Rscript --vanilla lira_ising_toy.R
```

Or copy and paste the code below into an R session:

```r
library(LIsegmentation)

set.seed(123)

# 1) Build a synthetic source image
n <- 16
img_true <- rectBase(
	bkg_param = 2,
	max_param = 8,
	n_img = n,
	width = 6,
	height = 6
)

# 2) Simulate LIRA-like draws (rows = draws, cols = pixels)
n_draws <- 100
lambda <- t(replicate(n_draws, c(img_true + matrix(rnorm(n^2, 0, 1), n, n))))

# 3) Load partition function for this n
G <- loadPartition(n = n, g_file = "path/to/beale_g_table_for_n16.txt")

# 4) Run Ising model MCMC (small iteration counts for speed)
ising <- isingGibbs(
	lambda = lambda,
	G = G,
	init_iter = 50,
	burn_iter = 10,
	beta_niter = 10,
	init_seed = 123,
	ncores = 1
)

# 5) Extract a boundary and plot
z_max <- getBound(Ziter = ising$ising_array, lambda = lambda, param = ising$param)
z_img <- array(z_max, dim = c(n, n))

img_mean <- array(colMeans(lambda), dim = c(n, n))
plotSource(img = img_mean, bound = z_img, title = "Toy example boundary")
```

This example uses a smaller image (n=16) and fewer iterations than typical production runs for speed. The complete version is saved in [lira_ising_toy.R](lira_ising_toy.R).

## Notes on inputs and outputs

- All images are assumed to be square `n x n`. The partition function file must match this `n`.
- `liraPost()` returns one row per LIRA draw and one column per pixel.
- `isingGibbs()` returns a list with `param` (tau and beta draws) and `ising_array` (pixel assignments).
- `getBound()` returns a vector of length `n^2` with values `-1` (background) and `1` (source).

## Suggested workflow

1) Run LIRA (pylira) on your data and save the `.out` and `.param` files.
2) Use `liraPost()` to reshape the LIRA draws.
3) Load the partition function for the image size using `loadPartition()`.
4) Run `isingGibbs()` to sample pixel assignments.
5) Extract a boundary with `getBound()` or use `gaBound()` for a genetic algorithm alternative.
6) Visualize with `plotSource()` or `plotGrid()`.
