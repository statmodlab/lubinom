# lubinom

`lubinom` is an R package for the Lambert-uniform binomial distribution,
a one-parameter model for bounded count data with support `{0, ..., n}`.

The package provides tools for probability evaluation, random generation,
threshold-based reliability measures, maximum likelihood estimation,
goodness-of-fit assessment, and comparison with alternative bounded-count
models.

## Installation

The package can be installed from source using:

```r
install.packages("lubinom_0.1.0.tar.gz", repos = NULL, type = "source")
```

After installation:

```r
library(lubinom)
```

## Basic usage

Evaluate the probability mass function:

```r
dLUB(
  x = 0:4,
  size = 4,
  beta = 3
)
```

Generate random observations:

```r
set.seed(123)

rLUB(
  n = 100,
  size = 4,
  beta = 3
)
```

## Maximum likelihood estimation

```r
x <- rep(
  0:4,
  times = c(69, 30, 9, 2, 0)
)

fit <- fitLUB_mle(
  x = x,
  size = 4
)

fit
```

## Model comparison

The package includes tools for comparison with the binomial,
beta-binomial, Kumaraswamy-binomial, and Conway-Maxwell-binomial models.

```r
compareLUBmodels(
  x = x,
  size = 4
)
```

## Reproducibility scripts

The scripts used to reproduce the figures, simulation study, and data
applications associated with the paper are included with the package.

After installation, their location can be obtained with:

```r
system.file("paper", package = "lubinom")
```

A script can be opened for inspection without executing it, for example:

```r
file.edit(
  system.file(
    "paper",
    "Application1.R",
    package = "lubinom"
  )
)
```

## License

This package is distributed under the GNU General Public License version 3.
