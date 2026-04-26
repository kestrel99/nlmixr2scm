# nlmixr2scm

`nlmixr2scm` provides stepwise covariate modeling (SCM) for `nlmixr2`
population PK/PD models.

SCM is a likelihood-ratio-based covariate selection workflow that proceeds in
two phases:

1. Forward inclusion: starting from a base model, add the candidate
   covariate-parameter relationship that gives the largest statistically
   significant OFV improvement.
2. Backward elimination: starting from the forward-final model, test each
   included relationship for removal using a stricter threshold.

The result is a reproducible covariate-search workflow that can summarize every
candidate tested, keep accepted and dropped relationships explicit, and return
the final forward and backward models for follow-up review.

## The details

`runSCM()` implements the SCM workflow for `nlmixr2` fits. It automatically
builds the covariate terms inside the model body, so continuous covariates can
be centered at their observed medians and categorical covariates can be turned
into indicator columns without pre-editing the model by hand.

The package supports:

* full Cartesian searches via `varsVec` and `covarsVec`
* exact candidate definitions via `pairsVec`
* multiple built-in continuous-covariate shapes (`"power"`, `"lin"`, `"log"`,
  and `"identity"`)
* automatic categorical preprocessing through `catvarsVec`
* optional forced backward-start relationships through `includedRelations`
* saved step tables and fitted candidate models for resume/review workflows
* optional parallel fitting of candidates through `workers`

`nlmixr2scm` is designed to work alongside `nlmixr2utils`, which provides 
shared worker-plan helpers and supporting infrastructure.

## Installation

The package is not on CRAN. Install it from GitHub together with
`nlmixr2utils`.

Using `pak`:

```r
pak::pkg_install(c(
  "kestrel99/nlmixr2utils",
  "kestrel99/nlmixr2scm"
))
```

Using `remotes`:

```r
remotes::install_github("kestrel99/nlmixr2utils")
remotes::install_github("kestrel99/nlmixr2scm")
```

## Basic Use

```r
library(nlmixr2)
library(nlmixr2data)
library(nlmixr2scm)

pkdata <- nlmixr2data::warfarin[nlmixr2data::warfarin$dvid == "cp", ]

warf_pk <- function() {
  ini({
    tka <- log(1.15)
    tcl <- log(0.135)
    tv <- log(7.0)
    eta.ka ~ 0.40
    eta.cl ~ 0.25
    eta.v ~ 0.10
    prop.err <- 0.10
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ prop(prop.err)
  })
}

fit_base <- nlmixr2(
  warf_pk,
  data = pkdata,
  est = "focei",
  control = nlmixr2est::foceiControl(print = 0),
  table = nlmixr2est::tableControl(cwres = TRUE)
)

scm <- runSCM(
  fit = fit_base,
  data = pkdata,
  varsVec = c("cl", "v"),
  covarsVec = "wt",
  catvarsVec = "sex",
  shapes = c("power", "lin"),
  pVal = list(fwd = 0.05, bck = 0.01),
  searchType = "scm",
  saveModels = FALSE,
  print = 0,
  workers = 1L,
  confirm = FALSE
)

scm$summaryTable
scm$resBck[[1]]$parFixedDf
```

## Requirements and Practical Notes

Pass the original, untransformed data to `runSCM()`. Centering, indicator
generation, and any other necessary covariate-shape transformations are generated inside the
SCM tool.

For practical use:

* Pass `data` explicitly whenever the base model does not already reference a
  covariate column that SCM needs to test.
* Use `catvarsVec` for categorical predictors instead of creating dummy columns
  manually.
* Use `pairsVec` when you want exact control over which relationships are
  tested, rather than the full parameter x covariate search space.
* Use `includedRelations` when backward elimination should start from a model
  that contains required relationships beyond what forward inclusion selected.
* Use `saveModels = FALSE` for quick exploratory runs; leave it at the default
  when you want saved step tables, fitted candidates, and resume support.
* Use `confirm = FALSE` in scripts, tests, or non-interactive workflows.

## How this differs from `nlmixr2extra`

The underlying forward/backward SCM idea is the same as in `nlmixr2extra`, but
the current package has a broader and more explicit workflow.

Compared with `covarSearchAuto()` in `nlmixr2extra`, `runSCM()`:

* is the sole public SCM entry point, with the legacy alias removed and the
  model-building helpers kept internal
* supports explicit `pairsVec` candidate definitions in addition to the older
  vector-based Cartesian search
* supports multiple built-in shapes, custom shape builders, and richer initial
  value / bounds control
* handles categorical covariates more explicitly through reference-level
  selection, indicator generation, and `catCutoff`
* supports `missingToken`, `includedRelations`, explicit `data`, and optional
  parallel candidate fitting through `workers`
* uses numbered `<fit>_scm_<N>` output directories with step summaries,
  candidate tables, log files, and resume support, while also allowing
  in-memory runs via `saveModels = FALSE`
* can run model fits in parallel via the `workers` argument

## References

* Jonsson, E.N. & Karlsson, M.O. (1998). Automated covariate model building
  within NONMEM. *Pharmaceutical Research*, 15(9), 1463-1468.
* Lindbom, L., Ribbing, J., & Jonsson, E.N. (2004). Perl-speaks-NONMEM (PsN):
  a Perl module for NONMEM related programming. *Computer Methods and Programs
  in Biomedicine*, 75, 85-94.
* Ribbing, J. & Jonsson, E.N. (2004). Power, selection bias and predictive
  performance of the population pharmacokinetic covariate model. *Journal of
  Pharmacokinetics and Pharmacodynamics*, 31(2), 109-134.

For a fuller worked example, see the package vignette:
`vignette("runSCM", package = "nlmixr2scm")`.
