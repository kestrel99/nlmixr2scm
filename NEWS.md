# nlmixr2scm 0.1

* `runSCM()` now keeps SCM model-building helpers internal; the former exported helpers were renamed to `scmAddCatCovariates()`, `scmAddOrRemoveCovariate()`, `scmBuildCovInfo()`, and `scmBuildUpdatedUi()`, and the deprecated `covarSearchAuto()` alias was removed.

* `runSCM()` provides a standalone home for nlmixr2 stepwise covariate modeling, including SCM model-building helpers, categorical covariate expansion, resume support, and the migrated SCM vignette/tests.
