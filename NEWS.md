# nlmixr2scm 0.2

* Auto-scaled default theta bounds for `lin`, `log`, and `identity` covariate shapes. Previously the shape-agnostic defaults of `(-5, 5)` could produce numerically catastrophic parameter spaces for physiologically-scaled covariates (e.g. weight in kg), causing both diverged optimizations and spurious minima. Bounds are now scaled by `1/|center|` for `lin` and `identity`, and by `1/|log(center)|` for `log`, giving results comparable to the dimensionless `power` shape.

* New retry mechanism for unrealistic OFV values in `runSCM()`. Six new arguments control the behaviour: `maxRetries` (default `3L`), `maxDeltaOFV` (default `Inf`), `retryPerturbSD` (default `0.5`), `retrySmallInit` (default `0.01`), `retryOFVTolerance` (default `NULL`, auto-detected), and `retryFailOnExhaustion` (default `FALSE`). When a candidate fit produces an OFV that is implausibly high or low, the fit is retried with perturbed or reduced initial estimates. Stochastic estimators (SAEM) get a 10-unit tolerance margin automatically to avoid spurious retries due to Monte Carlo noise.

# nlmixr2scm 0.1

* `runSCM()` provides a standalone home for nlmixr2 stepwise covariate modeling, including SCM model-building helpers, categorical covariate expansion, resume support, and the migrated SCM vignette/tests.
