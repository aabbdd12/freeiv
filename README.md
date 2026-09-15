# freeiv — instrument-free estimation of a linear structural model

`freeiv` estimates the coefficient on an endogenous regressor **without an
external instrument**, for the standard linear triangular model in which the
endogeneity comes from a single latent confounder.

Its organising idea is that the assumption-free result comes first. Under
scale consistency alone — the confounder's direct loading in the outcome
equation is proportional to its loading in the first stage — the coefficient
is bounded by

```
    [ gamma-tilde / 2 ,  gamma-tilde ]        gamma-tilde = the OLS slope
```

and that interval is *exactly* the region where the implied variances are
non-negative. No distributional assumption is used to obtain it. Every point
estimator then adds one assumption to close the system, so the package
reports them side by side against the interval, with the diagnostic that says
whether each one's assumption is tenable on these data.

An estimate that falls outside the interval is not a large estimate. It is an
infeasible one: it implies a negative variance under the maintained model.
`freeiv` says so rather than printing it without comment.

## Install

```stata
net install freeiv, from("https://raw.githubusercontent.com/aabbdd12/freeiv/main") replace
net get     freeiv, from("https://raw.githubusercontent.com/aabbdd12/freeiv/main") replace
```

`net get` is a separate step and it matters: the four example datasets are
ancillary files, so `net install` alone does not bring them down.

Stata 16 or later. No dependencies.

## Documentation

`help freeiv` is the reference. [`freeiv_paper.pdf`](freeiv_paper.pdf), in this
repository and retrieved by `net get freeiv`, presents the models, the command and worked
examples at article length, in the format of a Stata Journal article; the
output it prints is that of the installed version.

## Quick start

```stata
net get freeiv
use freeiv_card, clear

freeivmenu lwage exper black south smsa (educ)   // which routes can these data carry?
freeiv     lwage exper black south smsa (educ)   // the default route, with the interval
freeiv     lwage exper black south smsa (educ), method(all)
```

Read the menu before the estimate. It reports, from the data alone, whether
the third-order signal exists, whether the symmetry the HME needs is
plausible, and whether the loading is strong enough for the fourth-order
route to bind — that is, which of the assumptions below you are entitled to
make.

## The routes

Section 2 of [`freeiv_paper.pdf`](freeiv_paper.pdf) states each model — what it
assumes, what it identifies, and how the command estimates it.

| `method()` | what it adds to scale consistency |
|---|---|
| `bounds` | nothing — the interval alone |
| `ols` | the long-regression slope, reported as the upper bound it is |
| `qme` *(default)* | a non-zero third moment of the confounder |
| `sce` | equal variances of the two disturbances |
| `rre` | the scale-free counterpart of that restriction |
| `hme` | a symmetric first-stage disturbance |
| `qbe` | quantile bias extrapolation, calibrated |
| `gmm` | the joint nine-moment GMM, with a Hansen J on 1 df |
| `pgmm` | the same nine moments, profiled — a *different* estimator |
| `lsz` | Lewbel, Schennach and Zhang (2024), as `trigmm` implements it |
| `lewbel12` `copula` `rank` `rpiv` `ape` `oster` | the literature, on the same data, for comparison |
| `all` | every route above, in one table |

### The two GMM routes

`method(gmm)` is the estimator of Araar (2026c): nine moment conditions in
eight free parameters, so one over-identifying restriction. It is minimised
**deterministically**, and that is not housekeeping. Fix `gamma` *and*
`theta` and the only non-linearity left — the product `theta·sigma2_V2` —
disappears, so all nine residuals become linear in the six remaining
parameters. The eight-parameter problem is a two-dimensional surface, and a
surface can be gridded. No seed, no starting value, same data same number.
Multi-start over the eight parameters returns a *local* minimum on six of the
eight datasets of that paper.

Because the criterion is often nearly flat, the route also returns the range
of `gamma` over which the profiled Hansen statistic stays within 3.84 of its
minimum, together with that range's share of the identified interval
(`e(ar_lo)`, `e(ar_hi)`, `e(ar_frac)`). This inverts the J at each fixed
`gamma` rather than inverting a Wald statistic around the optimum, so it
remains valid on the boundary and on a flat criterion — the two situations in
which a standard error means nothing. **Read that share before the point
estimate.** When it fills the whole interval, the third- and fourth-order
moments have added nothing to the assumption-free bounds.

`method(pgmm)` forces the first five conditions to zero exactly, which is a
different estimator — it weights them infinitely rather than optimally. It is
kept because it is one-dimensional, hence gridded exhaustively and never
troubled by local minima, and because its fifth residual *is* the QME
quadratic, which exhibits the third-order route as the part the fourth
moments reweight. A wide gap between the two is a symptom of a flat
criterion, not an error in either.

## Two endogenous regressors

Two variables in the parentheses switch the command to the **two-indicator
model**, where both load on one latent confounder and the confounder's direct
loading in the outcome equation is left *free* and identified. That is what
makes scale consistency itself testable rather than maintained.

```stata
freeiv lwage exper black south smsa (educ motheduc)
```

No method name is needed: the syntax alone selects it.

## Companion commands

| command | what it does |
|---|---|
| `freeivmenu` | which routes the data can carry, **before** any estimation |
| `freeivdiag` | what a given value of gamma implies about the unobservables, whatever produced it |
| `freeivtest` | endogeneity, agreement between routes, and an outside estimate against the interval |
| `freeivreport` | all of it in one table, with an optional CSV |

`db freeiv` opens a dialog covering every route and option.

## The datasets

Four travel with the package and come down with `net get freeiv`. Three are
simulated and their generating equations are printed in `help freeiv`, so
none is a black box.

| dataset | what it is for |
|---|---|
| `freeiv_sim1` | the model holds, the truth is 0.40, every route can be judged against it |
| `freeiv_sim2` | identical but for a skewed first-stage disturbance, which breaks the symmetry route and leaves the others intact |
| `freeiv_proxy` | two endogenous regressors loading on one confounder |
| `freeiv_card` | the Card (1995) extract of Araar (2026d) — one endogenous regressor can only be bounded there, while a second indicator identifies it |

## Reference commands, not redistributed

The package was validated against three published commands. They are not in
this repository, because they are not ours to license. Install them yourself
to compare:

- `trigmm` — `ssc install trigmm`, or `st0797` from
  http://www.stata-journal.com/software/sj26-1: the reference for
  `method(lsz)`
- `ivreg2h` (SSC): the reference for `method(lewbel12)`
- `psacalc` (SSC): the reference for `method(oster)`

On the same data, `lewbel12` equals `ivreg2h` to six decimals, `oster`
equals `psacalc` to five, and `lsz` matches `trigmm` to ten.

What each comparison route maintains in place of the model, and why an
estimate from one of them can fall outside the identified interval, is set
out in Section 2.5 and Table 2 of [`freeiv_paper.pdf`](freeiv_paper.pdf); the
Card example in Section 4.3 shows four of the six doing exactly that.

## Papers

The methods implemented here are developed in four working papers, all on
Zenodo.

- Araar, A. (2026a). *Correcting Endogeneity without External Instruments:
  Four Closed-Form Estimators for the Linear Structural Model.*
  [10.5281/zenodo.22067980](https://doi.org/10.5281/zenodo.22067980)
- Araar, A. (2026b). *The Quadratic Moment Estimator for Instrument-Free
  Endogeneity Correction.*
  [10.5281/zenodo.22068143](https://doi.org/10.5281/zenodo.22068143)
- Araar, A. (2026c). *Instrument-Free Estimation under Linear
  Scale-Consistency: Closed-Form Identification, Partial Bounds, and
  Higher-Moment GMM.*
  [10.5281/zenodo.22753299](https://doi.org/10.5281/zenodo.22753299)
- Araar, A. (2026d). *Two Indicators of One Latent Confounder: Closed-Form
  Identification of the Triangular Model with a Free Proxy Effect.*
  [10.5281/zenodo.22764001](https://doi.org/10.5281/zenodo.22764001)

`method(lewbel12)` follows Lewbel (2012) as `ivreg2h` implements it,
`method(lsz)` follows Lewbel, Schennach and Zhang (2024) as `trigmm` does,
and `method(oster)` follows Oster (2019) as `psacalc` does.

## Citing

```
Araar, A. (2026). freeiv: Instrument-free estimation of a linear structural
model. Stata package, version 0.8.1. https://github.com/aabbdd12/freeiv
```

A `CITATION.cff` file is included for reference managers.

## Author

Abdelkrim Araar, Université Laval and the Partnership for Economic Policy —
<aabd@ecn.ulaval.ca>

## Licence

MIT, for the code in this repository. The external commands named above are
under their own licences and are not included.
