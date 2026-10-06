# freeiv — instrument-free estimation of a linear structural model

`freeiv` is a Stata command (Stata 16 or later) that estimates the coefficient
on an endogenous regressor **without an external instrument**, in the linear
triangular model in which a single latent confounder generates the
endogeneity. Version **1.0.0** (6 October 2026).

The assumption-free result comes first. Under scale consistency — the
confounder moves the outcome directly by exactly as much as it moves it
through the regressor — the second moments of the data leave one number free,
and the coefficient lies in

```
    [ gamma-tilde / 2 ,  gamma-tilde ]        gamma-tilde = the OLS slope
```

That interval assumes nothing more, and it is exactly the set of values for
which the implied variances are non-negative: an estimate outside it is not a
large estimate but an infeasible one, and the command says so. The routes of
the command then impose the same linear structure on the moments of the next
orders — the third order identifies the coefficient, the fourth
over-identifies it and makes the structure testable — and the estimators of
the literature are computed on the same data and judged by the same interval.
A second indicator of the confounder frees the confounder's direct loading and
turns scale consistency into a test. Standard errors come from the influence
functions of the moments, with sampling weights in sandwich form and
linearization over a survey design.

The models, the command and worked examples are presented in

> Araar, A. (2026). *freeiv: Instrument-free estimation of a linear structural
> model with an endogenous regressor* (Version 1.0.0). Zenodo.
> <https://doi.org/10.5281/zenodo.23190841>

which is [`paper/freeiv_paper.pdf`](paper/freeiv_paper.pdf) in this repository.
All versions of the paper: <https://doi.org/10.5281/zenodo.22770175> (always the
latest).

## Installation

From Stata:

```stata
net install freeiv, from("https://raw.githubusercontent.com/aabbdd12/freeiv/main/") replace
net get freeiv, from("https://raw.githubusercontent.com/aabbdd12/freeiv/main/") replace
help freeiv
```

`net get` copies the four example datasets and the paper into the current
folder. The examples of `help freeiv` run from their links — in the command
window, in the dialog box filled in, or as a do-file — and, without the copy,
read their data from the SSC archive or from this repository. The package will
also be submitted to the SSC archive (`ssc install freeiv, all`).

## Quick start

```stata
use freeiv_card, clear
freeivmenu lwage exper expersq black south smsa (educ)              // what can these data carry?
freeiv     lwage exper expersq black south smsa (educ)              // the default route, with the interval
freeiv     lwage exper expersq black south smsa (educ), method(all) // every route against the interval
```

The output is read in four steps: the menu, before any estimation, says
which routes the data can carry; the interval holds whatever route is chosen;
the routes follow, each against the interval; the diagnostics last — the
precision of the third-order signal, the share of the interval that the
fourth-order moments rule out, the tests of the over-identifying restrictions.

## The routes

| `method()` | what it adds to the linear model | the signal that says whether the data carry it |
|---|---|---|
| `bounds` | nothing: the interval | — |
| `ols` | nothing: the upper end of the interval | — |
| `qme` *(default)* | a third-order signal (the quadratic moment estimator) | z of m03, z of the discriminant |
| `hme` | V1 and V2 symmetric (a closed form) | B = E[V2^3] |
| `gmm` | the fourth-order identities: nine moments, a Hansen J on 1 df, and the region where the profiled J stays within 3.84 of its minimum | the region's share of the interval, J |
| `pgmm` | the same moments with orders 2 and 3 fitted exactly: the companion of `gmm` | J |
| `lsz` | Lewbel, Schennach and Zhang (2024), as `trigmm` implements it | convergence |
| `lewbel12` | Lewbel (2012): heteroskedasticity in the controls | F of eps2^2 on X |
| `copula` | a Gaussian copula (Park and Gupta 2012) | — |
| `rank` | the rank control function (Breitung, Mayer and Wied 2024) | — |
| `oster` | proportional selection (Oster 2019), delta and Rmax assumed | — |
| `all` | every route above, in one table | |

**Read the region's share before the point estimate.** When it fills the
interval, the moments of orders three and four have added nothing to the
second order, the point is the argmin of a flat criterion, and the interval is
what to report.

## Two indicators of one confounder

Two variables in the parentheses switch the command to the two-indicator
model: both load on the same latent confounder, the confounder's direct
loading in the outcome equation is free and identified in closed form, and
scale consistency becomes a test.

```stata
freeivmenu lwage exper expersq black south smsa (educ motheduc)   // relevance, guards, precision regime
freeiv     lwage exper expersq black south smsa (educ motheduc)
freeivtest                                                        // scale consistency, one factor, J
```

When the moments contradict one factor, a guard refuses the closed form, the
command ends with r(498) and `e(guard)` says which guard fired.

## Companion commands

| command | what it does |
|---|---|
| `freeivmenu` | which routes the data can carry, **before** any estimation; with two indicators, what the two-indicator model needs |
| `freeivdiag` | what a given value of the coefficient implies about the unobservables, whatever produced it |
| `freeivtest` | endogeneity itself, the agreement between routes, an outside estimate against the interval; after the two-indicator model, its tests |
| `freeivreport` | all of it in one table, optionally written to a file |

`db freeiv` opens the dialog box.

## The datasets

| dataset | what it is for |
|---|---|
| `freeiv_sim1` | simulated, the model holds, the truth is 0.40: the routes can be judged against it |
| `freeiv_sim2` | the same design with V2 drawn from the confounder's skewed law: the two roots of the quadratic merge and the symmetry route fails quietly |
| `freeiv_proxy` | simulated with two indicators of one confounder |
| `freeiv_card` | the Card (1995) extract distributed with Wooldridge's teaching datasets |

The generating equations of the simulated datasets are printed in
`help freeiv`.

## Reference commands, not redistributed

The comparison routes were validated against the commands that implement them;
they are not in this repository. `method(lsz)` reproduces `trigmm`
(`ssc install trigmm`) within its stopping tolerance, `method(lewbel12)`
reproduces `ivreg2h` (SSC) and `method(oster)` reproduces `psacalc` (SSC).

## Repository layout

| path | content |
|---|---|
| `freeiv.pkg`, `stata.toc` | the package description read by `net install` and `net get` |
| `src/` | the commands (`.ado`), their help files (`.sthlp`) and the dialog box (`.dlg`) |
| `examples/` | the four example datasets |
| `paper/freeiv_paper.pdf` | the paper of the command |
| `replication/` | `freeiv_paper/freeiv_note.do`, which reproduces every output and every number of the paper (see `replication/README.md`) |

## Papers of the series

The methods are developed in four working papers on Zenodo, cited by their
concept DOI (all versions, resolving to the latest):

- Araar, A. (2026a). *Correcting Endogeneity without External Instruments:
  Four Closed-Form Estimators for the Linear Structural Model.*
  [10.5281/zenodo.20312356](https://doi.org/10.5281/zenodo.20312356)
- Araar, A. (2026b). *The Quadratic Moment Estimator for Instrument-Free
  Endogeneity Correction.*
  [10.5281/zenodo.22068003](https://doi.org/10.5281/zenodo.22068003)
- Araar, A. (2026c). *Instrument-Free Estimation under Linear
  Scale-Consistency: Closed-Form Identification, Partial Bounds, and
  Higher-Moment GMM.*
  [10.5281/zenodo.22119231](https://doi.org/10.5281/zenodo.22119231)
- Araar, A. (2026d). *Two Indicators of One Latent Confounder: Closed-Form
  Identification of the Triangular Model with a Free Proxy Effect.*
  [10.5281/zenodo.22207331](https://doi.org/10.5281/zenodo.22207331)

## Citing

```
Araar, A. (2026). freeiv: Instrument-free estimation of a linear structural
model with an endogenous regressor (Version 1.0.0). Zenodo.
https://doi.org/10.5281/zenodo.23190841
```

To cite all versions, use <https://doi.org/10.5281/zenodo.22770175>, which
always resolves to the latest. A `CITATION.cff` file is included for reference
managers.

## Versions

Version 1.0.0 is on `main` and at the tag `v1.0.0`. The earlier version 0.8.2
(15 September 2026) stays installable from its tag:

```stata
net install freeiv, from("https://raw.githubusercontent.com/aabbdd12/freeiv/v0.8.2/") replace
```

## Author

Abdelkrim Araar, Université Laval and the Partnership for Economic Policy —
<aabd@ecn.ulaval.ca>

## Licence

MIT, for the code in this repository. The reference commands named above are
under their own licences and are not included.
