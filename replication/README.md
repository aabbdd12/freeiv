# Replication of the paper of freeiv 1.0.0

`freeiv_paper/freeiv_note.do` reproduces every output and every number of the
paper of the command, *freeiv: Instrument-free estimation of a linear
structural model with an endogenous regressor* (`paper/freeiv_paper.pdf`;
Zenodo, <https://doi.org/10.5281/zenodo.22770175>, all versions): each `stlog`
block of the paper is an excerpt of its log, commands included, at the paper's
line width, and each number quoted in the text is printed by it.

Run it from the root of this repository:

    do replication/freeiv_paper/freeiv_note.do

It loads the commands of `src/`, so it runs the code of this version whether or
not freeiv is installed, and reads the example datasets from `examples/`, as a
reader reads them after `net get freeiv`. About four minutes. nlsw88 comes
with Stata. The two parts that need an internet connection come last: nhanes2
(`webuse`) and wage1 (the Wooldridge datasets at Boston College), the second
with a check of the LSZ draws against `trigmm` (`ssc install trigmm`), skipped
when it is not installed. The do-file writes `freeiv_note.log` in the current
folder; `logs/freeiv_note.log` is the log of the run that produced the paper
(Stata 19.5, 6 October 2026).

| section of the paper | what the do-file runs |
|---|---|
| 3.6 The bootstrap | 500 replications of the qme on freeiv_sim1; LSZ on wage1, its 50 replications, and the same 50 draws given to `trigmm` |
| 4.1 Where the model holds | `freeivmenu` and `method(all)` on freeiv_sim1 |
| 4.2 The edge of identification | `method(all)` on freeiv_sim2 |
| 4.3 Real data with one regressor | the menu and `method(all)` on Card; the qme with fewer controls, and `freeivdiag` |
| 4.4 An estimate from elsewhere | `freeivtest` of Card's IV estimate against the interval |
| 4.5 A second indicator | the two-indicator menu, model B and `freeivtest` with motheduc; fatheduc; the guards with IQ and KWW |
| 4.6 A model the data refute | nlsw88: the qme outside the interval, the two J |
| 4.7 Survey data | nhanes2 under `vce(svy)`, against the sampling weights alone |
| 4.8 Reporting a range | Oster's range on Card; factor variables against hand-made terms |
| 5 Discussion | the outcome equation y1 = y2 + b y2^2 + e with no confounder, b = 0.10 and 0.02, n = 20,000 |

## The datasets

- `freeiv_sim1.dta`, `freeiv_sim2.dta`: n = 5,000, seed 20260913, generated
  in Python (numpy) from the equations printed in `help freeiv`: U
  standardized chi2(3), gamma = 0.40, alpha2 = 0.80, scale consistency
  holding; V2 normal in sim1, drawn from U's skewed law in sim2.
- `freeiv_proxy.dta`: simulated with two indicators of one confounder.
- `freeiv_card.dta`: the Card (1995) extract of Wooldridge's teaching
  datasets (`card`, through the Python package `wooldridge`), eleven
  variables.
