* freeiv_note.do -- every output and every number of the technical note of
* freeiv 1.0.0 (paper/freeiv_sj.tex).  The blocks of the note are excerpts of
* this log, commands included, at the note's line width.
*
* Run from the root of the repository, with freeiv installed (net install) or
* with src/ at hand:
*     do replication/freeiv_paper/freeiv_note.do
* The four datasets of the package are read from examples/; nlsw88 comes with
* Stata.  The two parts that need an internet connection come last: nhanes2
* (webuse) and wage1 (Boston College, the Wooldridge datasets), the second with
* a check of the LSZ draws against trigmm (ssc install trigmm), skipped when
* trigmm is not installed.  About four minutes, most of them in that check.
* Writes freeiv_note.log in the current folder.

version 16
clear all
set more off
set linesize 80
cap mata: mata clear

cap confirm file "examples/freeiv_sim1.dta"
if (_rc) {
    di as err "run this file from the root of the repository (examples/ not found)"
    exit 601
}
* the code of the repository when it is there, the installed package otherwise
cap confirm file "src/freeiv.ado"
if (_rc == 0) {
    foreach f in freeiv freeiv_engine freeiv_mata freeivmenu freeivdiag ///
                 freeivtest freeivreport {
        run "src/`f'.ado"
    }
}
cap log close fivnote
log using "freeiv_note.log", replace text name(fivnote)
di as txt "Stata " c(stata_version) ", " c(current_date)
* the datasets are read as a reader reads them after -net get freeiv-
cd examples

* ======================================================================
* Section 3.6 -- the bootstrap: the qme near the edge (z of D 1.49)
* ======================================================================
use freeiv_sim1, clear
set seed 20261006
bootstrap, reps(500) nodots: freeiv y1 x (y2)

* ======================================================================
* Section 4.1 -- where the model holds (truth 0.40)
* ======================================================================
use freeiv_sim1, clear
freeivmenu y1 x (y2)
freeiv y1 x (y2), method(all)

* ======================================================================
* Section 4.2 -- the edge of identification (truth 0.40, B close to 2A)
* ======================================================================
use freeiv_sim2, clear
freeiv y1 x (y2), method(all)

* ======================================================================
* Section 4.3 -- Card, one endogenous regressor: the interval is the result
* ======================================================================
use freeiv_card, clear
freeivmenu lwage exper expersq black south smsa (educ)
freeiv lwage exper expersq black south smsa (educ), method(all)
* the same model with fewer controls, and what its qme implies
freeiv lwage exper expersq (educ)
freeivdiag

* ======================================================================
* Section 4.4 -- an estimate from elsewhere: Card's IV against the interval
* ======================================================================
quietly freeiv lwage exper expersq black south smsa (educ)
freeivtest, gamma(0.132) segamma(0.049)

* ======================================================================
* Section 4.5 -- a second indicator of the confounder
* ======================================================================
freeivmenu lwage exper expersq black south smsa (educ motheduc)
freeiv lwage exper expersq black south smsa (educ motheduc)
freeivtest
quietly freeiv lwage exper expersq black south smsa (educ fatheduc)
display "with fatheduc: g2 " %5.3f e(g2) ", a1 - (g2 a2 + g3 a3) " %5.3f e(sc_d)
* two indicators the guards refuse: freeiv ends with r(498), e() posted
capture noisily freeiv lwage exper expersq black south smsa (educ IQ)
display "return code " _rc ", guard " e(guard)
capture noisily freeiv lwage exper expersq black south smsa (educ KWW)
display "return code " _rc ", guard " e(guard)

* ======================================================================
* Section 4.6 -- a model the data refute
* ======================================================================
sysuse nlsw88, clear
generate lwage = ln(wage)
freeiv lwage grade age i.race (tenure)
quietly freeiv lwage grade age i.race (tenure), method(gmm)
display "joint GMM: J " %4.1f e(J_gmm) ", p " %6.4f e(p_gmm)

* ======================================================================
* Section 4.8 -- Oster's range, and factor variables
* ======================================================================
use freeiv_card, clear
quietly freeiv lwage exper expersq black south smsa (educ), method(oster) rmax(0.5)
display "oster, rmax(0.5): " %5.3f e(oster) "   (OLS " %5.3f e(ols) ")"
quietly freeiv lwage exper expersq black south smsa (educ), method(oster) ///
    delta(2) rmax(0.8)
display "oster, delta(2) rmax(0.8): " %5.3f e(oster)
quietly freeiv lwage exper expersq black south smsa (educ), method(all)
local q1 = e(qme)
quietly freeiv lwage c.exper##c.exper i.black i.south i.smsa (educ), method(all)
display "factor variables: qme " %8.6f e(qme) " against " %8.6f `q1' ///
    " by hand (|difference| " %7.1e abs(e(qme) - `q1') ")"

* ======================================================================
* Section 5 -- what the third order cannot tell apart: no confounder,
* y1 = y2 + b y2^2 + e; the population value of the qme is 2/3 for every
* b != 0.  The same draws of y2 and e for both values of b.
* ======================================================================
foreach b in 0.10 0.02 {
    clear
    set seed 20260915
    quietly set obs 20000
    generate double y2 = rnormal()
    generate double y1 = y2 + `b'*y2^2 + rnormal()
    quietly freeiv y1 (y2), method(gmm)
    display "b = `b', n = " e(N) ": bounds [" %5.3f e(lo) ", " %5.3f e(hi) ///
        "], qme " %5.3f cond(e(at_vertex) == 1, e(vertex), e(qme)) ///
        cond(e(at_vertex) == 1, " (vertex)", "") ", z of D " %5.2f e(disc_z)
    display "    J of pgmm " %5.1f e(J_pgmm) " (p " %6.4f e(p_pgmm) ///
        "), J of gmm " %5.1f e(J_gmm) " (p " %6.4f e(p_gmm) ")"
}

* ======================================================================
* The parts that need an internet connection
* ======================================================================

* Section 4.7 -- survey data
webuse nhanes2, clear
freeiv bpsystol age female black (bmi), vce(svy)
display "the design: s.e. of the OLS slope " %5.3f e(se_gt)
quietly freeiv bpsystol age female black (bmi) [pweight=finalwgt]
display "pweights alone: s.e. of the OLS slope " %5.3f e(se_gt)
quietly freeiv bpsystol age female black (bmi), method(gmm) vce(svy)
display "joint GMM under the design: J " %4.1f e(J_gmm)

* Section 3.6 -- a route that often has no estimate: LSZ on wage1
use "http://fmwww.bc.edu/ec-p/data/wooldridge/wage1", clear
quietly freeiv lwage exper (educ), method(lsz)
display "wage1, full sample: lsz " %6.4f e(lsz) ", s.e. " %6.4f e(se_lsz)
set seed 20260913
bootstrap, reps(50) nodots: freeiv lwage exper (educ), method(lsz)
display "the draws kept: mean " %6.4f el(e(b_bs), 1, 1) ", s.e. " %6.4f _se[educ]
* The same 50 draws, redrawn with -bsample- (the RNG state is restored around
* trigmm, so the draws stay aligned with the prefix's): where freeiv's solver
* finds no root, does trigmm, given 1,000 iterations, find one (criterion
* below 1e-9), or not?  About two minutes.
cap which trigmm
if (_rc == 0) {
    tempfile wage1 res
    quietly save `wage1'
    tempname P
    postfile `P' byte conv double tQ int trc using `res'
    set seed 20260913
    forvalues r = 1/50 {
        quietly use `wage1', clear
        bsample
        cap quietly freeiv lwage exper (educ), method(lsz)
        local conv = e(lsz_conv)
        local tQ = .
        local trc = .
        if (`conv' != 1) {
            local rng = c(rngstate)
            cap quietly trigmm lwage educ, covariates(exper) p(0 1) iterate(1000)
            local trc = _rc
            if (_rc == 0) local tQ = e(Q)
            set rngstate `rng'
        }
        post `P' (`conv') (`tQ') (`trc')
    }
    postclose `P'
    use `res', clear
    quietly count if conv == 1
    display "freeiv finds a root in " r(N) " of the 50 draws"
    quietly count if conv != 1 & tQ < 1e-9
    display "  trigmm finds one in " r(N) " of the others"
    quietly count if conv != 1 & tQ >= 1e-9 & tQ < 1e-5
    display "  trigmm stops near zero (1e-9 to 1e-5) in " r(N)
    quietly count if conv != 1 & (tQ >= 1e-5 | trc != 0)
    display "  trigmm finds none either in " r(N)
}
else display "(trigmm is not installed: ssc install trigmm)"

cd ..
log close fivnote
