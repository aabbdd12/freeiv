*! freeivmenu 0.3.1  15sep2026  A. Araar (Universite Laval / PEP)
*! The identification menu: before any estimation, what these data can carry.
*! One line per family of strategies, with the signal it needs, its value, and
*! a verdict.
*!
*!   freeivmenu depvar [indepvars] (endogvar) [if] [in] [fw pw aw]
*!              [, QUANtile(#) BW(#) ]
*!
*! 0.3.0 adds the local slope profile of E[xi | eps2]: the kernel-weighted
*! slope of xi on eps2 at seven percentiles of eps2.  Flat means Gaussian-
*! like (the third order will find little); its minimum is an upper bound
*! on gamma when E[U | eps2] is increasing; the ratio of the two tail slopes
*! tends to 1 + alpha1/(gamma alpha2), i.e. 2 under scale consistency.
*! Reference: python/freeiv/profile.py, values in python/export_profile.txt.
*!
*! Every displayed statistic is returned in r().

cap program drop freeivmenu
cap program drop _fivm_verdict

program define freeivmenu, rclass
    version 16
    syntax anything(name=eqs equalok) [if] [in] [fw pw aw] ///
        [, QUANtile(real 0.25) BW(real -1) ]

    local p1 = strpos("`eqs'", "(")
    local p2 = strpos("`eqs'", ")")
    if (`p1' == 0 | `p2' == 0 | `p2' < `p1') {
        di as err "the endogenous variable must be given in parentheses"
        exit 198
    }
    local endog = trim(substr("`eqs'", `p1' + 1, `p2' - `p1' - 1))
    local rest  = trim(substr("`eqs'", 1, `p1' - 1)) + " " ///
                + trim(substr("`eqs'", `p2' + 1, .))
    local rest  = trim(stritrim("`rest'"))
    gettoken depvar exog : rest
    local nend : word count `endog'
    if (`nend' != 1) {
        di as err "freeivmenu expects exactly one endogenous variable"
        exit 198
    }
    confirm numeric variable `depvar' `endog'

    * factor variables among the controls, as in freeiv itself
    local exogfv   "`exog'"
    local exogbase ""
    if ("`exog'" != "") {
        cap fvrevar `exog', list
        if (_rc) {
            di as err "invalid varlist for the exogenous controls:"
            di as err "    `exog'"
            exit 198
        }
        local exogbase "`r(varlist)'"
        confirm numeric variable `exogbase'
    }
    local clash : list depvar & endog
    if ("`clash'" != "") {
        di as err "`clash' cannot be both the dependent and the endogenous variable"
        exit 198
    }

    freeiv_engine
    marksample touse
    markout `touse' `depvar' `endog' `exogbase'
    qui count if `touse'
    local N = r(N)

    tempvar wv xi e2 e2sq
    if ("`weight'" != "") {
        qui gen double `wv' `exp' if `touse'
        local wname "`wv'"
        local wexp2 "[aw = `wv']"
    }
    else {
        local wname ""
        local wexp2 ""
    }
    if ("`exogfv'" != "") {
        fvrevar `exogfv' if `touse', substitute
        local exog "`r(varlist)'"
    }
    qui gen double `xi' = .
    qui gen double `e2' = .
    mata: _freeiv_all("`depvar'", "`endog'", "`exog'", "`wname'", "`touse'", ///
                      `quantile', "`xi'", "`e2'")
    tempname M
    matrix `M' = __freeiv_M
    matrix drop __freeiv_M
    cap matrix drop __freeiv_V __freeiv_G
    local vals "n sum_w m02 m11 m20 m03 m12 m21 m30 m04 m13 m22 m02c m11c m20c gt skew2 lo hi disc disc_se disc_z vertex rstar root1 root2 nroots qme sce rre hme ols qbe qbe_R qbe_q1 qbe_skew qbe_exk qbe_nsub qbe_clip qbe_dom theta sV2 sV1 k kstar A B mu se_gt se_lo se_qme se_sce se_rre se_hme se_vertex"
    local j = 0
    foreach v of local vals {
        local ++j
        local `v' = `M'[1, `j']
    }

    * ---- heteroskedasticity signal for Lewbel (2012) ----------------------
    local F_lew = .
    local p_lew = .
    local df_lew = 0
    if ("`exog'" != "") {
        qui gen double `e2sq' = `e2'^2 if `touse'
        cap qui regress `e2sq' `exog' `wexp2' if `touse'
        if (_rc == 0) {
            local F_lew  = e(F)
            local df_lew = e(df_m)
            local p_lew  = Ftail(e(df_m), e(df_r), e(F))
        }
    }

    * ---- z statistic of the third moment of eps2 --------------------------
    tempvar c3
    qui summarize `e2' `wexp2' if `touse', meanonly
    qui gen double `c3' = (`e2' - r(mean))^3 if `touse'
    qui summarize `c3' `wexp2' if `touse'
    local se_m03 = r(sd) / sqrt(`N')
    local z_m03  = cond(`se_m03' > 0, `m03' / `se_m03', .)

    * ---- local slope profile of E[xi | eps2] ------------------------------
    * Gaussian kernel, h = 2 * 1.06 * sd(eps2) * N^(-1/5) unless bw() is given;
    * -regress [aw=kw], vce(robust)- is exactly the kernel-weighted local
    * linear slope with the HC1 sandwich the Python reference computes.
    qui summarize `e2' `wexp2' if `touse'
    local h = cond(`bw' > 0, `bw', 2 * 1.06 * r(sd) * `N'^(-0.2))
    local probs "5 10 25 50 75 90 95"
    _pctile `e2' `wexp2' if `touse', p(`probs')
    forvalues j = 1/7 {
        local q`j' = r(r`j')
    }
    tempname P
    matrix `P' = J(7, 3, .)
    tempvar kw
    forvalues j = 1/7 {
        qui gen double `kw' = exp(-0.5 * ((`e2' - `q`j'') / `h')^2) if `touse'
        if ("`wname'" != "") qui replace `kw' = `kw' * `wname' if `touse'
        cap qui regress `xi' `e2' [aw = `kw'] if `touse', vce(robust)
        if (_rc == 0) {
            matrix `P'[`j', 1] = `q`j''
            matrix `P'[`j', 2] = _b[`e2']
            matrix `P'[`j', 3] = _se[`e2']
        }
        drop `kw'
    }
    matrix colnames `P' = eps2 slope se
    matrix rownames `P' = p5 p10 p25 p50 p75 p90 p95
    local pmin = .
    local jmin = 0
    local jmax = 0
    local pmax = .
    forvalues j = 1/7 {
        local b = `P'[`j', 2]
        if (`b' < . & (`pmin' >= . | `b' < `pmin')) {
            local pmin = `b'
            local jmin = `j'
        }
        if (`b' < . & (`pmax' >= . | `b' > `pmax')) {
            local pmax = `b'
            local jmax = `j'
        }
    }
    local pmin_se = cond(`jmin' > 0, `P'[`jmin', 3], .)
    local pmax_se = cond(`jmax' > 0, `P'[`jmax', 3], .)
    local pratio  = cond(`P'[1, 2] != 0 & `P'[1, 2] < ., `P'[7, 2] / `P'[1, 2], .)
    local pz      = cond(`pmin_se' < . & `pmax_se' < ., ///
                    (`pmax' - `pmin') / sqrt(`pmax_se'^2 + `pmin_se'^2), .)
    local pnames "5 10 25 50 75 90 95"
    local pmin_p : word `jmin' of `pnames'

    * ================= display ============================================
    di
    di as txt "Identification menu: what these data can carry"
    di as txt "`depvar' on `endog'" cond("`exogfv'" != "", ", controls `exogfv'", "") ///
       _col(56) "n = " as res %8.0f `N'
    di as txt "{hline 76}"
    di as txt %-22s "route" %-28s "signal" %10s "value" "  verdict"
    di as txt "{hline 76}"

    * 1. bounds ------------------------------------------------------------
    local w = `hi' - `lo'
    di as txt %-22s "bounds" %-28s "width of the interval" ///
       as res %10.4f `w' as txt "  always"
    di as txt %-22s "" %-28s "[gamma-tilde/2, gamma-tilde]" ///
       as res "  [" %6.4f `lo' ", " %6.4f `hi' "]"

    * 2. third order -------------------------------------------------------
    _fivm_verdict `z_m03' 3 2
    local v1 "`s'"
    di as txt %-22s "third order: qme," %-28s "skewness of eps2" ///
       as res %10.4f `skew2'
    di as txt %-22s "  hme, lsz" %-28s "z of m03" as res %10.2f `z_m03' ///
       as txt "  `v1'"
    if (`disc' < .) {
        _fivm_verdict `disc_z' 2 1
        local v2 "`s'"
        di as txt %-22s "" %-28s "z of the discriminant D" ///
           as res %10.2f `disc_z' as txt "  `v2'"
    }
    if (`disc' < 0) {
        di as txt %-22s "" %-28s "D < 0: the roots merge" ///
           as res %10s "." as txt "  vertex used"
    }
    if (`se_qme' < .) {
        di as txt %-22s "" %-28s "implied s.e. of the QME" ///
           as res %10.4f `se_qme' as txt cond(`se_qme' > `w' / 4, ///
           "  wide", "  usable")
    }
    if (`A' < .) {
        di as txt %-22s "" %-28s "A = alpha2^3 E[U^3]" as res %10.4f `A'
        di as txt %-22s "" %-28s "mu = A/(A+B)" as res %10.4f `mu' ///
           as txt cond(abs(`mu' - 1/3) < 0.10, "  near 1/3", "")
    }

    * 3. equal variances, scale free (RRE) ----------------------------------
    if (`kstar' < .) {
        di as txt %-22s "equal variances: rre" %-28s "k* = sV1/(g^2 sV2), sets 1" ///
           as res %10.4f `kstar' as txt cond(abs(`kstar' - 1) < 0.5, "  plausible", "  doubtful")
        di as txt %-22s "" %-28s "rre = m20/(2 m11)" as res %10.4f `rre' ///
           as txt cond(`rre' < `lo' | `rre' > `hi', "  outside", "  inside")
    }
    else di as txt %-22s "equal variances: rre" %-28s "k* not identified here" ///
         as res %10s "." as txt "  undecidable"

    * 4. equal variances, unit dependent (SCE) ------------------------------
    if (`k' < .) {
        di as txt %-22s "  and its sce variant" %-28s "k = sV1/sV2, sce sets 1" ///
           as res %10.4f `k' as txt cond(abs(`k' - 1) < 0.5, "  plausible", "  doubtful")
        di as txt %-22s "" %-28s "k depends on units; k* not"
    }

    * 5. symmetry of V2 (HME) ------------------------------------------------
    if (`B' < .) {
        local rel = cond(`m03' != 0, abs(`B' / `m03'), .)
        di as txt %-22s "symmetry of V2: hme" %-28s "B = E[V2^3], hme sets 0" ///
           as res %10.4f `B' as txt cond(`rel' < 0.2, "  plausible", "  doubtful")
    }

    * 6. heteroskedasticity (Lewbel 2012) ------------------------------------
    if (`F_lew' < .) {
        di as txt %-22s "heteroskedasticity:" %-28s "F of eps2^2 on X" ///
           as res %10.2f `F_lew' as txt cond(`p_lew' < 0.05, "  present", "  absent")
        di as txt %-22s "  lewbel12" %-28s "p" as res %10.4f `p_lew'
    }
    else di as txt %-22s "heteroskedasticity:" %-28s "no exogenous control" ///
         as res %10s "." as txt "  unavailable"

    * 7. local slope profile --------------------------------------------------
    di as txt %-22s "local slope profile" %-28s "slope of xi on eps2 at p5" ///
       as res %10.4f `P'[1, 2] as txt "  h = " %5.3f `h'
    di as txt %-22s "  of E[xi | eps2]" %-28s "p10 / p25 / p50" ///
       as res %7.4f `P'[2, 2] " " %7.4f `P'[3, 2] " " %7.4f `P'[4, 2]
    di as txt %-22s "" %-28s "p75 / p90 / p95" ///
       as res %7.4f `P'[5, 2] " " %7.4f `P'[6, 2] " " %7.4f `P'[7, 2]
    di as txt %-22s "" %-28s "curvature z (max-min)" as res %10.2f `pz' ///
       as txt cond(`pz' >= 2 & `pz' < ., "  curved", cond(`pz' < ., "  flat", ""))
    * the bound reading needs a rising profile (m increasing): the minimum
    * must sit at p5 or p10.  An interior minimum is reported as such.
    local tight = cond(`jmin' > 2, "interior min", ///
        cond(`pmin' + 1.96 * `pmin_se' < `hi', "tightens", "no tightening"))
    di as txt %-22s "" %-28s "min slope, at p`pmin_p'" as res %10.4f `pmin' ///
       as txt "  `tight'"
    di as txt %-22s "" %-28s "tail ratio p95/p5" as res %10.2f `pratio' ///
       as txt cond(`pratio' < ., cond(abs(`pratio' - 2) < 0.5, "  near 2", "  not 2"), "")

    * 8. proxy ---------------------------------------------------------------
    di as txt %-22s "two indicators" %-28s "freeiv y1 x (y2 y3)" ///
       as res %10s "." as txt "  by syntax"

    di as txt "{hline 76}"
    di as txt "Reading: the first line is always true; the others say whether the"
    di as txt "assumption that identifies is carried by the data.  A weak signal does"
    di as txt "not make an estimator wrong, it makes it imprecise."
    di as txt "B = 0 and k* = 1 are restrictions on the model; k = 1 is a restriction"
    di as txt "on the units, which is why the sce is reported below the rre."
    di as txt "The slope profile is flat when the confounder is Gaussian-like; its"
    di as txt "minimum bounds gamma from above when E[U|eps2] is increasing, and the"
    di as txt "tail ratio tends to 2 under scale consistency.  A curved profile is"
    di as txt "a skewed confounder or a non-linear outcome equation: it does not say"
    di as txt "which."

    * ---- r() ---------------------------------------------------------------
    foreach v of local vals {
        return scalar `v' = ``v''
    }
    return scalar z_m03 = `z_m03'
    return scalar prof_h     = `h'
    return scalar prof_min   = `pmin'
    return scalar prof_minse = `pmin_se'
    return scalar prof_ratio = `pratio'
    return scalar prof_z     = `pz'
    return matrix profile = `P'
    return scalar F_lewbel = `F_lew'
    return scalar p_lewbel = `p_lew'
    return matrix moments = `M'
end


program define _fivm_verdict
    args z hi lo
    if ("`z'" == "" | `z' >= .) local s "undecidable"
    else if (abs(`z') >= `hi')  local s "strong"
    else if (abs(`z') >= `lo')  local s "weak"
    else                        local s "absent"
    c_local s "`s'"
end
