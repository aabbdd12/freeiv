*! freeivdiag 0.2.0  13sep2026  A. Araar (Universite Laval / PEP)
*! Post-estimation diagnostic: the nuisance parameters implied by a value of
*! gamma, and what they say about the model.
*!
*!   freeivdiag [, GAMma(#) ]
*!
*! The second-order system pins theta, sigma2_V1 and sigma2_V2 as functions of
*! gamma ALONE.  The diagnostic is therefore method free: hand it the estimate
*! from the QME, the SCE, the RRE, GMM, Lewbel (2012) or LSZ and it will say
*! what each of them implies.  gamma(#) evaluates an outside value.
*!
*! The partial identification bounds ARE the positivity constraints:
*!   theta >= 0      is equivalent to   gamma <= gamma-tilde
*!   sigma2_V2 >= 0  is equivalent to   gamma >= gamma-tilde / 2
*! An estimate outside the bounds therefore implies a negative variance, and
*! the command says so.

cap program drop freeivdiag

program define freeivdiag, rclass
    version 16
    syntax [, GAMma(real -99999) ]

    if ("`e(cmd)'" != "freeiv") {
        di as err "freeivdiag is used after freeiv"
        exit 301
    }

    tempname M
    matrix `M' = e(moments)
    local vals "n sum_w m02 m11 m20 m03 m12 m21 m30 m04 m13 m22 m02c m11c m20c gt skew2 lo hi disc disc_se disc_z vertex rstar root1 root2 nroots qme sce rre hme ols qbe qbe_R qbe_q1 qbe_skew qbe_exk qbe_nsub qbe_clip qbe_dom theta sV2 sV1 k kstar A B mu se_gt se_lo se_qme se_sce se_rre se_hme se_vertex"
    local j = 0
    foreach v of local vals {
        local ++j
        local `v' = `M'[1, `j']
    }

    if (`gamma' == -99999) {
        local g = e(gamma)
        local src "`e(method)'"
    }
    else {
        local g = `gamma'
        local src "supplied value"
    }
    if (`g' >= . | `g' == 0) {
        di as err "no value of gamma to diagnose"
        exit 198
    }

    * ---- nuisances implied by g -------------------------------------------
    local th  = `m02c' * (`gt' - `g') / `g'
    local s2  = `m02c' - `th'
    local s1  = `m20c' - `g'^2 * `m02c' - 3 * `g'^2 * `th'

    * At the exact bounds theta or sigma2_V2 is zero up to machine rounding;
    * without a tolerance we would print -0.000000 NEGATIVE.
    local tol = 1e-8 * max(1, abs(`m02c'), abs(`m20c'))
    if (abs(`th') < `tol') local th = 0
    if (abs(`s2') < `tol') local s2 = 0
    if (abs(`s1') < `tol') local s1 = 0

    local kk  = cond(`s2' != 0, `s1' / `s2', .)
    local kks = cond(`kk' < ., `kk' / `g'^2, .)
    local Ag  = `m12' / `g' - `m03'
    local Bg  = 2 * `m03' - `m12' / `g'
    local mug = cond(`m03' != 0, `Ag' / `m03', .)

    * ---- identified kurtosis of the confounder (ARAARP3) ------------------
    local kU = .
    if (`th' > 0) {
        local P = (`m13' - `g' * `m04') / `g'
        local Q = (`m22' - `g'^2 * `m04' - `s1' * `m02') / `g'^2
        local A4 = (3 * `Q' - 7 * `P') / 2
        local kU = `A4' / `th'^2
    }

    * ================= display =============================================
    di
    di as txt "freeiv diagnostic" _col(40) "gamma = " as res %9.6f `g' ///
       as txt "  (`src')"
    di as txt "{hline 72}"
    di as txt "Position inside the identification interval"
    di as txt "    bounds [" as res %8.6f `lo' as txt ", " as res %8.6f `hi' as txt "]"
    if (`g' < `lo' | `g' > `hi') {
        di as res "    OUTSIDE the bounds: under scale consistency this value"
        di as res "    implies a negative variance (see below)"
    }
    else {
        local pos = cond(`hi' > `lo', 100 * (`g' - `lo') / (`hi' - `lo'), .)
        di as txt "    relative position" _col(46) as res %9.1f `pos' as txt " %"
    }

    di as txt "{hline 72}"
    di as txt "Implied nuisances (second-order system, functions of gamma alone)"
    di as txt "    theta = alpha2^2 Var(U)" _col(46) as res %12.6f `th' ///
       as txt cond(`th' < 0, "   NEGATIVE", "")
    di as txt "    sigma2_V2" _col(46) as res %12.6f `s2' ///
       as txt cond(`s2' < 0, "   NEGATIVE", "")
    di as txt "    sigma2_V1" _col(46) as res %12.6f `s1' ///
       as txt cond(`s1' < 0, "   NEGATIVE", "")
    di as txt "    k = sigma2_V1 / sigma2_V2" _col(46) as res %12.6f `kk' ///
       as txt cond(`kk' >= ., "   not identified", ///
                   cond(abs(`kk' - 1) < 0.5, "   near 1", "   far from 1"))
    di as txt "    k* = k / gamma^2" _col(46) as res %12.6f `kks' ///
       as txt cond(`kks' >= ., "   not identified", ///
                   cond(abs(`kks' - 1) < 0.5, "   near 1", "   far from 1"))
    local sh = cond(`m02c' > 0, `th' / `m02c', .)
    di as txt "    share of the confounder in Var(eps2)" _col(46) as res %12.4f `sh'

    di as txt "{hline 72}"
    di as txt "Third moment"
    di as txt "    A = alpha2^3 E[U^3]" _col(46) as res %12.6f `Ag'
    di as txt "    B = E[V2^3]  (the hme sets it to 0)" _col(46) as res %12.6f `Bg'
    di as txt "    mu = A/(A+B)" _col(46) as res %12.6f `mug'
    if (`mug' < . & abs(`mug' - 1/3) < 0.10) {
        di as res "    mu near 1/3: this is where the two roots of the QME merge"
        di as res "    and where its standard error explodes"
    }
    if (`kU' < .) {
        di as txt "    identified kurtosis of the confounder" _col(46) as res %12.4f `kU' ///
           as txt cond(`kU' >= 3 & `kU' <= 10, "   plausible", "   outside 3-10")
    }

    di as txt "{hline 72}"
    local nw = 0
    if (`th' < 0 | `s2' < 0 | `s1' < 0) {
        di as res "warning: at least one implied variance is negative -- a model with"
        di as res "a single linear confounder cannot produce this value of gamma on"
        di as res "these data"
        local ++nw
    }
    if (`kks' < . & abs(`kks' - 1) > 0.5) {
        di as res "warning: k* departs from 1, so the scale-free equal-variance"
        di as res "restriction behind the RRE is doubtful here"
        local ++nw
    }
    if (`kk' < . & abs(`kk' - 1) > 0.5) {
        di as res "warning: k departs from 1, so the equal-variance assumption of the"
        di as res "SCE is doubtful -- and note that k is not invariant to the units of"
        di as res "`e(depvar)' and `e(endog)', which is why k* is the one to read"
        local ++nw
    }
    if (`Bg' < . & `m03' != 0 & abs(`Bg' / `m03') > 0.2) {
        di as res "warning: B departs from 0, so the symmetry assumption of V2 on which"
        di as res "the HME rests is doubtful"
        local ++nw
    }
    if (`nw' == 0) di as txt "no warning"

    return scalar gamma = `g'
    return scalar theta = `th'
    return scalar sV2   = `s2'
    return scalar sV1   = `s1'
    return scalar k     = `kk'
    return scalar kstar = `kks'
    return scalar A     = `Ag'
    return scalar B     = `Bg'
    return scalar mu    = `mug'
    return scalar kurtU = `kU'
    return scalar lo    = `lo'
    return scalar hi    = `hi'
    return scalar inbounds = (`g' >= `lo' & `g' <= `hi')
    return local  source "`src'"
end
