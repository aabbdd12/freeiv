*! freeivtest 0.1.0  13sep2026  A. Araar (Universite Laval / PEP)
*! A direct test of endogeneity, agreement tests between the instrument-free
*! routes, and the confrontation of an outside estimate with the identified set.
*!
*!   freeivtest [, GAMma(#) SEGamma(#) Level(#) ]
*!
*! FIRST BLOCK -- theta = 0, with no estimate of gamma.  Under theta = 0 the
*! OLS residual r = xi - gt eps2 is INDEPENDENT of eps2, not merely
*! uncorrelated, so E[r eps2^2] = m12 - gt m03 and E[r^2 eps2] = m21 - 2gt m12
*! + gt^2 m03 both vanish; the Wald statistic on the pair is chi2 with 2 df.
*! In simulation it holds its size (0.060 at n = 526, 0.047 at n = 5000 for a
*! nominal 5%) and it keeps power where the qme - ols comparison loses its own,
*! because it never estimates gamma and so escapes the 1/sqrt(D) fragility.
*!
*! Detecting endogeneity and measuring it need DIFFERENT signals.  Measuring
*! gamma needs A = alpha2^3 E[U^3] != 0, a skewed confounder.  Detecting
*! theta > 0 needs only m03 != 0, from U or from V2 indifferently.  A symmetric
*! confounder beside a skewed V2 is therefore detectable and not estimable
*! (power 0.34 at n = 526 in that configuration, against 0.24 with a skewed
*! confounder).  When U and V2 are BOTH symmetric the third-order block is
*! empty, the test has no power at all (0.03), and that case is visible before
*! any estimation as z(m03) ~ 0 in freeivmenu -- read z(m03) first, exactly as
*! one reads the first-stage F before an endogeneity test in IV.
*!
*! The test maintains the linear one-factor specification, so a rejection is
*! evidence of theta > 0 OR of a departure from that structure.
*!
*! Every closed form is a smooth function of the SAME moment vector, so one
*! covariance matrix V serves them all and the difference of any two has
*! variance (g_a - g_b)' V (g_a - g_b) / n.  Each pair is therefore a test of
*! the assumption that separates the two routes:
*!
*!   qme - ols   theta = 0        no confounding at all
*!   qme - rre   k* = 1           scale-free equal-variance restriction
*!   qme - sce   k  = 1           the same, in the units of the data
*!   qme - hme   B  = 0           symmetry of V2
*!
*! The QME is the reference because it is the route with the weakest
*! assumptions: scale consistency plus E[U^3] != 0, and nothing else.
*!
*! NOTE ON WHAT IS NOT TESTED HERE.  None of these pairs tests scale
*! consistency itself.  The map gamma -> (theta, sigma2_V1, sigma2_V2) is
*! DERIVED from alpha1 = gamma1 alpha2, so every statistic above presupposes
*! it.  At second order the restriction is in fact untestable: the interval
*! [gamma-tilde/2, gamma-tilde] is exactly the positivity region, so an
*! admissible gamma always exists.  It IS tested at order four, by the
*! over-identifying J of method(gmm): the unrestricted model's Jacobian has
*! determinant -(gamma - tau)^5 (B kurt_U - A kurt_V2), so the restriction
*! binds except where that factor vanishes -- and a normal V2 sits exactly
*! there.  Read e(idfac) before reading the J.  It is also refuted, without
*! any fourth moment, with a SECOND INDICATOR, where model B identifies a1
*! freely and it can be compared with g2 a2 + g3 a3, and by confronting the
*! identified set with a route that does not impose it -- LSZ, Lewbel (2012),
*! or a genuine instrument.  That last test is what gamma(#) performs.

cap program drop freeivtest
cap program drop _fivt_pair

program define freeivtest, rclass
    version 16
    syntax [, GAMma(real -99999) SEGamma(real -1) Level(cilevel) ]

    if ("`e(cmd)'" != "freeiv") {
        di as err "freeivtest is used after freeiv"
        exit 301
    }
    cap confirm matrix e(grad)
    if (_rc) {
        di as err "e(grad) not found: re-run freeiv with version 0.4.0 or later"
        exit 498
    }

    di
    di as txt "freeiv agreement tests" _col(52) "Number of obs = " ///
       as res %9.0f e(n)

    * ---- 1. endogeneity itself, with no estimate of gamma ------------------
    tempname VV GG gg SS WW
    matrix `VV' = e(Vmom)
    local gt  = e(gt)
    local m02 = e(m02)
    local m11 = e(m11)
    local m03 = e(m03)
    local m12 = e(m12)
    local m21 = e(m21)
    matrix `GG' = J(2, 10, 0)
    matrix `GG'[1, 1] = `m03' * `m11' / `m02'^2
    matrix `GG'[1, 2] = -`m03' / `m02'
    matrix `GG'[1, 4] = -`gt'
    matrix `GG'[1, 5] = 1
    local cc = -2 * `m12' + 2 * `gt' * `m03'
    matrix `GG'[2, 1] = `cc' * (-`m11' / `m02'^2)
    matrix `GG'[2, 2] = `cc' / `m02'
    matrix `GG'[2, 4] = `gt'^2
    matrix `GG'[2, 5] = -2 * `gt'
    matrix `GG'[2, 6] = 1
    local g1 = `m12' - `gt' * `m03'
    local g2 = `m21' - 2 * `gt' * `m12' + `gt'^2 * `m03'
    matrix `gg' = (`g1' \ `g2')
    local W = .
    local pW = .
    cap {
        matrix `SS' = `GG' * `VV' * `GG'' / e(n)
        matrix `WW' = `gg'' * syminv(`SS') * `gg'
        local W  = `WW'[1, 1]
        local pW = chi2tail(2, `W')
    }
    di as txt "{hline 76}"
    di as txt "Endogeneity itself: theta = 0, with no estimate of gamma"
    di as txt "    E[r eps2^2] = m12 - gt m03" _col(46) as res %12.6f `g1'
    di as txt "    E[r^2 eps2] = m21 - 2gt m12 + gt^2 m03" _col(46) as res %12.6f `g2'
    di as txt "    Wald, chi2(2)" _col(46) as res %12.4f `W' ///
       as txt "   P>chi2 " as res %6.4f `pW'
    di as txt "    r is the OLS residual of xi on eps2.  Under theta = 0 it is"
    di as txt "    independent of eps2, not merely uncorrelated, so both cross"
    di as txt "    moments vanish.  This never estimates gamma, so it keeps its"
    di as txt "    power where the qme loses its own (small D)."
    if (abs(e(skew2)) < 0.1) {
        di as res "    caution: eps2 is nearly symmetric, so the third-order block"
        di as res "    is close to empty and this test has little or no power"
    }
    di as txt "    Read it with z(m03) from -freeivmenu-: detecting theta > 0 needs"
    di as txt "    m03 != 0 from either U or V2, while ESTIMATING gamma needs"
    di as txt "    A = alpha2^3 E[U^3] != 0, a skewed confounder.  The two signals"
    di as txt "    are different, so detection can succeed where estimation fails."

    di as txt "{hline 76}"
    di as txt "Pairwise agreement between the closed forms (reference: qme)"
    di as txt "  " %-12s "pair" _col(18) %11s "difference" _col(31) %10s "s.e." ///
       _col(43) %8s "z" _col(53) %7s "P>|z|" "   tests"
    di as txt "{hline 76}"

    _fivt_pair qme ols "theta = 0  (no confounding)"
    local z_ols = r(z)
    _fivt_pair qme rre "k* = 1     (scale free)"
    local z_rre = r(z)
    _fivt_pair qme sce "k  = 1     (unit dependent)"
    local z_sce = r(z)
    _fivt_pair qme hme "B  = 0     (symmetry of V2)"
    local z_hme = r(z)

    di as txt "{hline 76}"
    di as txt "A large |z| says the two routes disagree by more than sampling noise,"
    di as txt "so the assumption that separates them is rejected.  A small |z| is"
    di as txt "mutual corroboration: two routes with different assumptions agree."
    di as txt "None of these tests scale consistency -- they all presuppose it."

    * ---- an outside estimate against the identified set --------------------
    local out_z = .
    local out_d = .
    if (`gamma' != -99999) {
        local lo = e(lo)
        local hi = e(hi)
        local selo = e(se_lo)
        local sehi = e(se_gt)
        di as txt "{hline 76}"
        di as txt "Outside estimate " as res %10.6f `gamma' as txt ///
           " against the identified set"
        di as txt "    bounds [" as res %8.6f `lo' as txt ", " as res %8.6f `hi' ///
           as txt "]   s.e. " as res %8.6f `selo' as txt " and " ///
           as res %8.6f `sehi'
        if (`gamma' >= `lo' & `gamma' <= `hi') {
            di as txt "    inside the set: these data do not refute the"
            di as txt "    scale-consistent one-factor model through this route"
        }
        else {
            if (`gamma' < `lo') {
                local out_d = `gamma' - `lo'
                local sb = `selo'
                local which "lower"
            }
            else {
                local out_d = `gamma' - `hi'
                local sb = `sehi'
                local which "upper"
            }
            local sv = `sb'^2
            if (`segamma' >= 0) local sv = `sv' + `segamma'^2
            local out_z = cond(`sv' > 0, `out_d' / sqrt(`sv'), .)
            di as txt "    distance to the `which' bound" _col(46) ///
               as res %12.6f `out_d'
            di as txt "    z" _col(46) as res %12.4f `out_z' ///
               as txt "   P>|z| " as res %6.4f 2 * normal(-abs(`out_z'))
            di as res "    => outside the set: under scale consistency no value of"
            di as res "    gamma can produce this estimate without a negative"
            di as res "    variance.  Run -freeivdiag, gamma(`gamma')- to see which."
            if (`segamma' < 0) {
                di as txt "    (the outside estimate is treated as fixed; supply"
                di as txt "    segamma(#) to add its own sampling error)"
            }
            else {
                di as txt "    (the two standard errors are added as if independent,"
                di as txt "    which they are not -- they come from the same sample)"
            }
        }
    }
    di as txt "{hline 76}"

    return scalar W_endo = `W'
    return scalar p_endo = `pW'
    return scalar g1     = `g1'
    return scalar g2     = `g2'
    return scalar z_ols = `z_ols'
    return scalar z_rre = `z_rre'
    return scalar z_sce = `z_sce'
    return scalar z_hme = `z_hme'
    return scalar out_z = `out_z'
    return scalar out_d = `out_d'
end


program define _fivt_pair, rclass
    args A B what
    local ga = e(`A')
    local gb = e(`B')
    local df = .
    local sd = .
    local z  = .
    local p  = .
    if (`ga' < . & `gb' < .) {
        local df = `ga' - `gb'
        tempname G V d q
        matrix `G' = e(grad)
        matrix `V' = e(Vmom)
        cap matrix `d' = `G'["`A'", 1...] - `G'["`B'", 1...]
        if (_rc == 0) {
            cap matrix `q' = `d' * `V' * `d''
            if (_rc == 0) {
                local vd = `q'[1, 1] / e(n)
                if (`vd' > 0 & `vd' < .) {
                    local sd = sqrt(`vd')
                    local z  = `df' / `sd'
                    local p  = 2 * normal(-abs(`z'))
                }
            }
        }
    }
    di as txt "  " %-12s "`A' - `B'" _col(18) as res %11.6f `df' ///
       _col(31) %10.6f `sd' _col(43) %8.2f `z' _col(53) %7.4f `p' ///
       as txt "   `what'"
    return scalar diff = `df'
    return scalar se   = `sd'
    return scalar z    = `z'
    return scalar p    = `p'
end
