*! freeivreport 0.2.0  14sep2026  A. Araar (Universite Laval / PEP)
*! One command for the whole pipeline: what the data can carry, every route
*! against the identified set, and the tests -- assembled into a single table
*! that can be written straight to a file for a paper.
*!
*!   freeivreport depvar [indepvars] (endogvar) [if] [in] [fw pw aw]
*!                [, SAVing(filename[, replace]) noMENU noTESTs
*!                   Level(#) QUANtile(#) DELta(#) RMAX(#) SIGN(#) ]
*!
*!   freeivreport depvar [indepvars] (endogvar1 endogvar2) ...
*!
*! The report is deliberately ordered the way the results should be read:
*! identification first, estimates second, tests last.  An estimate is never
*! shown without the interval that judges it, and never without the signal
*! that says whether its own route is available at all.
*!
*! saving() writes the same table as a comma-delimited file with one row per
*! line: block, label, value, standard error, note.

cap program drop freeivreport
cap program drop _fivr_row
cap program drop _fivr_open
cap program drop _fivr_close

program define freeivreport, rclass
    version 16
    syntax anything(name=eqs equalok) [if] [in] [fw pw aw] ///
        [, SAVing(string) noMENU noTESTs Level(cilevel) ///
           QUANtile(real 0.25) DELta(real 1) RMAX(real -1) SIGN(real 1) ]

    * ---- how many endogenous variables ------------------------------------
    local p1 = strpos("`eqs'", "(")
    local p2 = strpos("`eqs'", ")")
    if (`p1' == 0 | `p2' == 0 | `p2' < `p1') {
        di as err "the endogenous variable(s) must be given in parentheses"
        exit 198
    }
    local endog = trim(substr("`eqs'", `p1' + 1, `p2' - `p1' - 1))
    local nend : word count `endog'
    local rest  = trim(substr("`eqs'", 1, `p1' - 1)) + " " ///
                + trim(substr("`eqs'", `p2' + 1, .))
    gettoken depvar exog : rest
    local exog = trim("`exog'")

    local rmopt = ""
    if (`rmax' >= 0) local rmopt "rmax(`rmax')"
    local wt ""
    if ("`weight'" != "") local wt "[`weight'`exp']"

    _fivr_open, saving(`saving')
    local fh "`r(fh)'"

    * ======================= identification ================================
    if ("`menu'" == "" & `nend' == 1) {
        cap qui freeivmenu `eqs' `if' `in' `wt', quantile(`quantile')
        if (_rc == 0) {
            local m_zm03 = r(z_m03)
            local m_skew = r(skew2)
            local m_Flew = r(F_lewbel)
            local m_plew = r(p_lewbel)
            local m_dz   = r(disc_z)
            local m_k    = r(k)
            local m_ks   = r(kstar)
            local m_B    = r(B)
            local m_mu   = r(mu)
        }
    }

    * ======================= estimation =====================================
    * with two endogenous regressors freeiv now refuses the options that
    * steer a one-endogenous route, so they must not be forwarded
    if (`nend' == 2) {
        qui freeiv `eqs' `if' `in' `wt', level(`level')
    }
    else {
        qui freeiv `eqs' `if' `in' `wt', method(all) level(`level') ///
            quantile(`quantile') delta(`delta') `rmopt' sign(`sign')
    }

    di
    di as txt "freeiv report" _col(30) as res "`depvar'" as txt " on " ///
       as res "`endog'" as txt cond("`exog'" != "", ", controls " + "`exog'", "") ///
       _col(66) "n = " as res %8.0f e(n)
    di as txt "{hline 76}"

    if (`nend' == 1) {
        di as txt "IDENTIFICATION" _col(38) "value" _col(52) "verdict"
        _fivr_row "identification" "bounds, width" ///
            `=e(hi) - e(lo)' . "[`=string(e(lo),"%7.5f")', `=string(e(hi),"%7.5f")']" "`fh'"
        if ("`menu'" == "") {
            _fivr_row "identification" "z of m03 (third-order route)" ///
                `m_zm03' . "`=cond(abs(`m_zm03')>3,"strong",cond(abs(`m_zm03')>2,"weak","absent"))'" "`fh'"
            _fivr_row "identification" "z of the discriminant D" ///
                `m_dz' . "`=cond(abs(`m_dz')>2,"usable","fragile")'" "`fh'"
            _fivr_row "identification" "F of eps2^2 on X (lewbel12)" ///
                `m_Flew' . "`=cond(`m_plew'<0.05,"present","absent")'" "`fh'"
            _fivr_row "identification" "k* (rre sets 1)" ///
                `m_ks' . "`=cond(abs(`m_ks'-1)<0.5,"plausible","doubtful")'" "`fh'"
            _fivr_row "identification" "B = E[V2^3] (hme sets 0)" ///
                `m_B' . "" "`fh'"
        }

        di as txt "{hline 76}"
        di as txt "ESTIMATES" _col(32) "estimate" _col(46) "s.e." _col(60) "vs the set"
        foreach k in ols qme sce rre hme qbe gmm pgmm lsz lewbel12 copula ///
                     rank rpiv ape oster {
            local kk "`k'"
            if ("`k'" == "rank") local kk "g_rank"
            * since 0.2.0 gmm is the JOINT nine-moment route of Araar (2026c)
            * and pgmm the profiled one; they are different estimators
            if ("`k'" == "gmm")  local kk "g_jgmm"
            if ("`k'" == "pgmm") local kk "g_gmm"
            local v  = e(`kk')
            local se = .
            if ("`k'" == "ols") local se = e(se_gt)
            if ("`k'" == "qme") local se = e(se_qme)
            if ("`k'" == "sce") local se = e(se_sce)
            if ("`k'" == "rre") local se = e(se_rre)
            if ("`k'" == "hme") local se = e(se_hme)
            * suppressed, not absent: on the boundary there is nothing for the
            * delta method to expand around
            if ("`k'" == "gmm")  local se = cond(e(j_bound), ., e(se_jgmm))
            if ("`k'" == "pgmm") local se = e(se_gmm)
            * a non-converged solver does not get to show a standard error
            if ("`k'" == "lsz") local se = cond(e(lsz_conv) != 1, ., e(se_lsz))
            local nt ""
            if (`v' < . & e(lo) < .) {
                * the OLS slope IS the upper bound, so compare with a
                * relative tolerance or it is reported as outside its own set
                local tol = 1e-9 * max(1, abs(e(hi)))
                if (`v' < e(lo) - `tol' | `v' > e(hi) + `tol') local nt "outside"
                else if (abs(`v' - e(hi)) < `tol')             local nt "at the upper bound"
                else if (abs(`v' - e(lo)) < `tol')             local nt "at the lower bound"
                else                                          local nt "inside"
            }
            if ("`k'" == "lsz" & e(lsz_conv) != 1) local nt "`nt', DID NOT CONVERGE"
            _fivr_row "estimate" "`k'" `v' `se' "`nt'" "`fh'"
        }
        _fivr_row "estimate" "identified set" . . ///
            "[`=string(e(lo),"%7.5f")', `=string(e(hi),"%7.5f")']" "`fh'"

        if ("`tests'" == "") {
            di as txt "{hline 76}"
            di as txt "TESTS" _col(32) "statistic" _col(46) "p" _col(60) "tests"
            cap qui freeivtest, level(`level')
            if (_rc == 0) {
                _fivr_row "test" "endogeneity, chi2(2)" `=r(W_endo)' ///
                    `=r(p_endo)' "theta = 0" "`fh'"
                _fivr_row "test" "qme - ols, z" `=r(z_ols)' . "theta = 0" "`fh'"
                _fivr_row "test" "qme - rre, z" `=r(z_rre)' . "k* = 1" "`fh'"
                _fivr_row "test" "qme - sce, z" `=r(z_sce)' . "k = 1" "`fh'"
                _fivr_row "test" "qme - hme, z" `=r(z_hme)' . "B = 0" "`fh'"
            }
            _fivr_row "test" "GMM 2-3-4 joint, J chi2(1)" `=e(J9)' `=e(p_J9)' ///
                "scale consistency + linearity" "`fh'"
            * the region is valid where the standard error is not, so it is
            * reported whenever the minimum sits on a boundary -- and its
            * share of the identified interval says what the higher moments
            * actually added
            if (e(ar_frac) < .) {
                _fivr_row "test" "  J within 3.84 of its min" `=e(ar_lo)' ///
                    `=e(ar_hi)' "`=string(100*e(ar_frac),"%4.0f")'% of the set" "`fh'"
                if (e(ar_frac) > 0.95) {
                    di as txt "      the region IS the identified set: orders 3 and 4"
                    di as txt "      add nothing to the assumption-free bounds here"
                }
            }
            if (e(j_bound) == 1) {
                di as txt "      the joint minimum is on the boundary, so its s.e."
                di as txt "      and the chi2(1) law for its J are both invalid"
            }
            _fivr_row "test" "GMM 2-3-4 profiled, J chi2(1)" `=e(J)' `=e(p_J)' ///
                "h1..h5 forced to zero" "`fh'"
            _fivr_row "test" "  its identification factor" `=e(z_idfac)' . ///
                "`=cond(abs(e(z_idfac))<2,"near the knife edge","usable")'" "`fh'"
        }
    }

    * ======================= model B ========================================
    else {
        di as txt "MODEL B: TWO INDICATORS" _col(32) "estimate" _col(46) "s.e."
        _fivr_row "estimate" "g2 (`=word("`endog'",1)')" `=e(g2)' `=e(se_g2)' "" "`fh'"
        _fivr_row "estimate" "g3 (`=word("`endog'",2)')" `=e(g3)' `=e(se_g3)' "" "`fh'"
        _fivr_row "estimate" "a1 (free loading of U)" `=e(a1)' `=e(se_a1)' "" "`fh'"
        _fivr_row "estimate" "g2 a2 + g3 a3" `=e(sc)' . "model A would call this a1" "`fh'"
        _fivr_row "estimate" "condition number" `=e(cnum)' . ///
            "`=cond(e(cnum)>100,"ill conditioned","")'" "`fh'"
        _fivr_row "estimate" "`=e(endog3)' as an instrument" `=e(ivgap)' . ///
            "compare with g2 above" "`fh'"

        di as txt "{hline 76}"
        di as txt "TESTS" _col(32) "statistic" _col(46) "s.e. / p"
        _fivr_row "test" "scale consistency, a1 - sc" ///
            `=e(sc_d)' `=e(se_sc)' "z = `=string(e(z_sc),"%6.3f")'" "`fh'"
        _fivr_row "test" "one factor, R1 - R3" `=e(of_d)' `=e(se_of)' ///
            "z = `=string(e(z_of),"%6.3f")'; needs symmetric V" "`fh'"
        _fivr_row "test" "GMM16, J chi2(4)" `=e(q_J)' `=e(q_pJ)' ///
            "the one-factor linear structure" "`fh'"
        _fivr_row "test" "guard of Proposition 1" `=e(guard)' . ///
            "`=cond(e(guard)==0,"none fired","a guard fired")'" "`fh'"
    }

    di as txt "{hline 76}"
    di as txt "Read down, not across: the first block says which routes these"
    di as txt "data can carry, the second shows what each returns against the"
    di as txt "identified set, the third what the data say about the assumptions."
    _fivr_close, fh("`fh'") saving(`saving')

    return scalar n = e(n)
    if (`nend' == 1) {
        return scalar lo = e(lo)
        return scalar hi = e(hi)
        return scalar qme = e(qme)
    }
    return local endog "`endog'"
    return local depvar "`depvar'"
end


program define _fivr_row
    args block label value se note fh
    local sv "         ."
    local ss "         ."
    cap if (`value' < .) local sv = string(`value', "%10.6f")
    cap if (`se' < .)    local ss = string(`se', "%10.6f")
    if ("`block'" == "identification") {
        di as txt "  " %-34s "`label'" _col(38) as res "`sv'" ///
           _col(52) as txt "`note'"
    }
    else {
        di as txt "  " %-28s "`label'" _col(32) as res "`sv'" ///
           _col(46) "`ss'" _col(60) as txt "`note'"
    }
    if ("`fh'" != "") {
        * plain CSV, no quoting: commas inside a field become semicolons, so
        * no embedded double quote can reach the file and break a later read
        local lab = subinstr(`"`label'"', ",", ";", .)
        local nt2 = subinstr(`"`note'"', ",", ";", .)
        local lab = subinstr(`"`lab'"', `"""', "", .)
        local nt2 = subinstr(`"`nt2'"', `"""', "", .)
        file write `fh' `"`block',`lab',`=trim("`sv'")',`=trim("`ss'")',`nt2'"' _n
    }
end


program define _fivr_open, rclass
    syntax [, SAVing(string)]
    if ("`saving'" == "") {
        return local fh ""
        exit
    }
    gettoken fn rest : saving, parse(",")
    local fn = trim(`"`fn'"')
    local fn = subinstr(`"`fn'"', `"""', "", .)
    if (strpos(lower("`rest'"), "replace")) local rep "replace"
    * NOT a tempname: a tempname handle is closed when the program that
    * created it ends, so the handle must survive across _fivr_row calls
    cap file close __fivrpt
    file open __fivrpt using `"`fn'"', write text `rep'
    file write __fivrpt "block,label,value,se,note" _n
    return local fh "__fivrpt"
end


program define _fivr_close
    syntax [, fh(string) SAVing(string)]
    if ("`fh'" == "") exit
    file close `fh'
    gettoken fn rest : saving, parse(",")
    local fn = trim(subinstr(`"`fn'"', `"""', "", .))
    di as txt "table written to " as res `"`fn'"'
end
