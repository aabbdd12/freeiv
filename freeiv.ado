*! freeiv 0.8.2  15sep2026  A. Araar (Universite Laval / PEP)
*! Instrument-free estimation: linear structural model with one endogenous
*! regressor, identified through higher-order moments under the scale
*! consistency restriction alpha1 = gamma1 * alpha2.
*!
*!   freeiv depvar [indepvars] (endogvar) [if] [in] [fw pw aw]
*!          [, METHod(name) Level(#) QUANtile(#) noHEADer ]
*!
*! METHod : bounds | ols | qme | sce | rre | hme | qbe          (this family)
*!          gmm | pgmm | lsz                          (over-identified)
*!          lewbel12 | copula | rank | rpiv | ape | oster       (literature)
*!          all                                          (default: qme)
*!
*! 0.8.2 gates the one-factor ratio verdict of model B on the z of the two
*! third-order cross-moments (both >= 2), the same z as freeivmenu prints,
*! stored in e(z_m223), e(z_m233); R2 and R3 divide by those moments.
*! 0.8.0 splits the GMM in two, because they are two estimators and only one
*! of them is the paper's.
*!
*! method(gmm) is the JOINT GMM of Araar (2026c): nine moment conditions in
*! eight free parameters, one over-identifying restriction, a J with 1 df.
*! It is minimised deterministically rather than by multi-start, and the
*! reason is not tidiness.  Fix gamma AND theta and the only non-linearity
*! left, the product theta*sV2, disappears, so all nine residuals become
*! linear in the six remaining parameters; the eight-parameter problem is a
*! two-dimensional surface, and a surface can be gridded.  No seed, no
*! starting value.  Multi-start over the eight parameters returns a LOCAL
*! minimum on six of the eight datasets of that paper.
*!
*! Because the criterion is often nearly flat, the route also returns
*! e(ar_lo), e(ar_hi) and e(ar_frac): the values of gamma at which the
*! profiled J is within 3.84 of its minimum, and their share of the
*! identified interval.  That region inverts the J at each fixed gamma rather
*! than inverting a Wald statistic around the optimum, so it stays valid on
*! the boundary and on a flat criterion, where the standard error does not.
*! When it fills the whole interval the moments of orders three and four have
*! added nothing to the assumption-free bounds.  e(j_bound) marks a minimum
*! on the boundary, where no standard error is printed; e(j_weak) marks
*! theta < 0.05.
*!
*! method(pgmm) is the PROFILED variant, which is what method(gmm) meant
*! before 0.8.0: given gamma, five restrictions determine theta, sV2, sV1, A
*! and B exactly, leaving four residuals in three unknowns, hence again one
*! over-identifying restriction and a J with 1 df.  It weights the first five
*! conditions infinitely rather than optimally, so it is a different
*! estimator, not the same one solved differently.  It is kept for two
*! reasons: it is one-dimensional, hence gridded exhaustively and never
*! troubled by local minima, and its fifth residual is the QME quadratic
*! itself, which shows the third-order route as the part the fourth moments
*! reweight.  It therefore still returns a value when the discriminant is
*! negative and the QME has no real root.  A wide gap between e(g_jgmm) and
*! e(g_gmm) is a symptom of a flat criterion, not an error in either.
*!
*! WHAT THE J TESTS.  Scale consistency itself, alpha1 = gamma1 alpha2, plus
*! the linearity of the confounder's effect.  Write the unrestricted model as
*! eps2 = U + V2 and xi = tau U + gamma V2 + V1 with tau free; scale
*! consistency is tau = 2 gamma.  The Jacobian of the nine moments with
*! respect to the nine free parameters has determinant
*!
*!     -(gamma - tau)^5 * (B kurt_U - A kurt_V2)
*!
*! with kurt_U = A4 - 3 theta^2 and kurt_V2 = B4 - 3 sV2^2.  The unrestricted
*! model is therefore locally identified, and tau = 2 gamma is a genuine
*! testable restriction, everywhere except on two surfaces: tau = gamma, and
*! the vanishing of the second factor.  A NORMAL V2 sits exactly on the second
*! one, since it makes B = 0 and kurt_V2 = 0, and there the J has no power at
*! all.  Simulation confirms the algebra: against a violation by a factor two
*! at n = 4000, the rejection rate at 5% is 0.065 with V2 normal, 0.435 with
*! V2 symmetric but leptokurtic, and 1.000 with V2 skewed.
*!
*! The J also rejects a U entering non-linearly (power 1.00 at n = 8000) and a
*! loading that varies with X (0.38).  Read e(idfac) and its z before reading
*! the J: near zero, the J is uninformative about scale consistency.  The
*! factor's magnitude is not a calibrated measure of power -- a determinant is
*! not scale free -- only its vanishing is meaningful.
*!
*! method(lsz) is Lewbel, Schennach & Zhang (2024) with p(0 1), exactly as
*! the published command trigmm implements it: the coefficients of X are
*! parameters of the system rather than partialled out, and the start is
*! trigmm's own.  With k exogenous variables the system has 5 + 2(k+1)
*! parameters and as many moments, so it is JUST IDENTIFIED -- the estimate is
*! the root of the moment vector and does not depend on the optimiser, which
*! is why a damped Newton reproduces trigmm exactly.  sign(-1) flips the
*! orientation as trigmm's own option does.  p(0 1 2) over-identifies and its
*! criterion has several local minima; the ado does not attempt it, use trigmm
*! itself there.
*!
*! LSZ leaves beta, the loading of the confounder in the outcome equation,
*! entirely free.  That is the whole difference with scale consistency, and it
*! is why confronting e(lsz) with the identification bounds -- which the
*! command does -- is a test of the scale-consistent model rather than of LSZ.
*!
*! The literature routes are reported on the same data so that the
*! identification bounds judge them all together.  They are not peers: each
*! buys its point identification with an assumption of its own, and none has
*! an analytic standard error here -- use the bootstrap prefix.  oster()
*! takes delta(#) and rmax(#); the defaults are 1 and min(1.3 R1, 1).
*!
*! MODEL B -- two endogenous regressors sharing one latent confounder:
*!
*!   freeiv depvar [indepvars] (endogvar1 endogvar2) [if] [in] [fw pw aw]
*!
*! Two variables inside one parenthesis block switch to the proxy model of
*! ARAARP4.  Y2 and Y3 are two indicators of the same U, and the closed form
*! of Theorem 1 identifies g2, g3 and the free loading a1 without any
*! scale-consistency restriction, with the three guards of Proposition 1.
*! method() does not apply there.
*!
*! Bricks 1-2: closed forms and analytic standard errors by the delta method
*! on the stacked moments; the influence function accounts for the estimation
*! of the coefficients on X.  The Q-BE is not smooth and has no analytic
*! standard error: use the bootstrap prefix.
*!
*! Reference: Araar, A. (2026), Zenodo 10.5281/zenodo.22067980,
*! 22068143, 22753299, 22764001.

cap program drop freeiv
cap program drop _freeiv_display
cap program drop _freeiv_fmt
cap program drop _freeiv_ident
cap program drop _freeiv_pdisplay
cap program drop _freeiv_pident
cap program drop _fivp_line

program define freeiv, eclass
    version 16
    if replay() {
        if ("`e(cmd)'" != "freeiv") error 301
        syntax [, noHEADer Level(cilevel) *]
        _freeiv_display, `header'
        exit
    }

    syntax anything(name=eqs equalok) [if] [in] [fw pw aw] ///
        [, METHod(name) Level(cilevel) QUANtile(real 0.25) noHEADer ///
           DELta(real 1) RMAX(real -1) SIGN(real 1) ]

    * ---- option values that can be judged before the data ----------------
    * An option quietly absorbed is worse than one refused: sign() used to go
    * straight to Mata, where any positive value acted as 1 and any negative
    * as -1, so a typo survived intact into e(lsz_sign).
    if (`sign' != 1 & `sign' != -1) {
        di as err "sign() must be 1 or -1"
        exit 198
    }
    if (`delta' >= .) {
        di as err "delta() must be a number"
        exit 198
    }
    if (`rmax' != -1 & (`rmax' <= 0 | `rmax' > 1)) {
        di as err "rmax() is an R-squared: it must lie in (0, 1]"
        di as err "    omit it for the default, min(1.3 R1, 1)"
        exit 198
    }

    * ---- method ----------------------------------------------------------
    local usermeth "`method'"
    if ("`method'" == "") local method qme
    local method = lower("`method'")
    if (!inlist("`method'", "bounds", "ols", "qme", "sce", "rre") & ///
        !inlist("`method'", "hme", "qbe", "all") & ///
        !inlist("`method'", "lewbel12", "copula", "rank", "rpiv") & ///
        !inlist("`method'", "ape", "oster", "gmm", "pgmm", "lsz")) {
        di as err "method() must be one of:"
        di as err "    bounds ols qme sce rre hme qbe        (this family)"
        di as err "    gmm pgmm lsz                          (over-identified)"
        di as err "    lewbel12 copula rank rpiv ape oster   (the literature)"
        di as err "    all"
        exit 198
    }
    local rmaxopt = cond(`rmax' < 0, ., `rmax')
    if (`quantile' <= 0 | `quantile' >= 1) {
        di as err "quantile() must lie strictly between 0 and 1"
        exit 198
    }

    * ---- endogenous variables go in parentheses --------------------------
    local p1 = strpos("`eqs'", "(")
    local p2 = strpos("`eqs'", ")")
    if (`p1' == 0 | `p2' == 0 | `p2' < `p1') {
        di as err "the endogenous variable must be given in parentheses:"
        di as err "    freeiv depvar [indepvars] (endogvar) ..."
        exit 198
    }
    local endog = trim(substr("`eqs'", `p1' + 1, `p2' - `p1' - 1))
    local rest  = trim(substr("`eqs'", 1, `p1' - 1)) + " " ///
                + trim(substr("`eqs'", `p2' + 1, .))
    local rest  = trim(stritrim("`rest'"))
    gettoken depvar exog : rest

    if ("`depvar'" == "") {
        di as err "dependent variable missing"
        exit 198
    }
    confirm numeric variable `depvar'
    local nend : word count `endog'
    if (`nend' == 0) {
        di as err "no endogenous variable inside the parentheses"
        exit 198
    }
    if (`nend' > 2) {
        di as err "freeiv accepts one endogenous regressor, or two sharing a"
        di as err "single latent confounder:  freeiv y1 x (y2 y3)"
        exit 198
    }
    confirm numeric variable `endog'

    * ---- factor variables among the controls -----------------------------
    * i.sex, ib2.region, i.sex##c.age.  fvrevar expands them into temporary
    * indicators with the base level omitted; the projection that partials X
    * out is computed with invsym, a generalized inverse, so an expansion
    * that turns out collinear is absorbed rather than fatal.
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

    * ---- one variable, one role ------------------------------------------
    local clash : list depvar & endog
    if ("`clash'" != "") {
        di as err "`clash' cannot be both the dependent and an endogenous variable"
        exit 198
    }
    if ("`exogbase'" != "") {
        local clash : list depvar & exogbase
        if ("`clash'" != "") {
            di as err "`clash' cannot be both the dependent variable and a control"
            exit 198
        }
        local clash : list endog & exogbase
        if ("`clash'" != "") {
            di as err "`clash' cannot be both endogenous and a control"
            exit 198
        }
    }
    if (`nend' == 2) {
        local e2chk : word 1 of `endog'
        local e3chk : word 2 of `endog'
        if ("`e2chk'" == "`e3chk'") {
            di as err "the two endogenous variables must differ"
            exit 198
        }
    }

    * ---- a route the specification cannot carry --------------------------
    if ("`exogfv'" == "" & inlist("`method'", "lewbel12", "oster")) {
        di as err "method(`method') needs at least one exogenous control"
        if ("`method'" == "lewbel12") {
            di as err "    its instruments are (X - mean X) * eps2"
        }
        else {
            di as err "    it compares the R-squared with and without controls"
        }
        exit 198
    }
    if ("`exogfv'" == "" & "`method'" == "all") {
        di as txt "note: no exogenous control, so lewbel12 and oster are unavailable"
    }
    if (`nend' == 2) {
        * These steer a route of the one-endogenous model and are inert
        * here.  Accepting them silently would return a result that looks
        * like the result of what was asked for, which is the failure mode
        * the whole validation block above exists to prevent -- so they are
        * refused, exactly like method(lewbel12) with no exogenous control.
        local inert ""
        if ("`usermeth'" != "")  local inert "`inert' method()"
        if (`quantile' != 0.25)  local inert "`inert' quantile()"
        if (`delta' != 1)        local inert "`inert' delta()"
        if (`rmax' != -1)        local inert "`inert' rmax()"
        if (`sign' != 1)         local inert "`inert' sign()"
        if ("`inert'" != "") {
            di as err "these options do not apply with two endogenous" ///
                      " regressors:`inert'"
            di as err "    the two-indicator model has a single closed-form"
            di as err "    route; level() and noheader are the only options"
            di as err "    it takes"
            exit 198
        }
    }

    freeiv_engine

    marksample touse
    markout `touse' `depvar' `endog' `exogbase'
    qui count if `touse'
    local N = r(N)
    if (`N' < 20) {
        di as err "too few observations (`N')"
        exit 2001
    }

    * ---- weights ---------------------------------------------------------
    tempvar wv
    if ("`weight'" != "") {
        qui gen double `wv' `exp' if `touse'
        local wname "`wv'"
    }
    else local wname ""

    * the factor variables become temporary indicators, built on the
    * estimation sample so that an empty level outside it cannot create a
    * column of zeros
    if ("`exogfv'" != "") {
        fvrevar `exogfv' if `touse', substitute
        local exog "`r(varlist)'"
    }

    * ================= model B: two indicators of one confounder ==========
    if (`nend' == 2) {
        gettoken en2 en3 : endog
        local en3 = trim("`en3'")
        mata: _freeiv_proxy("`depvar'", "`en2'", "`en3'", "`exog'", ///
                            "`wname'", "`touse'")
        tempname P PV b V
        matrix `P'  = __freeiv_P
        matrix `PV' = __freeiv_PV
        matrix drop __freeiv_P __freeiv_PV

        local pvals "n sum_w m22 m33 m23 mx2 mx3 m223 m233 mx23 g2 g3 a1 a2 a3 mu3 s2 s3 cnum sc guard se_g2 se_g3 se_a1 R1 R2 R3 disc_R ivgap t_rho"
        local j = 0
        foreach v of local pvals {
            local ++j
            local q`v' = `P'[1, `j']
        }

        * the over-identified GMM on sixteen moments (brick 7b)
        mata: _freeiv_gmm16("`depvar'", "`en2'", "`en3'", "`exog'", ///
                            "`wname'", "`touse'")
        tempname Q16
        matrix `Q16' = __freeiv_Q16
        matrix drop __freeiv_Q16
        local qvals "n q_J q_pJ q_df q_conv q_iter q_g2 q_g3 q_a1 q_a2 q_a3 q_mu3 q_s1 q_s2 q_s3 q_k1 q_k2 q_k3 q_se_g2 q_se_g3 q_se_a1 q_se_a2 q_se_a3 q_se_mu3 q_se_s1 q_se_s2 q_se_s3 q_se_k1 q_se_k2 q_se_k3"
        local j = 0
        foreach v of local qvals {
            local ++j
            local Q_`v' = `Q16'[1, `j']
        }

        * the two tests model A cannot perform (brick 7a)
        mata: _freeiv_ptests("`depvar'", "`en2'", "`en3'", "`exog'", ///
                             "`wname'", "`touse'")
        tempname PT
        matrix `PT' = __freeiv_PT
        matrix drop __freeiv_PT
        local ptvals "n t_m22 t_m33 t_m23 t_mx2 t_mx3 t_m223 t_m233 t_mx23 t_mx22 t_mx33 t_mxx2 t_mxx3 sc_d se_sc z_sc of_d se_of z_of z_m223 z_m233"
        local j = 0
        foreach v of local ptvals {
            local ++j
            local T_`v' = `PT'[1, `j']
        }

        if (`qguard' == 0) {
            matrix `b' = (`qg2', `qg3')
            matrix colnames `b' = `en2' `en3'
            matrix rownames `b' = y1
            matrix `V' = `PV'[1..2, 1..2]
            matrix colnames `V' = `en2' `en3'
            matrix rownames `V' = `en2' `en3'
            cap ereturn post `b' `V', esample(`touse') obs(`N') depname(`depvar')
            if (_rc) ereturn post `b', esample(`touse') obs(`N') depname(`depvar')
        }
        else ereturn post, esample(`touse') obs(`N') depname(`depvar')

        ereturn local cmd     "freeiv"
        ereturn local cmdline "freeiv `0'"
        ereturn local method  "proxy"
        ereturn local model   "B"
        ereturn local depvar  "`depvar'"
        ereturn local endog   "`endog'"
        ereturn local endog2  "`en2'"
        ereturn local endog3  "`en3'"
        ereturn local exog    "`exogfv'"
        ereturn local wtype   "`weight'"
        ereturn local wexp    "`exp'"
        ereturn local title   "Instrument-free estimation, two indicators"
        foreach v of local pvals {
            ereturn scalar `v' = `q`v''
        }
        ereturn scalar level = `level'
        ereturn matrix pmoments = `P'
        cap ereturn matrix pcov = `PV'
        foreach v of local ptvals {
            if ("`v'" != "n") ereturn scalar `v' = `T_`v''
        }
        foreach v of local qvals {
            if ("`v'" != "n") ereturn scalar `v' = `Q_`v''
        }
        ereturn matrix ptests = `PT'
        ereturn matrix gmm16  = `Q16'
        _freeiv_pdisplay, `header'
        exit
    }

    * ---- engine ----------------------------------------------------------
    tempvar xi e2
    qui gen double `xi' = .
    qui gen double `e2' = .
    mata: _freeiv_all("`depvar'", "`endog'", "`exog'", "`wname'", "`touse'", ///
                      `quantile', "`xi'", "`e2'")
    tempname M VM GM
    matrix `M' = __freeiv_M
    matrix `VM' = __freeiv_V
    matrix `GM' = __freeiv_G
    matrix drop __freeiv_M __freeiv_V __freeiv_G

    local vals "n sum_w m02 m11 m20 m03 m12 m21 m30 m04 m13 m22 m02c m11c m20c gt skew2 lo hi disc disc_se disc_z vertex rstar root1 root2 nroots qme sce rre hme ols qbe qbe_R qbe_q1 qbe_skew qbe_exk qbe_nsub qbe_clip qbe_dom theta sV2 sV1 k kstar A B mu se_gt se_lo se_qme se_sce se_rre se_hme se_vertex"
    local j = 0
    foreach v of local vals {
        local ++j
        local `v' = `M'[1, `j']
    }

    * ---- the literature routes, for the same data --------------------------
    mata: _freeiv_lit("`depvar'", "`endog'", "`exog'", "`wname'", "`touse'", ///
                      `sce', `delta', `rmaxopt')
    tempname L
    matrix `L' = __freeiv_L
    matrix drop __freeiv_L
    local lvals "n ols lewbel12 copula rank rpiv ape ape_route ape_theta oster ost_b0 ost_b1 ost_R0 ost_R1 ost_rmax ost_delta"
    local j = 0
    foreach v of local lvals {
        local ++j
        local L_`v' = `L'[1, `j']
    }

    * ---- the PROFILED GMM, method(pgmm) -----------------------------------
    mata: _freeiv_gmm("`depvar'", "`endog'", "`exog'", "`wname'", "`touse'")
    tempname G4
    matrix `G4' = __freeiv_G4
    matrix drop __freeiv_G4
    local gvals "n g_gmm se_gmm J p_J df_J g_theta g_sV2 g_sV1 g_A g_B A4 B4 idfac se_idfac z_idfac se_A4 se_B4 g_lo g_hi"
    local j = 0
    foreach v of local gvals {
        local ++j
        local G_`v' = `G4'[1, `j']
    }

    * ---- the JOINT nine-moment GMM: the estimator of Araar (2026c) ---------
    * Nine residuals in eight free parameters, minimised on the two-
    * dimensional profile in (gamma, theta) -- deterministic, no seed.  This
    * is what the paper reports; method(pgmm) above is the profiled variant,
    * which sets h1..h5 to zero exactly and is a different estimator.
    local doar = cond(inlist("`method'", "gmm", "all"), 1, 0)
    mata: _freeiv_jgmm("`depvar'", "`endog'", "`exog'", "`wname'", "`touse'", `doar')
    tempname G9
    matrix `G9' = __freeiv_G9
    matrix drop __freeiv_G9
    local jvals "n g_jgmm se_jgmm J9 p_J9 df_J9 j_theta j_sV2 j_sV1 j_A j_B j_A4 j_B4 j_kurt ar_lo ar_hi ar_frac j_weak j_bound g_lo g_hi"
    local j = 0
    foreach v of local jvals {
        local ++j
        local G_`v' = `G9'[1, `j']
    }

    * ---- LSZ (2024), as trigmm implements it, p(0 1) -----------------------
    * This one is not free.  It is a Levenberg-Marquardt on 5 + 2(k+1)
    * parameters with a numeric Jacobian, so each iteration costs 5 + 2(k+1)
    * passes over the sample and the total grows with the number of controls
    * as well as with n -- unlike every other route here, which is closed form
    * or a one-dimensional search.  It therefore runs only when it is asked
    * for.  The other methods still post e(lsz*), as missing, so that the
    * shape of e() does not depend on method().
    local zvals "n lsz se_lsz lsz_beta lsz_var_u lsz_var_v lsz_var_r lsz_crit lsz_conv lsz_iter lsz_nmom lsz_npar lsz_sign"
    tempname LZ
    if (inlist("`method'", "lsz", "all")) {
        mata: _freeiv_lsz("`depvar'", "`endog'", "`exog'", "`wname'", "`touse'", `sign')
        matrix `LZ' = __freeiv_LSZ
        matrix drop __freeiv_LSZ
        local j = 0
        foreach v of local zvals {
            local ++j
            local Z_`v' = `LZ'[1, `j']
        }
    }
    else {
        foreach v of local zvals {
            local Z_`v' = .
        }
        local Z_n = `N'
        matrix `LZ' = J(1, `: word count `zvals'', .)
        matrix `LZ'[1, 1] = `N'
        matrix colnames `LZ' = `zvals'
    }

    * ---- retained value --------------------------------------------------
    * method(all) retains the same point as method(qme), vertex fallback
    * included, so that e(b) and e(V) mean the same thing in both
    local qme_flag ""
    if (inlist("`method'", "qme", "all") & `disc' < 0) {
        local gamma = `vertex'
        local qme_flag "vertex"
    }
    else if ("`method'" == "qme")    local gamma = `qme'
    else if ("`method'" == "ols")    local gamma = `ols'
    else if ("`method'" == "sce")    local gamma = `sce'
    else if ("`method'" == "rre")    local gamma = `rre'
    else if ("`method'" == "hme")    local gamma = `hme'
    else if ("`method'" == "qbe")    local gamma = `qbe'
    else if ("`method'" == "bounds") local gamma = .
    else if ("`method'" == "gmm")    local gamma = `G_g_jgmm'
    else if ("`method'" == "pgmm")   local gamma = `G_g_gmm'
    else if ("`method'" == "lsz")    local gamma = `Z_lsz'
    else if (inlist("`method'", "lewbel12", "copula", "rank", "rpiv", "ape") ///
           | "`method'" == "oster")  local gamma = `L_`method''
    else                             local gamma = `qme'

    * ---- standard error of the retained value ----------------------------
    if ("`qme_flag'" == "vertex")    local segam = `se_vertex'
    else if (inlist("`method'", "qme", "all")) local segam = `se_qme'
    else if ("`method'" == "ols")    local segam = `se_gt'
    else if ("`method'" == "sce")    local segam = `se_sce'
    else if ("`method'" == "rre")    local segam = `se_rre'
    else if ("`method'" == "hme")    local segam = `se_hme'
    else if ("`method'" == "gmm")    local segam = cond(`G_j_bound', ., `G_se_jgmm')
    else if ("`method'" == "pgmm")   local segam = `G_se_gmm'
    * a non-converged solver reports the most precise-looking number in the
    * whole table; the standard error is the first thing that must not be
    * shown, and the note below says why
    else if ("`method'" == "lsz")    local segam = cond(`Z_lsz_conv' != 1, ., `Z_se_lsz')
    else                             local segam = .

    * ---- e() -------------------------------------------------------------
    tempname b V
    if (`gamma' < .) {
        matrix `b' = (`gamma')
        matrix colnames `b' = `endog'
        matrix rownames `b' = y1
        if (`segam' < .) {
            matrix `V' = (`segam'^2)
            matrix colnames `V' = `endog'
            matrix rownames `V' = `endog'
            ereturn post `b' `V', esample(`touse') obs(`N') depname(`depvar')
        }
        else ereturn post `b', esample(`touse') obs(`N') depname(`depvar')
    }
    else {
        ereturn post, esample(`touse') obs(`N') depname(`depvar')
    }
    ereturn local cmd     "freeiv"
    ereturn local cmdline "freeiv `0'"
    ereturn local method  "`method'"
    ereturn local depvar  "`depvar'"
    ereturn local endog   "`endog'"
    ereturn local exog    "`exogfv'"
    ereturn local wtype   "`weight'"
    ereturn local wexp    "`exp'"
    ereturn local title   "Instrument-free estimation"
    ereturn local qme_flag "`qme_flag'"
    foreach v of local vals {
        ereturn scalar `v' = ``v''
    }
    foreach v of local gvals {
        if ("`v'" != "n") ereturn scalar `v' = `G_`v''
    }
    foreach v of local zvals {
        if ("`v'" != "n") ereturn scalar `v' = `Z_`v''
    }
    * e(rank) belongs to -ereturn post- (the rank of V), so the rank-based
    * control-function estimate is posted as e(g_rank) instead
    foreach v of local lvals {
        if ("`v'" == "n")      continue
        if ("`v'" == "rank")   ereturn scalar g_rank = `L_rank'
        else                   ereturn scalar `v' = `L_`v''
    }
    foreach v of local jvals {
        if ("`v'" == "n")      continue
        if ("`v'" == "g_lo")   continue
        if ("`v'" == "g_hi")   continue
        ereturn scalar `v' = `G_`v''
    }
    ereturn scalar gamma    = `gamma'
    ereturn scalar se       = `segam'
    ereturn scalar quantile = `quantile'
    ereturn scalar level    = `level'
    ereturn matrix moments  = `M'
    ereturn matrix Vmom     = `VM'
    ereturn matrix grad     = `GM'
    ereturn matrix lit      = `L'
    ereturn matrix gmm      = `G4'
    ereturn matrix jgmm     = `G9'
    ereturn matrix lszmat   = `LZ'

    _freeiv_display, `header'
end


program define _freeiv_fmt
    args v
    if ("`v'" == "") local v = .
    cap confirm number `v'
    if (_rc) local v = .
    if (`v' >= .) local s "        ."
    else          local s = string(`v', "%10.6f")
    c_local s "`s'"
end


program define _freeiv_display
    syntax [, noHEADer]
    if ("`e(cmd)'" != "freeiv") exit
    if ("`e(model)'" == "B") {
        _freeiv_pdisplay, `header'
        exit
    }
    local lev = e(level)
    cap confirm number `lev'
    if (_rc | "`lev'" == "") local lev = 95
    if (`lev' >= . | `lev' <= 0 | `lev' >= 100) local lev = 95
    * Under a prefix (bootstrap, jackknife, svy) the prefix owns the
    * coefficient table: let it print, and only add the identification block.
    local pref "`e(prefix)'"
    if ("`pref'" != "" | e(N_reps) < .) {
        cap noisily ereturn display, level(`lev')
        cap noisily _freeiv_ident
        exit
    }
    local m "`e(method)'"
    local z = invnormal(1 - (100 - `lev') / 200)
    di
    if ("`header'" == "") {
        di as txt "Instrument-free estimation" _col(47) "Number of obs = " ///
           as res %9.0f e(n)
        di as txt "model: " as res "`e(depvar)'" as txt " on " ///
           as res "`e(endog)'" as txt cond("`e(exog)'" != "", ", controls " + "`e(exog)'", "")
        di as txt "method: " as res "`m'"
    }
    di as txt "{hline 72}"

    if ("`m'" == "all") {
        di as txt "Estimates of gamma" _col(26) "estimate" _col(40) "std. err."
        foreach k in ols qme sce rre hme qbe {
            _freeiv_fmt `=e(`k')'
            local sg "`s'"
            local sek = .
            if ("`k'" == "ols") local sek = e(se_gt)
            if ("`k'" == "qme") local sek = e(se_qme)
            if ("`k'" == "sce") local sek = e(se_sce)
            if ("`k'" == "rre") local sek = e(se_rre)
            if ("`k'" == "hme") local sek = e(se_hme)
            _freeiv_fmt `sek'
            di as txt "    " %-10s "`k'" _col(26) as res "`sg'" _col(40) "`s'"
        }
        _freeiv_fmt `=e(g_jgmm)'
        local sgj "`s'"
        _freeiv_fmt `=cond(e(j_bound), ., e(se_jgmm))'
        di as txt "    " %-10s "gmm" _col(26) as res "`sgj'" _col(40) "`s'" ///
           as txt "   J " as res %6.3f e(J9) as txt " p " as res %5.3f e(p_J9)
        _freeiv_fmt `=e(g_gmm)'
        local sgm "`s'"
        _freeiv_fmt `=e(se_gmm)'
        di as txt "    " %-10s "pgmm" _col(26) as res "`sgm'" _col(40) "`s'" ///
           as txt "   J " as res %6.3f e(J) as txt " p " as res %5.3f e(p_J)
        _freeiv_fmt `=e(lsz)'
        local slz "`s'"
        * the same rule as the single-method display: a solver that did not
        * converge does not get to show the tightest standard error in the
        * table
        _freeiv_fmt `=cond(e(lsz_conv) != 1, ., e(se_lsz))'
        di as txt "    " %-10s "lsz" _col(26) as res "`slz'" _col(40) "`s'" ///
           as txt cond(e(lsz_conv) != 1, "   did not converge", "")
        _freeiv_fmt `=e(lo)'
        local slo "`s'"
        _freeiv_fmt `=e(hi)'
        di as txt "    " %-10s "bounds" as res "[`=trim("`slo'")', `=trim("`s'")']"
        di as txt "  the literature, on the same data"
        foreach k in lewbel12 copula rank rpiv ape oster {
            local kk "`k'"
            if ("`k'" == "rank") local kk "g_rank"
            _freeiv_fmt `=e(`kk')'
            local sg "`s'"
            local mark ""
            if (e(`kk') < . & e(lo) < .) {
                if (e(`kk') < e(lo) | e(`kk') > e(hi)) local mark "  outside the bounds"
            }
            di as txt "    " %-10s "`k'" _col(26) as res "`sg'" ///
               _col(40) as txt "`mark'"
        }
        di as txt "    an estimate outside the bounds implies a negative variance"
        di as txt "    under scale consistency; -freeivdiag, gamma(#)- says which"
    }
    else if ("`m'" == "bounds") {
        _freeiv_fmt `=e(lo)'
        local slo "`s'"
        _freeiv_fmt `=e(hi)'
        di as txt "Partial identification: gamma in [" as res ///
           "`=trim("`slo'")'" as txt ", " as res "`=trim("`s'")'" as txt "]"
    }
    else {
        _freeiv_fmt `=e(gamma)'
        local sg "`s'"
        _freeiv_fmt `=e(se)'
        if (e(se) < .) {
            local cl = e(gamma) - `z' * e(se)
            local cu = e(gamma) + `z' * e(se)
            di as txt "    gamma (" as res "`e(endog)'" as txt ")" _col(24) ///
               as res "`sg'" as txt "  s.e. " as res "`=trim("`s'")'" ///
               as txt "   [" as res %8.6f `cl' as txt ", " as res %8.6f `cu' as txt "]"
            di as txt "        z = " as res %7.3f e(gamma)/e(se) ///
               as txt "    P>|z| = " as res %6.4f 2*normal(-abs(e(gamma)/e(se)))
        }
        else {
            di as txt "    gamma (" as res "`e(endog)'" as txt ")" _col(24) as res "`sg'"
            if ("`m'" == "qbe") {
                di as res "note: no analytic standard error -- the Q-BE is not smooth;"
                di as res "      use the -bootstrap- prefix"
            }
            else if ("`m'" == "lsz" & e(lsz_conv) != 1) {
                di as res "note: the LSZ solver DID NOT CONVERGE, so no standard error"
                di as res "      is shown and the point above should not be read as an"
                di as res "      estimate.  The objective is non-convex and needs"
                di as res "      multi-start; -trigmm- with several starts is the"
                di as res "      reference implementation for this route."
            }
            else if ("`m'" == "gmm" & e(j_bound)) {
                * not a missing standard error but a suppressed one: at a
                * boundary the delta method has nothing to expand around
                di as res "note: no standard error is REPORTED here, not none exists:"
                di as res "      the minimum is on the boundary of the parameter space,"
                di as res "      where the delta method has nothing to expand around and"
                di as res "      the chi2(1) law for J does not hold either.  Read the"
                di as res "      region below the J instead; it stays valid there."
            }
            else {
                di as res "note: no analytic standard error for method(`m') -- the"
                di as res "      literature routes are reported here for comparison, not"
                di as res "      re-derived; use the -bootstrap- prefix"
            }
        }
        if ("`m'" == "qbe" & e(qbe_dom) != 1) {
            di as res "note: Q-BE outside its domain of validity -- skewness of the"
            di as res "      first-stage residual " %6.3f e(qbe_skew) ", calibrated R " ///
               %5.3f e(qbe_R) cond(e(qbe_clip) != 0, " (clipped)", "") "."
            di as res "      The ARAARP3 calibration was fitted on laws of U with positive"
            di as res "      skewness; the figure above is not reliable"
        }
        if ("`m'" == "rre" & (e(rre) < e(lo) | e(rre) > e(hi))) {
            di as res "note: the RRE falls outside the identification bounds, so the"
            di as res "      restriction sV1 = gamma^2 sV2 is rejected by these data"
        }
        if ("`e(qme_flag)'" == "vertex") {
            di as res "note: negative discriminant -- the value shown is the vertex"
            di as res "      3*m12/(4*m03), where the two roots merge; it coincides with"
            di as res "      gamma when B = 2A, the very configuration where D vanishes"
        }
    }

    cap noisily _freeiv_ident
end


program define _freeiv_ident
    di as txt "{hline 72}"
    di as txt "Identification"
    _freeiv_fmt `=e(gt)'
    di as txt "    implied OLS slope (upper bound)" _col(46) as res "`s'"
    _freeiv_fmt `=e(lo)'
    local slo "`s'"
    _freeiv_fmt `=e(hi)'
    di as txt "    partial identification bounds" _col(46) as res ///
       "[`=trim("`slo'")', `=trim("`s'")']"
    _freeiv_fmt `=e(skew2)'
    di as txt "    skewness of the first-stage residual" _col(46) as res "`s'"
    _freeiv_fmt `=e(disc)'
    local sd "`s'"
    _freeiv_fmt `=e(disc_se)'
    di as txt "    discriminant D = 9m12^2 - 8m03m21" _col(46) as res "`sd'" ///
       as txt "  (s.e. " as res "`=trim("`s'")'" as txt ")"
    _freeiv_fmt `=e(disc_z)'
    di as txt "        D = gamma^2 (2A - B)^2 : z against 0" _col(46) as res "`s'"
    if (e(disc) >= 0 & e(disc) < .) {
        _freeiv_fmt `=e(root1)'
        local s1 "`s'"
        _freeiv_fmt `=e(root2)'
        di as txt "    roots of the quadratic" _col(46) as res ///
           "`=trim("`s1'")' and `=trim("`s'")'"
        di as txt "        of which inside the bounds" _col(46) as res %10.0f e(nroots)
        if (e(nroots) == 0) {
            di as res "note: NEITHER root lies in [gamma-tilde/2, gamma-tilde]."
            di as res "      That interval is exactly the region where the implied"
            di as res "      variances are non-negative, so the retained value is"
            di as res "      not compatible with the model under scale consistency."
            di as res "      Read the J below, and freeivmenu, before using it."
        }
    }
    else {
        _freeiv_fmt `=e(vertex)'
        di as txt "    no real root; vertex" _col(46) as res "`s'"
    }
    _freeiv_fmt `=e(rstar)'
    di as txt "    R* = m12^2/(m03 m21)" _col(46) as res "`s'"
    di as txt "        its floor is 8/9 when m03 m21 > 0"

    if (e(k) < .) {
        * These are the nuisances the RETAINED gamma implies, so they must be
        * recomputed from it: e(theta) and friends are solved at the QME root,
        * which is the retained value only under method(qme).  Given gamma,
        * the first three moment conditions give them in closed form.
        tempname MM gr th s2 s1 kk ks
        matrix `MM' = e(moments)
        * by NAME, not by position: e(moments) opens with n and sum_w, and
        * its layout is not something this display should depend on
        local c02 = colnumb(`MM', "m02")
        local c11 = colnumb(`MM', "m11")
        local c20 = colnumb(`MM', "m20")
        scalar `gr' = e(gamma)
        if (`c02' < . & `c11' < . & `c20' < . & `gr' < . & `gr' != 0) {
            scalar `th' = `MM'[1,`c11'] / `gr' - `MM'[1,`c02']
            scalar `s2' = `MM'[1,`c02'] - `th'
            scalar `s1' = `MM'[1,`c20'] - `gr'^2 * `MM'[1,`c02'] ///
                        - 3 * `gr'^2 * `th'
            scalar `kk' = `s1' / `s2'
            scalar `ks' = `s1' / (`gr'^2 * `s2')
        }
        else {
            scalar `th' = e(theta)
            scalar `s2' = e(sV2)
            scalar `s1' = e(sV1)
            scalar `kk' = e(k)
            scalar `ks' = e(kstar)
        }
        di as txt "{hline 72}"
        di as txt "Nuisances implied by the retained value"
        _freeiv_fmt `=scalar(`th')'
        di as txt "    theta = alpha2^2 Var(U)" _col(46) as res "`s'"
        _freeiv_fmt `=scalar(`s2')'
        di as txt "    sigma2_V2" _col(46) as res "`s'"
        _freeiv_fmt `=scalar(`s1')'
        di as txt "    sigma2_V1" _col(46) as res "`s'"
        * theta, sV2 and sV1 are variances up to a positive factor.  A
        * negative one is not a small number to be read as approximately
        * zero: it says the retained gamma is outside the region the model
        * allows.  freeivdiag has always said so; freeiv itself did not.
        if (`th' < 0 | `s2' < 0 | `s1' < 0) {
            di as res "note: an implied variance is NEGATIVE, which no model can"
            di as res "      produce.  The retained value lies outside the region"
            di as res "      where scale consistency is feasible; the numbers just"
            di as res "      above are arithmetic, not estimates."
        }
        _freeiv_fmt `=scalar(`kk')'
        di as txt "    k = sV1 / sV2" _col(46) as res "`s'" ///
           as txt "  (sce sets k = 1)"
        _freeiv_fmt `=scalar(`ks')'
        di as txt "    k* = sV1 / (gamma^2 sV2)" _col(46) as res "`s'" ///
           as txt "  (rre sets k* = 1)"
        if (abs(`ks' - 1) > 0.5 & `ks' < .) {
            di as res "note: k* is far from 1 -- the scale-free equal-variance"
            di as res "      restriction behind the RRE is doubtful on these data"
        }
        di as txt "    k is not invariant to the units of `e(depvar)' or `e(endog)';"
        di as txt "    k* is.  Read k* when judging an equal-variance restriction."
    }
    if (e(J) < .) {
        di as txt "{hline 72}"
        di as txt "Over-identification at order four"
        di as txt "    joint nine-moment GMM, method(gmm)"
        _freeiv_fmt `=e(g_jgmm)'
        di as txt "        gamma" _col(46) as res "`s'"
        di as txt "        J, 1 df" _col(46) as res %10.4f e(J9) ///
           as txt "   P>chi2 " as res %6.4f e(p_J9)
        * The region costs a profile of the criterion, so it is computed only
        * when the joint route is the one being reported.  Guard on it being
        * there at all -- and note that in Stata a missing value is LARGER
        * than any number, so an unguarded -> 0.95- test fires on a missing.
        if (e(ar_frac) < .) {
            _freeiv_fmt `=e(ar_lo)'
            local arl "`s'"
            _freeiv_fmt `=e(ar_hi)'
            di as txt "        J within 3.84 of its minimum on [" as res "`arl'" ///
               as txt ", " as res "`s'" as txt "]"
            di as txt "        that is " as res %5.1f 100*e(ar_frac) ///
               as txt "% of [gamma-tilde/2, gamma-tilde]"
        }
        else {
            di as txt "        (the region where J stays within 3.84 of its"
            di as txt "        minimum is computed under method(gmm) and"
            di as txt "        method(all) only)"
        }
        * The region above inverts the J at each fixed gamma instead of
        * inverting a Wald statistic around the optimum, so it stays valid
        * where the standard error does not: on the boundary, and where the
        * criterion is flat.  When it fills the whole interval, the fourth-
        * order moments have added nothing to Proposition 1.
        if (e(ar_frac) > 0.95 & e(ar_frac) < .) {
            di as res "    note: that region is the WHOLE identified interval, so the"
            di as res "          moments of orders 3 and 4 add nothing here to the"
            di as res "          assumption-free bounds.  The point estimate is the"
            di as res "          argmin of a flat criterion; report the interval."
        }
        if (e(j_bound) == 1) {
            di as res "    note: the minimum is ON THE BOUNDARY -- gamma at an end of"
            di as res "          the interval, or theta at its floor.  Neither the"
            di as res "          standard error nor the chi2(1) law for J is valid"
            di as res "          there, so no standard error is printed."
        }
        if (e(j_weak) == 1) {
            di as res "    note: theta < 0.05 at the minimum: the implied confounder"
            di as res "          has almost no variance, so kurt_U is not defined and"
            di as res "          there is nothing for the higher moments to bind on."
        }
        di as txt "    profiled GMM, method(pgmm)"
        _freeiv_fmt `=e(g_gmm)'
        di as txt "        gamma" _col(46) as res "`s'"
        di as txt "        J, 1 df" _col(46) as res %10.4f e(J) ///
           as txt "   P>chi2 " as res %6.4f e(p_J)
        di as txt "        h1..h5 are set to zero exactly, which is a DIFFERENT"
        di as txt "        estimator, not the same one solved differently."
        di as txt "    B kurt_U - A kurt_V2" _col(46) as res %10.4f e(idfac) ///
           as txt "   z " as res %6.2f e(z_idfac)
        di as txt "    The J tests alpha1 = gamma1 alpha2.  The Jacobian of the"
        di as txt "    model with alpha1 FREE has determinant -(gamma - tau)^5"
        di as txt "    times the factor above, so the restriction is testable"
        di as txt "    everywhere except where that factor vanishes -- and a"
        di as txt "    normal V2 sits exactly there, since it makes both B and"
        di as txt "    kurt_V2 zero."
        if (e(p_J) < 0.05) {
            di as txt "    The J rejects, so the rejection stands on its own and the"
            di as txt "    factor need not be read: it is computed at the restricted"
            di as txt "    estimates, which are not consistent under the alternative."
        }
        else if (abs(e(z_idfac)) < 2) {
            di as res "    The J does not reject AND the factor is not distinguishable"
            di as res "    from zero, so this non-rejection carries no information"
            di as res "    about scale consistency: the design is near the surface"
            di as res "    where the restriction cannot be tested at all."
        }
        else {
            di as txt "    The J does not reject and the factor is clearly non-zero,"
            di as txt "    so the non-rejection is informative."
        }
        di as txt "    The factor's magnitude is NOT a calibrated measure of power;"
        di as txt "    only its vanishing, under the null, is meaningful."
    }
    di as txt "{hline 72}"
    di as txt "s.e.: delta method on the stacked moments; the estimation of the"
    di as txt "coefficients on X is accounted for.  The QME's own standard error is"
    di as txt "proportional to 1/sqrt(D): it explodes as the discriminant nears zero."
end


program define _freeiv_pdisplay
    syntax [, noHEADer]
    if ("`e(cmd)'" != "freeiv" | "`e(model)'" != "B") exit
    local lev = e(level)
    cap confirm number `lev'
    if (_rc | "`lev'" == "") local lev = 95
    if (`lev' >= . | `lev' <= 0 | `lev' >= 100) local lev = 95
    local z = invnormal(1 - (100 - `lev') / 200)

    * under a prefix the prefix owns the coefficient table
    local pref "`e(prefix)'"
    if ("`pref'" != "" | e(N_reps) < .) {
        cap noisily ereturn display, level(`lev')
        cap noisily _freeiv_pident
        exit
    }

    di
    if ("`header'" == "") {
        di as txt "Instrument-free estimation, two indicators" _col(47) ///
           "Number of obs = " as res %9.0f e(n)
        di as txt "model: " as res "`e(depvar)'" as txt " on " ///
           as res "`e(endog2)' `e(endog3)'" ///
           as txt cond("`e(exog)'" != "", ", controls " + "`e(exog)'", "")
        di as txt "method: " as res "proxy" as txt ///
           "  (Theorem 1, one latent confounder behind both)"
    }
    di as txt "{hline 72}"

    if (e(guard) != 0) {
        local g = e(guard)
        di as res "no estimate: guard " as res %1.0f `g' as res " of Proposition 1 fired"
        if (`g' == 1) {
            di as txt "    E[eps2^2 eps3] and E[eps2 eps3^2] do not share a sign, or the"
            di as txt "    second is numerically zero: the third moment of the confounder"
            di as txt "    is too weak to orient the loadings"
        }
        if (`g' == 2) {
            di as txt "    E[eps2 eps3] and the ratio of the two third moments disagree"
            di as txt "    in sign, so a2 a3 would be negative: the two indicators are"
            di as txt "    not loading on one common factor with the same sign"
        }
        if (`g' == 3) {
            di as txt "    an implied idiosyncratic variance is negative: the common"
            di as txt "    factor would have to explain more than the whole variance of"
            di as txt "    an indicator"
        }
        if (`g' == 4) di as txt "    the 3x3 system of Theorem 1 is singular"
        di as txt "{hline 72}"
        exit
    }

    di as txt "Structural coefficients"
    _fivp_line g2 `=e(g2)' `=e(se_g2)' "`e(endog2)'" `z'
    _fivp_line g3 `=e(g3)' `=e(se_g3)' "`e(endog3)'" `z'
    _fivp_line a1 `=e(a1)' `=e(se_a1)' "U in the outcome" `z'
    di as txt "        a1 is free here: the two-indicator route needs no"
    di as txt "        scale-consistency restriction to identify it"
    cap noisily _freeiv_pident
end


program define _fivp_line
    args nm est se lab z
    _freeiv_fmt `est'
    local sg "`s'"
    _freeiv_fmt `se'
    local ss "`s'"
    if (`se' < .) {
        di as txt "    `nm' (" as res "`lab'" as txt ")" _col(34) ///
           as res "`sg'" as txt "  s.e. " as res "`=trim("`ss'")'" ///
           as txt "   [" as res %8.5f `est' - `z' * `se' as txt ", " ///
           as res %8.5f `est' + `z' * `se' as txt "]"
    }
    else di as txt "    `nm' (" as res "`lab'" as txt ")" _col(34) as res "`sg'"
end


program define _freeiv_pident
    di as txt "{hline 72}"
    di as txt "The confounder, recovered from the two indicators"
    _freeiv_fmt `=e(a2)'
    local s2 "`s'"
    _freeiv_fmt `=e(a3)'
    di as txt "    loadings a2, a3" _col(40) as res "`=trim("`s2'")'   `=trim("`s'")'"
    _freeiv_fmt `=e(mu3)'
    di as txt "    E[U^3] (mu3)" _col(40) as res "`s'"
    _freeiv_fmt `=e(s2)'
    local s2 "`s'"
    _freeiv_fmt `=e(s3)'
    di as txt "    sigma2_V2, sigma2_V3" _col(40) as res "`=trim("`s2'")'   `=trim("`s'")'"
    _freeiv_fmt `=e(sc)'
    di as txt "    g2 a2 + g3 a3" _col(40) as res "`s'" ///
       as txt "   (what model A would call a1)"
    di as txt "    condition number of the 3x3 system" _col(40) as res ///
       %10.4f e(cnum) cond(e(cnum) > 100, "   ill conditioned", "")

    di as txt "{hline 72}"
    di as txt "One-factor check: three estimates of a2/a3 that must agree"
    _freeiv_fmt `=e(R1)'
    di as txt "    R1  from eps2, eps3 alone" _col(40) as res "`s'"
    _freeiv_fmt `=e(R2)'
    di as txt "    R2  from xi eps^2" _col(40) as res "`s'"
    _freeiv_fmt `=e(R3)'
    di as txt "    R3  from xi^2 eps" _col(40) as res "`s'"
    * the ratio is read only when both third-order cross-moments are
    * measured (z >= 2): R2 and R3 divide by them and are noise otherwise,
    * which is the same gate as the menu applies
    local weak3 = (abs(e(z_m223)) < 2 | abs(e(z_m233)) < 2 | e(z_m223) >= . | e(z_m233) >= .)
    di as txt "    |R3/R1 - 1|" _col(40) as res %10.4f e(disc_R) ///
       cond(`weak3', "   not informative", ///
       cond(e(disc_R) < 0.15, "   consistent", "   a second factor is likely"))
    di as txt "    z of eps2^2 eps3, eps2 eps3^2" _col(40) as res ///
       %10.2f e(z_m223) "  " %8.2f e(z_m233) ///
       as txt cond(`weak3', "   both must reach 2", "")
    di as txt "    With one factor all three equal a2/a3.  With two, R1 becomes"
    di as txt "    (a2^2 a3 E[U^3] + b2^2 b3 E[W^3]) / (a2 a3^2 E[U^3] +"
    di as txt "    b2 b3^2 E[W^3]), which is a2/a3 only if the second factor is"
    di as txt "    symmetric or loads proportionally, and R3 weights the two"
    di as txt "    factors differently again.  A gap is therefore a second"
    di as txt "    factor, by derivation, not a sampling artefact."

    di as txt "    R1 - R3" _col(40) as res %10.6f e(of_d) ///
       as txt "  s.e. " as res %8.6f e(se_of) as txt "  z " as res %7.3f e(z_of)
    di as txt "    R1 equals a2/a3 under one factor with no side condition;"
    di as txt "    R3 needs symmetric V2 and V3 as well, so a rejection here is"
    di as txt "    against one factor AND that symmetry, jointly."

    if (e(q_J) < .) {
        di as txt "{hline 72}"
        di as txt "Over-identified GMM, 16 moments for 12 parameters"
        di as txt "    J, 4 df" _col(40) as res %10.4f e(q_J) ///
           as txt "   P>chi2 " as res %6.4f e(q_pJ)
        di as txt "    g2, g3, a1 from the GMM" _col(40) as res %10.6f e(q_g2) ///
           "  " %10.6f e(q_g3) "  " %10.6f e(q_a1)
        di as txt "    E[U^3], E[V2^3], E[V3^3]" _col(40) as res %10.6f e(q_mu3) ///
           "  " %10.6f e(q_k2) "  " %10.6f e(q_k3)
        di as txt "    The closed form uses eight of these sixteen moments and"
        di as txt "    fits them exactly; the J asks whether the eight it leaves"
        di as txt "    out agree, so it tests the one-factor linear structure"
        di as txt "    itself.  k2 and k3 are E[V2^3] and E[V3^3], the symmetry"
        di as txt "    that R3 -- and only R3 -- needs."
    }

    di as txt "{hline 72}"
    di as txt "Scale consistency, which model A cannot test"
    di as txt "    a1 (free, identified here)" _col(40) as res %10.6f e(a1)
    di as txt "    g2 a2 + g3 a3" _col(40) as res %10.6f e(sc)
    di as txt "    difference" _col(40) as res %10.6f e(sc_d) ///
       as txt "  s.e. " as res %8.6f e(se_sc) as txt "  z " as res %7.3f e(z_sc)
    di as txt "    Model A imposes alpha1 = gamma1 alpha2, which with two"
    di as txt "    endogenous regressors reads a1 = g2 a2 + g3 a3.  Its own nine"
    di as txt "    moments pin the ratio of the two loadings only up to that"
    di as txt "    normalisation; the second indicator frees a1 and makes the"
    di as txt "    restriction testable.  No fourth moment is needed."
    if (abs(e(z_sc)) > 1.96) {
        di as res "    The restriction is rejected on these data, so the"
        di as res "    single-indicator route would not be valid here."
    }

    di as txt "{hline 72}"
    di as txt "What `e(endog3)' used as an instrument for `e(endog2)' would return"
    di as txt "    iv estimate" _col(40) as res %10.6f e(ivgap)
    di as txt "    An indicator of the confounder is not an instrument: it is"
    di as txt "    correlated with the very thing it is meant to purge.  The gap"
    di as txt "    against g2 above is the bias that route would carry."
    di as txt "    correlation of the two indicators, t" _col(40) as res %10.4f e(t_rho)
    di as txt "{hline 72}"
    di as txt "s.e.: delta method, numeric Jacobian of Theorem 1 on the eight"
    di as txt "moments, with the estimation of the coefficients on X accounted for."
    di as txt "Theorem 1 involves square roots, so the linearisation is only good"
    di as txt "locally: in the reference design it tracked a 600-replication"
    di as txt "bootstrap to within about ten percent, in either direction.  Prefer"
    di as txt "the -bootstrap- prefix when the inference matters."
end
