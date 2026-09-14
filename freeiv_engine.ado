*! freeiv_engine 0.2.0  14sep2026  A. Araar
*! Loader of the Mata engine.  An ado-file loaded automatically does not
*! execute its mata: block, so freeiv_mata.ado is run on demand, once.
*!
*! 0.2.0 checks EVERY entry point, not just the first.  Guarding on
*! _freeiv_all() alone was wrong in one ordinary case: a session that had
*! already used an older freeiv still holds that older engine in memory, the
*! guard reports it as loaded, and a newly added function -- _freeiv_jgmm()
*! in 0.8.0 -- is never defined.  The user then meets "_freeiv_jgmm() not
*! found" with no hint that the cure is to reload.  Naming all of them makes
*! an upgrade inside a live session reload, which is what the guard was for.
*! A new brick adds its entry point to the list below.

cap program drop _freeiv_engine_ck
cap program drop freeiv_engine

* sets s(fiv_missing) to the entry points not currently in memory, if any
program define _freeiv_engine_ck, sclass
    version 16
    args need
    sreturn clear
    local gone ""
    foreach f of local need {
        mata: st_local("ok", strofreal(findexternal("`f'()") != NULL))
        if ("`ok'" != "1") local gone "`gone' `f'()"
    }
    sreturn local fiv_missing = trim("`gone'")
end

program define freeiv_engine
    version 16
    local need "_freeiv_all _freeiv_proxy _freeiv_lit _freeiv_gmm"
    local need "`need' _freeiv_jgmm _freeiv_lsz _freeiv_ptests _freeiv_gmm16"

    _freeiv_engine_ck "`need'"
    if ("`s(fiv_missing)'" == "") exit

    cap findfile freeiv_mata.ado
    if (_rc) {
        cap findfile freeiv_mata.ado, path("stata")
    }
    if (_rc) {
        di as err "freeiv_mata.ado not found on the adopath"
        exit 601
    }
    local fn "`r(fn)'"
    cap noisily version `c(stata_version)': run "`fn'"

    _freeiv_engine_ck "`need'"
    if ("`s(fiv_missing)'" != "") {
        di as err "the freeiv Mata engine could not be loaded:"
        di as err "    `s(fiv_missing)' still undefined after running"
        di as err "    `fn'"
        di as err "the ado files and freeiv_mata.ado are out of step;"
        di as err "reinstall the package, or -clear all- and try again"
        exit 601
    }
end
