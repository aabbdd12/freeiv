{smcl}
{* *! version 0.8.0  13sep2026}{...}
{vieweralsosee "freeiv" "help freeiv"}{...}
{vieweralsosee "freeivmenu" "help freeivmenu"}{...}
{vieweralsosee "freeivdiag" "help freeivdiag"}{...}
{vieweralsosee "freeivtest" "help freeivtest"}{...}
{title:Title}

{phang}
{bf:freeivreport} {hline 2} The whole instrument-free pipeline in one table


{title:Syntax}

{p 8 17 2}
{cmd:freeivreport} {depvar} [{indepvars}]
{cmd:(}{it:endogvar}[ {it:endogvar2}]{cmd:)}
{ifin} {weight}
[{cmd:,} {it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt sav:ing(filename[, replace])}}also write the table as a
comma-delimited file{p_end}
{synopt:{opt nomenu}}skip the identification block{p_end}
{synopt:{opt notest:s}}skip the test block{p_end}
{synopt:{opt l:evel(#)}, {opt quan:tile(#)}, {opt del:ta(#)}, {opt rmax(#)}, {opt sign(#)}}as in {helpb freeiv}{p_end}
{synoptline}
{p2colreset}{...}


{title:Description}

{pstd}
{cmd:freeivreport} runs {helpb freeivmenu}, {helpb freeiv} with
{cmd:method(all)} and {helpb freeivtest} on the same sample and assembles them
into one table.  It exists because these three should never be read apart, and
because a table of this shape is what an empirical section needs.

{pstd}
The order is the order in which the results should be read, and it is not
negotiable:

{phang}{bf:1.  Identification.}  The width of the identified set, then the
signal each family of strategies needs -- the z of m03 for the third-order
route, the z of the discriminant, the heteroskedasticity F for Lewbel (2012),
k* for the equal-variance route, B for the symmetry route.  A route whose
signal is absent will still return a number; that number means nothing.{p_end}

{phang}{bf:2.  Estimates.}  Every route, with its standard error where one
exists, and a column saying whether it falls inside the identified set.  An
estimate outside implies a negative variance under scale consistency; run
{helpb freeivdiag}{cmd:, gamma(}{it:#}{cmd:)} to see which one.  The set
itself is printed on the last line, so no estimate is ever shown without
it.{p_end}

{phang}{bf:3.  Tests.}  Whether there is any endogeneity at all, whether the
routes agree with each other, and the over-identifying J with the factor that
says whether that J can test anything.{p_end}

{pstd}
With two endogenous variables the report switches to model B and shows g2, g3
and the free loading a1, then the three tests that model A cannot perform: the
scale-consistency test a1 - (g2 a2 + g3 a3), the one-factor test R1 - R3, and
the J of the sixteen-moment GMM.


{title:Options}

{phang}
{opt saving(filename[, replace])} writes the same table as a comma-delimited
file with a header line and one row per line: {it:block}, {it:label},
{it:value}, {it:se}, {it:note}.  The blocks are {cmd:identification},
{cmd:estimate} and {cmd:test}, so the file can be reshaped into a paper table
without editing.

{phang}
{opt nomenu} and {opt notests} trim the report to the estimates alone.  Use
them for a loop over many specifications, not for a result you intend to
report.

{phang}
The remaining options are passed through to {helpb freeiv} unchanged.


{title:Examples}

{pstd}These run on data that ship with Stata, so they can be executed as they
stand{p_end}
{phang2}{cmd:. sysuse nlsw88, clear}{p_end}
{phang2}{cmd:. generate lwage = ln(wage)}{p_end}

{phang2}{cmd:. freeivreport lwage grade age i.race (tenure)}{p_end}

{pstd}Writing the table out for a paper{p_end}
{phang2}{cmd:. freeivreport lwage grade age i.race (tenure), saving(tab1.csv, replace)}{p_end}

{pstd}Two endogenous regressors sharing one confounder{p_end}
{phang2}{cmd:. net get freeiv}{p_end}
{phang2}{cmd:. use freeiv_proxy, clear}{p_end}
{phang2}{cmd:. freeivreport y1 x (y2 y3)}{p_end}

{pstd}A loop over specifications, estimates only{p_end}
{phang2}{cmd:. sysuse nlsw88, clear}{p_end}
{phang2}{cmd:. generate lwage = ln(wage)}{p_end}
{phang2}{cmd:. foreach c of varlist grade ttl_exp {c -(}}{p_end}
{phang2}{cmd:.     freeivreport lwage `c' age (tenure), nomenu notests}{p_end}
{phang2}{cmd:. {c )-}}{p_end}


{title:Stored results}

{pstd}
{cmd:freeivreport} is {cmd:rclass} and returns {cmd:r(n)}, {cmd:r(lo)},
{cmd:r(hi)}, {cmd:r(qme)}, {cmd:r(depvar)} and {cmd:r(endog)}.  Everything
else remains in {cmd:e()} from the {helpb freeiv} call it made, so
{cmd:ereturn list} still works afterwards.


{title:Author}

{pstd}
Abdelkrim Araar, Universite Laval and PEP{break}
{browse "mailto:aabd@ecn.ulaval.ca":aabd@ecn.ulaval.ca}


{title:Also see}

{psee}
Online: {helpb freeiv}, {helpb freeivmenu}, {helpb freeivdiag},
{helpb freeivtest}
{p_end}
