{smcl}
{* *! version 0.9.1  15sep2026}{...}
{vieweralsosee "freeiv" "help freeiv"}{...}
{vieweralsosee "freeivdiag" "help freeivdiag"}{...}
{vieweralsosee "freeivtest" "help freeivtest"}{...}
{title:Title}

{phang}
{bf:freeivmenu} {hline 2} What these data can identify, before any estimation


{title:Syntax}

{p 8 17 2}
{cmd:freeivmenu} {depvar} [{indepvars}] {cmd:(}{it:endogvar}{cmd:)}
{ifin} {weight}
[{cmd:, quantile(}{it:#}{cmd:)} {cmd:bw(}{it:#}{cmd:)}]

{p 4 6 2}{it:fweight}s, {it:pweight}s and {it:aweight}s are allowed.{p_end}


{title:Description}

{pstd}
{cmd:freeivmenu} is run {it:before} {helpb freeiv}, in the same spirit as
reading a first-stage F before an instrumental-variables regression.  It
prints one line per family of identifying strategies, with the signal that
family needs, the value of that signal in these data, and a verdict.

{pstd}
The first line is always true: the identified set exists under the model
alone.  The others say whether the assumption that would deliver a {it:point}
is carried by the data.  A weak signal does not make an estimator wrong; it
makes it imprecise, and sometimes it makes it silent.

{pstd}
What each line reports:

{p2colset 5 26 28 2}{...}
{p2col:{bf:bounds}}the width of [gamma-tilde/2, gamma-tilde].  Always
available.{p_end}
{p2col:{bf:third order}}the skewness of the first-stage residual and the z of
m03.  The third-order route -- qme, hme, lsz -- lives entirely on this.  Also
the z of the discriminant D, the implied standard error of the QME, and
mu = A/(A+B), which near 1/3 is where the two roots merge and the standard
error explodes.{p_end}
{p2col:{bf:equal variances}}k* = sV1/(g^2 sV2), which the {cmd:rre} sets to 1,
and below it k = sV1/sV2, which the {cmd:sce} sets to 1.  Only k* is a
statement about the model: k moves when either variable is rescaled.{p_end}
{p2col:{bf:symmetry of V2}}B = E[V2^3], which the {cmd:hme} sets to 0.{p_end}
{p2col:{bf:heteroskedasticity}}the F of eps2^2 on X, which is what
{cmd:lewbel12} needs.  When it is absent, that route returns noise -- and
{cmd:ivreg2h} says the same thing in its own language, through a small
Cragg-Donald F.{p_end}
{p2col:{bf:local slope profile}}the kernel-weighted slope of xi on eps2 at
seven percentiles of eps2 -- the local-linear estimate of the derivative of
E[xi | eps2].  In the one-factor model that derivative is
gamma + alpha1 m'(s) with m(s) = E[U | eps2 = s].  If the confounder is
Gaussian-like, m is linear and the profile is {it:flat} at gamma-tilde: the
third-order routes will find little.  If it is skewed, m is curved and the
profile moves -- but a non-linear outcome equation curves it too, so a
curved profile says which of the two is in play only in conjunction with
the other lines.  Two readings are exact in the limit.  When m is
increasing (log-concave densities suffice) the {it:minimum} of the profile
is an upper bound on gamma tighter than gamma-tilde, with no assumption on
the law; "tightens" is printed when the minimum sits at p5 or p10 -- the
rising profile an increasing m implies -- and beats gamma-tilde by more
than 1.96 standard errors.  An interior minimum is reported as
"interior min": the profile is not monotone, the bound reading does not
apply, and a non-linear outcome equation is the first thing to suspect.  And the ratio of the two tail slopes tends to
1 + alpha1/(gamma alpha2), which is 2 under scale consistency; it is
informative only when one tail of eps2 is dominated by U and the other by
V2, as with a bounded or strongly skewed confounder, and it converges
slowly.  The bandwidth is h = 2 * 1.06 * sd(eps2) * N^(-1/5) unless
{opt bw()} is given.{p_end}
{p2col:{bf:two indicators}}a reminder that a second indicator of the same
confounder, if one exists, opens model B and with it the only direct test of
scale consistency.  Model B is selected by the syntax alone -- two variables
inside the parentheses -- and there is no {cmd:method()} name for it.{p_end}
{p2colreset}{...}

{pstd}
{bf:Read z of m03 first.}  Detecting endogeneity and measuring it need
different signals: measuring g1 needs A = a2^3 E[U^3] != 0, a {it:skewed
confounder}, while detecting theta > 0 needs only m03 != 0 from either U or
V2.  When U and V2 are both symmetric the third-order block is empty and
neither is possible; that case is visible here and nowhere else.


{title:Options}

{phang}
{opt quantile(#)} is passed through to the Q-BE diagnostics; default 0.25.

{phang}
{opt bw(#)} sets the kernel bandwidth of the local slope profile, in the
units of eps2; the default is 2 * 1.06 * sd(eps2) * N^(-1/5).  A larger
bandwidth gives a smoother, more precise and more biased profile.


{title:Examples}

{pstd}These run on data that ship with Stata, so they can be executed as they
stand{p_end}
{phang2}{cmd:. sysuse nlsw88, clear}{p_end}
{phang2}{cmd:. generate lwage = ln(wage)}{p_end}

{pstd}A specification these data can carry{p_end}
{phang2}{cmd:. freeivmenu lwage grade age i.race (tenure)}{p_end}
{phang2}{cmd:. freeiv lwage grade age i.race (tenure), method(all)}{p_end}

{pstd}And one they cannot, in the same file: the third-order signal is absent
and the discriminant is negative, so the QME falls back on the vertex and
lands outside its own identified set{p_end}
{phang2}{cmd:. freeivmenu lwage ttl_exp (grade)}{p_end}
{phang2}{cmd:. freeiv lwage ttl_exp (grade), method(all)}{p_end}


{title:Stored results}

{pstd}
{cmd:freeivmenu} is {cmd:rclass}.  Every displayed statistic is returned,
including {cmd:r(z_m03)}, {cmd:r(F_lewbel)}, {cmd:r(p_lewbel)},
{cmd:r(lo)}, {cmd:r(hi)}, {cmd:r(k)}, {cmd:r(kstar)}, {cmd:r(A)}, {cmd:r(B)},
{cmd:r(mu)}, {cmd:r(disc)}, {cmd:r(disc_z)}, and the matrix
{cmd:r(moments)}.  The local slope profile is returned as the 7 x 3 matrix
{cmd:r(profile)} (evaluation point, slope, standard error, rows p5 to p95),
with {cmd:r(prof_h)}, {cmd:r(prof_min)}, {cmd:r(prof_minse)},
{cmd:r(prof_ratio)} and {cmd:r(prof_z)}.


{title:Author}

{pstd}
Abdelkrim Araar, Universite Laval and PEP{break}
{browse "mailto:aabd@ecn.ulaval.ca":aabd@ecn.ulaval.ca}


{title:Also see}

{psee}
Online: {helpb freeiv}, {helpb freeivdiag}, {helpb freeivtest}
{p_end}
