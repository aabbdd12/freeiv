{smcl}
{* *! version 0.10.2  15sep2026}{...}
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

{p 8 17 2}
{cmd:freeivmenu} {depvar} [{indepvars}] {cmd:(}{it:endogvar1} {it:endogvar2}{cmd:)}
{ifin} {weight}
[{cmd:, bw(}{it:#}{cmd:)}]

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
the law; "tightens" is printed when the minimum sits at p5 or p10 and
beats gamma-tilde by more than 1.96 standard errors.  The minimum belongs
in a tail: on the side of eps2 that V2 dominates the slope tends to gamma,
on the side U dominates to gamma + alpha1/alpha2, and an increasing m
keeps it above gamma in between.  An interior minimum is reported as
"interior min": the profile does not have that shape, no bound is read
from it, and a non-linear outcome equation is the first thing to suspect.
And the ratio of the two tail slopes tends to
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
{bf:With two indicators.}  Two variables in the parentheses give the menu of
the two-indicator model, which says before estimation whether Theorem 1 of
Araar (2026d) has anything to work with:

{p2colset 5 26 28 2}{...}
{p2col:{bf:relevance of the pair}}the t of the residual correlation of the two
indicators, against the application rule t >= 10 below which the closed form
disperses steeply, and alpha2 alpha3 = E[eps2 eps3] with its z.{p_end}
{p2col:{bf:third order}}the two cross-moments E[eps2^2 eps3] and
E[eps2 eps3^2] with their z, and whether they share a sign -- guard (i) of
Proposition 1, which fires when the regressor affects the indicator or a
second factor is present.{p_end}
{p2col:{bf:loadings}}alpha2/alpha3 from the ratio of the two cross-moments,
guard (ii) on its coherence with the covariance, the implied loadings, and
the implied idiosyncratic variances, guard (iii).  A guard that fires here
fires in {helpb freeiv}.{p_end}
{p2col:{bf:confounder}}the implied skewness of U, which sets the precision
regime mapped in Araar (2026d) -- {it:low} (below 0.5): the causal
coefficient is recovered but the free direct effect a1 is not;
{it:moderate}; {it:full} (above 1.5): all parameters are precise and scale
consistency becomes testable.{p_end}
{p2col:{bf:one factor}}the three estimates R1, R2, R3 of alpha2/alpha3 and
|R3/R1 - 1|; read only when both third-order z exceed 2, since R2 and R3
are noise otherwise.{p_end}
{p2col:{bf:local slope profile}}of xi on each indicator, as above; "rising"
or "interior min" says whether the minimum sits in a tail, as the linear
one-factor model implies, or inside.  A guard that fires is explained on
the line below its verdict.{p_end}
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
{cmd:r(prof_ratio)} and {cmd:r(prof_z)}.  With two indicators the returned
set is that of the two-indicator engine ({cmd:r(t_rho)}, {cmd:r(m23)},
{cmd:r(m223)}, {cmd:r(m233)}, {cmd:r(guard)}, {cmd:r(R1)}, {cmd:r(R2)},
{cmd:r(R3)}, {cmd:r(disc_R)}, {cmd:r(a2)}, {cmd:r(a3)}, {cmd:r(mu3)},
{cmd:r(s2)}, {cmd:r(s3)}, ...), the z statistics {cmd:r(z_m23)},
{cmd:r(z_m223)}, {cmd:r(z_m233)}, and the two profiles {cmd:r(profile2)},
{cmd:r(profile3)} with their {cmd:r(prof2_*)} and {cmd:r(prof3_*)} scalars.


{title:Author}

{pstd}
Abdelkrim Araar, Universite Laval and PEP{break}
{browse "mailto:aabd@ecn.ulaval.ca":aabd@ecn.ulaval.ca}


{title:Also see}

{psee}
Online: {helpb freeiv}, {helpb freeivdiag}, {helpb freeivtest}
{p_end}
