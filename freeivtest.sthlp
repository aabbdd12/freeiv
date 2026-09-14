{smcl}
{* *! version 0.8.0  13sep2026}{...}
{vieweralsosee "freeiv" "help freeiv"}{...}
{vieweralsosee "freeivmenu" "help freeivmenu"}{...}
{vieweralsosee "freeivdiag" "help freeivdiag"}{...}
{title:Title}

{phang}
{bf:freeivtest} {hline 2} Endogeneity, agreement between routes, and an
outside estimate against the identified set


{title:Syntax}

{p 8 17 2}
{cmd:freeivtest} [{cmd:, gamma(}{it:#}{cmd:)} {cmd:segamma(}{it:#}{cmd:)}
{cmd:level(}{it:#}{cmd:)}]

{pstd}after {helpb freeiv}.


{title:Description}

{pstd}
Three blocks, in the order they should be read.


{title:1.  Is there any endogeneity at all}

{pstd}
Under theta = 0 the confounder is absent, so eps2 = V2 and the OLS residual
r = xi - gamma-tilde eps2 equals V1, which is {it:independent} of eps2 -- not
merely uncorrelated, which is true by construction.  Independence is
refutable, and two third-order cross moments must then vanish:

{p 8 8 2}g1 = E[r eps2^2] = m12 - gt m03 = gamma m03 (mu - lambda){p_end}
{p 8 8 2}g2 = E[r^2 eps2] = m21 - 2 gt m12 + gt^2 m03 = gamma^2 m03 (mu(1-2 lambda) + lambda^2){p_end}

{pstd}
The Wald statistic on the pair is chi2 with 2 df.  {bf:It never estimates
gamma}, so it keeps its power exactly where the qme loses its own: on
freeiv_sim2 the discriminant is negative and the QME has no real root, yet the
test rejects at p = 0.0004.

{pstd}
Both statistics vanish identically when m03 = 0, that is when U and V2 are
both symmetric.  There the test has no power {it:by construction}, and the
case is visible beforehand as z(m03) near zero in {helpb freeivmenu}.  Read
that first.

{pstd}
The test maintains the linear one-factor specification, so a rejection is
evidence of theta > 0 {it:or} of a departure from that structure.


{title:2.  Do the routes agree}

{pstd}
Every closed form is a smooth function of the same moment vector, so one
covariance matrix serves them all and the difference of any two has variance
(g_a - g_b)' V (g_a - g_b) / n.  Each pair is a test of the assumption that
separates the two routes:

{p2colset 8 24 26 2}{...}
{p2col:{bf:qme - ols}}theta = 0, no confounding{p_end}
{p2col:{bf:qme - rre}}k* = 1, the scale-free equal-variance restriction{p_end}
{p2col:{bf:qme - sce}}k = 1, the same in the units of the data{p_end}
{p2col:{bf:qme - hme}}B = 0, symmetry of V2{p_end}
{p2colreset}{...}

{pstd}
The QME is the reference because its assumptions are the weakest.  A large
|z| rejects the assumption that separates the pair; a small one is mutual
corroboration between two routes with different assumptions.

{pstd}
Note that {cmd:qme - ols} is {it:not} the right endogeneity test: it needs the
QME to exist and inherits its 1/sqrt(D) imprecision.  On wage1 it gives
p = 0.36 while block 1 gives p = 0.0002 on the same data.


{title:3.  An outside estimate against the identified set}

{pstd}
{cmd:gamma(}{it:#}{cmd:)} takes an estimate produced by any other method --
LSZ, Lewbel (2012), a published paper, a genuine instrument -- and asks
whether the scale-consistent model can produce it at all.  Outside the
interval it cannot, without a negative variance; {helpb freeivdiag} says which
one.  {cmd:segamma(}{it:#}{cmd:)} adds the outside estimate's own sampling
error, treating the two as independent, which they are not.


{title:What is NOT tested here}

{pstd}
No pair in block 2 tests scale consistency: the map from gamma to the
nuisances is derived from a1 = g1 a2, so every statistic presupposes it.  At
second order it is untestable outright, the bounds being exactly the
positivity region.

{pstd}
It {it:is} tested at fourth order, by the over-identifying J of
{cmd:method(gmm)} -- but only away from a knife edge.  The Jacobian of the
model with a1 free has determinant

{p 8 8 2}-(gamma - tau)^5 * (B kurt_U - A kurt_V2){p_end}

{pstd}
so the restriction binds except where that second factor vanishes, and a
{it:normal V2} sits exactly there, since it makes B and kurt_V2 both zero.
Read {cmd:e(idfac)} before the J.  In simulation, against a violation by a
factor two at n = 4000, the rejection rate at 5% is 0.065 with V2 normal,
0.435 with V2 symmetric but leptokurtic, and 1.000 with V2 skewed.

{pstd}
It is tested {it:directly}, without any fourth moment, by a second indicator.
In model B a1 is free and identified, and {cmd:freeiv} reports
{cmd:e(sc_d)} = a1 - (g2 a2 + g3 a3) with its standard error.  That is the
only unconditional test of scale consistency in this package.


{title:Examples}

{pstd}These run on data that ship with Stata, so they can be executed as they
stand{p_end}
{phang2}{cmd:. sysuse nlsw88, clear}{p_end}
{phang2}{cmd:. generate lwage = ln(wage)}{p_end}

{phang2}{cmd:. freeiv lwage grade age i.race (tenure)}{p_end}
{phang2}{cmd:. freeivtest}{p_end}

{pstd}Against an estimate obtained elsewhere, with its own standard error{p_end}
{phang2}{cmd:. freeivtest, gamma(0.03) segamma(0.01)}{p_end}


{title:Stored results}

{pstd}
{cmd:freeivtest} is {cmd:rclass} and returns {cmd:r(W_endo)},
{cmd:r(p_endo)}, {cmd:r(g1)}, {cmd:r(g2)}, the pairwise {cmd:r(z_ols)},
{cmd:r(z_rre)}, {cmd:r(z_sce)}, {cmd:r(z_hme)}, and {cmd:r(out_z)},
{cmd:r(out_d)} when {cmd:gamma()} is given.


{title:Author}

{pstd}
Abdelkrim Araar, Universite Laval and PEP{break}
{browse "mailto:aabd@ecn.ulaval.ca":aabd@ecn.ulaval.ca}


{title:Also see}

{psee}
Online: {helpb freeiv}, {helpb freeivmenu}, {helpb freeivdiag}
{p_end}
