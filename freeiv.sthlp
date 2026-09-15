{smcl}
{* *! version 0.16.5  15sep2026}{...}
{vieweralsosee "[R] ivregress" "help ivregress"}{...}
{vieweralsosee "freeivmenu" "help freeivmenu"}{...}
{vieweralsosee "freeivdiag" "help freeivdiag"}{...}
{vieweralsosee "freeivtest" "help freeivtest"}{...}
{viewerjumpto "Syntax" "freeiv##syntax"}{...}
{viewerjumpto "Description" "freeiv##description"}{...}
{viewerjumpto "Options" "freeiv##options"}{...}
{viewerjumpto "The routes, one by one" "freeiv##routes"}{...}
{viewerjumpto "How to read the output" "freeiv##reading"}{...}
{viewerjumpto "Which option acts on which route" "freeiv##whichopt"}{...}
{viewerjumpto "What freeiv refuses" "freeiv##refused"}{...}
{viewerjumpto "Examples" "freeiv##examples"}{...}
{viewerjumpto "Stored results" "freeiv##results"}{...}
{viewerjumpto "References" "freeiv##references"}{...}
{title:Title}

{phang}
{bf:freeiv} {hline 2} Instrument-free estimation of a linear structural model


{marker syntax}{...}
{title:Syntax}

{pstd}One endogenous regressor (model A){p_end}

{p 8 17 2}
{cmd:freeiv} {depvar} [{indepvars}] {cmd:(}{it:endogvar}{cmd:)}
{ifin} {weight}
[{cmd:,} {it:options}]

{pstd}Two endogenous regressors sharing one latent confounder (model B){p_end}

{p 8 17 2}
{cmd:freeiv} {depvar} [{indepvars}] {cmd:(}{it:endogvar1} {it:endogvar2}{cmd:)}
{ifin} {weight}
[{cmd:,} {cmd:level(}{it:#}{cmd:)} {cmd:noheader}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt meth:od(name)}}estimator; default is {cmd:method(qme)}{p_end}
{synopt:{opt quan:tile(#)}}sub-sample quantile used by {cmd:method(qbe)}; default {cmd:quantile(0.25)}{p_end}
{synopt:{opt del:ta(#)}}selection parameter of {cmd:method(oster)}; default {cmd:delta(1)}{p_end}
{synopt:{opt rmax(#)}}R-squared ceiling of {cmd:method(oster)}; default min(1.3 R1, 1){p_end}
{synopt:{opt sign(#)}}orientation of {cmd:method(lsz)}, 1 or -1; default {cmd:sign(1)}{p_end}

{syntab:Reporting}
{synopt:{opt l:evel(#)}}confidence level; default {cmd:level(95)}{p_end}
{synopt:{opt nohead:er}}suppress the header{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}{it:fweight}s, {it:pweight}s and {it:aweight}s are allowed.{p_end}
{p 4 6 2}{it:indepvars} may contain factor variables -- {cmd:i.sex},
{cmd:ib2.region}, {cmd:i.sex##c.age}.  They are expanded with the base level
omitted, and the estimate is identical to the one obtained from indicators
built by hand: dropping any one level of a full dummy set spans the same
column space, so the partialled-out residuals do not move.  {it:depvar} and
the endogenous variables must be plain numeric variables.{p_end}
{p 4 6 2}{cmd:bootstrap}, {cmd:jackknife} and {cmd:svy} prefixes are allowed;
{cmd:e(b)} and {cmd:e(V)} are posted whenever an estimate exists.{p_end}

{pstd}{it:name} in {cmd:method(}{it:name}{cmd:)} is one of{p_end}

{synoptset 14 tabbed}{...}
{synopt:{opt bounds}}the identified set, which needs no assumption beyond the model{p_end}
{synopt:{opt ols}}the naive slope, which is also the upper bound{p_end}
{synopt:{opt qme}}quadratic moment estimator; needs only E[U^3] != 0 {bf:(default)}{p_end}
{synopt:{opt sce}}cubic under Var(V1) = Var(V2){p_end}
{synopt:{opt rre}}reverse-regression estimator, under Var(V1) = g^2 Var(V2){p_end}
{synopt:{opt hme}}cube root under symmetric V1 and V2{p_end}
{synopt:{opt qbe}}quantile bias extrapolation{p_end}
{synopt:{opt gmm}}joint GMM on the nine moments of orders 2, 3 and 4{p_end}
{synopt:{opt pgmm}}the same nine moments, profiled{p_end}
{synopt:{opt lsz}}Lewbel, Schennach and Zhang (2024), as {cmd:trigmm} implements it{p_end}
{synopt:{opt lewbel12}}Lewbel (2012), heteroskedasticity-based{p_end}
{synopt:{opt copula}}Park and Gupta (2012){p_end}
{synopt:{opt rank}}Breitung, Mayer and Wied (2024){p_end}
{synopt:{opt rpiv}}residual purging IV; {it:proven inconsistent}, kept for comparison{p_end}
{synopt:{opt ape}}adaptive rule between OLS and RPIV{p_end}
{synopt:{opt oster}}Oster (2019), coefficient stability{p_end}
{synopt:{opt all}}every route above, side by side, against the bounds{p_end}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:freeiv} estimates the coefficient on an endogenous regressor without an
external instrument.  The model is

{p 8 8 2}Y2 = X'b2 + eps2,{space 6}eps2 = a2 U + V2{p_end}
{p 8 8 2}Y1 = X'b1 + g1 Y2 + eps1,{space 2}eps1 = a1 U + V1{p_end}

{pstd}
with U the unobserved confounder and V1, V2 independent of it and of each
other.  The identifying restriction of model A is {it:scale consistency},
a1 = g1 a2: the confounder's direct loading in the outcome equation is
proportional to its loading in the first stage, with the structural
coefficient as the factor.

{pstd}
Under that restriction alone, and with no distributional assumption at all,
g1 is bounded:

{p 8 8 2}g1 is in [ gamma-tilde / 2 , gamma-tilde ]{p_end}

{pstd}
where gamma-tilde is the OLS slope.  That interval is not a convenience, it is
exactly the region where the implied variances are non-negative, so it is the
most that can be said without a further assumption.  Every point estimator in
this command buys its point by adding one, and {cmd:method(all)} shows them
together with the interval so that each can be judged against it.

{pstd}
Two variables inside one parenthesis block switch to {bf:model B}, where Y2 and
Y3 are two indicators of the same U and a1 is free.  That is not merely a
second estimator: it is what makes scale consistency testable, since model A's
own moments pin the ratio of the two loadings only up to that normalisation.
See {helpb freeivtest}.


{marker options}{...}
{title:Options}

{phang}
{opt method(name)} selects the route.  Every route estimates the same g1; what
separates them is the assumption each adds to close the system, and those
assumptions are not comparable in strength.  The default, {cmd:qme}, adds the
weakest one in the family -- a non-zero third moment of the confounder -- and
{cmd:method(all)} computes every route and prints them against the identified
interval, which is how they are meant to be read.  The option does not apply
to the two-indicator model, where a single closed form is the whole method,
and it is refused there rather than ignored.

{phang}
{opt quantile(#)} acts on {cmd:method(qbe)} alone.  The Q-BE estimates the
bias where it is smallest and extrapolates: it runs OLS on the sub-sample with
the smallest {it:signed} first-stage residuals -- the observations least
contaminated by the confounder -- and projects back to the full sample.
{opt quantile()} is the share taken, the default being the lowest quartile.  A
smaller share is less contaminated but noisier; the estimate should be stable
across neighbouring values, and a large movement between 0.2 and 0.3 is itself
a finding about the design.  Values must lie strictly between 0 and 1.

{phang}
{opt delta(#)} and {opt rmax(#)} act on {cmd:method(oster)} alone, and neither
is estimated: both are {it:assumed}, and the answer moves with them.
{cmd:delta} is Oster's proportional-selection parameter -- how much the
unobservables matter relative to the observables already in the model, with
1 meaning "as much" -- and {cmd:rmax} is the R-squared a regression including
those unobservables would reach, defaulting to min(1.3 R1, 1) where R1 is the
R-squared of the controlled regression.  The honest use of this route is to
report a range obtained by varying both, not one number.  {cmd:rmax} must lie
in (0, 1]; a value below the observed R1 is a contradiction the formula cannot
absorb.  Both require at least one exogenous control, since the method
compares the R-squared with and without them.

{phang}
{opt sign(#)} acts on {cmd:method(lsz)} alone and takes 1 or -1.  The LSZ
moment system identifies the coefficient up to an orientation, and the option
selects which branch the solver follows -- the same choice {cmd:trigmm} offers.
It is not a preference: the two branches are different estimates, and the one
to keep is the one whose implied parameters are admissible.  Any other value
is refused.

{phang}
{opt level(#)} sets the confidence level.  It changes nothing for the seven
routes that carry no analytic standard error, since no interval is printed for
them -- but a {cmd:bootstrap:} or {cmd:jackknife:} prefix owns the coefficient
table when one is used, and reads {opt level()} itself, which is precisely the
case where it matters.

{phang}
{opt noheader} suppresses the header.  It applies to every route and to both
models.


{marker routes}{...}
{title:The routes, one by one}

{pstd}
Read this with {helpb freeivmenu} open: for each route the menu prints the
statistic that says whether these data can carry it.  A route whose signal is
absent still returns a number, and that number means nothing.

{pstd}
{bf:bounds} {hline 2} the identified set{p_end}

{pstd}
Not an estimator but the honest answer.  Under scale consistency alone, with
no distributional assumption whatever, g1 lies in [gamma-tilde/2,
gamma-tilde], gamma-tilde being the OLS slope m11/m02.  The interval is not a
convention: it is exactly the set of values for which the implied variances
theta, sigma2_V1 and sigma2_V2 are all non-negative, so nothing outside it is
compatible with the model.  Its width is gamma-tilde/2, which also says that in
magnitude OLS can overstate the coefficient by at most a factor of two under
this restriction.  Every other route buys a point by adding one assumption to this.
(Araar 2026a, 2026c.)

{pstd}
{bf:ols} {hline 2} the naive slope{p_end}

{pstd}
Reported because it {it:is} the upper bound, not as a competitor.  Seeing it
printed beside the interval makes the arithmetic visible: the confounder can
only inflate the slope, never deflate it, under scale consistency.

{pstd}
{bf:qme} {hline 2} quadratic moment estimator, the default{p_end}

{pstd}
The root of 2 m03 g^2 - 3 m12 g + m21 = 0, formed from the third-order moments
of the two residuals.  Its only requirement beyond the model is E[U^3] != 0: a
skewed confounder.  Nothing is assumed about variances, about symmetry of V,
or about the shape of anything else, which is what makes it the default.  Read
the z of m03 in {helpb freeivmenu} first -- a symmetric confounder leaves this
route empty, and the estimate is then noise.  The quadratic has two roots and the
command reports how many of them fall inside the identified set.  None need
do: when that count is zero the model is refuted on these data, and the
command says so rather than printing the number quietly.  When the discriminant
D = 9 m12^2 - 8 m03 m21 is negative the two roots merge and the vertex
3 m12/(4 m03) is returned instead, which is flagged in the output.  The
standard error is proportional to 1/sqrt(D) and therefore explodes as D nears
zero -- which is the same statement as "the two roots are about to become
complex".  (Araar 2026b, 2026c.)

{pstd}
{bf:hme} {hline 2} higher moments, symmetric V{p_end}

{pstd}
The closed form (1/2)(m30/m03)^(1/3), which needs both E[U^3] != 0 and V1 and
V2 symmetric, that is E[V2^3] = 0.  It was the central estimator of the early
work in this line before the QME replaced it, because the symmetry it requires
is a real restriction on the data rather than on the confounder alone.  The
command reports B = E[V2^3] so that the assumption can be read directly rather
than believed.  (Araar 2026a.)

{pstd}
{bf:sce} {hline 2} cubic under equal variances{p_end}

{pstd}
The admissible root of 2 g^3 - 3 gamma-tilde g^2 + (R - 2) g + gamma-tilde = 0
with R = m20/m02, under Var(V1) = Var(V2).  That restriction compares a
variance measured in the units of {it:depvar} with one measured in the units
of {it:endogvar}, so it is not scale free: multiply {it:depvar} by a thousand
and the estimate changes.  It is reported for continuity with the literature
and is not recommended as an estimator; its useful role is as a test, since
the QME identifies k = sigma2_V1/sigma2_V2 and k = 1 can then be read against
the data.  (Araar 2026a.)

{pstd}
{bf:rre} {hline 2} reverse regression, the scale-free counterpart{p_end}

{pstd}
The closed form m20/(2 m11), which imposes Var(V1) = g^2 Var(V2) instead --
the same idea as the SCE but stated in a way that survives a change of units.
The quantity to judge it is k* = sigma2_V1/(g^2 sigma2_V2), printed by both
{helpb freeivmenu} and {helpb freeivdiag}: k* is invariant to the scale of
either variable, k is not.  By Cauchy-Schwarz the RRE never falls below the
lower bound, and it breaches the upper bound exactly when the restriction is
rejected by the data, so a value outside the interval is informative rather
than a failure.  This route does not appear in the four papers; it was derived
for this package as the scale-free repair of the SCE.

{pstd}
{bf:qbe} {hline 2} quantile bias extrapolation{p_end}

{pstd}
See {opt quantile()} above for the mechanism.  It assumes only that the bias
varies smoothly with the signed first-stage residual, which is weaker than
what the closed forms assume, and it pays for that with no analytic standard
error: use the {cmd:bootstrap:} prefix.  (Araar 2026c.)

{pstd}
{bf:gmm} {hline 2} the joint nine-moment GMM of Araar (2026c){p_end}

{pstd}
Nine moment conditions in eight free parameters (g, theta, sV1, sV2, A, B,
A4, B4), so one over-identifying restriction and a J with 1 degree of
freedom.  This is the estimator the paper reports, and {cmd:freeiv}
reproduces its published table.  It is minimised deterministically, which
matters more here than it usually does: fix g {it:and} theta and the only
non-linearity left, the product theta*sV2, disappears, so all nine residuals
become linear in the six remaining parameters.  The eight-parameter problem
is therefore a two-dimensional surface, and a surface can be gridded.  There
is no seed and no starting value.  Multi-start over the eight parameters is
not adequate: on six of the eight datasets of Araar (2026c) it returns a
local minimum, and on one of them the Hansen statistic is twice its true
value.{p_end}

{pstd}
Because the criterion is often nearly flat, the route reports
{cmd:e(ar_lo)} and {cmd:e(ar_hi)}: the values of g at which the profiled J
is within 3.84 of its minimum.  That region inverts the J at each fixed g
rather than inverting a Wald statistic around the optimum, so it stays valid
where the standard error does not -- on the boundary, and on a flat
criterion.  When it fills the whole identified interval, as it does on seven
of the paper's eight applications, the moments of orders three and four have
added nothing to the assumption-free bounds, and the interval rather than
the point is what should be reported.  {cmd:e(j_bound)} marks a minimum on
the boundary of the parameter space, where no standard error is printed, and
{cmd:e(j_weak)} marks theta below 0.05.{p_end}

{pstd}
The J tests scale consistency {it:and} linearity jointly -- a rejection does
not say which failed.  Read {cmd:e(idfac)} before concluding from a
non-rejection: the Jacobian of the model with a1 {it:free} has determinant
-(g - tau)^5 times (B kurt_U - A kurt_V2), so the restriction is testable
everywhere except where that second factor vanishes, and a normal V2 sits
exactly there -- it makes both B and kurt_V2 zero.  The factor's
{it:magnitude} is not a calibrated power index; only its vanishing is
meaningful.  The search runs over [gamma-tilde/2, gamma-tilde], so the
estimate cannot leave the identified set.  (Araar 2026c.)

{pstd}
{bf:pgmm} {hline 2} the same nine moments, profiled{p_end}

{pstd}
The first five conditions are set to zero {it:exactly}, which solves theta,
sV1, sV2, A and B in closed form at each candidate g and leaves four
residuals in three unknowns -- again one over-identifying restriction and a
J with 1 degree of freedom.  This is a {it:different} estimator, not the
same one solved differently: it weights the first five conditions infinitely
rather than optimally.  Two properties recommend it as a companion rather
than a replacement.  It is one-dimensional, so it is gridded exhaustively
and has never had a local-minimum problem.  And substituting A into its
fifth residual returns m21 - 3 g m12 + 2 g^2 m03, which {it:is} the QME
quadratic, so it shows the third-order route as the part of the GMM that
the fourth moments reweight.  Compare the two: a large gap between
{cmd:e(g_jgmm)} and {cmd:e(g_gmm)} is a symptom of a flat criterion, not of
an error in either.  (Araar 2026c, section 6.3.)

{pstd}
{bf:lsz} {hline 2} Lewbel, Schennach and Zhang (2024){p_end}

{pstd}
A different model, not a variant of this one: the triangular system is the
same but the confounder's loading in the outcome equation is {it:free}, with
nothing tying it to g1.  The parameter spaces are nested; the identifying
assumptions are not.  As {cmd:trigmm} implements it with p(0 1) the system has
5 + 2(k+1) parameters and as many moments, so it is exactly identified and the
estimate is simply the root of the moment vector -- which is why a
deterministic solver reproduces {cmd:trigmm} to ten decimals rather than
approximately.  Read {cmd:e(lsz_crit)} and {cmd:e(lsz_conv)} before
{cmd:e(lsz)}: a criterion that has not reached its root makes the point
meaningless.  See {opt sign()} for the orientation.  (Lewbel, Schennach and
Zhang 2024.)

{pstd}
{bf:lewbel12} {hline 2} heteroskedasticity-generated instruments{p_end}

{pstd}
Builds instruments (X - mean X) * eps2 from the first-stage residual and runs
two-stage least squares with them.  The identifying condition is
Cov(X, eps2^2) != 0 -- heteroskedasticity of the first stage in the exogenous
controls -- and it therefore requires at least one such control; with none the
instrument set is empty and the route is refused.  {helpb freeivmenu} reports
the F that decides whether the condition holds; when it is absent the estimate
is noise with a standard error attached.  This route reproduces {cmd:ivreg2h}
to six decimals on the same data.  (Lewbel 2012.)

{pstd}
{bf:copula} {hline 2} Park and Gupta (2012){p_end}

{pstd}
Adds the control function Phi^-1(F-hat(Y2)) to the outcome equation, where
F-hat is the empirical distribution of the endogenous regressor.  It assumes a
Gaussian copula between Y2 and the structural error, and it needs Y2 to be
{it:non-normal}: if Y2 is Gaussian the control function is collinear with Y2
itself and nothing is identified.  (Park and Gupta 2012.)

{pstd}
{bf:rank} {hline 2} rank control function{p_end}

{pstd}
The same idea with a different transformation, using the ranks of the
residuals rather than a normal quantile.  It rests on non-Gaussianity in the
same way, and is reported beside the copula so that the two can be compared
where both are available.  (Breitung, Mayer and Wied 2024.)

{pstd}
{bf:oster} {hline 2} coefficient stability{p_end}

{pstd}
Not an instrument-free estimator in the sense of the others: it asks how far
the coefficient would move if selection on unobservables were proportional to
selection on the observables already included, and reports the bias-adjusted
value b1 - delta (b0 - b1)(rmax - R1)/(R1 - R0).  Everything hinges on delta
and rmax, which are assumed; see the options above.  It reproduces
{cmd:psacalc} to five decimals.  (Oster 2019.)

{pstd}
{bf:rpiv} and {bf:ape} {hline 2} kept for comparison{p_end}

{pstd}
{cmd:rpiv} purges the residual in two steps and uses the purged variable as an
instrument.  It is {it:inconsistent} -- this is established, not suspected --
and it is retained only because it appears in the earlier literature and
readers will ask.  {cmd:ape} is an adaptive rule that chooses between OLS and
RPIV according to the estimated relative bias.  Neither should be reported as
a preferred estimate.  (Araar 2026a.)

{pstd}
{bf:the two-indicator model} {hline 2} two endogenous regressors, one
confounder{p_end}

{pstd}
Given two endogenous regressors Y2 and Y3 that load on the same latent U, the
third-order cross-moments identify g2, g3 {it:and} the confounder's own
loading a1 in closed form, with no scale-consistency restriction at all.  That
is what makes the restriction testable rather than assumed: the command
reports a1 - (g2 a2 + g3 a3), which is identically zero under scale
consistency and has a standard error here.  Three further statistics come with
it -- the one-factor check R1 - R3 (its ratio verdict is read only when both
third-order cross-moments have z >= 2, since R2 and R3 divide by them), the J
of the sixteen-moment GMM, and what
Y3 would return if it were misused as an instrument for Y2, which it is not,
being correlated with the very thing it would have to purge.  Model B is
selected by the syntax alone: two variables inside the parentheses.  See
{helpb freeivtest}.  (Araar 2026d.)


{marker reading}{...}
{title:How to read the output}

{pstd}
In order, and the command prints them in this order:

{phang}1. {helpb freeivmenu} {it:before} estimating: it says which routes the
data can carry at all.  A weak signal does not make an estimator wrong, it
makes it imprecise.{p_end}

{phang}2. The identified set.  Any estimate outside it, from any route,
implies a negative variance under scale consistency.  {helpb freeivdiag} says
which one.{p_end}

{phang}3. The point estimate and its standard error.  The QME's is
proportional to 1/sqrt(D) and explodes as the discriminant nears zero; this is
not a numerical accident but the shape of the problem.{p_end}

{phang}4. {helpb freeivtest} for what the data say about the assumptions
themselves: whether there is any endogeneity at all, whether the routes agree,
and whether an outside estimate is compatible with the set.{p_end}


{marker whichopt}{...}
{title:Which option acts on which route}

{pstd}
An option is legal for every route but acts on one.  This is the table the
dialog box follows -- it shows the group box for the route chosen and nothing
else -- and it is the table to read before wondering why a value made no
difference.

{synoptset 14 tabbed}{...}
{synopthdr:route}
{synoptline}
{synopt:{opt qbe}}{opt quantile()}{p_end}
{synopt:{opt oster}}{opt delta()} and {opt rmax()}{p_end}
{synopt:{opt lsz}}{opt sign()}{p_end}
{synopt:{opt all}}all four: the table it prints contains the Q-BE, Oster and
LSZ rows, and each option steers its own row{p_end}
{synopt:every other}none.  {cmd:bounds ols qme sce rre hme gmm pgmm lewbel12 copula
rank rpiv ape} take no option of their own{p_end}
{synopt:two indicators}none but {opt level()} and {opt noheader}; the
two-indicator model has a single closed-form route.  The others are refused
rather than ignored{p_end}
{synoptline}
{p2colreset}{...}

{pstd}
{bf:Standard errors.}  Seven routes return a point with no analytic variance:
{cmd:qbe lewbel12 copula rank rpiv ape oster}; {cmd:bounds} returns no point at
all.  No confidence interval is printed for them, and {opt level()} changes
nothing in
{cmd:freeiv}'s own output -- but it is still read by a {cmd:bootstrap:} or
{cmd:jackknife:} prefix, which owns the coefficient table when one is used.
The routes that do carry an analytic standard error are {cmd:ols qme sce rre
hme gmm pgmm lsz}, computed by the delta method on the stacked moments, with the
estimation of the coefficients on X accounted for.

{pstd}
{bf:After estimation.}  {helpb freeivdiag} says what the retained value
implies about the unobservables, whatever produced it; {helpb freeivtest}
reports endogeneity, agreement between routes, and an outside estimate
against the identified set; {helpb freeivmenu} is run {it:before} estimating,
on the same variables; {helpb freeivreport} runs the three into one table.

{pstd}
{bf:Weights.}  Every moment is weighted, and the standard errors come from
the same weighted influence functions.


{marker refused}{...}
{title:What freeiv refuses}

{pstd}
An option that is quietly absorbed is worse than one that is refused, because
the result still looks like a result.  These are rejected with a message
naming the constraint, before any computation:

{phang}{cmd:sign()} other than 1 or -1.  It reaches the LSZ solver as a
direction, where any positive value would have behaved as 1 and any negative
as -1 -- and the typed value would have survived into {cmd:e(lsz_sign)}.{p_end}

{phang}{cmd:rmax()} outside (0, 1].  It is an R-squared.  Omit it for the
default, min(1.3 R1, 1).{p_end}

{phang}{cmd:quantile()} outside (0, 1), and a missing {cmd:delta()}.{p_end}

{phang}A variable in two roles: {it:depvar} also endogenous, {it:depvar} also
a control, an endogenous variable also a control, or the same variable twice
inside the parentheses.{p_end}

{phang}{cmd:method(lewbel12)} or {cmd:method(oster)} with no exogenous
control.  The first builds its instruments from (X - mean X) * eps2 and the
second compares the R-squared with and without the controls; with no X
neither exists.  Under {cmd:method(all)} they are reported as missing and a
note says why.{p_end}

{pstd}
Two situations produce a note rather than an error: {cmd:method()} given with
two endogenous regressors, which does not apply, and {cmd:method(all)} with no
exogenous control.

{pstd}
Numerical failures are reported, not hidden.  A negative discriminant makes
the QME return the vertex and say so; {cmd:e(lsz_conv)} and {cmd:e(lsz_crit)}
say whether the LSZ solver reached its root; {cmd:e(guard)} says whether a
guard of Proposition 1 fired in model B.  A collinear expansion of the
controls is absorbed rather than fatal, the projection being computed with a
generalized inverse.


{title:What the third order cannot tell apart}

{pstd}
One caveat applies to every instrument-free route, not only to these.  In
the single-indicator model the three third-order moments serve three
unknowns (gamma, A, B): the system is just-identified, so any triple
(m03, m12, m21) can be rationalised by a confounder, and a mechanism with no
confounder at all produces the same triple.  Take y1 = a*y2 + b*y2^2 + e with
y2 standard normal, e independent and no latent factor.  Then m03 = 0,
m12 = 2b, m21 = 4ab, and the QME quadratic reduces to -6b*gamma + 4ab = 0:
{cmd:method(qme)} returns 2a/3 for every b other than zero -- a correction of
one third of the OLS slope that does not depend on the size of the quadratic
term, lies inside the interval, has a positive discriminant, and implies a
confounder share of one half with positive variances.  Nothing at order
three can object, because nothing at order three is over-identified.

{pstd}
The fourth order can: on that design with b = 0.1 and n = 20,000 the J of
{cmd:method(pgmm)} is 14 and that of {cmd:method(gmm)} is 5.4, both
rejecting.  But the power of the test falls with b^2 and with n while the
bias does not: at b = 0.02 the J is 0.2 and the estimate is still 0.70.  So
read the z of the discriminant in {helpb freeivmenu} as the strength of a
third-order signal, not as its source: a curvature in E[xi | eps2] is what a
skewed confounder and a non-linear outcome equation both produce, and a
RESET-type check on eps2^2 rejects under either.  The two-indicator model is
the least exposed route, because its loading layer never reads y1.  And a
valid external instrument, where one exists, is the only thing that
identifies gamma with no assumption on the form of the outcome equation:
an instrument-free route complements one, it does not replace it.


{marker examples}{...}
{title:Examples}

{pstd}
Four datasets travel with the package.  {cmd:net get} copies them into the
current directory{p_end}
{phang2}{cmd:. net get freeiv}{p_end}

{pstd}
The order below is deliberate.  Simulated data first, because there the true
coefficient is known and each route can be judged against it.  Then the same
data with one assumption broken, to see which route notices.  Then real data,
where a single endogenous regressor turns out to carry no usable third-order
signal and only the interval survives -- and where a second indicator of the
same confounder then identifies the coefficient after all.


{pstd}{bf:1. The estimator, where the model holds}{p_end}

{pstd}
{cmd:freeiv_sim1.dta} is generated from{p_end}

{phang2}{cmd:u  = (chi2(3) - 3) / sqrt(6)}{space 5}{it:skewed, unit variance}{p_end}
{phang2}{cmd:y2 = 0.30 x + 0.80 u + v2}{space 6}{it:v1, v2, x standard normal}{p_end}
{phang2}{cmd:y1 = 0.60 x + 0.40 y2 + 0.32 u + v1}{p_end}

{pstd}
so the true coefficient is {bf:0.40}, and the confounder's loading in the
outcome equation is 0.40 x 0.80 = 0.32 -- which is scale consistency, holding
here by construction.  {cmd:u} is skewed, so the third-order routes are
available{p_end}

{phang2}{cmd:. use freeiv_sim1, clear}{p_end}
{phang2}{cmd:. freeivmenu y1 x (y2)}{p_end}
{phang2}{cmd:. freeiv y1 x (y2)}{p_end}
{phang2}{cmd:. freeiv y1 x (y2), method(all)}{p_end}

{pstd}
Read the menu first, then every route against the identified interval and
against 0.40.

{pstd}
The two over-identified routes are worth running side by side on these data,
because here they should agree{p_end}

{phang2}{cmd:. freeiv y1 x (y2), method(gmm)}{p_end}
{phang2}{cmd:. freeiv y1 x (y2), method(pgmm)}{p_end}
{phang2}{cmd:. display e(g_jgmm), e(ar_lo), e(ar_hi), e(ar_frac)}{p_end}

{pstd}
{cmd:method(gmm)} prints, under its Hansen J, the range of gamma at which that
J stays within 3.84 of its minimum, and what share of the identified interval
that range covers.  On {cmd:freeiv_sim1} it covers about 64% of the interval:
identification is strong here, so the moments of orders three and four really
do narrow what the second-order algebra gives on its own.  Keep that number in
mind -- example 3 shows the same statistic on real data, where it behaves
quite differently.

{pstd}
Read the share before the point estimate.  It inverts the J at each fixed
gamma instead of inverting a Wald statistic around the optimum, so it stays
valid whether or not the minimum is interior -- which the standard error does
not.


{pstd}{bf:2. The same data at the edge of identification}{p_end}

{pstd}
{cmd:freeiv_sim2.dta} differs in one respect only: {cmd:v2} is drawn from the
same skewed law as {cmd:u} instead of a normal.  Scale consistency still holds
and the truth is still 0.40 -- yet the QME returns nothing at all{p_end}

{phang2}{cmd:. use freeiv_sim2, clear}{p_end}
{phang2}{cmd:. freeivmenu y1 x (y2)}{p_end}
{phang2}{cmd:. freeiv y1 x (y2), method(all)}{p_end}

{pstd}
The reason is exact rather than accidental.  The discriminant of the quadratic
is D = gamma^2 (2A - B)^2, with A = alpha2^3 E[U^3] and B = E[V2^3], so it
vanishes when B = 2A.  A standardised chi2(3) has third moment
24/6^(3/2) = 4/sqrt(6) = 1.6330, so with alpha2 = 0.80 and U and V2 drawn from
that same law, A = 0.8^3 x 1.6330 = 0.8361 and B = 1.6330: 2A - B is 0.0392
rather than zero, and the population discriminant is 0.00025 against 0.447 in
{cmd:freeiv_sim1.dta}, where V2 is normal and B is zero.  This design sits, by
construction, almost exactly on
the surface where the two roots merge, and in any one sample D lands on either
side of zero.

{pstd}
What the output shows is worth more than an estimate.  The QME reports no real
root and falls back on the vertex, saying so; its standard error, proportional
to 1/sqrt(D), is enormous.  {cmd:method(hme)} meanwhile returns a
comfortable-looking number near 0.34 -- while the very assumption it rests on,
E[V2^3] = 0, is violated by construction.  A route that fails loudly is safer
than one that fails quietly, and the menu's B is what tells them apart.


{pstd}{bf:3. Real data: one endogenous regressor gives an interval}{p_end}

{pstd}
{cmd:freeiv_card.dta} is the Card (1995) extract used in the empirical section
of Araar (2026d).  With the controls of that paper{p_end}

{phang2}{cmd:. use freeiv_card, clear}{p_end}
{phang2}{cmd:. freeivmenu lwage exper expersq black south smsa (educ)}{p_end}
{phang2}{cmd:. freeiv lwage exper expersq black south smsa (educ)}{p_end}
{phang2}{cmd:. freeiv lwage exper expersq black south smsa (educ), method(all)}{p_end}

{pstd}
the retained value does fall inside the identified interval and the implied
variances are positive, so nothing is refuted -- but the menu shows the
discriminant indistinguishable from zero, and the QME's standard error is
proportional to 1/sqrt(D).  The honest reading on these data is the interval,
not the point.  This is the ordinary case with one endogenous regressor, and
it is why {helpb freeivmenu} exists.

{pstd}
Four independent diagnostics say it, which is the point of running
{cmd:method(all)} here.  The discriminant is indistinguishable from zero, so
the third-order route carries almost no signal.  The region where the joint
GMM's J stays within 3.84 of its minimum covers the identified interval
{it:entirely} -- every value the assumption-free bounds admit is a value the
nine moment conditions cannot reject -- against 64% on {cmd:freeiv_sim1}, so
the fourth-order moments add nothing here at all.  That GMM's minimum sits on the {it:boundary}, with the implied
confounder variance at its floor, and the command prints no standard error
because none is valid there.  And the LSZ solver does not converge on the same
data, which is reported rather than hidden behind the tightest-looking
interval in the table.  Four routes, four ways of saying that these data do
not identify a point.

{pstd}
Dropping a control changes the answer, and the command says why rather than
leaving it to be discovered{p_end}

{phang2}{cmd:. freeiv lwage exper expersq (educ)}{p_end}

{pstd}
Here the retained value leaves the interval and an implied variance turns
negative.  Both are printed: the interval is exactly the region where those
variances are non-negative, so a value outside it is not a large estimate but
an impossible one.


{pstd}{bf:4. The same data, with a second indicator of the confounder}{p_end}

{pstd}
Mother's education loads on the same family background that confounds
schooling, and it is itself endogenous in the wage equation -- it is not an
instrument, and using it as one is the comparison the command prints{p_end}

{phang2}{cmd:. freeiv lwage exper expersq black south smsa (educ motheduc)}{p_end}
{phang2}{cmd:. freeivreport lwage exper expersq black south smsa (educ motheduc)}{p_end}

{pstd}
Now the coefficient is identified as a point, with a standard error, and the
confounder's own loading a1 is free and estimated rather than restricted.  The
scale-consistency test that model A cannot perform is reported here, together
with what {cmd:motheduc} would have returned if it had been treated as an
instrument.  {cmd:fatheduc} can be used in its place.

{pstd}
These two lines reproduce Table 6b of Araar (2026d) on the data that paper
used: g2 = 0.066, g3 = -0.009, a1 = 0.067 and an IV gap of 0.102 with
{cmd:motheduc}, and g2 = 0.035 with {cmd:fatheduc}, the departure from scale
consistency being 0.032 and 0.563 respectively.  The agreement is checked
rather than asserted -- {cmd:freeiv_test10.do} in the replication material
runs the paper's own script and compares value by value, the guards included.

{pstd}
Two proxies of the same file are {it:refused} by the guards of Proposition 1,
which is the intended behaviour and worth seeing{p_end}

{phang2}{cmd:. freeiv lwage exper expersq black south smsa (educ IQ)}{p_end}
{phang2}{cmd:. freeiv lwage exper expersq black south smsa (educ KWW)}{p_end}

{pstd}
The guards fire when the third-order moments imply a negative variance or
incoherent signs -- conditions no one-factor structure can produce -- so the
closed form is not reported at all.  {cmd:e(guard)} records which fired: 1 for
opposite signs of E[e2^2 e3] and E[e2 e3^2], 2 for an incoherent product and
ratio, 3 for a negative implied variance.  {cmd:KWW} is refused by the first
and {cmd:IQ} by the third, which is Table 6a of Araar (2026d), reason by
reason.


{pstd}{bf:5. Two indicators with a known truth}{p_end}

{pstd}
{cmd:freeiv_proxy.dta} is simulated, n = 5000, from{p_end}

{phang2}{cmd:y2 = 0.3 x + 0.8 u + v2}{p_end}
{phang2}{cmd:y3 = 0.2 x + 0.6 u + v3}{p_end}
{phang2}{cmd:y1 = 0.6 x + 0.40 y2 + 0.20 y3 + 0.10 u + v1}{p_end}

{pstd}
with {cmd:v1 v2 v3 x} standard normal and {cmd:u} a centred chi-squared
rescaled to unit variance and skewness 1 -- the skewness is what the
third-order moments live on, and a symmetric {cmd:u} would leave the model
unidentified.  So g2 = 0.40, g3 = 0.20 and a1 = 0.10 are the truth.  Scale
consistency does {it:not} hold here, deliberately: it would require
a1 = g2 a2 + g3 a3 = 0.44 against 0.10, which is why the test rejects, and
should{p_end}

{phang2}{cmd:. use freeiv_proxy, clear}{p_end}
{phang2}{cmd:. freeiv y1 x (y2 y3)}{p_end}
{phang2}{cmd:. freeivtest}{p_end}


{pstd}{bf:Other things}{p_end}

{pstd}Inference for a route with no analytic standard error{p_end}
{phang2}{cmd:. use freeiv_sim1, clear}{p_end}
{phang2}{cmd:. bootstrap, reps(500): freeiv y1 x (y2), method(qbe)}{p_end}

{pstd}Oster's answer is a function of what you assume, so report a range{p_end}
{phang2}{cmd:. use freeiv_card, clear}{p_end}
{phang2}{cmd:. freeiv lwage exper expersq black south smsa (educ), method(oster) rmax(0.5)}{p_end}
{phang2}{cmd:. freeiv lwage exper expersq black south smsa (educ), method(oster) rmax(0.8) delta(2)}{p_end}

{pstd}Factor variables among the controls{p_end}
{phang2}{cmd:. sysuse nlsw88, clear}{p_end}
{phang2}{cmd:. generate lwage = ln(wage)}{p_end}
{phang2}{cmd:. freeiv lwage grade age i.race (tenure)}{p_end}

{pstd}Confronting an estimate obtained elsewhere with the identified set{p_end}
{phang2}{cmd:. use freeiv_card, clear}{p_end}
{phang2}{cmd:. freeiv lwage exper expersq black south smsa (educ)}{p_end}
{phang2}{cmd:. freeivtest, gamma(0.132) segamma(0.049)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:freeiv} is an {cmd:eclass} command.  The full list is long; {cmd:ereturn}
{cmd:list} shows it.  The ones most often wanted:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(gamma)}}the retained estimate{p_end}
{synopt:{cmd:e(se)}}its standard error{p_end}
{synopt:{cmd:e(lo)}, {cmd:e(hi)}}the identified set{p_end}
{synopt:{cmd:e(gt)}}the OLS slope, which is also {cmd:e(hi)}{p_end}
{synopt:{cmd:e(qme)}, {cmd:e(sce)}, {cmd:e(rre)}, {cmd:e(hme)}, {cmd:e(qbe)}}the closed forms{p_end}
{synopt:{cmd:e(g_jgmm)}, {cmd:e(J9)}, {cmd:e(p_J9)}}the joint nine-moment GMM and its over-identification test{p_end}
{synopt:{cmd:e(ar_lo)}, {cmd:e(ar_hi)}, {cmd:e(ar_frac)}}the region where the profiled J is within 3.84 of its minimum, and its share of the identified interval{p_end}
{synopt:{cmd:e(j_bound)}, {cmd:e(j_weak)}}the minimum is on the boundary; theta is below 0.05{p_end}
{synopt:{cmd:e(g_gmm)}, {cmd:e(J)}, {cmd:e(p_J)}}the profiled GMM and its over-identification test{p_end}
{synopt:{cmd:e(idfac)}, {cmd:e(z_idfac)}}the factor that says whether that J can test anything{p_end}
{synopt:{cmd:e(lsz)}, {cmd:e(lsz_crit)}, {cmd:e(lsz_conv)}}LSZ and its convergence{p_end}
{synopt:{cmd:e(lewbel12)}, {cmd:e(copula)}, {cmd:e(g_rank)}, {cmd:e(oster)}}the literature routes{p_end}
{synopt:{cmd:e(disc)}, {cmd:e(disc_z)}}the discriminant and its z{p_end}
{synopt:{cmd:e(theta)}, {cmd:e(sV1)}, {cmd:e(sV2)}, {cmd:e(k)}, {cmd:e(kstar)}}the implied nuisances{p_end}
{synopt:{cmd:e(A)}, {cmd:e(B)}, {cmd:e(mu)}}the third-moment decomposition{p_end}

{pstd}In model B{p_end}
{synopt:{cmd:e(g2)}, {cmd:e(g3)}, {cmd:e(a1)}}the structural coefficients and the free loading{p_end}
{synopt:{cmd:e(sc)}}g2 a2 + g3 a3, which is what model A would call a1{p_end}
{synopt:{cmd:e(sc_d)}, {cmd:e(z_sc)}}the scale-consistency test{p_end}
{synopt:{cmd:e(of_d)}, {cmd:e(z_of)}}the one-factor test{p_end}
{synopt:{cmd:e(z_m223)}, {cmd:e(z_m233)}}the z of the two third-order cross-moments, mean over sd/sqrt(n) as {helpb freeivmenu} prints them; the one-factor ratio |R3/R1 - 1| is read only when both exceed 2{p_end}
{synopt:{cmd:e(q_J)}, {cmd:e(q_pJ)}}the 16-moment GMM and its J with 4 df{p_end}
{synopt:{cmd:e(guard)}}0, or the guard of Proposition 1 that fired{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}, {cmd:e(V)}}the posted estimate and its variance{p_end}
{synopt:{cmd:e(moments)}}every moment and derived quantity, named{p_end}
{synopt:{cmd:e(Vmom)}, {cmd:e(grad)}}the moment covariance and the gradients, used by {helpb freeivtest}{p_end}
{p2colreset}{...}


{marker references}{...}
{title:Documentation}

{pstd}
{bf:freeiv_paper.pdf} presents the models, the command and worked examples at
article length, in the format of a Stata Journal article.  It comes down
with {cmd:net get freeiv}, next to the datasets, and is kept at
{browse "https://github.com/aabbdd12/freeiv"}.  The output it prints is
that of the installed version.


{title:References}

{phang}
Araar, A. 2026a.  Correcting endogeneity without external instruments: four
closed-form estimators.  Zenodo 10.5281/zenodo.22067980.

{phang}
Araar, A. 2026b.  The quadratic moment estimator.
Zenodo 10.5281/zenodo.22068143.

{phang}
Araar, A. 2026c.  Instrument-free estimation under linear scale-consistency.
Zenodo 10.5281/zenodo.22753299.

{phang}
Araar, A. 2026d.  Two indicators of one latent confounder: closed-form
identification of the triangular model with a free proxy effect.
Zenodo 10.5281/zenodo.22764001.

{phang}
Card, D. 1995.  Using geographic variation in college proximity to estimate
the return to schooling.  In {it:Aspects of Labour Market Behaviour}, ed.
L. N. Christofides et al., 201-222.  Toronto: University of Toronto Press.
The extract shipped as {cmd:freeiv_card.dta} comes from the teaching datasets
distributed with Wooldridge's {it:Introductory Econometrics}.

{phang}
Breitung, J., A. Mayer and D. Wied. 2024.  Asymptotic properties of
endogeneity corrections using nonlinear transformations.
{it:The Econometrics Journal} 27: 362-383.

{phang}
Lewbel, A. 2012.  Using heteroscedasticity to identify and estimate
mismeasured and endogenous regressor models.  {it:Journal of Business and
Economic Statistics} 30: 67-80.

{phang}
Lewbel, A., S. M. Schennach and L. Zhang. 2024.  Identification of a
triangular two equation system without instruments.  {it:Journal of Business
and Economic Statistics} 42: 14-25.

{phang}
Oster, E. 2019.  Unobservable selection and coefficient stability: theory and
evidence.  {it:Journal of Business and Economic Statistics} 37: 187-204.

{phang}
Park, S. and S. Gupta. 2012.  Handling endogenous regressors by joint
estimation using copulas.  {it:Marketing Science} 31: 567-586.


{title:Author}

{pstd}
Abdelkrim Araar, Universite Laval and PEP{break}
{browse "mailto:aabd@ecn.ulaval.ca":aabd@ecn.ulaval.ca}


{title:Also see}

{psee}
Online: {helpb freeivmenu}, {helpb freeivdiag}, {helpb freeivtest},
{helpb ivregress}
{p_end}
