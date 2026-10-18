# Q7 boundary-onset independent references (2026-10-02)

Scope: manufactured fixed-normal elastic tangent springs, testing only the new
boundary eligibility/event entry around the unchanged coupled cone rows. These
are not an AFSD, plasticity, thermal or continuum-space validation. Old Q7/Q7R
failures and the old onset RED are retained. All values below were derived before
running the new tests; no tested contact function produces reference answers.

## Units and tolerances

Use SI: tangent stiffness K in N/m, held normal N=2 N per contact, mu=.5,
capacity C=1 N, manufactured time in seconds. E=(.6,.8) is a unit direction in
the actual two-dimensional tangent plane. The force/slip/event limits are the
unchanged benchmark values: 9e-5 N, 1e-11 m, 1e-12 s; 40 bisections, 12 events,
16 active sets, probe fraction1e-4. Root fixture solves to1e-13 N, with32 bounded
Newton iterations; this is a small linear reference assembly, not a changed
27-particle root budget. Assertion atol1e-10 in the stated SI quantity and rtol0
allows floating-point/root and event-bracket error while easily resolving the
9e-5 N classification issue. Raw Q is never abs/max/clipped.

The shared old ElasticVectorOwner fixture assembles linear springs and invokes
the actual cone_response/cone_rows under test. New assertions use only the
formulas/literals here. Added guarded commit checks admissibility, raw work and
the existing no-unsplit-opposite-sliding rule; it does not create an independent
contact solver. Small assembly tests do not replace actual Q7 commit testing.

## Exact straight-path continuous references

1. Immediate outward onset: K=I N/m, initial u=0, tau=E N (prestress supplied
   by spring base=-E m), wall=.2*t*E m. For t=1, u=0, tau=E, Q=.2 J,
   Wwall=.2 J, Wmaterial=0. The mode probe increases the forced-stick force by
   only2e-5 N, below9e-5 N; the old path falsely retains forced stick and cannot
   bracket a strict feasible lower state. A dedicated expectedFailure retains
   this exact old RED instead of deleting it or claiming it passed.
2. Initial inward motion: same preload, wall=-.2*t*E. At t=1, u=-.2E m,
   tau=.8E N, Q=0, Wwall=Wmaterial=-.18 J. Boundary eligibility must not
   prescribe sliding: a complementarity row can return zero slip/stick.
3. True interior event: initial tau=.5E, wall=tE. The event is t=.5 s;
   final u=.5E m, tau=E N, Q=.5 J, Wwall=.875 J, Wmaterial=.375 J.
   An interior contact cannot be released merely because the target is overcone.
4. Coupled two-contact activation: K=[[2,-1],[-1,2]] tensor I N/m,
   initial tau=(1,.8) tensor E N, initial u=0, bases=-(14/15,13/15) tensor E m.
   Walls=(.4t,.2t) tensor E m. Contact1 is initially at capacity. While contact2
   sticks, tau1=1 gives u1=u2/2 and tau2=.8+1.5u2. Therefore its event is
   t=2/3, not the diagonal-contact prediction1/2. Final u=(1/15,2/15) tensor E,
   both tau=E; Q1=1/3 J, Q2=1/15 J, Qtotal=2/5 J, Wwall=44/75 J and
   Wmaterial=DeltaEs=14/75 J. This tests the same coupled stiffness, not separate
   scalar masses. Partitioning straight paths must preserve these totals.
5. Shared unconstrained vector DOFs, mixed slip/stick: use the independently
   specified K=[[2,0,-1],[0,2,-1],[-1,-1,3]] tensor I N/m. The first two vector
   nodes contact walls; the third has zero external force and is solved, not
   prescribed. Initial u=0 and base=-(.75,.75,.5) tensor E m give tau1=tau2=E
   and zero free-node force. Walls=(+.3t,-.3t) tensor E. For t=1, contact1 slides
   with tau1=E; contact2 sticks with u2=-.3E. The equilibrium equations
   2u1-u3=0 and -u1-u2+3u3=0 give u1=-.06E and u3=-.12E, hence tau2=.52E.
   Slip=(-.36E,0), Q=(.36,0) J, Wwall=.072 J, Wmaterial=DeltaEs=-.288 J.
   Both rows may remain eligible; eligibility is not physical sliding. This
   reference was added after the first9 tests PASS, before running this newcase.

## Signed initial roundoff-scale-to-tolerance perturbations

Repeat case1 with initial tau=(1+delta)E, delta=+4e-5 or -4e-5 N. These are
explicit finite-force-tolerance diagnostic states, not newly accepted physical
constitutive inputs. No projection of initial force is permitted. The endpoint
equations give tau=E, u=-delta E, slip=-(.2+delta)E. The actual *discrete midpoint*
ledger, derived independently, is

Q=(1+delta/2)(.2+delta), Wwall=.2+.1delta,
DeltaEs=-delta-delta^2/2.

For +4e-5: Q=.2000440008 J, Wwall=.200004 J, DeltaEs=-.0000400008 J.
For -4e-5: Q=.1999560008 J, Wwall=.199996 J, DeltaEs=.0000399992 J.

These values do NOT claim exact continuous work. The positive-offset initial
state has no exactly admissible Coulomb continuum solution. For the negative
offset, exact stick-then-slide continuous Q=.2+delta; the unsplit midpoint value
differs by delta/2*(.2+delta), about-4e-6 J here. Tolerance-based onset eligibility
therefore introduces a finite near-boundary temporal error even when the discrete
work identity closes. This limitation must not be hidden by ledger PASS.

## Transverse endpoint and remaining boundary

With K=I, initial tau=(1,0), wall increment(.2,1.6) m, the *discrete endpoint*
cone equation has independent root tau=(.6,.8) N: initialtau+dw=(1.2,1.6),
radius2; u=(-.4,.8), slip=(-.6,-.8). Midpoint Q=.8 J, Wwall=.8 J,
Wmaterial=0. This checks the two-tangent direction is not frozen to initial tau.
It is NOT an exact continuously rotating-slip trajectory or its integral.
Arbitrary curved path convergence remains NOT_RUN.

A directly committed full reversal with final force opposite the initial force
must still fail the existing unsplit-transition commit guard. For the resolved
straight reversal wall=-3tE from initial tau=E, the actual tiny mode probe moves
inward by3e-4 N, larger than9e-5 N. The onset driver must therefore retain stick,
locate the opposite boundary at t=2/3, and only then slide: final u=-2E,
tau=-E, Q=1 J, Wwall=1 J, Wmaterial=0. A probe that cannot resolve the inward
excursion is not covered by this continuous-reference assertion. Also reject an
initial overcone state beyond9e-5 N: do not normalize it silently. Such negative
guard tests are expected rejections, not successful physical paths.

All spring references keep strictly positive normal force. They do not certify
zero-capacity startup; the proposed seed excludes old capacity<=9e-5 N, and any
claim about the actual Q7 preload prefix requires a separate actual-body test.
