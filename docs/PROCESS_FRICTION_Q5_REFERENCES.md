# Q5 directional friction: independent manufactured references

These references are derived here from static force balance and piecewise Coulomb
inequalities, not by evaluating the implementation under test. All quantities are
SI. They are manufactured benchmark constants, not AA6061 or paper parameters.
The contact law has one frictional tangent direction; the other tangent is free.
This is not certification of isotropic three-dimensional Coulomb contact, rotating
features, thermal partition, inertia, or the earlier dynamic AFSD contact method.

## Sign and work convention

Let `u` be the material contact displacement, `w` the rigid-wall displacement,
`r = u - w` the relative displacement, and `f` the tangential force **on material**.
For normal compression `lambda >= 0`, the Coulomb capacity is `C = mu*lambda`.
Stick requires `delta r = 0` and `-C <= f <= C`. Sliding requires
`f = -C*sign(delta r)`. A mode value `+1/-1` denotes the sign of relative slip,
not the sign of force or wall motion.

On an accepted piecewise-smooth subinterval:

    W_wall     = f_mid * delta w
    W_material = f_mid * delta u
    Q_friction = -f_mid * delta r = W_wall - W_material

There is no absolute value, positive-part operation, residual heat, or inference
of correctness from Q alone. Event positions, final displacements, mode sequence,
and signed forces are separate reference observables. Stick can have nonzero
force, wall work, and material work while Q is exactly zero. `Q_friction` is a
mechanical dissipated-work diagnostic; temperature is not solved in this stage.

## A. Grounded elastic spring and reversing moving wall

Take `k = 1 N/m`, `C = 1 N`, initial `u = w = 0`, and prescribe the wall path
`w: 0 -> 3 -> 0 m`. Static balance gives `f = k*u`; there is no particle inertia.

1. Increasing w from 0 to 1: stick, `u=w`, `f=w`.
2. Increasing w from 1 to 3: negative relative slip, `u=1`, `f=1`.
3. Decreasing w from 3 to 1: stick, `u=w-2`, `f=w-2`.
4. Decreasing w from 1 to 0: positive relative slip, `u=-1`, `f=-1`.

Literal event endpoints and per-interval integrals:

| endpoint | w m | u m | f N | W_wall J | W_material J | Q J |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| initial | 0 | 0 | 0 | - | - | - |
| first slip threshold | 1 | 1 | 1 | .5 | .5 | 0 |
| wall reversal | 3 | 1 | 1 | 2 | 0 | 2 |
| reverse slip threshold | 1 | -1 | -1 | 0 | 0 | 0 |
| final | 0 | -1 | -1 | 1 | 0 | 1 |

Totals: `Q=3 J`, `W_wall=3.5 J`, `W_material=delta U=.5 J`.
After the wall returns to zero, the material remains at `u=-1 m`.
At wall reversal the force does **not** immediately reverse: the contact sticks
while f decreases from +1 through 0 to -1. Relative slip is zero throughout that
force reversal. A scheme that flips f solely with wall velocity is incorrect.

For arbitrary positive k,C, set `a=C/k`, wall path `0 -> 3a -> 0`, event values
`u=[0,a,a,-a,-a]`. Work scales with `C^2/k`: Q is 3 times this scale, W_wall 3.5
times, and delta U .5 times. These are exact piecewise-static values.

For the extended two-reversal wall program `0 -> +3a -> -3a -> 0`, literal
event arrays are `w=[0,a,3a,a,-3a,-a,0]`, `u=[0,a,a,-a,-a,a,a]`, and
`f=[0,C,C,-C,-C,C,C]`. Sliding dissipation by leg is `[2,4,1]*C^2/k`, giving
`Q=7*C^2/k`, `W_wall=7.5*C^2/k`, and final `U=.5*C^2/k` with `u=+a`.
The two-reversal case has a different final state and must not be checked
against the one-reversal `u=-a`, `Q=3*C^2/k` values.

This reference does not replace the earlier inertial test `m=1 kg, u0=1 m/s,
F_other=-2 N, C=1 N, dt=1 s`, whose independent values remain
`t_stop=1/3 s`, `u_final=-2/3 m/s`, `Q=7/18 J`, `W_other=1/9 J`, `delta K=-5/18 J`.
No Q5 static PASS can change the status of that dynamic test.

## B. Static wall and driven spring support

Let `k=C=1`, static wall `w=0`, and spring base displacement `b:0 -> 3 m`.
Static balance is `f=k*(u-b)`. For `0<=b<=1`, stick gives `u=0`, `f=-b`.
For `1<=b<=3`, positive relative slip gives `u=b-1`, `f=-1`.
The final values are `u=2`, `f=-1`, spring energy `.5 J`.

Exact integrals are `W_wall=0`, `W_material=-2 J`, `Q=2 J`, and prescribed-base
work `W_base=2.5 J`. Thus `W_base + W_wall = delta U + Q`.
This detects confusing laboratory material work with relative frictional work.

## C. Translating zero-slip preload

For `k=C=1`, choose a preloaded spring extension `.5 m`: `u=w=b+.5` and `f=.5 N`.
Translate wall, material, and spring base by the same `2 m`. Then
`W_wall=W_material=1 J`, `W_base=-1 J`, `Q=0`, and `delta U=0`.
This is a work-identity reference, not a claim that the fixed-base 3D block can
undergo such a rigid translation. It is also not a finite-rotation test.

## D. Two contacts sharing a free elastic degree of freedom

Use the positive-definite elastic stiffness (N/m)

    K = [[ 2,  0, -1],
         [ 0,  2, -1],
         [-1, -1,  3]]

Contacts act on `u1,u2`; the third displacement z is free. Both capacities are
`1 N`, and both walls have the same displacement w. The zero-force equation gives
`z=(u1+u2)/3`. Eliminating z gives the exact contact stiffness

    K_eff = [[5/3, -1/3], [-1/3, 5/3]]

For equal contact displacements the stiffness per contact is `4/3 N/m`, not the
isolated diagonal value 2. Define `a=.75 m` and prescribe `w:0 -> 2.25 -> 0 m`.
The common u follows reference A with this a. Literal event arrays are

    w       = [0, .75, 2.25, .75, 0]
    u1 = u2 = [0, .75,  .75, -.75, -.75]
    z       = [0,  .5,   .5,  -.5,  -.5]
    f1 = f2 = [0,   1,    1,   -1,   -1]

Each wall has `Q=2.25 J`, `W_wall=2.625 J`, and `W_material=.375 J`.
Totals are `Q=4.5 J`, `W_wall=5.25 J`, and final elastic `U=.75 J`.
Independent per-contact diagonal treatment would give the wrong onset `.5 m`;
this reference specifically tests coupling through a shared free DOF.

## E. Coupled normal reaction and tangential friction capacity

An optional stronger algebraic reference uses

    K = [[2, .5], [.5, 1]] N/m
    u_x = d = -1 m, held fixed; mu = .5

The force on material at the contact is `(-lambda, tau) = K*u`, hence
`lambda=2-.5*u_y` and `tau=-.5+u_y`. Starting at `u_y=0`, sticking requires
`lambda=2 N`, `tau=-.5 N`, which is feasible. Increasing wall displacement with
stick `u_y=w` reaches positive friction at

    u_y = w = 6/5 m; lambda=7/5 N; tau=7/10 N.

The opposite slip threshold is

    u_y = -2/3 m; lambda=7/3 N; tau=-7/6 N.

The different thresholds follow from coupled normal reaction, not fixed initial
capacity. Normal wall work is zero because d is held, even while lambda changes.
This algebraic reference has an initial tangential prestress; its preload must
not be silently confused with the zero-state spring in A.

## Numerical acceptance and use with the 3D block

For the unit-scale algebraic A/D/E references, absolute position/force and signed
work tolerance `1e-10 m`, `1e-10 N`, `1e-10 J` is declared before their execution.
It allows linear-solve roundoff while remaining many orders below the expected
slip and work. Event-position tolerance is `1e-10 m`. A reference has not passed
until actually executed; the execution record must name which sections ran.
The independent static fixture has monotone parameter-time intervals as long as
6 s and permits 48 bisections to reach the same `1e-12 s` event bracket width.
The actual block retains its declared 40 bisections for intervals within its
.2 s program. This is a bounded root-search iteration budget, not physical
mass/time scaling or a relaxation of the acceptance width.

The 3D companion retains its separately frozen force/gap/work checks. An
X-minus clamped cube with other free faces is not exactly a homogeneous simple
shear spring: end constraints, warping, and the actual preloaded Jaumann update
may couple shear and normal force. Therefore `k=G*A/L` must not be asserted as an
exact reference for that entire body. These algebraic references test the shared
contact algorithm; actual-body checks separately verify equilibrium, KKT,
event order, independent work readback, passive labels, and rejected-state
purity. They do not prove spatial convergence or remove the existing closed-
cycle constitutive and Q4 strict plastic-tail FAIL evidence.
