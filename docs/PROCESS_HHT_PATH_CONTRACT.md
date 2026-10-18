# HHT dense path and contact work — before tests

2026-10-03, candidate-only numerical interface. Original sources unchanged.
Weight a=.9, gamma=.6,beta=.3025 (a=1 midpoint control); d=beta-gamma/2.
Given actual Newmark endpoint states, choose the compatible quadratic path
v(s)=v0+s(v1-v0)+6h*d*(a1-a0)*s*(1-s), 0<=s<=1. Its integral equals x1-x0.
It is a declared temporal interpolation, NOT an imposed free-material trajectory.
Force-driven endpoint solving belongs to the joint caller.

For force linear between f0,f1, two-point Gauss exactly integrates degree-three
force dot velocity. Store each actual quadrature impulse and work velocity;
signed friction Qq=-Jq dot relative_vq, not abs/clipped. Also evaluate maximum
cubic friction power at both endpoints and ALL interior stationary points.
Nonnegative integrated Q alone is insufficient. Interior minimum slip is found
from the quartic squared-slip polynomial derivative; zero events require
splitting/restart by caller, never velocity snapping. Check original cone at
both endpoints; turning slip is not automatically exact Coulomb inside a step.

Contact directions fixed within this interface. Moving/rotating geometry must
be subdivided or separately verified by a caller; fixed normals are not general
AFSD finite-edge contact. Wall velocity affine; force/slip 3D. Whole-path power
bound256eps*force_scale*relative_velocity_scale; zero-speed256eps*velocity_scale.
These are floating-arithmetic certificates, not adjustable physical tolerances.
Explicit raw signed force-filter impulse and kinetic work terms are returned;
their sign is not assumed and they never feed thermal sources.

Independent preregistered tests: exact constant acceleration/displacement;
nonzero-a quadratic path integration; exact linear-force work (hand polynomial);
concealed two-zero reversal v(s)=1-6s+6s² whose integral Q=0 yet max friction
power=.5 W (must reject); 1D event-split reversal reference7/18 J and -2/3 m/s;
wall boost /3D tangential projection; kinematic/impulse-work certificates and
endpoint capacity rejection. Scalar tolerances1e-12 (SI fixture scales), polynomial
event positions1e-10. All runs external<=60s; first failure retained, max2 repairs.
No physical AFSD run or production permission follows from component PASS.
