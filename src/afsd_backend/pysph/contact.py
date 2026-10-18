"""Host adapter of pinned PySPH wall collision, not a new rigid-body solver.

Normal penalty: upstream RigidBodyWallCollision, restitution=1 (no hard-coded
effective-mass damping). Tangential trial is memoryless, impulse-limited sliding;
this is explicitly NOT a validated static-stick history or P44 Eq3-160 law.
"""
import contextlib
import io
import numpy as np
from pysph.sph.rigid_body import RigidBodyWallCollision
from .geometry import surfaces, assert_outside_tool


class WallContact:
    def __init__(self, stiffness, friction):
        with contextlib.redirect_stdout(io.StringIO()):
            self.law = RigidBodyWallCollision('solid', ['wall'], kn=stiffness, mu=friction, en=1.)

    def force(self, gap, normal, relative_velocity, radius, mass, dt):
        n = np.asarray(normal, float); v = np.asarray(relative_velocity, float)
        if abs(np.linalg.norm(n)-1) > 1e-12:
            raise ValueError('contact normal must be unit length')
        overlap = max(radius-gap, 0.)
        tangent = v-np.dot(v, n)*n
        speed = np.linalg.norm(tangent)
        fn = self.law.kn*overlap
        # Documented regularization: cannot reverse tangential slip in one kick.
        ft = min(self.law.mu*fn, mass*speed/dt)
        displacement = tangent/speed*ft/self.law.kt if speed > 0 else np.zeros(3)
        force = [np.zeros(1) for _ in range(3)]
        tang_v = [np.zeros(1) for _ in range(3)]
        ones = np.ones(1)
        self.law.loop(0, *force, ones, np.array([mass]), np.array([radius]),
                      *[np.array([d]) for d in displacement], *tang_v,
                      0, gap*n, abs(gap), gap*gap, v, *[np.array([d]) for d in n])
        f = np.array([q[0] for q in force])
        friction_power = -np.dot(f-fn*n, tangent)
        if friction_power < -1e-12:
            raise ValueError('negative friction dissipation')
        return f, friction_power, .5*self.law.kn*overlap**2, overlap

    def conductance_row(self, x, t, p):
        """Per-particle wall thermal row [W/K], to include in explicit dt bound."""
        assert_outside_tool(x, t, p)
        if p['contact_conductance'] < 0:
            raise ValueError('negative contact conductance')
        row = np.zeros(len(x))
        for ids, gaps, _, _, _ in surfaces(x, t, p):
            active = ids[gaps < p['dx']/2]
            np.add.at(row, active, p['contact_conductance']*p['dx']**2)
        return row

    def evaluate(self, x, v, mass, temperature, t, dt, p, heat_velocity=None):
        """Force at v; friction work at an independently supplied kick midpoint.

        With heat_velocity=None the returned power is instantaneous, as in force.
        A time integrator must pass (v_before+v_after)/2 AFTER all impulses and
        actuator constraints, not integrate old-velocity power as if exact heat.
        """
        assert_outside_tool(x, t, p)
        if heat_velocity is not None:
            heat_velocity = np.asarray(heat_velocity)
            if heat_velocity.shape != v.shape or not np.isfinite(heat_velocity).all():
                raise ValueError('invalid contact heat midpoint velocity')
        force = np.zeros_like(x); heat = np.zeros(len(x))
        ledger = dict(wall_work_rate=0., wall_heat_rate=0., friction_power=0., contact_heat_rate=0., contact_potential=0., max_penetration=0., contacts=0)
        for ids, gaps, normals, wall_v, wall_t in surfaces(x, t, p):
            for a, idx in enumerate(ids):
                f, dissipation, potential, penetration = self.force(gaps[a], normals[a], v[idx]-wall_v[a], p['dx']/2, mass[idx], dt)
                if heat_velocity is not None:
                    normal = normals[a]
                    slip_mid = heat_velocity[idx]-wall_v[a]
                    tangent_mid = slip_mid-np.dot(slip_mid, normal)*normal
                    friction_force = f-self.law.kn*penetration*normal
                    dissipation = -np.dot(friction_force, tangent_mid)
                    if dissipation < -1e-12:
                        raise ValueError('negative midpoint friction dissipation; coupled kick reverses contact slip')
                force[idx] += f
                # Contact conductance acts only on actually compressed contacts.
                flux = p['contact_conductance']*p['dx']**2*(wall_t-temperature[idx]) if penetration > 0 else 0.
                heat[idx] += p['heat_partition']*dissipation+flux
                ledger['wall_work_rate'] += np.dot(f, wall_v[a])
                ledger['wall_heat_rate'] += (1-p['heat_partition'])*dissipation-flux
                ledger['friction_power'] += dissipation
                ledger['contact_heat_rate'] += flux
                ledger['contact_potential'] += potential
                ledger['max_penetration'] = max(ledger['max_penetration'], penetration)
                ledger['contacts'] += int(penetration > 0)
        return force, heat, ledger
