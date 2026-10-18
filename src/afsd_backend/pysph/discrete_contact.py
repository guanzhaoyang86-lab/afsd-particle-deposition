"""Finite-interval normal contact with a literal discrete-gradient force.

For EACH finite patch, -f_bar . (y1-y0) = U1-U0.  This is a new contact
approximation, not a relabelling of endpoint work.  The force returned here
must enter mechanical equilibrium.  y is tool-body position; base positions
remain in the laboratory.  The inherited swept-solid certificate still uses
linear laboratory particle paths and shortest-geodesic tool rotation.

The effective Coulomb normal is the discrete force direction/magnitude, not
the endpoint distance gradient.  Tangential springs retain backward Euler:
algorithmic dissipation and projection/release energy remain unassigned and
must remain in the process energy gate.  Only plastic sliding work is heat.

The discrete-gradient identity follows the Gonzalez construction.  Its use
on our sum of finite-patch potentials is a project-specific approximation;
it does not transfer an upstream method's full stability/accuracy theorem.
"""
import numpy as np

from .quasistatic_contact import (QuasistaticContact, ContactGeometryError,
                                 FEATURES, _pose)


CLASSIFICATION = 'HOT_START_QUASISTATIC_DISCRETE_NORMAL_BE_FRICTION_CONTACT_APPROXIMATE_COMPONENT'
CHAIN_RULE_CLASSIFICATION = 'HOT_START_QUASISTATIC_CHAIN_RULE_DISCRETE_NORMAL_BE_FRICTION_CONTACT_APPROXIMATE_COMPONENT'


class DiscreteContact(QuasistaticContact):
    """Same SI constructor/evaluate API as QuasistaticContact; pure candidates.

    DG forces whose component along the endpoint geometric outward normal is
    negative beyond roundoff are rejected, not turned positive with abs().
    A body-coordinate midpoint on a nonsmooth/inside branch is also rejected.
    These are finite-step admissibility restrictions, not proofs that every
    swept-admissible path admits this discretization.
    """

    classification = CLASSIFICATION
    normal_model = 'GONZALEZ_SUM_OF_FINITE_PATCH_POTENTIALS_DUPLICATE_LIP_CONTRIBUTIONS_RETAINED'

    def _midpoint_geometry(self, local0, local1, x0, x1):
        middle = .5*(local0+local1)
        gaps = np.full((len(middle), 8), np.inf)
        normals = np.zeros((len(middle), 8, 3))
        for i, point in enumerate(middle):
            gaps[i, :7], normals[i, :7], _ = self._patches(point)
        if self.base_z is not None:
            gaps[:, 7] = .5*(x0[:, 2]+x1[:, 2])-self.base_z
            normals[:, 7, 2] = 1.
        return gaps, normals

    def _midpoint(self, local0, local1, x0, x1):
        gaps, normals = self._midpoint_geometry(local0, local1, x0, x1)
        overlap = np.where(gaps < self.particle_radius,
                           self.particle_radius-gaps, 0.)
        return self.kn*overlap[:, :, None]*normals

    def _normal_body(self, g0, g1, relative, x0, x1):
        """Original direct-potential DG retained for v1/v2 reproduction."""
        midforce = self._midpoint(g0['local_positions_m'], g1['local_positions_m'], x0, x1)
        potential_change = g1['potential_J']-g0['potential_J']
        squared_step = np.sum(relative**2, axis=2)
        length_scale = np.maximum(np.maximum(np.linalg.norm(g0['local_positions_m'], axis=1),
            np.linalg.norm(g1['local_positions_m'], axis=1)),
            max(self.radius, self.height, self.particle_radius))
        small = squared_step <= (64*np.finfo(float).eps*length_scale[:, None])**2
        correction = np.zeros_like(squared_step)
        defect = -potential_change-np.sum(midforce*relative, axis=2)
        np.divide(defect, squared_step, out=correction, where=~small)
        normal_body = midforce+correction[:, :, None]*relative
        return normal_body, midforce, small, length_scale, np.zeros_like(squared_step)

    def evaluate(self, x0, x1, pose0, pose1, dt, history=None, check_path=True):
        x0 = np.asarray(x0, float); x1 = np.asarray(x1, float)
        if (x0.ndim != 2 or x0.shape[1] != 3 or x1.shape != x0.shape or
                not np.isfinite(x0).all() or not np.isfinite(x1).all() or
                not np.isfinite(dt) or dt <= 0.):
            raise ValueError('finite equal (N,3) positions and positive dt required')
        origin0, rotation0 = _pose(pose0); origin1, rotation1 = _pose(pose1)
        old = np.zeros((len(x0), 8, 3)) if history is None else np.asarray(history, float)
        if old.shape != (len(x0), 8, 3) or not np.isfinite(old).all():
            raise ValueError('finite (N,8,3) history required')
        g0 = self.geometry(x0, pose0); g1 = self.geometry(x1, pose1)
        sweep = (self._sweep(x0, x1, origin0, origin1, rotation0, rotation1, g0, g1)
                 if check_path else dict(certified=False, status='NOT_RUN_TANGENT_PERTURBATION_ONLY'))
        relative = np.empty_like(old)
        relative[:, :7] = (g1['local_positions_m']-g0['local_positions_m'])[:, None, :]
        relative[:, 7] = x1-x0
        potential_change = g1['potential_J']-g0['potential_J']
        squared_step = np.sum(relative**2, axis=2)
        normal_body, midforce, small, length_scale, distance_roundoff = self._normal_body(
            g0, g1, relative, x0, x1)
        magnitude = np.linalg.norm(normal_body, axis=2)
        alignment = np.sum(normal_body*g1['normal_body'], axis=2)
        # Include the original force magnitude in the roundoff bound: a very
        # small corrected force can cancel several finite intermediate terms.
        alignment_tolerance = 256*np.finfo(float).eps*(magnitude+
            np.linalg.norm(midforce, axis=2)+g0['magnitude_N']+g1['magnitude_N'])
        if np.any(alignment < -alignment_tolerance):
            raise ContactGeometryError('discrete normal attracts along endpoint geometric normal; subdivide increment')
        effective_normal = np.zeros_like(normal_body)
        np.divide(normal_body, magnitude[:, :, None], out=effective_normal,
                  where=magnitude[:, :, None] > 0.)

        # Endpoint contact loss releases history even when the interval-average
        # normal force remains nonzero.  No invisible suspended spring remains.
        friction_active = g1['active'] & (magnitude > 0.)
        tangential_step = relative-np.sum(relative*effective_normal, axis=2)[:, :, None]*effective_normal
        projected_old = old-np.sum(old*effective_normal, axis=2)[:, :, None]*effective_normal
        projected_old = np.where(friction_active[:, :, None], projected_old, 0.)
        trial = np.where(friction_active[:, :, None], projected_old+tangential_step, 0.)
        trial_norm = np.linalg.norm(trial, axis=2)
        limit = self.mu*magnitude/self.kt
        sliding = trial_norm > limit
        factor = np.ones_like(trial_norm)
        factor[sliding] = limit[sliding]/trial_norm[sliding]
        updated = trial*factor[:, :, None]
        plastic_step = trial-updated
        tangential_body = -self.kt*updated
        friction = self.kt*np.sum(updated*plastic_step, axis=2)
        algorithmic = .5*self.kt*np.sum((updated-projected_old)**2, axis=2)
        energy0 = .5*self.kt*np.sum(old**2, axis=2)
        energy1 = .5*self.kt*np.sum(updated**2, axis=2)
        release = energy0-.5*self.kt*np.sum(projected_old**2, axis=2)
        tangential_relative_work = np.sum(tangential_body*relative, axis=2)
        tangential_balance = tangential_relative_work+energy1-energy0+friction+algorithmic+release
        tangential_scale = (np.abs(tangential_relative_work)+energy0+energy1+
                            np.abs(friction)+algorithmic+np.abs(release))
        # P = I-n*n^T is an orthogonal projector in real arithmetic only.
        # Its absolute forward error scales with the INPUT displacement and
        # old spring, even if exact P*relative is zero.  In that case the
        # computed spring can be O(eps*relative), its work O(eps*kt*relative^2),
        # while its stored energy is O(eps^2*kt*relative^2).  An energy-output-
        # only relative certificate would incorrectly reject roundoff here.
        # The modest operation-count bound below also covers normalization
        # and radial return.  It is a numerical identity certificate, not an
        # allowance in the process's physical energy budget.  Forces, history,
        # heat and the raw signed identity residual are left unmodified.
        projection_input_scale = self.kt*(np.linalg.norm(old, axis=2)+
            np.linalg.norm(relative, axis=2))*np.linalg.norm(relative, axis=2)
        projection_roundoff = 256*np.finfo(float).eps*projection_input_scale
        dissipation_tolerance = 256*np.finfo(float).eps*tangential_scale
        tangent_tolerance = dissipation_tolerance+projection_roundoff
        # Keep the pre-existing sign certificate separate: no new allowance
        # for negative physical sliding heat or signed release is introduced.
        if np.any(friction < -dissipation_tolerance) or np.any(release < -dissipation_tolerance):
            raise ArithmeticError('negative raw friction/release energy; candidate rejected')
        if np.any(np.abs(tangential_balance) > tangent_tolerance):
            raise ArithmeticError('tangential spring work identity failed')

        normal = normal_body.copy(); normal[:, :7] = normal_body[:, :7]@rotation1.T
        tangential = tangential_body.copy(); tangential[:, :7] = tangential_body[:, :7]@rotation1.T
        relative_world = relative.copy(); relative_world[:, :7] = relative[:, :7]@rotation1.T
        normals_world = effective_normal.copy(); normals_world[:, :7] = effective_normal[:, :7]@rotation1.T
        wall_step = (x1-x0)[:, None, :]-relative_world
        total = normal+tangential
        normal_work = np.sum(normal*relative_world, axis=2)
        identity = potential_change+normal_work
        identity_scale = (g0['potential_J']+g1['potential_J']+np.abs(normal_work)+
                          self.kn*squared_step)
        identity_tolerance = 1024*np.finfo(float).eps*identity_scale
        # The derivative-limit branch can contain endpoint roundoff in a zero
        # potential contact activation; retain its absolute arithmetic bound.
        identity_tolerance += np.where(small, 256*np.finfo(float).eps*
            length_scale[:, None]*(g0['magnitude_N']+g1['magnitude_N']), 0.)
        identity_tolerance += distance_roundoff
        if np.any(np.abs(identity) > identity_tolerance):
            raise ArithmeticError('discrete normal work identity failed')
        endpoint_error = potential_change+np.sum(g1['force_by_feature_N']*relative_world, axis=2)
        particle_work = np.sum(total*(x1-x0)[:, None, :], axis=2)
        wall_work = np.sum(total*wall_step, axis=2)
        return dict(classification=self.classification, features=FEATURES,
            normal_model=self.normal_model,
            friction_model='BE_BODY_SPRING_COULOMB_WITH_DISCRETE_NORMAL_DIRECTION_AND_LOAD',
            work_quadrature='DISCRETE_NORMAL_AND_BE_TANGENTIAL_FORCE_FINITE_BODY_COORDINATE_WALL_DISPLACEMENT',
            force_N=np.sum(total, axis=1), force_by_feature_N=total,
            normal_force_by_feature_N=normal, tangential_force_by_feature_N=tangential,
            history=updated, active=g1['active'], friction_active=friction_active,
            sliding=sliding & friction_active, gaps_m=g1['gaps_m'], normals=normals_world,
            geometric_normals=g1['normal_world'], normal_load_N=magnitude,
            endpoint_normal_load_N=g1['magnitude_N'], normal_alignment_N=alignment,
            closest_points_m=g1['closest_points_m'], overlap_m=g1['overlap_m'],
            penetration_m=g1['overlap_m'], normal_potential0_J=g0['potential_J'],
            normal_potential1_J=g1['potential_J'], tangential_potential0_J=energy0,
            tangential_potential1_J=energy1, friction_dissipation_J=friction,
            algorithmic_dissipation_J=algorithmic, release_energy_J=release,
            tangential_work_identity_error_J=tangential_balance,
            tangential_work_identity_tolerance_J=tangent_tolerance,
            tangential_projection_roundoff_bound_J=projection_roundoff,
            # Compatibility key now describes the returned DG force.  The
            # hypothetical old endpoint quadrature has its own separate key.
            normal_endpoint_work_error_J=identity, dg_identity_error_J=identity,
            dg_identity_tolerance_J=identity_tolerance,
            dg_distance_roundoff_bound_J=distance_roundoff,
            normal_endpoint_quadrature_error_J=endpoint_error,
            normal_derivative_limit_branch=small,
            normal_particle_work_J=np.sum(normal*(x1-x0)[:, None, :], axis=2),
            normal_wall_work_J=np.sum(normal*wall_step, axis=2),
            particle_work_J=particle_work, wall_work_J=wall_work,
            relative_step_body_m=relative, wall_step_m=wall_step,
            particle_reaction_on_tool_N=-np.sum(total[:, :7], axis=(0, 1)),
            particle_reaction_torque_on_tool_Nm=-np.sum(np.cross(
                (x1-origin1)[:, None, :], total[:, :7]), axis=(0, 1)),
            particle_reaction_on_base_N=-np.sum(total[:, 7], axis=0), sweep=sweep,
            dt_s=float(dt), release_energy_policy='UNASSIGNED_ACCOUNT_NOT_AUTOMATICALLY_HEAT',
            algorithmic_dissipation_policy='REPORTED_NOT_AUTOMATICALLY_HEAT')


class ChainRuleDiscreteContact(DiscreteContact):
    """DG of distance times the scalar penalty's distance secant.

    U=phi(d), phi(d)=kn/2*max(r-d,0)^2.  Applying the discrete chain rule
    keeps an affine plane-distance gradient normal even across activation,
    where a direct U-Gonzalez correction can create a large tangential force.
    Curved/edge distance branches still use a Gonzalez distance gradient and
    retain the outward-alignment rejection.  This is not contact event
    integration: shell contact occurring entirely between inactive endpoints
    can be missed; a time-refinement/process-domain check remains required.
    """
    classification = CHAIN_RULE_CLASSIFICATION
    normal_model = 'CHAIN_RULE_DISTANCE_GONZALEZ_SCALAR_PENALTY_SECANT_DUPLICATE_LIP_CONTRIBUTIONS_RETAINED'

    def _normal_body(self, g0, g1, relative, x0, x1):
        _, midnormal = self._midpoint_geometry(g0['local_positions_m'],
            g1['local_positions_m'], x0, x1)
        squared_step = np.sum(relative**2, axis=2)
        length_scale = np.maximum(np.maximum(np.linalg.norm(g0['local_positions_m'], axis=1),
            np.linalg.norm(g1['local_positions_m'], axis=1)),
            max(self.radius, self.height, self.particle_radius))
        small = squared_step <= (64*np.finfo(float).eps*length_scale[:, None])**2
        finite = np.isfinite(g0['gaps_m']) & np.isfinite(g1['gaps_m'])
        distance_change = np.zeros_like(squared_step)
        np.subtract(g1['gaps_m'], g0['gaps_m'], out=distance_change, where=finite)
        correction = np.zeros_like(squared_step)
        defect = distance_change-np.sum(midnormal*relative, axis=2)
        np.divide(defect, squared_step, out=correction, where=finite & ~small)
        gradient = midnormal+correction[:, :, None]*relative
        affine = (np.all(g0['normal_body'] == g1['normal_body'], axis=2) &
                  np.all(g0['normal_body'] == midnormal, axis=2))
        # An exactly unchanged plane normal defines the affine distance
        # derivative; keep it exactly instead of promoting gap-subtraction
        # roundoff into a spurious transverse force.  Retain its raw work
        # residual and certify endpoint distance arithmetic below.
        gradient = np.where(affine[:, :, None], midnormal, gradient)
        q0, q1 = g0['overlap_m'], g1['overlap_m']
        coefficient = np.zeros_like(squared_step)  # minus the phi secant
        both_active = g0['active'] & g1['active']
        coefficient[both_active] = .5*self.kn*(q0[both_active]+q1[both_active])
        crossing = g0['active'] ^ g1['active']
        potential_change = .5*self.kn*(q1-q0)*(q1+q0)
        np.divide(-potential_change, distance_change, out=coefficient, where=crossing)
        if not np.isfinite(coefficient).all() or np.any(coefficient < 0.):
            raise ArithmeticError('nonpositive-distance penalty secant invariant failed')
        force = coefficient[:, :, None]*gradient
        reference_force = coefficient[:, :, None]*midnormal
        distance_roundoff = np.where(affine, 1024*np.finfo(float).eps*
            length_scale[:, None]*coefficient, 0.)
        return force, reference_force, small, length_scale, distance_roundoff
