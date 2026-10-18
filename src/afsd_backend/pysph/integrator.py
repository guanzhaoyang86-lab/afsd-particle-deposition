"""Execute pinned PySPH LeapFrogIntegrator and LeapFrogStep on host arrays.

Unmodified upstream drift/kick/drift; constitutive/thermal split is explicit and
first-order. This adapter does not claim the author's unavailable integrator code.
"""
import numpy as np
from pysph.sph.integrator import LeapFrogIntegrator
from pysph.sph.integrator_step import LeapFrogStep


class HostLeapfrog(LeapFrogIntegrator):
    def __init__(self, pa, acceleration_callback):
        self.stepper = LeapFrogStep()
        super().__init__(solid=self.stepper)
        self.pa = pa
        self.callback = acceleration_callback

    def advance(self, time, dt):
        self.time, self.dt = time, dt
        self.constraint_work = 0.
        self.one_timestep(time, dt)

    def stage1(self):
        p = self.pa
        self.stepper.stage1(slice(None), p.x, p.y, p.z, p.u, p.v, p.w, p.ax, p.ay, p.az, self.dt)

    def compute_accelerations(self):
        self.callback(self.time+self.dt/2, self.dt)

    def stage2(self):
        p = self.pa
        self.stepper.stage2(slice(None), p.x, p.y, p.z, p.u, p.au, p.v, p.av, p.w, p.aw, p.ax, p.ay, p.az, p.rho, p.arho, p.e, p.ae, self.dt)

    def update_domain(self):
        # The callback constructs PySPH NNPS after stage1 from current positions.
        pass

    def do_post_stage(self, dt, stage):
        pass
