import math
import sys


class Reaction:
    def __init__(self, p_binding = 1, p_dissociation = 1, c0_ligand = 1.0, c0_receptor = 1.0, c0_complex = 0.0):
        self.p_dissociation = p_binding
        self.p_binding = p_dissociation
        self._c0_ligands = c0_ligand
        self._c0_recptrs = c0_receptor
        self._c0_complex = c0_complex
        self._c_ligands = self._c0_ligands
        self._c_recptrs = self._c0_recptrs
        self._c_complex = self._c0_complex
        self._reset()
        # print("# initial concentrations: %d %d %d" % (self.c0_receptor, self.c0_ligand, self.c0_complex), file=sys.stderr)

    @property
    def c_receptor(self):
        return self._c_recptrs
    @property
    def c_ligand(self):
        return self._c_ligands
    @property
    def c_complex(self):
        return self._c_complex
    @property
    def c0_complex(self): return self._c0_complex
    @property
    def c0_ligand(self): return self._c0_ligands
    @property
    def c0_receptor(self): return self._c0_recptrs

    @c0_complex.setter
    def c0_complex(self, c0_rl):
        self._c0_complex = c0_rl
        self._reset()

    @c0_ligand.setter
    def c0_ligand(self, c0_l):
        self._c0_ligands = c0_l
        self._reset()

    @c0_receptor.setter
    def c0_receptor(self, c0_r):
        self._c0_recptrs = c0_r
        self._reset()

    def _reset(self):
        self._c_ligands = self._c0_ligands
        self._c_recptrs = self._c0_recptrs
        self._c_complex = self._c0_complex
        # print(f"# actual concentrations after reset: {self._c_ligands} {self._c_recptrs} {self._c_complex}", file=sys.stderr)

    def Kd(self):
        if self._c_complex == 0: return 0.0
        return self._c_recptrs*self._c_ligands/self._c_complex

    def self_test(self):
        nl = self.c_ligand
        nr = self.c_receptor
        nc = self.c_complex
        assert nl+nc == self._c0_ligands + self._c0_complex
        assert nr + nc == self._c0_recptrs + self._c0_complex

    def equilibrate(self, dt=0.001, epsilon=0.0001):
        delta = 1
        i_step = 1
        while delta/dt > epsilon:
            d_a = self.p_binding * self._c_ligands * self._c_recptrs * dt
            d_d = self.p_dissociation * self._c_complex * dt
            delta = d_a - d_d
            self._c_complex += delta
            self._c_ligands -= delta
            self._c_recptrs -= delta
            i_step += 1
        return i_step

def log_space(x_from, n_steps, step):
    return [x_from * pow(10,i*step) for i in range(n_steps)]

if __name__ == "__main__":

    c0_lig = log_space(0.1, 22, 0.1)
    for cl in c0_lig:
        sim = Reaction(p_binding=1, p_dissociation=1, c0_ligand=cl)
        sim.equilibrate()
        cL, cR, cRL = sim.c_ligand, sim.c_receptor, sim.c_complex
        print("%.6f %.6f %.6f %f %.3f" % (cL, cR, cRL, cl, -math.log10(sim.Kd())))

