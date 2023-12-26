import math
import sys

from agonist import Reaction, log_space


class ReactionCompetitive(Reaction):
    def __init__(self, p_binding = 1, p_dissociation = 1, p_a_binding = 1, p_a_dissociation = 1, c0_antagonist = 0.0, c0_ligand = 1.0, c0_receptor = 1.0, c0_complex = 0.0):
        super().__init__(p_binding, p_dissociation, c0_ligand, c0_receptor, c0_complex)
        self.p_a_dissociation = p_a_binding
        self.p_a_binding = p_a_dissociation
        self.__c0_antagonist = c0_antagonist
        self.__c_antagonist = self.__c0_antagonist
        self._c_competitive = 0.0
        self._reset()
        # print("# initial concentrations: %d %d %d" % (self.c0_receptor, self.c0_ligand, self.c0_complex), file=sys.stderr)

    @property
    def c_antagonist(self):
        return self.__c_antagonist
    @property
    def c0_antagonist(self): return self.__c0_antagonist

    def _reset_comp(self):
        super()._reset()
        self.__c_antagonist = self.__c0_antagonist

    def equilibrate(self, dt=0.001, epsilon=0.0001):
        delta = 1
        i_step = 1
        while delta/dt > epsilon:
            # L + R -> RL
            d_rl = self.p_binding * self._c_ligands * self._c_recptrs * dt
            # RL -> R + L
            d_r = self.p_dissociation * self._c_complex * dt
            # A + R -> RA
            d_ra = self.p_a_binding * self.__c_antagonist * self._c_recptrs * dt
            # RA -> R + A
            d_a = self.p_a_dissociation * self._c_competitive * dt
            self._c_complex += d_rl - d_r
            self._c_ligands -= d_rl - d_r
            self._c_recptrs -= d_rl - d_r + d_ra + d_a
            self.__c_antagonist -= d_ra - d_r
            self._c_competitive += d_ra - d_a
            i_step += 1
            delta = d_rl - d_r
        return i_step


if __name__ == "__main__":

    c0_lig = log_space(0.1, 22, 0.1)
    for cl in c0_lig:
        sim = ReactionCompetitive(p_binding=1, p_dissociation=1, c0_ligand=cl, c0_antagonist = 0.75)
        sim.equilibrate()
        cL, cR, cRL = sim.c_ligand, sim.c_receptor, sim.c_complex
        print("%.6f %.6f %.6f %f %.3f" % (cL, cR, cRL, cl, -math.log10(sim.Kd())))

