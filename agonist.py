import math
import random
import sys


class Reaction:
    def __init__(self, p_binding = 1, p_dissociation = 1, c0_receptor = 1000, c0_ligand = 1000, c0_complex = 0):
        self.p_dissociation = p_binding
        self.p_binding = p_dissociation
        self.__c0_ligands = c0_ligand
        self.__c0_recptrs = c0_receptor
        self.__c0_complex = c0_complex
        self.__c_ligands = self.__c0_ligands
        self.__c_recptrs = self.__c0_recptrs
        self.__c_complex = self.__c0_complex
        self.__reset()
        self.__c_avg_ligands = 0
        self.__c_avg_recptrs = 0
        self.__c_avg_complex = 0
        self.__n_avg = 0
        print("# initial concentrations: %d %d %d" % (self.c0_receptor, self.c0_ligand, self.c0_complex), file=sys.stderr)

    @property
    def c_receptor(self):
        return self.__c_recptrs
    @property
    def c_ligand(self):
        return self.__c_ligands
    @property
    def c_complex(self):
        return self.__c_complex

    @property
    def c_avg_complex(self): return self.__c_avg_complex / self.__n_avg
    @property
    def c_avg_ligand(self): return self.__c_avg_ligands / self.__n_avg
    @property
    def c_avg_receptor(self): return self.__c_avg_recptrs / self.__n_avg

    @property
    def c0_complex(self): return self.__c0_complex
    @property
    def c0_ligand(self): return self.__c0_ligands
    @property
    def c0_receptor(self): return self.__c0_recptrs

    @c0_complex.setter
    def c0_complex(self, c0_rl):
        self.__c0_complex = c0_rl
        self.__reset()
    @c0_ligand.setter
    def c0_ligand(self, c0_l):
        self.__c0_ligands = c0_l
        self.__reset()

    @c0_receptor.setter
    def c0_receptor(self, c0_r):
        self.__c0_recptrs = c0_r
        self.__reset()

    def __reset(self):
        self.__c_ligands = self.__c0_ligands
        self.__c_recptrs = self.__c0_recptrs
        self.__c_complex = self.__c0_complex
        print(f"# actual concentrations after reset: {self.__c_ligands} {self.__c_recptrs} {self.__c_complex}", file=sys.stderr)

    def Kd(self):
        if self.c_avg_complex == 0: return 0.0
        return self.c_avg_ligand*self.c_avg_receptor/self.c_avg_complex

    def step(self, delta=1):
        p_a = self.p_binding * self.__c_ligands * self.__c_recptrs
        p_d = self.p_dissociation * self.__c_complex
        total = p_a + p_d
        p_a, p_d = p_a / total, p_d/total
        delta = delta if random.random() < p_a else -delta
        self.__c_complex += delta
        self.__c_recptrs -= delta
        self.__c_ligands -= delta

    def self_test(self):
        nl = self.c_ligand
        nr = self.c_receptor
        nc = self.c_complex
        assert nl+nc == self.__c0_ligands + self.__c0_complex
        assert nr + nc == self.__c0_recptrs + self.__c0_complex

    def simulate(self, n_repeats, n_step_factor=1):

        total = self.c0_complex+self.c0_ligand+self.c0_receptor
        for r in range(n_repeats):
            for k in range(int(total*n_step_factor)):
                self.step()
            self.__n_avg += 1
            self.__c_avg_complex += self.c_complex
            self.__c_avg_ligands += self.c_ligand
            self.__c_avg_recptrs += self.c_receptor
            self.self_test()

    def resets_avg(self):
        self.__n_avg = 0
        self.__c_avg_complex = 0
        self.__c_avg_ligands = 0
        self.__c_avg_recptrs = 0

if __name__ == "__main__":
    c0_lig = [100, 130, 160, 200, 250, 320, 400, 500, 630, 790, 1000, 1300, 1600,2000, 2500, 3200, 4000, 5000, 6300, 7900, 10000, 13000, 16000,20000, 25000]
    # c0_lig = [c*0.1 for c in c0_lig]
    for cl in c0_lig:
        sim = Reaction(p_binding=1, p_dissociation=1, c0_receptor=10000, c0_ligand=cl)
        ref = sim.c0_receptor
        sim.simulate(1000)
        sim.resets_avg()
        sim.simulate(100000)
        cL, cR, cRL = sim.c_avg_ligand / ref, sim.c_avg_receptor / ref, sim.c_avg_complex / ref
        print("%.6f %.6f %.6f %f %.3f" % (cL, cR, cRL, cl, -math.log10(sim.Kd())))

