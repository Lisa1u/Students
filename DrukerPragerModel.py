import numpy as np

class DPModel:
    def __init__(self, fi, c, E, vu, sx, sy, tau_xy):
        self.fi = np.radians(fi)
        self.c = c
        self.E = E
        self.vu = vu
        self.sx = sx
        self.sy = sy
        self.tau_xy = tau_xy

    def _main_stresses(self):
        s1 = 0.5 * (self.sx + self.sy)
        s2 = np.sqrt(0.5 * (self.sx - self.sy)**2 + self.tau_xy**2)
        sigma_1 = s1 + s2
        sigma_2 = s1 - s2
        return sigma_1, sigma_2

    def _invariants(self):
        sigma_1, sigma_2 = self._main_stresses()



mat1 = DPModel(fi=30)
print(mat1.sin())
