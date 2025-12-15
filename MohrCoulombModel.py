import numpy as np

class SmoothMohrCoulombPlaneStrain:
    """
    Сглаженный Мора–Кулон для плоской деформации.
    σzz исключено явно и учитывается неявно через εzz = 0.
    """

    def __init__(self, sx, sy, tau_xy, nu, phi_deg, c):
        # напряжения (атрибуты)
        self.sx = sx
        self.sy = sy
        self.tau_xy = tau_xy

        # параметры материала
        self.nu = nu
        self.phi = np.radians(phi_deg)
        self.c = c

    # --------------------------------------------------
    # Неявное σzz из условия plane strain
    # --------------------------------------------------
    def _sigma_zz(self):
        return self.nu * (self.sx + self.sy)

    # --------------------------------------------------
    # Главные напряжения в плоскости XY
    # --------------------------------------------------
    def _principal_xy(self):
        c = 0.5 * (self.sx + self.sy)
        r = np.sqrt((0.5 * (self.sx - self.sy))**2 + self.tau_xy**2)
        s1 = c + r
        s2 = c - r
        return s1, s2

    # --------------------------------------------------
    # Инварианты p, q (3D, но σzz подставлено неявно)
    # --------------------------------------------------
    def _invariants(self):
        s1, s2 = self._principal_xy()
        s3 = self._sigma_zz()

        # среднее напряжение
        p = (s1 + s2 + s3) / 3.0

        # девиаторные главные
        s1d = s1 - p
        s2d = s2 - p
        s3d = s3 - p

        # интенсивность девиатора
        q = np.sqrt(1.5 * (s1d**2 + s2d**2 + s3d**2))

        return p, q, np.array([s1d, s2d, s3d])

    # --------------------------------------------------
    # Функция текучести (сглаженная)
    # --------------------------------------------------
    def _yield_function(self, p, q):
        return q + p * np.sin(self.phi) - self.c * np.cos(self.phi)

    # --------------------------------------------------
    # Угол главных напряжений (из trial-состояния)
    # --------------------------------------------------
    def _principal_angle(self):
        if abs(self.sx - self.sy) < 1e-14 and abs(self.tau_xy) < 1e-14:
            return 0.0
        return 0.5 * np.arctan2(2.0 * self.tau_xy, self.sx - self.sy)

    # --------------------------------------------------
    # Восстановление осевых напряжений из главных в плоскости
    # --------------------------------------------------
    @staticmethod
    def _axial_from_principal(s1, s2, theta):
        c = np.cos(theta)
        s = np.sin(theta)

        sx = s1 * c**2 + s2 * s**2
        sy = s1 * s**2 + s2 * c**2
        tau = (s1 - s2) * s * c

        return sx, sy, tau

    # --------------------------------------------------
    # Return mapping (геометрический)
    # --------------------------------------------------
    def return_mapping(self):
        # trial инварианты
        p_tr, q_tr, s_tr = self._invariants()
        f_tr = self._yield_function(p_tr, q_tr)

        # упругое поведение
        if f_tr <= 0.0:
            return self.sx, self.sy, self.tau_xy

        # направление девиатора
        norm_s = np.linalg.norm(s_tr)
        if norm_s < 1e-14:
            alpha = 0.0
        else:
            # новое q на поверхности текучести
            q_new = self.c * np.cos(self.phi) - p_tr * np.sin(self.phi)
            alpha = q_new / q_tr if q_tr > 1e-14 else 0.0

        # масштабирование девиатора
        s_new = alpha * s_tr

        # новые главные напряжения
        s1_tr, s2_tr = self._principal_xy()
        s3_tr = self._sigma_zz()

        s1 = p_tr + s_new[0]
        s2 = p_tr + s_new[1]
        s3 = p_tr + s_new[2]  # неявно, но согласовано

        # восстановление в осях XY
        theta = self._principal_angle()
        sx_new, sy_new, tau_new = self._axial_from_principal(s1, s2, theta)

        return sx_new, sy_new, tau_new

sx = 1
sy = 1
tau = 1
fi = 30
c = 0
nu = 0.3

model = SmoothMohrCoulombPlaneStrain(sx, sy, tau, nu, fi, c)
print(model.return_mapping())
