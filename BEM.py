import numpy as np
from math import pi, sin, cos, atan, radians, degrees, tan, sqrt
import matplotlib.pyplot as plt
import pandas as pd
from naca16_509_python.naca16_509_m06_clcd import naca16_509_m06
from stdatm import stdatm
import scipy.optimize as optimize


class BEM:
    def __init__(self, diameter, hub_diameter, num_blades, chord, beta_ref, beta_pitch):
        self.diameter = diameter

        self.num_blades = num_blades
        self.chord = chord
        self.beta_ref = radians(beta_ref)
        self.beta_pitch = radians(beta_pitch)
        self.p_ref = 2*pi*0.75*(diameter/2) * tan(self.beta_ref)
        
        self.epsilon = 1e-6
        self.relax = 0.3
        self.max_iter = 1000

        self.r = np.arange(hub_diameter/2, diameter/2, 0.01)
        self.a_list = np.zeros_like(self.r)
        self.aprime_list = np.zeros_like(self.r)
    
    def get_C(self, alpha):
        cl, cd, _ = naca16_509_m06(alpha)
        return cl, cd
    
    def get_pitch(self, r):
        beta = atan(self.p_ref / (2 * pi * r)) + (self.beta_pitch - self.beta_ref)
        return beta
    
    def iterate(self, a0, aprim0, J, r):
        a = a0
        aprim = aprim0
        iter = 0

        while True:
            if iter > self.max_iter: break

            phi = atan( ((J*self.diameter) / (2*pi*r)) * ((1+a) / (1-aprim)) )

            sigma = (self.num_blades * self.chord)/(2 * pi * r)
            alpha = self.get_pitch(r) - phi

            Cl, Cd = self.get_C(alpha)
            Cn = Cl * cos(phi) - Cd * sin(phi)
            Ct = Cl * sin(phi) + Cd * cos(phi)

            a_new = ((sigma*Cn)*(1+a)) / (2*(1 - cos(2*phi)))
            aprim_new = ((sigma*Ct)*(1 - aprim)) / (2*sin(2*phi))

            if abs(a - a_new) < self.epsilon and abs(aprim - aprim_new) < self.epsilon:
                break

            a = (1 - self.relax) * a + self.relax * a_new
            aprim = (1 - self.relax) * aprim + self.relax * aprim_new

            iter += 1
        
        return a, aprim
        
    
    def solve(self, J):
        a_curent = 0
        aprime_curent = 0
        for i in range(self.r.size):
            a_new, aprim_new = self.iterate(a_curent, aprime_curent, J, self.r[i])
            self.a_list[i] = a_new
            self.aprime_list[i] = aprim_new
            a_curent = a_new
            aprime_curent = aprim_new
        return
    
    def get_ks(self, J):
        kT = np.zeros_like(J)
        kQ = np.zeros_like(J)
        kP = np.zeros_like(J)
        etaP = np.zeros_like(J)
        
        for j in range(J.size):
            self.solve(J[j])
            kT[j] = ((4 * pi * (J[j]**2))/(self.diameter**2)) * np.trapz( self.r * (1 + self.a_list) * self.a_list, self.r)
            kQ[j] = ((J[j] * 8 * pi**2)/(self.diameter**4)) * np.trapz( self.r**3 * (1 + self.a_list) * self.aprime_list, self.r)
            kP[j] = 2*pi * kQ[j]
            etaP[j] = (kT[j] * J[j]) / kP[j]

            if abs(etaP[j]) >= 1 or np.isnan(etaP[j]):
                etaP[j:] = np.nan
                kT[j:] = np.nan
                kQ[j:] = np.nan
                kP[j:] = np.nan
                break

            print(f"kT: {kT[j]}, kQ: {kQ[j]}, kP: {kP[j]}, etaP: {etaP[j]}")
            

        return kT, kQ, kP, etaP
    


if __name__ == "__main__":
    settings = True
    optimisation = False

    if settings:
        blade10 = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=10)
        blade20 = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=20)
        blade30 = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=30)
        blade40 = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=40)
        blade50 = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=50)
        blade60 = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=60)

        J = np.linspace(0.1, 5, 200)
        kT10, kQ10, kP10, etaP10 = blade10.get_ks(J)
        kT20, kQ20, kP20, etaP20 = blade20.get_ks(J)
        kT30, kQ30, kP30, etaP30 = blade30.get_ks(J)
        kT40, kQ40, kP40, etaP40 = blade40.get_ks(J)
        kT50, kQ50, kP50, etaP50 = blade50.get_ks(J)
        kT60, kQ60, kP60, etaP60 = blade60.get_ks(J)

        plt.subplot(2,2,1)
        plt.plot(J, kT10, label="Beta pitch 10°")
        plt.plot(J, kT20, label="Beta pitch 20°")
        plt.plot(J, kT30, label="Beta pitch 30°")
        plt.plot(J, kT40, label="Beta pitch 40°")
        plt.plot(J, kT50, label="Beta pitch 50°")
        plt.plot(J, kT60, label="Beta pitch 60°")
        plt.xlabel("J")
        plt.ylabel("kT")
        # plt.legend()
        plt.grid()
        plt.subplot(2,2,2)
        plt.plot(J, kQ10, label="Beta pitch 10°")
        plt.plot(J, kQ20, label="Beta pitch 20°")
        plt.plot(J, kQ30, label="Beta pitch 30°")
        plt.plot(J, kQ40, label="Beta pitch 40°")
        plt.plot(J, kQ50, label="Beta pitch 50°")
        plt.plot(J, kQ60, label="Beta pitch 60°")
        plt.xlabel("J")
        plt.ylabel("kQ")
        # plt.legend()
        plt.grid()
        plt.subplot(2,2,3)
        plt.plot(J, kP10, label="Beta pitch 10°")
        plt.plot(J, kP20, label="Beta pitch 20°")
        plt.plot(J, kP30, label="Beta pitch 30°")
        plt.plot(J, kP40, label="Beta pitch 40°")
        plt.plot(J, kP50, label="Beta pitch 50°")
        plt.plot(J, kP60, label="Beta pitch 60°")
        plt.xlabel("J")
        plt.ylabel("kP")
        # plt.legend()
        plt.grid()
        plt.subplot(2,2,4)
        plt.plot(J, etaP10, label="Beta pitch 10°")
        plt.plot(J, etaP20, label="Beta pitch 20°")
        plt.plot(J, etaP30, label="Beta pitch 30°")
        plt.plot(J, etaP40, label="Beta pitch 40°")
        plt.plot(J, etaP50, label="Beta pitch 50°")
        plt.plot(J, etaP60, label="Beta pitch 60°")
        plt.xlabel("J")
        plt.ylabel("etaP")
        # plt.legend()
        plt.tight_layout()
        plt.grid()
        plt.savefig("BEM_results_beta_pitch_10.png")
        plt.show()

    if optimisation:
        blade = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=0)
        n = (3000/60) *0.477
        omega = n * 2 * pi
        p, t, rho = stdatm(20e3*0.3048)  # Standard atmosphere at 20000 ft
        M = 0.5
        u0 = M * sqrt(1.4 * 287 * t)
        J = 2*pi*u0 / (blade.diameter * omega)  # for 3000 RPM
        print(f"J: {J}, rho: {rho}, t: {t} K")

        def objective(beta_pitch):
            blade.beta_pitch = radians(beta_pitch)
            _, _, _, etaP = blade.get_ks(np.array([J]))
            if np.isnan(etaP[0]): return 1e6
            etaP[0] = min(etaP[0], 1.0)
            return -etaP[0]
        
        result = optimize.minimize_scalar(objective, bounds=(0, 90), method='bounded')
        optimal_beta_pitch = result.x
        blade.beta_pitch = radians(optimal_beta_pitch)
        kT_opt, kQ_opt, kP_opt, etaP_opt = blade.get_ks(np.array([J]))
        T = kT_opt[0] * rho * n**2 * blade.diameter**4
        Q = kQ_opt[0] * rho * n**2 * blade.diameter**5
        P = kP_opt[0] * rho * n**3 * blade.diameter**5
        
        print(f"Optimal beta pitch: {optimal_beta_pitch} degrees")
        print(f"At optimal beta pitch, kT: {kT_opt[0]}, kQ: {kQ_opt[0]}, kP: {kP_opt[0]}, etaP: {etaP_opt[0]}")
        print(f"Thrust: {T} N, Torque: {Q} Nm, Power: {P} W")



 








