import numpy as np
from math import pi, sin, cos, atan, radians, tan, degrees, atan2
from naca16_509_python.naca16_509_m06_clcd import naca16_509_m06


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

        self.r = np.arange(hub_diameter/2, diameter/2, 0.005)
        self.a_list = np.zeros_like(self.r)
        self.aprime_list = np.zeros_like(self.r)
    
    def get_C(self, alpha):
        cl, cd, _ = naca16_509_m06(alpha)
        return cl, cd
    
    def get_pitch(self, r):
        # beta = atan(self.p_ref / (2 * pi * r)) + (self.beta_pitch - self.beta_ref)
        beta = atan2(self.p_ref, (2 * pi * r)) + (self.beta_pitch - self.beta_ref)
        return beta
    
    def iterate(self, a0, aprim0, J, r):
        a = a0
        aprim = aprim0
        iter = 0

        while True:
            if iter > self.max_iter: break

            # phi = atan( ((J*self.diameter) / (2*pi*r)) * ((1+a) / (1-aprim)) )
            phi = atan2( ((J*self.diameter) / (2*pi*r)) * (1 + a), (1 - aprim) )

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
            print(j + 1, "/", J.size)
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

            # print(f"kT: {kT[j]}, kQ: {kQ[j]}, kP: {kP[j]}, etaP: {etaP[j]}")
            

        return kT, kQ, kP, etaP