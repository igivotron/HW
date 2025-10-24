import numpy as np
from math import pi, sin, cos, atan, radians, degrees
import matplotlib.pyplot as plt
import pandas as pd
import time 

class BEM:
    def __init__(self, diameter, hub_diameter, num_blades, chord, pitch):
        self.diameter = diameter

        self.num_blades = num_blades
        self.chord = chord
        self.pitch = radians(pitch)
        
        self.epsilon = 1e-6
        self.relax = 0.3
        self.max_iter = 1000

        self.r = np.arange(hub_diameter/2, diameter/2, 0.001)
        self.a_list = np.zeros_like(self.r)
        self.aprime_list = np.zeros_like(self.r)

    
    def get_Cl(self, alpha):
        return 2 * pi * alpha
    
    def get_Cd(self, alpha):
        return 0
    
    def get_pitch(self):
        return self.pitch
    
    def iterate(self, a0, aprim0, J, r):
        a = a0
        aprim = aprim0
        iter = 0

        while True:
            if iter > self.max_iter: break

            phi = atan( ((J*self.diameter) / (2*pi*r)) * ((1+a) / (1-aprim)) )

            sigma = (self.num_blades * self.chord)/(2 * pi * r)
            alpha = self.get_pitch() - phi

            Cl = self.get_Cl(alpha)
            Cd = self.get_Cd(alpha)
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
        
        for j in range(J.size):
            self.solve(J[j])
            kT[j] = ((4 * pi * (J[j]**2))/(self.diameter**2)) * np.trapz( self.r * (1 + self.a_list) * self.a_list, self.r)
            kQ[j] = ((J[j] * 8 * pi**2)/(self.diameter**4)) * np.trapz( self.r**3 * (1 + self.a_list) * self.aprime_list, self.r)
            kP[j] = 2*pi * kQ[j]
            if kT[j] < 0: kT[j] = None
            if kQ[j] < 0: kQ[j] = None
            if kP[j] < 0: kP[j] = None
            
            if kT[j] is None and kQ[j] is None and kP[j] is None: break

        etaP = (kT * J) / kP

        return kT, kQ, kP, etaP
    
    

def get_solution():
    data = pd.read_csv("Verification.txt", delim_whitespace=True)
    J = data["J"].values
    kT = data["kT"].values
    kQ = data["kQ"].values
    kP = data["kP"].values
    etaP = data["etaP"].values
    return J, kT, kQ, kP, etaP

if __name__ == "__main__":
    J, kT_ref, kQ_ref, kP_ref, etaP_ref = get_solution()
    bem = BEM(diameter=1, hub_diameter=0.25, num_blades=2, chord=0.15, pitch=25)
    kT, kQ, kP, etaP = bem.get_ks(J)








