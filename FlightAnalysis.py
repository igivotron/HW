from BEM import BEM
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import pi, cos

class Aircraft:
    def __init__(self, Mass, Awing, WingSpan, Cd0, e):
        """
        Mass in lb
        Awing in m^2
        WingSpan in m
        """
        self.M = Mass/2.205
        self.A = Awing
        self.b = WingSpan
        self.Cd0 = Cd0
        self.e = e
        self.AR = self.b**2 / self.A
        self.K = 1 / (pi * self.e * self.AR)

    def getDrag(self, rho, V, theta):
        return 0.5 * rho * V**2 * self.A * (self.Cd0 + self.K * ( (self.M * 9.81 * cos(theta)) / (0.5 * rho * V**2 * self.A) )**2 )
    
    def getCl(self, rho, V, theta):
        return (self.M * 9.81 * cos(theta)) / (0.5 * rho * V**2 * self.A)
        

hs = pd.read_csv("FlightData/P51D-hs.csv")
climb = pd.read_csv("FlightData/P51D-climb.csv")

P51D = Aircraft(Mass=8430, Awing=21.83, WingSpan=11.28, Cd0=0.0163, e=0.8)
bem = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=0)
