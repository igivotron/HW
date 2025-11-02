from BEM import BEM
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import pi, cos, asin, sqrt, sin
import scipy as sp
from stdatm import stdatm

class Aircraft:
    def __init__(self, Mass, Awing, WingSpan, Cd0, e, bem):
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
        self.theta = 0
        self.bem = bem

    def getDrag(self, rho, V, theta):
        return 0.5 * rho * V**2 * self.A * (self.Cd0 + self.K * ( (self.M * 9.81 * cos(theta)) / (0.5 * rho * V**2 * self.A) )**2 )
    
    def getCl(self, rho, V, theta):
        return (self.M * 9.81 * cos(theta)) / (0.5 * rho * V**2 * self.A)

    def fixTheta(self, Vz, V):
        self.theta = asin(Vz / V)
        return self.theta

class FlightData:
    def __init__(self, RPM, MP, power, altitude, TAS_measured, Vz, reduce):
        self.altitude = altitude*0.3048 # ft to m
        self.V = TAS_measured *0.44704 # mph to m/s
        self.Vz = Vz *0.3048/60  # ft/min to m/s
        self.power = power * 745.7  # bhp to W
        self.MP = MP * 3386.388 # inHg to Pa
        self.RPM = RPM

        self.n = (RPM / 60) * reduce  # rev/s
        self.omega = self.n * 2 * pi
        
        self.getAtmosphere()

    def getAtmosphere(self):
        p, T, rho = stdatm(self.altitude)
        self.p = p
        self.T = T
        self.rho = rho


class FlightAnalysis:
    def __init__(self, aircraft, flightdata):
        self.aircraft = aircraft
        self.flightdata = flightdata
    
    def function(self, kQ, kT):
        P = self.flightdata.power
        rho  = self.flightdata.rho
        n = self.flightdata.n
        D = self.aircraft.bem.diameter
        V = self.flightdata.V
        Drag = self.aircraft.getDrag(rho, V, self.aircraft.theta)
        M = self.aircraft.M


        return [kQ - (P / (rho * 2 * pi * n**3 * D**5)), \
                kT - Drag / (rho * n**2 * D**4) - M*9.81*sin(self.aircraft.theta) / (rho * n**2 * D**4)]


    def find_kT_kQ(self):
        kT_init = 0.1
        kQ_init = 0.01

        def equations(vars):
            kT, kQ = vars
            return self.function(kQ, kT)

        kT_sol, kQ_sol = sp.optimize.fsolve(equations, (kT_init, kQ_init))

        return kT_sol, kQ_sol
    
    
bem = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=20)
aircraft = Aircraft(Mass=3600, Awing=16.2, WingSpan=10.97, Cd0=0.025, e=0.8, bem=bem)

flight = FlightData(RPM=3000, MP=60.5, power=1450, altitude=5000, TAS=363, Vz=0,  reduce=0.477)
analysis = FlightAnalysis(aircraft, flight)

kT, kQ = analysis.find_kT_kQ()
print(f"Calculated kT: {kT}, kQ: {kQ}")
