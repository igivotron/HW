from BEM import BEM
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import pi, cos, asin, sin, radians, degrees
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
    
    def printAircraft(self):
        print(f"Mass: {self.M:.2f} kg")
        print(f"Wing Area: {self.A:.2f} m^2")
        print(f"Wing Span: {self.b:.2f} m")
        print(f"C_d0: {self.Cd0:.4f}")
        print(f"Oswald Efficiency Factor: {self.e:.4f}")
        print(f"Aspect Ratio: {self.AR:.2f}")
        print(f"Induced Drag Factor K: {self.K:.4f}")

class FlightData:
    def __init__(self, RPM, MP, power, altitude, TAS_measured, Vz, reduce):
        self.altitude = altitude*0.3048 # ft to m
        self.TAS_measured = TAS_measured *0.44704 # mph to m/s
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

    def printData(self):
        print(f"Altitude: {self.altitude:.2f} m")
        print(f"TAS measured: {self.TAS_measured:.2f} m/s")
        print(f"Vertical Speed: {self.Vz:.2f} m/s")
        print(f"Power: {self.power:.2f} W")
        print(f"Manifold Pressure: {self.MP/1e3:.2f} kPa")
        print(f"RPM: {self.RPM:.2f} rpm")
        print(f"Air density: {self.rho:.4f} kg/m^3")
        print(f"Air temperature: {self.T:.2f} K")
        print(f"Atmospheric pressure: {self.p/1000:.2f} kPa")


class FlightAnalysis:
    def __init__(self, aircraft, flightdata):
        self.aircraft = aircraft
        self.flightdata = flightdata

    def function(self, u0, beta_pitch):
        rho = self.flightdata.rho
        P = self.flightdata.power
        n = self.flightdata.n
        D = self.aircraft.bem.diameter
        M = self.aircraft.M
        A = self.aircraft.A
        Cd0 = self.aircraft.Cd0
        K = self.aircraft.K
        theta = self.aircraft.theta
        J = u0 / (n * D)
        if J < 0.1: J = 0.1
        self.aircraft.bem.beta_pitch = beta_pitch

        kT, kQ, _, _ = self.aircraft.bem.get_k(J)

        rep1 = kQ - (P / (rho * 2 * pi * n**3 * D**5))
        rep2 = kT - 0.5 * ((A * J**2)/(D**2)) * (Cd0 + K * ( (M * 9.81 * cos(theta)) / (0.5 * rho * (n * D * J)**2 * A) )**2 ) - M * 9.81 * sin(theta) / (rho * n**2 * D**4)
        return [rep1, rep2]
    

    def solve(self):
        u0_guess = self.flightdata.TAS_measured
        theta_pitch_guess = radians(40)

        def equations(vars):
            u0, theta_pitch = vars
            if theta_pitch < 0: theta_pitch = 0
            if theta_pitch > radians(60): theta_pitch = radians(60)
            return self.function(u0, theta_pitch)
        
        TAS_measured = self.flightdata.TAS_measured
        bounds = ((0.9*TAS_measured, radians(30)), (1.1*TAS_measured, radians(50)))
        sol = sp.optimize.least_squares(equations, (u0_guess, theta_pitch_guess), bounds=bounds)
        u0_sol, theta_pitch_sol = sol.x
        

        return u0_sol, theta_pitch_sol
    
    
    
bem = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=20)
aircraft = Aircraft(Mass=8430, Awing=21.83, WingSpan=11.28, Cd0=0.0163, e=0.8, bem=bem)

flight = FlightData(RPM=3000, MP=40.7, power=845, altitude=38000, TAS_measured=403, Vz=0,  reduce=0.477)
analysis = FlightAnalysis(aircraft, flight)
u0_sol, theta_pitch_sol = analysis.solve()

print(f"Solved True Air Speed: {u0_sol:.2f} m/s, {u0_sol*1.94384:.2f} knots, {u0_sol*2.23694:.2f} mph")
print(f"Solved Pitch Angle: {degrees(theta_pitch_sol):.2f} degrees")