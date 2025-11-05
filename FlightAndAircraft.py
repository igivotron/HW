from BEM import BEM
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import pi, cos, asin, sin, radians
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

        self.C1LB = -36.12
        self.C2LB = 1.785e-4
        self.C1HB = -20.13
        self.C2HB = 1.849e-4

    def getDrag(self, rho, V, theta):
        return 0.5 * rho * V**2 * self.A * (self.Cd0 + self.K * ( (self.M * 9.81 * cos(theta)) / (0.5 * rho * V**2 * self.A) )**2 )
    
    def getCl(self, rho, V, theta):
        return (self.M * 9.81 * cos(theta)) / (0.5 * rho * V**2 * self.A)

    def getTheta(self, Vz, V):
        self.theta = asin(Vz / V)
        return self.theta
    
    def getConsomption(self, power, SuperchargerMode):
        if SuperchargerMode == "Low blower":
            C1 = self.C1LB
            C2 = self.C2LB
        else:
            C1 = self.C1HB
            C2 = self.C2HB
        
        SFC = C1 + C2 * power  # in lb/h per hp
        SFC *= 3.78541
        return SFC
    
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
        theta = self.aircraft.getTheta(self.flightdata.Vz, u0)
        J = u0 / (n * D)
        if J < 0.1: J = 0.1
        self.aircraft.bem.beta_pitch = beta_pitch

        kT, kQ, kP, etaP = self.aircraft.bem.get_k(J)

        rep1 = kQ - (P / (rho * 2 * pi * n**3 * D**5))
        rep2 = kT - 0.5 * ((A * J**2)/(D**2)) * (Cd0 + K * ( (M * 9.81 * cos(theta)) / (0.5 * rho * (n * D * J)**2 * A) )**2 ) - M * 9.81 * sin(theta) / (rho * n**2 * D**4)
        return [rep1, rep2]
    

    def solve(self):
        u0_guess = self.flightdata.TAS_measured
        beta_pitch_guess = radians(40)

        def equations(vars):
            u0, beta_pitch = vars
            if beta_pitch < 0: beta_pitch = 0
            if beta_pitch > radians(60): beta_pitch = radians(60)
            if self.aircraft.Vz > u0: u0 = self.aircraft.Vz - 1
            return self.function(u0, beta_pitch)
        TAS_measured = self.flightdata.TAS_measured
        bounds = ((0.9*TAS_measured, radians(30)), (1.1*TAS_measured, radians(50)))
        sol = sp.optimize.least_squares(equations, (u0_guess, beta_pitch_guess), bounds=bounds)
        u0_sol, beta_pitch_sol = sol.x
        return u0_sol, beta_pitch_sol

    def solve2(self, u0_guess, beta_pitch_guess):
        def equations(vars):
            u0, beta_pitch = vars
            if beta_pitch < 0: beta_pitch = 0
            if beta_pitch > radians(60): beta_pitch = radians(60)
            return self.function(u0, beta_pitch)
        Vz = self.flightdata.Vz
        
        bounds = ((Vz, radians(0)), (2*u0_guess, radians(60)))
        sol = sp.optimize.least_squares(equations, (u0_guess, beta_pitch_guess), bounds=bounds)
        u0_sol, beta_pitch_sol = sol.x
        

        return u0_sol, beta_pitch_sol
    

