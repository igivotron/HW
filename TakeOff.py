import numpy as np
from BEM import BEM
from FlightAndAircraft import Aircraft, FlightData, FlightAnalysis
from stdatm import stdatm
from math import cos, pi, radians
from scipy.optimize import fsolve, least_squares




bem = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=20)
aircraft = Aircraft(Mass=8430, Awing=21.83, WingSpan=11.28, Cd0=0.0163, e=0.8, bem=bem)
flight = FlightData(RPM=3000,
                    MP=61,
                    power=1400,
                    altitude=0,
                    TAS_measured=150,
                    Vz=0,
                    reduce=0.477)

analysis = FlightAnalysis(aircraft, flight)

p, T, rho = stdatm(flight.altitude)

Cl = aircraft.getCl(rho, flight.TAS_measured, 0)
Cd = aircraft.getCd(rho, flight.TAS_measured, 0)
print(f"Cl = {Cl:.4f}, Cl_max = 2")
print(f"Cd = {Cd:.4f}")

## TakeOff thrust and pitch angle
n = flight.n
D = bem.diameter
J = flight.TAS_measured / (n * D)

def function_takeoff(beta_pitch):
    aircraft.bem.beta_pitch = beta_pitch
    K = aircraft.K
    Cd0 = aircraft.Cd0
    Pengine = flight.power
    u0 = flight.TAS_measured
    A = aircraft.A

    kT, kQ, kP, etaP = aircraft.bem.get_k(J)
    rep1 = kQ - Pengine / (rho * 2 * pi * n**3 * D**5)
    rep2 = kT - 0.5 * (A*J**2)/D**2 * (Cd0 + K * Cl**2)
    return rep1

solution = least_squares(function_takeoff, radians(25))
print(solution.x)
beta_pitch_solution = solution.x[0]
aircraft.bem.beta_pitch = beta_pitch_solution



kT, kQ, kP, etaP = aircraft.bem.get_k(J)
K = aircraft.K
Cd0 = aircraft.Cd0
T = kT * rho * n**2 * D**4
D = 0.5 * rho * (flight.TAS_measured**2) * aircraft.A * (Cd0 + K * Cl**2)
R = T - D
a = R / aircraft.M
G = a / 9.81

print(f"TakeOff pitch angle (degrees): {np.degrees(beta_pitch_solution):.2f}")
print(f"TakeOff Thrust (N): {T:.2f}")
print(f"TakeOff Drag (N): {D:.2f}")
print(f"TakeOff Excess Thrust (N): {R:.2f}")
print(f"TakeOff Acceleration (m/s²): {a:.2f}")
print(f"TakeOff Load Factor: {G:.2f}")






