import numpy as np
from BEM import BEM
from FlightAndAircraft import Aircraft, FlightData, FlightAnalysis
from stdatm import stdatm
from math import cos, pi




blade = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=23)
aircraft = Aircraft(Mass=8430, Awing=21.83, WingSpan=11.28, Cd0=0.0163, e=0.8, bem=blade)
flight = FlightData(RPM=3000,
                    MP=61,
                    power=1400,
                    altitude=0,
                    TAS_measured=150,
                    Vz=0,
                    reduce=0.477)

analysis = FlightAnalysis(aircraft, flight)

p, T, rho = stdatm(flight.altitude)

beta = analysis.solve3()[0]
blade.beta_pitch = beta

Cl = aircraft.getCl(rho, flight.TAS_measured, 0)
Cd = aircraft.getCd(rho, flight.TAS_measured, 0)
print(f"Cl = {Cl:.4f}, Cl_max = 2")
print(f"Cd = {Cd:.4f}")

## TakeOff thrust and pitch angle
n = flight.n
D = blade.diameter
J = flight.TAS_measured / (n * D)

kT, kQ, kP, etaP = blade.get_k(J)
T = kT * rho * n**2 * D**4
D = aircraft.getDrag(rho, flight.TAS_measured, 0)

print(f"Advance ratio J: {J:.4f}")
print(f"Beta pitch angle: {np.degrees(blade.beta_pitch):.2f} deg")
print(f"Thrust at take-off: {T:.2f} N")
print(f"Drag at take-off: {D:.2f} N")

ax = (T - D) / aircraft.M
ay = 9.81
a = np.sqrt(ax**2 + ay**2)
g = a/9.81

print(f"Acceleration: {g:.2f} g")


