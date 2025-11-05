import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


data = pd.read_csv("output/P51D-climb-analysis.csv")

consumptions = data["Consumption [L/h]"]
altitude = data["Altitude [ft]"]
Vz = data["Vz [ft/min]"]
TAS = data["TAS solved [mph]"]

fuelDensity = 0.755  # kg/L for aviation gasoline

t = np.zeros(len(Vz))
t[0] = 0

for i in range(1, len(Vz)):
    delta_h = altitude[i] - altitude[i-1]
    avg_Vz = (Vz[i] + Vz[i-1]) / 2
    t[i] = (delta_h / avg_Vz) / 60  # hours

CONSO = np.zeros_like(altitude)
for i in range(1, len(Vz)):
    avg_Conso = (consumptions[i] + consumptions[i-1]) / 2  # L/h
    CONSO[i] = CONSO[i-1] + avg_Conso * t[i]  # L

T = np.cumsum(t)
CONSO_tot = CONSO[-1]
M_fuel = CONSO_tot * fuelDensity  # kg
M_tot = 8430*0.453592 + M_fuel  # kg
print(f"Total fuel consumed during climb: {CONSO_tot:.2f} L, {M_fuel:.2f} kg")
print(f"Initial aircraft mass: {M_tot:.2f} kg")

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(T, altitude, 'g-', marker='o', label='Altitude [ft]')
ax2.plot(T, CONSO, 'b-', marker='x', label='Cumulative Consumption [L]')
ax1.set_xlabel('Time [h]')
ax1.set_ylabel('Altitude [ft]')
ax2.set_ylabel('Cumulative Consumption [L]')


fig.legend(loc="upper left", bbox_to_anchor=(0.15,0.85))
plt.title("Altitude and Cumulative Fuel Consumption vs Time for P-51D Mustang (Climb)")
plt.grid()
plt.show()