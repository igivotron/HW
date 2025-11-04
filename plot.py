import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

HS = pd.read_csv("output/P51D-HS-analysis.csv")

TAS_measured = HS["TAS measured [mph]"]
TAS_solved = HS["TAS solved [mph]"]
pitch_angle = HS["Pitch angle [deg]"]
J = HS["J"]
etaP = HS["Etap"]
altitude = HS["Altitude [ft]"]

alpha = np.loadtxt("output/P51D-HS-alpha.txt")

plt.plot(altitude, TAS_measured, label="TAS measured", marker='o')
plt.plot(altitude, TAS_solved, label="TAS solved", marker='x')
plt.xlabel("Altitude [ft]")
plt.ylabel("True Airspeed [mph]")
plt.title("True Airspeed vs Altitude for P-51D Mustang (High Speed)")
plt.legend()
plt.grid()
plt.savefig("figures/P51D-HS-TAS-vs-Altitude.png")
plt.clf()

plt.plot(altitude, pitch_angle, label="Pitch angle", marker='s', color='orange')
plt.xlabel("Altitude [ft]")
plt.ylabel("Pitch Angle [deg]")
plt.title("Pitch Angle vs Altitude for P-51D Mustang (High Speed)")
plt.legend()
plt.grid()
plt.savefig("figures/P51D-HS-PitchAngle-vs-Altitude.png")
plt.clf()

plt.plot(altitude, J, label="Advance Ratio J", marker='^', color='green')
plt.xlabel("Altitude [ft]")
plt.ylabel("Advance Ratio J")
plt.title("Advance Ratio vs Altitude for P-51D Mustang (High Speed)")
plt.legend()
plt.grid()
plt.savefig("figures/P51D-HS-AdvanceRatio-vs-Altitude.png")
plt.clf()

plt.plot(altitude, etaP, label="Propulsive Efficiency ηP", marker='d', color='red')
plt.xlabel("Altitude [ft]")
plt.ylabel("Propulsive Efficiency ηP")
plt.title("Propulsive Efficiency vs Altitude for P-51D Mustang (High Speed)")
plt.legend()
plt.grid()
plt.savefig("figures/P51D-HS-PropulsiveEfficiency-vs-Altitude.png")
plt.clf()

for i in range(alpha.shape[0]):
    plt.plot(np.linspace(0, 1, alpha.shape[1]), np.degrees(alpha[i]), label=f"Point {i+1}")
plt.xlabel("r/R")
plt.ylabel("Angle of Attack α [deg]")
plt.title("Angle of Attack Distribution along the Blade for P-51D Mustang (High Speed)")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
plt.grid()
plt.savefig("figures/P51D-HS-AoA-distribution.png", bbox_inches='tight')
plt.clf()