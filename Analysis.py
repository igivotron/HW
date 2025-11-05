from FlightAndAircraft import FlightData, Aircraft, FlightAnalysis
from BEM import BEM

import numpy as np
import pandas as pd

flightData = pd.read_csv("FlightData/P51D-hs.csv")

bem = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=20)
aircraft = Aircraft(Mass=8430, Awing=21.83, WingSpan=11.28, Cd0=0.0163, e=0.8, bem=bem)

N = flightData.shape[0]
solutions = np.zeros((N, 2))
consumptions = np.zeros(N)
MachBlade = np.zeros(N)
alphaList = np.zeros((N, bem.r.size))
J = np.zeros(N)
etaP = np.zeros(N)
for i in range(N):
    flight = FlightData(RPM=flightData["RPM"][i],
                        MP=flightData["MP"][i],
                        power=flightData["Pengine"][i],
                        altitude=flightData["Altitude"][i],
                        TAS_measured=flightData["TAS"][i],
                        Vz=0,
                        reduce=0.477)
    analysis = FlightAnalysis(aircraft, flight)
    u0_sol, theta_pitch_sol = analysis.solve()
    solutions[i, 0] = u0_sol
    solutions[i, 1] = theta_pitch_sol
    J[i] = u0_sol / ( (flight.RPM / 60 * 0.477) * bem.diameter )
    etaP[i] = bem.get_k(J[i])[3]
    SFC = aircraft.getConsomption(flight.power * 745.7, flightData["Supercharger Mode"][i])  # W to bhp

    a = np.sqrt(1.4 * 287 * flight.T)
    # MachBlade[i] = bem.BladeTipSpeed(flight.RPM, 0.477) / a
    MachBlade[i] = np.sqrt( (bem.BladeTipSpeed(flight.RPM, 0.477))**2 + (u0_sol)**2) / a
    consumptions[i] = SFC/1000 # in L/h
    alphaList[i] = bem.alpha
    print(f"Point {i+1}/{N} solved.")


table = pd.DataFrame(columns=["Supercharger Mode", "RPM", "MP [inHg]", "Pengine [bhp]", "Altitude [ft]", "TAS measured [mph]", "TAS solved [mph]", "Speed difference","Pitch angle [deg]", "J", "Etap", "Consumption [L/h]", "Mach Blade Tip"])
table["Supercharger Mode"] = flightData["Supercharger Mode"]
table["RPM"] = flightData["RPM"]
table["MP [inHg]"] = flightData["MP"]
table["Pengine [bhp]"] = flightData["Pengine"]
table["Altitude [ft]"] = flightData["Altitude"]
table["TAS measured [mph]"] = flightData["TAS"]
table["TAS solved [mph]"] = solutions[:, 0] / 0.44704
table["Speed difference in %"] = (table["TAS solved [mph]"] - table["TAS measured [mph]"]) / table["TAS measured [mph]"] * 100
table["Pitch angle [deg]"] = np.degrees(solutions[:, 1])
table["J"] = J
table["Etap"] = etaP
table["Consumption [L/h]"] = consumptions
table["Mach Blade Tip"] = MachBlade

table.to_csv("output/P51D-HS-analysis.csv", index=False)
print(table)

np.savetxt("output/P51D-HS-alpha.txt", alphaList)