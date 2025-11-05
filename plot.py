import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



climb = True

if not climb:
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

if climb:
    data = pd.read_csv("output/P51D-climb-analysis.csv")
    altitude = data["Altitude [ft]"]
    Vz = data["Vz [ft/min]"]
    TAS_solved = data["TAS solved [mph]"]
    pitch_angle = data["Pitch angle [deg]"]
    J = data["J"]
    etaP = data["Etap"]
    consumptions = data["Consumption [L/h]"]
    MachBlade = data["Mach Blade Tip"]
    SuperChargerMode = data["Supercharger Mode"]
    alphas = np.loadtxt("output/P51D-climb-alpha.txt")
    low = data["Supercharger Mode"] == "Low blower"
    high = data["Supercharger Mode"] == "High blower"

    plt.plot(altitude, Vz, label="Vertical Speed Vz", marker='o')
    plt.xlabel("Altitude [ft]")
    plt.ylabel("Vertical Speed Vz [ft/min]")
    plt.title("Vertical Speed vs Altitude for P-51D Mustang (Climb)")
    plt.legend()
    plt.grid()
    plt.savefig("figures/P51D-climb-Vz-vs-Altitude.png")
    plt.clf()

    plt.scatter(data["Altitude [ft]"][low],  data["TAS solved [mph]"][low],  label="Low Blower")
    plt.scatter(data["Altitude [ft]"][high], data["TAS solved [mph]"][high], label="High Blower")
    plt.plot(data["Altitude [ft]"], data["TAS solved [mph]"], label="TAS solved", color='gray', alpha=0.5)
    plt.xlabel("Altitude [ft]")
    plt.ylabel("True Airspeed [mph]")
    plt.title("True Airspeed vs Altitude for P-51D Mustang (Climb)")
    plt.legend()
    plt.grid()
    plt.savefig("figures/P51D-climb-TAS-vs-Altitude.png")
    plt.clf()


    plt.scatter(data["Vz [ft/min]"][low],  data["TAS solved [mph]"][low],  label="Low Blower")
    plt.scatter(data["Vz [ft/min]"][high], data["TAS solved [mph]"][high], label="High Blower")
    plt.plot(data["Vz [ft/min]"], data["TAS solved [mph]"], label="TAS solved", color='gray', alpha=0.5)
    plt.xlabel("Vertical Speed Vz [ft/min]")
    plt.ylabel("True Airspeed [mph]")
    plt.title("True Airspeed vs Vertical Speed for P-51D Mustang (Climb)")
    plt.legend()
    plt.grid()
    plt.savefig("figures/P51D-climb-TAS-vs-Vz.png")
    plt.clf()


    for i in range(alphas.shape[0]):
        plt.plot(np.linspace(0, 1, alphas.shape[1]), np.degrees(alphas[i]), label=f"Point {i+1}")
    plt.xlabel("r/R")
    plt.ylabel("Angle of Attack α [deg]")
    plt.title("Angle of Attack Distribution along the Blade for P-51D Mustang (High Speed)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    plt.grid()
    plt.savefig("figures/P51D-climbing-AoA-distribution.png", bbox_inches='tight')
    plt.clf()

