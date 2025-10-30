import numpy as np
from math import pi, radians, sqrt
import matplotlib.pyplot as plt
from BEM import BEM
from stdatm import stdatm
import scipy.optimize as optimize


settings = True
optimisation = 0

if settings:
    blade10 = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=10)
    blade20 = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=20)
    blade30 = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=30)
    blade40 = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=40)
    blade50 = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=50)
    blade60 = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=60)

    J = np.linspace(0.1, 5, 200)

    kT10, kQ10, kP10, etaP10 = blade10.get_ks(J)
    kT20, kQ20, kP20, etaP20 = blade20.get_ks(J)
    kT30, kQ30, kP30, etaP30 = blade30.get_ks(J)
    kT40, kQ40, kP40, etaP40 = blade40.get_ks(J)
    kT50, kQ50, kP50, etaP50 = blade50.get_ks(J)
    kT60, kQ60, kP60, etaP60 = blade60.get_ks(J)

    plt.figure(figsize=(10,8))
    plt.subplot(2,2,1)
    plt.plot(J, kT10, label="Beta pitch 10°")
    plt.plot(J, kT20, label="Beta pitch 20°")
    plt.plot(J, kT30, label="Beta pitch 30°")
    plt.plot(J, kT40, label="Beta pitch 40°")
    plt.plot(J, kT50, label="Beta pitch 50°")
    plt.plot(J, kT60, label="Beta pitch 60°")
    plt.xlabel("J")
    plt.ylabel("kT")
    plt.legend()
    plt.grid()
    plt.subplot(2,2,2)
    plt.plot(J, kQ10, label="Beta pitch 10°")
    plt.plot(J, kQ20, label="Beta pitch 20°")
    plt.plot(J, kQ30, label="Beta pitch 30°")
    plt.plot(J, kQ40, label="Beta pitch 40°")
    plt.plot(J, kQ50, label="Beta pitch 50°")
    plt.plot(J, kQ60, label="Beta pitch 60°")
    plt.xlabel("J")
    plt.ylabel("kQ")
    plt.legend()
    plt.grid()
    plt.subplot(2,2,3)
    plt.plot(J, kP10, label="Beta pitch 10°")
    plt.plot(J, kP20, label="Beta pitch 20°")
    plt.plot(J, kP30, label="Beta pitch 30°")
    plt.plot(J, kP40, label="Beta pitch 40°")
    plt.plot(J, kP50, label="Beta pitch 50°")
    plt.plot(J, kP60, label="Beta pitch 60°")
    plt.xlabel("J")
    plt.ylabel("kP")
    plt.legend()
    plt.grid()
    plt.subplot(2,2,4)
    plt.plot(J, etaP10, label="Beta pitch 10°")
    plt.plot(J, etaP20, label="Beta pitch 20°")
    plt.plot(J, etaP30, label="Beta pitch 30°")
    plt.plot(J, etaP40, label="Beta pitch 40°")
    plt.plot(J, etaP50, label="Beta pitch 50°")
    plt.plot(J, etaP60, label="Beta pitch 60°")
    plt.xlabel("J")
    plt.ylabel("etaP")
    plt.legend()
    plt.tight_layout()
    plt.grid()
    plt.savefig("BEM_results_beta_pitch_10.png")
    plt.show()

if optimisation:
    blade = BEM(diameter=3.4, hub_diameter=0.45, num_blades=4, chord=0.25, beta_ref=15, beta_pitch=0)
    n = (3000/60) *0.477
    omega = n * 2 * pi
    p, t, rho = stdatm(20e3*0.3048)  # Standard atmosphere at 20000 ft
    M = 0.5
    u0 = M * sqrt(1.4 * 287 * t)
    J = 2*pi*u0 / (blade.diameter * omega)  # for 3000 RPM
    print(f"J: {J}, rho: {rho}, t: {t} K")

    def objective(beta_pitch):
        blade.beta_pitch = radians(beta_pitch)
        _, _, _, etaP = blade.get_ks(np.array([J]))
        if np.isnan(etaP[0]): return 1e6
        etaP[0] = min(etaP[0], 1.0)
        return -etaP[0]
    
    result = optimize.minimize_scalar(objective, bounds=(0, 90), method='bounded')
    optimal_beta_pitch = result.x
    blade.beta_pitch = radians(optimal_beta_pitch)
    kT_opt, kQ_opt, kP_opt, etaP_opt = blade.get_ks(np.array([J]))
    T = kT_opt[0] * rho * n**2 * blade.diameter**4
    Q = kQ_opt[0] * rho * n**2 * blade.diameter**5
    P = kP_opt[0] * rho * n**3 * blade.diameter**5
    
    print(f"Optimal beta pitch: {optimal_beta_pitch} degrees")
    print(f"At optimal beta pitch, kT: {kT_opt[0]}, kQ: {kQ_opt[0]}, kP: {kP_opt[0]}, etaP: {etaP_opt[0]}")
    print(f"Thrust: {T} N, Torque: {Q} Nm, Power: {P} W")