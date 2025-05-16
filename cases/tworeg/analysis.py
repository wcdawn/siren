import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.optimize

if __name__ == "__main__":

    extension = "png"
    dpi = 600

    nusf = 0.02
    sigma_r_fuel = 0.02
    D_fuel = 1.2

    sigma_r_refl = 0.015
    D_refl = 0.7

    kappa_refl = np.sqrt(sigma_r_refl / D_refl)

    LF = 80.0
    LR = 100.0

    f = lambda bf: (D_refl * kappa_refl) / np.tanh(
        kappa_refl * (LR - LF)
    ) - D_fuel * bf * np.tan(bf * LF)
    BF = scipy.optimize.fsolve(f, 0.017)[0]
    print("BF={:.20f}".format(BF))

    keff = nusf / (D_fuel * BF**2 + sigma_r_fuel)
    print("keff={:.20f}".format(keff))

    phi0 = 1.0

    phi_fuel = lambda x: phi0 * np.cos(BF * x)
    phi_refl = (
        lambda x: phi0
        * np.cos(BF * LF)
        * (
            np.cosh(kappa_refl * (x - LF))
            - np.sinh(kappa_refl * (x - LF)) / np.tanh(kappa_refl * (LR - LF))
        )
    )

    def phi(x):
        p = np.zeros_like(x)
        for i in range(len(x)):
            if x[i] < LF:
                p[i] = phi_fuel(x[i])
            else:
                p[i] = phi_refl(x[i])
        return p

    x = np.linspace(0.0, LR, 1024)

    plt.figure()
    plt.plot(x, phi(x))
    plt.title("Flux")
    plt.xlabel("x [cm]")
    plt.ylabel("$\\phi(x)$")
    plt.tight_layout()
    plt.savefig("analytic_phi." + extension, dpi=dpi)

    plt.show()
