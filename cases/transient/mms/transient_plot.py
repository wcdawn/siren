import numpy as np
import matplotlib.pyplot as plt
import sys

# import plot_settings
import matplotlib

matplotlib.rcParams["lines.linewidth"] = 2
matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams["font.size"] = 14

L = 100.0
alpha = 1e9
kappa_nu = 1.5e-11
nusf = 0.156
kappasf = kappa_nu * nusf
p0 = 1e-3
integral = 163.0 / 225.0 * L
phi0 = p0 / kappasf * np.pi * 0.5 / L


def relpow_exact(t):
    return 1.0 + np.pi * alpha * integral * 0.5 / L / phi0 * t**2


if __name__ == "__main__":

    # fname = sys.argv[1]
    fname = "mms_transient.csv"

    extension = "png"
    dpi = 600

    dat = np.loadtxt(fname, delimiter=",", skiprows=1)

    elapt = dat[:, 0]
    relpow = dat[:, 2]

    plt.figure()
    plt.plot(elapt, relpow_exact(elapt), label="Exact")
    plt.plot(elapt, relpow, label="Siren")
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel("Relative Power [-]")
    plt.title("Siren Transient")
    plt.tight_layout()
    plt.savefig("transient." + extension, dpi=600)

    plt.show()
