import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys

matplotlib.rcParams["lines.linewidth"] = 2
matplotlib.rcParams["font.size"] = 14

if __name__ == "__main__":

    fname = "anl_slab_transient.csv"
    extension = "png"
    dpi = 600

    ref_time = np.array((0.000, 0.001, 0.005, 0.010, 0.012, 0.015, 0.018, 0.020))
    ref_powr = np.array(
        (
            1.000e0,
            1.022e0,
            1.659e0,
            1.565e1,
            7.019e1,
            6.803e2,
            6.612e3,
            3.012e4,
        )
    )

    dat = np.loadtxt(fname, delimiter=",", skiprows=1)

    elapt = dat[:, 0]
    relpow = dat[:, 2]

    plt.figure()
    plt.semilogy(elapt, relpow, label="Siren")
    plt.plot(ref_time, ref_powr, "or", label="Ref.")
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel("Relative Power [-]")
    plt.title("Siren Transient")
    plt.tight_layout()
    plt.savefig("transient." + extension, dpi=600)

    plt.show()
