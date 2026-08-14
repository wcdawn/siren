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

    ref_time = np.array((0.000, 0.001, 0.002, 0.003, 0.0035, 0.004, 0.0045, 0.005))
    ref_powr = np.array(
        (
            1.000,
            1.178,
            1.558,
            2.797,
            5.284,
            20.72,
            471.6,
            153500.0,
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
