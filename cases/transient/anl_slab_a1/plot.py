import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams["lines.linewidth"] = 2
matplotlib.rcParams["font.size"] = 14

if __name__ == "__main__":

    fname = "anl_slab_transient.csv"
    extension = "png"
    dpi = 600

    ref_time = np.array((0.0, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0))
    ref_powr = np.array(
        (
            1.0000,
            0.9299,
            0.8733,
            0.7597,
            0.6588,
            0.6432,
            0.6307,
        )
    )

    dat = np.loadtxt(fname, delimiter=",", skiprows=1)

    elapt = dat[:, 0]
    relpow = dat[:, 2]

    plt.figure()
    plt.plot(elapt, relpow, label="Siren")
    plt.plot(ref_time, ref_powr, "or", label="Ref.")
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel("Relative Power [-]")
    plt.title("Siren Transient")
    plt.tight_layout()
    plt.savefig("transient." + extension, dpi=600)

    plt.show()
