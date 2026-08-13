import numpy as np
import matplotlib.pyplot as plt
import sys

import plot_settings

if __name__ == "__main__":

    fname = sys.argv[1]
    extension = "png"
    dpi = 600

    dat = np.loadtxt(fname, delimiter=",", skiprows=1)

    elapt = dat[:, 0]
    relpow = dat[:, 2]

    plt.figure()
    plt.plot(elapt, relpow)
    plt.xlabel("Time [s]")
    plt.ylabel("Relative Power [-]")
    plt.title("Siren Transient")
    plt.tight_layout()
    plt.savefig("transient." + extension, dpi=600)

    plt.show()
