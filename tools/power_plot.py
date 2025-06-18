import numpy as np
import matplotlib.pyplot as plt
import sys

if __name__ == "__main__":

    fname = sys.argv[1]
    extension = "png"
    dpi = 600

    dat = np.loadtxt(fname, delimiter=",", skiprows=1)

    x = dat[:, 0]
    power = dat[:, 1:]

    plt.figure()
    plt.plot(x, power)
    plt.xlabel("x [cm]")
    plt.ylabel("Power (arb. units)")
    plt.title("SIREN Power")
    plt.tight_layout()
    plt.savefig("power." + extension, dpi=dpi)

    plt.show()
