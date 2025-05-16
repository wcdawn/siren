import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    fname = "flux.csv"
    extension = "png"
    dpi = 600

    dat = np.loadtxt(fname, delimiter=",", skiprows=1)

    x = dat[:, 0]
    flux = dat[:, 1:]
    ngroup = flux.shape[1]

    plt.figure()
    for g in range(ngroup):
        plt.plot(x, flux[:, g], label="g={:d}".format(g))
    if ngroup <= 10:
        plt.legend()
    plt.xlabel("x [cm]")
    plt.ylabel("Flux (arb. units)")
    plt.title("SIREN Flux")
    plt.tight_layout()
    plt.savefig("flux." + extension, dpi=dpi)

    plt.show()
