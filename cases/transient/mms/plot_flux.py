import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    extension = "png"
    resolution = 600

    flux = []
    for i in range(5):
        dat = np.loadtxt("mms_flux_t{:d}.csv".format(i + 1), delimiter=",", skiprows=1)
        x = dat[:, 0]
        flux.append(dat[:, 1])

    tedit = [0.0, 0.25, 0.5, 0.75, 1.0]
    plt.figure()
    for i in range(5):
        plt.plot(x, flux[i], label="t={:.2f}s".format(tedit[i]))
    plt.legend()
    plt.xlabel("x [cm]")
    plt.ylabel("Flux [arb. units]")
    plt.title("Thermal Flux During Transient")
    plt.tight_layout()
    plt.savefig("flux." + extension, dpi=resolution)

    plt.show()
