import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    extension = "png"
    resolution = 600

    dat = np.loadtxt("anl_slab_flux_t1.csv", delimiter=",", skiprows=1)
    x = dat[:, 0]
    thermal_flux_t0 = dat[:, 2]

    dat = np.loadtxt("anl_slab_flux_t2.csv", delimiter=",", skiprows=1)
    thermal_flux_t1 = dat[:, 2]

    dat = np.loadtxt("anl_slab_flux_t3.csv", delimiter=",", skiprows=1)
    thermal_flux_t2 = dat[:, 2]

    plt.figure()
    plt.plot(x, thermal_flux_t0, label="t=0s")
    plt.plot(x, thermal_flux_t1, label="t=1s")
    plt.plot(x, thermal_flux_t2, label="t=2s")
    plt.legend()
    plt.xlabel("x [cm]")
    plt.ylabel("Flux [arb. units]")
    plt.title("Thermal Flux During Transient")
    plt.tight_layout()
    plt.savefig("flux." + extension, dpi=resolution)

    plt.show()
