import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    extension = "png"
    resolution = 600

    dat = np.loadtxt("albedo_power_t1.csv", delimiter=",", skiprows=1)
    x = dat[:, 0]
    power_t0 = dat[:, 1]

    dat = np.loadtxt("albedo_power_t2.csv", delimiter=",", skiprows=1)
    power_t1 = dat[:, 1]

    dat = np.loadtxt("albedo_power_t3.csv", delimiter=",", skiprows=1)
    power_t2 = dat[:, 1]

    dat = np.loadtxt("albedo_power_t4.csv", delimiter=",", skiprows=1)
    power_t3 = dat[:, 1]

    dat = np.loadtxt("albedo_power_t5.csv", delimiter=",", skiprows=1)
    power_t4 = dat[:, 1]

    plt.figure()
    plt.plot(x, power_t0, label="t=0s")
    plt.plot(x, power_t1, label="t=1s")
    plt.plot(x, power_t2, label="t=2s")
    plt.plot(x, power_t3, label="t=3s")
    plt.plot(x, power_t4, label="t=4s")
    plt.legend()
    plt.xlabel("x [cm]")
    plt.ylabel("Power [arb. units]")
    plt.title("Power Distribution During Transient")
    plt.tight_layout()
    plt.savefig("power." + extension, dpi=resolution)

    plt.show()
