import numpy as np
import matplotlib.pyplot as plt
import sys

if __name__ == "__main__":

    fname = sys.argv[1]
    extension = "png"
    resolution = 600

    dat = np.loadtxt(fname, delimiter=",", skiprows=1)

    with open(fname, "r") as f:
        line = f.readline()
    line = line[:-1]  # newline character
    s = line.split(",")
    s = s[1:]  # remove x [cm] label
    ngroup = -1
    pnorder = -1
    for title in s:
        t = title.split("_")

        g = int(t[-1].replace("g", ""))
        ngroup = np.max((ngroup, g))

        n = int(t[1].replace("n", ""))
        pnorder = np.max((pnorder, n))
    pnorder += 1

    arr_per_phi = int((dat.shape[1] - 1) / 3)

    x = dat[:, 0]
    phi_siren = dat[:, arr_per_phi * 0 + 1 : arr_per_phi * 1 + 1]
    phi_exact = dat[:, arr_per_phi * 1 + 1 : arr_per_phi * 2 + 1]
    phi_diff = dat[:, arr_per_phi * 2 + 1 : arr_per_phi * 3 + 1]

    for n in range(pnorder):
        for g in range(ngroup):
            plt.figure()
            plt.plot(x, phi_siren[:, g + n * ngroup], label="Siren")
            plt.plot(x, phi_exact[:, g + n * ngroup], label="Exact")
            plt.legend()
            plt.xlabel("x [cm]")
            plt.ylabel("$\\phi(x)$")
            plt.title("$\\phi(x)$ -- n={:d} -- g={:d}".format(n, g + 1))
            plt.tight_layout()
            plt.savefig("analytic_n{:d}_g{:d}.".format(n, g + 1), dpi=resolution)

            plt.figure()
            plt.plot(x, phi_diff[:, g + n * ngroup], label="Difference")
            plt.legend()
            plt.xlabel("x [cm]")
            plt.ylabel("$\\phi(x)$")
            plt.title("$\\phi(x)$ Diff. -- n={:d} -- g={:d}".format(n, g + 1))
            plt.tight_layout()
            plt.savefig("analytic_diff_n{:d}_g{:d}.".format(n, g + 1), dpi=resolution)

    plt.show()
