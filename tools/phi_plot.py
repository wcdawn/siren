import numpy as np
import matplotlib.pyplot as plt
import sys


def load(fname):

    dat = np.loadtxt(fname, delimiter=",", skiprows=1)

    x = dat[:, 0]
    phi = dat[:, 1:]

    with open(fname, "r") as f:
        line = f.readline()
    line = line[:-1]  # newline character
    s = line.split(",")
    s = s[1:]  # remove "x [cm]" label
    ngroup = -1
    pnorder = -1
    for title in s:
        t = title.split("_")

        g = int(t[-1].replace("g", ""))
        ngroup = np.max((ngroup, g))

        n = int(t[1].replace("n", ""))
        pnorder = np.max((pnorder, n))
    pnorder += 1

    phi_reshape = np.zeros((pnorder, ngroup, len(x)))
    for n in range(pnorder):
        for g in range(ngroup):
            phi_reshape[n, g, :] = phi[:, g + n * ngroup]

    return x, phi_reshape


if __name__ == "__main__":

    fname = sys.argv[1]
    extension = "png"
    resolution = 600

    x, phi = load(fname)
    pnorder = phi.shape[0]
    ngroup = phi.shape[1]

    dphi = np.zeros_like(phi)
    for g in range(ngroup):
        for n in range(pnorder):
            dphi[n, g, :] = np.gradient(phi[n, g, :], x)

    for n in range(pnorder):

        plt.figure()
        for g in range(ngroup):
            plt.plot(x, phi[n, g, :], label="g={:d}".format(g + 1))
        if ngroup <= 10:
            plt.legend()
        plt.xlabel("x [cm]")
        plt.ylabel("$\\phi(x)$ (arb. units)")
        plt.title("Siren $\\phi$ {:d}".format(n))
        plt.tight_layout()
        plt.savefig("phi_{:d}".format(n) + "." + extension, dpi=resolution)

    plt.show()
