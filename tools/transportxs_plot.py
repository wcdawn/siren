import numpy as np
import matplotlib.pyplot as plt
import sys


def load(fname):

    dat = np.loadtxt(fname, delimiter=",", skiprows=1)

    x = dat[:, 0]
    sigma_tr = dat[:, 1:]

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

        n = int(t[2].replace("n", ""))
        pnorder = np.max((pnorder, n))
    pnorder += 1

    sigma_tr_reshape = np.zeros((pnorder, ngroup, len(x)))
    for n in range(pnorder):
        for g in range(ngroup):
            sigma_tr_reshape[n, g, :] = sigma_tr[:, g + n * ngroup]

    return x, sigma_tr_reshape


if __name__ == "__main__":

    fname = sys.argv[1]
    extension = "png"
    dpi = 600

    x, sigma_tr = load(fname)
    pnorder = sigma_tr.shape[0]
    ngroup = sigma_tr.shape[1]

    for n in range(pnorder):

        plt.figure()
        for g in range(ngroup):
            plt.plot(x, sigma_tr[n, g, :], label="g={:d}".format(g + 1))
        if ngroup <= 10:
            plt.legend()
        plt.xlabel("x [cm]")
        plt.ylabel("$\\sigma_{{tr}}(x)$ [1/cm]")
        plt.title("Siren $\\sigma_{{tr}}$ {:d}".format(n))
        plt.tight_layout()
        plt.savefig("sigma_tr_n{:d}".format(n) + extension, dpi=dpi)

    plt.show()
