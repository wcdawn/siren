import numpy as np
import matplotlib.pyplot as plt
import sys

import plot_settings
import transportxs_plot
import xslib
from private_ebound import c5_586
from ebound import anl425


def compute_dx(x):
    dx = np.zeros_like(x)
    xleft = 0.0
    for i in range(len(dx)):
        dx[i] = 2.0 * (x[i] - xleft)
        xleft += dx[i]
    return dx


def compute_groupwise_average(dx, table):
    ngroup = table.shape[0]
    nx = table.shape[1]
    avg = np.zeros(ngroup)
    for i in range(nx):
        avg += dx[i] * table[:, i]
    avg /= np.sum(dx)
    return avg


if __name__ == "__main__":

    extension = "png"
    resolution = 600

    fname = sys.argv[1]  # transportxs
    fname_xslib = sys.argv[2]  # xs library
    hydrogen_name = sys.argv[3]  # name in xslib for hydrogen

    x, sigma_tr = transportxs_plot.load(fname)
    pnorder = sigma_tr.shape[0]
    ngroup = sigma_tr.shape[1]
    nx = sigma_tr.shape[2]

    dx = compute_dx(x)
    sigtr = compute_groupwise_average(dx, sigma_tr[1, :, :])

    egrid = anl425
    ebins_midpoint = np.zeros(ngroup)
    for g in range(ngroup):
        ebins_midpoint[g] = 0.5 * (egrid[g] + egrid[g + 1])

    xs = xslib.load(fname_xslib)
    sigt = xs[hydrogen_name]["sigma_t"].copy()

    plt.figure()
    plt.loglog(ebins_midpoint, sigtr, "-x")
    plt.xlabel("Energy [eV]")
    plt.ylabel("Cross Section [barn]")
    plt.title("Transport Cross Section")
    plt.tight_layout()
    plt.savefig("sigma_tr_spectrum." + extension, dpi=resolution)

    plt.figure()
    plt.loglog(ebins_midpoint, sigt, "-x")
    plt.xlabel("Energy [eV]")
    plt.ylabel("Cross Section [barn]")
    plt.title("Total Cross Section")
    plt.tight_layout()
    plt.savefig("sigma_t_spectrum." + extension, dpi=resolution)

    plt.figure()
    plt.semilogx(ebins_midpoint, sigtr / sigt, "-x", label="Siren")
    plt.xlabel("Energy [eV]")
    plt.ylabel("$\\Sigma_{tr} / \\Sigma_{t}$")
    plt.title("Transport-to-Total Cross Section Ratio")
    plt.tight_layout()
    plt.savefig("sigma_tr_ratio." + extension, dpi=resolution)

    plt.show()
