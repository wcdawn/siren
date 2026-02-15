import numpy as np
import matplotlib.pyplot as plt
import sys

import plot_settings
import phi_plot
import xslib
from private_ebound import c5_586


def energy_expand(energy_edge, xs):
    ngroup = len(xs)
    xx = np.zeros(ngroup * 2)
    yy = np.zeros(ngroup * 2)
    for i in range(ngroup):
        xx[2 * i + 0] = energy_edge[i]
        xx[2 * i + 1] = energy_edge[i + 1]
        dE = energy_edge[i] - energy_edge[i + 1]
        yy[2 * i + 0] = xs[i] / dE
        yy[2 * i + 1] = xs[i] / dE
    return xx, yy


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

    fname = sys.argv[1]

    x, phi = phi_plot.load(fname)
    pnorder = phi.shape[0]
    ngroup = phi.shape[1]
    nx = phi.shape[2]

    dx = compute_dx(x)
    phi0_spectrum = compute_groupwise_average(dx, phi[0, :, :])
    phi1_spectrum = compute_groupwise_average(dx, phi[1, :, :])

    c5_586_midpoint = np.zeros(ngroup)
    for g in range(ngroup):
        c5_586_midpoint[g] = 0.5 * (c5_586[g] + c5_586[g + 1])

    xs = xslib.load("c5xs_speng.xs")
    chi = xs["FUE_U235"]["chi"].copy()

    plt.figure()
    xx, yy = energy_expand(c5_586, phi0_spectrum)
    plt.loglog(xx, yy, label="$\\phi_0$")
    xx, yy = energy_expand(c5_586, phi1_spectrum)
    plt.loglog(xx, yy, label="$\\phi_1$")
    plt.legend()
    plt.xlabel("Energy [eV]")
    plt.ylabel("Arb. Units")
    plt.title("Transport Spectra")
    plt.tight_layout()
    plt.savefig("transport_spectra." + extension, dpi=resolution)

    plt.figure()
    xx, yy = energy_expand(c5_586, chi)
    plt.loglog(xx, yy, label="Source")
    plt.xlabel("Energy [eV]")
    plt.ylabel("Arb. Units")
    plt.title("Fixed-Source Spectra")
    plt.tight_layout()
    plt.savefig("transport_spectra." + extension, dpi=resolution)

    plt.show()
