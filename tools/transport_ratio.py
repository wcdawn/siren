import numpy as np
import matplotlib.pyplot as plt
import sys

import plot_settings
import transportxs_plot
import xslib
from ebound import c5_586

if __name__ == "__main__":

    extension = "png"
    resolution = 600

    fname = sys.argv[1]

    x, sigma_tr = transportxs_plot.load(fname)

    dx = np.zeros_like(x)
    xleft = 0.0
    for i in range(len(dx)):
        dx[i] = 2.0 * (x[i] - xleft)
        xleft += dx[i]

    length = np.sum(dx)

    ngroup = sigma_tr.shape[1]
    nx = sigma_tr.shape[2]

    sigtr = np.zeros(ngroup)  # just first one...
    for i in range(nx):
        sigtr += dx[i] * sigma_tr[1, :, i]
    sigtr /= length

    c5_586_midpoint = np.zeros(ngroup)
    for g in range(ngroup):
        c5_586_midpoint[g] = 0.5 * (c5_586[g] + c5_586[g + 1])

    xs = xslib.load("c5xs_speng.xs")
    sigt = xs["COO_H1"]["sigma_t"].copy()

    plt.figure()
    plt.loglog(c5_586_midpoint, sigtr, "-x")
    plt.xlabel("Energy [eV]")
    plt.ylabel("Cross Section [barn]")
    plt.title("Transport Cross Section")
    plt.tight_layout()
    plt.savefig("sigma_tr_spectrum." + extension, dpi=resolution)

    plt.figure()
    plt.loglog(c5_586_midpoint, sigt, "-x")
    plt.xlabel("Energy [eV]")
    plt.ylabel("Cross Section [barn]")
    plt.title("Total Cross Section")
    plt.tight_layout()
    plt.savefig("sigma_t_spectrum." + extension, dpi=resolution)

    sigtr_ref = np.loadtxt("/Users/williamdawn/work/siren/cases/pno.slab/sig_tr.txt")

    plt.figure()
    plt.semilogx(c5_586_midpoint, sigtr / sigt, "-x", label="Siren")
    plt.semilogx(c5_586_midpoint, sigtr_ref / sigt, "-x", label="SPENG")
    plt.legend()
    plt.xlabel("Energy [eV]")
    plt.ylabel("$\\Sigma_{tr} / \\Sigma_{t}$")
    plt.title("Transport-to-Total Cross Section Ratio")
    plt.tight_layout()
    plt.savefig("sigma_tr_ratio." + extension, dpi=resolution)

    plt.show()
