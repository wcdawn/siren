import numpy as np
import scipy
import scipy.linalg
import matplotlib.pyplot as plt
import sys

import xslib
import phi_plot
import material_plot
import plot_settings

if __name__ == "__main__":

    fname_phi = sys.argv[1]
    fname_xs = (
        "/Users/williamdawn/work/siren/cases/pin_slab/c5xs.xs"  # TODO command line?
    )

    extension = "png"
    resolution = 600

    x, phi = phi_plot.load(fname_phi)
    pnorder = phi.shape[0]
    ngroup = phi.shape[1]
    nx = phi.shape[2]

    xs = xslib.load(fname_xs)
    nmoment = xs[list(xs.keys())[0]]["scatter"].shape[0]

    fname_matmap = fname_phi.replace("phi", "matmap")
    xstart, xend, mat_map, mat_list = material_plot.load(fname_matmap)

    dphi = np.zeros_like(phi)
    for g in range(ngroup):
        for n in range(pnorder):
            dphi[n, g, :] = np.gradient(phi[n, g, :], x)

    for i in range(nx):
        material_name = mat_list[mat_map[i]]
        xsmat = xs[material_name]
        for n in range(1, pnorder, 2):

            if n < nmoment:
                trans = np.diag(xsmat["sigma_t"]) - xsmat["scatter"][n, :, :]
                itrans = scipy.linalg.inv(trans)
            else:
                trans = np.diag(xsmat["sigma_t"])
                itrans = np.diag(1.0 / xsmat["sigma_t"])

            xmul_next = (n + 1.0) / (2.0 * n + 1.0)
            xmul_prev = n / (2.0 * n + 1.0)

            dphi_prev = dphi[n - 1, :, i]
            if n + 1 < pnorder:
                dphi_next = dphi[n + 1, :, i]
            else:
                dphi_next = np.zeros(ngroup)

            phi[n, :, i] = -itrans @ (xmul_next * dphi_next + xmul_prev * dphi_prev)

    for n in range(pnorder):

        plt.figure()
        for g in range(ngroup):
            plt.plot(x, phi[n, g, :], label="g={:d}".format(g + 1))
        if ngroup <= 10:
            plt.legend()
        plt.ylim((-0.6, 0.4))
        plt.xlabel("x [cm]")
        plt.ylabel("$\\phi(x)$ (arb. units)")
        plt.title("Siren $\\phi$ {:d}".format(n))
        plt.tight_layout()
        plt.savefig("phi_sharp_{:d}".format(n) + "." + extension, dpi=resolution)

    plt.show()
