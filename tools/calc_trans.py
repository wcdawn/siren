import numpy as np
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

    sigma_t_max = {}
    sigma_t_min = {}
    for m in xs:
        sigma_t_max[m] = xs[m]["sigma_t"].max()
        sigma_t_min[m] = xs[m]["sigma_t"].min()
    print("Maximum Sigma_t")
    for m in sigma_t_max:
        print(m, "{:.6e}".format(sigma_t_max[m]))
    print("Minimum Sigma_t")
    for m in sigma_t_min:
        print(m, "{:.6e}".format(sigma_t_min[m]))

    fname_matmap = fname_phi.replace("phi", "matmap")
    xstart, xend, mat_map, mat_list = material_plot.load(fname_matmap)

    sigma_tr = np.zeros_like(phi)
    sigma_tr_clip = np.zeros_like(phi)

    for i in range(nx):
        material_name = mat_list[mat_map[i]]
        xsmat = xs[material_name]
        for n in range(pnorder):
            trans = xsmat["sigma_t"].copy()
            if n < nmoment:
                for g in range(ngroup):
                    ratio = phi[n, :, i] / phi[n, g, i]
                    trans[g] -= np.sum(xsmat["scatter"][n, g, :] * ratio)

            sigma_tr[n, :, i] = trans.copy()

            sigma_tr_clip[n, :, i] = np.clip(
                trans, sigma_t_min[material_name], sigma_t_max[material_name]
            )

    for n in range(pnorder):

        plt.figure()
        for g in range(ngroup):
            plt.plot(x, sigma_tr[n, g, :], label="g={:d}".format(g + 1))
        if ngroup <= 10:
            plt.legend()
        plt.xlim((0.0, 0.4))
        plt.ylim((-150, 150))
        plt.xlabel("x [cm]")
        plt.ylabel("$\\sigma_{{tr,g}}(x)$ [1/cm]")
        plt.title("CALC $\\sigma_{{tr,g}}$ {:d}".format(n))
        plt.tight_layout()
        plt.savefig("calc_trans_{:d}".format(n) + "." + extension, dpi=resolution)

        plt.figure()
        for g in range(ngroup):
            plt.plot(x, sigma_tr_clip[n, g, :], label="g={:d}".format(g + 1))
        if ngroup <= 10:
            plt.legend()
        plt.xlim((0.0, 0.4))
        plt.ylim((-150, 150))
        plt.xlabel("x [cm]")
        plt.ylabel("$\\sigma_{{tr,g}}(x)$ [1/cm]")
        plt.title("CALC CLIP $\\sigma_{{tr,g}}$ {:d}".format(n))
        plt.tight_layout()
        plt.savefig("calc_trans_clip_{:d}".format(n) + "." + extension, dpi=resolution)

    plt.show()
