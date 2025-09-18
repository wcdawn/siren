import numpy as np
import matplotlib.pyplot as plt
import phi_plot
import material_plot
import xslib


def calc_phi_odd(x, dx, phi, mat_map, xs):

    # TODO this uses second-order everywhere
    # may be too numerically noisy,
    # consider using first-order at material interfaces
    dphi = np.zeros_like(phi)
    for n in range(pnorder):
        for g in range(ngroup):
            dphi[n, g, :] = np.gradient(phi[n, g, :], x)


def calc_transportxs(x, dx, phi, mat_map, xs):

    pnorder = phi.shape[0]
    ngroup = phi.shape[1]
    nx = phi.shape[2]

    nmoment = xs[list(xs.keys())[0]]["scatter"].shape[0]

    sigma_tr = np.zeros_like(phi)

    for n in range(pnorder):
        if n >= nmoment:
            for g in range(ngroup):
                for i in range(nx):
                    xsmat = xs[mat_map[i]]
                    sigma_tr[n, g, i] = xsmat["sigma_t"][g]
        else:
            for g in range(ngroup):
                for i in range(nx):
                    xsmat = xs[mat_map[i]]
                    sigma_tr[n, g, i] = (
                        xsmat["sigma_t"][g]
                        - np.sum(xsmat["scatter"][n, g, :] * phi[n, :, i])
                        / phi[n, g, i]
                    )

    return sigma_tr


if __name__ == "__main__":

    extension = "png"
    resolution = 600

    fname_phi = "../cases/pin_slab/pin_slab_phi.csv"
    fname_xs = "../cases/pin_slab/c5xs_p0.xs"  # TODO
    fname_matmap = "../cases/pin_slab/pin_slab_matmap.csv"

    x, phi = phi_plot.load(fname_phi)
    xs = xslib.load(fname_xs)

    # read the mat_map and convert to strings that we can use in the xslib dictionary
    xstart, xend, mat_map, mat_list = material_plot.load(fname_matmap)
    mat_map_name = []
    for m in mat_map:
        mat_map_name.append(mat_list[m])
    dx = xend.copy() - xstart

    sigma_tr = calc_transportxs(x, dx, phi, mat_map_name, xs)

    pnorder = sigma_tr.shape[0]
    ngroup = sigma_tr.shape[1]

    for n in range(pnorder):
        plt.figure()
        for g in range(ngroup):
            plt.plot(x, sigma_tr[n, g, :], label="g={:d}".format(g + 1))
        if ngroup <= 10:
            plt.legend()
        plt.xlabel("x [cm]")
        plt.ylabel("$\\Sigma_{{tr}}(x)$ (arb. units)")
        plt.title("Siren $\\Sigma_{{tr}}$ {:d}".format(n))
        plt.tight_layout()
        plt.savefig("sigma_tr_n{:d}.".format(n) + extension, dpi=resolution)
    plt.show()
