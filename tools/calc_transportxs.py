import numpy as np
import matplotlib.pyplot as plt
import phi_plot
import material_plot
import xslib


def deriv(x1, x2, x3, f1, f2, f3):
    h_left = x2 - x1
    h_right = x3 - x2

    a = -(h_right**2) / (h_left * h_right * (h_left + h_right))
    c = h_left**2 / (h_left * h_right * (h_left + h_right))
    b = -(a + c)

    return a * f1 + b * f2 + c * f3


# evaluate all groups
# require a scattering matrix, it may be all zeros
def eval_transportxs(sigma_t, scatter, phi_g):
    sigma_tr = np.zeros_like(sigma_t)
    ngroup = len(sigma_tr)
    for g in range(ngroup):
        sigma_tr[g] = sigma_t[g] - np.sum(scatter[g, :] * phi_g) / phi_g[g]
    return sigma_tr


# evaluate all groups
# requires dphi_g_next, it may be all zeros
def eval_phi_odd(n, sigma_tr, dphi_g_prev, dphi_g_next):
    xmul_prev = n / (2.0 * n + 1.0)
    xmul_next = (n + 1.0) / (2.0 * n + 1.0)
    ngroup = len(sigma_tr)
    phi_g = np.zeros(ngroup)
    for g in range(ngroup):
        phi_g[g] = (
            -(xmul_prev * dphi_g_prev[g] + xmul_next * dphi_g_next[g]) / sigma_tr[g]
        )
    return phi_g


def nonlin_picard(n, sigma_t, sigma_tr, scatter, phi_g, dphi_g_prev, dphi_g_next):

    sigma_tr = eval_transportxs(sigma_t, scatter, phi_g)
    phi_g = eval_phi_odd(n, sigma_tr, dphi_g_prev, dphi_g_next)

    return sigma_tr, phi_g


def calc_transportxs(x, dx, phi, mat_map, xs):

    pnorder = phi.shape[0]
    ngroup = phi.shape[1]
    nx = phi.shape[2]

    nmoment = xs[list(xs.keys())[0]]["scatter"].shape[0]

    sigma_tr = np.zeros_like(phi)

    for i in range(1, nx - 1):
        xsmat = xs[mat_map[i]]

        for n in range(pnorder):

            if n >= nmoment:
                scat = np.zeros((ngroup, ngroup))
            else:
                scat = xsmat["scatter"][n, :, :]

            if n % 2 == 0:
                # simple update for even moments, they are not really used
                sigma_tr[n, :, i] = eval_transportxs(
                    xsmat["sigma_t"], scat, phi[n, :, i]
                )
            else:
                if n == pnorder - 1:
                    dphi_g_next = np.zeros(ngroup)
                else:
                    if (mat_map[i] == mat_map[i + 1]) and (
                        mat_map[i] == mat_map[i - 1]
                    ):
                        dphi_g_next = deriv(
                            -0.5 * (dx[i - 1] + dx[i]),
                            0.0,
                            +0.5 * (dx[i] + dx[i + 1]),
                            phi[n + 1, :, i - 1],
                            phi[n + 1, :, i],
                            phi[n + 1, :, i + 1],
                        )
                    elif mat_map[i] == mat_map[i + 1]:
                        dphi_g_next = (phi[n + 1, :, i + 1] - phi[n + 1, :, i]) / (
                            0.5 * (dx[i] + dx[i + 1])
                        )
                    elif mat_map[i] == mat_map[i - 1]:
                        dphi_g_next = (phi[n + 1, :, i] - phi[n + 1, :, i - 1]) / (
                            0.5 * (dx[i - 1] + dx[i])
                        )
                    else:
                        # fall back, try second-order
                        dphi_g_next = deriv(
                            -0.5 * (dx[i - 1] + dx[i]),
                            0.0,
                            +0.5 * (dx[i] + dx[i + 1]),
                            phi[n + 1, :, i - 1],
                            phi[n + 1, :, i],
                            phi[n + 1, :, i + 1],
                        )

                if (mat_map[i] == mat_map[i + 1]) and (mat_map[i] == mat_map[i - 1]):
                    dphi_g_prev = deriv(
                        -0.5 * (dx[i - 1] + dx[i]),
                        0.0,
                        +0.5 * (dx[i] + dx[i + 1]),
                        phi[n - 1, :, i - 1],
                        phi[n - 1, :, i],
                        phi[n - 1, :, i + 1],
                    )
                elif mat_map[i] == mat_map[i + 1]:
                    dphi_g_prev = (phi[n - 1, :, i + 1] - phi[n - 1, :, i]) / (
                        0.5 * (dx[i] + dx[i + 1])
                    )
                elif mat_map[i] == mat_map[i - 1]:
                    dphi_g_prev = (phi[n - 1, :, i] - phi[n - 1, :, i - 1]) / (
                        0.5 * (dx[i - 1] + dx[i])
                    )
                else:
                    # fall back, try second-order
                    dphi_g_prev = deriv(
                        -0.5 * (dx[i - 1] + dx[i]),
                        0.0,
                        +0.5 * (dx[i] + dx[i + 1]),
                        phi[n - 1, :, i - 1],
                        phi[n - 1, :, i],
                        phi[n - 1, :, i + 1],
                    )

                maxiter = 1
                reltol = 1e-6

                sigma_tr[n, :, i], phi[n, :, i] = nonlin_picard(
                    n,
                    xsmat["sigma_t"],
                    sigma_tr[n, :, i],
                    scat,
                    phi[n, :, i],
                    dphi_g_prev,
                    dphi_g_next,
                )

    for n in range(pnorder):

        xsmat = xs[mat_map[0]]
        if n >= nmoment:
            scat = np.zeros((ngroup, ngroup))
        else:
            scat = xsmat["scatter"][n, :, :]
        sigma_tr[n, :, 0] = eval_transportxs(xsmat["sigma_t"], scat, phi[n, :, 0])

        xsmat = xs[mat_map[-1]]
        if n >= nmoment:
            scat = np.zeros((ngroup, ngroup))
        else:
            scat = xsmat["scatter"][n, :, :]
        sigma_tr[n, :, -1] = eval_transportxs(xsmat["sigma_t"], scat, phi[n, :, -1])

        if n % 2 == 1:
            # these are only valid for mirror bcs
            phi[n, :, 0] = phi[n, :, 1] * 0.5 * dx[0] / (dx[0] + 0.5 * dx[1])
            phi[n, :, -1] = phi[n, :, -2] * 0.5 * dx[-1] / (dx[-1] + 0.5 * dx[-2])

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
