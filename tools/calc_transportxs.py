import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.optimize
import phi_plot
import material_plot
import xslib
import plot_settings
import ebound


def energy_expand(energy_edge, xs):
    ngroup = len(xs)
    xx = np.zeros(ngroup * 2)
    yy = np.zeros(ngroup * 2)
    for i in range(ngroup):
        xx[2 * i + 0] = energy_edge[i]
        xx[2 * i + 1] = energy_edge[i + 1]
        dE = np.abs(energy_edge[i] - energy_edge[i + 1])  # in case backward order
        yy[2 * i + 0] = xs[i] / dE
        yy[2 * i + 1] = xs[i] / dE
    return xx, yy


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
        ratio = phi_g / phi_g[g]
        sigma_tr[g] = sigma_t[g] - np.sum(scatter[g, :] * ratio)
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

    maxiter = 1
    reltol = 1e-6

    for i in range(maxiter):
        sigma_tr_old = sigma_tr.copy()
        sigma_tr = eval_transportxs(sigma_t, scatter, phi_g)
        phi_g_old = phi_g.copy()
        phi_g = eval_phi_odd(n, sigma_tr, dphi_g_prev, dphi_g_next)
        conv = np.max(
            ((sigma_tr - sigma_tr_old) / sigma_tr, (phi_g - phi_g_old) / phi_g)
        )
        if conv < reltol:
            return sigma_tr, phi_g

    print("failed to converge conv={:.2e}".format(conv))
    return sigma_tr, phi_g


def nonlin_richardson(n, sigma_t, sigma_tr, scatter, phi_g, dphi_g_prev, dphi_g_next):

    maxiter = 10
    reltol = 1e-6

    omega_sigma_tr = 0.8
    omega_phi = 0.75

    for i in range(maxiter):

        sigma_tr_old = sigma_tr.copy()
        sigma_tr = eval_transportxs(sigma_t, scatter, phi_g)
        sigma_tr = sigma_tr_old + omega_sigma_tr * (sigma_tr - sigma_tr_old)

        phi_g_old = phi_g.copy()
        phi_g = eval_phi_odd(n, sigma_tr, dphi_g_prev, dphi_g_next)
        phi_g = phi_g_old + omega_phi * (phi_g - phi_g_old)

        conv = np.max(
            ((sigma_tr - sigma_tr_old) / sigma_tr, (phi_g - phi_g_old) / phi_g)
        )
        if conv < reltol:
            return sigma_tr, phi_g

    print("failed to converge conv={:.2e}".format(conv))
    return sigma_tr, phi_g


def nonlin_experiment(n, sigma_t, sigma_tr, scatter, phi_g, dphi_g_prev, dphi_g_next):

    maxiter = 1
    reltol = 1e-6

    omega_sigma_tr = 1.0
    omega_phi = 1.0

    sigma_tr = sigma_t - np.diag(scatter)

    for i in range(maxiter):

        phi_g_old = phi_g.copy()
        phi_g = eval_phi_odd(n, sigma_tr, dphi_g_prev, dphi_g_next)
        phi_g = phi_g_old + omega_phi * (phi_g - phi_g_old)

        sigma_tr_old = sigma_tr.copy()
        sigma_tr = eval_transportxs(sigma_t, scatter, phi_g)
        sigma_tr = sigma_tr_old + omega_sigma_tr * (sigma_tr - sigma_tr_old)

        conv = np.max(
            ((sigma_tr - sigma_tr_old) / sigma_tr, (phi_g - phi_g_old) / phi_g)
        )
        if conv < reltol:
            return sigma_tr, phi_g

    print("failed to converge conv={:.2e}".format(conv))
    return sigma_tr, phi_g


def func(n, sigma_t, scatter, dphi_g_prev, dphi_g_next, xvec):

    ngroup = len(sigma_t)
    sigma_tr = xvec[:ngroup]
    phi_g = xvec[ngroup:]

    sigma_tr_update = eval_transportxs(sigma_t, scatter, phi_g)
    phi_g_update = eval_phi_odd(n, sigma_tr, dphi_g_prev, dphi_g_next)

    xvec = np.zeros_like(xvec)
    xvec[:ngroup] = sigma_tr - sigma_tr_update
    xvec[ngroup:] = phi_g - phi_g_update
    return xvec


def nonlin_solve(n, sigma_t, sigma_tr, scatter, phi_g, dphi_g_prev, dphi_g_next):

    # initial guess
    sigma_tr = sigma_t
    phi_g = eval_phi_odd(n, sigma_tr, dphi_g_prev, dphi_g_next)
    x0 = np.append(sigma_tr, phi_g)

    f = lambda x: func(n, sigma_t, scatter, dphi_g_prev, dphi_g_next, x)
    xvec = scipy.optimize.fsolve(f, x0)

    ngroup = len(sigma_t)
    return xvec[:ngroup], xvec[ngroup:]


def nonlin_jfnk(n, sigma_t, sigma_tr, scatter, phi_g, dphi_g_prev, dphi_g_next):

    # initial guess
    sigma_tr = sigma_t
    phi_g = eval_phi_odd(n, sigma_tr, dphi_g_prev, dphi_g_next)
    x0 = np.append(sigma_tr, phi_g)

    f = lambda x: func(n, sigma_t, scatter, dphi_g_prev, dphi_g_next, x)
    # sol = scipy.optimize.root(f, x0, method="krylov", options={"disp": True, "maxiter": int(2e4)})
    sol = scipy.optimize.root(f, x0, method="krylov", options={"maxiter": int(2e4)})
    xvec = sol.x
    if not sol.success:
        print("failed to converge")

    ngroup = len(sigma_t)
    return xvec[:ngroup], xvec[ngroup:]


def calc_transportxs(x, dx, phi_in, mat_map, xs):

    phi = phi_in.copy()
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

                sigma_tr[n, :, i] = xsmat["sigma_t"] - np.diag(scat)
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

        if n % 2 == 1:
            # these are only valid for mirror bcs
            phi[n, :, 0] = phi[n, :, 1] * 0.5 * dx[0] / (dx[0] + 0.5 * dx[1])
            phi[n, :, -1] = phi[n, :, -2] * 0.5 * dx[-1] / (dx[-1] + 0.5 * dx[-2])

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

    return sigma_tr, phi


def calc_transportxs_noisy(x, dx, phi_in, mat_map, xs):

    phi = phi_in.copy()
    pnorder = phi.shape[0]
    ngroup = phi.shape[1]
    nx = phi.shape[2]

    nmoment = xs[list(xs.keys())[0]]["scatter"].shape[0]

    sigma_tr = np.zeros_like(phi)
    dphi = np.zeros_like(phi)
    for n in range(pnorder):
        for g in range(ngroup):
            dphi[n, g, :] = np.gradient(phi[n, g, :], x)

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
                    dphi_g_next = dphi[n + 1, :, i]

                dphi_g_prev = dphi[n - 1, :, i]

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

        if n % 2 == 1:
            # these are only valid for mirror bcs
            phi[n, :, 0] = phi[n, :, 1] * 0.5 * dx[0] / (dx[0] + 0.5 * dx[1])
            phi[n, :, -1] = phi[n, :, -2] * 0.5 * dx[-1] / (dx[-1] + 0.5 * dx[-2])

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

    return sigma_tr, phi


if __name__ == "__main__":

    extension = "pdf"
    resolution = 600

    fname_phi = "../cases/pin_slab/pin_slab_phi.csv"
    fname_xs = "../cases/pin_slab/c5xs.xs"  # TODO
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

    interfaces = []
    for i in range(1, len(mat_map)):
        if mat_map[i - 1] != mat_map[i]:
            interfaces.append(xstart[i])

    sigma_tr, phi_update = calc_transportxs(x, dx, phi, mat_map_name, xs)
    sigma_tr_noisy, phi_update_noisy = calc_transportxs_noisy(
        x, dx, phi, mat_map_name, xs
    )

    pnorder = sigma_tr.shape[0]
    ngroup = sigma_tr.shape[1]

    for n in range(pnorder):
        plt.figure()
        for g in range(ngroup):
            plt.plot(x, sigma_tr[n, g, :], label="g={:d}".format(g + 1))
        yl = plt.ylim()
        for xx in interfaces:
            plt.plot((xx, xx), yl, "-k", lw=1.0, label="_hide")
        plt.ylim(yl)
        if ngroup <= 10:
            plt.legend()
        plt.xlim((x[1], 0.4))
        plt.ylim((-1000, 1000))
        plt.xlabel("x [cm]")
        plt.ylabel("$\\Sigma_{{n,g}}(x)$ [1/cm]")
        plt.title("Transport Cross Section -- n={:d}".format(n))
        plt.tight_layout()
        plt.savefig("sigma_tr_n{:d}.".format(n) + extension, dpi=resolution)

    for n in range(pnorder):
        plt.figure()
        for g in range(ngroup):
            plt.plot(x, phi_update[n, g, :], label="g={:d}".format(g + 1))
        # yl = plt.ylim()
        # for xx in interfaces:
        #    plt.plot((xx, xx), yl, "-k", lw=1.0, label="_hide")
        # plt.ylim(yl)
        if ngroup <= 10:
            plt.legend()
        plt.xlabel("x [cm]")
        plt.ylabel("$\\phi_{{n,g}}(x)$ (arb. units)")
        plt.title("Flux Moment -- n={:d}".format(n))
        plt.tight_layout()
        plt.savefig("phi_update_n{:d}".format(n) + "." + extension, dpi=resolution)

    for n in range(pnorder):
        plt.figure()
        for g in range(ngroup):
            plt.plot(x, phi_update_noisy[n, g, :], label="g={:d}".format(g + 1))
        # yl = plt.ylim()
        # for xx in interfaces:
        #    plt.plot((xx, xx), yl, "-k", lw=1.0, label="_hide")
        # plt.ylim(yl)
        if ngroup <= 10:
            plt.legend()
        plt.xlabel("x [cm]")
        plt.ylabel("$\\phi_{{n,g}}(x)$ (arb. units)")
        plt.title("Flux Moment -- n={:d}".format(n))
        plt.tight_layout()
        plt.savefig(
            "phi_update_noisy_n{:d}".format(n) + "." + extension, dpi=resolution
        )

    for n in range(pnorder):
        plt.figure()
        for g in range(ngroup):
            plt.plot(x, phi[n, g, :] / phi[n, g, 0], label="g={:d}".format(g + 1))
        yl = plt.ylim()
        for xx in interfaces:
            plt.plot((xx, xx), yl, "-k", lw=1.0, label="_hide")
        plt.ylim(yl)
        if ngroup <= 10:
            plt.legend()
        plt.xlim((x[1], 0.4))
        plt.ylim((-500, 500))
        plt.plot((x[1], 0.4), (0, 0), "-k", lw=1, label="_hide")
        plt.xlabel("x [cm]")
        plt.ylabel("$\\phi_{{n,g}}(x) / \\phi_{{n,1}}(x)$")
        plt.title("Flux Ratios -- n={:d}".format(n))
        plt.tight_layout()
        plt.savefig("ratio_n{:d}.".format(n) + extension, dpi=resolution)

    for n in range(0, pnorder, 2):
        plt.figure()
        for g in range(ngroup):
            xnext = (n + 1) * (n + 1) / ((2 * n + 1) * (2 * n + 3))
            xprev = n * n / ((2 * n + 1) * (2 * n - 1))
            if n == 0:
                d = xnext / sigma_tr[n + 1, g, :]
            else:
                d = xnext / sigma_tr[n + 1, g, :] + xprev / sigma_tr[n - 1, g, :]
            plt.plot(x, d, label="g={:d}".format(g + 1))
        yl = plt.ylim()
        for xx in interfaces:
            plt.plot((xx, xx), yl, "-k", lw=1.0, label="_hide")
        plt.ylim(yl)
        plt.xlim((x[1], 0.4))
        plt.ylim((-0.1, 0.1))
        plt.plot((x[1], 0.4), (0, 0), "-k", lw=1.0, label="_hide")
        if ngroup <= 10:
            plt.legend()
        plt.xlabel("x [cm]")
        plt.ylabel("$\\hat{{D}}_{{n,g}}(x)$")
        plt.title("Pseudo-Diffusion Coefficient -- n={:d}".format(n))
        plt.tight_layout()
        plt.savefig("diff_n{:d}.".format(n) + extension, dpi=resolution)

    nx = phi.shape[2]
    spectra = {}
    width = {}
    for i in range(nx):
        mname = mat_map_name[i]
        if mname not in spectra:
            spectra[mname] = np.zeros((pnorder, ngroup))
            width[mname] = 0.0
        spectra[mname] += phi[:, :, i] * dx[i]
        width[mname] += dx[i]
    for m in spectra:
        spectra[m] /= width[m]

    for n in range(pnorder):
        plt.figure()
        for m in spectra:
            xx, yy = energy_expand(ebound.c5_586, spectra[m][n, :])
            plt.semilogx(xx, yy, label=m)
        xl = plt.xlim()
        plt.plot(xl, (0, 0), "-k", lw=1, label="_hide")
        plt.xlim(xl)
        plt.legend()
        plt.xlabel("Energy [eV]")
        plt.ylabel("Spectrum (arb. units)")
        plt.title("Flux Spectra -- n={:d}".format(n))
        plt.tight_layout()
        plt.savefig("spectra_n{:d}.".format(n) + extension, dpi=resolution)
        plt.close()

    plt.show()
