import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    fname = "phi.csv"
    extension = "png"
    dpi = 600

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

    print("g=", ngroup)
    print("pnorder=", pnorder)

    dphi = np.zeros_like(phi)
    hx = x[1]-x[0]
    for i in range(phi.shape[1]):
        dphi[:,i] = np.gradient(phi[:,i], hx)

    # post-compute the odd moments as well
    # TODO remove this
    dat = np.loadtxt("sigma_tr.csv", delimiter=",", skiprows=1)
    sigma_tr = dat[:,1:]
    phi_odd = np.zeros_like(phi)
    for n in range(1,pnorder,2):
        print("nodd=", n, "(n-1)=", n-1)
        for g in range(ngroup):
            if (n < pnorder-1):
                phi_odd[:,g + n*ngroup] = - (
                    (n+1)/(2*n+1) * dphi[:, g + (n+1)*ngroup]
                    + n/(2*n+1) * dphi[:, g + (n-1)*ngroup]
                )/sigma_tr[:,g+n*ngroup]
            else:
                phi_odd[:,g+n*ngroup] = -(
                    n/(2*n+1)*dphi[:,g+(n-1)*ngroup]
                )/sigma_tr[:,g+n*ngroup]

    for n in range(pnorder):

        plt.figure()
        for g in range(ngroup):
            plt.plot(x, phi[:, g + n * ngroup], label="g={:d}".format(g + 1))
        if ngroup <= 10:
            plt.legend()
        plt.xlabel("x [cm]")
        plt.ylabel("$\\phi(x)$ (arb. units)")
        plt.title("SIREN $\\phi$ {:d}".format(n))
        plt.tight_layout()
        plt.savefig("phi_{:d}".format(n) + extension, dpi=dpi)

        if (n%2 == 1):
            plt.figure()
            for g in range(ngroup):
                plt.plot(x, phi_odd[:, g + n * ngroup], label="g={:d}".format(g + 1))
            if ngroup <= 10:
                plt.legend()
            plt.xlabel("x [cm]")
            plt.ylabel("$\\phi(x)$ (arb. units)")
            plt.title("SIREN ODD COMPUTE $\\phi$ {:d}".format(n))
            plt.tight_layout()
            plt.savefig("phi_odd_{:d}".format(n) + extension, dpi=dpi)


    plt.show()
