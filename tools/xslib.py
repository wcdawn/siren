import numpy as np
import matplotlib.pyplot as plt


def read(fname):

    xs = {}
    xsthis = None

    with open(fname, "r") as f:
        while True:

            line = f.readline()

            if line == "":
                # EOF
                break
            elif line == "\n":
                # empty line
                continue

            s = line.split()

            if s[0] == "ngroup":
                ngroup = int(s[1])
            elif s[0] == "nmoment":
                nmoment = int(s[1])
            elif s[0] == "name":
                if xsthis is not None:
                    # insert into dictionary
                    xs[name] = xsthis
                xsthis = {}
                name = s[1]
            elif s[0] == "sigma_a":
                sigma_a = np.zeros(ngroup)
                for g in range(ngroup):
                    sigma_a[g] = float(f.readline())
                xsthis["sigma_a"] = sigma_a.copy()
            elif s[0] == "diffusion":
                diffusion = np.zeros(ngroup)
                for g in range(ngroup):
                    diffusion[g] = float(f.readline())
                xsthis["diffusion"] = diffusion.copy()
            elif s[0] == "sigma_t":
                sigma_t = np.zeros(ngroup)
                for g in range(ngroup):
                    sigma_t[g] = float(f.readline())
                xsthis["sigma_t"] = sigma_t.copy()
            elif s[0] == "sigma_f":
                sigma_f = np.zeros(ngroup)
                for g in range(ngroup):
                    sigma_f[g] = float(f.readline())
                xsthis["sigma_f"] = sigma_f.copy()
            elif s[0] == "nusf":
                nusf = np.zeros(ngroup)
                for g in range(ngroup):
                    nusf[g] = float(f.readline())
                xsthis["nusf"] = nusf.copy()
            elif s[0] == "chi":
                chi = np.zeros(ngroup)
                for g in range(ngroup):
                    chi[g] = float(f.readline())
                xsthis["chi"] = chi.copy()
            elif s[0] == "scatter":
                if "scatter" not in xsthis:
                    xsthis["scatter"] = np.zeros((nmoment + 1, ngroup, ngroup))
                ell = int(s[1])
                for g in range(ngroup):
                    xsthis["scatter"][ell, g, :] = np.array(
                        f.readline().split()
                    ).astype(float)

    if xsthis is not None:
        xs[name] = xsthis

    return xs


def eigen(xsmat):
    ngroup = len(xsmat["sigma_t"])
    if "nusf" not in xsmat:
        return 0.0, np.zeros(ngroup)

    keff = 1.0
    phi = np.ones(ngroup)

    A = np.diag(xsmat["sigma_t"])
    A -= xsmat["scatter"][0, :, :]

    F = np.outer(xsmat["chi"], xsmat["nusf"])

    M = np.linalg.inv(A) @ F

    eigval, eigvec = np.linalg.eig(M)

    kinf = np.max(np.real(eigval))

    # find the matching eigenvector
    for i in range(len(eigval)):
        if np.real(eigval[i]) == keff:
            phi = eigvec[:i]
            break

    phi = np.real(phi)
    phi = np.abs(phi)
    phi /= np.max(phi)  # bit of arbitrary normalization

    return kinf, phi


def summary(lib):

    materials = list(lib.keys())

    niso = len(materials)
    print("nmaterial: {:d}".format(niso))

    ngroup = len(lib[materials[0]]["sigma_t"])
    print("ngroup: {:d}".format(ngroup))

    nmoment = lib[materials[0]]["scatter"].shape[0] - 1
    print("nmoment: {:d}".format(nmoment))

    s = "materials:"
    for mat in materials:
        s += " " + mat
    print(s)

    for mat in lib:
        if "nusf" not in lib[mat]:
            print("Material " + mat + " is not fissile.")
        else:
            kinf, phi = eigen(lib[mat])
            print("Material " + mat + " kinf = {:.8f}".format(kinf))


if __name__ == "__main__":

    fname = "../cases/pin_slab/c5xs.xs"
    xs = read(fname)
    summary(xs)
