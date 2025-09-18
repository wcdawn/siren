import numpy as np
import matplotlib.pyplot as plt
import sys


def load(fname):

    with open(fname, "r") as f:

        line = f.readline()
        s = line.split()
        if s[0] != "niso":
            print("Failure to parse mat_map")
            sys.exit(1)
        niso = int(s[1])

        mat_list = []
        for i in range(niso):
            line = f.readline().strip()
            mat_list.append(line)

    dat = np.loadtxt(fname, delimiter=",", skiprows=niso + 2)
    xstart = dat[:, 0]
    xend = dat[:, 1]
    mat_map = dat[:, 2].astype(int) - 1

    return xstart, xend, mat_map, mat_list


if __name__ == "__main__":

    fname = sys.argv[1]
    extension = "png"
    resolution = 600

    xstart, xend, mat_map, mat_list = load(fname)

    print("mat_list: ", mat_list)

    xx = np.zeros(2 * len(xstart))
    yy = np.zeros(2 * len(xstart), dtype=int)
    for i in range(len(xstart)):
        xx[2 * i] = xstart[i]
        xx[2 * i + 1] = xend[i]
        yy[2 * i] = mat_map[i]
        yy[2 * i + 1] = mat_map[i]

    plt.figure()
    plt.plot(xx, yy)
    plt.xlabel("x [cm]")
    plt.ylabel("Material ID")
    plt.title("Siren Material Map")
    plt.tight_layout()
    plt.savefig(fname.replace(".csv", "." + extension), dpi=resolution)
    plt.show()
