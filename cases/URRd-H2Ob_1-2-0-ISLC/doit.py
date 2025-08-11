import numpy as np
import subprocess
import matplotlib
import matplotlib.pyplot as plt
import sys

matplotlib.rcParams["lines.linewidth"] = 2
# matplotlib.rcParams["mathtext.fontset"] = "stix"
# matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams["font.size"] = 16


def set_input(txt, pnorder, refine):
    out = []
    for line in txt:
        if "refine" in line:
            out.append("refine {:d}\n".format(refine))
        elif "pnorder" in line:
            out.append("pnorder {:d}\n".format(pnorder))
        else:
            out.append(line)
    return out


def run(exe, inp):
    result = subprocess.run([exe, inp], stdout=subprocess.PIPE)
    return result.stdout.decode("utf-8")


def get_keff(lines):
    for line in lines:
        if "keff = " in line:
            line = line.split()
            return float(line[2])


def get_nx(lines):
    ready = False
    for line in lines:
        if "after refinement" in line:
            ready = True
        elif "nx =" in line:
            nx = int(line.split()[2])
            if ready:
                return nx
    # in case refinement not performed
    return nx


if __name__ == "__main__":

    executable = "/Users/williamdawn/work/siren/src/siren.x"
    fname_base = "URRd-H2Ob_1-2-0-ISLC.inp"
    max_pnorder = 13
    max_refine = 12

    fname_run = fname_base.replace(".inp", "_run.inp")

    runtxt = open(fname_base, "r").readlines()

    table = np.zeros((int((max_pnorder + 1) / 2), max_refine))
    nx = np.zeros(max_refine, dtype=int)

    pidx = 0
    for p in range(1, max_pnorder + 1, 2):
        ridx = 0
        for r in range(0, max_refine):
            inp = set_input(runtxt, p, r)
            open(fname_run, "w").writelines(inp)

            out = run(executable, fname_run)
            if "WARNING" in out:
                print("Encountered warning p=", p, "r=", r)
                sys.exit(1)
            out = out.split("\n")
            keff = get_keff(out)
            table[pidx, ridx] = keff

            # note: overwritten several times
            nx[ridx] = get_nx(out)

            print("p=", p, "r=", r, "keff= {:.16f}".format(keff))
            ridx += 1
        pidx += 1

    with open("result.csv", "w") as f:
        for pidx in range(table.shape[0]):
            f.write(" , p{:d}".format(pidx * 2 + 1))
        f.write("\n")
        for ridx in range(table.shape[1]):
            f.write("r{:d}".format(ridx))
            for pidx in range(table.shape[0]):
                f.write(" , {:.16f}".format(table[pidx, ridx]))
            f.write("\n")

    diff = 1.0 - table

    plt.figure()
    for pidx in range(table.shape[0]):
        plt.semilogx(
            nx, diff[pidx, :] * 1e5, "-o", label="$P_{{{:d}}}$".format(2 * pidx + 1)
        )
    plt.legend()
    plt.xlabel("NX")
    plt.ylabel("Error [pcm]")
    plt.tight_layout()
    plt.title("Spatial Refinement")

    plt.figure()
    for pidx in range(table.shape[0]):
        plt.semilogx(
            nx, table[pidx, :], "-o", label="$P_{{{:d}}}$".format(2 * pidx + 1)
        )
    plt.legend()
    plt.xlabel("NX")
    plt.ylabel("keff")
    plt.tight_layout()
    plt.title("Spatial Refinement")

    pnorder = np.zeros(table.shape[0])
    for pidx in range(table.shape[0]):
        pnorder[pidx] = 2 * pidx + 1

    plt.figure()
    for ridx in range(table.shape[1]):
        plt.plot(pnorder, diff[:, ridx] * 1e5, "-o", label="NX={:d}".format(nx[ridx]))
    plt.legend()
    plt.xlabel("$P_N$")
    plt.ylabel("Error [pcm]")
    plt.tight_layout()
    plt.title("Moment Refinement")

    plt.figure()
    for ridx in range(table.shape[1]):
        plt.plot(pnorder, table[:, ridx], "-o", label="NX={:d}".format(nx[ridx]))
    plt.legend()
    plt.xlabel("$P_N$")
    plt.ylabel("keff")
    plt.tight_layout()
    plt.title("Moment Refinement")

    plt.show()
