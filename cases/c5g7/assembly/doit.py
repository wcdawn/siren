import numpy as np
import subprocess
import matplotlib
import matplotlib.pyplot as plt
import sys

matplotlib.rcParams["lines.linewidth"] = 2
matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams["font.size"] = 11.5


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

    extension = "pdf"
    resolution = 600

    if len(sys.argv) > 1:
        fname_result = sys.argv[1]
        quick = True
    else:
        quick = False

    executable = "/Users/williamdawn/work/siren/src/siren.x"
    fname_base = "assembly_uo2_uniform.inp"
    arr_pnorder = [1, 3, 5, 9, 17, 33, 65, 129, 257, 513]
    max_refine = 7

    fname_run = fname_base.replace(".inp", "_run.inp")

    runtxt = open(fname_base, "r").readlines()

    table = np.zeros((len(arr_pnorder), max_refine))

    if quick:
        table = np.loadtxt(fname_result, delimiter=",", dtype=str)
        max_refine = int(table[-1, 0].replace("r", "")) + 1
        max_pnorder = int(table[0, -1].replace("p", ""))
        table = table[1:, 1:].astype(float)
        nx = 7 * np.ones(max_refine, dtype=int)
        for i in range(len(nx)):
            nx[i] *= 2**i
        table = table.transpose()
    else:
        nx = np.zeros(max_refine, dtype=int)
        for pidx in range(len(arr_pnorder)):
            p = arr_pnorder[pidx]
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

        with open("result.csv", "w") as f:
            for pidx in range(table.shape[0]):
                f.write(" , p{:d}".format(arr_pnorder[pidx]))
            f.write("\n")
            for ridx in range(table.shape[1]):
                f.write("r{:d}".format(ridx))
                for pidx in range(table.shape[0]):
                    f.write(" , {:.16f}".format(table[pidx, ridx]))
                f.write("\n")

    # perform Richardson extrapolation
    extrap = np.zeros(table.shape[0])
    for i in range(table.shape[0]):
        extrap[i] = table[i, -1] + (table[i, -1] - table[i, -2]) / 3.0
    # diff is the difference
    diff = np.zeros_like(table)
    for i in range(table.shape[0]):
        diff[i, :] = extrap[i] - table[i, :]

    # dump table in LaTeX format
    with open("pin_slab_result.tex", "w") as f:
        for i in range(table.shape[0]):
            f.write(" & {{$P_{{{:d}}}$}}".format(2*i+1))
        f.write("\\\\\n")
        f.write("\\midrule\n")
        for j in range(table.shape[1]):
            f.write("r{:d}".format(j))
            for i in range(table.shape[0]):
                f.write(" & {:.6f}".format(table[i,j]))
            f.write("\\\\\n")
        f.write("\\midrule\n")
        f.write("Extrap.")
        for e in extrap:
            f.write(" & {:.6f}".format(e))
        f.write("\\\\\n")

    plt.figure()
    for pidx in range(table.shape[0]):
        plt.semilogx(
            nx, table[pidx, :], "-o", label="$P_{{{:d}}}$".format(arr_pnorder[pidx])
        )
    plt.legend()
    plt.xlabel("NX")
    plt.ylabel("keff")
    plt.title("Spatial Refinement")
    plt.tight_layout()
    plt.savefig("refinement_space." + extension, dpi = resolution)

    plt.figure()
    for pidx in range(table.shape[0]):
        plt.loglog(
            nx,
            np.abs(diff[pidx, :]) * 1e5,
            "-o",
            label="$P_{{{:d}}}$".format(arr_pnorder[pidx]),
        )
    plt.loglog(nx, 1e5 * nx.astype(float) ** (-2), "-k", lw=1, label="_hide")
    plt.legend()
    plt.xlabel("NX")
    plt.ylabel("|extrapolation - keff| [pcm]")
    plt.title("Spatial Refinement to Extrapolation")
    plt.tight_layout()
    plt.savefig("refinement_conv." + extension, dpi = resolution)

    plt.figure()
    for ridx in range(table.shape[1]):
        plt.plot(arr_pnorder, table[:, ridx], "-o", label="NX={:d}".format(nx[ridx]))
    plt.legend()
    plt.xlabel("$P_N$")
    plt.gca().set_xticks(arr_pnorder)
    plt.ylabel("keff")
    plt.title("Moment Refinement")
    plt.tight_layout()
    plt.savefig("refinement_moment." + extension, dpi = resolution)

    plt.show()
