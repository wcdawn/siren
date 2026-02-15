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


def get_keyword(lines, keyword):
    for line in lines:
        if keyword + " = " in line:
            line = line.split()
            return float(line[2])


def get_keff(lines):
    return get_keyword(lines, "keff")


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
    fname_base = "tworeg.inp"
    fname_base = "tworeg_nonuniform.inp"
    max_refine = 11

    fname_run = fname_base.replace(".inp", "_run.inp")

    runtxt = open(fname_base, "r").readlines()

    linferr = np.zeros(max_refine)
    keff_diff = np.zeros(max_refine)
    nx = np.zeros(max_refine, dtype=int)

    p = 0
    for r in range(0, max_refine):
        inp = set_input(runtxt, p, r)
        open(fname_run, "w").writelines(inp)

        out = run(executable, fname_run)
        if "CONVERGENCE!" not in out:
            print("failed to converge p=", p, "r=", r)
            sys.exit(1)
        out = out.split("\n")

        keff = get_keff(out)
        keff_diff[r] = get_keyword(out, "keff_diff")

        linferr[r] = get_keyword(out, "linferr_n0_g1")

        nx[r] = get_nx(out)

        print("p=", p, "r=", r, "keff= {:.16f}".format(keff))

    with open("result.csv", "w") as f:
        f.write(" , keff_diff [pcm]")
        f.write(" , linferr")
        f.write("\n")
        for ridx in range(max_refine):
            f.write("r{:d}".format(ridx))
            f.write(" , {:.16f}".format(keff_diff[ridx]))
            f.write(" , {:.16e}".format(linferr[ridx]))
            f.write("\n")

    plt.figure()
    plt.loglog(nx, np.abs(keff_diff), "-o")
    plt.xlabel("NX")
    plt.ylabel("keff error [pcm]")
    plt.title("Spatial Refinement")
    plt.tight_layout()

    plt.figure()
    plt.loglog(nx, linferr, "-o")
    plt.xlabel("NX")
    plt.ylabel("$\\| \\phi \\|$")
    plt.title("Spatial Refinement")
    plt.tight_layout()

    plt.show()
