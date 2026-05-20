import numpy as np
import subprocess
import sys

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

    executable = "/Users/williamdawn/work/siren/src/siren.x"
    fname_base = "assembly_mix_uniform.inp"
    start_pnorder = 7
    start_refine = 0
    max_refine = 7

    fname_run = fname_base.replace(".inp", "_run.inp")

    runtxt = open(fname_base, "r").readlines()

    keff = np.zeros(max_refine)
    nx = np.zeros(max_refine, dtype=int)

    start_neq = int((start_pnorder+1)/2)

    for i in range(max_refine):
        neq = start_neq * 2 **i
        p = neq*2 - 1
        r = start_refine + i

        inp = set_input(runtxt, p, r)
        open(fname_run, "w").writelines(inp)

        out = run(executable, fname_run)
        if "WARNING" in out:
            print("Encountered warning p=", p, "r=", r)
            sys.exit(1)
        out = out.split("\n")
        keff[i] = get_keff(out)

        # note: overwritten several times
        nx[i] = get_nx(out)

        print("p=", p, "r=", r, "keff= {:.16f}".format(keff[i]))

    for k in keff:
        print("{:.16f}".format(k))
