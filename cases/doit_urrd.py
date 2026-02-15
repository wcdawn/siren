import numpy as np
import subprocess
import sys
import re

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
    pnorder = 17
    refine = 12

    cases = [
        "./URRd-H2Ob_1-2-0-ISLC/URRd-H2Ob_1-2-0-ISLC.inp",
        "./URRd-H2Ob_10-2-0-ISLC/URRd-H2Ob_10-2-0-ISLC.inp",
        "./URRd-H2Oc_1-2-0-ISLC/URRd-H2Oc_1-2-0-ISLC.inp",
        "./URRd-H2Oc_10-2-0-ISLC/URRd-H2Oc_10-2-0-ISLC.inp",
    ]

    ncases = len(cases)

    keff = {}

    for c in cases:

        fname_run = c.replace(".inp", "_run.inp")
        runtxt = open(c, "r").readlines()
        print(fname_run)
        inp = set_input(runtxt, pnorder, refine)
        open(fname_run, "w").writelines(inp)
        out = run(executable, fname_run)
        if "WARNING" in out:
            print("Encountered warning: ", c)
        out = out.split("\n")
        keff[c] = get_keff(out)

    print(keff)

    for c in keff:
        stub = re.sub("^.*/", "", c).replace(".inp", "")
        print(stub + " & " + " {:.6f}".format(keff[c]) + " & " + "{:.2f}".format((1.0-keff[c])*1e5) + " \\\\")
