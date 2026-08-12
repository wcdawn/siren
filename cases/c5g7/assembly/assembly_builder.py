import sys

if __name__ == "__main__":

    mat_uo2 = ["moderator", "uo2", "moderator"]
    dx_uo2 = [0.09, 1.08, 0.09]

    mat_mox43 = ["moderator", "mox43", "moderator"]
    dx_mox43 = [0.09, 1.08, 0.09]

    mat_mox70 = ["moderator", "mox70", "moderator"]
    dx_mox70 = [0.09, 1.08, 0.09]

    mat_mox87 = ["moderator", "mox87", "moderator"]
    dx_mox87 = [0.09, 1.08, 0.09]

    mat_gdt = ["moderator", "guide_tube", "moderator"]
    dx_gdt = [0.09, 1.08, 0.09]

    assembly_uo2 = [
        "uo2",
        "uo2",
        "gdt",
        "uo2",
        "uo2",
        "gdt",
        "uo2",
        "uo2",
        "gdt",
        "uo2",
        "uo2",
        "gdt",
        "uo2",
        "uo2",
        "gdt",
        "uo2",
        "uo2",
    ]

    assembly_mox = [
        "mox43",
        "mox70",
        "gdt",
        "mox87",
        "mox87",
        "gdt",
        "mox87",
        "mox87",
        "gdt",
        "mox87",
        "mox87",
        "gdt",
        "mox87",
        "mox87",
        "gdt",
        "mox70",
        "mox43",
    ]

    dx = []
    mat = []
    for pin in assembly_uo2:
        if pin == "uo2":
            dx += dx_uo2
            mat += mat_uo2
        elif pin == "mox43":
            dx += dx_mox43
            mat += mat_mox43
        elif pin == "mox70":
            dx += dx_mox70
            mat += mat_mox70
        elif pin == "mox87":
            dx += dx_mox87
            mat += mat_mox87
        elif pin == "gdt":
            dx += dx_gdt
            mat += mat_gdt
        else:
            print("unknown pin: ", pin)
            sys.exit(1)

    print("nx ", len(dx))

    s = "dx"
    for x in dx:
        s += " {:.2f}".format(x)
    print(s)

    s = "mat_map"
    for m in mat:
        s += " " + m
    print(s)
