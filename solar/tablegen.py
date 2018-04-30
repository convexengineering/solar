from solar import Mission, Aircraft
from gpkit.repr_conventions import unitstr
from numpy import array

def var_table(sols, colns, varks, latns, title, label):
    assert len(sols) == len(colns) - 1

    filename = "./figs/%s.generated.tex" % label
    with open(filename, "w") as f:
        if title:
            f.write("\\caption{%s}\n" % title)
        f.write("\\begin{tabular}{lccccccccccccc}\n")
        f.write("\\toprule\n")
        f.write("\\toprule\n")
        f.write("\\label{t:%s}\n" % label)
        f.write(" & ".join(colns) + " \\\\\n")
        f.write("\\midrule\n")

        nvars = max(varks.shape)
        for i in range(nvars):
            vrs = varks[:, i]
            uts = unitstr(vrs[0].units) if vrs[0].units else ""
            if "**" in uts:
                uts = "$^".join(uts.split("**")) + "$"
            vals = [s(v).magnitude for s, v in zip(sols, vrs)]
            line = "$" + latns[i] + "$" + " & " + " & ".join(["%.3g " % x + uts for x in vals])
            f.write(line + "\\\\\n")
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}")

if __name__ == "__main__":

    V = Aircraft(sp=False)
    M = Mission(V, latitude=[20])
    M.cost = M[M.aircraft.Wtotal]

    # Baseline table
    sol = M.solve("mosek")
    CLNS = ["Design Variable", "Optimum Value"]
    VARS = array([[M.aircraft.Wtotal, M.aircraft.battery.W,
                  M.aircraft.wing.W, M.aircraft.wing.planform.AR,
                  M.aircraft.wing.planform.b, M.mission[1].fs.V,
                  M.mission[1].aircraftPerf.drag.wing.CL]])
    LATNS = ["W_{\\mathrm{MTO}}", "W_{\\mathrm{batt}}",
             "W_{\\mathrm{wing}}", "AR", "b", "V", "C_L"]
    var_table([sol], CLNS, VARS, LATNS, None, "vars")

    # twist trade
    V1 = Aircraft(sp=True)
    M1 = Mission(V1, latitude=[20])
    M1.cost = M1.aircraft.Wtotal
    M1.substitutions[M1.aircraft.emp.htail.Vh] = 0.45
    sol1 = M1.localsolve("mosek")

    CLNS = ["Variable", "Baseline Value", "With Twist Constraint"]
    VARS = array([[M.aircraft.Wtotal, M.aircraft.wing.W,
                   M.aircraft.wing.spar.W, M.aircraft.wing.planform.AR,
                   M.aircraft.wing.planform.b],
                  [M1.aircraft.Wtotal, M1.aircraft.wing.W,
                   M1.aircraft.wing.spar.W, M1.aircraft.wing.planform.AR,
                   M1.aircraft.wing.planform.b]])
    LATNS = ["W_{\\mathrm{MTO}}", "W_{\\mathrm{wing}}",
             "W_{\\mathrm{spar}}", "AR", "b"]
    var_table([sol, sol1], CLNS, VARS, LATNS,
              "Effect of Twist Constraint", "twist")

    V1 = Aircraft(sp=True)
    M1 = Mission(V1, latitude=[20])
    M1.cost = M1.aircraft.Wtotal
    M1.substitutions[M1.mission[1].wingg.twmax] = 10
    M1.substitutions[M1.mission[1].winggust.twmax] = 10
    sol1 = M1.localsolve("mosek")

    CLNS = ["Variable", "Baseline Value", "With Tail Flex"]
    VARS = array([[M.aircraft.Wtotal, M.aircraft.emp.htail.Vh,
                   M.aircraft.emp.htail.lh, M.aircraft.emp.htail.planform.S,
                   M.aircraft.emp.htail.W],
                  [M1.aircraft.Wtotal, M1.aircraft.emp.htail.Vh,
                   M1.aircraft.emp.htail.lh, M1.aircraft.emp.htail.planform.S,
                   M1.aircraft.emp.htail.W]])
    LATNS = ["W_{\\mathrm{MTOW}}", "V_{\\mathrm{h}}",
             "l_{\\mathrm{h}}", "S_{\\mathrm{h}}",
             "W_{\\mathrm{h}}"]
    var_table([sol, sol1], CLNS, VARS, LATNS,
              "Effect of tail boom flexibility on design variables", "tbflex")
