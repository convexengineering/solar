" contour plots "
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
from solar import Mission, Aircraft
from plotting import labelLines
from gpkit.tools.autosweep import autosweep_1d

N = 100
plt.rcParams.update({'font.size':15})

def battsolarcon(latitudes):
    time = 0
    numsolves = 0
    etamax = []
    lines = []
    midx = []
    data = {}
    passing = True
    hmax = 500
    etamax = 0.4
    for lat in latitudes:
        if passing:
            V = Aircraft(sp=True)
            M = Mission(V, latitude=[lat])
            del M.substitutions["hbatt"]
            M.substitutions.update({"etasolar": etamax})
            M.cost = M["hbatt"]
            sol = M.localsolve("mosek")
            time += sol["soltime"]
            numsolves += 1
            hmin = sol["cost"].magnitude + 1e-3
            if hmin >= hmax:
                break
            del M.substitutions["etasolar"]
            M.cost = M["etasolar"]
            hbatts = np.linspace(hmin, hmax, 10)
            tol = 0.01
            notpassing = True
            etas = []
            for h in hbatts:
                M.substitutions["hbatt"] = h
                try:
                    sol = M.localsolve("mosek")
                    time += sol["soltime"]
                    numsolves += 1
                    etas.append(sol("etasolar").magnitude)
                except RuntimeWarning:
                    etas.append(np.NaN)
            data["%d-batt" % lat] = hbatts
            data["%d-cell" % lat] = etas

    df = pd.DataFrame(data)
    return df

def plot_battsolar(df, lats):

    fig, ax = plt.subplots()
    lines = []
    midx = []
    maxetas = []
    for la in lats:
        hbatts, etas = df["%d-batt" % la], df["%d-cell" % la]
        maxetas.append(max(etas))
        if la % 4 == 0:
            l = ax.plot(hbatts, etas, "k", label="%d$^{\\circ}$ Lat" % la)
            if la is not 20:
                lines.append(l[0])
            midx.append(np.median(hbatts))
        elif la % 2 == 0:
            l = ax.plot(hbatts, etas, "--", c="0.5", label="%d$^{\\circ}$ Lat" % la)

    labelLines(lines, align=False, xvals=midx, zorder=[10]*len(lines))
    ax.set_ylabel("Solar Cell Efficiency")
    ax.set_xlabel("Battery Specific Energy [Whr/kg]")
    ax.set_xlim([250, 500])
    ax.set_ylim([0.05, max(maxetas)])
    ax.grid()
    return fig, ax

def test():
    _, _ = plot_battsolarcon()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        path = sys.argv[1]
        GENERATE = True if sys.argv[2] == "GEN" else False
    else:
        path = ""
        GENERATE = False

    Lats = range(20, 44, 2)
    if GENERATE:
        DF = battsolarcon(Lats)
        DF.to_csv("battsolarcon.generated.csv")
    else:
        DF = pd.read_csv("battsolarcon.generated.csv")
        del DF["Unnamed: 0"]

    fig, ax = plot_battsolar(DF, Lats)
    fig.savefig(path + "battsolarcon.pdf", bbox_inches="tight")
