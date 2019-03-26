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

def payloadcon(latitudes):
    time = 0
    numsolves = 0
    etamax = []
    lines = []
    midx = []
    data = {}
    passing = True
    wmin = 0.01
    pmin = 0.01
    for lat in latitudes:
        if passing:
            V = Aircraft(sp=True)
            M = Mission(V, latitude=[lat])
            del M.substitutions["Wpay"]
            M.substitutions.update({"Ppay": pmin})
            M.cost = 1./M["Wpay"]
            try:
                sol = M.localsolve("mosek")
            except RuntimeWarning:
                break
            time += sol["soltime"]
            numsolves += 1
            wmax = sol("Wpay").magnitude - 1e-3
            if wmin >= wmax:
                break
            del M.substitutions["Ppay"]
            M.cost = 1./M["Ppay"]
            wpays = np.linspace(wmin, wmax, 10)
            tol = 0.01
            notpassing = True
            ppays = []
            for w in wpays:
                M.substitutions["Wpay"] = w
                try:
                    sol = M.localsolve("mosek")
                    time += sol["soltime"]
                    numsolves += 1
                    ppays.append(sol("Ppay").magnitude)
                except RuntimeWarning:
                    ppays.append(np.NaN)
            data["%d-W" % lat] = wpays
            data["%d-P" % lat] = ppays

    df = pd.DataFrame(data)
    return df

def plot_payloadcon(df, lats):

    lats = np.array(lats)
    maxlat = int(df.columns[-1][0:2])
    lats = lats[lats <= 30]
    fig, ax = plt.subplots()
    lines = []
    midx = []
    maxps = []
    for la in lats:
        wpays, ppays= df["%d-W" % la], df["%d-P" % la]
        maxps.append(max(ppays))
        if la % 4 == 0:
            l = ax.plot(wpays, ppays, "k", label="%d$^{\\circ}$ Lat" % la)
            if la is not 20:
                lines.append(l[0])
            midx.append(np.median(wpays))
        elif la % 2 == 0:
            l = ax.plot(wpays, ppays, "--", c="0.5", label="%d$^{\\circ}$ Lat" % la)

    labelLines(lines, align=False, xvals=midx, zorder=[10]*len(lines))
    ax.set_ylabel("Payload Power [W]")
    ax.set_xlabel("Payload Weight [lbf]")
    ax.set_xlim([0, 500])
    ax.set_ylim([0, max(maxps)])
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

    Lats = range(16, 30, 2)
    if GENERATE:
        DF = payloadcon(Lats)
        DF.to_csv("payloadcon.generated.csv")
    else:
        DF = pd.read_csv("payloadcon.generated.csv")
        del DF["Unnamed: 0"]

    fig, ax = plot_payloadcon(DF, Lats)
    fig.savefig(path + "paycon.pdf", bbox_inches="tight")
