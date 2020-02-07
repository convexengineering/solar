" run seasonal trade "
from builtins import zip
from builtins import range
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from solar import Mission, Aircraft

#pylint: disable=invalid-name

def season(lats):
    " trade seasonal availibility with weight "
    days = [80, 52, 21, 355]

    data = {}
    for l in lats:
        failed = False
        mtows = []
        for d in days:
            if failed:
                mtows = mtows + [np.nan]*(4-len(mtows))
                break
            V = Aircraft(sp=False)
            M = Mission(V, latitude=list(range(1, l+1, 1)), day=d)
            M.cost = M[M.aircraft.Wtotal]
            try:
                sol = M.solve()
                mtow = sol(M.aircraft.Wtotal).magnitude
                mtows.append(mtow)
            except RuntimeWarning:
                mtows.append(np.nan)
                failed = True

        mtows = np.hstack(mtows)
        data[l] = mtows

    df = pd.DataFrame(data)
    return df

def plot_season(df):
    " plot season "
    fig, ax = plt.subplots()
    colors = ["#014636", "#016c59", "#02818a", "#3690c0", "#67a9cf"]
    mrks = ["o", "v", "D", "^", "s"]
    for d, cl, mk in zip(df, colors, mrks):
        ax.plot(list(range(1, 5)), df[d], ms=7, lw=2, color=cl, ls="dashed",
                marker=mk, label="%s$^{\\circ}$ Lat" % d)

    ax.set_xlim([0.5, 4.5])
    ax.set_ylim([0, 300])
    ax.set_ylabel("Max Take-off Weight [lbf]")
    ax.set_xlabel("Seasonal Requirement")
    ax.set_xticklabels(["", "6-months", "", "8-months", "", "10-months", "",
                        "12-months"], rotation=-45, ha="left")
    ax.grid()
    hnd, lbls = ax.get_legend_handles_labels()
    ax.legend(np.flip(hnd, 0), np.flip(lbls, 0), loc=2, fontsize=15,
              numpoints=1)
    return fig, ax

def test():
    "setup for integrated testing"
    _ = season([20])

if __name__ == "__main__":

    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        path = ""

    GENERATE = False

    if GENERATE:
        DF = season(list(range(20, 30, 2)))
        DF.to_csv("season.generated.csv")
    else:
        DF = pd.read_csv("season.generated.csv")
        del DF["Unnamed: 0"]

    f, a = plot_season(DF)
    f.savefig(path + "season.pdf", bbox_inches="tight")
