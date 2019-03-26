"sweep latitude for gas and solar"
from solar import Mission, Aircraft
from gassolar.environment.wind_speeds import get_windspeed
from relaxed_constants import relaxed_constants, post_process
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys

FIGPATH = (os.path.abspath(__file__).replace(os.path.basename(__file__), "")
           + "figs" + os.sep)

plt.rcParams.update({'font.size':15})

def lats():
    fig, ax = plt.subplots()
    lat = np.arange(10, 35, 1)
    data = {"l": lat}
    faillat = []
    stime = 0.0
    snsols = 0.0
    for a in [80, 90, 95]:
        ws = []
        runagains = True
        for l in lat:
            if runagains:
                V = Aircraft(sp=True)
                M = Mission(V, latitude=[l])
                M.substitutions[M.mission[1].fs.pct] = a/100.
                M.cost = M[M.aircraft.Wtotal]
                try:
                    sol = M.localsolve("mosek")
                    stime += sol["soltime"]
                    snsols += 1
                    ws.append(sol("Wtotal").magnitude)
                except RuntimeWarning:
                    ws.append(np.nan)
                    faillat.append(l)
                    runagains = False
                    # feas = relaxed_constants(M)
                    # sol = feas.localsolve("mosek")
                    # vks = post_process(sol)
                    # if vks:
                    #     ws.append(np.nan)
                    #     faillat.append(l)
                    #     runagains = False
                    # else:
                    #     ws.append(sol("Wtotal").magnitude)
            else:
                ws.append(np.nan)
        data[a] = np.hstack(ws)

    print "Solar: %d solves in %.4f seconds" % (snsols, stime)
    df = pd.DataFrame(data)
    return df

def plot_lats(df):

    lat = df["l"]
    psolar = [list(df["80"].values), list(df["90"].values), list(df["95"].values)]

    fig, ax = plt.subplots()
    indl = psolar[0].index(max(psolar[0]))
    indh = psolar[2].index(max(psolar[2]))
    a = (psolar[2][indh]-psolar[0][indl])/(lat[indh]-lat[indl])
    b = psolar[0][indl]-a*lat[indl]
    c = a*lat[indh+1:indl+1] + b
    ax.fill_between(lat[0:indl+1], psolar[0][0:indl+1],
                    np.append(np.array(psolar[-1][0:indh+1]), c), alpha=0.3,
                    facecolor="b", edgecolor="None")
    psolar[1][np.where(lat == 31)[0][0]] = np.nan
    ax.plot(lat, psolar[1], "b", lw=2)
    ax.plot(lat, psolar[0], "b")
    ax.plot(lat, psolar[2], "b")

    ax.annotate("80%", xy=(30,psolar[0][np.where(lat==30)[0][0]]), xytext=(15,15), textcoords="offset points", arrowprops=dict(arrowstyle="-"), fontsize=12)
    ax.annotate("90%", xy=(28, psolar[1][np.where(lat==28)[0][0]]), xytext=(15,15), textcoords="offset points", arrowprops=dict(arrowstyle="-"), fontsize=12)
    ax.annotate("95%", xy=(27,psolar[2][np.where(lat==27)[0][0]]), xytext=(-30,15), textcoords="offset points", arrowprops=dict(arrowstyle="-"), fontsize=12)

    ax.set_ylim([0, max(np.hstack(psolar))*1.1])
    ax.set_xlim([10, 40])
    ax.grid()
    ax.set_xlabel("Latitude Requirement [deg]")
    ax.set_ylabel("Max Take Off Weight [lbf]")
    labels = ["$\\pm$" + item.get_text() for item in ax.get_xticklabels()]
    labels = ["$\\pm$%d" % l for l in np.linspace(20, 40, len(labels))]
    ax.set_xticklabels(labels)
    # ax.legend(["Solar-electric Powered", "Gas Powered"], fontsize=15, loc=2)
    return fig, ax

if __name__ == "__main__":
    if len(sys.argv) > 1:
        path = sys.argv[1]
        GENERATE = True if sys.argv[2] == "GEN" else False
    else:
        path = ""
        GENERATE = False

    if GENERATE:
        DF = lats()
        DF.to_csv("lats.generated.csv")
    else:
        DF = pd.read_csv("lats.generated.csv")
        del DF["Unnamed: 0"]

    fig, ax = plot_lats(DF)
    fig.savefig(path + "mtowvslat.pdf", bbox_inches="tight")
