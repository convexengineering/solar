from solar import Mission, Aircraft
from numpy import linspace, array
import pandas as pd
import matplotlib.pyplot as plt
import sys

def taus():

    V = Aircraft(sp=True)
    M = Mission(V, latitude=[20])
    M.cost = M.aircraft.Wtotal

    taumax = linspace(0.11, 0.145, 10)
    data = {}
    for t in taumax:
        M.substitutions[M.aircraft.maxtau] = t
        sol = M.localsolve("mosek")
        data[t] = sol("Wtotal").magnitude

    df = pd.DataFrame(data, index=[0])
    return df

def plot_taus(df):
    fig, ax = plt.subplots()

    ts = list(array(df.columns.values, float))
    ws = list(df.loc[0].values)
    ax.plot(ts, ws)
    ax.grid()
    ax.set_ylim([0, max(ws)*1.1])
    ax.set_xlim([min(ts), max(ts)])
    ax.set_xlabel("t/c - thickness to chord ratio")
    ax.set_ylabel("Max Take off Weight [lbf]")
    return fig, ax


if __name__ == "__main__":
    if len(sys.argv) > 1:
        path = sys.argv[1]
        GENERATE = True if sys.argv[2] == "GEN" else False
    else:
        path = ""
        GENERATE = False

    if GENERATE:
        DF = taus()
        DF.to_csv("taus.generated.csv")
    else:
        DF = pd.read_csv("taus.generated.csv")
        del DF["Unnamed: 0"]

    f, a = plot_taus(DF)
    f.savefig(path + "maxtau.pdf", bbox_inches="tight")

