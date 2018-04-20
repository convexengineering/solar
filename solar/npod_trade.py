" number of pods trade study "
from solar import Mission, Aircraft
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

def pods():
    "trade number of pods"
    SP = True
    data = {}
    N = [1, 3, 5, 7, 9, 0]
    for i in N:
        Vehicle = Aircraft(Npod=i, sp=SP)
        M = Mission(Vehicle, latitude=[15])
        M.cost = M[M.aircraft.Wtotal]
        sol = M.localsolve("mosek")
        data[i] = sol("Wtotal").magnitude

    df = pd.DataFrame(data, index=[0])
    return df

def plot_pods(df):
    " plot pod trade "

    N = list(np.array(df.columns.values, float))
    del N[0]
    N.append(11)
    Wtot = list(df.loc[0].values)
    Wtot.append(Wtot[0])
    del Wtot[0]
    fig, ax = plt.subplots()
    ax.plot(N, Wtot)
    ax.set_xlabel("Number of battery pods")
    ax.set_ylabel("Max Take off Weight [lbf]")
    ax.grid()
    ax.set_ylim([0, max(Wtot)*1.1])
    ax.set_xlim([1, 11])
    ax.set_xticks(N)
    del N[-1]
    labels = ["%d" % n for n in N] + ["battery\nin wing"]
    ax.set_xticklabels(labels)
    fig.savefig("npod_trade.pdf", bbox_inches="tight")
    return fig, ax

def test():
    " for unit testing "
    pods()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        path = ""

    GENERATE = False

    if GENERATE:
        DF = pods()
        DF.to_csv("pds.generated.csv")
    else:
        DF = pd.read_csv("pds.generated.csv")
        del DF["Unnamed: 0"]

    f, a = plot_pods(DF)
    f.savefig(path + "npod_trade.pdf", bbox_inches="tight")
