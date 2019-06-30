" number of pods trade study "
# from solar import Mission, Aircraft
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
from relaxed_constants import relaxed_constants, post_process
from subprocess import CalledProcessError

def pods(N=[1, 3, 5, 7, 9, 0], Nplot=5):
    "trade number of pods"
    SP = True
    data = {}
    for i in N:
        print("\nN=%i" % i)
        from solar import Mission, Aircraft
        Vehicle = Aircraft(Npod=i, sp=SP)
        M = Mission(Vehicle, latitude=[20])
        M.cost = M[M.aircraft.Wtotal]
        try:
            sol = M.localsolve()
            data[i] = sol("Wtotal").magnitude
        except RuntimeWarning as e:
            print(e)
            V2 = Aircraft(Npod=i, sp=SP)
            M2 = Mission(V2, latitude=[20])
            M2.cost = M2[M2.aircraft.Wtotal]
            feas = relaxed_constants(M2)
            sol2 = feas.localsolve()
            vks = post_process(sol2)
            data[i] = np.NaN if vks else sol2("Wtotal").magnitude
            M, sol = M2, sol2
        except CalledProcessError as e:
            print(e)
            data[i] = np.NaN  # mosek_cli can't process certain Ns

        if Nplot == i:
            plot_shear(M, sol)

    df = pd.DataFrame(data, index=[0])
    return df

def plot_shear(model, result):
    " plot shear and moment diagrams "

    S = result(model.mission[1].winggust.S)
    m = result(model.mission[1].winggust.M)
    fig, ax = plt.subplots(2)
    ax[0].plot(range(20), S)
    ax[1].plot(range(20), m)
    ax[0].grid(); ax[1].grid()
    fig.savefig("shearandmoment.pdf")

    S = result(model.mission[1].wingg.S)
    m = result(model.mission[1].wingg.M)
    fig, ax = plt.subplots(2)
    ax[0].plot(range(20), S)
    ax[1].plot(range(20), m)
    ax[0].grid(); ax[1].grid()
    fig.savefig("shearandmoment2.pdf")

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
    pods(Nplot=100)

if __name__ == "__main__":
    test()
    # if len(sys.argv) > 1:
    #     path = sys.argv[1]
    #     GENERATE = True if sys.argv[2] == "GEN" else False
    # else:
    #     path = ""
    #     GENERATE = True
    #
    # if GENERATE:
    #     DF = pods()
    #     DF.to_csv("pds.generated.csv")
    # else:
    #     DF = pd.read_csv("pds.generated.csv")
    #     del DF["Unnamed: 0"]
    #
    # f, a = plot_pods(DF)
    # f.savefig(path + "npod_trade.pdf", bbox_inches="tight")
