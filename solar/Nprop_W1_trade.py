" number of motors trade study "
# from solar import Mission, Aircraft
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
from relaxed_constants import relaxed_constants, post_process

def props(Nprops=[1, 2, 3, 4, 5, 7, 9], W1 = [300,687, 1200], Nplot=5):
    "trade number of motors"
    SP = True
    data = {}
    N = range(len(Nprops))
    W = range(len(W1))
    wmto = []
    Kv   = []
    wmotor = []
    etam = []
    etap = []
    W1s = []
    Npropss= []
    omega = []
    Rprop = []

    for j in W:
        for i in N:
            from solar import Mission, Aircraft
            Vehicle = Aircraft(Npod=1, sp=SP)
            Vehicle.substitutions[Vehicle.Nprop] = Nprops[i]
            Vehicle.substitutions[Vehicle.motor.W1] = W1[j]
            M = Mission(Vehicle, latitude=[20])
            M.cost = M[M.aircraft.Wtotal]
            try:
                sol = M.localsolve("mosek")
                wmto.append(sol("Wtotal").magnitude)
                etam.append(sol("etam_Mission/FlightSegment/AircraftPerf/AircraftDrag/MotorPerf").magnitude)
                etap.append(sol("eta_Mission/FlightSegment/AircraftPerf/AircraftDrag/BladeElementProp").magnitude)
                omega.append(sol("omega_Mission/FlightSegment/AircraftPerf/AircraftDrag/BladeElementProp").magnitude)
                Rprop.append(sol("R_Aircraft/Propeller").magnitude)
                Kv.append(sol("Kv").magnitude)
                wmotor.append(float(sol("W_Aircraft/Motor").magnitude*Nprops[i]))
                W1s.append(W1[j])
                Npropss.append(Nprops[i])

            except RuntimeWarning:
                pass
                #V2 = Aircraft(Npod=1, sp=SP)
                #Vehicle.substitutions[Vehicle.Nprop] = i
                #M2 = Mission(V2, latitude=[20])
                #M2.cost = M2[M2.aircraft.Wtotal]
                #feas = relaxed_constants(M2)
                #sol2 = feas.localsolve("mosek")
                #vks = post_process(sol2)
                #data[i] = np.NaN if vks else sol2("Wtotal").magnitude
                #M, sol = M2, sol2
    data = {"Wtotal":wmto, "etam":etam, "etap":etap, "wmotor":wmotor, "Nprops":Npropss, "Kv":Kv, "W1s":W1s, "omega":omega, "Rprop":Rprop}
    df = pd.DataFrame(data)
    return df

def plot_props(df, C = 7, K = 3):
    " plot Nprops trade "

    Wtot = df['Wtotal']
    etam = df['etam']
    etap = df['etap']
    wmotor = df['wmotor']
    Nprops = df['Nprops']
    Kv = df['Kv']
    W1s = df['W1s']
    wmtot = []
    fmotor = []

    for i in range(len(Nprops)):
        fmotor.append(wmotor[i]/Wtot[i]*100.)

    colors = ["#084081", "#0868ac", "#2b8cbe", "#4eb3d3"]
    ls = ["-", "--", "-."]

    fig, ax = plt.subplots(3,sharex = True)
    for j in range(K):
        ax[0].plot(Nprops[0:C],Wtot[0+j*C:C+j*C], label=('%i $\mathrm{kg W}/\mathrm{rpm}^2$'%(float(W1s[j*C])/2.20462)), c= colors[0], ls=ls[j])
        ax[0].set_ylabel("Max Takeoff\nWeight [lbf]")
        ax[0].set_ylim([300, max(Wtot)*1.2])
        ax[0].set_xlim([1, 9])
        ax[1].plot(Nprops[0:C],etam[0+j*C:C+j*C], c = colors[1], ls=ls[j])
        ax[1].plot(Nprops[0:C],etap[0+j*C:C+j*C], c = 'r', ls=ls[j])
        ax[1].set_ylabel("Efficiency")
        ax[2].plot(Nprops[0:C],fmotor[0+j*C:C+j*C], c = colors[1], ls=ls[j])
        ax[2].set_ylabel("Motor weight\nfraction (%MTOW)")
        ax[2].set_xlabel("Number of propellers")
        ax[2].set_ylim([2.5,4.5])


    ax[0].legend(loc=2,ncol=3, title="Motor Weight Parameter W1", fontsize = 10)
    ax[1].legend(["Motor","Propeller"], fontsize = 10)
    fig.savefig("nprop_w1_trade.pdf", bbox_inches="tight")
    return fig, ax

def test():
    " for unit testing "
    props(Nplot=100)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        path = sys.argv[1]
        GENERATE = True if sys.argv[2] == "GEN" else False
    else:
        path = ""
        GENERATE = False

    if GENERATE:
        N = [1,2, 3, 4, 5, 7, 9]
        DF = props(N)
        DF.to_csv("pds.generated.csv")
    else:
        DF = pd.read_csv("pds.generated.csv")
        del DF["Unnamed: 0"]

    f, a = plot_props(DF)
    f.savefig(path + "nprops_w1_trade.pdf", bbox_inches="tight")