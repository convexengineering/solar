from solar import Mission, Aircraft
import matplotlib.pyplot as plt

SP = True
Wtot = []
N = [1, 3, 5, 7, 9, 0]
for i in N:
    Vehicle = Aircraft(Npod=i, sp=SP)
    M = Mission(Vehicle, latitude=[28])
    M.cost = M[M.aircraft.Wtotal]
    sol = M.localsolve("mosek")
    Wtot.append(sol("Wtotal").magnitude)

fig, ax = plt.subplots()
N[-1] = 11
ax.plot(N, Wtot)
ax.set_xlabel("Number of pods")
ax.set_ylabel("Max Take off Weight [lbf]")
ax.grid()
ax.set_ylim([0, max(Wtot)*1.1])
ax.set_xlim([1, 11])
ax.set_xticks(N)
N[-1] = 0
ax.set_xticklabels(["%d" % n for n in N])
fig.savefig("npod_trade.pdf", bbox_inches="tight")


