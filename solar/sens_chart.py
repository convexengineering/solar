from solar import Mission
import matplotlib.pyplot as plt
import numpy as np
import sys
from gpkit.repr_conventions import unitstr

#pylint: disable=invalid-name, anomalous-backslash-in-string

def plot_sens(model, sol, varnames=None):
    fig, ax = plt.subplots()
    pss = []
    ngs = []
    sens = {}
    if varnames:
        for vname in varnames:
            sen = sol["sensitivities"]["constants"][vname]
            if hasattr(sen, "__len__"):
                val = max(np.abs(sen.values()))
                vk = [svk for svk in sen if abs(sen[svk])==val][0]
                # sen = sol["sensitivities"]["constants"][vk]
                sen = sum(sen.values())
            else:
                vk = model[vname].key
            sens[vk] = sen
    else:
        for s in sol["sensitivities"]["constants"]:
            sens[s] = sum(np.hstack([sol["sensitivities"]["constants"][s]]))

    labels = []
    i = 0
    sorted_sens = sorted(sens.items(), key=lambda x: np.absolute(x[1]),
                         reverse=True)

    for s in sorted_sens:
        i += 1
        if i > 10:
            break
        if hasattr(model[s[0]], "__len__"):
            vk = model[s[0]][0]
        else:
            vk = s[0]
        val = sum(np.hstack([model.substitutions[vk]]))
        if "units" in vk.descr:
            uts = unitstr(vk.descr["units"])
        else:
            uts = "-"

        if varnames:
            for vn in varnames:
                if vk.descr["name"] == vn:
                    lbl = varnames[vk.descr["name"]]
                else:
                    allname = (vk.descr["name"] + "_" +
                               "/".join([x for x in vk.descr["models"]]))
                    if allname == vn:
                        lbl = varnames[allname]
        else:
            lbl = vk.descr["label"]
        labels.append(lbl + "$ =%.2f$ [%s]" % (val, uts))
        if s[1] > 0:
            pss.append(s[1])
            ngs.append(0)
        else:
            ngs.append(abs(s[1]))
            pss.append(0)

    ind = np.arange(0.5, i + 0.5, 1)
    ax.bar(ind, pss, 0.5, color="#4D606E")
    ax.bar(ind, ngs, 0.5, color="#3FBAC2")
    ax.set_xlim([0.0, ind[-1]+1])
    ax.set_xticks(ind)
    ax.set_xticklabels(labels, rotation=-45, ha="left")
    ax.legend(["Positive", "Negative"])
    ax.set_ylabel("sensitivities")
    return fig, ax

if __name__ == "__main__":
    M = Mission(latitude=25)
    M.cost = M["W_{total}"]
    sol = M.solve("mosek")

    vns = {"W_{pay}": "$W_{\\mathrm{pay}}$",
           "\\eta_{charge}": "$\\eta_{\\mathrm{charge}}$",
           "\\eta_{discharge}": "$\\eta_{\\mathrm{discharge}}$",
           "h_{batt}": "$h_{\\mathrm{batt}}$",
           "\\eta_Mission/Aircraft/SolarCells": "$\\eta_{\\mathrm{solar}}$",
           "N_{max}": "$N_{\\mathrm{max}}$",
           "e": "$e$", "\\eta_{prop}": "$\\eta_{\\mathrm{prop}}$"}

    fig, ax = plot_sens(M, sol, vns)
    fig.savefig("sensbar.pdf", bbox_inches="tight")
