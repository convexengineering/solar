" create bar chart "
import sys
import matplotlib.pyplot as plt
import numpy as np
from gpkit.repr_conventions import unitstr
from solar import Mission

#pylint: disable=invalid-name, anomalous-backslash-in-string

def get_highestsens(model, res, varnames=None):
    " plot bar chart of sensitivities "
    pss = []
    ngs = []
    sens = {}
    if varnames:
        for vname in varnames:
            sen = res["sensitivities"]["constants"][vname]
            if hasattr(sen, "__len__"):
                val = max(np.abs(sen.values()))
                vk = [svk for svk in sen if abs(sen[svk]) == val][0]
                sen = sum(sen.values())
            else:
                vk = model[vname].key
            sens[vk] = sen
    else:
        for s in res["sensitivities"]["constants"]:
            sens[model[s].key] = sum(
                np.hstack([res["sensitivities"]["constants"][s]]))

    labels = []
    i = 0
    sorted_sens = sorted(sens.items(), key=lambda x: np.absolute(x[1]),
                         reverse=True)

    for s in sorted_sens:
        if i > 10:
            break
        i += 1
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
    sensdict = {"positives": pss, "negatives": ngs, "indicies": ind,
                "labels": labels}
    return sensdict

def plot_chart(sensdict):
    "plot sensitivities on bar chart"
    fig, ax = plt.subplots()
    ax.bar(sensdict["indicies"], sensdict["positives"], 0.5, color="#4D606E")
    ax.bar(sensdict["indicies"], sensdict["negatives"], 0.5, color="#3FBAC2")
    ax.set_xlim([0.0, sensdict["indicies"][-1]+1])
    ax.set_xticks(sensdict["indicies"])
    ax.set_xticklabels(sensdict["labels"], rotation=-45, ha="left")
    ax.legend(["Positive", "Negative"])
    ax.set_ylabel("sensitivities")
    return fig, ax

def test():
    " test for integrated testing "
    model = Mission(latitude=[10])
    model.cost = model["W_{total}"]
    result = model.solve("mosek")
    _ = get_highestsens(model, result)

    vn = {"W_{pay}": "$W_{\\mathrm{pay}}$",
          "\\eta_{charge}": "$\\eta_{\\mathrm{charge}}$"}
    _ = get_highestsens(model, result, vn)

if __name__ == "__main__":

    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        path = ""

    M = Mission(latitude=[25])
    M.cost = M["W_{total}"]
    sol = M.solve("mosek")

    vns = {"W_{pay}": "$W_{\\mathrm{pay}}$",
           "\\eta_{charge}": "$\\eta_{\\mathrm{charge}}$",
           "\\eta_{discharge}": "$\\eta_{\\mathrm{discharge}}$",
           "h_{batt}": "$h_{\\mathrm{batt}}$",
           "\\eta_Mission/Aircraft/SolarCells": "$\\eta_{\\mathrm{solar}}$",
           "N_{max}": "$N_{\\mathrm{max}}$",
           "e": "$e$", "\\eta_{prop}": "$\\eta_{\\mathrm{prop}}$"}

    sd = get_highestsens(M, sol, vns)
    f, a = plot_chart(sd)
    f.savefig(path + "sensbar.pdf", bbox_inches="tight")
