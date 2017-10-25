" Simple Solar-Electric Powered Aircraft Model "
#pylint: disable=attribute-defined-outside-init, invalid-name, unused-variable
#pylint: disable=too-many-locals, redefined-variable-type
#pylint: disable=too-many-instance-attributes
import os
import pandas as pd
import numpy as np
import gassolar.environment
from gassolar.environment.solar_irradiance import get_Eirr, twi_fits
from gassolar.environment.wind_speeds import get_month
from gpkit import Model, Variable
from gpkit.tests.helpers import StdoutCaptured
from gpkitmodels.GP.aircraft.wing.wing import Wing as WingGP
from gpkitmodels.SP.aircraft.wing.wing import Wing as WingSP
from gpkitmodels.GP.aircraft.wing.boxspar import BoxSpar
from gpkitmodels.GP.aircraft.wing.capspar import CapSpar
from gpkitmodels.GP.aircraft.tail.empennage import Empennage
from gpkitmodels.GP.aircraft.tail.tail_boom import TailBoomState
from gpkitmodels.SP.aircraft.tail.tail_boom_flex import TailBoomFlexibility
from gpkitmodels.GP.aircraft.fuselage.elliptical_fuselage import Fuselage
from gpkitmodels.tools.summing_constraintset import summing_vars
from gpfit.fit_constraintset import FitCS as FCS

path = os.path.dirname(gassolar.environment.__file__)

class Aircraft(Model):
    "vehicle"
    fuseModel = None
    def setup(self, sp=False):

        self.sp = sp
        self.empennage = Empennage()
        self.solarcells = SolarCells()
        if sp:
            WingSP.sparModel = CapSpar
            WingSP.fillModel = None
            self.wing = WingSP()
        else:
            WingGP.sparModel = CapSpar
            WingGP.fillModel = None
            self.wing = WingGP()
        self.battery = Battery()
        self.motor = Motor()

        self.components = [self.solarcells, self.wing, self.battery,
                           self.empennage, self.motor]
        loading = []
        if sp:
            tbstate = TailBoomState()
            loading = TailBoomFlexibility(self.empennage.htail,
                                          self.empennage.tailboom,
                                          self.wing, tbstate)

        Wpay = Variable("W_{pay}", 10, "lbf", "payload weight")
        Wavn = Variable("W_{avn}", 8, "lbf", "avionics weight")
        Wtotal = Variable("W_{total}", "lbf", "aircraft weight")
        Wwing = Variable("W_{wing}", "lbf", "wing weight")
        Wcent = Variable("W_{cent}", "lbf", "center weight")

        self.empennage.substitutions["V_v"] = 0.04
        self.wing.substitutions[self.wing.topvar("m_{fac}")] = 1.1

        if not sp:
            self.empennage.substitutions["V_h"] = 0.45
            self.empennage.substitutions["m_h"] = 0.1

        constraints = [
            self.solarcells["S"]/self.solarcells["m_{fac}"] <= self.wing["S"],
            self.empennage.htail["V_h"] <= (
                self.empennage.htail["S"]
                * self.empennage.htail["l_h"]/self.wing["S"]**2
                * self.wing["b"]),
            self.empennage.vtail["V_v"] <= (
                self.empennage.vtail["S"]
                * self.empennage.vtail["l_v"]/self.wing["S"]
                / self.wing["b"]),
            self.empennage.tailboom["d_0"] <= (
                self.wing["\\tau"]*self.wing["c_{root}"])
            ]

        if self.fuseModel:
            self.fuselage = self.fuseModel()
            self.components.extend([self.fuselage])
            constraints.extend([
                self.battery["\\mathcal{V}"] <= self.fuselage["\\mathcal{V}"],
                Wwing >= self.wing.topvar("W") + self.solarcells["W"],
                Wcent >= (Wpay + Wavn + self.empennage.topvar("W")
                          + self.motor["W"] + self.fuselage["W"]
                          + self.battery["W"]),
                ])
        else:
            constraints.extend([
                Wwing >= (sum(summing_vars([self.wing, self.battery,
                                            self.solarcells], "W"))),
                Wcent >= (Wpay + Wavn +
                          sum(summing_vars([self.empennage, self.motor], "W"))),
                self.battery["\\mathcal{V}"] <= (
                    self.wing["c_{MAC}"]**2*0.5*self.wing["\\tau"]
                    * self.wing["b"])
                ])

        constraints.extend([Wtotal >= (
            Wpay + Wavn + sum(summing_vars(self.components, "W")))])

        return constraints, self.components, loading

    def flight_model(self, state):
        " what happens during flight "
        return AircraftPerf(self, state)

class Motor(Model):
    "the thing that provides power"
    def setup(self):

        W = Variable("W", "lbf", "motor weight")
        Pmax = Variable("P_{max}", "W", "max power")
        Bpm = Variable("B_{PM}", 4140.8, "W/kg", "power mass ratio")
        m = Variable("m", "kg", "motor mass")
        g = Variable("g", 9.81, "m/s**2", "gravitational constant")
        eta = Variable("\\eta", 0.95, "-", "motor efficiency")

        constraints = [Pmax == Bpm*m,
                       W >= m*g]

        return constraints

class Battery(Model):
    "battery model"
    def setup(self):

        W = Variable("W", "lbf", "battery weight")
        eta_charge = Variable("\\eta_{charge}", 0.98, "-",
                              "charging efficiency")
        eta_discharge = Variable("\\eta_{discharge}", 0.98, "-",
                                 "discharging efficiency")
        E = Variable("E", "J", "total battery energy")
        g = Variable("g", 9.81, "m/s**2", "gravitational constant")
        hbatt = Variable("h_{batt}", 350, "W*hr/kg", "battery specific energy")
        vbatt = Variable("(E/\\mathcal{V})", 800, "W*hr/l",
                         "volume battery energy density")
        Volbatt = Variable("\\mathcal{V}", "m**3", "battery volume")

        constraints = [W >= E/hbatt*g,
                       Volbatt >= E/vbatt,
                       eta_charge == eta_charge,
                       eta_discharge == eta_discharge]

        return constraints

class SolarCells(Model):
    "solar cell model"
    def setup(self):

        rhosolar = Variable("\\rho_{solar}", 0.27, "kg/m^2",
                            "solar cell area density")
        g = Variable("g", 9.81, "m/s**2", "gravitational constant")
        S = Variable("S", "ft**2", "solar cell area")
        W = Variable("W", "lbf", "solar cell weight")
        etasolar = Variable("\\eta", 0.22, "-", "solar cell efficiency")
        mfac = Variable("m_{fac}", 1.0, "-", "solar cell area margin")

        constraints = [W >= rhosolar*S*g]

        return constraints

class AircraftPerf(Model):
    "aircraft performance"
    def setup(self, static, state):

        self.wing = static.wing.flight_model(static.wing, state)
        self.htail = static.empennage.htail.flight_model(static.empennage.htail,
                                                         state)
        self.vtail = static.empennage.vtail.flight_model(static.empennage.vtail,
                                                         state)
        self.tailboom = static.empennage.tailboom.flight_model(state)

        self.flight_models = [self.wing, self.htail, self.vtail,
                              self.tailboom]
        areadragmodel = [self.htail, self.vtail, self.tailboom]
        areadragcomps = [static.empennage.htail,
                         static.empennage.vtail,
                         static.empennage.tailboom]

        self.wing.substitutions["e"] = 0.95

        if static.fuseModel:
            self.fuse = static.fuselage.flight_model(state)
            areadragmodel.extend([self.fuse])
            areadragcomps.extend([static.fuselage])
            self.flight_models.extend([self.fuse])

        CD = Variable("C_D", "-", "aircraft drag coefficient")
        cda = Variable("CDA", "-", "non-wing drag coefficient")
        Pshaft = Variable("P_{shaft}", "hp", "shaft power")
        Pavn = Variable("P_{avn}", 0.0, "W", "Accessory power draw")
        Poper = Variable("P_{oper}", "W", "operating power")
        mfac = Variable("m_{fac}", 1.05, "-", "drag margin factor")

        dvars = []
        for dc, dm in zip(areadragcomps, areadragmodel):
            if "C_d" in dm.varkeys:
                dvars.append(dm["C_d"]*dc["S"]/static.wing["S"])
            elif "C_f" in dm.varkeys:
                dvars.append(dm["C_f"]*dc["S"]/static.wing["S"])

        constraints = [
            state["(E/S)_{irr}"] >= (
                state["(E/S)_{day}"] + static.battery["E"]
                / static.battery["\\eta_{charge}"]
                / static.solarcells["\\eta"]/static.solarcells["S"]),
            static.battery["E"]*static.battery["\\eta_{discharge}"] >= (
                Poper*state["t_{night}"]
                + state["(E/S)_C"]*static.solarcells["\\eta"]
                * static.solarcells["S"]),
            Poper >= Pavn + Pshaft/static.motor["\\eta"],
            Poper == (state["(P/S)_{min}"]*static.solarcells["S"]
                      * static.solarcells["\\eta"]),
            cda >= sum(dvars),
            CD/mfac >= cda + self.wing["C_d"],
            Poper <= static.motor["P_{max}"]
            ]

        return self.flight_models, constraints

class FlightState(Model):
    """
    environmental state of aircraft

    inputs
    ------
    latitude: earth latitude [deg]
    altitude: flight altitude [ft]
    percent: percentile wind speeds [%]
    day: day of the year [Jan 1st = 1]
    """
    def setup(self, latitude=45, day=355):

        month = get_month(day)
        df = pd.read_csv(path + os.sep + "windfits" + month +
                         "/windaltfit_lat%d.csv" % latitude).to_dict(
                             orient="records")[0]
        with StdoutCaptured(None):
            dft, dfd = twi_fits(latitude, day, gen=True)
        esirr, _, tn, _ = get_Eirr(latitude, day)

        Vwind = Variable("V_{wind}", "m/s", "wind velocity")
        V = Variable("V", "m/s", "true airspeed")
        rho = Variable("\\rho", "kg/m**3", "air density")
        mu = Variable("\\mu", 1.42e-5, "N*s/m**2", "viscosity")
        ESirr = Variable("(E/S)_{irr}", esirr, "W*hr/m^2",
                         "solar energy")
        PSmin = Variable("(P/S)_{min}", "W/m^2",
                         "minimum necessary solar power")
        ESday = Variable("(E/S)_{day}", "W*hr/m^2",
                         "solar cells energy during daytime")
        ESc = Variable("(E/S)_C", "W*hr/m^2",
                       "energy for batteries during sunrise/set")
        ESvar = Variable("(E/S)_{ref}", 1, "W*hr/m^2", "energy units variable")
        PSvar = Variable("(P/S)_{ref}", 1, "W/m^2", "power units variable")
        tnight = Variable("t_{night}", tn, "hr", "night duration")
        pct = Variable("p_{wind}", 0.9, "-", "percentile wind speeds")
        Vwindref = Variable("V_{wind-ref}", 100.0, "m/s",
                            "reference wind speed")
        rhoref = Variable("\\rho_{ref}", 1.0, "kg/m**3",
                          "reference air density")
        mfac = Variable("m_{fac}", 1.0, "-", "wind speed margin factor")

        constraints = [
            V/mfac >= Vwind,
            FCS(df, Vwind/Vwindref, [rho/rhoref, pct]),
            FCS(dfd, ESday/ESvar, [PSmin/PSvar]),
            FCS(dft, ESc/ESvar, [PSmin/PSvar]),
            ]

        return constraints

def altitude(density):
    " find air density "
    g = 9.80665 # m/s^2
    R = 287.04 # m^2/K/s^2
    T11 = 216.65 # K
    p11 = 22532 # Pa
    p = density*R*T11
    h = (11000 - R*T11/g*np.log(p/p11))/0.3048
    return h

class FlightSegment(Model):
    "flight segment"
    def setup(self, aircraft, latitude=35, day=355):

        self.latitude = latitude
        self.day = day

        self.aircraft = aircraft
        self.fs = FlightState(latitude=latitude, day=day)
        self.aircraftPerf = self.aircraft.flight_model(self.fs)
        self.slf = SteadyLevelFlight(self.fs, self.aircraft,
                                     self.aircraftPerf)

        self.loading = self.aircraft.wing.loading(
            self.aircraft.wing, self.aircraft["W_{cent}"],
            self.aircraft["W_{wing}"], self.aircraftPerf["V"],
            self.aircraftPerf["C_L"])

        for vk in self.loading.varkeys["N_{max}"]:
            if "ChordSparL" in vk.descr["models"]:
                self.loading.substitutions.update({vk: 5})
            if "GustL" in vk.descr["models"]:
                self.loading.substitutions.update({vk: 2})

        self.loading.substitutions["V_{gust}"] = 5

        self.submodels = [self.fs, self.aircraftPerf, self.slf, self.loading]

        return self.submodels

class SteadyLevelFlight(Model):
    "steady level flight model"
    def setup(self, state, aircraft, perf):

        T = Variable("T", "N", "thrust")
        etaprop = Variable("\\eta_{prop}", 0.8, "-", "propeller efficiency")

        constraints = [
            aircraft["W_{total}"] <= (
                0.5*state["\\rho"]*state["V"]**2*perf["C_L"]
                * aircraft.wing["S"]),
            T >= (0.5*state["\\rho"]*state["V"]**2*perf["C_D"]
                  *aircraft.wing["S"]),
            perf["P_{shaft}"] >= T*state["V"]/etaprop]

        return constraints

class Mission(Model):
    "define mission for aircraft"
    def setup(self, latitude=range(1, 21, 1), day=355, sp=False):

        self.solar = Aircraft(sp=sp)
        self.mission = []
        if day == 355 or day == 172:
            for l in latitude:
                self.mission.append(FlightSegment(self.solar, l, day))
        else:
            assert day < 172
            for l in latitude:
                self.mission.append(FlightSegment(self.solar, l, day))
                self.mission.append(FlightSegment(self.solar, l,
                                                  355 - 10 - day))

        return self.mission, self.solar

    def process_result(self, result):
        super(Mission, self).process_result(result)
        result["latitude"] = []
        result["day"] = []
        sens = result["sensitivities"]["constants"]
        for f in self.mission:
            const = False
            for sub in f.substitutions:
                if sub not in f.varkeys:
                    continue
                if sub in self.solar.varkeys:
                    continue
                if any(s > 1e-5 for s in np.hstack([abs(sens[sub])])):
                    const = True
                    break
            if const:
                print "%d is a constraining latitude" % f.latitude
                result["latitude"].extend([f.latitude])
                result["day"].extend([f.day])
                continue
            for vk in f.varkeys:
                if vk in self.solar.varkeys:
                    continue
                del result["variables"][vk]
                if vk in result["freevariables"]:
                    del result["freevariables"][vk]
                else:
                    del result["constants"][vk]
                    del result["sensitivities"]["constants"][vk]
        self.setup(latitude=result["latitude"])

def test():
    " test model for continuous integration "
    m = Mission()
    m.cost = m["W_{total}"]
    m.solve()
    m = Mission(sp=True)
    m.cost = m["W_{total}"]
    m.localsolve()

if __name__ == "__main__":
    M = Mission(sp=False)
    M.cost = M["W_{total}"]
    sol = M.solve("mosek")
