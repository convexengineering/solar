" Simple Solar-Electric Powered Aircraft Model "
#pylint: disable=attribute-defined-outside-init, invalid-name, unused-variable
#pylint: disable=too-many-locals, redefined-variable-type
#pylint: disable=too-many-instance-attributes
from os.path import abspath, dirname
from os import sep
import os
import pandas as pd
import numpy as np
import gassolar.environment
from gassolar.environment.solar_irradiance import get_Eirr, twi_fits
from gassolar.environment.wind_speeds import get_month
from gpkit import Model, Variable, parse_variables
from gpkit.tests.helpers import StdoutCaptured
from gpkitmodels.GP.aircraft.wing.wing import Wing as WingGP
from gpkitmodels.SP.aircraft.wing.wing import Wing as WingSP
from gpkitmodels.GP.aircraft.wing.boxspar import BoxSpar
from gpkitmodels.GP.aircraft.wing.capspar import CapSpar
from gpkitmodels.GP.aircraft.wing.wing_skin import WingSkin
from gpkitmodels.GP.aircraft.tail.empennage import Empennage
from gpkitmodels.GP.aircraft.tail.horizontal_tail import HorizontalTail
from gpkitmodels.GP.aircraft.tail.vertical_tail import VerticalTail
from gpkitmodels.GP.aircraft.tail.tail_boom import TailBoomState
from gpkitmodels.SP.aircraft.tail.tail_boom_flex import TailBoomFlexibility
from gpkitmodels.GP.aircraft.fuselage.elliptical_fuselage import Fuselage
from gpkitmodels import g
from gpfit.fit_constraintset import FitCS as FCS

path = dirname(gassolar.environment.__file__)

class Aircraft(Model):
    """ Aircraft Model

    Variables
    ---------
    Wpay        10      [lbf]   payload weight
    Wavn        8       [lbf]   avionics weight
    Wtotal              [lbf]   aircraft weight
    Wwing               [lbf]   wing weight
    Wcent               [lbf]   center weight

    Upper Unbounded
    ---------------
    Wtotal, Wwing, Wcent

    LaTex Strings
    -------------
    Wpay        W_{\\mathrm{pay}}
    Wavn        W_{\\mathrm{avn}}
    Wtotal      W_{\\mathrm{total}}
    Wwing       W_{\\mathrm{wing}}
    Wcent       W_{\\mathrm{cent}}

    """
    fuseModel = None
    def setup(self, sp=False):
        exec parse_variables(Aircraft.__doc__)

        self.sp = sp

        HorizontalTail.sparModel = BoxSpar
        HorizontalTail.fillModel = None
        VerticalTail.sparModel = BoxSpar
        VerticalTail.fillModel = None
        self.emp = Empennage(N=5, tailboomSpar=BoxSpar)
        self.solarcells = SolarCells()
        if sp:
            WingSP.sparModel = BoxSpar
            WingSP.fillModel = None
            self.wing = WingSP()
        else:
            WingGP.sparModel = BoxSpar
            WingGP.fillModel = None
            self.wing = WingGP()
        self.battery = Battery()
        self.motor = Motor()

        self.components = [self.solarcells, self.wing, self.battery,
                           self.emp, self.motor]
        loading = []
        if sp:
            tbstate = TailBoomState()
            loading = TailBoomFlexibility(self.emp.htail, self.emp.hbend,
                                          self.wing)

        Sw = self.Sw = self.wing.planform.S
        cmac = self.cmac = self.wing.planform.cmac
        tau = self.tau = self.wing.planform.tau
        croot = self.croot = self.wing.planform.croot
        b = self.b = self.wing.planform.b
        Vh = self.Vh = self.emp.htail.Vh
        lh = self.lh = self.emp.htail.lh
        Sh = self.Sh = self.emp.htail.planform.S
        Vv = self.Vv = self.emp.vtail.Vv
        Sv = self.Sv = self.emp.vtail.planform.S
        lv = self.lv = self.emp.vtail.lv
        d0 = self.d0 = self.emp.tailboom.d0
        Volbatt = self.Volbatt = self.battery.Volbatt
        Ssolar = self.Ssolar = self.solarcells.S
        mfsolar = self.mfsolar = self.solarcells.mfac

        self.emp.substitutions[Vv] = 0.04
        self.wing.substitutions[self.wing.mfac] = 1.1
        if not sp:
            self.emp.substitutions[Vh] = 0.45
            self.emp.substitutions[self.emp.htail.mh] = 0.1


        constraints = [Ssolar*mfsolar <= Sw,
                       Vh <= Sh*lh/Sw/cmac,
                       Vv <= Sv*lv/Sw/b,
                       d0 <= tau*croot,
                      ]

        if self.fuseModel:
            self.fuselage = self.fuseModel()
            self.components.extend([self.fuselage])
            constraints.extend([
                Volbatt <= self.fuselage.Vol,
                Wwing >= self.wing.W + self.solarcells.W,
                Wcent >= (Wpay + Wavn + self.emp.W + self.motor.W
                          + self.fuselage.W + self.battery.W),
                ])
        else:
            constraints.extend([
                Wwing >= sum([c.W for c in [self.wing, self.battery,
                                            self.solarcells]]),
                Wcent >= Wpay + Wavn + self.emp.W + self.motor.W,
                Volbatt <= cmac**2*0.5*tau*b
                ])

        constraints.extend([Wtotal >= (
            Wpay + Wavn + sum([c.W for c in self.components]))])

        return constraints, self.components, loading

    def flight_model(self, state):
        " what happens during flight "
        return AircraftPerf(self, state)

class Motor(Model):
    """ Motor Model

    Variables
    ---------
    W                   [lbf]       motor weight
    Pmax                [W]         max power
    Bpm     4140.8      [W/kg]      power mass ratio
    m                   [kg]        motor mass
    eta     0.95        [-]         motor efficiency

    Upper Unbounded
    ---------------
    W

    Lower Unbounded
    ---------------
    Pmax

    LaTex Strings
    -------------
    Pmax        P_{\\mathrm{max}}
    Bpm         B_{\\mathrm{PM}}
    eta         \\eta

    """
    def setup(self):
        exec parse_variables(Motor.__doc__)
        return [Pmax <= Bpm*m, W >= m*g]

class Battery(Model):
    """ Battery Model

    Variables
    ---------
    W                               [lbf]        battery weight
    etacharge           0.98        [-]          charging efficiency
    etadischarge        0.98        [-]          discharging efficiency
    E                               [J]          total battery energy
    hbatt               350         [W*hr/kg]    battery specific energy
    vbatt               800         [W*hr/l]     volume battery energy density
    Volbatt                         [m**3]       battery volume

    Upper Unbounded
    ---------------
    W, Volbatt

    Lower Unbounded
    ---------------
    E

    LaTex Strings
    -------------
    eta_charge          \\eta_{\\mathrm{charge}}
    eta_discharge       \\eta_{\\mathrm{discharge}}
    hbatt               h_{\\mathrm{batt}}
    vbatt               (EV)_{\\mathrm{batt}}
    Volbatt             \\mathcal{V}_{\\mathrm{batt}}

    """
    def setup(self):
        exec parse_variables(Battery.__doc__)
        return [W >= E/hbatt*g, Volbatt >= E/vbatt]

class SolarCells(Model):
    """solar cell model

    Variables
    ---------
    rhosolar            0.27    [kg/m^2]        solar cell area density
    S                           [ft**2]         solar cell area
    W                           [lbf]           solar cell weight
    etasolar            0.22    [-]             solar cell efficiency
    mfac                1.0     [-]             solar cell area margin

    Upper Unbounded
    ---------------
    W

    Lower Unbounded
    ---------------
    S

    LaTex Strings
    -------------
    rhosolar        \\rho_{\\mathrm{solar}}
    etasolar        \\eta_{\\mathrm{solar}}
    mfac            m_{\\mathrm{fac}}

    """
    def setup(self):
        exec parse_variables(SolarCells.__doc__)
        return [W >= rhosolar*S*g]

class AircraftPerf(Model):
    """ Aircaft Performance

    Variables
    ---------
    CD                  [-]     aircraft drag coefficient
    cda                 [-]     non-wing drag coefficient
    Pshaft              [hp]    shaft power
    Pavn        0.0     [W]     accessory power draw
    Poper               [W]     operating power
    mfac        1.05    [-]     drag margin factor


    LaTex Strings
    -------------
    CD          C_D
    cda         CDA
    Pshaft      P_{\\mathrm{shaft}}
    Pavn        P_{\\mathrm{avn}}
    Poper       P_{\\mathrm{oper}}
    mfac        m_{\\mathrm{fac}}
    """
    def setup(self, static, state):
        exec parse_variables(AircraftPerf.__doc__)

        fd = dirname(abspath(__file__)) + sep + "dai1336a.csv"

        self.wing = static.wing.flight_model(static.wing, state, fitdata=fd)
        self.htail = static.emp.htail.flight_model(static.emp.htail, state)
        self.vtail = static.emp.vtail.flight_model(static.emp.vtail, state)
        self.tailboom = static.emp.tailboom.flight_model(static.emp.tailboom,
                                                         state)

        self.flight_models = [self.wing, self.htail, self.vtail, self.tailboom]

        e = self.e = self.wing.e
        cdht = self.cdht = self.htail.Cd
        cdvt = self.cdvt = self.vtail.Cd
        Sh = self.Sh = static.Sh
        Sv = self.Sv = static.Sv
        Sw = self.Sw = static.Sw
        cftb = self.cftb = self.tailboom.Cf
        Stb = self.Stb = static.emp.tailboom.S
        E = self.E = static.battery.E
        etacharge = self.etacharge = static.battery.etacharge
        etadischarge = self.etadischarge = static.battery.etadischarge
        etasolar = self.etasolar = static.solarcells.etasolar
        Ssolar = self.Ssolar = static.Ssolar
        etamotor = self.etamotor = static.motor.eta
        Pmax = self.Pmax = static.motor.Pmax
        cdw = self.cdw = self.wing.Cd
        ESirr = self.ESirr = state.ESirr
        ESday = self.ESday = state.ESday
        EStwi = self.EStwi = state.EStwi
        tnight = self.tnight = state.tnight
        PSmin = self.PSmin = state.PSmin

        self.wing.substitutions[e] = 0.95

        dvars = [cdht*Sh/Sw, cdvt*Sv/Sw, cftb*Stb/Sw]

        if static.fuseModel:
            self.fuse = static.fuselage.flight_model(static.fuselage, state)
            self.flight_models.extend([self.fuse])

            cdfuse = self.fuse.Cd
            Sfuse = static.fuselage.S
            dvars.append(cdfuse*Sfuse/Sw)

        constraints = [
            ESirr >= (ESday + E/etacharge/etasolar/Ssolar),
            E*etadischarge >= (Poper*tnight + EStwi*etasolar*Ssolar),
            Poper >= Pavn + Pshaft/etamotor,
            Poper == PSmin*Ssolar*etasolar,
            cda >= sum(dvars),
            CD/mfac >= cda + cdw,
            Poper <= Pmax
            ]

        return self.flight_models, constraints

class FlightState(Model):
    """Flight State (wind speed, solar irradiance, atmosphere)

    Arguments
    ------
    latitude        [deg]       earth latitude
    day                         day of the year [Jan 1st = 1]

    Variables
    ---------
    Vwind                   [m/s]       wind velocity
    V                       [m/s]       true airspeed
    rho                     [kg/m^3]    air density
    mu          1.42e-5     [N*s/m^2]   viscosity
    ESirr       self.esirr  [W*hr/m^2]  solar energy
    PSmin                   [W/m^2]     minimum necessary solar power
    ESday                   [W*hr/m^2]  solar cells energy during daytime
    EStwi                   [W*hr/m^2]  twilight required battery energy
    ESvar       1           [W*hr/m^2]  energy units variable
    PSvar       1           [W/m^2]     power units variable
    tnight      self.tn     [hr]        night duration
    pct         0.9         [-]         percentile wind speeds
    Vwindref    100.0       [m/s]       reference wind speed
    rhoref      1.0         [kg/m^3]    reference air density
    mfac        1.0         [-]         wind speed margin factor
    rhosl       1.225       [kg/m^3]    sea level air density

    LaTex Strings
    -------------
    Vwind       V_{\\mathrm{wind}}}
    V           V
    rho         \\rho
    mu          \\mu
    ESirr       (E/S)_{\\mathrm{irr}}
    PSmin       (P/S)_{\\mathrm{min}}
    ESday       (E/S)_{\\mathrm{day}}
    EStwi       (E/S)_{\\mathrm{twi}}
    ESvar       (E/S)_{\\mathrm{ref}}
    PSvar       (P/S)_{\\mathrm{ref}}
    tnight      t_{\\mathrm{night}}
    pct         p_{\\mathrm{wind}}
    Vwindref    V_{\\mathrm{wind-ref}}
    rhoref      \\rho_{\\mathrm{ref}}
    mfac        m_{\\mathrm{fac}}

    """
    def setup(self, latitude=45, day=355):
        self.esirr, _, self.tn, _ = get_Eirr(latitude, day)
        exec parse_variables(FlightState.__doc__)

        month = get_month(day)
        df = pd.read_csv(path + sep + "windfits" + month +
                         "/windaltfit_lat%d.csv" % latitude).to_dict(
                             orient="records")[0]
        with StdoutCaptured(None):
            dft, dfd = twi_fits(latitude, day, gen=True)

        return [V/mfac >= Vwind,
                FCS(df, Vwind/Vwindref, [rho/rhoref, pct]),
                FCS(dfd, ESday/ESvar, [PSmin/PSvar]),
                FCS(dft, EStwi/ESvar, [PSmin/PSvar]),
               ]

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
    """ Flight Segment

    Variables
    ---------
    Vne         30      [m/s]       never exceed velocity

    """
    def setup(self, aircraft, latitude=35, day=355):
        exec parse_variables(FlightSegment.__doc__)

        self.latitude = latitude
        self.day = day

        self.aircraft = aircraft
        self.fs = FlightState(latitude=latitude, day=day)
        self.aircraftPerf = self.aircraft.flight_model(self.fs)
        self.slf = SteadyLevelFlight(self.fs, self.aircraft,
                                     self.aircraftPerf)

        self.wingg = self.aircraft.wing.spar.loading(self.aircraft.wing)
        self.winggust = self.aircraft.wing.spar.gustloading(self.aircraft.wing)
        self.htailg = self.aircraft.emp.htail.spar.loading(
            self.aircraft.emp.htail)
        self.vtailg = self.aircraft.emp.vtail.spar.loading(
            self.aircraft.emp.vtail)

        self.loading = [self.wingg, self.winggust, self.htailg, self.vtailg]

        self.wingg.substitutions[self.wingg.Nmax] = 5
        self.winggust.substitutions[self.winggust.vgust] = 5
        self.winggust.substitutions[self.winggust.Nmax] = 2

        rhosl = self.fs.rhosl
        Sh = self.aircraft.emp.htail.planform.S
        CLhmax = self.aircraft.emp.htail.planform.CLmax
        Sv = self.aircraft.emp.vtail.planform.S
        CLvmax = self.aircraft.emp.vtail.planform.CLmax

        constraints = [self.aircraft.Wcent == self.wingg.W,
                       self.aircraft.Wcent == self.winggust.W,
                       self.aircraft.Wwing == self.winggust.Ww,
                       self.fs.V == self.winggust.v,
                       self.aircraftPerf.wing.CL == self.winggust.cl,
                       self.htailg.W == 0.5*Vne**2*rhosl*Sh*CLhmax,
                       self.vtailg.W == 0.5*Vne**2*rhosl*Sv*CLvmax
                      ]

        self.submodels = [self.fs, self.aircraftPerf, self.slf, self.loading]

        return constraints, self.submodels

class SteadyLevelFlight(Model):
    """ steady level flight model

    Variables
    ---------
    T                       [N]     thrust
    etaprop         0.8     [-]     propeller efficiency

    """
    def setup(self, state, aircraft, perf):
        exec parse_variables(SteadyLevelFlight.__doc__)

        Wtotal = self.Wtotal = aircraft.Wtotal
        CL = self.CL = perf.wing.CL
        CD = self.CD = perf.CD
        Pshaft = self.Pshaft = perf.Pshaft
        S = self.S = aircraft.wing.planform.S
        rho = self.rho = state.rho
        V = self.V = state.V

        return [Wtotal <= (0.5*rho*V**2*CL*S),
                T >= 0.5*rho*V**2*CD*S,
                Pshaft >= T*V/etaprop]

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

def test():
    " test model for continuous integration "
    m = Mission()
    m.cost = m[m.solar.Wtotal]
    m.solve()
    m = Mission(sp=True)
    m.cost = m[m.solar.Wtotal]
    m.localsolve()

if __name__ == "__main__":
    SP = True
    M = Mission(latitude=[20], sp=SP)
    del M.substitutions[M.solar.wing.planform.tau]
    M.cost = M[M.solar.Wtotal]
    sol = M.localsolve("mosek") if SP else M.solve("mosek")
