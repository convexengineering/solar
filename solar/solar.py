" Simple Solar-Electric Powered Aircraft Model "
#pylint: disable=invalid-name, too-many-instance-attributes, too-many-locals
#pylint: disable=redefined-variable-type, too-many-statements, not-callable
from os.path import abspath, dirname
from os import sep
import pandas as pd
import numpy as np
import gassolar.environment
from gassolar.environment.solar_irradiance import get_Eirr, twi_fits
from gassolar.environment.wind_speeds import get_month
from gpkit import Model, parse_variables, Variable
from gpkit.tests.helpers import StdoutCaptured
from gpkitmodels.GP.aircraft.wing.wing import Wing as WingGP
from gpkitmodels.SP.aircraft.wing.wing import Wing as WingSP
from gpkitmodels.GP.aircraft.wing.boxspar import BoxSpar
from gpkitmodels.GP.aircraft.wing.wing_skin import WingSecondStruct
from gpkitmodels.GP.aircraft.tail.empennage import Empennage
from gpkitmodels.GP.aircraft.tail.horizontal_tail import HorizontalTail
from gpkitmodels.GP.aircraft.tail.vertical_tail import VerticalTail
from gpkitmodels.GP.aircraft.tail.tail_boom import TailBoom
from gpkitmodels.SP.aircraft.tail.tail_boom_flex import TailBoomFlexibility
from gpkitmodels import g
from gpfit.fit_constraintset import FitCS as FCS
from gpkitmodels.GP.materials import cfrpud, cfrpfabric, foamhd


path = dirname(gassolar.environment.__file__)

class AircraftPerf(Model):
    """ Aircaft Performance
    """
    def setup(self, static, state):
        exec parse_variables(AircraftPerf.__doc__)

        self.drag = AircraftDrag(static, state)
        self.CD = self.drag.CD
        self.CL = self.drag.CL
        self.Pshaft = self.drag.Pshaft
        Poper = self.drag.Poper
        E = self.E = static.battery.E
        etacharge = self.etacharge = static.battery.etacharge
        etadischarge = self.etadischarge = static.battery.etadischarge
        etasolar = self.etasolar = static.solarcells.etasolar
        Ssolar = self.Ssolar = static.Ssolar
        ESirr = self.ESirr = state.ESirr
        ESday = self.ESday = state.ESday
        EStwi = self.EStwi = state.EStwi
        tnight = self.tnight = state.tnight
        PSmin = self.PSmin = state.PSmin

        constraints = [
            ESirr >= (ESday + E/etacharge/etasolar/Ssolar),
            E*etadischarge >= (Poper*tnight + EStwi*etasolar*Ssolar),
            Poper == PSmin*Ssolar*etasolar,
            ]

        return self.drag, constraints

class AircraftDrag(Model):
    """ Aircaft Performance

    Variables
    ---------
    CD                  [-]     aircraft drag coefficient
    cda                 [-]     non-wing drag coefficient
    mfac        1.05    [-]     drag margin factor
    Pshaft              [hp]    shaft power
    Pavn        200     [W]     avionics power draw
    Ppay        100     [W]     payload power draw
    Poper               [W]     operating power
    mpower      1.05    [-]     power margin

    LaTex Strings
    -------------
    CD          C_D
    cda         CDA
    mfac        m_{\\mathrm{fac}}
    """
    def setup(self, static, state):
        exec parse_variables(AircraftDrag.__doc__)

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
        cdw = self.cdw = self.wing.Cd
        self.CL = self.wing.CL
        etamotor = self.etamotor = static.motor.eta
        Pmax = self.Pmax = static.motor.Pmax
        self.wing.substitutions[e] = 0.95

        dvars = [cdht*Sh/Sw, cdvt*Sv/Sw, cftb*Stb/Sw]

        if static.fuseModel:
            self.fuse = static.fuselage.flight_model(static.fuselage, state)
            self.flight_models.extend([self.fuse])

            cdfuse = self.fuse.Cd
            Sfuse = static.fuselage.S
            dvars.append(cdfuse*Sfuse/Sw)

        constraints = [cda >= sum(dvars),
                       CD/mfac >= cda + cdw,
                       Poper/mpower >= Pavn + Ppay + Pshaft/etamotor,
                       Poper <= Pmax]

        return self.flight_models, constraints

class Aircraft(Model):
    """ Aircraft Model

    Variables
    ---------
    Wpay        11      [lbf]   payload weight
    Wavn        22      [lbf]   avionics weight
    Wtotal              [lbf]   aircraft weight
    Wwing               [lbf]   wing weight
    Wcent               [lbf]   center weight
    mfac        1.05    [-]     total weight margin
    fland       0.02    [-]     fractional landing gear weight
    Wland               [lbf]   landing gear weight

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
    flight_model = AircraftPerf

    def setup(self, sp=False):
        exec parse_variables(Aircraft.__doc__)

        self.sp = sp

        cfrpud.substitutions.update({cfrpud.rho: 1.5,
                                     cfrpud.E: 200,
                                     cfrpud.tmin: 0.1,
                                     cfrpud.sigma: 1500})
        cfrpfabric.substitutions.update({cfrpfabric.rho: 1.3,
                                         cfrpfabric.E: 40,
                                         cfrpfabric.tmin: 0.1,
                                         cfrpfabric.sigma: 300,
                                         cfrpfabric.tau: 80})
        foamhd.substitutions.update({foamhd.rho: 0.03})
        materials = [cfrpud, cfrpfabric, foamhd]

        HorizontalTail.sparModel = BoxSpar
        HorizontalTail.fillModel = None
        HorizontalTail.skinModel = WingSecondStruct
        VerticalTail.sparModel = BoxSpar
        VerticalTail.fillModel = None
        VerticalTail.skinModel = WingSecondStruct
        TailBoom.__bases__ = (BoxSpar,)
        TailBoom.secondaryWeight = True
        self.emp = Empennage(N=5)
        self.solarcells = SolarCells()
        if sp:
            WingSP.sparModel = BoxSpar
            WingSP.fillModel = None
            WingSP.skinModel = WingSecondStruct
            self.wing = WingSP()
        else:
            WingGP.sparModel = BoxSpar
            WingGP.fillModel = None
            WingGP.skinModel = WingSecondStruct
            self.wing = WingGP()
        self.battery = Battery()
        self.motor = Motor()

        self.components = [self.solarcells, self.wing, self.battery,
                           self.emp, self.motor]

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
        self.emp.substitutions[self.emp.htail.skin.rhoA] = 0.4
        self.emp.substitutions[self.emp.vtail.skin.rhoA] = 0.4
        self.emp.substitutions[self.emp.tailboom.wlim] = 1.0
        self.wing.substitutions[self.wing.mfac] = 1.0
        if not sp:
            self.emp.substitutions[Vh] = 0.45
            self.emp.substitutions[self.emp.htail.mh] = 0.1


        constraints = [Ssolar*mfsolar <= Sw,
                       Vh <= Sh*lh/Sw/cmac,
                       Vv <= Sv*lv/Sw/b,
                       d0 <= tau*croot,
                       Wland >= fland*Wtotal,
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

        constraints.extend([Wtotal/mfac >= (
            Wpay + Wavn + Wland + sum([c.W for c in self.components]))])

        return constraints, self.components, materials

class Motor(Model):
    """ Motor Model

    Variables
    ---------
    W                   [lbf]       motor weight
    Pmax                [W]         max power
    Bpm     1/0.0003    [W/kg]      power mass ratio
    m                   [kg]        motor mass
    eta                 [-]         motor system efficiency
    etam    0.95        [-]         motor efficiency
    etac    0.97        [-]         controller efficiency

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
        return [Pmax <= Bpm*m, W >= m*g, eta <= etam*etac]

class Battery(Model):
    """ Battery Model

    Variables
    ---------
    W                               [lbf]        battery weight
    etacharge           0.98        [-]          charging efficiency
    etadischarge        0.98        [-]          discharging efficiency
    E                               [J]          total battery energy
    hbatt               350         [W*hr/kg]    battery specific energy
    vbatt               800         [W*hr/l]     battery energy density
    Volbatt                         [m**3]       battery volume
    etapack             0.85        [-]          packing efficiency
    etaRTE              0.95        [-]          battery RTE
    minSOC              1.03        [-]          minimum state of charge
    rhomppt             0.4223      [kg/kW]      power system mass density
    etamppt             0.975       [-]          power system efficiency

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
        return [W >= E*minSOC/hbatt/etaRTE/etapack*g,
                Volbatt >= E/vbatt]

class SolarCells(Model):
    """solar cell model

    Variables
    ---------
    rhosolar            0.3     [kg/m^2]        solar cell area density
    S                           [ft**2]         solar cell area
    W                           [lbf]           solar cell weight
    etasolar            0.2     [-]             solar cell efficiency
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
    Vne                     [m/s]       never exceed speed at altitude
    qne                     [kg/s^2/m]  never exceed dynamic pressure
    N           1.4         [-]         factor on Vne

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
                FCS(df, Vwind/Vwindref, [rho/rhoref, pct], name="wind"),
                FCS(dfd, ESday/ESvar, [PSmin/PSvar]),
                FCS(dft, EStwi/ESvar, [PSmin/PSvar]),
                Vne == N*V,
                qne == 0.5*rho*Vne**2
               ]

def altitude(density):
    " find air density "
    g1 = 9.80665 # m/s^2
    R = 287.04 # m^2/K/s^2
    T11 = 216.65 # K
    p11 = 22532 # Pa
    p = density*R*T11
    h = (11000 - R*T11/g1*np.log(p/p11))/0.3048
    return h

class FlightSegment(Model):
    """ Flight Segment
    """
    def setup(self, aircraft, latitude=35, day=355):
        exec parse_variables(FlightSegment.__doc__)

        self.latitude = latitude
        self.day = day

        self.aircraft = aircraft
        self.fs = FlightState(latitude=latitude, day=day)
        self.aircraftPerf = self.aircraft.flight_model(aircraft, self.fs)
        self.slf = SteadyLevelFlight(self.fs, self.aircraft,
                                     self.aircraftPerf)

        self.wingg = self.aircraft.wing.spar.loading(
            self.aircraft.wing, self.fs)
        self.winggust = self.aircraft.wing.spar.gustloading(
            self.aircraft.wing, self.fs)
        self.htailg = self.aircraft.emp.htail.spar.loading(
            self.aircraft.emp.htail, self.fs)
        self.vtailg = self.aircraft.emp.vtail.spar.loading(
            self.aircraft.emp.vtail, self.fs)

        self.tbhbend = self.aircraft.emp.tailboom.tailLoad(
            self.aircraft.emp.tailboom, self.aircraft.emp.htail,
            self.fs)
        self.tbvbend = self.aircraft.emp.tailboom.tailLoad(
            self.aircraft.emp.tailboom, self.aircraft.emp.vtail,
            self.fs)

        self.loading = [self.wingg, self.winggust, self.htailg, self.vtailg,
                        self.tbhbend, self.tbvbend]

        if self.aircraft.sp:
            self.tbflex = TailBoomFlexibility(self.aircraft.emp.htail,
                                              self.tbhbend, self.aircraft.wing)
            self.tbflex.substitutions[self.tbflex.SMcorr] = 0.05
            self.loading.append(self.tbflex)

        self.wingg.substitutions[self.wingg.Nmax] = 2
        self.wingg.substitutions[self.wingg.Nsafety] = 1.5
        self.winggust.substitutions[self.winggust.vgust] = 5
        self.winggust.substitutions[self.winggust.Nmax] = 2
        self.winggust.substitutions[self.winggust.Nsafety] = 1.5
        self.tbhbend.substitutions[self.tbhbend.Nsafety] = 1.5
        self.tbvbend.substitutions[self.tbvbend.Nsafety] = 1.5

        Sh = self.aircraft.emp.htail.planform.S
        CLhmax = self.aircraft.emp.htail.planform.CLmax
        Sv = self.aircraft.emp.vtail.planform.S
        CLvmax = self.aircraft.emp.vtail.planform.CLmax
        qne = self.fs.qne

        constraints = [self.aircraft.Wcent == self.wingg.W,
                       self.aircraft.Wcent == self.winggust.W,
                       self.aircraft.Wwing == self.winggust.Ww,
                       self.fs.V == self.winggust.v,
                       self.aircraftPerf.CL == self.winggust.cl,
                       self.htailg.W == qne*Sh*CLhmax,
                       self.vtailg.W == qne*Sv*CLvmax
                      ]

        self.submodels = [self.fs, self.aircraftPerf, self.slf, self.loading]

        return constraints, self.submodels

class Climb(Model):
    """ Climb model

    Variables
    ---------
    h           60000       [ft]            climb altitude
    t           500         [min]           time to climb
    hdotmin                 [ft/min]        minimum climb rate
    hdot                    [ft/min]        climb rate
    T                       [N]             thrust to climb
    rho         0.003097    [kg/m^3]        air density
    V                       [m/s]           vehicle speed
    mu          1.42e-5     [N*s/m^2]       viscosity
    etaprop     0.85        [-]             propeller efficiency

    """

    def setup(self, aircraft):
        exec parse_variables(Climb.__doc__)

        self.drag = AircraftDrag(aircraft, self)
        Wtotal = self.Wtotal = aircraft.Wtotal
        CD = self.CD = self.drag.CD
        CL = self.CL = self.drag.CL
        S = self.S = aircraft.wing.planform.S
        Pshaft = self.drag.Pshaft

        constraints = [
            Wtotal <= 0.5*rho*V**2*CL*S,
            T >= 0.5*rho*V**2*CD*S + Wtotal*hdot/V,
            hdot >= hdotmin,
            hdotmin == h/t,
            Pshaft >= T*V/etaprop]

        return self.drag, constraints

class SteadyLevelFlight(Model):
    """ steady level flight model

    Variables
    ---------
    T                       [N]     thrust
    etaprop         0.85    [-]     propeller efficiency

    """
    def setup(self, state, aircraft, perf):
        exec parse_variables(SteadyLevelFlight.__doc__)

        Wtotal = self.Wtotal = aircraft.Wtotal
        CL = self.CL = perf.CL
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
        self.mission.append(Climb(self.solar))
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
    M.cost = M[M.solar.Wtotal]
    sol = M.localsolve("mosek") if SP else M.solve("mosek")
