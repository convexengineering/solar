" Simple Solar-Electric Powered Aircraft Model "
#pylint: disable=invalid-name, too-many-instance-attributes, too-many-locals
#pylint: disable=redefined-variable-type, too-many-statements, not-callable
from os.path import abspath, dirname
from os import sep
from numpy import hstack
import pandas as pd
import gassolar.environment
from ad.admath import exp
from gassolar.environment.solar_irradiance import get_Eirr, twi_fits
from gassolar.environment.wind_speeds import get_month
from gpkit import Model, parse_variables, Vectorize, SignomialsEnabled
from gpkit.tests.helpers import StdoutCaptured
from gpkitmodels.GP.aircraft.wing.wing import Wing as WingGP
from gpkitmodels.SP.aircraft.wing.wing import Wing as WingSP
from gpkitmodels.GP.aircraft.wing.boxspar import BoxSpar as BoxSparGP
from gpkitmodels.SP.aircraft.wing.boxspar import BoxSpar as BoxSparSP
from gpkitmodels.GP.aircraft.wing.wing_skin import WingSecondStruct
from gpkitmodels.GP.aircraft.tail.empennage import Empennage
from gpkitmodels.GP.aircraft.tail.horizontal_tail import HorizontalTail
from gpkitmodels.GP.aircraft.tail.vertical_tail import VerticalTail
from gpkitmodels.GP.aircraft.tail.tail_boom import TailBoom
from gpkitmodels.SP.aircraft.tail.tail_boom_flex import TailBoomFlexibility
from gpkitmodels.GP.materials import cfrpud, cfrpfabric, foamhd
from gpkitmodels.GP.aircraft.fuselage.elliptical_fuselage import Fuselage
from gpkitmodels.GP.aircraft.prop.propeller import Propeller, ActuatorProp
from gpkitmodels.SP.aircraft.prop.propeller import BladeElementProp
from gpkitmodels.GP.aircraft.motor.motor import Motor
from gpkitmodels import g
from gpfit.fit_constraintset import FitCS as FCS
from relaxed_constants import relaxed_constants, post_process

path = dirname(gassolar.environment.__file__)

class AircraftPerf(Model):
    "Aircaft Performance"

    def setup(self, static, state, onDesign=False):
        self.drag = AircraftDrag(static, state, onDesign)
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
            Poper == PSmin*Ssolar*etasolar]

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
    T                   [lbf]     thrust

    LaTex Strings
    -------------
    CD          C_D
    cda         CDA
    mfac        m_{\\mathrm{fac}}
    """
    @parse_variables(__doc__, globals())
    def setup(self, static, state, onDesign=False):
        fd = dirname(abspath(__file__)) + sep + "dai1336a.csv"

        self.wing = static.wing.flight_model(static.wing, state, fitdata=fd)
        self.htail = static.emp.htail.flight_model(static.emp.htail, state)
        self.vtail = static.emp.vtail.flight_model(static.emp.vtail, state)
        self.tailboom = static.emp.tailboom.flight_model(static.emp.tailboom,
                                                         state)
        self.motor = static.motor.flight_model(static.motor, state)
        if static.sp:
            if onDesign:
                static.propeller.flight_model = BladeElementProp

        self.propeller = static.propeller.flight_model(static.propeller, state)

        self.flight_models = [self.wing, self.htail, self.vtail,
                              self.tailboom, self.motor, self.propeller]

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
        Nprop = static.Nprop
        Tprop = self.propeller.T
        Qprop = self.propeller.Q
        RPMprop = self.propeller.omega
        Qmotor = self.motor.Q
        RPMmotor = self.motor.omega
        Pelec = self.motor.Pelec


        self.wing.substitutions[e] = 0.95
        self.wing.substitutions[self.wing.CLstall] = 4

        self.wing.substitutions[e] = 0.95
        dvars = [cdht*Sh/Sw, cdvt*Sv/Sw, cftb*Stb/Sw]

        if static.Npod is not 0:
            with Vectorize(static.Npod):
                self.fuse = static.fuselage.flight_model(static.fuselage,
                                                         state)
            self.flight_models.extend([self.fuse])
            cdfuse = self.fuse.Cd
            Sfuse = static.fuselage.S
            dvars.extend(cdfuse*Sfuse/Sw)
            self.fuse.substitutions[self.fuse.mfac] = 1.1

        constraints = [cda >= sum(dvars),
                       Tprop == T/Nprop,
                       Qmotor == Qprop,
                       RPMmotor == RPMprop,
                       CD/mfac >= cda + cdw,
                       Poper/mpower >= Pavn + Ppay + (Pelec*Nprop),
                      ]

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
    Nprop       4       [-]     Number of propulsors
    minvttau    0.09    [-]     minimum vertical tail tau ratio
    minhttau    0.06    [-]     minimum horizontal tail tau ratio
    maxtau      0.144   [-]     maximum wing tau ratio

    SKIP VERIFICATION

    Upper Unbounded
    ---------------
    Wwing, Wcent, wing.mw (if sp), propeller.c

    Lower Unbounded
    ---------------
    wing.spar.J, wing.spar.Sy
    emp.htail.spar.J, emp.htail.spar.Sy
    emp.vtail.spar.J, emp.vtail.spar.Sy
    emp.tailboom.J, emp.tailboom.Sy
    motor.Qmax
    battery.E, solarcells.S
    emp.htail.mh (if sp), emp.htail.Vh (if sp)
    propeller.R, propeller.T_m, propeller.c

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

    @parse_variables(__doc__, globals())
    def setup(self, Npod=0, sp=False):
        self.Npod = Npod
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

        HorizontalTail.sparModel = BoxSparGP
        HorizontalTail.fillModel = None
        HorizontalTail.skinModel = WingSecondStruct
        VerticalTail.sparModel = BoxSparGP
        VerticalTail.fillModel = None
        VerticalTail.skinModel = WingSecondStruct
        TailBoom.__bases__ = (BoxSparGP,)
        TailBoom.secondaryWeight = True
        self.emp = Empennage(N=5)
        self.solarcells = SolarCells()
        self.battery = Battery()
        if sp:
            WingSP.sparModel = BoxSparSP
            WingSP.fillModel = None
            WingSP.skinModel = WingSecondStruct
            self.wing = WingSP(N=20)
        else:
            WingGP.sparModel = BoxSparGP
            WingGP.fillModel = None
            WingGP.skinModel = WingSecondStruct
            self.wing = WingGP(N=20)
        self.motor = Motor()
        Propeller.flight_model = ActuatorProp
        self.propeller = Propeller()
        self.components = [self.solarcells, self.wing, self.battery,
                           self.emp]
        self.propulsor = [self.motor, self.propeller]

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
        Ssolar = self.Ssolar = self.solarcells.S
        mfsolar = self.mfsolar = self.solarcells.mfac
        Volbatt = self.battery.Volbatt
        vttau = self.emp.vtail.planform.tau
        httau = self.emp.htail.planform.tau

        self.emp.substitutions[Vv] = 0.02
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
                       vttau >= minvttau,
                       httau >= minhttau,
                       tau <= maxtau
                      ]

        if self.Npod is not 0:
            with Vectorize(1):
                with Vectorize(self.Npod):
                    self.fuselage = Fuselage()

            self.k = self.fuselage.k
            Volfuse = self.Volfuse = self.fuselage.Vol[:, 0]
            Wbatt = self.battery.W
            Wfuse = sum(self.fuselage.W)
            self.fuselage.substitutions[self.fuselage.nply] = 5

            constraints.extend([
                Volbatt <= Volfuse,
                Wwing >= self.wing.W + self.solarcells.W,
                Wcent >= (Wpay + Wavn + self.emp.W + self.motor.W*Nprop
                          + self.fuselage.W[0] + Wbatt/self.Npod),
                Wtotal/mfac >= (Wpay + Wavn + Wland + Wfuse +
                                sum([c.W for c in self.components]) +
                                (Nprop)*sum([c.W for c in self.propulsor]))
                ])

            self.components.append(self.fuselage)
        else:
            constraints.extend([
                Wwing >= sum([c.W for c in [self.wing, self.battery,
                                            self.solarcells]]),
                Wcent >= Wpay + Wavn + self.emp.W + self.motor.W*Nprop,
                Volbatt <= cmac**2*0.5*tau*b,
                Wtotal/mfac >= (Wpay + Wavn + Wland
                                + sum([c.W for c in self.components])
                                + Nprop*sum([c.W for c in self.propulsor]))
                ])

        return constraints, self.components, materials, self.propulsor


class Battery(Model):
    """ Battery Model

    Variables
    ---------
    W                               [lbf]        battery weight
    etacharge           0.98        [-]          charging efficiency
    etadischarge        0.98        [-]          discharging efficiency
    E                               [kJ]         total battery energy
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
    @parse_variables(__doc__, globals())
    def setup(self):
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
    @parse_variables(__doc__, globals())
    def setup(self):
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
    ESirr       esirr       [W*hr/m^2]  solar energy
    PSmin                   [W/m^2]     minimum necessary solar power
    ESday                   [W*hr/m^2]  solar cells energy during daytime
    EStwi                   [W*hr/m^2]  twilight required battery energy
    ESvar       1           [W*hr/m^2]  energy units variable
    PSvar       1           [W/m^2]     power units variable
    tnight      tn          [hr]        night duration
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
    @parse_variables(__doc__, globals())
    def setup(self, latitude, day, esirr, tn):
        self.esirr = esirr
        self.tn = tn

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


class FlightSegment(Model):
    """ Flight Segment
    """
    def setup(self, aircraft, latitude=35, day=355):
        # exec parse_variables(FlightSegment.__doc__)

        self.latitude = latitude
        self.day = day
        esirr, _, tn, _ = get_Eirr(latitude, day)

        self.aircraft = aircraft
        self.fs = FlightState(latitude, day, esirr, tn)
        self.aircraftPerf = self.aircraft.flight_model(aircraft, self.fs, False)
        self.slf = SteadyLevelFlight(self.fs, self.aircraft,
                                     self.aircraftPerf)

        if aircraft.Npod is not 0 and aircraft.Npod is not 1:
            assert self.aircraft.sp
            loadsp = self.aircraft.sp
        else:
            loadsp = False

        self.wingg = self.aircraft.wing.spar.loading(
            self.aircraft.wing, self.fs, out=loadsp)
        self.winggust = self.aircraft.wing.spar.gustloading(
            self.aircraft.wing, self.fs, out=loadsp)
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

        constraints = [
            self.aircraft.Wcent == self.wingg.W,
            self.aircraft.Wcent == self.winggust.W,
            self.aircraft.Wwing == self.winggust.Ww,
            self.fs.V == self.winggust.v,
            self.aircraftPerf.CL == self.winggust.cl,
            self.htailg.W == qne*Sh*CLhmax,
            self.vtailg.W == qne*Sv*CLvmax,
            ]

        if self.aircraft.Npod is not 0 and self.aircraft.Npod is not 1:
            Nwing, Npod = self.aircraft.wing.N, self.aircraft.Npod
            ypod = Nwing/((Npod-1)/2 + 1)
            ypods = [ypod*n for n in range(1, (Npod-1)/2+1)]
            Sgust, Mgust = self.winggust.S, self.winggust.M
            qgust, Sg, Mg = self.winggust.q, self.wingg.S, self.wingg.M
            qg = self.wingg.q
            deta = self.aircraft.wing.planform.deta
            b = self.aircraft.wing.planform.b
            weight = self.aircraft.battery.W/Npod*self.wingg.N
            for i in range(Nwing-1):
                if i in ypods:
                    with SignomialsEnabled():
                        constraints.extend([
                            Sgust[i] >= (Sgust[i+1] + 0.5*deta[i]*(b/2)
                                         * (qgust[i] + qgust[i+1]) - weight),
                            Sg[i] >= (Sg[i+1] + 0.5*deta[i]*(b/2)
                                      * (qg[i] + qg[i+1]) - weight),
                            Mgust[i] >= (Mgust[i+1] + 0.5*deta[i]*(b/2)
                                         * (Sgust[i] + Sgust[i+1])),
                            Mg[i] >= (Mg[i+1] + 0.5*deta[i]*(b/2)
                                      * (Sg[i] + Sg[i+1]))
                            ])
                else:
                    constraints.extend([
                        Sgust[i] >= (Sgust[i+1] + 0.5*deta[i]*(b/2)
                                     * (qgust[i] + qgust[i+1])),
                        Sg[i] >= Sg[i+1] + 0.5*deta[i]*(b/2)*(qg[i] + qg[i+1]),
                        Mgust[i] >= (Mgust[i+1] + 0.5*deta[i]*(b/2)
                                     * (Sgust[i] + Sgust[i+1])),
                        Mg[i] >= Mg[i+1] + 0.5*deta[i]*(b/2)*(Sg[i] + Sg[i+1])
                        ])

        self.submodels = [self.fs, self.aircraftPerf, self.slf, self.loading]

        return constraints, self.submodels

class Climb(Model):
    """ Climb model

    Variables
    ---------
    h           60000           [ft]            climb altitude
    t           500             [min]           time to climb
    hdotmin                     [ft/min]        minimum climb rate
    mu          1.42e-5         [N*s/m^2]       viscosity
    dh          self.hstep      [ft]            change in altitude

    Variables of length [1,N]
    ---------------------
    dt                          [min]           time step
    V                           [m/s]           vehicle speed
    hdot                        [ft/min]        climb rate
    rho         self.density    [kg/m^3]        air density
    """
    def density(self, c):
        " find air density "
        alpha = 0.0065 # K/m
        h11k, T11k, p11k, rhosl = 11019, 216.483, 22532, 1.225 #m, K, Pa, kg/m^3
        T0, R, gms, n = 288.16, 287.04, 9.81, 5.2561 #K, m^2/K/s^2, m/s^2, -
        hrange = [c(self.h).to("m").magnitude*i/(self.N+1)
                  for i in range(1, self.N+1)]
        rho = []
        for al in hrange:
            if al < h11k:
                T = T0 - alpha*al
                rho.append(rhosl*(T/T0)**(n-1))
            else:
                p = p11k*exp((h11k - al)*gms/R/T11k)
                rho.append(p/R/T11k)
        return [rho]

    def hstep(self, c):
        " find delta altitude "
        return c[self.h]/self.N

    @parse_variables(__doc__, globals())
    def setup(self, N, aircraft):
        self.N = N

        with Vectorize(self.N):
            self.drag = AircraftDrag(aircraft, self)


        Wtotal = self.Wtotal = aircraft.Wtotal
        CD = self.CD = self.drag.CD
        CL = self.CL = self.drag.CL
        S = self.S = aircraft.wing.planform.S
        E = aircraft.battery.E
        Poper = self.drag.Poper
        T = self.drag.T
        self.rho = rho

        constraints = [
            Wtotal <= 0.5*rho*V**2*CL*S,
            T >= 0.5*rho*V**2*CD*S + Wtotal*hdot/V,
            hdot >= dh/dt,
            t >= sum(hstack(dt)),
            E >= sum(hstack(Poper*dt))]

        return self.drag, constraints

class SteadyLevelFlight(Model):
    """ steady level flight model

    """
    def setup(self, state, aircraft, perf):
        Wtotal = self.Wtotal = aircraft.Wtotal
        CL = self.CL = perf.CL
        CD = self.CD = perf.CD
        S = self.S = aircraft.wing.planform.S
        rho = self.rho = state.rho
        V = self.V = state.V
        T = perf.drag.T

        return [Wtotal <= (0.5*rho*V**2*CL*S),
                T >= 0.5*rho*V**2*CD*S]

class Mission(Model):
    "define mission for aircraft"
    def setup(self, aircraft, latitude=range(1, 21, 1), day=355):

        self.aircraft = aircraft
        self.mission = []
        self.mission.append(Climb(5, self.aircraft))
        if day == 355 or day == 172:
            for l in latitude:
                self.mission.append(FlightSegment(self.aircraft, l, day))
        else:
            assert day < 172
            for l in latitude:
                self.mission.append(FlightSegment(self.aircraft, l, day))
                self.mission.append(FlightSegment(self.aircraft, l,
                                                  355 - 10 - day))

        return self.mission, self.aircraft

def test():
    " test model for continuous integration "
    v = Aircraft(sp=False)
    m = Mission(v, latitude=[20])
    m.cost = m[m.aircraft.Wtotal]
    m.solve()
    v = Aircraft(sp=True)
    m = Mission(v, latitude=[20])
    m.cost = m[m.aircraft.Wtotal]
    m.localsolve()
    v = Aircraft(Npod=3, sp=True)
    m = Mission(v, latitude=[20])
    m.cost = m[m.aircraft.Wtotal]
    f = relaxed_constants(M)
    s = f.localsolve("mosek")
    post_process(s)

if __name__ == "__main__":
    SP = True
    Vehicle = Aircraft(Npod=3, sp=SP)
    M = Mission(Vehicle, latitude=[20])
    M.cost = M[M.aircraft.Wtotal]
    try:
        sol = (M.localsolve("mosek") if SP else M.solve("mosek"))
    except RuntimeWarning:
        V2 = Aircraft(Npod=3, sp=SP)
        M2 = Mission(V2, latitude=[20])
        M2.cost = M2[M2.aircraft.Wtotal]
        feas = relaxed_constants(M2)
        sol = feas.localsolve("mosek")
        vks = post_process(sol)
