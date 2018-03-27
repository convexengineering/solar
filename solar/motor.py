"Electric motor model "
from numpy import pi
from gpkit import Model, parse_variables, SignomialsEnabled, SignomialEquality, units
from gpkitmodels.GP.aircraft.prop.propeller import Propeller
from gpkitmodels import g
from gpkitmodels.GP.aircraft.wing.wing_test import FlightState
import matplotlib.pyplot as plt
class ElecMotor(Model):
    """ Electric Motor Model

    Variables
    ---------
    Qstar       1           [kg/(N*m)]         motor specific torque
    W                       [lbf]              motor weight
    Qmax                    [N*m]              motor max. torque
    V_max       44          [V]                motor max voltage

    """
    def setup(self):
        exec parse_variables(ElecMotor.__doc__)

        constraints = [W >= Qstar*Qmax*g]

        return constraints

    def flight_model(self,state):
        return ElecMotor_Performance(self, state)

class ElecMotor_Performance(Model):
    """ Electric Motor Performance Model

    Variables
    ---------
    Pshaft                  [kW]            motor output shaft power
    Pelec                   [kW]            motor input shaft power
    etam                    [-]             motor efficiency
    Q                       [N*m]           torque
    omega                   [rad/s]         propeller rotation rate 
    i                       [amps]          current
    v                       [V]             woltage
    i0         .001           [amps]          zero-load current
    Kv         20           [rad/s/V]       motor voltage constant
    R          .05             [ohms]          internal resistance
    """
    def setup(self,parent,  state):
        exec parse_variables(ElecMotor_Performance.__doc__)


        constraints = [Pshaft == Q*omega,
                Pelec == v*i,
                etam == Pshaft/Pelec, 
                parent.Qmax >= Q,
                v <= parent.V_max,
                i >= Q*Kv+i0,
                v >= omega/Kv + i*R
                ]
        return constraints

class Propulsor(Model):
    """Propulsor model

    Variables
    ---------
    W                       [lbf]              propulsor weight


    """
    def setup(self):
        exec parse_variables(Propulsor.__doc__)

        self.prop = Propeller()
        self.motor = ElecMotor()

        components = [self.prop, self.motor]

        return [self.W >= self.prop.W + self.motor.W], components

    def flight_model(self,state):
        return Propulsor_Performance(self,state)

class Propulsor_Performance(Model):
    """Propulsor Performance Model

    """
    def setup(self, parent,state):

        self.prop    = parent.prop.flight_model(state)
        self.motor   = parent.motor.flight_model(state)

        self.components = [self.prop, self.motor]

        constraints = [self.prop.Q == self.motor.Q,
                        self.prop.omega == self.motor.omega
                        ]

        return constraints, self.components

class Propulsor_Test(Model):
    """Propulsor Test Model
    """

    def setup(self):
        fs = FlightState()
        p = Propulsor()
        pp = p.flight_model(fs)
        pp.substitutions[pp.prop.T] = 20
        self.cost = 1./pp.motor.etam + p.W/(100*units('lbf'))

        return fs,p,pp

def propulsor_test():

    test = Propulsor_Test()
    #sol = test.debug()
    sol = test.solve()
    print sol.table()

class Motor_P_Test(Model):
    def setup(self):
        fs = FlightState()
        m  = ElecMotor()
        mp = ElecMotor_Performance(m,fs)
        self.mp = mp
        mp.substitutions[m.Qmax] = 100
        mp.substitutions[mp.Q]    = 10
        self.cost = 1./mp.etam
        return self.mp, fs

def motor_test():
    test = Motor_P_Test()
    sol = test.solve()
    #sol = test.debug()

    print sol.table()
    
def motor_eta_speed():
    test = Motor_P_Test()
    omega = [.01, 10, 100, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 3000, 4000, 10000, 100000]
    eta = []
    for o in omega:
        test = Motor_P_Test()
        test.substitutions[test.mp.omega] = o
        sol = test.solve()
        eta.append(sol["freevariables"]["etam"])

    print omega
    print eta
    #sol = test.debug()
    plt.plot(omega, eta)
    plot.show()
    print sol.table()
    
if __name__ == "__main__":
    #motor_eta_speed()
    #motor_test()
    propulsor_test()
