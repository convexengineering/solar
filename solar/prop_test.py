from gpkitmodels.GP.aircraft.prop.propeller import Propeller
from gpkitmodels.GP.aircraft.wing.wing_test import FlightState


fs = FlightState()
p = Propeller(fs)
p.substitutions[p.T] = 1000
p.cost = 1/p.eta
sol = p.solve()
print sol.table()