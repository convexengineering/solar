from solar import Aircraft, Mission
from gpkitmodels.GP.aircraft.fuselage.elliptical_fuselage import Fuselage

Aircraft.fuseModel = Fuselage
m = Mission(latitude=[20])
sol = m.solve()
