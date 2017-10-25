" tail boom flexibility trade "
from solar import Mission
from print_sens import sol_table

Mg = Mission(latitude=25, sp=False)
Mg.cost = Mg["W_{total}"]
solg = Mg.solve("mosek")

Ms = Mission(latitude=25, sp=True)
Ms.cost = Ms["W_{total}"]
sols = Ms.localsolve("mosek")

