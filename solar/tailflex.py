" tail boom flexibility trade "
from solar import Mission
from print_sens import sol_table

Mg = Mission(latitude=25, sp=False)
solg = Mg.solve("mosek")

Ms = Mission(latitude=25, sp=True)
sols = Ms.localsolve("mosek")

