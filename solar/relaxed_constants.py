from __future__ import print_function
from builtins import str
from gpkit.constraints.relax import ConstantsRelaxed
from gpkit.constraints.bounded import Bounded
from gpkit import Model

"""
Methods to precondition an SP so that it solves with a relaxed constants algorithim
and postcondition an SP to ensure all relax values are 1
"""

def relaxed_constants(model, include_only=None, exclude=None):
    """
    Method to precondition an SP so it solves with a relaxed constants
        algorithim

    ARGUMENTS
    ---------
    model: the model to solve with relaxed constants

    RETURNS
    -------
    feas: the input model but with relaxed constants and a new objective
    """

    if model.substitutions:
        constsrelaxed = ConstantsRelaxed(Bounded(model))
        feas = Model(constsrelaxed.relaxvars.prod()**20 * model.cost,
                     constsrelaxed)
        # NOTE: It hasn't yet been seen but might be possible that
        #       the model.cost component above could cause infeasibility
    else:
        feas = Model(model.cost, model)

    return feas

def post_process(sol):
    """
    Model to print relevant info for a solved model with relaxed constants

    ARGUMENTS
    --------
    sol: the solution to the solved model
    """

    bdvars = [d for d in sol.program.varkeys if "Relax" in str(d)
              and "before" not in str(d) and sol(d).magnitude >= 1.001]
    if bdvars:
        print("GP iteration has relaxed constants")
        print(sol.program.result.table(varkeys))

    return bdvars

