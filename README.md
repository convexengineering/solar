# Solar GP Model

This repo contains a solar aircraft sizing model that leverages geometric programming (GP).  Provisions are added to allow signomial programming (SP) for higher fidelity results.  All changes are unit tested on continuous integration software.  Code is python dependent. 

[![Build Status](https://acdl.mit.edu/csi/buildStatus/icon?job=gpkit_ResearchModel_solar_Push)](https://acdl.mit.edu/csi/job/gpkit_ResearchModel_solar_Push/)

## Installation

Currently, this model depends on 3 other repos that need to be installed for this repo to run: [gpkit](https://github.com/convexengineering/gpkit), [gplibrary](https://github.com/convexengineering/gplibrary) and [gassolar](https://github.com/convexengineering/gassolar).  All dependent repos are unit tested on continuous intregration testing software.  A summary of all 3 repos is included below: 

### GPkit 

GPkit is a python package for defining and manuipulating geometric programming models.  More information, including installation instructions can be found at [gpkit.readthedocs.io](gpkit.readthedocs.io).  GPkit uses the information in the constraints defined in the solar model to create the A-matrix, solve the GP, and parse and return the solution to the user.  

NOTE: GPkit supports off the shelf solvers.  Given the complexity of this model, it is advisable that the user install [Mosek](https://www.mosek.com/downloads/) as their solver.  Mosek's presolve and solver algorithms are more reliable and faster than the alternatives.  Covergence to a solution is not guaranteed with the open-source solver, cvxopt.

### GPLibrary

GP library is a collection of GP and SP compatible models.  The solar model imports various models from this library including Wing, Empennage, Fuselage. 

### Gassolar

The gassolar repository contains simple gas and solar GP models used to evaluate a trade study between gas and solar powered long endurance aircraft.  This solar repository relies on the extensive data on wind speeds and solar flux in the gassolar repository. This may be phased out in later versions.

## Usage

To use this model navigate to `./solar` and run `solar.py`.  This will return a dict, `sol` that contains every variable and it's corresponding optimized value.  The python variable, `M`, is the GP model.  `M.substitutions` will return every fixed variable, and it's value.  Substitutions, or fixed values, can be changed using the command: `M.substitutions.update({'Wpay': 15})`. To resolve with an updated substitutions dict enter `sol = M.solve()`.  The objective function for this model is total weight.  The default latitude is 20.    In most cases, the worst case latitude is the highest latitude.  However, to be conservative the user can define a latitude range by defining an new `Mission` model: `M = Mission(latitude=range(1, 20, 1))` and resolving.  The default mission day is the winter solstice, `day=355`.

