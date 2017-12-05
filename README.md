# Solar GP Model

This repo contains a solar aircraft sizing model that leverages geometric programming (GP).  Provisions are added to allow signomial programming (SP) for higher fidelity results.  All changes are unit tested on continuous integration software.  Code is python dependent. 

[![Build Status](https://acdl.mit.edu/csi/buildStatus/icon?job=gpkit_ResearchModel_solar_Push)](https://acdl.mit.edu/csi/job/gpkit_ResearchModel_solar_Push/)

## Installation

Currently, this model depends on 3 other repos that need to be installed for this repo to run: [gpkit](https://github.com/convexengineering/gpkit), [gplibrary](https://github.com/convexengineering/gplibrary) and [gassolar](https://github.com/convexengineering/gassolar).  All dependent repos are unit tested on continuous intregration testing software.  A summary of all 3 repos is included below. 

To install GPkit follow the instructions here: [gpkit.readthedocs.installation](https://gpkit.readthedocs.io/en/latest/installation.html)
NOTE: GPkit supports off the shelf solvers.  Given the complexity of this model, it is advisable that the user install [Mosek](https://www.mosek.com/downloads/) as their solver.  Mosek's presolve and solver algorithms are more reliable and faster than the alternatives.  Covergence to a solution is not guaranteed with the open-source solver, cvxopt.

To install gplibrary and gassolar run the following in the command line in your home directory: 

```
pip install git+git://github.com/convexengineering/gplibrary.git
pip install git+git://github.com/convexengineering/gassolar.git
```

Finally, to install this repo run: 

`pip install git+git://github.com/convexengineering/solar.git`

Please report any installation issues in the respective repositories. 

### GPkit 

GPkit is a python package for defining and manuipulating geometric programming models.  More information can be found at [gpkit.readthedocs.io](gpkit.readthedocs.io).  GPkit uses the information in the constraints defined in the solar model to create the A-matrix, solve the GP, and parse and return the solution to the user.  

### GPLibrary

GP library is a collection of GP and SP compatible models.  The solar model imports various models from this library including Wing, Empennage, Fuselage. 

### Gassolar

The gassolar repository contains simple gas and solar GP models used to evaluate a trade study between gas and solar powered long endurance aircraft.  This solar repository relies on the extensive data on wind speeds and solar flux in the gassolar repository. This may be phased out in later versions.

## Usage

### Understanding the Model

The actual model is in `./solar/solar.py`.  Running this file will return a dict, `sol`, that constains every variable and it's corresponding optimized value.  The python variable, `M`, is the GP model and contains all of the submodels used in the optimization. 

Taking a line by line overview of `solar.py`, there are multiple submodels imported from gplibrary and gassolar.  The first created model is `Aircraft`.  This models creates subcomponents such as the wing, empennage, solar cells, many of which are imported from gplibrary. Each submodel, contains a list of constraints specific to that model.  The `Aircraft` model has an aircraft performance model `AircraftPerf`, that is created when a flight segment is created.  `AircraftPerf` sums the drag of each component and has constraints relating the solar flux, motor power draw and battery energy.  The `Mission` model creates multiple flight segments depending on how many latitudes are being evaluated. 

Each model has a list of Variables that are created when the mission is created.  Variables with no number are free variables that are optimized during a solve.  Variables with a number are parameters and are fixed.  These can be changed to evaluate the model at different parameter values. A common way to think about variables in a GP is that every variable needs an upper and lower bound.  A variable `x`, can be upper bounded by minimizing its value in the objective function or including it in a constraint such that it is on the right side of a greater than inequality (i.e. `>= x`).  A lower bound is created by maximizing the value in the objective function or placing it on the left side of a greater than inequality in a constraint (i.e. `x >=`).  Free variables that are not upper or lower bounded by the constraints unique to the model in which they are created are identified in that model's docstring.  These variables should be bounded by constraints in a different model or by the objective function. 

### Creating and Solving 

To create and run the primary GP model run the following in a python notebook: 

```python
from solar.solar import Mission
M = Mission()
sol = M.solve()
```
The equivalent is also done by running `./solar/solar.py`. To view the entire solution enter 

```python 
print sol.table()
```

To view any individual variable enter the variable as an argument into the solution dict. For example, to see the total weight variable `Wtotal`.  

```python
sol(M.solar.Wtotal)
```

### Changing Parameters and Values

To change any fixed variable and resolve, enter 
```python
M.substitutions.update({M.solar.batteries.hbatt: 400})
sol = M.solve()
```
Notice that each variable can be accessed using the dot notation.  Latitude is not a variable that can be updated in this way.  Because of the descrete wind speed and solar flux equations that vary with latitude, different latitudes are evaluated by recreating the model

```python
M = Mission(latitude=[25])
sol = M.solve()
```

The default latitude is 20 degrees.  In most cases, the worst case latitude is the highest latitude. However to be conservative the user can define a latitude range

```python
M = Mission(latitude=range(1, 25, 0))
sol = M.solve()
```

Changing parameters may result in an infeasible solution.  If this happens a `PRIM_INFEAS` error will be thrown when trying to solve. This will happen if one or more of the constraints cannot be satisfied for the given input parameters. Running `M.debug()` may help identify which variables are causing the solution to be infeasible. 

### Changing Configurations

For more experienced users, you can change the configuration by altering the submodels during creation in the `Aircraft` model. For example, to change spar configuration from a `BoxSpar` to a `CapSpar`, change the line in the `Aircraft` model from 

```python
WingGP.sparModel = BoxSpar
```
to
```python
WingGP.sparModel = CapSpar
```

## Reporting Issues

Please report any issues on the github repository by making any issue. In the issue, copy the error message and include a brief description of what you were trying to do.  

