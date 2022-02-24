# FortranNOAu111

[![][ci-img]][ci-url]

[ci-img]: https://github.com/nqcd/NQCDynamics.jl/actions/workflows/CI.yml/badge.svg
[ci-url]: https://github.com/nqcd/NQCDynamics.jl/actions/workflows/CI.yml

This interfaces the Fortran version of the NO on Au(111) model from Roy, Shenvi and Tully in J. Chem. Phys. 130, 174716 (2009)
with NQCModels.jl.

To initialise the model:
```julia
FortranNOAu111Model("path/to/library", initial_geometry; Ms=number_of_bath_states)
```
* `"path/to/library/"` must point to the Fortran shared library containing the relevant functions.
* `initial_geometry` is used to initialise the neighbours of the metal atoms.
* `Ms` determines how many bath states are added.
