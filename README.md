# FortranNOAu111

[![][ci-img]][ci-url]

[ci-img]: https://github.com/NQCD/FortranNOAu111.jl/actions/workflows/CI.yml/badge.svg
[ci-url]: https://github.com/NQCD/FortranNOAu111.jl/actions/workflows/CI.yml

This interfaces the Fortran version of the NO on Au(111) model from Roy, Shenvi and Tully in J. Chem. Phys. 130, 174716 (2009)
with NQCModels.jl.

To initialise the model,
```julia
using Libdl: DL_LOAD_PATH
push!(DL_LOAD_PATH, "path/to/library/directory")
FortranNOAu111Model(initial_geometry; Ms=number_of_bath_states)
```
* `"path/to/library/directory"` must point to the Fortran shared library that can be compiled in the `lib` directory.
* `initial_geometry` is used to initialise the neighbours of the metal atoms.
* `Ms` determines how many bath states are added.
