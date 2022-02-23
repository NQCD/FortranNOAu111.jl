using NQCModels
using FortranNOAu111
using Test
using DelimitedFiles
using Unitful, UnitfulAtomic
using FiniteDiff

r = permutedims(readdlm("$(@__DIR__)/surface_Au111.dat"; skipstart=4))
no = [
    0.0 2.0;
    0.0 0.0;
    2.0 2.0;
]

model = FortranNOAu111Model("$(@__DIR__)/tullynoau111.dylib", r; Ms=5, VAuFlag=0)
potential(model, r)

function V_neutral(x)
    potential(model, x)
    e = model.Hp[1]
    return austrip(e * u"J")
end

F = FiniteDiff.finite_difference_gradient(V_neutral, r)
@info "Finite difference"
display(F)

potential(model, r)
F2 = austrip.(reshape(model.dHp[1,:], 3, :) .* u"J/m")

F ./ F2
