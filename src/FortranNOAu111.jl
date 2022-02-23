module FortranNOAu111

export FortranNOAu111Model

using NQCModels: NQCModels
using Unitful, UnitfulAtomic
using Libdl: Libdl
using LinearAlgebra: Hermitian

struct FortranNOAu111Model <: NQCModels.DiabaticModels.LargeDiabaticModel
    energy_force_func::Ptr{Nothing}
    x::Matrix{Cdouble}
    nn::Matrix{Cint}
    N::Cint
    dnn::Cdouble
    aPBC::Vector{Cdouble}
    Hp::Vector{Cdouble}
    dHp::Matrix{Cdouble}
    VAuFlag::Cint
    r0::Matrix{Cdouble}
    mass::Vector{Cdouble}
    Ms::Cint
end

function FortranNOAu111Model(library_path, r; Ms, VAuFlag=1)

    lib = Libdl.dlopen(library_path)
    energy_force_func = Libdl.dlsym(lib, :get_energies_and_forces)

    x = copy(r)
    x = ustrip.(auconvert.(u"m", x))
    N = size(r, 2) - 2
    dnn = 2.95217081u"Å"
    aPBC = [32.47387893u"Å", 30.67985848u"Å", sqrt(6) * dnn * 100]
    dnn = ustrip(uconvert(u"m", dnn))
    aPBC = ustrip.(uconvert.(u"m", aPBC))
    nn = get_nearest_neighbours(lib, N, x, dnn * 1.01, aPBC)
    Hp = zeros(3)
    dHp = zeros(3,3*(N+2))

    a = sqrt(2) * dnn
    U = reshape([-1/sqrt(2),0,1/sqrt(2), 1/sqrt(2),-2/sqrt(6),1/sqrt(6), -1/sqrt(3),-1/sqrt(3),-1/sqrt(3) ], (3,3))
    r0 = transpose(U) * reshape([0,1,1,0,1,-1,1,0,1,-1,0,1,1,1,0,1,-1,0,0,-1,-1,0,-1,1,-1,0,-1,1,0,-1,-1,-1,0,-1,1,0], (3,12)) ./ 2a

    mass = zeros(N+2)
    mass[1] = ustrip(uconvert(u"kg", 14.00307440u"u"))
    mass[2] = ustrip(uconvert(u"kg", 15.99491502u"u"))
    mass[3:N+2] .= ustrip(uconvert(u"kg", 196.966548u"u"))

    FortranNOAu111Model(energy_force_func, x, nn, N, dnn, aPBC, Hp, dHp, VAuFlag, r0, mass, Ms)
end

NQCModels.ndofs(::FortranNOAu111Model) = 3
NQCModels.nstates(model::FortranNOAu111Model) = model.Ms+1

function get_nearest_neighbours(lib, N, r, dnn, aPBC)
    nn_func = Libdl.dlsym(lib, :get_nn)
    nn = zeros(Cint, N, 12)
    ccall(
        nn_func, Cvoid,
        (Ref{Cint}, Ptr{Cdouble}, Ref{Cdouble}, Ptr{Cdouble}, Ptr{Cint}),
        N, r, dnn * 1.01, aPBC, nn
    )
    return nn
end

function set_coordinates!(model::FortranNOAu111Model, r)
    model.x .= ustrip.(auconvert.(u"m", r))
end

function NQCModels.potential!(model::FortranNOAu111Model, V::Hermitian, r::AbstractMatrix)
    evaluate_energy_force_func!(model, r)
end

function NQCModels.derivative!(model::FortranNOAu111Model, D::AbstractMatrix{<:Hermitian}, r::AbstractMatrix)
    evaluate_energy_force_func!(model, r)
end

function evaluate_energy_force_func!(model::FortranNOAu111Model, r::AbstractMatrix)
    set_coordinates!(model, r)
    (;N, x, nn, r0, aPBC, mass, VAuFlag, Hp, dHp) = model
    ccall(
        model.energy_force_func, Cvoid,
        (Ref{Cint},
        Ptr{Cdouble},
        Ptr{Cint},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ref{Cint},
        Ptr{Cdouble}, Ptr{Cdouble}),
        N, x, nn, r0, aPBC, mass, VAuFlag, Hp, dHp
    )
end

end
