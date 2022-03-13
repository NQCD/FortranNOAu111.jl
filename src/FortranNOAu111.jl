module FortranNOAu111

export FortranNOAu111Model

using NQCModels: NQCModels
using Unitful, UnitfulAtomic
using Libdl: Libdl
using LinearAlgebra: Hermitian, eigen, diagind
using FastGaussQuadrature: gausslegendre

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
    DeltaE::Float64
    nelectrons::Int
    freeze::Vector{Int}

    x_gauss::Vector{Float64}
    w_gauss::Vector{Float64}
end

function FortranNOAu111Model(library_path, r; Ms, VAuFlag=1, freeze=Int[], freeze_layers=0)

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

    U = [
       -1/sqrt(2) 1/sqrt(6) -1/sqrt(3)
        0.0      -2/sqrt(6) -1/sqrt(3)
        1/sqrt(2) 1/sqrt(6) -1/sqrt(3)
    ]
    a = sqrt(2) * dnn
    r0 = transpose(U) * reshape([0,1,1,0,1,-1,1,0,1,-1,0,1,1,1,0,1,-1,0,0,-1,-1,0,-1,1,-1,0,-1,1,0,-1,-1,-1,0,-1,1,0], (3,12)) ./ 2 .* a

    mass = zeros(N+2)
    mass[1] = ustrip(uconvert(u"kg", 14.00307440u"u"))
    mass[2] = ustrip(uconvert(u"kg", 15.99491502u"u"))
    mass[3:N+2] .= ustrip(uconvert(u"kg", 196.966548u"u"))

    DeltaE = austrip(7u"eV")
    nelectrons = fld(Ms, 2)

    x_gauss, w_gauss = gausslegendre(nelectrons)

    if freeze_layers != 0
        freeze = find_layer_indices(r, freeze_layers)
    end

    FortranNOAu111Model(energy_force_func, x, nn, N, dnn, aPBC, Hp, dHp, VAuFlag, r0, mass, Ms,
        DeltaE, nelectrons, freeze, x_gauss, w_gauss
    )
end

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

function find_layer_indices(r, layers)
    freeze = Int[]
    layers == 0 && return freeze

    z_coordinates = @view r[3,:]
    permutation = sortperm(z_coordinates)
    ordered_z = z_coordinates[permutation]
    current_z = ordered_z[begin]
    current_layer = 1

    for (i, z) in enumerate(ordered_z)
        if !isapprox(z, current_z)
            current_layer += 1
            current_layer > layers && break
            current_z = z
        end
        push!(freeze, permutation[i])
    end
    
    return freeze
end

NQCModels.ndofs(::FortranNOAu111Model) = 3
NQCModels.nstates(model::FortranNOAu111Model) = model.Ms+1

convert_energy(x) = austrip(x * u"J")
convert_force(x) = austrip(x * u"J/m")

get_neutral_element(model::FortranNOAu111Model) = convert_energy(model.Hp[1])
get_ion_element(model::FortranNOAu111Model) = convert_energy(model.Hp[2])
get_coupling_element(model::FortranNOAu111Model) = convert_energy(model.Hp[3])

get_neutral_force(model::FortranNOAu111Model) = @views convert_force.(reshape(model.dHp[1,:], 3, :))
get_ion_force(model::FortranNOAu111Model) = @views convert_force.(reshape(model.dHp[2,:], 3, :))
get_coupling_force(model::FortranNOAu111Model) = @views convert_force.(reshape(model.dHp[3,:], 3, :))

function set_coordinates!(model::FortranNOAu111Model, r)
    model.x .= ustrip.(auconvert.(u"m", r))
end

function NQCModels.state_independent_potential(model::FortranNOAu111Model, r::AbstractMatrix)
    evaluate_energy_force_func!(model, r)
    return get_neutral_element(model)
end

function NQCModels.state_independent_derivative!(model::FortranNOAu111Model, D::AbstractMatrix, r::AbstractMatrix)
    evaluate_energy_force_func!(model, r)

    copyto!(D, get_neutral_force(model))
    freeze_atoms!(D, model.freeze)
    return D
end

function NQCModels.potential!(model::FortranNOAu111Model, V::Hermitian, r::AbstractMatrix)
    evaluate_energy_force_func!(model, r)

    (;DeltaE, x_gauss, w_gauss) = model

    bath = @view V[diagind(V)[2:end]]
    set_bath_energies!(bath, x_gauss, DeltaE)

    E_coup = get_coupling_element(model)
    couplings = @view V.data[2:end,1]
    set_coupling_elements!(couplings, w_gauss, DeltaE, E_coup)
    couplings = @view V.data[1,2:end]
    set_coupling_elements!(couplings, w_gauss, DeltaE, E_coup)

    V[1,1] = get_ion_element(model) - get_neutral_element(model)

    return V
end

function NQCModels.derivative!(model::FortranNOAu111Model, D::AbstractMatrix{<:Hermitian}, r::AbstractMatrix)
    evaluate_energy_force_func!(model, r)

    (;w_gauss, DeltaE) = model

    F_coup = get_coupling_force(model)
    for I in eachindex(D)
        coupling = @view D[I].data[2:end,1]
        set_coupling_elements!(coupling, w_gauss, DeltaE, F_coup[I])
        coupling = @view D[I].data[1,2:end]
        set_coupling_elements!(coupling, w_gauss, DeltaE, F_coup[I])
        D[I][1,1] = get_ion_force(model)[I] - get_neutral_force(model)[I]
    end
    freeze_atoms!(D, model.freeze)

    return D
end

function freeze_atoms!(D::AbstractMatrix, freeze_indices)
    @views for i in freeze_indices
        fill!(D[:,i], zero(eltype(D)))
    end
end

function freeze_atoms!(D::AbstractMatrix{<:Hermitian}, freeze_indices)
    for i in freeze_indices
        for j in axes(D, 1)
            fill!(D[j,i], zero(eltype(D[j,i])))
        end
    end
end

function set_bath_energies!(bath, x_gauss, DeltaE)
    for i in eachindex(x_gauss)
        bath[i] = DeltaE * (-1 + x_gauss[i]) / 4
    end
    n = length(x_gauss)
    for i in eachindex(x_gauss)
        bath[i+n] = DeltaE * (1 + x_gauss[i]) / 4
    end
end

function set_coupling_elements!(coupling, w_gauss, DeltaE, E_coup)
    nelectrons = length(w_gauss)
    for i in eachindex(w_gauss)
        coupling[i] = sqrt(DeltaE * w_gauss[i]) / 2 * E_coup / sqrt(DeltaE)
    end
    for i in eachindex(w_gauss)
        coupling[i+nelectrons] = coupling[i]
    end
end

function evaluate_energy_force_func!(model::FortranNOAu111Model, r::AbstractMatrix)
    set_coordinates!(model, r)
    (;N, x, nn, r0, aPBC, mass, VAuFlag, Hp, dHp) = model
    ccall(
        model.energy_force_func, Cvoid,
        (Ref{Cint}, Ref{Cdouble}, Ref{Cint}, Ref{Cdouble}, Ref{Cdouble},
        Ref{Cdouble}, Ref{Cint}, Ref{Cdouble}, Ref{Cdouble}),
        N, x, nn, r0, aPBC, mass, VAuFlag, Hp, dHp
    )
end

end
