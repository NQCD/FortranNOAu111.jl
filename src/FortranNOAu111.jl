module FortranNOAu111

export FortranNOAu111Model

using NQCModels: NQCModels
using Unitful, UnitfulAtomic
using LinearAlgebra: Hermitian, eigen, diagind
using StaticArrays: SMatrix

const libtullynoau111 = "tullynoau111"

struct FortranNOAu111Model <: NQCModels.DiabaticModels.DiabaticModel
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
    freeze::Vector{Int}
    mobile_atoms::Vector{Int}

    tmp_neutral_force::Matrix{Float64}
    tmp_ion_force::Matrix{Float64}
    tmp_coupling_force::Matrix{Float64}

    cache_positions::Matrix{Float64}
end

function FortranNOAu111Model(r; VAuFlag=1, freeze=Int[], freeze_layers=0)

    x = copy(r)
    x = ustrip.(auconvert.(u"m", x))
    N = size(r, 2) - 2
    dnn = 2.95217081u"Å"
    aPBC = [32.47387893u"Å", 30.67985848u"Å", sqrt(6) * dnn * 100]
    dnn = ustrip(uconvert(u"m", dnn))
    aPBC = ustrip.(uconvert.(u"m", aPBC))
    nn = get_nearest_neighbours(N, x, dnn * 1.01, aPBC)
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

    if freeze_layers != 0
        freeze = find_layer_indices(r, freeze_layers)
    end
    mobile_atoms = setdiff(axes(r, 2), freeze)

    tmp_neutral_force = zero(r)
    tmp_ion_force = zero(r)
    tmp_coupling_force = zero(r)

    cache_positions = zero(r)

    FortranNOAu111Model(x, nn, N, dnn, aPBC, Hp, dHp, VAuFlag, r0, mass,
        freeze, mobile_atoms,
        tmp_neutral_force, tmp_ion_force, tmp_coupling_force,
        cache_positions
    )
end

NQCModels.mobileatoms(model::FortranNOAu111Model, ::Int) = model.mobile_atoms
NQCModels.mobileatoms(model::FortranNOAu111Model) = model.mobile_atoms
NQCModels.nelectrons(model::FortranNOAu111Model) = model.nelectrons

function get_nearest_neighbours(N, r, dnn, aPBC)
    nn = zeros(Cint, N, 12)
    ccall(
        (:get_nn, libtullynoau111), Cvoid,
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
NQCModels.nstates(model::FortranNOAu111Model) = 2

convert_energy(x) = austrip(x * u"J")
convert_force(x) = austrip(x * u"J/m")

get_neutral_element(model::FortranNOAu111Model) = convert_energy(model.Hp[1])
get_ion_element(model::FortranNOAu111Model) = convert_energy(model.Hp[2])
get_coupling_element(model::FortranNOAu111Model) = convert_energy(model.Hp[3])

function get_neutral_force!(model::FortranNOAu111Model)
    source = @view model.dHp[1,:]
    output = model.tmp_neutral_force
    convert_force_in_buffer!(source, output)
    return output
end

function get_ion_force!(model::FortranNOAu111Model)
    source = @view model.dHp[2,:]
    output = model.tmp_ion_force
    convert_force_in_buffer!(source, output)
    return output
end

function get_coupling_force!(model::FortranNOAu111Model)
    source = @view model.dHp[3,:]
    output = model.tmp_coupling_force
    convert_force_in_buffer!(source, output)
    return output
end

function convert_force_in_buffer!(source, output)
    for I in eachindex(output, source)
        output[I] = convert_force(source[I])
    end
end

function set_coordinates!(model::FortranNOAu111Model, r)
    for I in eachindex(model.x, r)
        model.x[I] = ustrip(auconvert(u"m", r[I]))
    end
end

function NQCModels.potential(model::FortranNOAu111Model, r::AbstractMatrix)
    evaluate_energy_force_func!(model, r)

    v11 = get_neutral_element(model)
    v12 = get_coupling_element(model)
    v22 = get_ion_element(model)

    return Hermitian(SMatrix{2,2}(v11, v12, v12, v22))

    return V
end

function NQCModels.derivative!(model::FortranNOAu111Model, D::AbstractMatrix{<:Hermitian}, r::AbstractMatrix)
    evaluate_energy_force_func!(model, r)

    F_neutral = get_neutral_force!(model)
    F_ion = get_ion_force!(model)
    F_coup = get_coupling_force!(model)

    @inbounds for i in NQCModels.mobileatoms(model)
        for j in axes(r, 1)
            d11 = F_neutral[j,i]
            d12 = F_coup[j,i]
            d22 = F_ion[j,i]
            D[j,i] = Hermitian(SMatrix{2,2}(d11, d12, d12, d22))
        end
    end

    return D
end
function evaluate_energy_force_func!(model::FortranNOAu111Model, r::AbstractMatrix)
    if r != model.cache_positions
        set_coordinates!(model, r)
        (;N, x, nn, r0, aPBC, mass, VAuFlag, Hp, dHp) = model
        ccall(
            (:get_energies_and_forces, libtullynoau111), Cvoid,
            (Ref{Cint}, Ref{Cdouble}, Ref{Cint}, Ref{Cdouble}, Ref{Cdouble},
            Ref{Cdouble}, Ref{Cint}, Ref{Cdouble}, Ref{Cdouble}),
            N, x, nn, r0, aPBC, mass, VAuFlag, Hp, dHp
        )
        copy!(model.cache_positions, r)
    end
end

end
