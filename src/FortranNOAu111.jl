module FortranNOAu111

export FortranNOAu111Model

using NQCModels: NQCModels
using Unitful, UnitfulAtomic
using Libdl: Libdl
using LinearAlgebra: Hermitian, eigen

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
    dA::Vector{Float64}
    vA::Matrix{Float64}
    DeltaE::Float64
    nelectrons::Int
    freeze::Vector{Int}
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

    a = sqrt(2) * dnn
    U = reshape([-1/sqrt(2),0,1/sqrt(2), 1/sqrt(2),-2/sqrt(6),1/sqrt(6), -1/sqrt(3),-1/sqrt(3),-1/sqrt(3) ], (3,3))
    r0 = transpose(U) * reshape([0,1,1,0,1,-1,1,0,1,-1,0,1,1,1,0,1,-1,0,0,-1,-1,0,-1,1,-1,0,-1,1,0,-1,-1,-1,0,-1,1,0], (3,12)) ./ 2 .* a

    mass = zeros(N+2)
    mass[1] = ustrip(uconvert(u"kg", 14.00307440u"u"))
    mass[2] = ustrip(uconvert(u"kg", 15.99491502u"u"))
    mass[3:N+2] .= ustrip(uconvert(u"kg", 196.966548u"u"))

    DeltaE = austrip(7u"eV")
    nelectrons = fld(Ms, 2)
    eigs = eigen(get_Agauss(nelectrons))
    dA = eigs.values
    vA = eigs.vectors

    if freeze_layers != 0
        freeze = find_layer_indices(r, freeze_layers)
    end

    FortranNOAu111Model(energy_force_func, x, nn, N, dnn, aPBC, Hp, dHp, VAuFlag, r0, mass, Ms,
        dA, vA, DeltaE, nelectrons, freeze
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
    D .= get_neutral_force(model)
    D[:,model.freeze] .= zero(eltype(r))
    return D
end

function NQCModels.potential!(model::FortranNOAu111Model, V::Hermitian, r::AbstractMatrix)
    evaluate_energy_force_func!(model, r)

    burkey_cantrell_continuum!(model, V)
    V[1,1] = get_ion_element(model) - get_neutral_element(model)

    return V
end

function NQCModels.derivative!(model::FortranNOAu111Model, D::AbstractMatrix{<:Hermitian}, r::AbstractMatrix)
    evaluate_energy_force_func!(model, r)

    for I in eachindex(D)
        burkey_cantrell_continuum_derivative!(model, D, I)
        D[I][1,1] = get_ion_force(model)[I] - get_neutral_force(model)[I]
    end
    for i in model.freeze
        for j in axes(D, 1)
            fill!(D[j,i], zero(eltype(r)))
        end
    end

    return D
end

function get_Agauss(nelectrons)
    mtemp = zeros(nelectrons,nelectrons)
    for i in range(1,nelectrons-1)
        mtemp[i,i+1] = i/sqrt((2i+1) * (2i-1))
    end
    mtemp = Hermitian(mtemp)
    return mtemp
end

function burkey_cantrell_continuum!(model::FortranNOAu111Model, V::Hermitian)

    (;nelectrons, dA, vA, DeltaE) = model
    for i in range(2, nelectrons+1) # Set diagonal entries
        V.data[i,i] =  (DeltaE/4.0)*dA[i-1] - DeltaE/4.0
        V.data[i+nelectrons,i+nelectrons] =  (DeltaE/4.0)*dA[i-1] + DeltaE/4.0
    end

    E_coup = get_coupling_element(model)
    vtemp = E_coup*sqrt.((DeltaE*2*vA[1,:].*vA[1,:])/4)./sqrt(DeltaE)

    V.data[1,2:nelectrons+1] .= vtemp
    V.data[2:nelectrons+1,1] .= vtemp
    V.data[1,2+nelectrons:end] .= vtemp
    V.data[2+nelectrons:end,1] .= vtemp
end

function burkey_cantrell_continuum_derivative!(model::FortranNOAu111Model, derivative::AbstractMatrix, I)

    (;nelectrons, vA, DeltaE) = model
    FF = get_coupling_force(model)[I]
    # The expression here is without the division by the square root,
    # bcs the square-root comes into the energy expression from the integration
    # of the Schroedinger equation. The forces are not subject to Schroedinger.
    # James: I'm using the sqrt to be consistent with the energy. Maybe we don't need is here or in the energy?
    dvtemp =FF*sqrt.((DeltaE*2*vA[1,:].*vA[1,:])/4.0)./sqrt(DeltaE)
    # dvtemp = FF*sqrt.((DeltaE*2*vA[1,:].*vA[1,:])/4.0)

    derivative[I].data[1,2:nelectrons+1] .= dvtemp
    derivative[I].data[2:nelectrons+1,1] .= dvtemp
    derivative[I].data[1,2+nelectrons:end] .= dvtemp
    derivative[I].data[2+nelectrons:end,1] .= dvtemp
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
