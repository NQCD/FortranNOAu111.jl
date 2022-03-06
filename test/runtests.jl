using NQCModels
using FortranNOAu111
using Test
using DelimitedFiles
using Unitful, UnitfulAtomic
using FiniteDiff
using LinearAlgebra

function finite_difference_gradient(model, R)
    f(x, j, i) = potential(model, x)[j,i]
    grad = [Hermitian(zeros(nstates(model), nstates(model))) for _ in CartesianIndices(R)]
    for i=1:nstates(model)
        for j=1:nstates(model)
            gradient = FiniteDiff.finite_difference_gradient(x->f(x,j,i), R)
            for k in eachindex(R)
                grad[k].data[j,i] = gradient[k]
            end
        end
    end
    grad
end

r = permutedims(readdlm("$(@__DIR__)/surface_Au111.dat"; skipstart=4))
no = [
    0.0 2.0;
    0.0 0.0;
    2.0 2.0;
]
r = austrip.(hcat(no, r) .* u"Å")

model = FortranNOAu111Model("$(@__DIR__)/../lib/tullynoau111", r; Ms=6)

@test length(FortranNOAu111.find_layer_indices(r, 5)) == 530
@test length(FortranNOAu111.find_layer_indices(r, 4)) == 528
@test length(FortranNOAu111.find_layer_indices(r, 3)) == div(528, 4) * 3
@test length(FortranNOAu111.find_layer_indices(r, 2)) == div(528, 2)
@test length(FortranNOAu111.find_layer_indices(r, 1)) == div(528, 4)

@time @testset "Finite difference gradient for individual elements" begin
    function V_neutral(x)
        FortranNOAu111.evaluate_energy_force_func!(model, x)
        return FortranNOAu111.get_neutral_element(model)
    end

    function V_ion(x)
        FortranNOAu111.evaluate_energy_force_func!(model, x)
        return FortranNOAu111.get_ion_element(model)
    end

    function V_coupling(x)
        FortranNOAu111.evaluate_energy_force_func!(model, x)
        return FortranNOAu111.get_coupling_element(model)
    end

    F1 = FiniteDiff.finite_difference_gradient(V_neutral, r)
    F2 = FiniteDiff.finite_difference_gradient(V_ion, r)
    F3 = FiniteDiff.finite_difference_gradient(V_coupling, r)

    FortranNOAu111.evaluate_energy_force_func!(model, r)
    analytic_neutral = FortranNOAu111.get_neutral_force(model)
    analytic_ion = FortranNOAu111.get_ion_force(model)
    analytic_coupling = FortranNOAu111.get_coupling_force(model)
    @test F1 ≈ analytic_neutral 
    @test F2 ≈ analytic_ion
    @test F3 ≈ analytic_coupling
end

@time @testset "Finite difference state independent" begin
    V(x) = NQCModels.state_independent_potential(model, x)
    F = FiniteDiff.finite_difference_gradient(V, r)
    D = zero(r)
    NQCModels.state_independent_derivative!(model, D, r)
    @test D ≈ F
end

@time @testset "Finite difference state dependent" begin
    finite_diff = finite_difference_gradient(model, r)
    analytic = derivative(model, r)
    @test finite_diff ≈ analytic
end
