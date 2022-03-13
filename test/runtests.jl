using NQCModels
using FortranNOAu111
using Test
using DelimitedFiles
using Unitful, UnitfulAtomic
using FiniteDiff
using LinearAlgebra

@time @testset "find_layer_indices" begin
    r = permutedims(readdlm("$(@__DIR__)/surface_Au111.dat"; skipstart=4))
    no = [
        0.0 2.0;
        0.0 0.0;
        2.0 2.0;
    ]
    r = austrip.(hcat(no, r) .* u"Å")
    @test length(FortranNOAu111.find_layer_indices(r, 5)) == 530
    @test length(FortranNOAu111.find_layer_indices(r, 4)) == 528
    @test length(FortranNOAu111.find_layer_indices(r, 3)) == div(528, 4) * 3
    @test length(FortranNOAu111.find_layer_indices(r, 2)) == div(528, 2)
    @test length(FortranNOAu111.find_layer_indices(r, 1)) == div(528, 4)
end

@time @testset "Numerical comparison to reference" begin

    r = permutedims(readdlm("$(@__DIR__)/surface_Au111.dat"; skipstart=4))
    no = ustrip.(uconvert.(u"Å", [
        0.0 0.0;
        0.0 0.0;
        1.5e-10 2.65e-10;
    ] .* u"m"))
    r = austrip.(hcat(no, r) .* u"Å")
    model = FortranNOAu111Model("$(@__DIR__)/../lib/tullynoau111", r; Ms=10)
    FortranNOAu111.evaluate_energy_force_func!(model, r)

    @testset "Hamiltonian" begin
        reference = readdlm(joinpath(@__DIR__, "reference_data", "hamiltonian.txt"); skipstart=3)
        V = potential(model, r)
        @test reference ≈ V rtol=1e-6
    end

    @testset "Neighbours" begin
        reference = readdlm(joinpath(@__DIR__, "reference_data", "neighbours.txt"), Int; skipstart=0)
        @test reference == model.nn
    end

    @testset "Hp" begin
        reference = [6.4635724443288317E-019, 1.3379130168294729E-018, -4.3176728318301907E-019]
        @test model.Hp ≈ reference
    end

    @testset "dHp1" begin

        reference = permutedims(readdlm(joinpath(@__DIR__, "reference_data", "dHp1.txt"); skipstart=0))
        F = FortranNOAu111.get_neutral_force(model)
        @test reference ≈ F rtol=1e-6
    end

    @time @testset "state_independent_force" begin
        F = zero(r)
        model = FortranNOAu111Model("$(@__DIR__)/../lib/tullynoau111", r; Ms=10, freeze_layers=1)
        NQCModels.state_independent_derivative!(model, F, r)
        LinearAlgebra.lmul!(-1, F)
        reference = permutedims(readdlm(joinpath(@__DIR__, "reference_data", "state_independent_force.txt"); skipstart=0))
        @test reference ≈ F rtol=1e-6
    end

    @time @testset "complete_force" begin
        F = zero(r)
        model = FortranNOAu111Model("$(@__DIR__)/../lib/tullynoau111", r; Ms=10, freeze_layers=1)
        NQCModels.state_independent_derivative!(model, F, r)
        LinearAlgebra.lmul!(-1, F)

        V = potential(model, r)
        eig = eigen(V)
        U = eig.vectors
        D = derivative(model, r)

        state = 1:5
        for I in eachindex(F)
            for k in state
                adiabatic_derivative = U' * D[I] * U
                F[I] -= adiabatic_derivative[k, k]
            end
        end

        reference = permutedims(readdlm(joinpath(@__DIR__, "reference_data", "force.txt"); skipstart=2))
        @test reference ≈ F rtol=1e-6
    end

end

@time @testset "Finite difference" begin

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
    model = FortranNOAu111Model("$(@__DIR__)/../lib/tullynoau111", r; Ms=10)

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
end
