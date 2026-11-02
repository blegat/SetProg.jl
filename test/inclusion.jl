using Test
using LinearAlgebra
import GLPK
using SetProg, Polyhedra

function _test_preprocess(h, nhyperplanes, nhalfspaces)
    for solver in (GLPK.Optimizer, nothing)
        lib = Polyhedra.DefaultLibrary{Float64}(solver)
        domain = SetProg._preprocess_domain(polyhedron(h, lib))
        @test Polyhedra.nhyperplanes(domain) == nhyperplanes
        @test Polyhedra.nhalfspaces(domain) == nhalfspaces
        B = reduce(vcat, [hp.a' for hp in hyperplanes(domain)], init = zeros(0, fulldim(h)))
        @test B * B' ≈ I
        for hs in halfspaces(domain)
            @test norm(hs.a) ≈ 1
            @test norm(B * hs.a) < 1e-8
        end
    end
end

@testset "Preprocess domain" begin
    @testset "Opposite halfspaces with noise" begin
        h = HalfSpace([-0.5], 0.0) ∩ HalfSpace([5.551115123125783e-17], 0.0) ∩
            HalfSpace([0.5000000000000002], 0.0)
        _test_preprocess(h, 1, 0)
    end
    @testset "Duplicate halfspaces" begin
        h = HalfSpace([0.5], 0.0) ∩ HalfSpace([1.5], 0.0) ∩ HalfSpace([0.0], 0.0)
        _test_preprocess(h, 0, 1)
    end
    @testset "Halfspace implied by hyperplane" begin
        h = HyperPlane([1.5], 0.0) ∩ HalfSpace([-0.5], 0.0)
        _test_preprocess(h, 1, 0)
        h = HyperPlane([1.0, 1e-17], 0.0) ∩ HalfSpace([1.0, 1.0], 0.0) ∩
            HalfSpace([2.0, 2.0], 0.0) ∩ HalfSpace([-3.0, 0.0], 0.0)
        _test_preprocess(h, 1, 1)
    end
    @testset "Full-dimensional cone" begin
        h = HalfSpace([-1.0, 0.0], 0.0) ∩ HalfSpace([0.0, -2.0], 0.0)
        _test_preprocess(h, 0, 2)
    end
    @testset "Only cones are supported" begin
        lib = Polyhedra.DefaultLibrary{Float64}(GLPK.Optimizer)
        @test_throws ErrorException("only cones are supported") SetProg._preprocess_domain(
            polyhedron(hrep([HalfSpace([1.0], 1.0)]), lib),
        )
    end
end
