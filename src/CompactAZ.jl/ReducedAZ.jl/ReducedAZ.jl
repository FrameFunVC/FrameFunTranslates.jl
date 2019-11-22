module ReducedAZ

using FrameFun.Platforms, BasisFunctions, FrameFun.ApproximationProblems, FrameFun.FrameFunInterface, LinearAlgebra
import FrameFun.FrameFunInterface: solver

export ReducedAZStyle
"""
    struct ReducedAZStyle <: SolverStyle end

Use the reduced AZ algorithm to solve the discretized approximation problem.
"""
struct ReducedAZStyle <: SolverStyle end

include("ReductionSolvers.jl")
using .ReductionSolvers

solver(style::ReducedAZStyle, ap::ApproximationProblem, A::DictionaryOperator; options...) =
    solver(style, ap, A, AZ_Zt(DictionaryOperatorStyle(), ap; (options)...); options...)

function solver(::ReducedAZStyle, ap::ApproximationProblem, A::AbstractOperator, Zt::AbstractOperator;
        threshold=default_threshold(A), verbose=false, crop_tol=threshold, options...)
    plunge_op = I-A*Zt

    verbose && @info "ReducedAZStyle: use `ReductionSolver` as solver for first AZ step"
    psolver = ReductionSolver(plunge_op*A; verbose=verbose, crop_tol=crop_tol, options...)

    AZSolver(A, Zt, plunge_op, psolver)
end

end
