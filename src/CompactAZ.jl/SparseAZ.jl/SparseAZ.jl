# make number of nonzerocolumns smaller
# check correctness of difference_indices
# write nonzero_cols if domain is a grid

module SparseAZ



using FrameFun.Platforms
using FrameFun.BasisFunctions, FrameFun.ApproximationProblems, FrameFun.FrameFunInterface, LinearAlgebra
import FrameFun.FrameFunInterface: solver
using ...FrameFunTranslates.SPQR_Solvers: SPQR_solver, sparseQR_solver

export SparseAZStyle
"""
    struct SparseAZStyle <: SolverStyle end

Use the sparse AZ algorithm to solve the discretized approximation problem.
"""
struct SparseAZStyle <: SolverStyle end

include("sparseAAZA.jl")

solver(style::SparseAZStyle, ap::ApproximationProblem, A::DictionaryOperator; options...) =
    solver(style, ap, A, AZ_Zt(DictionaryOperatorStyle(), ap; (options)...); options...)

function solver(::SparseAZStyle, ap::ApproximationProblem, A::AbstractOperator, Zt::AbstractOperator;
        verbose=false, REG=sparseQR_solver, threshold=default_threshold(A), options...)
    plunge_op = I-A*Zt
    dict1 = src(A)
    dict2 = dest(Zt)
    g = grid(dest(A))

    colix = nonzero_cols(dict1,g)
    verbose && @info "SparseAZStyle: create sparse A-AZ^*A"
    AAZA = sparseAAZAmatrix(dict1,dict2,g; atol=threshold)
    verbose && @info "SparseAZStyle: A-AZ^*A has size $(size(AAZA)) and $(nnz(AAZA)) nonzero elements ($(100nnz(AAZA)/prod(size(AAZA)))% fill)"
    verbose && @info "SparseAZStyle: use $(REG) as solver for first AZ step"
    psolver = IndexExtensionOperator(dict1,colix)*REG(ArrayOperator(sparseAAZAmatrix(dict1,dict2,g; atol=threshold), dict1[colix], dest(A));
        verbose=verbose, threshold=threshold, options...)
    verbose && @info "SparseAZStyle: $(REG) created"
    AZSolver(A, Zt, plunge_op, psolver)
end

end
