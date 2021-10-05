module SparseAZ

using FrameFun, BasisFunctions
using LinearAlgebra, ..CompactFrameFunExtension
import FrameFun: solver

using FrameFun: ApproximationProblem, default_threshold

solver(style::SparseAZStyle, ap::ApproximationProblem, A::DictionaryOperator; samplingstyle=SamplingStyle(ap), options...) =
    solver(samplingstyle, style, ap, A, AZ_Zt(ap; options)...; options...)

function solver(samplingstyle::SamplingStyle, solverstyle::SparseAZStyle, ap::ApproximationProblem, A::AbstractOperator, Zt::AbstractOperator;
        threshold=default_threshold(A), nz_tol=threshold, verbose=false, options...)
    plunge_op = I-A*Zt
    dict1 = src(A)
    dict2 = dest(Zt)
    os_grid = grid(dest(A))

    psolver = sparseAZ_AAZAreductionsolver(samplingstyle, ap; threshold=threshold, nz_tol=nz_tol, solverstyle=solverstyle, os_grid=os_grid, verbose=verbose, options...)
    AZSolver(A, Zt, plunge_op, psolver)
end
end
