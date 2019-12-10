module ReducedAZ

using FrameFun.Platforms, BasisFunctions, FrameFun.ApproximationProblems, FrameFun.FrameFunInterface, LinearAlgebra, ..CompactFrameFunExtension
import FrameFun.FrameFunInterface: solver

solver(style::ReducedAZStyle, ap::ApproximationProblem, A::DictionaryOperator; samplingstyle = SamplingStyle(ap), options...) =
    solver(samplingstyle, style, ap, A, AZ_Zt(DictionaryOperatorStyle(), ap; (options)...); options...)

function solver(samplingstyle::OversamplingStyle, solverstyle::ReducedAZStyle, ap::ApproximationProblem, A::AbstractOperator, Zt::AbstractOperator;
        threshold=default_threshold(A), nz_tol=threshold, verbose=false, options...)
    plunge_op = I-A*Zt
    os_grid = grid(dest(A))
    verbose && @info "ReducedAZStyle: use `ReductionSolver` as solver for first AZ step"
    psolver = reducedAZ_AAZAreductionsolver(samplingstyle, ap; threshold=threshold, nz_tol=nz_tol, solverstyle=solverstyle, L=samplingparameter(ap), os_grid=os_grid, verbose=verbose, options...)

    AZSolver(A, Zt, plunge_op, psolver)
end

end
