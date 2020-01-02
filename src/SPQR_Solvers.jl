module SPQR_Solvers
using SuiteSparse.CHOLMOD, SparseArrays, LinearAlgebra
using FrameFun: default_threshold
using SuiteSparse.SPQR: _default_tol, ORDERING_DEFAULT
using BasisFunctions: VectorizingSolverOperator, ArrayOperator, GenericSolverOperator, DictionaryOperator
import BasisFunctions: linearized_apply!


export sparseQR_solver
"""
    sparseQR_solver(op::DictionaryOperator; options...)

A solver that stores the sparse QR factorization of a DictionaryOperator.

Note: the method relies on  `sparse(op)`. It should be efficient.
"""
sparseQR_solver(op::DictionaryOperator; options...) =
    GenericSolverOperator(op, sparseqr_factorization(op; options...))

sparseqr_factorization(op::DictionaryOperator; threshold=default_threshold(op), options...) = qr(sparse(op).A;tol=threshold)


"""
    struct SPQR_Solver{T <: CHOLMOD.VTypes} <: VectorizingSolverOperator{T}

This is no true Solver object. It does not store a factorization during creation.
It applies `spqr_solve` (see MATLAB)` at every application.
"""
struct SPQR_Solver{T <: CHOLMOD.VTypes} <: VectorizingSolverOperator{T}
    op   :: ArrayOperator{T}

    src_linear::Vector{T}
    dest_linear::Vector{T}

    threshold::T

    function SPQR_Solver(op::DictionaryOperator{T}; verbose=false, threshold=default_threshold(op)) where T
        A = sparse(op)
        verbose && @info "SPQR_Solver: create sparse matrix with $(nnz(A.A)) entries"
        droptol!(A.A,threshold)
        verbose && @info "SPQR_Solver: drop entries smaller than $(threshold): $(nnz(A.A)) entries"
        new{T}(A, Vector{T}(undef, size(op, 1)), Vector{T}(undef, size(op, 2)), threshold)
    end
end

export SPQR_solver
SPQR_solver(op::DictionaryOperator; threshold=default_threshold(op), verbose=false, options...) =
    SPQR_Solver(op; threshold=threshold, verbose=verbose)

function linearized_apply!(op::SPQR_Solver, dest, src)
    x = spqr_solve(op.op.A, src; tol=op.threshold)
    copyto!(dest, x)
end

function spqr_solve(A::SparseMatrixCSC{Tv}, b::Vector{Tv}; tol = _default_tol(A)) where {Tv <: CHOLMOD.VTypes}
    X = Ref{Ptr{CHOLMOD.C_Dense{Tv}}}()
    AA = Sparse(A,0)
    B = Dense(b)

    r = _spqr_solve!(ORDERING_DEFAULT, tol, 0, 2, AA,
        C_NULL, B, C_NULL, X,
        C_NULL, C_NULL, C_NULL, C_NULL, C_NULL)
    Vector(Dense(X[]))
end

function _spqr_solve!(ordering::Integer, tol::Real, econ::Integer, getCTX::Integer,
        A::Sparse{Tv},
        Bsparse::Union{Sparse{Tv}                      , Ptr{Cvoid}} = C_NULL,
        Bdense::Union{Dense{Tv}                        , Ptr{Cvoid}} = C_NULL,
        Zsparse::Union{Ref{Ptr{CHOLMOD.C_Sparse{Tv}}}  , Ptr{Cvoid}} = C_NULL,
        Zdense::Union{Ref{Ptr{CHOLMOD.C_Dense{Tv}}}    , Ptr{Cvoid}} = C_NULL,
        R::Union{Ref{Ptr{CHOLMOD.C_Sparse{Tv}}}        , Ptr{Cvoid}} = C_NULL,
        E::Union{Ref{Ptr{CHOLMOD.SuiteSparse_long}}    , Ptr{Cvoid}} = C_NULL,
        H::Union{Ref{Ptr{CHOLMOD.C_Sparse{Tv}}}        , Ptr{Cvoid}} = C_NULL,
        HPinv::Union{Ref{Ptr{CHOLMOD.SuiteSparse_long}}, Ptr{Cvoid}} = C_NULL,
        HTau::Union{Ref{Ptr{CHOLMOD.C_Dense{Tv}}}      , Ptr{Cvoid}} = C_NULL) where {Tv<:CHOLMOD.VTypes}

    AA   = unsafe_load(pointer(A))
    m, n = AA.nrow, AA.ncol
    rnk  = ccall((:SuiteSparseQR_C, :libspqr), CHOLMOD.SuiteSparse_long,
        (Cint, Cdouble, CHOLMOD.SuiteSparse_long, Cint,
         Ptr{CHOLMOD.C_Sparse{Tv}}, Ptr{CHOLMOD.C_Sparse{Tv}}, Ptr{CHOLMOD.C_Dense{Tv}},
         Ptr{Ptr{CHOLMOD.C_Sparse{Tv}}}, Ptr{Ptr{CHOLMOD.C_Dense{Tv}}}, Ptr{Ptr{CHOLMOD.C_Sparse{Tv}}},
         Ptr{Ptr{CHOLMOD.SuiteSparse_long}}, Ptr{Ptr{CHOLMOD.C_Sparse{Tv}}}, Ptr{Ptr{CHOLMOD.SuiteSparse_long}},
         Ptr{Ptr{CHOLMOD.C_Dense{Tv}}}, Ptr{Cvoid}),
        ordering,       # all, except 3:given treated as 0:fixed
        tol,            # columns with 2-norm <= tol treated as 0
        econ,           # e = max(min(m,econ),rank(A))
        getCTX,         # 0: Z=C (e-by-k), 1: Z=C', 2: Z=X (e-by-k)
        A,              # m-by-n sparse matrix to factorize
        Bsparse,        # sparse m-by-k B
        Bdense,         # dense  m-by-k B
        # /* outputs: */
        Zsparse,        # sparse Z
        Zdense,         # dense Z
        R,              # e-by-n sparse matrix */
        E,              # size n column perm, NULL if identity */
        H,              # m-by-nh Householder vectors
        HPinv,          # size m row permutation
        HTau,           # 1-by-nh Householder coefficients
        CHOLMOD.common_struct) # /* workspace and parameters */

    if rnk < 0
        error("Sparse QR factorization failed")
    end

    return rnk
end

end
