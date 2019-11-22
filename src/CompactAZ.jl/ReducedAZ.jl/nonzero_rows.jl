
function nonzero_rows(m::AbstractMatrix, s=size(m,1);nonzero_tol=0)
    mask = falses(s)
    if nonzero_tol == 0
        nonzero_rows!(mask, m)
    else
        nonzero_rows!(mask, m, nonzero_tol)
    end
    mask
end

function nonzero_rows!(mask::BitArray, m::Matrix)
    for i in 1:size(m,1)
        for j in 1:size(m,2)
            if m[i,j] != 0
                mask[i] = true
                break
            end
        end
    end
end

function nonzero_rows!(mask::BitArray, m::Matrix, tol)
    for i in 1:size(m,1)
        for j in 1:size(m,2)
            if abs(m[i,j]) > tol
                mask[i] = true
                break
            end
        end
    end
end
using SparseArrays: SparseMatrixCSC
function nonzero_rows!(mask::BitArray, m::SparseMatrixCSC)
    mask[m.rowval] .= true
end

function nonzero_rows!(mask::BitArray, m::SparseMatrixCSC, tol)
    for (x,y,v) in zip(findnz(m)...)
        if abs(v) > tol
            mask[x] = true
        end
   end
end
