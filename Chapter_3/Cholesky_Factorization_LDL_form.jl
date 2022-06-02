using LinearAlgebra

function choleskyFactorize(A_Matrix::Matrix)
    if isposdef(A_Matrix) == false
        println("A_Matrix is not positive definite. \nReturned Value is A_Matrix")
        return A_Matrix
    end
    D_Matrix = Diagonal(A_Matrix)
    L_Matrix = UnitLowerTriangular(A_Matrix)
    C_Matrix = A_Matrix
    for j in range(1, length=length(L_Matrix[:,1]))
        temp_value_j = 0
        for s in range(1, length=j-1)
            temp_value_j += (D_Matrix[s,s]*L_Matrix[j,s]^2)
        end
        C_Matrix[j,j] = A_Matrix[j,j] - temp_value_j
        D_Matrix[j,j] = C_Matrix[j,j]
        for i in range(j+1, length=length(L_Matrix[1,:])-j)
            temp_value_i = 0
            for s in range(1, length=j-1)
                temp_value_i += (D_Matrix[s,s]*L_Matrix[i,s]*L_Matrix[j,s])
            end
            C_Matrix[i,j] = A_Matrix[i,j] - temp_value_i
            L_Matrix[i,j] = C_Matrix[i,j] / D_Matrix[j,j]
        end
    end
    return L_Matrix
end

A = [1. 0. 2.; 1. 2. 0.; 5. 0. 2.]
test_A = A*Transpose(A)
println(choleskyFactorize(test_A))