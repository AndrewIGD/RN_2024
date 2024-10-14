def getCoefficient(coeff):
    if coeff == '-':
        return -1
    elif coeff == '+' or coeff == '':
        return 1
    else:
        return int(coeff)

def parseEquation(equation):
    equation = equation.replace(" ", "")
    left, right = equation.split("=")

    coefficients = []
    constant = int(right)

    prevIndex = 0
    for i in range(len(left)):
        if "+-0123456789".find(left[i]) == -1:
            coefficients.append(getCoefficient(left[prevIndex:i]))
            prevIndex = i + 1

    return [coefficients, constant]

def getMatrixFromFile(file):
    A = []
    B = []

    f = open(file, "r")
    for line in f:
        coefficients, constant = parseEquation(line)
        A.append(coefficients)
        B.append(constant)

    return [A, B]

cached_determinants = {}
def getDeterminant(matrix):
    matrixStr = str(matrix)
    if matrixStr in cached_determinants:
        return cached_determinants[matrixStr]

    def calculateDeterminant(matrix):
        if len(matrix) == 2:
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]

        if len(matrix) == 3:
            return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) - matrix[0][1] * (
                    matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) + matrix[0][2] * (
                    matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])

        determinant = 0
        for i in range(len(matrix)):
            determinant += (-1) ** i * matrix[0][i] * getDeterminant(minor(matrix, 0, i))

        return determinant

    determinant = calculateDeterminant(matrix)
    cached_determinants[matrixStr] = determinant
    return determinant

def trace(matrix):
    sum = 0
    for i in range(len(matrix)):
        sum += matrix[i][i]

    return sum

def norm(vector):
    sum = 0
    for i in range(len(vector)):
        sum += vector[i] ** 2

    return sum ** 0.5

def transpose(matrix):
    transposed = []
    for i in range(len(matrix)):
        vector = []
        for j in range(len(matrix)):
            vector.append(matrix[j][i])

        transposed.append(vector)

    return transposed

def multiply(matrix, vector):
    newVector = []
    for i in range(len(matrix)):
        sum = 0
        for j in range(len(matrix)):
            sum = sum + matrix[i][j] * vector[j]

        newVector.append(sum)

    return newVector

def solveByCramer():
    A, B = getMatrixFromFile("equations.txt")
    determinant = getDeterminant(A)

    sol_vector = []
    for i in range(len(A)):
        new_matrix = []
        for j in range(len(B)):
            new_matrix.append(A[j][0:i] + [B[j]] + A[j][i+1:])

        new_determinant = getDeterminant(new_matrix)
        sol_vector.append(new_determinant/determinant)

    return sol_vector

def minor(matrix, i, j):
    new_matrix = matrix[0:i] + matrix[i+1:]

    for k in range(len(new_matrix)):
        new_matrix[k] = new_matrix[k][0:j] + new_matrix[k][j+1:]

    return new_matrix

def scalarMultiply(matrix, scalar):
    new_matrix = []
    for i in range(len(matrix)):
        vector = []
        for j in range(len(matrix)):
            vector.append(matrix[i][j] * scalar)

        new_matrix.append(vector)

    return new_matrix

def cofactorMatrix(matrix):
    new_matrix = []
    for i in range(len(matrix)):
        vector = []
        for j in range(len(matrix)):
            vector.append((-1) ** ((i+1) + (j+1)) * getDeterminant(minor(matrix, i, j)))

        new_matrix.append(vector)

    return new_matrix

def solveByInversion():
    A, B = getMatrixFromFile("equations.txt")
    adjointMatrix = transpose(cofactorMatrix(A))
    inverseMatrix = scalarMultiply(adjointMatrix, 1 / getDeterminant(A))

    return multiply(inverseMatrix, B)

print(solveByCramer())
print(solveByInversion())
