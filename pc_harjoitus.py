from copy import deepcopy
from itertools import product

omega = {
    0: (1, 2, -3, 4),
    1: (1, 3, 4),
    2: (3, -3)
}
alpha = {
    0: True,
    1: False,
    2: True
}

def xro2Matrix(omega, alpha):
    '''
    This function takes omega sets and the alpha
    values and produces an augmented binary matrix
    '''

    # Collect all of the variables in omega
    allVariables = {
        var for key, variables in omega.items() for var in variables
    }

    # Sort the variables in the correct column order. The negations should
    # appear after the corresponding literal, i.e. 1, 2, 3, -3, 4, 5, -5, 6, 7
    allVariablesSorted = sorted(allVariables, key=abs)

    # Store the column indexes in a dictionary for faster access
    # (constant time). Using the list function index() would be O(n).
    columns = {
        var: col for col, var in enumerate(allVariablesSorted)
    }

    columnCount = len(columns) + 1 # Add one for the alpha column
    A = [] # Augmented binary matrix
    alphaColumn = columnCount - 1

    for rowIdx, variables in omega.items():
        row = [0] * columnCount
               
        for var in variables:
            row[columns[var]] = 1
        
        row[alphaColumn] = 1 if alpha[rowIdx] else 0
        A.append(row)
        
    return A
        
def replaceRowValues(A, cols, targetRow, otherRow):
    '''
    This function replaces the target row values by the
    sum of the target row values and some other row values.
    '''
    for col in range(cols):
        A[targetRow][col] = (A[otherRow][col] + A[targetRow][col]) % 2
        
def rredMod2(A):
    '''
    Transforms a matrix to row reduced echelon form by performing Gaussian
    elimination. This function only considers matrixes with binary values.
    '''

    if A == None or A == []:
        return

    B = deepcopy(A)
    rows = len(B)
    cols = len(B[0])
    
    # Traverse the matrix diagonally. In other words the row index and the
    # columns index for the pivot is always the same.
    for pivotIndex in range(rows):
        pivot = B[pivotIndex][pivotIndex]
        
        # Make all of the values below the pivot 0
        for row in range(pivotIndex + 1, rows):
            if B[row][pivotIndex] == 1:
                if pivot == 1:
                    replaceRowValues(B, cols, row, pivotIndex)
                else:
                    # If pivot is 0, swap the values of
                    # the pivot row and the current row
                    for col in range(cols):
                        temp = A[pivotIndex][col]
                        A[pivotIndex][col] = A[row][col]
                        A[row][col] = temp
    
    # Traverse the matrix again but this make
    # all of the values above the pivot 0.
    for pivotIndex in range(rows):        
        for row in range(pivotIndex - 1, -1, -1):
            if B[row][pivotIndex] == 1:
                # This time ignore a pivot with the value 0
                replaceRowValues(B, cols, row, pivotIndex)
    return B
    
def xSols(B):
    '''
    This function takes in a row-reduced augmented matrix and checks whether it
    has a solution and if so, returns True and a matrix with all the possible
    solutions. If there is no solution, False and an empty matrix is returned.

    The rows of the matrix represent a variable (literal) and the columns
    different solutions (1 or 0). For example, consider the following matrix:

    [[1, 0, 0, 1], [0, 0, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1], [0, 1, 0, 1]]

    This means that the first solution is:
    x1 = 1, x2 = 0, x3 = 1, x4 = 0, x5 = 0

    and the second solutions is:
    x1 = 0, x2 = 0, x3 = 1, x4 = 0, x5 = 1

    and so on.
    '''

    if B == None or B == []:
        return [], False

    pivotCols = []
    freeVariables = {} # Variables corresponding to the non-pivot columns

    for rowIdx, row in enumerate(B):
        freeVariables[rowIdx] = set()

        # Scan each row for the pivot columns and free variables
        for colIdx, variable in enumerate(row):
            if variable == 1:
                # Pivot found

                if colIdx == len(row) - 1:
                    return [], False # 0 = 1 -> contradiction

                pivotCols.append(colIdx)

                # Move to the next column and collect the free variables
                for col in range(colIdx + 1, len(row) - 1):
                    if row[col] == 1:
                        # Row has the variable
                        freeVariables[rowIdx].add(col)
                break

    variableCount = len(B[0]) - 1
    nonPivotCols = [
        col for col in range(variableCount) if col not in pivotCols
    ]

    # Every possible combination of t binary parameters
    tCombinations = list(map(list, product([0, 1], repeat = len(nonPivotCols))))

    # Output matrix
    X = [[] for _ in range(variableCount)]

    for tParams in tCombinations:
        variableValues = [] # x̂ + t1v1 + t2v2 + ... + tnvn

        for rowIdx, row in enumerate(B):
            result = row[-1]

            for i, t in enumerate(tParams):
                if t == 1 and nonPivotCols[i] in freeVariables[rowIdx]:
                    # x + t2 + t3 + ... + tn = α - mod(t1, 2)
                    result = (result - 1) % 2

            variableValues.append(result)

        variableValues += tParams

        # The result matrix is presented so that the variables
        # are the rows and their results the columns
        for row, value in enumerate(variableValues):
            X[row].append(value)

    return X, True

print(xSols(rredMod2(xro2Matrix(omega, alpha))))