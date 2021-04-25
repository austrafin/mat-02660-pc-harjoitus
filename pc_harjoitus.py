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

    allVarsUnsorted = set()
    
    # Collect and order all of the variables in omega
    # The negations should appear after the corresponding
    # literal, i.e. 1, 2, 3, -3, 4, 5, -5, 6, 7
    for key, variables in omega.items():
        for var in variables:
            allVarsUnsorted.add(var)
    allVarsSorted = sorted(allVarsUnsorted, key=abs)

    columns = dict()

    for col, var in enumerate(allVarsSorted):
        columns[var] = col

    columnCount = len(allVarsSorted) + 1
    A = []

    for key, variables in omega.items():
        row = [0] * columnCount
               
        for var in variables:
            row[columns[var]] = 1
        
        row[columnCount - 1] = 1 if alpha[key] else 0
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
    Transforms a matrix to row reduced echelon form by performing
    Gaussian elimination. This function only considers matrixes
    with binary values.
    '''
    if A == None or A == []:
        return

    B = deepcopy(A)
    rows = len(B)
    cols = len(B[0])
    
    # Traverse the matrix diagonally. In other words
    # the row index and the columns index for the
    # pivot is always the same.
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
    if B == None or B == []:
        return [], False

    pivotCols = list()
    freeVariables = dict()

    for row_idx, row in enumerate(B):
        freeVariables[row_idx] = set()

        for col_idx, col in enumerate(row):
            if col == 1:
                if col_idx == len(row) - 1:
                    return [], False # contradiction
                pivotCols.append(col_idx)

                # Collect the free variables
                for j in range(col_idx + 1, len(row) - 1):
                    if row[j] == 1:
                        freeVariables[row_idx].add(j)
                break

    nonPivotCols = [item for item in range(len(B[0]) - 1) if item not in pivotCols]
    tCombinations = list(map(list, product([0, 1], repeat = len(nonPivotCols))))
    X = [[] for _ in range((len(B[0]) - 1))]

    for tParams in tCombinations:
        xRow = []

        for row_idx, row in enumerate(B):
            result = row[-1]

            for i, t in enumerate(tParams):
                if t == 1 and nonPivotCols[i] in freeVariables[row_idx]:
                    result = (result - 1) % 2
            xRow.append(result)
        xRow += tParams
        print(xRow)

        for i, value in enumerate(xRow):
            X[i].append(value)

    return X, True

print(xSols(rredMod2(xro2Matrix(omega, alpha))))