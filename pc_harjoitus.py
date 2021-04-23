from copy import deepcopy

omega = {
    0: [1, 2, -3, 4],
    1: [1, 3, 4],
    2: [3, -3]
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
    
    A = []
    all_variables = []
    
    # Collect and order all of the variables in omega
    # The negations should appear after the corresponding
    # literal, i.e. 1, 2, 3, -3, 4, 5, -5 6, 7
    for key, variables in omega.items():
        for var in variables:
            # The variable should only appear once in the list
            if var not in all_variables:
                if not all_variables:
                    all_variables.append(var)
                elif abs(var) > abs(all_variables[-1]):
                    all_variables.append(var)
                else:
                    for i, value in enumerate(all_variables):
                        if abs(var) < abs(value):
                            all_variables.insert(i, var)
                            break
                        elif abs(var) == abs(value): # if negation of a literal
                            all_variables.insert(i if var > value else i + 1, var)
                            break
    columns = dict()
    
    for col, var in enumerate(all_variables):
        columns[var] = col
        
    column_count = len(all_variables) + 1
    
    for key, variables in omega.items():
        row = [0] * column_count
               
        for var in variables:
            row[columns[var]] = 1
        
        row[column_count - 1] = 1 if alpha[key] else 0
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
    
print(rredMod2(xro2Matrix(omega, alpha)))