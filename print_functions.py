import math

def getNumberOfChars(number):
    '''
    Helper function which returns the number of characters in a number
    '''

    if number > 0:
        return int(math.log10(number)) + 1
    if number < 0:
        return int(math.log10(-number)) + 2
    return 1 # number is 0

def printSolutions(solutionsMatrix, columns):
    '''
    Helper function which prints the solution in a clearer form
    '''

    if not solutionsMatrix:
        print("No solutions")
        return

    variables = []
    maxSpacing = 1

    for var, index in columns.items():
        variables.insert(index, var)
        numberOfChars = getNumberOfChars(var)

        if numberOfChars > maxSpacing:
            maxSpacing = numberOfChars

    variablesCount = len(variables)
    solutionsCount = len(solutionsMatrix[0])

    for var in variables:
        print(var, end=' ' * (maxSpacing - getNumberOfChars(var) + 1))
    print('\n')

    for i, var in enumerate(variables):
        numberOfChars = getNumberOfChars(var)
        spaceLine = '' if i == len(variables) - 1 else '-' * (maxSpacing - numberOfChars + 1)
        print('-' * numberOfChars, end=spaceLine)
    print('\n')

    spacing = ' ' * maxSpacing

    for i in range(solutionsCount):
        for j in range(variablesCount):
            print(solutionsMatrix[j][i], end=spacing)
        print('\n')