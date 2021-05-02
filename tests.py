from pc_harjoitus import xro2Matrix, rredMod2, xSols, printSolutions

# (a)
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

print("(a)\n")
augmentedBinaryMatrix, columns = xro2Matrix(omega, alpha)
solution = xSols(rredMod2(augmentedBinaryMatrix))
printSolutions(solution[0], columns)

# (b)
omega = {
    0: (1, 2, 5),
    1: (2, 3),
    2: (3, 4),
    3: (1, 4, 5)
}
alpha = {
    0: True,
    1: True,
    2: True,
    3: True
}

print("(b)\n")
augmentedBinaryMatrix, columns = xro2Matrix(omega, alpha)
solution = xSols(rredMod2(augmentedBinaryMatrix))
printSolutions(solution[0], columns)

# (c)
omega = {
    0: (1, 2, 6),
    1: (2, 3),
    2: (3, 4),
    3: (4, 5),
    4: (1, 5, 6)
}
alpha = {
    0: True,
    1: True,
    2: True,
    3: True,
    4: True
}

print("(c)\n")
augmentedBinaryMatrix, columns = xro2Matrix(omega, alpha)
solution = xSols(rredMod2(augmentedBinaryMatrix))
printSolutions(solution[0], columns)

# (d)
omega = {
    0: (1, -2, 5),
    1: (2, -3),
    2: (3, -4),
    3: (1, 4, -5)
}
alpha = {
    0: True,
    1: True,
    2: True,
    3: True
}

print("(d)\n")
augmentedBinaryMatrix, columns = xro2Matrix(omega, alpha)
solution = xSols(rredMod2(augmentedBinaryMatrix))
printSolutions(solution[0], columns)