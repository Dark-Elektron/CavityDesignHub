

M = []

p = 29

if 0 < p <= 6:
    for j in range(2):
        for i in range(3):
            M.append((i, j))

if 6 < p <= 12:
    for j in range(2):
        for i in range(3):
            M.append((i, j))

    for i in range(3):
        for j in range(2):
            M.append((i, j))

if 12 < p <= 18:
    for j in range(2):
        for i in range(3):
            M.append((i, j))

    for i in range(3, 6):
        for j in range(2):
            M.append((i, j))

    for i in range(6):
        M.append((i, 2))

if 18 < p <= 27:
    for j in range(2):
        for i in range(3):
            M.append((i, j))

    for i in range(3, 6): # n, n+3
        for j in range(2):
            M.append((i, j))

    for i in range(6): # n
        M.append((i, 2))

    for i in range(6, 9): # n, n+3
        for j in range(3):
            M.append((i, j))

if 27 < p <= 36:
    for j in range(2):
        for i in range(3):
            M.append((i, j))

    for i in range(3, 6): # n, n+3
        for j in range(2):
            M.append((i, j))

    for i in range(6): # n
        M.append((i, 2))

    for i in range(6, 9): # n, n+3
        for j in range(3):
            M.append((i, j))

    for i in range(9):
        M.append((i, 3))



print(M)