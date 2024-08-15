# Read the data in Domjudge
def read_data() -> list:
    # Length of data
    n = int(input())

    # Read all lines
    lines = []
    for _ in range(2 * n):
        data = input()
        line = data.rstrip("\n")
        lines += [[line]]
    
    return n, lines


# Read the data
n, lines = read_data()

# Get matrix A
Astr = lines[0:n]
A = []
for line in Astr:
    A += [list(map(int, line[0].split(" ")))]

    
# Get matrix T
Tstr = lines[n:]
T = []
for line in Tstr:
    T += [list(map(int, line[0].split(" ")))]

# Get the minor of a matrix
def minor(M, i, j):
    Mcopy = M[:i] + M[i+1:]
    for i, item in enumerate(Mcopy):
        Mcopy[i] = item[:j] + item[j+1:]
    
    return Mcopy


def firstminor(M):
    Mcopy = M[1:]
    for i, item in enumerate(Mcopy):
        Mcopy[i] = item[1:]
    
    return Mcopy


# Special double minor
def doubleminor(M):
    Mcopy = M[2:]
    for i, item in enumerate(Mcopy):
        Mcopy[i] = item[2:]
    
    return Mcopy




# Calculate determinant of integer matrix
def det(M):
    total = 0
    if len(M) == 2:
        return M[0][0] * M[1][1] - M[1][0] * M[0][1]

    row = M[0]
    for j, item in enumerate(row):
        if item != 0:
            total += (-1)**(j+1) * item * det(minor(M, 0, j))
    
    return total



# Calculate determinant of integer tridiagonal matrix
def tridet(M):
    if len(M) == 3:
        return M[0][0] * (M[1][1] * M[2][2] - M[2][1] * M[1][2]) - M[0][1] * M[1][0] * M[2][2]

    row = M[0]
    item1 = row[0]
    item2 = row[1]

    return item1 * tridet(minor(M, 0, 0)) - item2 * tridet(minor(M, 0, 1))



# Calculate determinant of integer tridiagonal matrix
def tridetopt(M):
    if len(M) == 2:
        return M[0][0] * M[1][1] - M[1][0] * M[0][1]
    if len(M) == 3:
        return M[0][0] * (M[1][1] * M[2][2] - M[2][1] * M[1][2]) - M[0][1] * M[1][0] * M[2][2]

    row = M[0]
    item1 = row[0]
    item2 = row[1]
    item3 = M[1][0]

    return item1 * tridetopt(firstminor(M)) - item2 * item3 * tridetopt(doubleminor(M))


temp = dict()
def tridetit(M, k, n):
    global temp
    maxk = len(temp)
    if k <= maxk:
        return temp[k-1]
    elif k == maxk + 1:
        item1 = M[n-k][n-k]
        item2 = M[n-k][n-k+1]
        item3 = M[n-k+1][n-k]

        temp += [item1 * tridetit(M, k-1, n) - item2 * item3 * tridetit(M, k-2, n)]
        return temp[-1]

    item1 = M[n-k][n-k]
    item2 = M[n-k][n-k+1]
    item3 = M[n-k+1][n-k]

    return item1 * tridetit(M, k-1, n) - item2 * item3 * tridetit(M, k-2, n)


temp += [T[n-1][n-1]]
temp += [T[n-2][n-2] * T[n-1][n-1] - T[n-2][n-1] * T[n-1][n-2]]
print(abs(tridetit(T, len(T), len(T))))
