# Read the data in Domjudge
def read_data() -> list:
    # Length of data
    n = 3

    # Read all lines
    lines = []
    for _ in range(n):
        data = input()
        line = data.rstrip("\n")
        lines += [[line]]
    
    return n, lines


def paths(vertex, steps):
    global n
    if steps == 1:
        return 1 if vertex == 1 or vertex == n else 2
    else:
        if vertex == 1:
            return paths(vertex + 1, steps-1)
        elif vertex == n:
            return paths(vertex - 1, steps-1)
        else:
            return paths(vertex + 1, steps-1) + paths(vertex - 1, steps-1) 


print(paths(i, k))