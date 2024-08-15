# First we read the testset.
n = int(input())


# Define what we need.
def dayornight(n):
    if n == 1:
        return 'Day'
    else:
        return 'Night'
    

# Finally, we print the answer of the testset.
answer = 'Hello ' + dayornight(n)
print(answer)