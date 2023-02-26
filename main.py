#LU Gauss

A = [ [2, -7, 8, -4],
     [0, -1, 4, -1],
     [3, -4, 2, -1],
     [-9, 1, -4, 6] ]
E = [ [1, 0, 0, 0],
      [0, 1, 0, 0],
      [0, 0, 1, 0],
      [0, 0, 0, 1]]

b = [57, 24, 18, 12]

# A = [ [10, 1, 1],
#       [2, 10, 1],
#       [2, 2, 10]]
# E = [ [1, 0, 0],
#       [0, 1, 0],
#       [0, 0, 1]]
# b = [12, 13, 14]


z = [0] * len(b)
x = [0] * len(b)
n = len(A)
for i in range(n):
    elem = A[i][i]

    for j in range(i+1, n):
        coeff = A[j][i]/elem
        for k in range(0, n):
            if k >= i:
                A[j][k] = A[j][k] - A[i][k] * coeff #Условие здесь, чтобы не залезать на матрицу L
            E[j][k] = E[j][k] - E[i][k] * coeff
        #A[j][i:] = [a-b for a, b in zip(A[j][i:], [t * coeff for t in A[i][i:]])]
        A[j][i] = coeff

det = 1
for i in range(n):
    det *= A[i][i]



for i in range(n):
    s = 0
    for k in range(0, i):
        s += z[k]*A[i][k]
    z[i] = b[i] - s
    #z.append(b[i] - sum([z[k]*A[i][k] for k in range(0, i)]))

for i in range(n-1, -1, -1):
    ui = A[i][i]
    s = 0
    for k in range(i+1, n):
        s += x[k] * A[i][k]
    x[i] = (z[i] - s) / ui
    #x.insert(0, (z[i] - sum([x[k]*A[i][k] for k in range(i+1, n)]))/ui)

print("X =", x)
print("Det =", det)
if det > 0:
    print("Det > 0 => reversable")

for i in range(n):
    for j in range(n-1, -1, -1):
        ui = A[j][j]
        s = 0
        for k in range(j+1, n):
            s += E[k][i] * A[j][k]
        E[j][i] = round((E[j][i] - s) / ui, 4)

print("reverse A matrix:")
for i in E:
    print(i)