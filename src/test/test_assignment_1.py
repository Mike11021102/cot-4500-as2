import numpy as np

#Question 1
x = [ 3.6, 3.8, 3.9]
val = [ 1.675, 1.436, 1.318]
w = 3.7

lenx = len(x)
neville = np.zeros((lenx, lenx))

for i in range(lenx):
    neville[i][0] = val[i]

for j in range(1, lenx):
    for i in range(j, lenx):
        term1 = (w - x[i - j]) * neville[i][j - 1]
        term2 = (w - x[i]) * neville[i - 1][j - 1]

        neville[round(i , 7)][round(j, 7)] = (term1 - term2) / (x[i] - x[i - j])
        

print(neville[lenx - 1][lenx - 1])
print(" ")

#Question 2
xi = [ 7.2, 7.4, 7.5, 7.6]
fxi = [23.5492, 25.3913, 26.8224, 27.4589]
lim = len(xi)


diffs = np.zeros((lim ,lim))


for i in range(lim):
    diffs[i][0] = fxi[i]

     
for j in range(1, lim):
    for i in range(j, lim):
        diffs[i][j] = (diffs[i][j - 1] - diffs[i - 1][j - 1]) / (xi[i] - xi[i - j])


# Extract coefficients from the first row of each column
coefficients = [diffs[i][i] for i in range(lim)]

# Print coefficients of Newton’s Interpolation Polynomial
print(coefficients[1])
print(coefficients[2])
print(coefficients[3])
print()

#Question 3
def newton_approx(x, degree):
    result = coefficients[0]
    term = 1

    for i in range(1, degree + 1):
        term *= (x - xi[i - 1])
        result += coefficients[i] * term
    
    return result
print(newton_approx(7.3, 3))
print()

#Question 4

def hermite_divided_differences(x, f, f_prime):
   
    
    n = len(x)
    size = 2 * n  
    table = np.zeros((size, 5))  # Adjust table size to 6x5
    x_hermite = np.zeros(size)  # Expanded x-values with repetitions
    
    # Step 1: Populate x-values with repeated entries
    for i in range(n):
        x_hermite[2 * i] = x[i]
        x_hermite[2 * i + 1] = x[i]

    # Step 2: Populate f(x) values and first-order divided differences
    for i in range(n):
        table[2 * i, 0] = x_hermite[2 * i]  # Set x value
        table[2 * i + 1, 0] = x_hermite[2 * i + 1]  # Set x value
        table[2 * i, 1] = f[i]  # f(x)
        table[2 * i + 1, 1] = f[i]  # f(x)
        table[2 * i + 1, 2] = f_prime[i]  # f'(x)

        if i > 0:
            table[2 * i, 2] = (f[i] - f[i - 1]) / (x[i] - x[i - 1])  # Regular divided difference

    # Step 3: Compute higher-order divided differences (up to d.d2)
    for j in range(3, 5):  # Compute only up to d.d2 (columns 3 and 4)
        for i in range(size - (j - 1)):
            table[i, j] = (table[i + 1, j - 1] - table[i, j - 1]) / (x_hermite[i + (j - 1)] - x_hermite[i])

    
    return table

# Example Data
x_values = np.array([3.6, 3.8, 3.9])
f_values = np.array([1.675, 1.436, 1.318])
f_prime_values = np.array([-1.195, -1.188, -1.182])

# Compute the Hermite divided difference table
hermite_table = hermite_divided_differences(x_values, f_values, f_prime_values)

# Print the result in the specified format
formatted_table = []

# Add each row of the computed table
for row in hermite_table:
    formatted_row = [f"{value:.10e}" for value in row]  # Format each value to 10 decimal places
    formatted_table.append(formatted_row)

# Print the formatted table row by row
for row in formatted_table:
    print(row)
    
print()

       
#Question 5
print()
x = np.array([2, 5, 8, 10])
fx = np.array([3, 5, 7, 9])


h = np.diff(x)  
n = len(x) - 1  


A = np.zeros((n + 1, n + 1))  

A[0][0] = 1  
A[n][n] = 1  


for i in range(1, n):
    A[i][i - 1] = h[i - 1]  
    A[i][i] = 2 * (h[i - 1] + h[i]) 
    A[i][i + 1] = h[i] 


b = np.zeros(n + 1)  
b[0] = 0
b[n] = 0
b[1] = ((fx[2] - fx[1]) * (3 / h[1])) - ((3 / h[0]) * (fx[1] - fx[0]))
b[2] = ((fx[3] - fx[2]) * (3 / h[2])) - ((3 / h[1]) * (fx[2] - fx[1]))

    
print(A)

# Print Vector b

print(b)


M = np.linalg.solve(A, b)

M_values = np.array([0, M[1], M[2], 0])
print(M_values)



