import numpy as np

k = 17
p = 22
s = 0.02 * k
B_val = 0.02 * p

A = np.array([
    [8.3,  2.62 + s, 4.1,      1.9],
    [3.92, 8.45,     7.78 - s, 2.46],
    [3.77, 7.21 + s, 8.04,     2.28],
    [2.21, 3.65 - s, 1.69,     6.69]
], dtype=float)

B = np.array([-10.65 + B_val, 12.21, 15.45 - B_val, -8.35], dtype=float)

V = A.copy()
P = B.copy()
n = len(P)


for i in range(n):            
    diag = V[i][i]
    
    for k_col in range(n):
        V[i][k_col] = V[i][k_col] / diag
    P[i] = P[i] / diag
    
    for j in range(i + 1, n):
        factor = V[j][i]
        
        for k in range(n): 
            V[j][k] = V[j][k] - factor * V[i][k]
            
        P[j] = P[j] - factor * P[i]

print("Верхня трикутна матриця C")
print(np.round(V, 4))
print("Змінений вектор правих частин")
print(np.round(P, 4))

X = np.zeros(n)
X[n-1] = P[n-1]

for i in range(n - 2, -1, -1):
    sum_cx = 0
    for j in range(i + 1, n):
        sum_cx += V[i][j] * X[j]
    X[i] = P[i] - sum_cx

print("Знайдені корені")
for i, x in enumerate(X):
    print(f"x{i+1} = {x:.4f}")

print("Перевірка вірності роботи алгоритму")

B_check = np.zeros(n)
for i in range(n):
    for j in range(n):
        B_check[i] += A[i][j] * X[j]

print("Обчислені ліві частини рівнянь:  ", np.round(B_check, 4))
print("Праві значення елементів вектора B:", np.round(B, 4))