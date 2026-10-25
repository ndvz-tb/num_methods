import math

def f(x):
    return 1 / (x**2 * math.sqrt(6 - x**2))

def F(x):
    return -math.sqrt(6 - x**2) / (6 * x)

a = 1
b = 2
m = 10
n = 2 * m
h = (b - a) / n

integral = 0
fa = f(a)
fb = f(b)

for i in range(1, m + 1):
    x = a + (2 * i - 1) * h
    integral = integral + 4 * f(x)

for i in range(1, m):
    x = a + 2 * i * h
    integral = integral + 2 * f(x)

integral = (h / 3) * (fa + fb + integral)
exact = F(b) - F(a)

print("Метод Сімпсона:", integral)
print("Точне значення (Ньютона-Лейбніца):", exact)
