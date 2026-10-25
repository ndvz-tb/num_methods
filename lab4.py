import math

def f(x):
    return 1 / (x**2 * math.sqrt(6 - x**2))

def F(x):
    return -math.sqrt(6 - x**2) / (6 * x)

def simpson(a, b, m):
    integral = 0.0
    n = 2 * m
    h = (b - a) / n
    
    fa = f(a)
    fb = f(b)
    
    for i in range(1, m + 1):
        x = a + (2 * i - 1) * h
        integral += 4 * f(x)
        
    for i in range(1, m):
        x = a + 2 * i * h
        integral += 2 * f(x)
        
    integral = (h / 3) * (fa + fb + integral)
    return integral

a = 1.0
b = 2.0
m = 10
n = 2 * m

print(f"Межі інтегрування: [{a}, {b}]")
print(f"Кількість розбиттів: n = 2m = {n} (m = {m})")
print(f"Крок інтегрування: h = {(b - a) / n}")

approx_val = simpson(a, b, m)
exact_val = F(b) - F(a)
abs_err = abs(exact_val - approx_val)
rel_err = (abs_err / abs(exact_val)) * 100

print("\nРезультати обчислення:")
print(f"Наближене значення (метод Сімпсона):   {approx_val:.10f}")
print(f"Точне значення (Ньютона-Лейбніца):     {exact_val:.10f}")
print(f"Абсолютна похибка:                     {abs_err:.6e}")
print(f"Відносна похибка (%):                  {rel_err:.6e} %")
