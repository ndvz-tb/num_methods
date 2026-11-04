import math

def f(x):
    return x**2 - math.cos(x)

def df(x):
    return 2*x + math.sin(x)

def simplified_newton(x0, eps_percent):
    df_x0 = df(x0)
    
    if df_x0 == 0:
        print("Помилка: похідна дорівнює нулю.")
        return None
        
    x_old = x0
    i = 1
    
    for i in range(100):
        fx = f(x_old)
        
        x_new = x_old - fx / df_x0
        
        if x_new != 0:
            error_percent = abs((x_new - x_old) / x_new) * 100
        else:
            error_percent = float('inf')
            
        print("Крок:", i, "| x =", x_new, "| Похибка (%) =", error_percent)
        
        if error_percent < eps_percent:
            break
            
        x_old = x_new
        i += 1
        
    else: 
        print("Кількість ітерацій більша за 100")
            
    return x_new

x_start = -4.0 
epsilon = 0.01 

print("Початкове наближення:", x_start)
print("Задана похибка:", epsilon, "%")
print("Починаємо рахувати...")

root = simplified_newton(x_start, epsilon)

if root is not None:
    print("Знайдений корінь:", root)
    print("Перевірка (підставляємо корінь у функцію):", f(root))