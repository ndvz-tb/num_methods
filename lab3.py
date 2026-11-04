import numpy as np

def f1(x1, x2):
    return x1 / (x1**2 + x2**2) - 0.4 - x1

def f2(x1, x2):
    return 1.4 - x2 / (x1**2 + x2**2) - x2


def g1(x1, x2):
    return x1 / (x1**2 + x2**2) - 0.4

def g2(x1, x2):
    return 1.4 - x2 / (x1**2 + x2**2)

def G(x):
    return np.array([g1(x[0], x[1]), g2(x[0], x[1])], dtype=float)

def epsilon_algorithm(x0, eps_percent=1e-5, q=2, p=2, max_stages=100):
    m = len(x0)
    n = 2 * q + 1
    e = np.zeros((n + 1, n + 1, m), dtype=float)
    
    x = x0.copy()
    
    for stage in range(1, max_stages + 1):
        for _ in range(p):
            x = G(x)
            
        e[0, 1] = x.copy()
        cond = False
        
        for j in range(2 * q):
            x = G(x)
            e[j + 1, 1] = x.copy()
            
            if j == 0:
                rel_err = np.abs((e[1, 1] - e[0, 1]) / e[1, 1]) * 100
                if np.any(rel_err > eps_percent):
                    cond = True
                    
        cur_err = np.max(rel_err)
        print(f"Стадія: {stage} | x1 = {x[0]:.6f}, x2 = {x[1]:.6f} | Похибка (%) = {cur_err:.6e}")
        
        if not cond:
            print(f"Збіжність досягнута на стадії {stage}")
            return e[0, 1]
            
        for k in range(1, n):
            for j in range(n - k):
                V = e[j + 1, k] - e[j, k]
                e[j, k + 1] = e[j + 1, k - 1] + V / np.sum(V * V)
                
        x = e[0, n].copy()
        
    print("Кількість стадій більша за", max_stages)
    return x


x_start = np.array([0.5, 0.5])
epsilon = 1e-5  

print("Початкове наближення:", x_start)
print("Задана похибка:", epsilon, "%")
print("Починаємо рахувати...")

root = epsilon_algorithm(x_start, epsilon)

if root is not None:
    print("\nЗнайдений розв'язок:")
    print(f"x1 = {root[0]:.8f}")
    print(f"x2 = {root[1]:.8f}")
    print("\nПеревірка (підставляємо корінь у систему):")
    print(f"f1 = {f1(root[0], root[1]):.6e}")
    print(f"f2 = {f2(root[0], root[1]):.6e}")
