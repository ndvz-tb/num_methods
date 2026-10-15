import sys
import numpy as np

# Забезпечення коректного виведення в консолі Windows
if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

# =====================================================================
# Лабораторна робота №3 з дисципліни "Чисельні методи"
# Тема: Системи нелінійних рівнянь. Метод Ньютона та е-алгоритм
# Варіант: №17 (Група 2)
# Метод: Векторний е-алгоритм (епсилон-алгоритм)
# =====================================================================

# 1. Початкова система нелінійних рівнянь F(X) = 0 (Варіант 17)
#    f1(x1, x2) = x1 / (x1^2 + x2^2) - 0.4 - x1 = 0
#    f2(x1, x2) = 1.4 - x2 / (x1^2 + x2^2) - x2 = 0
def F(x):
    x1, x2 = x[0], x[1]
    r2 = x1**2 + x2**2
    f1 = (x1 / r2) - 0.4 - x1
    f2 = 1.4 - (x2 / r2) - x2
    return np.array([f1, f2], dtype=float)

# 2. Система, приведена до вигляду зручного для простої ітерації: X = G(X)
#    x1 = g1(x1, x2) = x1 / (x1^2 + x2^2) - 0.4
#    x2 = g2(x1, x2) = 1.4 - x2 / (x1^2 + x2^2)
def G(x):
    x1, x2 = x[0], x[1]
    r2 = x1**2 + x2**2
    g1 = (x1 / r2) - 0.4
    g2 = 1.4 - (x2 / r2)
    return np.array([g1, g2], dtype=float)

def vector_epsilon_algorithm(x0, eps_percent=1e-5, q=2, p=2, max_stages=100):
    """
    Векторний е-алгоритм для розв'язування системи нелінійних рівнянь.
    x0: початкове наближення
    eps_percent: відносна похибка у відсотках (за умовою 10^-5 %)
    q: параметр алгоритму (розмірність q=2)
    p: кількість початкових ітерацій (p=2)
    """
    m = len(x0)          # розмірність системи (m = 2)
    n = 2 * q + 1        # розмір Sn-послідовності (n = 5)
    
    # 3D-матриця розмірністю (n+1) x (n+1) x m
    # Перший індекс: j (від 0 до n)
    # Другий індекс: k+1 (від 0 до n, де k_idx=0 відповідає k=-1, k_idx=1 відповідає k=0, тощо)
    # Третій індекс: i (компоненти вектора, від 0 до m-1)
    e = np.zeros((n + 1, n + 1, m), dtype=float)
    
    # Завдання нулів для k = -1 (другий індекс = 0)
    for j in range(n + 1):
        for i in range(m):
            e[j, 0, i] = 0.0
            
    x = x0.copy()
    total_evals = 0
    converged = False
    
    print("=" * 65)
    print("ПОЧАТОК ОБЧИСЛЕНЬ ВЕКТОРНИМ ЕПСИЛОН-АЛГОРИТМОМ")
    print(f"Початкове наближення: x0 = {x0}")
    print(f"Параметри: m = {m}, q = {q}, p = {p}, n = {n}, eps = {eps_percent}%")
    print("=" * 65)
    
    for stage in range(1, max_stages + 1):
        print(f"\n>>> Стадія екстраполяції №{stage}")
        print(f"    Вхідне наближення стадії: x = [{x[0]:.8f}, {x[1]:.8f}]")
        
        # 1. Початкові ітерації: j = 1..p
        for j in range(p):
            x = G(x)
            total_evals += 1
        print(f"    Після p={p} початкових ітерацій: x = [{x[0]:.8f}, {x[1]:.8f}]")
        
        # 2. Встановлення s0 = x у e[0, 1]
        for i in range(m):
            e[0, 1, i] = x[i]
            
        cond = 0
        # 3. Генерування S-послідовності: j = 0..(2q-1)
        for j in range(2 * q):
            x = G(x)
            total_evals += 1
            for i in range(m):
                e[j + 1, 1, i] = x[i]
                
            # Перевірка збіжності на першому кроці (j == 0)
            if j == 0:
                rel_errors = []
                for i in range(m):
                    denom = e[1, 1, i]
                    if denom != 0:
                        rel_err = abs((e[1, 1, i] - e[0, 1, i]) / denom) * 100.0
                    else:
                        rel_err = abs(e[1, 1, i] - e[0, 1, i]) * 100.0
                    rel_errors.append(rel_err)
                    if rel_err > eps_percent:
                        cond = 1
                print(f"    Відносна похибка кроку: e1 = {rel_errors[0]:.8e}%, e2 = {rel_errors[1]:.8e}%")
                
        # Перевірка умови виходу: cond == 0
        if cond == 0:
            print(f"\n[+] Умова збіжності досягнута на стадії №{stage}!")
            converged = True
            root = e[0, 1].copy()
            break
            
        # 4. Екстраполяція (обчислення е-таблиці з обертанням Самельсона)
        for k in range(1, n):
            for j in range(n - k):
                # Різниця векторів: V = e[j+1, k] - e[j, k]
                V = e[j + 1, k] - e[j, k]
                
                # Квадрат норми вектора (сума квадратів компонент)
                norm_sq = np.sum(V * V)
                if norm_sq == 0:
                    raise ZeroDivisionError("Помилка: нульова норма вектора в обертанні Самельсона!")
                    
                # Обернений вектор за Самельсоном: V / norm_sq
                V_inv = V / norm_sq
                
                # е-формула: e[j, k+1] = e[j+1, k-1] + V_inv
                e[j, k + 1] = e[j + 1, k - 1] + V_inv
                
        # Нове наближення береться з верхівки е-ромба: e[0, n]
        x = e[0, n].copy()
        print(f"    Результат екстраполяції e[0, {n}]: x = [{x[0]:.8f}, {x[1]:.8f}]")
    else:
        print("\n[-] Перевищено максимальну кількість стадій!")
        root = x.copy()
        
    return root, stage, total_evals, converged

# --- Головна частина програми ---
if __name__ == "__main__":
    # Вхідні дані з Таблиці 1 та Таблиці 2 (варіант №17, Група 2)
    x_start = np.array([0.5, 0.5], dtype=float)
    epsilon_perc = 1e-5  # 10^-5 %
    q_param = 2
    p_param = 2
    
    root, stages_count, total_iters, is_converged = vector_epsilon_algorithm(
        x0=x_start,
        eps_percent=epsilon_perc,
        q=q_param,
        p=p_param
    )
    
    # Перевірка вірності роботи алгоритму
    f_residuals = F(root)
    
    print("\n" + "=" * 65)
    print("РЕЗУЛЬТАТИ РОБОТИ АЛГОРИТМУ:")
    print("Знайдений розв'язок системи:")
    print(f"  x1 = {root[0]:.10f}")
    print(f"  x2 = {root[1]:.10f}")
    print(f"Кількість виконаних стадій екстраполяції: {stages_count}")
    print(f"Загальна кількість обчислень вектора функцій: {total_iters}")
    print("\nПеревірка (підстановка розв'язку у початкову систему F(X) = 0):")
    print(f"  f1(x1, x2) = {f_residuals[0]:.10e}")
    print(f"  f2(x1, x2) = {f_residuals[1]:.10e}")
    print("=" * 65)
