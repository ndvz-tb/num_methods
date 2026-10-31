import math
import matplotlib.pyplot as plt

Umax = 100
f = 50
R1 = 5
R2 = 4
R3 = 7
L1 = 0.01
L2 = 0.02
C1 = 300e-6

t_end = 0.2
h = 0.00001
steps = int(t_end / h)

x = [0.0, 0.0, 0.0]
t = 0.0

t_list = []
u2_list = []

for step in range(steps):
    u1 = Umax * math.sin(2 * math.pi * f * t)

    dx0 = (u1 - x[0] * (R1 + R2) + x[2] * R2) / L1
    dx1 = x[2] / C1
    dx2 = (x[0] * R2 - x[1] - x[2] * (R2 + R3)) / L2

    x[0] = x[0] + h * dx0
    x[1] = x[1] + h * dx1
    x[2] = x[2] + h * dx2

    t = t + h
    u2 = x[2] * R3

    t_list.append(t)
    u2_list.append(u2)

print("Кінцевий час t =", round(t, 4))
print("x[0] (i1) =", x[0])
print("x[1] (uC1) =", x[1])
print("x[2] (i3) =", x[2])
print("Вихідна напруга U2 =", u2)

plt.plot(t_list, u2_list)
plt.title("Перехідний процес вихідної напруги U2(t)")
plt.xlabel("t, с")
plt.ylabel("U2, В")
plt.grid(True)
plt.show()
