import numpy as np
import matplotlib.pyplot as plt

C1 = 6.2e-3
C2 = 55e-3
Lmax = 1.2
Lmin = 0.12
imin = 1.0
imax = 2.0
R1 = 3.0
R2 = 1.1
R3 = 630.0
R4 = 23.0
a = 0.04
T = 4 * a

def u1_func(t):
    tau = t % T
    if tau < a:
        return 10.0 * (tau / a)
    elif tau < 2 * a:
        return 10.0 * (1.0 - (tau - a) / a)
    else:
        return -10.0 * ((tau - 2 * a) / (2 * a))

def L2_func(i2):
    x = abs(i2)
    if x <= imin:
        return Lmax
    elif x >= imax:
        return Lmin
    else:
        hs = imax - imin
        m1 = 0.0
        m2 = 0.0
        b1 = (2 * (x - imin) + hs) * (imax - x)**2
        b2 = (2 * (imax - x) + hs) * (x - imin)**2
        b3 = (x - imin) * (imax - x)**2
        b4 = (x - imax) * (x - imin)**2
        return (b1 * Lmax + b2 * Lmin) / (hs**3) + (b3 * m1 + b4 * m2) / (hs**2)

def f_sys(t, y):
    uc1, i2, uc2 = y
    u1 = u1_func(t)
    i1 = (u1 - uc1 + i2 * R3) / (R1 + R3)
    ur3 = (i1 - i2) * R3

    duc1 = i1 / C1
    di2 = (ur3 - i2 * (R2 + R4) - uc2) / L2_func(i2)
    duc2 = i2 / C2
    return np.array([duc1, di2, duc2], dtype=float)

h = T / 400.0
t_end = 5 * T
steps = int(round(t_end / h))

t = 0.0
y = np.array([0.0, 0.0, 0.0], dtype=float)

t_list = []
u1_list = []
uc1_list = []
i2_list = []
uc2_list = []
u2_list = []

for step in range(steps + 1):
    if step % 4 == 0:
        u1 = u1_func(t)
        u2 = y[1] * R4 + y[2]
        t_list.append(t)
        u1_list.append(u1)
        uc1_list.append(y[0])
        i2_list.append(y[1])
        uc2_list.append(y[2])
        u2_list.append(u2)

    if step == steps:
        break

    K1 = h * f_sys(t, y)
    K2 = h * f_sys(t + h / 3.0, y + K1 / 3.0)
    K3 = h * f_sys(t + 2.0 * h / 3.0, y - K1 / 3.0 + K2)
    K4 = h * f_sys(t + h, y + K1 - K2 + K3)

    y = y + (K1 + 3.0 * K2 + 3.0 * K3 + K4) / 8.0
    t = t + h

print("Кінцевий час t =", round(t, 4))
print("uC1 =", y[0])
print("i2 =", y[1])
print("uC2 =", y[2])
print("Вихідна напруга U2 =", u2_list[-1])

plt.figure(figsize=(10, 8))

plt.subplot(3, 1, 1)
plt.plot(t_list, u1_list, label="u1(t)")
plt.plot(t_list, u2_list, label="u2(t)")
plt.ylabel("U, В")
plt.legend()
plt.grid(True)

plt.subplot(3, 1, 2)
plt.plot(t_list, uc1_list, label="uC1(t)")
plt.plot(t_list, uc2_list, label="uC2(t)")
plt.ylabel("uC, В")
plt.legend()
plt.grid(True)

plt.subplot(3, 1, 3)
plt.plot(t_list, i2_list, label="i2(t)", color="green")
plt.xlabel("t, с")
plt.ylabel("i2, А")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
