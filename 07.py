# Импорт модулей
import matplotlib.pyplot as plt  # модуль для построения графиков
import numpy as np  # модуль для общих математических функций, используя сокращенное название np


def J2eV(E):
    return E * 6.242 * 10 ** 18


m0 = 9.1e-31  # Масса частицы
hbar = 1.0546e-34  # Постоянная Дирака
L = 101e-10  # Длина ямы
Np = 100
# Нахождение Гамильтонианна
Lrr = np.linspace(0, L, Np)
Lr = Lrr[:]
Lan = np.linspace(0, L, 10000)
a = Lr[1] - Lr[0]
t0 = hbar ** 2 / 2 / m0 / a ** 2
H = np.diag(np.ones(len(Lr)) * 2 * t0)
H = H + np.diag(np.diag(H, 1) - t0, 1)
H = H + np.diag(np.diag(H, -1) - t0, -1)
# Получение собственныйх векторов гамильтониана (матрица Psy)
En, Psy = np.linalg.eigh(H)
n = 1
plt.figure(figsize=(20, 10))  # Задание размеров области построения графиков
plt.subplot(2, 2, 2)  # Расположение графика
# Построение графика плотности вероятности от x при n = 1
plt.plot(Lan * pow(10, 9), 2 / L * np.abs(np.sin(np.pi * n * Lan / L)) ** 2, color='red', label="n = 1")
n = 25
# Построение графика плотности вероятности от x при n = 25
plt.plot(Lan * pow(10, 9), 2 / L * np.abs(np.sin(np.pi * n * Lan / L)) ** 2, label="n = 25")
plt.xlabel("x, нм")  # Подпись оси абцисс
plt.ylabel("Плотность вероятности")  # Подпись оси ординат
plt.title("Распределение вероятностей")  # Подпись графика
plt.grid()  # Включение сетки
plt.legend()  # Включение легенды
plt.subplot(2, 2, 1)  # Расположение графика
Eanalytic = []
k = [(i + 1) * np.pi / L for i in range(Np)]
# Вычисление энергии численным методом
for i in range(len(En)):
    En[i] = J2eV(En[i])
# Вычисление энергии аналитическим методом
for i in k:
    Eanalytic.append(J2eV(hbar ** 2 * i ** 2 / 2 / m0))
k1 = [i * 2 for i in k]
plt.plot(k, En, "x", label='numeric', color='red')  # Построение численного графика
plt.plot(k, Eanalytic, label='analitic')  # Построение аналитического графика
plt.xlabel("Собственное число, α")  # Подпись оси абцисс
plt.ylabel("E (эВ)")  # Подпись оси ординат
plt.title("Дисперсионная зависимость")  # Подпись графика
plt.grid()  # Включение сетки
plt.legend()  # Включение легенды
plt.subplot(2, 2, 3)  # Расположение графика
# Вычисление численным методом дисперсионной зависимости
T = np.diag(np.ones(Np) * 7.6265)
t01 = 3.8132
for i in range(len(T)):
    for j in range(len(T[i])):
        if i + 1 == j:
            T[i][j] = -t01
        if i - 1 == j:
            T[i][j] = -t01

T[0][Np - 1] = -t01
T[Np - 1][0] = -t01
D, V = np.linalg.eigh(T)
plt.plot(D, "x", color='red')  # Построение графика
plt.title("Собственные значения энергии")  # Подпись оси абцисс
plt.xlabel("Собственное число, α")  # Подпись оси ординат
plt.ylabel("E (эВ)")  # Подпись графика
plt.grid()  # Включение сетки
plt.show()  # Вывод всех графиков
