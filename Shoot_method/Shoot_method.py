from sympy import *
import numpy as np
import matplotlib.pyplot as plt

def input_data(file_name):   # как обычно.  Это функция для ввода информации
    with open(file_name) as file:
        file.readline()
        n = int(file.readline())   # число уравнений и начальных условий
        file.readline()
        equations = []   #  это правые части
        for i in range(n):
            equations.append(simplify(file.readline().split("=")[-1]))
        file.readline()
        start_conditions = []  # это начальные условия
        for i in range(n):
            start_conditions.append(float(file.readline().split("=")[-1]))
        file.readline()
        segment = list(map(float, file.readline().split()))   #  это отрезок
        file.readline()
        step = float(file.readline())   #  это шаг
        file.readline()
        eps = float(file.readline())   # это точность
        file.readline();
        target = float(file.readline())
        return target, n, equations, start_conditions, segment, step, eps


def find_value(equations, value_x, values_Y): # эта функция подставляет value_x и values_Y в equations
    result = [eq for eq in equations]
    for i in range(len(equations)):
        result[i] = result[i].subs("s", value_x)
        for j in range(len(equations)):
            result[i] = result[i].subs("u" + str(j + 1), values_Y[j])
    return result
    

def solve_system_runge_kutta(n, equations, start_conditions, segment, step):  # это решает систему методом Рунге-Кутты 4 порядка точности и заданным шагом
    values_x = [segment[0] + i * step for i in range(int((segment[1] - segment[0]) / step) + 1)]  # задаем значения аргумента
    values_Y = [[start_conditions[i] for i in range(n)]]   # создаем список значений функции.  Пока известно только начальное
                                                           # состояние.
    curr_value_x = values_x[0] # задаем текущее значение аргумента
    curr_value_Y = [start_conditions[i] for i in range(n)]  # задаем текущие значения функций
    for value in values_x[1:]:   # цикл по всем значениям аргумента
        # вычисляем все k
        k1 = list(map(lambda X: step * X, find_value(equations, curr_value_x, curr_value_Y)))
        var = [curr_value_Y[i] + k1[i] / 3 for i in range(n)]
        k2 = list(map(lambda X: step * X, find_value(equations, curr_value_x + step / 3, var)))
        var = [curr_value_Y[i] - (k1[i] / 3) + k2[i] for i in range(n)]
        k3 = list(map(lambda X: step * X, find_value(equations, curr_value_x + (2 * step) / 3, var)))
        var = [curr_value_Y[i] + k1[i] - k2[i] + k3[i] for i in range(n)]
        k4 = list(map(lambda X: step * X, find_value(equations, curr_value_x + step, var)))
        values_Y.append([])
        # вычисляем новые значения функций
        for i in range(n):
            values_Y[-1].append(curr_value_Y[i] + (1 / 8) * (k1[i] + 3 * k2[i] + 3 * k3[i] + k4[i])) 
        #  задаем новые знаяения аргумента и функций
        curr_value_x = value
        curr_value_Y = [values_Y[-1][i] for i in range(n)]
    return values_x, values_Y  # возвращяем список из значений аргумента и списка из списков значений
                               # функций



def show_graph(values_x, values_Y):
    """
    values_U1 = [exp(-30 * x) for x in values_x]
    Y1 = [values_Y[j][0] for j in range(len(values_Y))]
    m1 = [abs(values_U1[i] - Y1[i]) for i in range(len(values_U1))]
    err = max(m1)   
    print("Погрешность: " + str(err))
    """
    for i in range(n):
        value = [values_Y[j][i] for j in range(len(values_Y))]
        plt.plot(values_x, value, label = "u" + str(i + 1) + "(x)")
    plt.grid()
    plt.legend()
    plt.show()
    

if __name__ == "__main__":
    target, n, equations, start_conditions, segment, step, eps = input_data("input.txt")
    last_change = 0
    left = 0
    right = 1
    # теперь нужно найти отрезок, на котором будет искомое решение.
    while True:
        values_x, values_Y = solve_system_runge_kutta(n, equations, start_conditions, segment, step)
        if (values_Y[-1][0] < target):         
            if  last_change == 0:
                last_change = 1
            if  last_change == -1:
                left = start_conditions[-1] 
                right = start_conditions[-1] + 0.1
                break
            start_conditions[-1] += 0.1

        if (values_Y[-1][0] > target):
            if  last_change == 0:
                last_change = -1
            if  last_change == 1:
                left = start_conditions[-1]  - 0.1
                right = start_conditions[-1]
                break
            start_conditions[-1] -= 0.1
    # теперь есть отрезок [left, right],  чтобы найти решение более точно.
    while True:
        middle = (left + right) / 2
        left_segment = [start_conditions[0], left]
        middle_segment = [start_conditions[0], middle]
        right_segment = [start_conditions[0], right]
        values_x, values_middle = solve_system_runge_kutta(n, equations, middle_segment, segment, step)
        if abs(values_middle[-1][0] - target) < eps:
            start_conditions = [start_conditions[0], middle]
            break
        if (values_middle[-1][0] > target):
            right = middle
        if (values_middle[-1][0] < target):
            left = middle
    x1, Y1 = solve_system_runge_kutta(n, equations, start_conditions, segment, step)
    show_graph(x1, Y1)