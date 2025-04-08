from os import system

import numpy as np
import math
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['figure.figsize'] = (10.0, 8.0)  # set default size of plots
plt.rcParams['image.interpolation'] = 'nearest'
plt.rcParams['image.cmap'] = 'gray'

#需要用到的气动函数
def pi_ma_equation(ma, k=1.4):
    pi_ma = (1 + (k - 1) / 2 * ma ** 2) ** (-k / (k - 1))
    return pi_ma

def q_ma_equation(ma, k=1.4):
    q_ma = ma * (2 / (k + 1) *(1 + (k - 1) / 2 * ma ** 2)) ** (-(k + 1) / (2 * k - 2))
    return q_ma

def q_lamda_equation(lamda, k=1.4):
    q_lamda = ((k + 1) / 2) ** (1 / (k - 1)) * lamda * (1 - (k - 1) / (k + 1) * lamda ** 2) **(1 / (k - 1))
    return q_lamda

def pi_lamda_equation(lamda, k=1.4):
    pi_lamda = (1 - (k - 1) / (k + 1) * lamda ** 2) ** (k / (k - 1))
    return pi_lamda

def lamda_ma_equation(ma, k=1.4):
    lamda = math.sqrt(((k + 1) / 2 * ma ** 2) / (1 + (k - 1) / 2 * ma ** 2))
    return lamda

def ma_lamda_equation(lamda, k=1.4):
    ma = math.sqrt((2 / (k + 1) * lamda ** 2) / (1 - (k - 1) / (k + 1) * lamda ** 2))
    return ma

def y_lamda_equation(lamda):
    y_lamda = q_lamda_equation(lamda) / pi_lamda_equation(lamda)
    return y_lamda

def t_ma_equation(ma, k=1.4):
    t_ma = 1 / (1 + (k - 1) / 2 * ma **2)
    return t_ma

# 激波角方程
def shock_angle_equation(beta, ma, delta, k=1.4):
    delta_rad = np.radians(delta)
    beta_rad = np.radians(beta)
    term1 = ma ** 2 * np.sin(beta_rad) ** 2 - 1
    term2 = np.tan(beta_rad) * (((k + 1) / 2) * ma ** 2 - ma ** 2 * np.sin(beta_rad) ** 2 + 1)
    return np.tan(delta_rad) - (term1 / term2)

# 斜激波总压恢复系数方程
def sigma_equation(ma, beta, k=1.4):
    beta_rad = np.radians(beta)
    numerator = 1 / (((2 * k) / (k + 1)) * ma ** 2 * np.sin(beta_rad) ** 2 - ((k - 1) / (k + 1)))
    denominator = 1 / ((2 / (k + 1)) / (ma ** 2 * np.sin(beta_rad) ** 2) + ((k - 1) / (k + 1)))
    return numerator ** (1 / (k - 1)) * denominator ** (k / (k - 1))

# 正激波总压恢复系数方程
def sigma_z_equation(ma, k=1.4):
    numerator = 1 / (((2 * k) / (k + 1)) * ma ** 2 - ((k - 1) / (k + 1)))
    denominator = ((2 / (k + 1)) / ma ** 2 + ((k - 1) / (k + 1)))
    return numerator ** (1 / (k - 1)) * denominator ** (-k / (k - 1))

# 激波前后后马赫数关系
def Ma_after_shock_wave_equation(Ma, beta, k=1.4):
    part1 = (Ma ** 2 + (2 / (k - 1))) / ((2 * k / (k - 1)) * Ma ** 2 * np.sin(np.radians(beta)) ** 2 - 1)
    part2 = (Ma ** 2 - Ma ** 2 * np.sin(np.radians(beta)) ** 2) / ((k - 1) / 2 * Ma ** 2 * np.sin(np.radians(beta)) ** 2 + 1)
    Ma_after = math.sqrt(part1 + part2)
    return Ma_after

# 求解激波角β
def shock_angle(Ma, delta):
    if delta >= 15:
        guesses = np.linspace(delta + 10, 70, 5)
    else:
        guesses = np.linspace(delta + 10, 80, 10)
    solutions = []

    for guess in guesses:
        try:
            sol = fsolve(shock_angle_equation, guess, args=(Ma, delta))
            if 0 <= sol[0] <= 90:
                if not any(np.isclose(sol[0], s, atol=1e-2) for s in solutions):
                    solutions.append(sol[0])
        except:
            continue

    if solutions:
        return min(solutions)  # 返回弱解
    else:
        raise ValueError(f"No valid solution found for delta={delta}°")

#两直线交点求解函数
def find_intersection_with_angles(point1, angle1, point2, angle2):

    x1, y1 = point1
    x2, y2 = point2

    # 将角度转换为斜率
    slope1 = math.tan(math.radians(angle1))
    slope2 = math.tan(math.radians(angle2))

    # 检查斜率是否相等（平行或重合）
    if slope1 == slope2:
        if y1 - slope1 * x1 != y2 - slope2 * x2:
            return "Lines are parallel and do not intersect."
        else:
            return "Lines overlap entirely."

    # 计算交点
    x = (slope1 * x1 - y1 - slope2 * x2 + y2) / (slope1 - slope2)
    y = slope1 * (x - x1) + y1

    return x, y

#遍历求解δ1和δ2
def find_best_deltas(Ma, a, max_delta=30):
    best_result = None
    max_sigma_all = -float('inf')
    results = []  # 用于存储调试信息

    # 遍历 delta1 和 delta2
    for delta1 in np.linspace(5, max_delta, 100):  # 更精细的步长
        for delta2 in np.linspace(5, max_delta, 100):
            try:
                beta1 = shock_angle(Ma, delta1)
                Ma2 = Ma_after_shock_wave_equation(Ma, beta1)
                beta2 = shock_angle(Ma2, delta2)
                Ma3 = Ma_after_shock_wave_equation(Ma2, beta2)
                beta3 = shock_angle(Ma3, delta1 + delta2)
                Ma4 = Ma_after_shock_wave_equation(Ma3, beta3)

                sigma1 = sigma_equation(Ma, beta1)
                sigma2 = sigma_equation(Ma2, beta2)
                sigma3 = sigma_equation(Ma3, beta3)
                sigma4 = sigma_z_equation(Ma4)
                sigma_front = sigma1 * sigma2 * sigma3
                sigma_all = sigma1 * sigma2 * sigma3 * sigma4

                # 调试信息
                result_line = f"delta1={delta1:.2f}, delta2={delta2:.2f}, Ma4={Ma4:.4f}, sigma_all={sigma_all:.4f}, sigma_front={sigma_front:.4f}"
                results.append(result_line)

                if Ma4 > 1 and sigma_all > max_sigma_all:
                    max_sigma_all = sigma_all
                    best_result = (delta1, delta2, Ma4, sigma_all, sigma_front)

            except ValueError:
                continue

    return best_result

#初始参数
Ma = 4.0                                                #初始马赫数
a = 5                                                   #坐标度量标准（点D横坐标）
total_pressure = 5529.3 / pi_ma_equation(Ma)            #初始总压
total_temperature = 216.65 / t_ma_equation(Ma)          #初始总温
pi = 3.1415926535                                       #数学常数

#初始变量清零
delta1 = 0.0
delta2 = 0.0

# 找到最佳的 delta1 和 delta2
best_deltas = find_best_deltas(Ma, a)

if best_deltas:
    delta1, delta2, Ma4, sigma_all, sigma_front = best_deltas
    print(f"找到最佳解：delta1 = {delta1:.2f}°, delta2 = {delta2:.2f}°, Ma4 = {Ma4:.4f}, sigma_all = {sigma_all:.4f}, sigma_front={sigma_front:.4f}")
else:
    print("未找到符合条件的 delta1 和 delta2。")

#各种中间过程参数计算
beta1 = shock_angle(Ma, delta1)
Ma2 = Ma_after_shock_wave_equation(Ma, beta1)
beta2 = shock_angle(Ma2, delta2)
Ma3 = Ma_after_shock_wave_equation(Ma2, beta2)
beta3 = shock_angle(Ma3, delta1 + delta2)

#坐标求解中间量
angle_OA = beta1
angle_DA = beta2 + delta1
angle_AB = beta3 - (delta1 + delta2)
angle_DB = delta1 + delta2

point_OA = (0.0, 0.0)
xd = np.tan(np.radians(delta1)) * a
point_DA = a, float(xd)
xa, ya = find_intersection_with_angles(point_OA, angle_OA, point_DA, angle_DA)
point_AB = (xa, ya)
point_DB = a, float(xd)
xb, yb = find_intersection_with_angles(point_DB, angle_DB, point_AB, -1 * angle_AB)

#第二部分计算赋予初始值
delta1, delta2, Ma4, sigma_all, sigma_front = best_deltas

#第二部分总压
p_z4 = total_pressure * sigma_front

#给定初始外界反压
p_e = 220000

#进口面积
A_4 = (3.3493 - 2.7794) ** 2 / 4 * pi

recovery_coefficients = []
out_mach_numbers = []
shock_positions = []
pressure_ratios = []
inner_recovery_coefficients = []  # 内压部分总压恢复系数


gamma_values = np.linspace(2, 15, 50)

for gamma in gamma_values:
    A_e = (gamma * A_4)

    # 进口流量函数
    q_lamda_in = q_lamda_equation(lamda_ma_equation(Ma4))

    # 流量计算
    def flow_equation(A, p_z, T_z, q_lamda_in, K=0.0404):
        m = K * (p_z / math.sqrt(T_z)) * A * q_lamda_in
        return m

    # 非线性方程组
    def flow_equation_e(lamda_e, A_e):
        equation = q_lamda_equation(lamda_e) - (A_4 / A_e) * q_lamda_equation(lamda_ma_equation(Ma4))
        return equation

    # 求解 lamda_e（不考虑正激波时出口lamda数）
    def solve_flow_equation(A_e):
        if gamma <= 2:
            guesses = 1.2, 1.5, 1.8
        elif 2 < gamma <= 3.5:
            guesses = 1.8, 1.9
        else:
            guesses = 2.1, 2.3
        solutions = []

        for guess in guesses:
            try:
                sol = fsolve(flow_equation_e, guess, args=A_e)
                if sol[0] <= math.sqrt(6):  # 确保解在限制范围内
                    solutions.append(sol[0])
            except:
                continue

        solutions.sort()
        return solutions[-1]

    # 考虑正激波后流动方程
    def new_flow_equation(lamda_e, gamma):
        equation = y_lamda_equation(lamda_e) - (p_z4 / p_e) * (1 / gamma) * q_lamda_equation(q_ma_equation(Ma4))
        return equation

    # 求解lamda_e
    def solve_new_flow_equation(gamma):
        guesses = [0.4, 0.8, 1.2, 1.6]  # 初始猜测值的列表
        solutions = []

        for guess in guesses:
            try:
                sol = fsolve(new_flow_equation, guess, args=(gamma,), maxfev=1000, xtol=1e-8)
                # 确保返回的解是唯一的
                if len(solutions) == 0 or not any(np.isclose(sol[0], s, atol=1e-5) for s in solutions):
                    solutions.append(sol[0])  # 保存唯一解
            except Exception as e:
                print(f"错误: {e}")
                continue

        # 返回唯一解
        return solutions[0] if solutions else None

    #print(solve_new_flow_equation(gamma))
    # 不考虑正激波出口lamda数
    lamda_e = solve_flow_equation(A_e)
    #print(f"解得 lamda_e(不考虑正激波) = {lamda_e:.6f}")

    # 进气道内压部分总压恢复系数求解
    sigma_neiya = (1 / gamma) * (
                q_lamda_equation(lamda_ma_equation(Ma4)) / q_lamda_equation(solve_new_flow_equation(gamma)))

    # 不考虑正激波出口压力求解
    p_out = p_z4 * pi_lamda_equation(lamda_e)
    #print(p_out)
    # 求解正激波前马赫数的方程
    def solve_sigma_z_equation(Ma, sigma_neiya, k=1.4):
        numerator = 1 / (((2 * k) / (k + 1)) * Ma ** 2 - ((k - 1) / (k + 1)))
        denominator = ((2 / (k + 1)) / Ma ** 2 + ((k - 1) / (k + 1)))
        return sigma_neiya - numerator ** (1 / (k - 1)) * denominator ** (- k / (k - 1))

    # 求解正激波前马赫数
    def solve_Ma_for_sigma(sigma_neiya, k=1.4):
        # 初始猜测值，可以根据实际情况调整
        guesses = [1.5, 2.0, 2.5, 3.0]
        solutions = []

        for guess in guesses:
            try:
                sol = fsolve(solve_sigma_z_equation, guess, args=(sigma_neiya, k), maxfev=1000, xtol=1e-8)

                # 检查解是否唯一，如果是则添加到 solutions 列表中
                if not any(np.isclose(sol[0], s, atol=1e-6) for s in solutions):
                    solutions.append(sol[0])
            except Exception as e:
                print(f"错误: {e}")
                continue

        # 返回唯一解
        return solutions[0] if solutions else None

    # 正激波前马赫数
    Ma_z_front = solve_Ma_for_sigma(sigma_neiya)
    # 出口 lamda 计算（考虑和不考虑正激波）
    lamda_out = solve_new_flow_equation(gamma)  # 出口 lamda
    Ma_out = ma_lamda_equation(lamda_out)  # 出口马赫数
    sigma_inner = (1 / gamma) * (q_lamda_equation(lamda_ma_equation(Ma4)) / q_lamda_equation(lamda_out))  # 内压部分
    #print(sigma_inner)
    sigma = sigma_inner * sigma_front
    p_out = p_z4 * pi_lamda_equation(lamda_out)
    R = p_out / p_e  # 增压比
    # 激波位置计算
    shock_position = q_lamda_equation(lamda_ma_equation(Ma4)) / q_lamda_equation(lamda_ma_equation(Ma_z_front))

    # 保存结果
    recovery_coefficients.append(sigma)
    out_mach_numbers.append(Ma_out)
    shock_positions.append(shock_position)
    pressure_ratios.append(R)
    inner_recovery_coefficients.append(sigma_inner)

plt.figure(figsize=(8, 6))
plt.plot(pressure_ratios, recovery_coefficients , label="总压恢复系数", color='b')
plt.xlabel("增压比 (R = $p_{out} / p_e$)")
plt.ylabel("总压恢复系数")
plt.title("总压恢复系数随增压比变化")
plt.legend()
plt.grid(True)
plt.show()

# 图 2: 出口马赫数 vs 增压比
plt.figure(figsize=(8, 6))
plt.plot(pressure_ratios, out_mach_numbers, label="出口马赫数", color='r')
plt.xlabel("增压比 (R = $p_{out} / p_e$)")
plt.ylabel("出口马赫数")
plt.title("出口马赫数随增压比变化")
plt.legend()
plt.grid(True)
plt.show()

# 图 3: 正激波位置 vs 增压比
plt.figure(figsize=(8, 6))
plt.plot(pressure_ratios, shock_positions, label="正激波位置", color='g')
plt.xlabel("增压比 (R = $p_{out} / p_e$)")
plt.ylabel("正激波位置")
plt.title("正激波位置随增压比变化")
plt.legend()
plt.grid(True)
plt.show()

system("pause")