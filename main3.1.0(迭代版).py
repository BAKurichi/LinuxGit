import math
from scipy.optimize import minimize

# --- 实用函数 ---
def pi_ma_equation(ma, k):
    return (1 + (k - 1) / 2 * ma ** 2) ** (-k / (k - 1))

def t_ma_equation(ma, k):
    return 1 / (1 + (k - 1) / 2 * ma ** 2)

# --- 发动机模型 ---
def engine_model(x, log=False):
    # 解包设计变量
    bypass_ratio, pi_cl, pi_ch, t_t4 = x

    # --- 固定参数（环境、效率、气体常数等） ---
    ma = 0.8; t0 = 216.77; p0 = 22700
    sigma_i = 0.99; eta_cl = 0.89; eta_ch = 0.86
    delta1 = 0.01; delta2 = 0.01
    eta_b = 0.99; sigma_b = 0.97
    eta_th = 0.90; eta_tl = 0.90
    eta_mh = 0.99; eta_ml = 0.99
    sigma_e = 0.99
    c_p = 1005; k = 1.4; c_pg = 1244; k_g = 1.33
    h_u = 42900000; R = 287
    w_in = 14.5
    w_out = w_in * bypass_ratio

    # --- 各段参数计算 ---
    c = math.sqrt(k * R * t0)
    p_t0 = p0 / pi_ma_equation(ma, k)
    t_t0 = t0 / t_ma_equation(ma, k)
    p2 = p_t0 * sigma_i
    t_t2 = t_t0
    p_t22 = p2 * pi_cl
    t_t22 = t_t2 * (1 + (pi_cl ** ((k - 1) / k) - 1) / eta_cl)
    l_cl = c_p * (t_t22 - t_t2)
    p_t3 = p_t22 * pi_ch
    t_t3 = t_t22 * (1 + (pi_ch ** ((k - 1) / k) - 1) / eta_ch)
    l_ch = c_p * (t_t3 - t_t22)
    w_3a = w_in * (1 - delta1 - delta2)
    f = (c_pg * t_t4 - c_p * t_t3) / (eta_b * h_u - c_pg * t_t4)
    p_t4 = p_t3 * sigma_b
    f_all = 1 + f
    w_4 = w_3a * f_all
    w_4a = w_4 + w_in * delta1
    tau_m1 = c_p * w_in * delta1 * t_t3 / (c_pg * w_4a * t_t4) + w_4 / w_4a
    t_t4a = t_t4 * tau_m1
    eta_h = eta_mh * eta_th
    pi_h = (1 - l_ch / (eta_h * c_pg * t_t4a * f_all)) ** (-(k_g / (k_g - 1)))
    p_t5 = p_t4 / pi_h
    t_t5 = t_t4a - t_t4a * eta_th * (1 - 1 / (pi_h ** ((k_g -1) / k_g)))
    w_4c = w_4a + w_in * delta2
    tau_m2 = ((1 - delta1 - delta2) * f_all * c_pg * t_t5 + delta1 * delta2 * c_p * t_t3) / ((f_all + delta1 + delta2) * c_pg * t_t5)
    t_t4c = t_t5 * tau_m2
    eta_l = eta_ml * eta_tl
    pi_l = (1 - l_cl * bypass_ratio / (eta_l * c_pg * t_t4c * f_all)) ** (-(k_g / (k_g - 1)))
    p_t6 = p_t5 / pi_l
    t_t6 = t_t4c - t_t4c * eta_th * (1 - 1 / (pi_l ** ((k_g -1) / k_g)))
    t_t7 = t_t6
    p_t7 = sigma_e * p_t6
    ma1 = math.sqrt((2 / (k_g - 1)) * ((p_t7 / p_t0) ** ((k_g -1) / k_g) - 1))
    t_out = t_t7 / t_ma_equation(ma1, k_g)
    c_out = math.sqrt(k_g * R * t_out) * ma1
    t_t8 = t_t22
    p_t8 = sigma_e * p_t22
    ma2 = math.sqrt((2 / (k_g - 1)) * ((p_t8 / p_t0) ** ((k_g -1) / k_g) - 1))
    t_out1 = t_t8 / t_ma_equation(ma2, k_g)
    c_out1 = math.sqrt(k_g * R * t_out1) * ma2

    # --- 推力与耗油率 ---
    f_g = w_in * c_out + w_out * c_out1 - (w_in + w_out) * c * ma
    w_f = f * w_in * (1 - delta2 - delta1)
    sfc = 3600 * w_f / f_g

    # --- 日志输出 ---
    if log:
        print(f"迭代参数: 涵道比={bypass_ratio:.3f}, 风扇压比={pi_cl:.3f}, 压气机压比={pi_ch:.3f}, Tt4={t_t4:.1f}")
        print(f"    → 推力 = {f_g:.2f} N, 耗油率 = {sfc:.4f} kg/(N·h)")

    # --- 目标函数 ---
    target_thrust = 25000
    target_sfc = 0.6
    thrust_error = (f_g - target_thrust) / target_thrust  # 推力误差
    sfc_error = (sfc - target_sfc) / target_sfc  # 耗油率误差

    # 目标函数加权求和，使得推力误差和耗油率误差对结果影响均衡
    loss = (thrust_error ** 2) + (sfc_error ** 2)
    return loss

# --- 初始设计点与边界 ---
x0 = [5.0, 3.0, 8.0, 1500.0]  # 初始猜测
bounds = [
    (4.5, 5.5),    # 涵道比
    (2.8, 3.2),    # 风扇压比
    (7.5, 8.5),    # 压气机压比
    (1450, 1550)   # 涡轮前温度
]

# 包装器，用于带日志地运行模型
def objective_with_log(x):
    return engine_model(x, log=True)

# --- 开始优化 ---
result = minimize(objective_with_log, x0, bounds=bounds, method='L-BFGS-B')

# --- 打印最终结果 ---
print("\n=== 优化完成 ===")
print(f"最优参数:")
print(f"  涵道比 = {result.x[0]:.3f}")
print(f"  风扇压比 = {result.x[1]:.3f}")
print(f"  压气机压比 = {result.x[2]:.3f}")
print(f"  涡轮前温度 = {result.x[3]:.2f} K")

# 最终性能评估
print("\n=== 最终性能评估 ===")
engine_model(result.x, log=True)
