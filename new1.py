import math
from scipy.optimize import minimize

# --- 实用函数 ---
def pi_ma_equation(ma, k):
    return (1 + (k - 1) / 2 * ma ** 2) ** (-k / (k - 1))

def t_ma_equation(ma, k):
    return 1 / (1 + (k - 1) / 2 * ma ** 2)

# --- 发动机模型 ---
def engine_model(x, log=False, iteration_count=None):
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

    #计算0-0截面的总温和总压
    c = math.sqrt(k * R * t0)          #声速
    p_t0 = p0 / pi_ma_equation(ma,1.4)     #截面总压
    #print("总压p_t0 = ", p_t0)
    t_t0 = t0 / t_ma_equation(ma,1.4)      #截面总温
    #print("总温t_t0 = ", t_t0)

    #进气道出口的总压和总温
    p2 = p_t0 * sigma_i #总压
    t_t2 = t_t0         #总温

    #风扇出口参数
    p_t22 = p2 * pi_cl                                              #出口总压
    #print("风扇出口总压p_t22 = ", p_t22)
    t_t22 = t_t2 * (1 + (pi_cl ** ((k - 1) / k) - 1) / eta_cl)     #出口总温
    #print("风扇出口总温t_t22 = ", t_t22)

    #风扇每千克空气消耗的功
    l_cl = c_p * (t_t22 - t_t2)
    #print("风扇每千克空气消耗的功l_cl = ", l_cl)

    #高压压气机出口总压和总温
    p_t3 = p_t22 * pi_ch
    t_t3 = t_t22 * (1 + (pi_ch ** ((k - 1) / k) - 1) / eta_ch)
    #print("高压压气机出口总温t_t3 = ", t_t3)

    #压气机每千克空气所消耗的功
    l_ch = c_p * (t_t3 - t_t22)
    #print("压气机每千克空气所消耗的功l_ch = ", l_ch)

    #主燃烧室出口参数
    w_3a = w_in * (1 - delta1 - delta2)                            #主燃烧室进口流量
    f = (c_pg * t_t4 - c_p * t_t3) / (eta_b * h_u - c_pg * t_t4)   #燃烧室油气比
    #print("燃烧室油气比f = ", f)
    p_t4 = p_t3 * sigma_b                                          #燃烧室出口总压
    f_all = 1 + f                                                  #内涵道总流量系数
    w_4 = w_3a *(1 + f)

    #高压涡轮进口参数
    w_4a = w_4 + w_in * delta1
    tau_m1 = c_p * w_in * delta1 * t_t3 / (c_pg * w_4a * t_t4) + w_4 / w_4a
    t_t4a = t_t4 * tau_m1
    p_t4a = p_t4

    #高压涡轮的出口参数
    eta_h = eta_mh * eta_th
    pi_h = (1 - l_ch / (eta_h * c_pg * t_t4a * f_all)) ** (-(k_g / (k_g - 1)))
    #print(1 - l_ch/ (eta_h * c_pg * t_t4a * f_all))
    #print("落压比pi_h = ", pi_h)

    p_t5 = p_t4a / pi_h #高压涡轮出口总压
    t_t5 = t_t4a - t_t4a * eta_th * (1 - 1 / (pi_h ** ((k_g -1) / k_g)))#高压涡轮出口总温
    #print("高压涡轮出口总温t_t5 = ", t_t5)

    #低压涡轮进口参数
    w_4c = w_4a + w_in * delta2
    tau_m2 = ((1 - delta1 - delta2) * (1 + f) * c_pg * t_t5 + delta1 * delta2 * c_p * t_t3) / (((1 - delta1 - delta2) * (1 + f) + delta1 + delta2) * c_pg * t_t5)
    t_t4c = t_t5 * tau_m2
    p_t4c = p_t5



    #低压涡轮的出口参数
    eta_l = eta_ml * eta_tl
    pi_l = (1 - l_cl * bypass_ratio /(eta_l * c_pg * t_t4c * f_all)) ** (-(k_g / (k_g - 1)))
    #print("落压比pi_l = ", pi_l)
    p_t6 = p_t4c / pi_l #低压涡轮出口总压
    #print("低压涡轮出口总压p_t6 = ", p_t6)
    t_t6 = t_t4c - t_t4c * eta_th * (1 - 1 / (pi_l ** ((k_g -1) / k_g)))#低压涡轮出口总温
    #print("低压涡轮出口总温t_t6 = ", t_t6)

    #尾喷管
    t_t7 = t_t6             #出口总温
    p_t7 = sigma_e * p_t6   #出口总压
    #print("尾喷管出口总压p_t7 = ",p_t7)
    #print((2 / (k_g - 1)) * ((p_t7 / p_t0) ** ((k_g -1) / k_g) - 1))
    ma1 = math.sqrt((2 / (k_g - 1)) * ((p_t7 / p_t0) ** ((k_g -1) / k_g) - 1))
    #print("尾喷管出口马赫数ma1 = ", ma1)
    t_out = t_t7 / t_ma_equation(ma1,1.33) #出口静温
    c_out = math.sqrt(k_g * R * t_out) * ma1 #出口速度

    #外涵道尾喷s
    t_t8 = t_t22
    p_t8 = sigma_e * p_t22

    ma2 = math.sqrt((2 / (k_g - 1)) * ((p_t8 / p_t0) ** ((k_g -1) / k_g) - 1))
    #print("外涵道排气出口马赫数ma2 = ", ma2)
    t_out1 = t_t8 / t_ma_equation(ma2,1.4) #出口静温
    c_out1 = math.sqrt(k * R* t_out1) * ma2 #出口速度

    #推力
    f_g = w_in * c_out + w_out * c_out1 - (w_in + w_out) * c * 0.8
    w_f = f * w_in * (1 - delta2 - delta1)
    sfc = 3600 * w_f / f_g  

    # --- 日志输出 ---
    if log:
        print(f"迭代次数: {iteration_count}, 涵道比={bypass_ratio:.4f}, 风扇压比={pi_cl:.4f}, 压气机压比={pi_ch:.4f}, Tt4={t_t4:.1f}")
        print(f"    → 推力 = {f_g:.4f} N, 耗油率 = {sfc:.4f} kg/(N·h)")

    # --- 目标函数 ---
    target_thrust = 25000
    target_sfc = 0.6
    thrust_error = (f_g - target_thrust) / target_thrust  # 推力误差
    sfc_error = (sfc - target_sfc) / target_sfc  # 耗油率误差

    # 目标函数加权求和，使得推力误差和耗油率误差对结果影响均衡
    loss = (thrust_error ** 2) + (sfc_error ** 2)
    return loss

# --- 初始设计点与边界 ---
x0 = [5, 3, 8, 1600]  # 初始猜测
bounds = [
    (3.8, 5.5),    # 涵道比
    (2.5, 3.4),    # 风扇压比
    (7.0, 8.9),    # 压气机压比
    (1400, 1600)   # 涡轮前温度
]

# 包装器，用于带日志地运行模型
def objective_with_log(x, iteration_count):
    return engine_model(x, log=True, iteration_count=iteration_count)

# --- 强制进行2000次迭代的优化 ---
max_iterations = 2000
iteration_count = 0

# 定义优化函数
def optimization_iteration(x):
    global iteration_count
    iteration_count += 1
    return engine_model(x, log=True, iteration_count=iteration_count)

# 使用L-BFGS-B算法进行优化，确保至少进行2000次迭代
result = minimize(optimization_iteration, x0, bounds=bounds, method='L-BFGS-B', options={'maxiter': max_iterations})

# --- 打印最终结果 ---
print("\n=== 优化完成 ===")
print(f"最优参数:")
print(f"  涵道比 = {result.x[0]:.3f}")
print(f"  风扇压比 = {result.x[1]:.3f}")
print(f"  压气机压比 = {result.x[2]:.3f}")
print(f"  涡轮前温度 = {result.x[3]:.2f} K")

# 最终性能评估
print("\n=== 最终性能评估 ===")
engine_model(result.x, log=True, iteration_count=iteration_count)
