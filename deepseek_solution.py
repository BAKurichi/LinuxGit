import math
import numpy as np
from scipy.optimize import minimize

# 原始计算模型封装
def engine_calculation(params):
    bypass_ratio, pi_cl, pi_ch, t_t4 = params
    
    # 固定参数保持不变
    ma = 0.8; h = 11000; t0 = 216.77; p0 = 22700
    sigma_i = 0.99; eta_cl = 0.89; eta_ch = 0.86
    delta1 = 0.01; delta2 = 0.01; eta_b = 0.99; sigma_b = 0.97
    eta_th = 0.90; eta_tl = 0.90; eta_mh = 0.99; eta_ml = 0.99
    sigma_e = 0.99; c_p = 1005; k = 1.4; c_pg = 1244; k_g = 1.33
    h_u = 42900000; R = 287; w_in = 14.5
    
    try:
        # ------ 核心计算流程 ------
        # 0-0截面
        c = math.sqrt(k * R * t0)
        p_t0 = p0 / ((1 + (k-1)/2 * ma**2)**(k/(k-1)))
        t_t0 = t0 * (1 + (k-1)/2 * ma**2)
        
        # 进气道出口
        p2 = p_t0 * sigma_i
        t_t2 = t_t0
        
        # 风扇计算
        p_t22 = p2 * pi_cl
        t_t22 = t_t2 * (1 + (pi_cl**((k-1)/k) -1)/eta_cl)
        l_cl = c_p * (t_t22 - t_t2)
        
        # 高压压气机
        p_t3 = p_t22 * pi_ch
        t_t3 = t_t22 * (1 + (pi_ch**((k-1)/k) -1)/eta_ch)
        l_ch = c_p * (t_t3 - t_t22)
        
        # 燃烧室
        w_3a = w_in * (1 - delta1 - delta2)
        f = (c_pg*t_t4 - c_p*t_t3)/(eta_b*h_u - c_pg*t_t4)
        p_t4_val = p_t3 * sigma_b
        w_4 = w_3a * (1 + f)
        
        # 高压涡轮
        w_4a = w_4 + w_in*delta1
        tau_m1 = (c_p*w_in*delta1*t_t3)/(c_pg*w_4a*t_t4) + w_4/w_4a
        t_t4a = t_t4 * tau_m1
        eta_h = eta_mh * eta_th
        denominator = eta_h * c_pg * t_t4a * (1 + f)
        
        if denominator <= 1e-6:  # 防止除零
            return np.inf, np.inf
            
        pi_h = (1 - l_ch/denominator) ** (-k_g/(k_g-1))
        p_t5 = p_t4_val / pi_h
        t_t5 = t_t4a - t_t4a*eta_th*(1 - 1/pi_h**((k_g-1)/k_g))
        
        # 低压涡轮
        w_4c = w_4a + w_in*delta2
        tau_m2_numerator = (1-delta1-delta2)*(1+f)*c_pg*t_t5 + (delta1+delta2)*c_p*t_t3
        tau_m2_denominator = ((1-delta1-delta2)*(1+f) + delta1+delta2)*c_pg*t_t5
        tau_m2 = tau_m2_numerator / tau_m2_denominator
        t_t4c = t_t5 * tau_m2
        
        eta_l = eta_ml * eta_tl
        denominator_l = eta_l * c_pg * t_t4c * (1 + f)
        if denominator_l <= 1e-6:  # 防止除零
            return np.inf, np.inf
            
        pi_l = (1 - l_cl*bypass_ratio/denominator_l) ** (-k_g/(k_g-1))
        p_t6 = p_t5 / pi_l
        t_t6 = t_t4c - t_t4c*eta_tl*(1 - 1/pi_l**((k_g-1)/k_g))
        
        # 尾喷管
        p_t7 = sigma_e * p_t6
        try:
            ma1 = math.sqrt(max(0, (2/(k_g-1))*((p_t7/p_t0)**((k_g-1)/k_g) -1)))
        except:
            ma1 = 0.01
        t_out = t_t6 / (1 + (k_g-1)/2*ma1**2)
        c_out = math.sqrt(k_g*R*t_out) * ma1
        
        # 外涵道
        w_out = w_in * bypass_ratio
        p_t8 = sigma_e * p_t22
        try:
            ma2 = math.sqrt(max(0, (2/(k_g-1))*((p_t8/p_t0)**((k_g-1)/k_g) -1)))
        except:
            ma2 = 0.01
        t_out1 = t_t22 / (1 + (k_g-1)/2*ma2**2)
        c_out1 = math.sqrt(k_g*R*t_out1) * ma2
        
        # 推力和耗油率
        f_g = w_in*c_out + w_out*c_out1 - (w_in + w_out)*c*ma
        w_f = f * w_in * (1 - delta1 - delta2)
        sfc = 3600 * w_f / f_g if f_g > 100 else np.inf  # 过滤不合理结果
        
        return f_g, sfc
    
    except Exception as e:
        return np.inf, np.inf

# 优化目标函数
def objective(x):
    # 参数边界约束
    x[0] = np.clip(x[0], 3, 8)      # bypass_ratio
    x[1] = np.clip(x[1], 2.5, 4)    # pi_cl
    x[2] = np.clip(x[2], 6, 12)     # pi_ch
    x[3] = np.clip(x[3], 1400, 1600)# t_t4
    
    f_g, sfc = engine_calculation(x)
    
    # 目标权重
    thrust_error = (f_g - 25000)/25000
    sfc_error = (sfc - 0.6)/0.6
    
    # 参数偏离惩罚项
    init_params = [5, 3, 8, 1500]
    param_penalty = 0.1 * sum((x[i]-init_params[i])**2 for i in range(4))
    
    # 综合损失函数
    loss = thrust_error**2 + 10*sfc_error**2 + param_penalty
    
    return loss

# 优化执行
if __name__ == "__main__":
    # 初始猜测值
    x0 = np.array([5.0, 3.0, 8.0, 1500.0])
    
    # 设置优化器参数
    options = {
        'maxiter': 500,
        'xatol': 1e-4,
        'fatol': 1e-5,
        'adaptive': True
    }
    
    # 执行优化
    result = minimize(
        objective,
        x0,
        method='Nelder-Mead',
        options=options
    )
    
    # 提取最佳参数
    optimal_params = np.clip(result.x, [3,2.5,6,1400], [8,4,12,1600])
    final_thrust, final_sfc = engine_calculation(optimal_params)
    
    # 结果输出
    print("\n=== 优化结果 ===")
    print(f"涵道比: {optimal_params[0]:.2f} (初始5.0)")
    print(f"风扇压比: {optimal_params[1]:.2f} (初始3.0)")
    print(f"压气机压比: {optimal_params[2]:.2f} (初始8.0)")
    print(f"涡轮前温度: {optimal_params[3]:.1f} K (初始1500)")
    print(f"推力: {final_thrust:.1f} N (目标25000 N)")
    print(f"耗油率: {final_sfc:.4f} kg/(N·h) (目标0.6000)")
    print(f"优化状态: {result.message}")