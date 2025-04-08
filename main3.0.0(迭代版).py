import math
import random

def pi_ma_equation(ma, k):
    return (1 + (k - 1)/2 * ma**2) ** (-k/(k-1))

def t_ma_equation(ma, k):
    return 1 / (1 + (k-1)/2 * ma**2)

def engine_model(bypass_ratio, pi_cl, pi_ch, t_t4):
    # 环境变量
    ma = 0.8
    h = 11000
    t0 = 216.77
    p0 = 22700
    
    # 压气机参数设定
    sigma_i = 0.99
    eta_cl = 0.89
    eta_ch = 0.86
    delta1 = 0.01
    delta2 = 0.01
    
    # 燃烧室参数
    eta_b = 0.99
    sigma_b = 0.97
    
    # 涡轮参数
    eta_th = 0.90
    eta_tl = 0.90
    eta_mh = 0.99
    eta_ml = 0.99
    
    # 尾喷管
    sigma_e = 0.99
    
    # 物理常数
    c_p = 1005
    k = 1.4
    c_pg = 1244
    k_g = 1.33
    h_u = 42900000
    R = 287
    
    # 流量参数
    w_in = 14.5
    w_out = w_in * bypass_ratio
    
    try:
        # 0-0截面计算
        c = math.sqrt(k * R * t0)
        p_t0 = p0 / pi_ma_equation(ma, 1.4)
        t_t0 = t0 / t_ma_equation(ma, 1.4)
        
        # 进气道出口
        p2 = p_t0 * sigma_i
        t_t2 = t_t0
        
        # 风扇出口
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
        f_all = 1 + f
        w_4 = w_3a * (1 + f)
        
        # 高压涡轮
        w_4a = w_4 + w_in * delta1
        tau_m1 = (c_p * w_in * delta1 * t_t3)/(c_pg * w_4a * t_t4) + w_4/w_4a
        t_t4a = t_t4 * tau_m1
        p_t4a = p_t4_val
        
        eta_h = eta_mh * eta_th
        denominator = eta_h * c_pg * t_t4a * f_all
        if denominator <= 0:
            raise ValueError("Negative denominator in pi_h calculation")
        
        pi_h = (1 - l_ch/denominator) ** (-k_g/(k_g-1))
        p_t5 = p_t4a / pi_h
        t_t5 = t_t4a - t_t4a*eta_th*(1 - 1/(pi_h**((k_g-1)/k_g)))
        
        # 低压涡轮
        w_4c = w_4a + w_in * delta2
        numerator = (1-delta1-delta2)*(1+f)*c_pg*t_t5 + (delta1+delta2)*c_p*t_t3
        denominator = ((1-delta1-delta2)*(1+f) + delta1+delta2)*c_pg*t_t5
        tau_m2 = numerator / denominator
        t_t4c = t_t5 * tau_m2
        p_t4c = p_t5
        
        eta_l = eta_ml * eta_tl
        denominator_l = eta_l * c_pg * t_t4c * f_all
        if denominator_l <= 0:
            raise ValueError("Negative denominator in pi_l calculation")
        
        pi_l = (1 - l_cl*bypass_ratio/denominator_l) ** (-k_g/(k_g-1))
        p_t6 = p_t4c / pi_l
        t_t6 = t_t4c - t_t4c*eta_tl*(1 - 1/(pi_l**((k_g-1)/k_g)))
        
        # 尾喷管
        t_t7 = t_t6
        p_t7 = sigma_e * p_t6
        try:
            ma1 = math.sqrt((2/(k_g-1)) * ((p_t7/p_t0)**((k_g-1)/k_g) -1))
        except:
            ma1 = 0.01  # 防止负值
        
        t_out = t_t7 / t_ma_equation(ma1, k_g)
        c_out = math.sqrt(k_g*R*t_out) * ma1
        
        # 外涵道
        t_t8 = t_t22
        p_t8 = sigma_e * p_t22
        try:
            ma2 = math.sqrt((2/(k_g-1)) * ((p_t8/p_t0)**((k_g-1)/k_g) -1))
        except:
            ma2 = 0.01
        
        t_out1 = t_t8 / t_ma_equation(ma2, k_g)
        c_out1 = math.sqrt(k_g*R*t_out1) * ma2
        
        # 推力和耗油率
        f_g = w_in*c_out + w_out*c_out1 - (w_in + w_out)*c*ma
        w_f = f * w_in * (1 - delta2 - delta1)
        sfc = 3600 * w_f / f_g if f_g != 0 else float('inf')
        
        return f_g, sfc
    
    except Exception as e:
        return float('inf'), float('inf')

# 优化参数配置
params = {
    'bypass_ratio': 8.0,
    'pi_cl': 4.5,
    'pi_ch':14.0,
    't_t4': 1700.0
}

bounds = {
    'bypass_ratio': (5.0, 9.0),
    'pi_cl': (3.0, 5.0),
    'pi_ch': (6.0, 15.0),
    't_t4': (1400.0, 1700.0)
}

# 优化控制参数
target_thrust = 25000.0
target_sfc = 0.6
learning_rate = 0.1
max_iterations = 1000
perturb = 0.01
tolerance = 1e-5

best_loss = float('inf')
best_params = params.copy()
history = []

for iter in range(max_iterations):
    try:
        # 当前参数计算
        thrust, sfc = engine_model(**params)
        
        # 损失计算（双目标加权）
        thrust_error = (thrust - target_thrust)/target_thrust
        sfc_error = (sfc - target_sfc)/target_sfc
        loss = (thrust_error**2) + 8*(sfc_error**2)  # 耗油率权重更大
        
        # 保存最佳参数
        if loss < best_loss:
            best_loss = loss
            best_params = params.copy()
            stagnation = 0
        else:
            stagnation += 1
        
        # 收敛检查
        if stagnation > 20 or loss < tolerance:
            break
        
        # 有限差分法计算梯度
        gradients = {}
        for key in params:
            # 正向扰动
            params_plus = params.copy()
            params_plus[key] += perturb
            for k in bounds:  # 边界约束
                params_plus[k] = max(bounds[k][0], min(params_plus[k], bounds[k][1]))
            
            # 计算扰动后的损失
            t_p, s_p = engine_model(**params_plus)
            l_p = ((t_p-target_thrust)/target_thrust)**2 + 10*((s_p-target_sfc)/target_sfc)**2
            
            # 梯度计算
            gradients[key] = (l_p - loss) / perturb
        
        # 参数更新（带动量项）
        momentum = 0.3
        for key in params:
            delta = -learning_rate * gradients[key] + momentum * (params[key] - best_params[key])
            new_val = params[key] + delta
            
            # 应用边界约束
            new_val = max(bounds[key][0], min(new_val, bounds[key][1]))
            params[key] = new_val
        
        # 学习率衰减
        if iter % 50 == 0:
            learning_rate *= 0.95
        
        # 打印进度
        if iter % 10 == 0:
            print(f"Iter {iter:3d} | Thrust: {thrust:8.1f} | SFC: {sfc:.4f} | Loss: {loss:.4e}")
    
    except Exception as e:
        print(f"Error at iteration {iter}: {str(e)}")
        params = best_params.copy()
        learning_rate *= 0.5

# 最终结果输出
final_thrust, final_sfc = engine_model(**best_params)
print("\n=== 优化结果 ===")
print(f"涵道比: {best_params['bypass_ratio']:.2f} (5.0-8.0)")
print(f"风扇压比: {best_params['pi_cl']:.2f} (3.0-4.0)")
print(f"压气机压比: {best_params['pi_ch']:.2f} (6.0-12.0)")
print(f"涡轮前温度: {best_params['t_t4']:.1f} K (1400-1600)")
print(f"推力: {final_thrust:.1f} N (目标25000 N)")
print(f"耗油率: {final_sfc:.4f} kg/(N·h) (目标0.6000)")