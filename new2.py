import math
import numpy as np
from multiprocessing import Pool, cpu_count
import itertools
import os, psutil
from tqdm import tqdm

# --- 基本气动函数 ---
def pi_ma_equation(ma, k):
    return (1 + (k - 1) / 2 * ma ** 2) ** (-k / (k - 1))
def t_ma_equation(ma, k):
    return 1 / (1 + (k - 1) / 2 * ma **2)

# --- 发动机模型 ---
def engine_model(x):
    bypass_ratio, pi_cl, pi_ch, t_t4 = x
    
    #环境变量和各部件基本设计参数
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
    c = math.sqrt(k * R * t0)              #声速
    p_t0 = p0 / pi_ma_equation(ma,1.4)     #截面总压
    t_t0 = t0 / t_ma_equation(ma,1.4)      #截面总温

    #进气道出口的总压和总温
    p2 = p_t0 * sigma_i #总压
    t_t2 = t_t0         #总温

    #风扇出口参数
    p_t22 = p2 * pi_cl                                              #出口总压
    t_t22 = t_t2 * (1 + (pi_cl ** ((k - 1) / k) - 1) / eta_cl)      #出口总温

    #风扇每千克空气消耗的功
    l_cl = c_p * (t_t22 - t_t2)

    #高压压气机出口总压和总温
    p_t3 = p_t22 * pi_ch
    t_t3 = t_t22 * (1 + (pi_ch ** ((k - 1) / k) - 1) / eta_ch)

    #压气机每千克空气所消耗的功
    l_ch = c_p * (t_t3 - t_t22)

    #主燃烧室出口参数
    w_3a = w_in * (1 - delta1 - delta2)                            #主燃烧室进口流量
    f = (c_pg * t_t4 - c_p * t_t3) / (eta_b * h_u - c_pg * t_t4)   #燃烧室油气比
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
    p_t5 = p_t4a / pi_h #高压涡轮出口总压
    t_t5 = t_t4a - t_t4a * eta_th * (1 - 1 / (pi_h ** ((k_g -1) / k_g)))#高压涡轮出口总温

    #低压涡轮进口参数
    w_4c = w_4a + w_in * delta2
    tau_m2 = ((1 - delta1 - delta2) * (1 + f) * c_pg * t_t5 + delta1 * delta2 * c_p * t_t3) / (((1 - delta1 - delta2) * (1 + f) + delta1 + delta2) * c_pg * t_t5)
    t_t4c = t_t5 * tau_m2
    p_t4c = p_t5

    #低压涡轮的出口参数
    eta_l = eta_ml * eta_tl
    pi_l = (1 - l_cl * (1 + bypass_ratio) /(eta_l * c_pg * t_t4c * f_all)) ** (-(k_g / (k_g - 1)))
    p_t6 = p_t4c / pi_l #低压涡轮出口总压
    t_t6 = t_t4c - t_t4c * eta_th * (1 - 1 / (pi_l ** ((k_g -1) / k_g)))#低压涡轮出口总温

    #尾喷管
    t_t7 = t_t6             #出口总温
    p_t7 = sigma_e * p_t6   #出口总压
    ma1 = math.sqrt((2 / (k_g - 1)) * ((p_t7 / p_t0) ** ((k_g -1) / k_g) - 1))
    t_out = t_t7 / t_ma_equation(ma1,1.33) #出口静温
    c_out = math.sqrt(k_g * R * t_out) * ma1 #出口速度

    #外涵道尾喷s
    t_t8 = t_t22
    p_t8 = sigma_e * p_t22
    ma2 = math.sqrt((2 / (k - 1)) * ((p_t8 / p_t0) ** ((k -1) / k) - 1))
    t_out1 = t_t8 / t_ma_equation(ma2,1.4) #出口静温
    c_out1 = math.sqrt(k * R* t_out1) * ma2 #出口速度

    #推力
    f_g = w_in * c_out + w_out * c_out1 - (w_in + w_out) * c * 0.8
    w_f = f * w_in * (1 - delta2 - delta1)
    sfc = 3600 * w_f / f_g
    
    #函数返回值
    return f_g, sfc

# --- 单点评估函数 ---
def evaluate_combo(params):
    target_thrust = 25000
    target_sfc = 0.06
    try:
        thrust, sfc = engine_model(params)
        error = (0.01*(thrust - target_thrust)) ** 2 + (1000 * (sfc - target_sfc)) ** 2 #给耗油率权重更高，因为耗油率的误差对总误差的影响太小，容易造成计算结果偏差较大但是总误差较小的情况
        return (error, params, thrust, sfc)
    except Exception:
        return (float('inf'), params, None, None)

# --- 分块生成器 ---
def generate_chunks(iterable, chunk_size):
    iterable = iter(iterable)
    while True:
        chunk = list(itertools.islice(iterable, chunk_size))
        if not chunk:
            break
        yield chunk

# --- 搜索空间定义 ---
bypass_range = np.arange(3.8, 5.5 + 0.01, 0.01)
pi_cl_range = np.arange(2.5, 3.4 + 0.01, 0.01)
pi_ch_range = np.arange(7.0, 8.9 + 0.01, 0.01)
t_t4_range = np.arange(1400, 1600 + 1, 1)

all_combos = itertools.product(bypass_range, pi_cl_range, pi_ch_range, t_t4_range)
total_combos = len(bypass_range) * len(pi_cl_range) * len(pi_ch_range) * len(t_t4_range)

# --- 暴力搜索 + 并行 + 分块 + 进度条 ---
if __name__ == '__main__':
    chunk_size = 3000000
    best_error = float('inf')
    best_combo = None
    best_result = None

    processed = 0
    with tqdm(total=total_combos, desc="搜索进度") as pbar:
        for chunk in generate_chunks(all_combos, chunk_size):
            with Pool(processes=cpu_count()) as pool:
                results = pool.map(evaluate_combo, chunk)

            for error, params, thrust, sfc in results:
                if error < best_error:
                    best_error = error
                    best_combo = params
                    best_result = (thrust, sfc)

            processed += len(chunk)
            pbar.update(len(chunk))

            # 打印当前内存使用
            mem = psutil.Process(os.getpid()).memory_info().rss / (1024 ** 3)
            print(f"当前内存占用: {mem:.2f} GB")

    print("\n--- 搜索完成 ---")
    print("最优参数组合:", best_combo)
    print("对应输出: 推力 = {:.2f} N, SFC = {:.4f} kg/(N·h)".format(*best_result))
    print("总误差（平方差和）:", best_error)
