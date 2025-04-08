import math

#基本气动函数
def pi_ma_equation(ma, k):
    pi_ma = (1 + (k - 1) / 2 * ma ** 2) ** (-k / (k - 1))
    return pi_ma

def t_ma_equation(ma, k):
    t_ma = 1 / (1 + (k - 1) / 2 * ma **2)
    return t_ma


#环境变量
ma = 0.8;h = 11000;t0 = 216.77; p0 = 22700

#压气机参数设定
bypass_ratio = 3.94      #涵道比
pi_cl = 2.55             #风扇的压比
pi_ch = 7.25             #压气机的压比
sigma_i = 0.99          #进气道总压恢复系数
eta_cl = 0.89           #风扇绝热效率
eta_ch = 0.86           #高压压气机绝热效率
delta1 = 0.01           #高压涡轮冷却气指数
delta2 = 0.01           #低压涡轮冷却气指数

#燃烧室设定参数
eta_b = 0.99            #主燃烧室燃烧效率
sigma_b = 0.97          #主燃烧室总压恢复系数
t_t4 = 1459             #涡轮前温度

#高压和低压涡轮设定参数
eta_th = 0.90           #高压涡轮效率
eta_tl = 0.90           #低压涡轮效率
eta_mh = 0.99           #高压轴机械效率
eta_ml = 0.99           #低压轴机械效率

#尾喷管设定参数
sigma_e = 0.99          #尾喷管总压恢复系数

#基本物理参数
c_p = 1005              #空气的定压比热容
k = 1.4                 #空气比热比
c_pg = 1244             #燃气的定压比热容
k_g = 1.33              #燃气比热比
h_u = 42900000          #燃料热值
R = 287

w_in = 14.5                     #内涵道进口流量
w_out = w_in * bypass_ratio     #外涵道进口流量
c_fg = 1.0

#计算0-0截面的总温和总压
c = math.sqrt(k * R * t0)               #声速
p_t0 = p0 / pi_ma_equation(ma,1.4)      #进气道进口截面总压
print("进气道进口截面总压p_t0 = ", p_t0)
t_t0 = t0 / t_ma_equation(ma,1.4)       #进气道出口截面总温
print("进气道进口截面总温t_t0 = ", t_t0)

#进气道出口的总压和总温
p2 = p_t0 * sigma_i                     #进气道出口截面总压
t_t2 = t_t0                             #进气道出口截面总温

#风扇(低压压气机)进口截面参数
w_fin = w_in + w_out                    #风扇(低压压气机)进口截面流量
p_fin = p2                              #风扇(低压压气机)进口截面总压
t_fin = t_t2                            #风扇(低压压气机)进口截面总温

#风扇(低压压气机)出口截面参数
w_fout = w_fin                                                  #风扇(低压压气机)出口截面流量
print("风扇(低压压气机)出口截面流量 = ", w_fout)
p_t22 = p2 * pi_cl                                              #风扇(低压压气机)出口截面总压
print("风扇(低压压气机)出口截面总压p_t22 = ", p_t22)
t_t22 = t_t2 * (1 + (pi_cl ** ((k - 1) / k) - 1) / eta_cl)      #风扇(低压压气机)出口截面总温
print("风扇(低压压气机)出口截面总温t_t22 = ", t_t22)

#风扇每千克空气消耗的功
l_cl = c_p * (t_t22 - t_t2)
print("风扇每千克空气消耗的功l_cl = ", l_cl)

#高压压气机进口截面参数
w_hin = w_in                                                    #高压压气机进口截面流量
p_hin = p_t22                                                   #高压压气机进口截面总压
t_hin = t_t22                                                   #高压压气机进口截面总温

#高压压气机出口截面参数
w_hout = w_hin                                                  #高压压气机出口截面流量
p_t3 = p_t22 * pi_ch                                            #高压压气机出口截面总压
t_t3 = t_t22 * (1 + (pi_ch ** ((k - 1) / k) - 1) / eta_ch)      #高压压气机出口截面总温
print("高压压气机出口总温t_t3 = ", t_t3)

#压气机每千克空气所消耗的功
l_ch = c_p * (t_t3 - t_t22)
print("压气机每千克空气所消耗的功l_ch = ", l_ch)

#主燃烧室进口截面参数
w_3a = w_in * (1 - delta1 - delta2)                             #主燃烧室进口截面流量
p_rin = p_t3

#主燃烧室出口参数
f = (c_pg * t_t4 - c_p * t_t3) / (eta_b * h_u - c_pg * t_t4)    #燃烧室油气比
print("燃烧室油气比f = ", f)
p_t4 = p_t3 * sigma_b                                           #燃烧室出口截面总压
f_all = 1 + f                                                   #内涵道总流量系数
w_4 = w_3a *(1 + f)                                             #燃烧室出口截面流量

#高压涡轮进口参数
w_4a = w_4 + w_in * delta1                                                  #高压涡轮进口截面流量
tau_m1 = c_p * w_in * delta1 * t_t3 / (c_pg * w_4a * t_t4) + w_4 / w_4a     #高压涡轮进口截面处冷却气与高温燃气混合前后总温比
t_t4a = t_t4 * tau_m1                                                       #高压涡轮进口截面总温
p_t4a = p_t4                                                                #高压涡轮进口截面总压

#高压涡轮的出口参数
eta_h = eta_mh * eta_th
pi_h = (1 - l_ch / (eta_h * c_pg * t_t4a * f_all)) ** (-(k_g / (k_g - 1)))
print(1 - l_ch/ (eta_h * c_pg * t_t4a * f_all))
print("落压比pi_h = ", pi_h)

p_t5 = p_t4a / pi_h                                                         #高压涡轮出口截面总压
t_t5 = t_t4a - t_t4a * eta_th * (1 - 1 / (pi_h ** ((k_g -1) / k_g)))        #高压涡轮出口截面总温
print("高压涡轮出口总温t_t5 = ", t_t5)

#低压涡轮进口参数
w_4c = w_4a + w_in * delta2                             #低压涡轮进口截面处流量
tau_m2 = ((1 - delta1 - delta2) * (1 + f) * c_pg * t_t5 + delta1 * delta2 * c_p * t_t3) / (((1 - delta1 - delta2) * (1 + f) + delta1 + delta2) * c_pg * t_t5)                                #低压涡轮进口处冷却气与高压涡轮出口气体混合前后总温比
t_t4c = t_t5 * tau_m2                                   #低压涡轮进口截面总温
p_t4c = p_t5                                            #低压涡轮进口截面总压



#低压涡轮的出口参数
eta_l = eta_ml * eta_tl                                                                         #低压涡轮总效率
pi_l = (1 - l_cl * bypass_ratio /(eta_l * c_pg * t_t4c * f_all)) ** (-(k_g / (k_g - 1)))        #低压涡轮落压比
print("落压比pi_l = ", pi_l)
p_t6 = p_t4c / pi_l                                                                             #低压涡轮出口截面总压
print("低压涡轮出口总压p_t6 = ", p_t6)
t_t6 = t_t4c - t_t4c * eta_th * (1 - 1 / (pi_l ** ((k_g - 1) / k_g)))                           #低压涡轮出口截面总温
print("低压涡轮出口总温t_t6 = ", t_t6)
w_5c = w_4c                                                                                     #低压涡轮出口截面流量 

#尾喷管出口截面参数(假定尾喷管完全膨胀)
t_t7 = t_t6                                     #尾喷管进出口截面总温
p_t7 = sigma_e * p_t6                           #尾喷管进口截面总压
print("尾喷管出口总压p_t7 = ",p_t7)
ma1 = math.sqrt((2 / (k_g - 1)) * ((p_t7 / p_t0) ** ((k_g - 1) / k_g) - 1))
print("尾喷管出口马赫数ma1 = ", ma1)
t_out = t_t7 / t_ma_equation(ma1, k_g)          #尾喷管出口截面静温
c_out = math.sqrt(k_g * R * t_out) * ma1        #尾喷管出口截面速度

#外涵道出口截面参数
t_t8 = t_t22
p_t8 = sigma_e * p_t22
ma2 = math.sqrt((2 / (k_g - 1)) * ((p_t8 / p_t0) ** ((k_g -1) / k_g) - 1))
print("外涵道排气出口马赫数ma2 = ", ma2)
t_out1 = t_t8 / t_ma_equation(ma2,1.4)              #外涵道出口截面静温
c_out1 = math.sqrt(k * R* t_out1) * ma2             #外涵道出口截面速度

#推力
f_g = w_in * c_out + w_out * c_out1 - (w_in + w_out) * c * 0.8
w_f = f * w_in * (1 - delta2 - delta1)
sfc = 3600 * w_f / f_g


print("推力f_g = ", f_g, "N")
print("耗油率sfc = ", sfc, "kg/N·h")

