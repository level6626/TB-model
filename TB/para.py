"""
Parameter definitions for the model
"""
# time-tick: 10 min
k = 600 / 6

# Diffusion calc constant
Dk = 1

#Chemokine diffusion coefficient /0.1 min
diffC = 0.64 * Dk / 4

#Chemokine degradation coefficient /0.1 min
decayC = 0.001 * Dk

#Prob. of bacteria being killed within M_R
p_k = 8.51 * 0.01

#Prob. T cell kills a macrophage
pT_k = 6.31 * 0.01

#Percentage of B_I being destroyed when killed
P_kill = 50 * 0.01

#Capacity of extracellular bacteria
K_BE = 200

#Intracellular growth rate of bacteria /min
alpha_BI = 0.00021 / (10.0) * Dk

#Extracellular growth rate of bacteria /min
alpha_BE = 0.00015 / (10.0) * Dk

#No. of intracellular bacteria when -> chronically infected MP
N_c = 10

#No. of bacteria causes a MP to burst
K_BI = 20

#No. of T needed to activate MP
N_tact = 4

#Prob. of T recruitment
T_recr = 32 * 0.01

#Prob. of T movement
T_move = 4.97 * 0.01

#Prob. of a T to activate a MP
T_actm = 6 * 0.01

#T lifespan day / 10 min
T_ls = 3 * 24 * 6

#T decay day
T_decay = 10

#Prob. of MP recruitment
M_recr = 2.11 * 0.01

#Resting MP lifespan day / 10 min
M_rls = 100 * 24 * 6

#Activated MP lifespan day / 10 min
M_als = 10 * 24 * 6

#No. of burstings to become necrotic
N_necr = 8

#T cell speed mium/min
T_sp = 10

#Resting MP speed mium/min
M_rsp = 1

#Activated MP speed mium/min
M_asp = 0.025

#Infected MP speed mium/min
M_isp = 0.0007

#Initial No. of resting MP
M_init = 105

#No. of bacteria killed by MP
N_RK = 2

#No. of bacteria killed by MP
N_phag = 10

#Initial No. of extracellular bacteria
BE_init = 16

#Chemokine released by MP
c_I = 5000

#Running time
t_total = round(100 * 24 * 600 / k)