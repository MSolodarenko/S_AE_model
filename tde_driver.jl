include("Functions/print_sameline.jl")
############################
######### Libraries ########
print_sameline("Initialising libraries")
using Plots
using JLD2

print_sameline("Loading functions for transitional_dynamics procedure")
include("Functions/transitional_dynamics_equilibrium.jl")

LOCAL_DIR = "$(@__DIR__)/Results/Stationary/GE_lambda/General/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Stationary\\GE_lambda\\General\\"
end
@load "$(LOCAL_DIR)SSS.jld2" SSS

# global parameters of the model's code
#                   1           2           3       4
#                gen_tol_x, gen_tol_f, distr_tol, val_tol
GLOBAL_PARAMS = SSS[10][1][46]
# global parameters of the approximation objects
#                               1               2               3               4                       5                   6
#                       number_a_nodes, number_u_m_nodes, number_u_w_nodes, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes
GLOBAL_APPROX_PARAMS = SSS[10][1][47]

# parameters of the model's economy (Italy)
MODEL_PARAMS_INIT = SSS[10][1][48]
LAMBDA =            MODEL_PARAMS_INIT[1]
BETA =              MODEL_PARAMS_INIT[2]
DELTA =             MODEL_PARAMS_INIT[3]
GAMMA =             MODEL_PARAMS_INIT[4]
ETA =               MODEL_PARAMS_INIT[5]
THETA =             MODEL_PARAMS_INIT[6]

C_E =               MODEL_PARAMS_INIT[7]
RHO_M =             MODEL_PARAMS_INIT[8]
RHO_W =             MODEL_PARAMS_INIT[9]
SIGMA_EPS_M =       MODEL_PARAMS_INIT[10]
RHO_EPS_M_W =       MODEL_PARAMS_INIT[11]
SIGMA_ZETA =        MODEL_PARAMS_INIT[12]
P_ALPHA =           MODEL_PARAMS_INIT[13]
ETA_ALPHA =         MODEL_PARAMS_INIT[14]
PROB_NODE1_ALPHA =  MODEL_PARAMS_INIT[15]
MU_M_ALPHA =        MODEL_PARAMS_INIT[16]
RHO_ALPHA_M_W =     MODEL_PARAMS_INIT[17]
SIGMA_ALPHA_W =     MODEL_PARAMS_INIT[18]
SIGMA_EPS_W =       MODEL_PARAMS_INIT[19]

CRRA =              MODEL_PARAMS_INIT[20]

R_ =                SSS[10][2]
W_ =                SSS[10][3]

gen_tol_x = GLOBAL_PARAMS[1]
gen_tol_f = GLOBAL_PARAMS[2]

# interest rate bounds
r_min = -DELTA#-delta
r_max = 1/BETA-1#1/beta-1

# wage bounds #### NOT DONE!!! ########
W_MIN = 0.01
W_MAX = 0.47

optimal_r = R_#(r_min+r_max)/2#r#
optimal_w = W_#(w_min+w_max)/2#w#

text_output = false
fig_output = false
calc_add_results = false

LAMBDA_GE_DIR = "$(@__DIR__)/Results/Stationary/GE_lambda/"
if Sys.iswindows()
   LAMBDA_GE_DIR = "$(@__DIR__)\\Results\\Stationary\\GE_lambda\\"
end
# Step 2
# ss_star = [res, r, w, approx_object, params]
@load "$(LAMBDA_GE_DIR)SS_lambda_1.67.jld2" SS
lambda_1 = SS[5][1]
ss_1 = copy(SS)
#@time ss_star = steady_state(lambda_star)

#@load "$(LAMBDA_GE_DIR)SS_lambda_2.6.jld2" SS
@load "$(LAMBDA_GE_DIR)SS_lambda_10.0.jld2" SS
lambda_2 = SS[5][1]
ss_2 = copy(SS)
#@time ss_starstar = steady_state(lambda_2)
#=
@load "$(LAMBDA_GE_DIR)SS_lambda_1.67.jld2" SS
lambda_3 = SS[5][1]
ss_3 = copy(SS)
#@time ss_starstar = steady_state(lambda_2)
=#
MAXITERS = 100#50#75#111#500#100#
############ change to 100
TIME_PERIODS = 50#
############
SMOOTHING = false#true
RUNWAY = 0# Int64(round(TIME_PERIODS/2; digits=0))
lambda_i = 2#=,3=#
lambda_local = lambda_2
ss_local = copy(ss_2)
if lambda_i == 3
    lambda_local = lambda_3
    ss_local = copy(ss_3)
end

LOCAL_DIR = "$(@__DIR__)/Results/Transitionary/TDE/"
if Sys.iswindows()
   LOCAL_DIR = "$(@__DIR__)\\Results\\Transitionary\\TDE\\"
end
mkpath(LOCAL_DIR)

# one-time change
LAMBDA_S = ones(TIME_PERIODS).*lambda_local
LAMBDA_S[1] = lambda_1
MODEL_PARAMS = [lambda_1, BETA, DELTA, GAMMA, ETA, THETA, C_E, RHO_M, RHO_W, SIGMA_EPS_M, RHO_EPS_M_W, SIGMA_ZETA, P_ALPHA, ETA_ALPHA, PROB_NODE1_ALPHA, MU_M_ALPHA, RHO_ALPHA_M_W, SIGMA_ALPHA_W, SIGMA_EPS_W, CRRA]

trans_res = open("$(LOCAL_DIR)trans_$(TIME_PERIODS)_results.txt", "w") do F
    #1         2        3            4              5                     6             7          8         9
    #lambda_s, ss_star, ss_starstar, GLOBAL_PARAMS, GLOBAL_APPROX_PARAMS, model_params, file_name, GUESS_RS, maxiters
    #                                 1         2    3         4             5                     6             7  8     9         10
    @time res = transitional_dynamics(LAMBDA_S, ss_1,ss_local, GLOBAL_PARAMS,GLOBAL_APPROX_PARAMS, MODEL_PARAMS, F, true, MAXITERS, LOCAL_DIR)
    (res)
end

p1 = Plots.plot([trans_res[3], ones(TIME_PERIODS).*ss_1[2],ones(TIME_PERIODS).*ss_2[2]], legend=false)
p2 = Plots.plot([trans_res[4], ones(TIME_PERIODS).*ss_1[3],ones(TIME_PERIODS).*ss_2[3]], legend=false)
display(Plots.plot(p1,p2, layout=(1,2)))


ss_2_temp = copy(ss_2)
ss_2 = copy(ss_local)
@save "$(LOCAL_DIR)trans_$(TIME_PERIODS)_results.jld2" trans_res ss_1 ss_2
ss_2 = copy(ss_2_temp)
