include("Functions/print_sameline.jl")
############################
######### Libraries ########
print_sameline("Initialising libraries")
using Plots
using JLD2

print_sameline("Loading functions for transitional_dynamics procedure")
include("Functions/transitional_dynamics_equilibrium.jl")

# global parameters of the model's code
#                   1           2           3       4
#                gen_tol_x, gen_tol_f, distr_tol, val_tol
GLOBAL_PARAMS = [1e-6, 1e-4, 1e-9, 1e-7]#[1e-8, 1e-4, 1e-9, 1e-7]#[1e-8, 1e-4, 1e-12, 1e-9]#[1e-8, 1e-4, 1e-7, 1e-5]#

# global parameters of the approximation objects
#                               1               2               3               4                       5                   6
#                       number_a_nodes, number_u_m_nodes, number_u_w_nodes, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes
GLOBAL_APPROX_PARAMS = [15,3,3,3,6,3]#[49,3,3,3,6,3]#[25,5,5,3,6,3]#

# parameters of the model's economy (Italy)
LAMBDA = 1.513028#1.548387
BETA = 0.917506#0.859906
DELTA = 0.1
GAMMA = 0.16
ETA = 0.3
THETA = 0.54

C_E = 0.0#0.007001#0.036321
RHO_M = 0.710266#0.849951
RHO_W = 0.96
SIGMA_EPS_M = 1.022857#0.942809
RHO_EPS_M_W = 0.361605#0.384691
SIGMA_ZETA = 0.091089#0.151709
P_ALPHA = 0.024445#0.043692
ETA_ALPHA = 5.896035#5.202754
PROB_NODE1_ALPHA = 0.39
MU_M_ALPHA = -4.817675#-3.909137
RHO_ALPHA_M_W = 0.181454#0.095903
SIGMA_ALPHA_W = 0.12
SIGMA_EPS_W = 0.021947#0.285519

# if you change CRRA dynamically, check the code beforehand
CRRA = 1.0 # >0.0

R_ = -0.005338#-0.069128
W_ = 0.230294#0.393906

gen_tol_x = GLOBAL_PARAMS[1]
gen_tol_f = GLOBAL_PARAMS[2]

# interest rate bounds
r_min = -DELTA#-delta
r_max = 1/BETA-1#1/beta-1

# wage bounds #### NOT DONE!!! ########
w_min = 0.18#0.25
w_max = 0.3#0.5

optimal_r = R_#(r_min+r_max)/2#r#
optimal_w = W_#(w_min+w_max)/2#w#

text_output = false
fig_output = false
calc_add_results = false

country = "Italy"#"Brazil"
LAMBDA_GE_DIR = "$(@__DIR__)/Results/Lambda_grid_c_e_0_test/$(country)/"
if Sys.iswindows()
   LAMBDA_GE_DIR = "$(@__DIR__)\\Results\\Lambda_grid_c_e_0_test\\$(country)\\"
end
# Step 2
# ss_star = [res, r, w, approx_object, params]
@load "$(LAMBDA_GE_DIR)SS_lambda_1.51.jld2" SS
lambda_1 = SS[5][1]
ss_1 = copy(SS)
#@time ss_star = steady_state(lambda_star)

@load "$(LAMBDA_GE_DIR)SS_lambda_2.6.jld2" SS
lambda_2 = SS[5][1]
ss_2 = copy(SS)
#@time ss_starstar = steady_state(lambda_2)

@load "$(LAMBDA_GE_DIR)SS_lambda_1.67.jld2" SS
lambda_3 = SS[5][1]
ss_3 = copy(SS)
#@time ss_starstar = steady_state(lambda_2)

MAXITERS = 75#111#50#500#100#
TIME_PERIODS = 100#150#20#50#200#
SMOOTHING = false#true
RUNWAY = 0# Int64(round(TIME_PERIODS/2; digits=0))
lambda_i = 2#=,3=#
lambda_local = lambda_2
ss_local = copy(ss_2)
if lambda_i == 3
    lambda_local = lambda_3
    ss_local = copy(ss_3)
end

LOCAL_DIR = "$(@__DIR__)/Results/Transitional_dynamics_c_e_0_test/$(country)/"
if Sys.iswindows()
   LOCAL_DIR = "$(@__DIR__)\\Results\\Transitional_dynamics_c_e_0_test\\$(country)\\"
end
mkpath(LOCAL_DIR)

# one-time change
LAMBDA_S = ones(TIME_PERIODS).*lambda_local
LAMBDA_S[1] = lambda_1
MODEL_PARAMS = [lambda_1, BETA, DELTA, GAMMA, ETA, THETA, C_E, RHO_M, RHO_W, SIGMA_EPS_M, RHO_EPS_M_W, SIGMA_ZETA, P_ALPHA, ETA_ALPHA, PROB_NODE1_ALPHA, MU_M_ALPHA, RHO_ALPHA_M_W, SIGMA_ALPHA_W, SIGMA_EPS_W, CRRA]

trans_res = open("$(LOCAL_DIR)trans_$(TIME_PERIODS)_results.txt", "w") do F
    #1         2        3            4              5                     6             7          8         9
    #lambda_s, ss_star, ss_starstar, GLOBAL_PARAMS, GLOBAL_APPROX_PARAMS, model_params, file_name, GUESS_RS, maxiters
    #                                 1         2    3         4             5                     6             7  8     9
    @time res = transitional_dynamics(LAMBDA_S, ss_1,ss_local, GLOBAL_PARAMS,GLOBAL_APPROX_PARAMS, MODEL_PARAMS, F, true, MAXITERS)
    (res)
end

p1 = Plots.plot([trans_res[3], ones(TIME_PERIODS).*ss_1[2],ones(TIME_PERIODS).*ss_2[2]], legend=false)
p2 = Plots.plot([trans_res[4], ones(TIME_PERIODS).*ss_1[3],ones(TIME_PERIODS).*ss_2[3]], legend=false)
display(Plots.plot(p1,p2, layout=(1,2)))


ss_2_temp = copy(ss_2)
ss_2 = copy(ss_local)
@save "$(LOCAL_DIR)trans_$(TIME_PERIODS)_results.jld2" trans_res ss_1 ss_2
ss_2 = copy(ss_2_temp)
