include("Functions/print_sameline.jl")
############################
######### Libraries ########
print_sameline("Initialising libraries")
using Plots
using JLD2

print_sameline("Loading functions for steady_state procedure")
include("Functions/steady_state.jl")

country = "Italy"
LOCAL_DIR = "$(@__DIR__)/Results/LAMBDA_grid_big_grid/$(country)_SS_2092_69/General/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\LAMBDA_grid_big_grid\\$(country)_SS_2092_69\\General\\"
end
@load "$(LOCAL_DIR)SSS.jld2" SSS

# global parameters of the model's code
#                   1           2           3       4
#                gen_tol_x, gen_tol_f, distr_tol, val_tol
GLOBAL_PARAMS = SSS[1][1][46]
# global parameters of the approximation objects
#                               1               2               3               4                       5                   6
#                       number_a_nodes, number_u_m_nodes, number_u_w_nodes, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes
GLOBAL_APPROX_PARAMS = SSS[1][1][47]

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

R_MIN = r_min
R_MAX = r_max

# wage bounds #### NOT DONE!!! ########
w_min = 0.01#0.18
w_max = 0.47#0.3

W_MIN = w_min
W_MAX = w_max

optimal_r = R_#(r_min+r_max)/2#r#
optimal_w = W_#(w_min+w_max)/2#w#

text_output = false
fig_output = false
calc_add_results = false

country = "Italy"
guess_R = -6.323809959525245/100#(r_min+r_max)/2#
guess_W = 0.19983052719342714#(w_min+w_max)/2#

# Generate gen eq'm for different lambdas
lambda = 1.513028#,2.6,1.6733333333333333

MODEL_PARAMS = zeros(19)
if country == "Brazil"
    #           1       2     3     4       5   6       7   8       9       10          11              12      13          14          15              16          17              18              19       20
    #         lambda, beta, delta, gamma, eta, theta, c_e, rho_m, rho_w, sigma_eps_m, rho_eps_m_w, sigma_zeta, p_alpha, eta_alpha, prob_node1_alpha, mu_m_alpha, rho_alpha_m_w, sigma_alpha_w, sigma_eps_w, crra
    # brazil MODEL_PARAMS
    MODEL_PARAMS = [lambda, 0.886291098471, 0.06, 0.197999, 0.325612, 0.476388, 0.080459934034, 0.78839752660496, 0.96, 1.14553769575612, 0.3351447009, 0.211613239161, 0.041225789925, 5.403585195245, 0.22878505546, -2.975991405139, 0.147638317002, 0.42, 0.07352027977526, CRRA]

elseif country == "Italy"
    # italy MODEL_PARAMS
    MODEL_PARAMS = [lambda, BETA, DELTA, GAMMA, ETA, THETA, C_E, RHO_M, RHO_W, SIGMA_EPS_M, RHO_EPS_M_W, SIGMA_ZETA, P_ALPHA, ETA_ALPHA, PROB_NODE1_ALPHA, MU_M_ALPHA, RHO_ALPHA_M_W, SIGMA_ALPHA_W, SIGMA_EPS_W, CRRA]
else
    throw(error)
end

# ss_star = [res, r, w, approx_object, MODEL_PARAMS]
@time ss_star = steady_state(guess_R, guess_W, GLOBAL_PARAMS,GLOBAL_APPROX_PARAMS, MODEL_PARAMS)
SS = copy(ss_star)

LOCAL_DIR = "$(@__DIR__)/Results/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\"
end

open("$(LOCAL_DIR)res_$(round(lambda;digits=2)).txt", "a") do f
    write(f, "Lambda: $(lambda)\n")

    write(f, "Interest rate: $(ss_star[2])\n")
    write(f, "Wage: $(ss_star[3])\n")
    write(f, "\n")
    write(f, "\nCapital-to-output: $(ss_star[1][12])\n")
    write(f, "Credit-to-output: $(ss_star[1][13])\n")

    write(f, "Share of workers: $(ss_star[1][14])\n")
    write(f, "Share of self-employed: $(ss_star[1][15])\n")
    write(f, "Share of employers: $(ss_star[1][16])\n")

    write(f, "Variance of log-consumption: $(ss_star[1][17])\n")

    write(f, "Occupational transition probabilities:\n")
    write(f, "$(ss_star[1][18])\n")

    write(f, "Variance of log-earnings: $(ss_star[1][19])\n")

    write(f, "Gini for worker's income: $(ss_star[1][20])\n")
    write(f, "Gini for entrepreneurial's income: $(ss_star[1][21])\n")

    write(f, "\n------Productivity results---------\n")
    write(f, "TFP_ideal: $(ss_star[1][43][end-4]*100)%\n")
    write(f, "TFP_data: $(ss_star[1][43][end-3]*100)%\n")
    write(f, "Labour productivity: $(ss_star[1][43][end-2])\n")
    write(f, "Capital productivity (Output-to-capital): $(ss_star[1][43][end-1])\n")
    write(f, "Capital productivity (Expenses-to-revenue = (r+delta)*Capital-earnings): $(ss_star[1][43][end])\n")
end

@save "$(LOCAL_DIR)SS_lambda_$(round(lambda;digits=2)).jld2" SS
