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
GLOBAL_APPROX_PARAMS = [69,3,3,3,6,3]#[35,3,3,3,6,3]

CODENAME = #="SS_2642"#"SS_2065"=#"SS_2092"
CODENAME = "$(CODENAME)_$(GLOBAL_APPROX_PARAMS[1])"
# parameters of the model's economy (Italy)
                    #SS_2642 #SS_2065 #SS_2092 #prev.calibration
LAMBDA =            #=1.633951#1.405096=#1.665907 #1.513028
BETA =              #=0.923514#0.910782=#0.934172 #0.917506
DELTA =             0.1
GAMMA =             0.16
ETA =               0.3
THETA =             0.54

C_E =               #=0.034618#0.011300=#0.057572 #0.007001
RHO_M =             #=0.961272#0.963735=#0.951032 #0.710266
RHO_W =             0.96
SIGMA_EPS_M =       #=0.996198#0.877067=#0.754644 #1.022857
RHO_EPS_M_W =       #=0.296954#0.309361=#0.216037 #0.361605
SIGMA_ZETA =        #=0.256800#0.247935=#0.211565 #0.091089
P_ALPHA =           #=0.039732#0.034184=#0.022080 #0.024445
ETA_ALPHA =         #=5.860434#5.567285=#5.811008 #5.896035
PROB_NODE1_ALPHA =  0.39
MU_M_ALPHA =       #=-4.950012#-2.089921=#-3.266806#-4.817675
RHO_ALPHA_M_W =     #=0.103225#0.166184=#0.147661 #0.181454
SIGMA_ALPHA_W =     0.12
SIGMA_EPS_W =       #=0.119669#0.152209=#0.048689 #0.021947

CRRA = 1.0 # >0.0

R_ =                #=-2.627421/100#-5.147463/100=#1.434096/100#-0.005338
W_ =                #=0.136814#0.287390=#0.294537 #0.230294

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

country = "Italy"#"Brazil"
LAMBDA_GE_DIR = "$(@__DIR__)/Results/Lambda_grid_big_grid/$(country)_$(CODENAME)/"
if Sys.iswindows()
   LAMBDA_GE_DIR = "$(@__DIR__)\\Results\\Lambda_grid_big_grid\\$(country)_$(CODENAME)\\"
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

LOCAL_DIR = "$(@__DIR__)/Results/Transitional_dynamics_big_grid/$(country)_$(CODENAME)/"
if Sys.iswindows()
   LOCAL_DIR = "$(@__DIR__)\\Results\\Transitional_dynamics_big_grid\\$(country)_$(CODENAME)\\"
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
