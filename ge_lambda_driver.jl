include("Functions/print_sameline.jl")
############################
######### Libraries ########
print_sameline("Initialising libraries")
using Plots
using JLD2

print_sameline("Loading functions for steady_state procedure")
include("Functions/steady_state.jl")

# global parameters of the model's code
#                   1           2           3       4
#                gen_tol_x, gen_tol_f, distr_tol, val_tol
GLOBAL_PARAMS = [1e-6,      1e-4,      1e-9,      1e-7]

# global parameters of the approximation objects
#                               1               2               3               4                       5                   6
#                       number_a_nodes, number_u_m_nodes, number_u_w_nodes, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes
GLOBAL_APPROX_PARAMS = [69,             3,                3,                3,                 6,                    3]

# parameters of the model's economy (Italy)
LAMBDA =            1.665907
BETA =              0.934172
DELTA =             0.1
GAMMA =             0.16
ETA =               0.3
THETA =             0.54

C_E =               0.057572
RHO_M =             0.951032
RHO_W =             0.96
SIGMA_EPS_M =       0.754644
RHO_EPS_M_W =       0.216037
SIGMA_ZETA =        0.211565
P_ALPHA =           0.022080
ETA_ALPHA =         5.811008
PROB_NODE1_ALPHA =  0.39
MU_M_ALPHA =       -3.266806
RHO_ALPHA_M_W =     0.147661
SIGMA_ALPHA_W =     0.12
SIGMA_EPS_W =       0.048689

CRRA = 1.0 # >0.0

R_ =                1.434096/100
W_ =                0.294537

gen_tol_x = GLOBAL_PARAMS[1]
gen_tol_f = GLOBAL_PARAMS[2]

# interest rate bounds
r_min = -DELTA#-delta
r_max = 1/BETA-1#1/beta-1

R_MIN = r_min
R_MAX = r_max

# wage bounds #### NOT DONE!!! ########
W_MIN = 0.01#0.23#
W_MAX = 0.47#0.33#

optimal_r = R_#(r_min+r_max)/2#r#
optimal_w = W_#(w_min+w_max)/2#w#

text_output = false
fig_output = false
calc_add_results = false

# Generate gen eq'm for different lambdas
lambdas = [collect(range(1.01;  length = 5, stop = 1.270437995455));
           collect(range(1.27;  length = 6, stop = LAMBDA)        )[2:end];
           collect(range(LAMBDA;length = 4, stop = 2.0)           )[2:end];
           collect(range(2.0;   length = 6, stop = 5.0)           )[2:end];
           collect(range(5.0;   length = 6, stop = 10.0)          )[2:end]#=;
           collect(range(10.0;  length = 5, stop = 25.0)          )[2:end];
           collect(range(25.0;  length = 4, stop = 100.0)         )[2:end];
           collect(range(100.0; length = 4, stop = 2500.0)        )[2:end]=#]

l1_i = findfirst(x -> x==LAMBDA, lambdas)

LOCAL_DIR = "$(@__DIR__)/Results/Stationary/GE_lambda/"
if Sys.iswindows()
   LOCAL_DIR = "$(@__DIR__)\\Results\\Stationary\\GE_lambda\\"
end
mkpath(LOCAL_DIR)

guess_R = R_#-9.88337/100#(r_min+r_max)/2#
guess_W = W_#0.2357213#(w_min+w_max)/2#

Rs = zeros(length(lambdas))
Ws = zeros(length(lambdas))
K_Ys = zeros(length(lambdas))
C_Ys = zeros(length(lambdas))
occWs = zeros(length(lambdas))
occSEs = zeros(length(lambdas))
occEMP = zeros(length(lambdas))
occENT = zeros(length(lambdas))
logcs = zeros(length(lambdas))
loges = zeros(length(lambdas))
giniWs = zeros(length(lambdas))
giniEnts = zeros(length(lambdas))

TFPis = zeros(length(lambdas))
TFPds = zeros(length(lambdas))
Labour_productivity_s = zeros(length(lambdas))
Capital_productivity_s = zeros(length(lambdas))
Investment_productivity_s = zeros(length(lambdas))

SSS = Array{Any}(undef,length(lambdas))

for i = [l1_i; l1_i-1:-1:1; l1_i+1:1:length(lambdas)]#1:length(lambdas)#[[l1_i, l3_i, l2_i]; deleteat!(collect(1:length(lambdas)),[l1_i, l2_i, l3_i])]#
    println("\n$(country) - $(i)/$(length(lambdas)) - $(lambdas[i])")
    if i >= l1_i-1 && i <= l1_i+1
        guess_R = R_
        guess_W = W_
    end
    MODEL_PARAMS = [lambdas[i], BETA, DELTA, GAMMA, ETA, THETA, C_E, RHO_M, RHO_W, SIGMA_EPS_M, RHO_EPS_M_W, SIGMA_ZETA, P_ALPHA, ETA_ALPHA, PROB_NODE1_ALPHA, MU_M_ALPHA, RHO_ALPHA_M_W, SIGMA_ALPHA_W, SIGMA_EPS_W, CRRA]
    
    ss_star = []
    try
        # ss_star = [res, r, w, approx_object, model_params]
        @time ss_star = steady_state(guess_R, guess_W, GLOBAL_PARAMS,GLOBAL_APPROX_PARAMS, MODEL_PARAMS)
    catch e
        println_sameline("Calculation has failed")
        try
            try
                ss_star = copy(SSS[i-1])
            catch eee
                ss_star = copy(SSS[i+1])
            end
        catch ee
            nothing
        end
    end
    SSS[i] = copy(ss_star)

    SS = copy(ss_star)
    @save "$(LOCAL_DIR)SS_lambda_$(round(lambdas[i];digits=2)).jld2" SS

    open("$(LOCAL_DIR)res_$(round(lambdas[i];digits=3)).txt", "a") do f
        write(f, "Lambda: $(lambdas[i])\n")

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

    Rs[i] = ss_star[2]
    Ws[i] = ss_star[3]
    K_Ys[i] = ss_star[1][12]
    C_Ys[i] = ss_star[1][13]
    occWs[i] = ss_star[1][14]
    occSEs[i] = ss_star[1][15]
    occEMP[i] = ss_star[1][16]
    occENT[i] = occSEs[i] + occEMP[i]
    logcs[i] = ss_star[1][17]
    loges[i] = ss_star[1][19]
    giniWs[i] = ss_star[1][20]
    giniEnts[i] = ss_star[1][21]

    TFPis[i] = ss_star[1][43][end-4]
    TFPds[i] = ss_star[1][43][end-3]
    Labour_productivity_s[i] = ss_star[1][43][end-2]
    Capital_productivity_s[i] = ss_star[1][43][end-1]
    Investment_productivity_s[i] = ss_star[1][43][end]

    guess_R = ss_star[2]
    guess_W = ss_star[3]
end

@save "$(LOCAL_DIR)SSS.jld2" GLOBAL_PARAMS GLOBAL_APPROX_PARAMS SSS
