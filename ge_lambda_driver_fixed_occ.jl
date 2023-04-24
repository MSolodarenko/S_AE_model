include("Functions/print_sameline.jl")
############################
######### Libraries ########
print_sameline("Initialising libraries")
using Plots
using JLD2
using ProgressMeter
using Suppressor

print_sameline("Loading functions for steady_state procedure")
include("Functions/steady_state.jl")

country = "Italy"
LOCAL_DIR = "$(@__DIR__)/Results/Lambda_grid_big_grid/$(country)_SS_2092_69/General/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Lambda_grid_big_grid\\$(country)_SS_2092_69\\General\\"
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

optimal_r = R_
optimal_w = W_

text_output = false
fig_output = false
calc_add_results = false

include("Functions/AllubErosa_fixed_occ.jl")
OCC_SHARES = [SSS[10][1][14], SSS[10][1][15], SSS[10][1][16]]
function AllubErosa(r,w,global_params,global_approx_params,model_params,approx_object)
    return AllubErosa(r,w,global_params,global_approx_params,model_params,approx_object,OCC_SHARES)
end

# Generate gen eq'm for different lambdas
lambdas = zeros(length(SSS))
@showprogress for i = 1:length(SSS)
    lambdas[i] = SSS[i][5][1]
end
l1_i = findfirst(x -> x==LAMBDA, lambdas)

#reuse the calculated results
country = "Italy"
LOCAL_DIR = "$(@__DIR__)/Results/Fixed_occ_shares/Lambda_grid/Italy/General/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Fixed_occ_shares\\Lambda_grid\\Italy\\General\\"
end
print_sameline("Loading data from Fixed_occ_shares/Lambda_grid")
@load "$(LOCAL_DIR)SSS_fixed.jld2" SSS_fixed_occ
lambdas_old = zeros(length(SSS_fixed_occ))
@showprogress for i = 1:length(SSS_fixed_occ)
    lambdas_old[i] = SSS_fixed_occ[i][5][1]
end

LOCAL_DIR_ITALY = "$(@__DIR__)/Results/Fixed_occ_shares/Lambda_grid/Italy_updated/"
if Sys.iswindows()
   LOCAL_DIR_ITALY = "$(@__DIR__)\\Results\\Fixed_occ_shares\\Lambda_grid\\Italy_updated\\"
end
for country = ["Italy"]
    if country == "Italy"
        LOCAL_DIR = LOCAL_DIR_ITALY
    else
        throw("No directory set for country $country")
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

    #list_of_iterators = [l1_i; l1_i-1:-1:1; l1_i+1:1:length(lambdas)]
    list_of_iterators = [l1_i; l1_i+1:1:length(lambdas); 1:l1_i-1]
    if Sys.iswindows()
        list_of_iterators = [l1_i; l1_i-1:-1:1; length(lambdas):-1:l1_i+1]
    end

    for i = list_of_iterators
        global R_MIN, W_MIN, R_MAX, W_MAX
        println("\n$(country) - $(i)/$(length(lambdas)) - $(lambdas[i])")

        MODEL_PARAMS = zeros(19)
        if country == "Italy"
            # italy MODEL_PARAMS
            MODEL_PARAMS = [lambdas[i], BETA, DELTA, GAMMA, ETA, THETA, C_E, RHO_M, RHO_W, SIGMA_EPS_M, RHO_EPS_M_W, SIGMA_ZETA, P_ALPHA, ETA_ALPHA, PROB_NODE1_ALPHA, MU_M_ALPHA, RHO_ALPHA_M_W, SIGMA_ALPHA_W, SIGMA_EPS_W, CRRA]
        else
            throw(error)
        end
        ss_star = []
        SS = []
        try
            @suppress_err @load "$(LOCAL_DIR)SS_lambda_$(round(lambdas[i];digits=2)).jld2" SS
            println("Already calculated")
            ss_star = copy(SS)
            println("The interest rate - $(ss_star[2]*100)% and the wage - $(ss_star[3])")
        catch e
            # iter = findmin(abs.(lambdas_old.-lambdas[i]))[2]
            # guess_R = SSS_fixed_occ[iter][2]
            # guess_W = SSS_fixed_occ[iter][3]

            # try
            #     # throw(error)
            #     @load "$(LOCAL_DIR)SS_lambda_$(round(lambdas[i-1];digits=2)).jld2" SS
            #     R_MIN = SS[2]
            #     W_MIN = SS[3]
            # catch e
            #     if i > l1_i
            #         iter = l1_i
            #         try
            #             while iter < i
            #                 @load "$(LOCAL_DIR)SS_lambda_$(round(lambdas[iter];digits=2)).jld2" SS
            #                 iter+=1
            #             end
            #         catch e
            #             @load "$(LOCAL_DIR)SS_lambda_$(round(lambdas[iter-1];digits=2)).jld2" SS
            #             R_MIN = SS[2]
            #             W_MIN = SS[3]
            #         end
            #     else
            #         R_MIN = r_min
            #         W_MIN = w_min
            #     end
            # end
            min_is_found = false
            iter = i-1
            while iter > 0 && !min_is_found
                @suppress_err try
                    @load "$(LOCAL_DIR)SS_lambda_$(round(lambdas[iter];digits=2)).jld2" SS
                    min_is_found = true
                    R_MIN = SS[2]
                    W_MIN = SS[3]
                catch e
                    iter -= 1
                end
            end
            if !min_is_found
                R_MIN = r_min
                W_MIN = w_min
            end

            # try
            #     # throw(error)
            #     @load "$(LOCAL_DIR)SS_lambda_$(round(lambdas[i+1];digits=2)).jld2" SS
            #     R_MAX = SS[2]
            #     W_MAX = SS[3]
            # catch e
            #     if i < l1_i
            #         iter = l1_i
            #         try
            #             while iter > i
            #                 @load "$(LOCAL_DIR)SS_lambda_$(round(lambdas[iter];digits=2)).jld2" SS
            #                 iter-=1
            #             end
            #         catch e
            #             @load "$(LOCAL_DIR)SS_lambda_$(round(lambdas[iter+1];digits=2)).jld2" SS
            #             R_MAX = SS[2]
            #             W_MAX = SS[3]
            #         end
            #     else
            #         R_MAX = r_max
            #         W_MAX = w_max
            #     end
            # end
            max_is_found = false
            iter = i+1
            while iter <= length(lambdas) && !max_is_found
                @suppress_err try
                    @load "$(LOCAL_DIR)SS_lambda_$(round(lambdas[iter];digits=2)).jld2" SS
                    max_is_found = true
                    R_MAX = SS[2]
                    W_MAX = SS[3]
                catch e
                    iter += 1
                end
            end
            if !max_is_found
                R_MAX = r_max
                W_MAX = w_max
            end

            if i == 10
                guess_R = 1.7929539540888162/100#R_
                guess_W = 0.15997672649914652#W_
            elseif i == 8
                guess_R = 1.169/100#R_
                guess_W = 0.159688#W_
            elseif i == 9
                guess_R = 1.5995147810553203/100#R_
                guess_W = 0.15997522673188153#W_
            elseif i == 11
                guess_R = 2.078678700816376/100#R_
                guess_W = 0.1619210016023443#W_
            elseif i == 12
                guess_R = 2.3116713029092075/100#R_
                guess_W = 0.16229482981762994#W_
            elseif i == 13
                guess_R = 2.4904744014401006/100#R_
                guess_W = 0.1642226728241757#W_
            elseif i == 14
                guess_R = 2.9391361991653033/100#R_
                guess_W = 0.16422267358861903#W_
            elseif i == 15
                guess_R = 2.9869593599518462/100#R_
                guess_W = 0.18189517491443577#W_
            elseif i == 16
                guess_R = 2.9945/100#R_
                guess_W = 0.18309#W_
            elseif i > l1_i
                guess_R = 0.95*R_MIN+0.05*R_MAX
                guess_W = 0.99*W_MIN+0.01*W_MAX
            elseif i < l1_i
                guess_R = 0.05*R_MIN+0.95*R_MAX
                guess_W = 0.01*W_MIN+0.99*W_MAX
            end
            println("R: $(R_MIN)<=$(guess_R)<=$(R_MAX)")
            println("W: $(W_MIN)<=$(guess_W)<=$(W_MAX)")

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
            if SS[6] <= GLOBAL_PARAMS[2]
                @save "$(LOCAL_DIR)SS_lambda_$(round(lambdas[i];digits=2)).jld2" SS
            else
                @save "$(LOCAL_DIR)SS_lambda_$(round(lambdas[i];digits=2))_failed.jld2" SS
            end
        end

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

    #=
    asset_grid = SSS[1][1][3]
    number_non_zero_asset_grid = SSS[1][1][43][3]
    #policy = SSS[l][1][4]
    number_u_nodes = SSS[1][4][7]
    number_zeta_nodes = global_approx_params[4]
    number_alpha_m_nodes = global_approx_params[5]
    number_alpha_w_nodes = global_approx_params[6]
    labels = permutedims(["l$(l),u$(u),z$(zeta),am$(alpha_m),aw$(alpha_w)" for l in 1:length(lambdas) for u in [Int64(round((1+number_u_nodes)/2;digits=0))] for zeta in [Int64(round((1+number_zeta_nodes)/2;digits=0))] for alpha_m in [1,number_alpha_m_nodes] for alpha_w in [1,number_alpha_w_nodes]])
    p2 = plot(asset_grid[1:number_non_zero_asset_grid],[SSS[l][1][4][1:number_non_zero_asset_grid,u,zeta,alpha_m,alpha_w] for l in 1:length(lambdas) for u in [Int64(round((1+number_u_nodes)/2;digits=0))] for zeta in [Int64(round((1+number_zeta_nodes)/2;digits=0))] for alpha_m in [1,number_alpha_m_nodes] for alpha_w in [1,number_alpha_w_nodes]], label=labels, legend=:outertopleft)
    display(plot(p2, title="Policy functions (high vs low fixed effects)"))
    savefig("$(LOCAL_DIR)_policy_hvsl_fixed_lambda.png")
    =#

end
